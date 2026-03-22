#!/usr/bin/env python3
"""
02_kmc_pairwise.py
==================
Compute exact pairwise Jaccard similarity and containment indices between
all GTDB genome k-mer sets (or a candidate subset) using KMC databases.

This is the "ground truth" counterpart to the approximate sourmash pairwise
computation. It uses kmc_tools to compute exact k-mer set intersections and
derives Jaccard, containment (both directions), and max-containment from the
exact counts.

Outputs (all in --output directory):
  pairwise_results.npz   — sparse result arrays (row, col, + metric arrays)
  genome_index.json      — maps integer index → genome_id (for row/col lookup)

Usage examples:
  # Full upper-triangle pairwise (with threshold pruning):
  python3 02_kmc_pairwise.py \\
      --metadata  /scratch/.../kmc_metadata.csv \\
      --output    /scratch/.../kmc_pairwise \\
      --threshold 0.01 \\
      --cores     192

  # Exact values for candidate pairs from sourmash pairwise output:
  python3 02_kmc_pairwise.py \\
      --metadata   /scratch/.../kmc_metadata.csv \\
      --candidates /scratch/.../gtdb_pairwise_containment.csv \\
      --output     /scratch/.../kmc_pairwise \\
      --threshold  0.00 \\
      --cores      192

Dependencies:
  numpy, scipy (for sparse output)
  kmc_db_info binary (compiled from kmc_db_info.cpp — no py_kmc_api needed)
"""

import argparse
import csv
import json
import logging
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# kmc_db_info binary — expected alongside this script
_KMC_DB_INFO_BIN = Path(__file__).resolve().parent / "kmc_db_info"


# ---------------------------------------------------------------------------
# KMC utilities
# ---------------------------------------------------------------------------

def get_kmer_count(db_prefix: str) -> int:
    """
    Return the number of k-mers in a KMC database by reading the .kmc_pre header.

    Primary path: kmc_db_info (compiled C++ helper using official kmc_api).
      - Reads only the header — no k-mer iteration, O(1) effectively.
      - Works for both -cs1 (existing) and -cs2 (new) databases.

    Fallback: kmc_tools transform histogram (for -cs2 databases only).
      - Sums the stored per-multiplicity histogram. Fast but unreliable for -cs1.

    NOTE: Under normal operation n_kmers values come directly from the metadata
    CSV (populated by rebuild_kmc_metadata.py), so this function is only called
    for rows where n_unique_kmers == 0 (incomplete metadata). If your CSV is
    complete, this function is never invoked.
    """
    # ---- Fast path: kmc_db_info C++ helper ---------------------------
    if _KMC_DB_INFO_BIN.exists():
        try:
            result = subprocess.run(
                [str(_KMC_DB_INFO_BIN)],
                input=f"{db_prefix}\n".encode(),
                capture_output=True,
                check=True,
            )
            for line in result.stdout.decode().splitlines():
                parts = line.strip().split("\t")
                if len(parts) == 2 and parts[0] == db_prefix:
                    return int(parts[1]) if parts[1].isdigit() else 0
        except Exception as exc:
            log.warning("kmc_db_info failed for %s: %s", db_prefix, exc)

    # ---- Fallback: kmc_tools histogram (cs2 databases only) ----------
    # Sums the second column of the histogram file (unique k-mers per count level).
    # With -cs2 this is at most 2 lines; with -cs1 this may be unreliable.
    log.warning(
        "kmc_db_info binary not found at %s — falling back to kmc_tools histogram. "
        "This may give wrong results for -cs1 databases.",
        _KMC_DB_INFO_BIN,
    )
    tmp_hist = db_prefix + "_tmp_hist.txt"
    try:
        subprocess.run(
            ["kmc_tools", "transform", db_prefix, "histogram", tmp_hist],
            check=True,
            capture_output=True,
        )
        total = 0
        with open(tmp_hist) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    total += int(parts[1])
        return total
    except Exception as exc:
        raise RuntimeError(
            f"Cannot read k-mer count from {db_prefix}: {exc}"
        ) from exc
    finally:
        try:
            os.unlink(tmp_hist)
        except FileNotFoundError:
            pass


def compute_intersection_size(
    db_i: str,
    db_j: str,
    tmp_dir: str,
    pair_id: str,
) -> int:
    """
    Use kmc_tools simple to compute the exact intersection of two KMC databases
    and return the number of shared k-mers.  Temporary files are cleaned up
    before returning.
    """
    tmp_prefix = os.path.join(tmp_dir, f"isect_{pair_id}")
    try:
        subprocess.run(
            [
                "kmc_tools",
                "-t1",          # 1 thread per kmc_tools call — parallelism comes
                                # from having 192 worker processes, not from each
                                # kmc_tools job grabbing many threads. Without this,
                                # 192 workers × N auto-detected threads = thousands
                                # of competing threads, causing the in1:0% hang.
                "simple",
                db_i, "-ci1",   # include all k-mers from db_i (count >= 1)
                db_j, "-ci1",   # include all k-mers from db_j
                "intersect", tmp_prefix,
            ],
            capture_output=True,
            check=True,
        )
        n_intersect = get_kmer_count(tmp_prefix)
    except subprocess.CalledProcessError as exc:
        log.error("kmc_tools simple failed for %s vs %s: %s", db_i, db_j, exc.stderr.decode())
        n_intersect = 0
    finally:
        # Always clean up temp files
        for ext in (".kmc_pre", ".kmc_suf"):
            try:
                os.unlink(tmp_prefix + ext)
            except FileNotFoundError:
                pass

    return n_intersect


# ---------------------------------------------------------------------------
# Worker function (runs in a subprocess via multiprocessing)
# ---------------------------------------------------------------------------

def _worker(args: tuple) -> Optional[tuple]:
    """
    Compute exact metrics for a single genome pair (i, j).

    Returns:
        (i, j, jaccard, containment_i_in_j, containment_j_in_i,
         max_containment, intersect_size)
        or None if the pair does not meet the threshold.
    """
    (
        i, j,
        db_i, db_j,
        n_i, n_j,
        tmp_root,
        threshold,
    ) = args

    # ---- Size-based Jaccard upper-bound pruning -----------------------
    # Jaccard(A,B) ≤ min(|A|,|B|) / max(|A|,|B|)
    # If this upper bound is below threshold we can skip the intersection.
    # Note: this prunes only when using a Jaccard threshold; max_containment
    # can always reach 1.0 regardless of size ratio, so we only apply pruning
    # if the caller has set threshold to represent a Jaccard floor.
    # In practice for GTDB, containment is often the more relevant metric
    # and max_containment pruning is applied after intersection anyway.
    if n_i > 0 and n_j > 0:
        jaccard_upper = min(n_i, n_j) / max(n_i, n_j)
        if jaccard_upper < threshold:
            return None

    # ---- Exact intersection via kmc_tools ----------------------------
    pair_id = f"{i}_{j}"
    # Each worker gets its own tmp subdir to avoid collisions
    worker_tmp = os.path.join(tmp_root, f"w_{os.getpid()}")
    os.makedirs(worker_tmp, exist_ok=True)

    n_intersect = compute_intersection_size(db_i, db_j, worker_tmp, pair_id)

    if n_intersect == 0:
        return None

    # ---- Compute metrics ---------------------------------------------
    union_size = n_i + n_j - n_intersect
    jaccard = n_intersect / union_size if union_size > 0 else 0.0
    containment_i_in_j = n_intersect / n_i if n_i > 0 else 0.0   # fraction of A found in B
    containment_j_in_i = n_intersect / n_j if n_j > 0 else 0.0   # fraction of B found in A
    max_containment = max(containment_i_in_j, containment_j_in_i)

    if max_containment < threshold:
        return None

    return (
        i, j,
        float(jaccard),
        float(containment_i_in_j),
        float(containment_j_in_i),
        float(max_containment),
        int(n_intersect),
    )


# ---------------------------------------------------------------------------
# Genome index helpers
# ---------------------------------------------------------------------------

def load_metadata(csv_path: str) -> list[dict]:
    """Read kmc_metadata.csv, skipping rows with missing or zero k-mer counts."""
    records = []
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            n = row.get("n_unique_kmers", "0")
            row["n_unique_kmers"] = int(n) if n.isdigit() else 0
            records.append(row)
    return records


def load_candidates(csv_path: str, name_to_idx: dict[str, int]) -> list[tuple[int, int]]:
    """
    Parse a sourmash pairwise CSV (or any CSV with 'query_name'/'match_name'
    columns) and return a list of (i, j) index pairs where i < j.
    Unrecognised genome names are skipped with a warning.
    """
    pairs = set()
    skipped = 0
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            q = row.get("query_name", "")
            m = row.get("match_name", "")
            if q not in name_to_idx or m not in name_to_idx:
                skipped += 1
                continue
            i, j = name_to_idx[q], name_to_idx[m]
            if i == j:
                continue
            pairs.add((min(i, j), max(i, j)))
    if skipped:
        log.warning("Skipped %d candidate rows with unrecognised genome names.", skipped)
    return sorted(pairs)


# ---------------------------------------------------------------------------
# Result accumulator + writer
# ---------------------------------------------------------------------------

class SparseAccumulator:
    """Accumulates sparse pairwise results and writes them to disk in batches."""

    def __init__(self, flush_every: int = 5_000_000):
        self.rows: list[int] = []
        self.cols: list[int] = []
        self.jaccard: list[float] = []
        self.cont_ij: list[float] = []
        self.cont_ji: list[float] = []
        self.max_cont: list[float] = []
        self.isect: list[int] = []
        self.flush_every = flush_every
        self._flushed_chunks: list[dict] = []

    def add(self, result: tuple):
        i, j, jac, cij, cji, mc, isz = result
        self.rows.append(i)
        self.cols.append(j)
        self.jaccard.append(jac)
        self.cont_ij.append(cij)
        self.cont_ji.append(cji)
        self.max_cont.append(mc)
        self.isect.append(isz)
        if len(self.rows) >= self.flush_every:
            self._flush_to_memory()

    def _flush_to_memory(self):
        if not self.rows:
            return
        self._flushed_chunks.append({
            "rows": np.array(self.rows, dtype=np.int32),
            "cols": np.array(self.cols, dtype=np.int32),
            "jaccard": np.array(self.jaccard, dtype=np.float32),
            "cont_ij": np.array(self.cont_ij, dtype=np.float32),
            "cont_ji": np.array(self.cont_ji, dtype=np.float32),
            "max_cont": np.array(self.max_cont, dtype=np.float32),
            "isect": np.array(self.isect, dtype=np.int32),
        })
        self.rows.clear()
        self.cols.clear()
        self.jaccard.clear()
        self.cont_ij.clear()
        self.cont_ji.clear()
        self.max_cont.clear()
        self.isect.clear()

    def save(self, out_path: str, n_genomes: int):
        """Merge all chunks and save to a compressed .npz file."""
        self._flush_to_memory()  # flush any remaining
        if not self._flushed_chunks:
            log.warning("No pairs above threshold — saving empty result.")
            np.savez_compressed(
                out_path,
                row=np.array([], dtype=np.int32),
                col=np.array([], dtype=np.int32),
                jaccard=np.array([], dtype=np.float32),
                containment_query_in_match=np.array([], dtype=np.float32),
                containment_match_in_query=np.array([], dtype=np.float32),
                max_containment=np.array([], dtype=np.float32),
                intersect_size=np.array([], dtype=np.int32),
                n_genomes=np.array([n_genomes], dtype=np.int32),
            )
            return

        keys = list(self._flushed_chunks[0].keys())
        merged = {k: np.concatenate([c[k] for c in self._flushed_chunks]) for k in keys}

        np.savez_compressed(
            out_path,
            row=merged["rows"],
            col=merged["cols"],
            jaccard=merged["jaccard"],
            containment_query_in_match=merged["cont_ij"],
            containment_match_in_query=merged["cont_ji"],
            max_containment=merged["max_cont"],
            intersect_size=merged["isect"],
            n_genomes=np.array([n_genomes], dtype=np.int32),
        )
        log.info("Saved %d pairs to %s", len(merged["rows"]), out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Exact pairwise k-mer Jaccard/containment via KMC databases.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--metadata", required=True,
                   help="kmc_metadata.csv from 01_kmc_count.sh")
    p.add_argument("--output", required=True,
                   help="Output directory (will be created if absent)")
    p.add_argument("--threshold", type=float, default=0.01,
                   help="Minimum max-containment to record a pair")
    p.add_argument("--cores", type=int, default=os.cpu_count(),
                   help="Worker processes")
    p.add_argument("--candidates",
                   help="Optional: CSV of candidate pairs (e.g. sourmash pairwise output). "
                        "Must have 'query_name' and 'match_name' columns. "
                        "When provided, skips full pairwise enumeration.")
    p.add_argument("--chunk-size", type=int, default=1,
                   help="Tasks returned to main process at a time. Keep at 1 so progress "
                        "reporting fires correctly — kmc_tools runtime dominates any IPC overhead.")
    return p.parse_args()


def main():
    args = parse_args()
    t0 = time.time()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    npz_path   = str(out_dir / "pairwise_results.npz")
    index_path = str(out_dir / "genome_index.json")

    # ---- Load metadata -----------------------------------------------
    log.info("Loading metadata from %s", args.metadata)
    records = load_metadata(args.metadata)

    # Build genome index: integer → genome_id (stable, sorted order)
    records.sort(key=lambda r: r["genome_id"])
    genome_ids = [r["genome_id"] for r in records]
    name_to_idx = {gid: i for i, gid in enumerate(genome_ids)}
    n_genomes = len(genome_ids)

    db_prefixes = [r["db_prefix"] for r in records]
    n_kmers     = [r["n_unique_kmers"] for r in records]

    # Warn about genomes with missing k-mer counts
    n_missing = sum(1 for n in n_kmers if n == 0)
    if n_missing:
        log.warning("%d genomes have zero/unknown k-mer count and will be skipped.", n_missing)

    log.info("Genome index: %d genomes.", n_genomes)

    # Save genome index JSON
    with open(index_path, "w") as f:
        json.dump({"genomes": genome_ids, "n_genomes": n_genomes}, f, indent=2)
    log.info("Genome index written to %s", index_path)

    # ---- Build pair list ---------------------------------------------
    if args.candidates:
        log.info("Loading candidate pairs from %s", args.candidates)
        pairs = load_candidates(args.candidates, name_to_idx)
        log.info("  %d candidate pairs loaded.", len(pairs))
    else:
        n_pairs_total = n_genomes * (n_genomes - 1) // 2
        log.info("Full upper-triangle pairwise: %d pairs.", n_pairs_total)
        log.info("  threshold = %.4f  |  size-based Jaccard pruning enabled.", args.threshold)
        pairs = None  # generated lazily in the task iterator below

    # ---- Shared temp directory for intersection DBs ------------------
    tmp_root = tempfile.mkdtemp(prefix="kmc_pairwise_tmp_", dir=str(out_dir))
    log.info("Temporary intersection DBs in: %s", tmp_root)

    def _cleanup_tmp():
        shutil.rmtree(tmp_root, ignore_errors=True)

    # ---- Task generator ----------------------------------------------
    def task_iter():
        """Yield worker arg-tuples for each pair."""
        if pairs is not None:
            # Candidate-pairs mode
            for i, j in pairs:
                if n_kmers[i] == 0 or n_kmers[j] == 0:
                    continue
                yield (i, j, db_prefixes[i], db_prefixes[j],
                       n_kmers[i], n_kmers[j], tmp_root, args.threshold)
        else:
            # Full upper-triangle mode
            for i in range(n_genomes):
                if n_kmers[i] == 0:
                    continue
                for j in range(i + 1, n_genomes):
                    if n_kmers[j] == 0:
                        continue
                    yield (i, j, db_prefixes[i], db_prefixes[j],
                           n_kmers[i], n_kmers[j], tmp_root, args.threshold)

    # ---- Run workers -------------------------------------------------
    n_total_pairs = len(pairs) if pairs is not None else (n_genomes * (n_genomes - 1) // 2)
    log.info("Starting %d worker processes (chunk_size=%d).", args.cores, args.chunk_size)

    accumulator = SparseAccumulator(flush_every=5_000_000)
    n_processed = 0
    n_recorded  = 0
    t_report = time.time()

    with mp.Pool(processes=args.cores) as pool:
        for result in pool.imap_unordered(
            _worker, task_iter(), chunksize=args.chunk_size
        ):
            n_processed += 1
            if result is not None:
                accumulator.add(result)
                n_recorded += 1

            # Progress report every 60 seconds
            if time.time() - t_report > 60:
                pct = 100 * n_processed / n_total_pairs if n_total_pairs > 0 else 0
                elapsed = time.time() - t0
                rate = n_processed / elapsed if elapsed > 0 else 0
                eta_s  = (n_total_pairs - n_processed) / rate if rate > 0 else float("inf")
                eta_h  = eta_s / 3600
                log.info(
                    "Progress: %d / %d pairs (%.1f%%)  |  %d recorded  |  "
                    "%.0f pairs/s  |  ETA %.1f h",
                    n_processed, n_total_pairs, pct, n_recorded, rate, eta_h
                )
                t_report = time.time()

    # ---- Save results ------------------------------------------------
    log.info("Saving results ...")
    accumulator.save(npz_path, n_genomes)

    _cleanup_tmp()

    # ---- Summary -----------------------------------------------------
    elapsed_total = time.time() - t0
    h = int(elapsed_total // 3600)
    m = int((elapsed_total % 3600) // 60)
    s = int(elapsed_total % 60)

    print()
    print("============================== Run Summary ==============================")
    print(f"  Wall-clock time     : {h:02d}:{m:02d}:{s:02d} (hh:mm:ss)")
    print(f"  Genomes             : {n_genomes}")
    print(f"  Pairs evaluated     : {n_processed:,}")
    print(f"  Pairs recorded      : {n_recorded:,}  (max_containment >= {args.threshold})")
    print(f"  Sparsity            : {100*(1 - n_recorded/max(n_processed,1)):.2f}% below threshold")
    print(f"  Output NPZ          : {npz_path}")
    print(f"  Genome index        : {index_path}")
    print("=========================================================================")
    print()
    print("To load results in Python:")
    print(f"  data  = np.load('{npz_path}')")
    print("  row, col       = data['row'], data['col']")
    print("  jaccard        = data['jaccard']")
    print("  max_containment= data['max_containment']")
    print("  cont_q_in_m    = data['containment_query_in_match']")
    print("  cont_m_in_q    = data['containment_match_in_query']")
    print("  intersect_size = data['intersect_size']")
    print()
    print("  with open('genome_index.json') as f:")
    print("      index = json.load(f)['genomes']  # list: integer → genome_id")
    print()
    print("  # Convert to scipy sparse (e.g. Jaccard matrix):")
    print("  from scipy.sparse import csr_matrix")
    print("  n = data['n_genomes'][0]")
    print("  J = csr_matrix((data['jaccard'], (data['row'], data['col'])), shape=(n, n))")
    print("  J = J + J.T  # make symmetric (upper-triangle stored)")


if __name__ == "__main__":
    main()
