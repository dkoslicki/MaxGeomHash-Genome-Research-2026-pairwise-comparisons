#!/usr/bin/env python3
"""
06_minhash_pairwise.py
=======================
Compute pairwise Jaccard similarity and containment for candidate GTDB genome
pairs using MinHash sketches (kmer-sketch binary, --algo minhash).

Only genome pairs that passed the FracMinHash 0.01 max-containment threshold
(listed in --candidates CSV) are evaluated.  Both genomes in a pair must have
a sketch in --sketch-dir; pairs missing either sketch are skipped automatically.
This means the same script works whether you sketched 200 genomes (test mode)
or all 143,614 (full run) — it adapts to however many sketches are present.

For each candidate pair (A, B) the script collects two similarity estimates:

  (A as query, B as reference):
    - Jaccard(A, B) via `filter --metric jaccard`
    - containment(B in A) via `filter --metric containment`
        = fraction of B's sketch elements also found in A

  (B as query, A as reference):
    - containment(A in B) via `filter --metric containment`
        = fraction of A's sketch elements also found in B

Combined:
    max_containment = max(containment(A in B), containment(B in A))

The `filter` binary is invoked once per unique genome (as query) with all its
candidate neighbors as references — one call per metric, so two filter
invocations per query genome total.  All calls run in parallel.

Outputs (in --output directory):
  pairwise_results.npz    — sparse arrays matching the kmc_pairwise format:
                             row, col, jaccard, containment_query_in_match,
                             containment_match_in_query, max_containment,
                             size_query_sketch, size_ref_sketch, n_genomes
  genome_index.json       — maps integer index → genome_id (same convention
                             as kmc_pairwise/genome_index.json)
  pairwise_run_stats.json — timing, counts, disk usage for later comparison

Usage:
  # Test run (after TEST_N=200 bash 03_minhash_sketch.sh):
  python3 06_minhash_pairwise.py \\
      --sketch-dir /scratch/.../minhash_sketches \\
      --candidates /scratch/.../gtdb_pairwise_containment.csv \\
      --output     /scratch/.../minhash_pairwise \\
      --cores      192

  # Full run — identical command; the script automatically uses all available
  # sketches, regardless of how many were created in the sketch step.

Dependencies: numpy (standard); scipy optional (for downstream sparse-matrix use)
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
from collections import defaultdict
from pathlib import Path

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# Locate binaries relative to this script file
_SCRIPT_DIR = Path(__file__).resolve().parent
_FILTER_BIN  = _SCRIPT_DIR / "kmer-sketch" / "bin" / "filter"

_SKETCH_EXT = ".minhash.sketch"


# ---------------------------------------------------------------------------
# Candidate pairs
# ---------------------------------------------------------------------------

def discover_sketches(sketch_dir: str) -> set:
    """Return the set of genome IDs that have a sketch file in sketch_dir."""
    return {
        p.name.replace(_SKETCH_EXT, "")
        for p in Path(sketch_dir).glob(f"*{_SKETCH_EXT}")
    }


def load_candidates(csv_path: str, available: set) -> list:
    """
    Read the FracMinHash pairwise CSV and return canonical (A < B) pairs where
    both genomes have been sketched.  Pairs missing a sketch are quietly skipped
    — this is how test mode works without any extra flags.
    """
    pairs = set()
    n_missing = 0
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            q = row.get("query_name", "")
            m = row.get("match_name", "")
            if q not in available or m not in available:
                n_missing += 1
                continue
            if q == m:
                continue
            pairs.add((min(q, m), max(q, m)))   # canonical: A < B lexicographically
    if n_missing:
        log.info(
            "  Skipped %d candidate rows: sketch missing for at least one genome.",
            n_missing,
        )
    return sorted(pairs)


# ---------------------------------------------------------------------------
# Per-query worker (called via multiprocessing pool)
# ---------------------------------------------------------------------------

def _run_filter_for_query(args: tuple):
    """
    Run the `filter` binary for one query genome against all its candidate
    neighbors, collecting both Jaccard and containment estimates.

    Returns:
        (query_id, neighbor_results) or None on fatal error.

        neighbor_results: dict[neighbor_id -> dict with keys:
            'jaccard'      — Jaccard estimate (float)
            'containment'  — containment of neighbor in query (float)
            'size_query'   — #elements in query sketch (int)
            'size_ref'     — #elements in neighbor sketch (int)
        ]
    """
    query_id, neighbor_ids, sketch_dir, filter_bin, tmp_root = args

    query_sketch = os.path.join(sketch_dir, f"{query_id}{_SKETCH_EXT}")
    if not os.path.exists(query_sketch):
        return None

    # Validate which neighbors actually have sketch files
    ref_paths = {
        nid: os.path.join(sketch_dir, f"{nid}{_SKETCH_EXT}")
        for nid in neighbor_ids
        if os.path.exists(os.path.join(sketch_dir, f"{nid}{_SKETCH_EXT}"))
    }
    if not ref_paths:
        return None

    # Unique temp directory per call
    worker_tmp = tempfile.mkdtemp(prefix=f"mh_{query_id[:12]}_", dir=tmp_root)
    try:
        refs_file = os.path.join(worker_tmp, "refs.txt")
        with open(refs_file, "w") as f:
            for npath in ref_paths.values():
                f.write(npath + "\n")

        path_to_id = {v: k for k, v in ref_paths.items()}
        neighbor_results = defaultdict(dict)

        for metric in ("jaccard", "containment"):
            out_file = os.path.join(worker_tmp, f"{metric}.tsv")
            cmd = [
                filter_bin,
                "--query",         query_sketch,
                "--refs-filelist", refs_file,
                "--metric",        metric,
                "--threshold",     "0.0",
                "--output",        out_file,
            ]
            try:
                subprocess.run(cmd, capture_output=True, check=True)
            except subprocess.CalledProcessError as exc:
                log.warning(
                    "filter failed (query=%s, metric=%s): %s",
                    query_id, metric, exc.stderr.decode()[:300],
                )
                continue

            if not os.path.exists(out_file):
                continue

            with open(out_file) as f:
                next(f, None)   # skip header line
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    ref_path = parts[0].strip()
                    try:
                        score  = float(parts[1])
                        size_q = int(parts[3])
                        size_r = int(parts[4])
                    except (ValueError, IndexError):
                        continue

                    nid = path_to_id.get(ref_path)
                    if nid is None:
                        continue

                    neighbor_results[nid][metric] = score
                    if "size_query" not in neighbor_results[nid]:
                        neighbor_results[nid]["size_query"] = size_q
                        neighbor_results[nid]["size_ref"]   = size_r

        return query_id, dict(neighbor_results)

    finally:
        shutil.rmtree(worker_tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="MinHash pairwise similarity for GTDB candidate pairs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sketch-dir", required=True,
                   help="Directory of *.minhash.sketch files "
                        "(from 03_minhash_sketch.sh)")
    p.add_argument("--candidates", required=True,
                   help="FracMinHash pairwise CSV with query_name/match_name columns "
                        "(gtdb_pairwise_containment.csv)")
    p.add_argument("--output", required=True,
                   help="Output directory (created if absent)")
    p.add_argument("--cores", type=int, default=os.cpu_count(),
                   help="Parallel worker processes")
    p.add_argument("--filter-bin", default=str(_FILTER_BIN),
                   help="Path to the `filter` binary")
    return p.parse_args()


def main():
    args = parse_args()
    t0 = time.time()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    filter_bin = args.filter_bin
    if not Path(filter_bin).exists():
        log.error("filter binary not found: %s", filter_bin)
        sys.exit(1)

    # ---- Discover sketches ----------------------------------------------------
    log.info("Scanning sketch directory: %s", args.sketch_dir)
    available = discover_sketches(args.sketch_dir)
    log.info("  %d sketches found.", len(available))

    # ---- Load candidate pairs -------------------------------------------------
    log.info("Loading candidate pairs from: %s", args.candidates)
    pairs = load_candidates(args.candidates, available)
    log.info("  %d candidate pairs where both genomes are sketched.", len(pairs))

    if not pairs:
        log.error(
            "No pairs to process. Run 03_minhash_sketch.sh first "
            "(with TEST_N set if testing)."
        )
        sys.exit(1)

    # ---- Build adjacency list -------------------------------------------------
    neighbors = defaultdict(set)
    for a, b in pairs:
        neighbors[a].add(b)
        neighbors[b].add(a)

    unique_queries = sorted(neighbors.keys())
    log.info("  %d unique genomes will be processed as queries.", len(unique_queries))

    # ---- Genome index ---------------------------------------------------------
    genome_ids  = sorted(available)
    name_to_idx = {gid: i for i, gid in enumerate(genome_ids)}
    n_genomes   = len(genome_ids)

    index_path = str(out_dir / "genome_index.json")
    with open(index_path, "w") as f:
        json.dump({"genomes": genome_ids, "n_genomes": n_genomes}, f, indent=2)
    log.info("Genome index written to %s  (%d entries)", index_path, n_genomes)

    # ---- Shared temp directory ------------------------------------------------
    tmp_root = tempfile.mkdtemp(prefix="mh_pairwise_tmp_", dir=str(out_dir))

    # ---- Build task list ------------------------------------------------------
    tasks = [
        (qid, list(neighbors[qid]), args.sketch_dir, filter_bin, tmp_root)
        for qid in unique_queries
    ]
    n_tasks = len(tasks)
    log.info("Submitting %d query tasks (2 filter calls each) across %d cores ...",
             n_tasks, args.cores)

    # ---- Run workers ----------------------------------------------------------
    raw_data = {}
    n_done   = 0
    t_report = time.time()

    with mp.Pool(processes=args.cores) as pool:
        for result in pool.imap_unordered(_run_filter_for_query, tasks, chunksize=1):
            n_done += 1
            if result is None:
                continue
            query_id, neighbor_results = result
            for nid, metrics in neighbor_results.items():
                raw_data[(query_id, nid)] = metrics

            if time.time() - t_report > 60:
                elapsed = time.time() - t0
                rate    = n_done / elapsed if elapsed > 0 else 0
                eta_h   = (n_tasks - n_done) / rate / 3600 if rate > 0 else float("inf")
                log.info(
                    "Progress: %d / %d queries | %.1f q/s | ETA %.1f h",
                    n_done, n_tasks, rate, eta_h,
                )
                t_report = time.time()

    shutil.rmtree(tmp_root, ignore_errors=True)

    # ---- Aggregate results per canonical pair ---------------------------------
    log.info("Aggregating results for %d canonical pairs ...", len(pairs))

    rows_, cols_           = [], []
    jac_, cqm_, cmq_, mc_ = [], [], [], []
    sq_, sr_               = [], []
    n_incomplete           = 0

    for a, b in pairs:
        ab = raw_data.get((a, b), {})
        ba = raw_data.get((b, a), {})

        jac_val     = ab.get("jaccard") if ab.get("jaccard") is not None \
                      else ba.get("jaccard")
        cont_b_in_a = ab.get("containment")
        cont_a_in_b = ba.get("containment")

        if jac_val is None or cont_b_in_a is None or cont_a_in_b is None:
            n_incomplete += 1
            continue

        size_a = ab.get("size_query", 0) or ba.get("size_ref", 0)
        size_b = ab.get("size_ref",   0) or ba.get("size_query", 0)

        i = name_to_idx[a]
        j = name_to_idx[b]
        rows_.append(i)
        cols_.append(j)
        jac_.append(jac_val)
        cqm_.append(cont_a_in_b)
        cmq_.append(cont_b_in_a)
        mc_.append(max(cont_a_in_b, cont_b_in_a))
        sq_.append(size_a)
        sr_.append(size_b)

    if n_incomplete:
        log.warning(
            "%d pairs had incomplete filter results and were skipped.",
            n_incomplete,
        )
    log.info("Recording %d complete pairs.", len(rows_))

    # ---- Save NPZ -------------------------------------------------------------
    npz_path = str(out_dir / "pairwise_results.npz")
    np.savez_compressed(
        npz_path,
        row                        = np.array(rows_, dtype=np.int32),
        col                        = np.array(cols_, dtype=np.int32),
        jaccard                    = np.array(jac_,  dtype=np.float32),
        containment_query_in_match = np.array(cqm_,  dtype=np.float32),
        containment_match_in_query = np.array(cmq_,  dtype=np.float32),
        max_containment            = np.array(mc_,   dtype=np.float32),
        size_query_sketch          = np.array(sq_,   dtype=np.int32),
        size_ref_sketch            = np.array(sr_,   dtype=np.int32),
        n_genomes                  = np.array([n_genomes], dtype=np.int32),
    )
    log.info("Results saved to %s", npz_path)

    # ---- Run stats -----------------------------------------------------------
    elapsed_total = time.time() - t0
    npz_size_mb   = os.path.getsize(npz_path) / 1024 ** 2
    sketch_dir_size_bytes = sum(
        p.stat().st_size for p in Path(args.sketch_dir).glob(f"*{_SKETCH_EXT}")
    )

    stats = {
        "script":                   "06_minhash_pairwise.py",
        "algo":                     "minhash",
        "sketch_dir":               args.sketch_dir,
        "candidates_csv":           args.candidates,
        "n_sketches_available":     len(available),
        "n_candidate_pairs_input":  len(pairs),
        "n_pairs_recorded":         len(rows_),
        "n_pairs_skipped":          n_incomplete,
        "n_genomes_in_index":       n_genomes,
        "cores":                    args.cores,
        "wall_clock_seconds":       round(elapsed_total, 1),
        "output_npz_mb":            round(npz_size_mb, 2),
        "sketch_dir_total_bytes":   sketch_dir_size_bytes,
        "sketch_dir_total_mb":      round(sketch_dir_size_bytes / 1024 ** 2, 1),
        "output_dir":               str(out_dir),
    }
    stats_path = str(out_dir / "pairwise_run_stats.json")
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)

    # ---- Summary -------------------------------------------------------------
    h = int(elapsed_total // 3600)
    m = int((elapsed_total % 3600) // 60)
    s = int(elapsed_total % 60)

    print()
    print("============================== Run Summary ==============================")
    print(f"  Wall-clock time         : {h:02d}:{m:02d}:{s:02d} (hh:mm:ss)")
    print(f"  Sketches available      : {len(available)}")
    print(f"  Candidate pairs         : {len(pairs)}")
    print(f"  Pairs recorded          : {len(rows_)}")
    print(f"  Pairs skipped           : {n_incomplete}  (incomplete filter output)")
    print(f"  Output NPZ              : {npz_path}  ({npz_size_mb:.1f} MB)")
    print(f"  Genome index            : {index_path}")
    print(f"  Run stats               : {stats_path}")
    print("=========================================================================")
    print()
    print("Next step (sanity check):")
    print(f"  python3 {_SCRIPT_DIR}/09_sanity_check.py \\")
    print(f"      --minhash-pairwise   {out_dir} \\")
    print(f"      --kmc-pairwise       <kmc_pairwise_dir> \\")
    print(f"      --output             {out_dir}/sanity_check")


if __name__ == "__main__":
    main()
