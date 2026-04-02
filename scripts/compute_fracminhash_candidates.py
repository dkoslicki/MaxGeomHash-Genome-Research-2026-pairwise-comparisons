#!/usr/bin/env python3
"""
compute_fracminhash_candidates.py
===================================
All-vs-all FracMinHash pairwise to find candidate genome pairs above a
max-containment threshold. Replaces the sourmash branchwater pairwise step,
using the kmer-sketch FracMinHash sketches already computed by
04_fracminhash_sketch.sh (scale=0.01, scaled=100).

Algorithm
---------
Phase 1 — Load sketches (parallel):
    Read all .fracminhash.sketch files using a multiprocessing pool. Each file
    is plain text with '#' header lines and one uint64 hash per data line
    (pre-sorted by the sketch binary). Stored as sorted numpy uint64 arrays.

Phase 2 — Sharded inverted index → pair counts (parallel):
    Divide the hash space [0, max_hash] into NUM_SHARDS equal shards and
    process each shard in a separate worker process.

    Worker: for one shard, binary-search each sketch to extract hashes in
    [lo, hi), concatenate all (hash, genome_idx) pairs, numpy-argsort by hash,
    then sweep runs of equal hash values to emit (i, j, count) triplets for
    genome pairs sharing ≥1 hash.  Returns three compact int32 numpy arrays
    (rows, cols, counts).

    Parallelism: workers share the ~35 GB sketch data read-only via Linux's
    fork copy-on-write (no serialisation of sketch arrays). Only the small
    result arrays travel over the IPC pipe.

    Merge: all per-shard (rows, cols, counts) arrays are concatenated, then
    a single numpy lexsort + np.add.at sums the counts per unique (i, j) pair.

Phase 3 — Filter and write CSV:
    Emit pairs where max_containment = shared / min(size_A, size_B) ≥ threshold.

Output CSV columns (compatible with all downstream scripts)
-----------------------------------------------------------
  query_name, match_name,
  containment_query_in_match, containment_match_in_query,
  max_containment, jaccard,
  intersect_hashes, size_query, size_match

Usage
-----
  python3 compute_fracminhash_candidates.py \\
      --sketch-dir /scratch/.../fracminhash_sketches \\
      --output     /scratch/.../gtdb_pairwise_containment_thr0001.csv \\
      --threshold  0.001 \\
      --cores      384 \\
      --num-shards 128

Performance notes (143 k genomes, ~32 k hashes each, scale=0.01)
-----------------------------------------------------------------
  Phase 1 : I/O-bound; ~35 GB loaded in ~2 min with 384 parallel readers.
  Phase 2 : CPU-bound; 128 shards × ~5 s each / 128 workers ≈ 5–15 s compute.
            Merge of per-shard numpy arrays adds ~1–5 min depending on pair count.
  Phase 3 : Fast (CSV write).
  Peak RSS: ~50–70 GB (sketches + shard result arrays during merge).
"""

import argparse
import csv
import json
import logging
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

_SKETCH_EXT = ".fracminhash.sketch"

# ---------------------------------------------------------------------------
# Module-level global: populated in the parent process BEFORE pool creation
# so that forked workers inherit it via copy-on-write (Linux default).
# Workers never write to this; refcount writes on Python objects cause only
# ~15 MB of CoW overhead per worker, not the 35 GB of underlying C data.
# ---------------------------------------------------------------------------
_SKETCHES_GLOBAL = None   # list[np.ndarray[uint64]]


# ---------------------------------------------------------------------------
# Phase 1: sketch parsing (pool worker)
# ---------------------------------------------------------------------------

def _parse_one_sketch(args: tuple):
    """Load one .fracminhash.sketch file; return (idx, genome_id, sorted uint64 array)."""
    idx, path = args
    genome_id = Path(path).name.replace(_SKETCH_EXT, "")
    hashes = []
    try:
        with open(path, "r") as fh:
            for line in fh:
                if line and line[0] != "#":
                    s = line.strip()
                    if s:
                        hashes.append(int(s))
    except Exception as exc:
        log.warning("Failed to parse %s: %s", path, exc)
        return idx, genome_id, np.array([], dtype=np.uint64)
    arr = np.array(hashes, dtype=np.uint64)
    if len(arr) > 1 and not np.all(arr[:-1] <= arr[1:]):
        arr.sort()
    return idx, genome_id, arr


# ---------------------------------------------------------------------------
# Phase 2: shard processing (pool worker — uses _SKETCHES_GLOBAL via fork)
# ---------------------------------------------------------------------------

def _process_shard(args: tuple):
    """
    Process one hash-space shard.

    Extracts all (hash, genome_idx) entries for hashes in [lo, hi) across
    every sketch, sorts by hash, sweeps runs of equal hashes, and returns
    three int32 numpy arrays: rows, cols, shared_counts — one entry per
    (i, j) genome pair (i < j) that shares ≥1 hash in this shard.
    Returns None if no pairs found.
    """
    shard_idx, lo, hi = args
    sketches = _SKETCHES_GLOBAL   # read-only access; no CoW on the C data

    lo = np.uint64(lo)
    hi = np.uint64(hi)

    hash_chunks = []
    gidx_chunks = []

    for gidx, sketch in enumerate(sketches):
        if len(sketch) == 0:
            continue
        start = int(np.searchsorted(sketch, lo))
        end   = int(np.searchsorted(sketch, hi))
        if start >= end:
            continue
        hash_chunks.append(sketch[start:end])
        gidx_chunks.append(np.full(end - start, gidx, dtype=np.int32))

    if not hash_chunks:
        return None

    all_hashes = np.concatenate(hash_chunks)
    all_gidxs  = np.concatenate(gidx_chunks)

    # Sort by hash value
    order      = np.argsort(all_hashes, kind="stable")
    all_hashes = all_hashes[order]
    all_gidxs  = all_gidxs[order]

    # Sweep runs of identical hashes and collect (i, j, 1) events
    diffs       = np.diff(all_hashes)
    boundaries  = np.where(diffs != 0)[0] + 1   # start of each new hash group
    prev        = 0

    shard_rows   = []
    shard_cols   = []

    for boundary in boundaries.tolist() + [len(all_hashes)]:
        run_len = boundary - prev
        if run_len >= 2:
            run = sorted(all_gidxs[prev:boundary].tolist())
            m   = len(run)
            if m == 2:
                shard_rows.append(run[0])
                shard_cols.append(run[1])
            else:
                for a in range(m):
                    for b in range(a + 1, m):
                        shard_rows.append(run[a])
                        shard_cols.append(run[b])
        prev = boundary

    if not shard_rows:
        return None

    rows_arr = np.array(shard_rows, dtype=np.int32)
    cols_arr = np.array(shard_cols, dtype=np.int32)
    return (rows_arr, cols_arr)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="All-vs-all FracMinHash candidate pair detection.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sketch-dir", required=True,
                   help="Directory of *.fracminhash.sketch files")
    p.add_argument("--output", required=True,
                   help="Output CSV file path")
    p.add_argument("--threshold", type=float, default=0.001,
                   help="Min max-containment to emit")
    p.add_argument("--cores", type=int, default=os.cpu_count(),
                   help="Parallel workers (used for both Phase 1 and Phase 2)")
    p.add_argument("--num-shards", type=int, default=128,
                   help="Number of hash-space shards (≥ --cores recommended so "
                        "each core gets ≥1 shard)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args    = parse_args()
    t_start = time.time()

    sketch_dir = Path(args.sketch_dir)
    paths = sorted(sketch_dir.glob(f"*{_SKETCH_EXT}"))
    n = len(paths)

    log.info("Found %d sketch files in %s", n, sketch_dir)
    if n < 2:
        log.error("Need at least 2 sketches. Run 04_fracminhash_sketch.sh first.")
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Phase 1: load all sketches in parallel
    # -----------------------------------------------------------------------
    log.info("Phase 1: loading %d sketch files (%d workers) ...", n, args.cores)
    t1 = time.time()

    sketches   = [None] * n
    genome_ids = [None] * n
    tasks      = [(i, str(p)) for i, p in enumerate(paths)]
    chunk      = max(1, n // (args.cores * 4))

    with mp.Pool(processes=min(args.cores, n)) as pool:
        for idx, genome_id, arr in pool.imap_unordered(
            _parse_one_sketch, tasks, chunksize=chunk
        ):
            sketches[idx]   = arr
            genome_ids[idx] = genome_id

    sizes        = np.array([len(s) for s in sketches], dtype=np.int64)
    total_hashes = int(sizes.sum())
    log.info(
        "  Loaded in %.1f s.  Genomes: %d  |  Total hashes: %d  |  "
        "Avg sketch size: %.0f  |  Min: %d  Max: %d",
        time.time() - t1, n, total_hashes,
        total_hashes / n, int(sizes.min()), int(sizes.max()),
    )

    # -----------------------------------------------------------------------
    # Phase 2: parallel sharded inverted index
    #
    # Critical: set the module-level global BEFORE creating the Pool.
    # On Linux, Pool() uses fork, so child processes inherit this reference
    # via copy-on-write. The 35 GB of C-allocated numpy data is never copied.
    # -----------------------------------------------------------------------
    max_hash_val = int(max(
        int(s[-1]) if len(s) > 0 else 0 for s in sketches
    ))
    if max_hash_val == 0:
        log.error("All sketches are empty.")
        sys.exit(1)

    NUM_SHARDS  = args.num_shards
    shard_width = (max_hash_val + NUM_SHARDS) // NUM_SHARDS  # integer, plain Python

    # Expose sketches to forked workers via module global (no serialisation).
    # We explicitly request the 'fork' start method so workers are cloned from
    # the parent mid-execution and inherit _SKETCHES_GLOBAL directly.
    # Python 3.14 on Linux may no longer default to fork (e.g. it switches to
    # spawn/forkserver when numpy's internal threads are detected), which would
    # leave workers with the module-initialisation value of None.
    global _SKETCHES_GLOBAL
    _SKETCHES_GLOBAL = sketches

    fork_ctx = mp.get_context("fork")

    n_workers = min(args.cores, NUM_SHARDS)
    log.info(
        "Phase 2: %d shards × %d parallel workers  |  shard width ≈ %d",
        NUM_SHARDS, n_workers, shard_width,
    )
    t2 = time.time()

    shard_tasks = [
        (s, s * shard_width, min((s + 1) * shard_width, max_hash_val + 1))
        for s in range(NUM_SHARDS)
    ]

    # Collect per-shard (rows, cols) arrays; counts are implicit (each row = 1 event)
    all_rows_parts = []
    all_cols_parts = []
    n_done = 0

    with fork_ctx.Pool(processes=n_workers) as pool:
        for result in pool.imap_unordered(_process_shard, shard_tasks, chunksize=1):
            n_done += 1
            if result is not None:
                rows_arr, cols_arr = result
                all_rows_parts.append(rows_arr)
                all_cols_parts.append(cols_arr)
            if n_done % max(1, NUM_SHARDS // 16) == 0 or n_done == NUM_SHARDS:
                elapsed_s2 = time.time() - t2
                rate = n_done / elapsed_s2
                eta  = (NUM_SHARDS - n_done) / rate if rate > 0 else 0
                n_events = sum(len(r) for r in all_rows_parts)
                log.info(
                    "  Shard %4d / %d  |  pair-events so far: %d  |  "
                    "%.1f shards/s  |  ETA %.0f s",
                    n_done, NUM_SHARDS, n_events, rate, eta,
                )

    log.info("  Shard processing done in %.1f s.", time.time() - t2)

    # -----------------------------------------------------------------------
    # Merge: concatenate all (row, col) events, sum counts per unique pair
    # -----------------------------------------------------------------------
    log.info("  Merging shard results ...")
    t_merge = time.time()

    if not all_rows_parts:
        log.warning("No pairs found at all — output will be empty.")
        all_rows   = np.array([], dtype=np.int32)
        all_cols   = np.array([], dtype=np.int32)
        pair_rows  = np.array([], dtype=np.int32)
        pair_cols  = np.array([], dtype=np.int32)
        pair_counts = np.array([], dtype=np.int32)
    else:
        all_rows = np.concatenate(all_rows_parts)
        all_cols = np.concatenate(all_cols_parts)
        del all_rows_parts, all_cols_parts   # free memory

        # Encode each (row, col) pair as a single int64 key: row * n + col
        # (row < n and col < n, so no overflow for n ≤ 143 614 < 2^17)
        composite = all_rows.astype(np.int64) * n + all_cols.astype(np.int64)
        del all_rows, all_cols

        # np.unique sorts internally (pdqsort on int64 — faster than a
        # separate stable mergesort pass + random-access fancy-index reorder).
        uniq_comp, uniq_counts = np.unique(composite, return_counts=True)
        del composite

        pair_rows   = (uniq_comp // n).astype(np.int32)
        pair_cols   = (uniq_comp  % n).astype(np.int32)
        pair_counts = uniq_counts.astype(np.int32)
        del uniq_comp, uniq_counts

    n_pairs_any = len(pair_rows)
    log.info(
        "  Merge done in %.1f s.  Unique pairs sharing ≥1 hash: %d",
        time.time() - t_merge, n_pairs_any,
    )

    # -----------------------------------------------------------------------
    # Phase 3: filter by max-containment and write CSV
    # -----------------------------------------------------------------------
    log.info(
        "Phase 3: filtering at max-containment ≥ %.4f and writing CSV ...",
        args.threshold,
    )
    t3 = time.time()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sizes_i = sizes[pair_rows]
    sizes_j = sizes[pair_cols]
    min_sizes = np.minimum(sizes_i, sizes_j)

    # Avoid division by zero (shouldn't happen with valid sketches)
    valid = min_sizes > 0
    max_cont_arr = np.where(valid, pair_counts / min_sizes, 0.0)
    keep = max_cont_arr >= args.threshold
    keep_indices = np.where(keep)[0]

    n_written = 0
    with open(out_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "query_name", "match_name",
            "containment_query_in_match", "containment_match_in_query",
            "max_containment", "jaccard",
            "intersect_hashes", "size_query", "size_match",
        ])
        for k in keep_indices:
            i       = int(pair_rows[k])
            j       = int(pair_cols[k])
            shared  = int(pair_counts[k])
            si      = int(sizes[i])
            sj      = int(sizes[j])
            union   = si + sj - shared
            jaccard = shared / union if union > 0 else 1.0
            cont_q_in_m = shared / sj    # |i∩j| / |j|
            cont_m_in_q = shared / si    # |i∩j| / |i|
            max_cont    = shared / min(si, sj)
            writer.writerow([
                genome_ids[i], genome_ids[j],
                f"{cont_q_in_m:.8f}", f"{cont_m_in_q:.8f}",
                f"{max_cont:.8f}",   f"{jaccard:.8f}",
                shared, si, sj,
            ])
            n_written += 1

    log.info(
        "  Filter done in %.1f s.  Pairs written: %d  |  Below threshold: %d",
        time.time() - t3, n_written, n_pairs_any - n_written,
    )

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    elapsed_total = time.time() - t_start
    h  = int(elapsed_total // 3600)
    m_ = int((elapsed_total % 3600) // 60)
    s_ = int(elapsed_total % 60)

    out_size = out_path.stat().st_size / 1024 ** 2 if out_path.exists() else 0.0

    stats = {
        "script":                   "compute_fracminhash_candidates.py",
        "sketch_dir":               str(sketch_dir),
        "n_genomes":                n,
        "threshold":                args.threshold,
        "num_shards":               NUM_SHARDS,
        "n_workers":                n_workers,
        "n_pairs_sharing_any_hash": n_pairs_any,
        "n_pairs_above_threshold":  n_written,
        "wall_clock_seconds":       round(elapsed_total, 1),
        "output_csv":               str(out_path),
        "output_csv_mb":            round(out_size, 2),
    }
    stats_path = out_path.with_suffix(".stats.json")
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)

    print()
    print("============================== Run Summary ==============================")
    print(f"  Wall-clock time           : {h:02d}:{m_:02d}:{s_:02d} (hh:mm:ss)")
    print(f"  Genomes processed         : {n}")
    print(f"  Max-containment threshold : {args.threshold}")
    print(f"  Pairs sharing ≥1 hash     : {n_pairs_any}")
    print(f"  Pairs above threshold     : {n_written}")
    print(f"  Output CSV                : {out_path}  ({out_size:.1f} MB)")
    print(f"  Run stats JSON            : {stats_path}")
    print("=========================================================================")
    print()
    print("Next step:")
    print("  bash scripts/02_kmc_pairwise.sh")


if __name__ == "__main__":
    main()
