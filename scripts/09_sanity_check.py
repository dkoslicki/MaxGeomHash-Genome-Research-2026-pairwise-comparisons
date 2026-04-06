#!/usr/bin/env python3
"""
09_sanity_check.py
==================
Compare sketching-method pairwise estimates against each other and against
exact KMC values.  Supports all three kmer-sketch algorithms studied here:

    MaxGeomHash       (--mg-pairwise)
    AlphaMaxGeomHash  (--amg-pairwise)
    BottomK           (--bottomk-pairwise)
    FracMinHash       (--fracminhash-pairwise)    <- kmer-sketch binary

and, optionally, the Sourmash FracMinHash baseline:

    Sourmash FMH      (--sourmash-csv)            <- kept for backward-compat
                                                     but excluded from default
                                                     output; use --include-sourmash
                                                     to add it back

For each shared pair, computes per-metric error statistics:
  - Mean Absolute Error (MAE)
  - Root Mean Square Error (RMSE)
  - Pearson correlation (r)
  - Spearman rank correlation (rho)
  - Median absolute error

Produces a JSON summary and, if matplotlib is available, scatter plots.

Performance design
------------------
All string genome IDs are mapped to int64 indices immediately after loading.
Every pair is encoded as a single int64 key:

    key = min(genome_idx_A, genome_idx_B) * N_genomes + max(idx_A, idx_B)

Pair intersection is done with np.intersect1d (not Python set operations).
This allows the full 2M-pair GTDB dataset to be processed in ~20 seconds.

Usage:
  # Compare all three kmer-sketch methods against KMC (most common case):
  python3 09_sanity_check.py \\
      --mg-pairwise           /scratch/.../maxgeom_pairwise \\    
      --amg-pairwise          /scratch/.../alphamaxgeom_pairwise \\
      --bottomk-pairwise      /scratch/.../bottomk_pairwise \\
      --fracminhash-pairwise  /scratch/.../fracminhash_pairwise \\
      --kmc-pairwise          /scratch/.../kmc_pairwise \\
      --output                /scratch/.../sanity_check

  # Include Sourmash FracMinHash in output as well:
  python3 09_sanity_check.py ... --include-sourmash \\
      --sourmash-csv /scratch/.../gtdb_pairwise_containment.csv

  # Minimal: only AMG vs KMC (backwards-compatible with old pipeline):
  python3 09_sanity_check.py \\
      --mg-pairwise   /scratch/.../maxgeom_pairwise \\  
      --amg-pairwise  /scratch/.../alphamaxgeom_pairwise \\
      --kmc-pairwise  /scratch/.../kmc_pairwise \\
      --output        /scratch/.../sanity_check
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def encode_pair_keys(rows, cols, N):
    """
    Encode (row, col) integer pairs as a single int64 key.
    key = min(row, col) * N + max(row, col)
    """
    rows = rows.astype(np.int64)
    cols = cols.astype(np.int64)
    lo = np.minimum(rows, cols)
    hi = np.maximum(rows, cols)
    return lo * N + hi


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_npz_pairwise(pairwise_dir, label, genome_id_to_int=None, N=None):
    """
    Load a pairwise NPZ produced by any of the *_pairwise.py scripts.

    If genome_id_to_int is None (first call), builds the index from the
    accompanying genome_index.json and returns it alongside the data.

    Returns (when genome_id_to_int is None):
        genome_ids, genome_id_to_int, N, keys, jac, mc, cqm, cmq

    Returns (when genome_id_to_int is already provided):
        keys, jac, mc, cqm, cmq
    """
    pairwise_dir = Path(pairwise_dir)
    npz_path   = pairwise_dir / "pairwise_results.npz"
    index_path = pairwise_dir / "genome_index.json"

    for p in (npz_path, index_path):
        if not p.exists():
            log.error("%s file not found: %s", label, p)
            sys.exit(1)

    with open(index_path) as f:
        index_data = json.load(f)
    local_genome_ids = index_data["genomes"]

    build_index = genome_id_to_int is None
    if build_index:
        genome_ids = local_genome_ids
        genome_id_to_int = {gid: i for i, gid in enumerate(genome_ids)}
        N = len(genome_ids)

    # Remap local integer indices to the shared integer space
    local_to_shared = np.array(
        [genome_id_to_int.get(gid, -1) for gid in local_genome_ids],
        dtype=np.int64,
    )

    data = np.load(npz_path)
    raw_rows = local_to_shared[data["row"].astype(np.int64)]
    raw_cols = local_to_shared[data["col"].astype(np.int64)]

    valid = (raw_rows >= 0) & (raw_cols >= 0)
    raw_rows = raw_rows[valid]
    raw_cols = raw_cols[valid]

    keys  = encode_pair_keys(raw_rows, raw_cols, N)
    order = np.argsort(keys, kind="stable")
    keys  = keys[order]

    jac = data["jaccard"][valid][order].astype(np.float64)
    mc  = data["max_containment"][valid][order].astype(np.float64)
    cqm = data["containment_query_in_match"][valid][order].astype(np.float64)
    cmq = data["containment_match_in_query"][valid][order].astype(np.float64)

    log.info("%s: %d pairs over %d genomes (from %s)", label, len(keys), N, npz_path)

    if build_index:
        return genome_ids, genome_id_to_int, N, keys, jac, mc, cqm, cmq
    return keys, jac, mc, cqm, cmq


def load_sourmash_csv(csv_path, genome_id_to_int, N):
    """
    Load the Sourmash FracMinHash pairwise CSV with pandas.
    Kept for backward-compatibility; excluded from default output.

    Returns: (keys, jac, mc, cont) -- all sorted by keys.
    """
    csv_path = Path(csv_path)
    log.info("Sourmash FMH: reading %s ...", csv_path)
    df = pd.read_csv(
        csv_path,
        usecols=["query_name", "match_name", "jaccard",
                 "max_containment", "containment"],
        dtype={"jaccard": "float32", "max_containment": "float32",
               "containment": "float32"},
    )
    log.info("  %d rows loaded", len(df))

    q_idx = df["query_name"].map(genome_id_to_int)
    m_idx = df["match_name"].map(genome_id_to_int)
    valid = q_idx.notna() & m_idx.notna() & (q_idx != m_idx)
    df    = df[valid].copy()
    q_idx = q_idx[valid].astype(np.int64).values
    m_idx = m_idx[valid].astype(np.int64).values

    lo = np.minimum(q_idx, m_idx)
    hi = np.maximum(q_idx, m_idx)
    raw_keys = lo * N + hi

    jac  = df["jaccard"].values.astype(np.float64)
    mc   = df["max_containment"].values.astype(np.float64)
    cont = df["containment"].values.astype(np.float64)

    # Deduplicate: keep the row with the highest max_containment per canonical pair
    order2 = np.lexsort((-mc, raw_keys))
    raw_keys, jac, mc, cont = (raw_keys[order2], jac[order2],
                                mc[order2],      cont[order2])
    _, first = np.unique(raw_keys, return_index=True)
    keys = raw_keys[first]
    jac  = jac[first]
    mc   = mc[first]
    cont = cont[first]

    log.info("Sourmash FMH: %d canonical pairs retained", len(keys))
    return keys, jac, mc, cont


# ---------------------------------------------------------------------------
# Pair alignment
# ---------------------------------------------------------------------------

def align_pairs(keys_a, keys_b):
    """Return index arrays (ia, ib) such that keys_a[ia] == keys_b[ib]."""
    _, ia, ib = np.intersect1d(keys_a, keys_b, return_indices=True)
    return ia, ib


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def error_stats(true_vals, pred_vals, label):
    """Compute standard accuracy metrics between reference and predicted arrays."""
    diff = pred_vals - true_vals
    mae  = float(np.mean(np.abs(diff)))
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    med  = float(np.median(np.abs(diff)))

    pearson_r = float(np.corrcoef(true_vals, pred_vals)[0, 1]) \
        if true_vals.std() > 0 and pred_vals.std() > 0 else float("nan")

    rk_t = np.argsort(np.argsort(true_vals)).astype(float)
    rk_p = np.argsort(np.argsort(pred_vals)).astype(float)
    spearman_rho = float(np.corrcoef(rk_t, rk_p)[0, 1]) \
        if rk_t.std() > 0 and rk_p.std() > 0 else float("nan")

    result = {
        "metric":       label,
        "n_pairs":      len(true_vals),
        "mae":          round(mae,  6),
        "rmse":         round(rmse, 6),
        "median_ae":    round(med,  6),
        "pearson_r":    round(pearson_r,    6) if not np.isnan(pearson_r)    else None,
        "spearman_rho": round(spearman_rho, 6) if not np.isnan(spearman_rho) else None,
    }
    log.info(
        "  %-55s  n=%7d  MAE=%.5f  RMSE=%.5f  r=%.4f  rho=%.4f",
        label, len(true_vals), mae, rmse,
        pearson_r    if not np.isnan(pearson_r)    else -99,
        spearman_rho if not np.isnan(spearman_rho) else -99,
    )
    return result


# ---------------------------------------------------------------------------
# Plotting (optional)
# ---------------------------------------------------------------------------

def make_scatter_plot(x, y, xlabel, ylabel, title, outpath):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.info("  matplotlib not available -- skipping: %s", outpath)
        return

    fig, ax = plt.subplots(figsize=(5, 5))
    n = len(x)
    if n > 10_000:
        rng = np.random.default_rng(42)
        idx = rng.choice(n, 10_000, replace=False)
        x, y = x[idx], y[idx]

    ax.scatter(x, y, s=4, alpha=0.4, linewidths=0)
    ax.plot([0, 1], [0, 1], "r--", linewidth=1, label="y = x")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=120)
    plt.close(fig)
    log.info("  Scatter plot saved: %s", outpath)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Sanity check: compare kmer-sketch methods against KMC exact values.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # --- method pairwise directories ------------------------------------------
    p.add_argument("--mg-pairwise", default="",
                   help="MaxGeomHash pairwise directory "
                        "(pairwise_results.npz + genome_index.json)")
    p.add_argument("--amg-pairwise", default="",
                   help="AlphaMaxGeomHash pairwise directory "
                        "(pairwise_results.npz + genome_index.json)")
    p.add_argument("--bottomk-pairwise", default="",
                   help="BottomK pairwise directory "
                        "(pairwise_results.npz + genome_index.json)")
    p.add_argument("--fracminhash-pairwise", default="",
                   help="FracMinHash (kmer-sketch) pairwise directory "
                        "(pairwise_results.npz + genome_index.json)")
    # --- ground truth ----------------------------------------------------------
    p.add_argument("--kmc-pairwise", default="",
                   help="KMC exact pairwise directory (optional)")
    # --- Sourmash baseline (kept but not shown by default) --------------------
    p.add_argument("--sourmash-csv", default="",
                   help="Sourmash FracMinHash pairwise CSV "
                        "(gtdb_pairwise_containment.csv). "
                        "Loaded only if --include-sourmash is set.")
    p.add_argument("--include-sourmash", action="store_true",
                   help="Include Sourmash FracMinHash in comparisons and output "
                        "(requires --sourmash-csv). Off by default.")
    # --- output ---------------------------------------------------------------
    p.add_argument("--output", required=True,
                   help="Output directory for JSON summary and scatter plots")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Determine which NPZ datasets are available
    # ------------------------------------------------------------------
    method_dirs = {}
    if args.amg_pairwise and Path(args.amg_pairwise).exists():
        method_dirs["AMG"] = Path(args.amg_pairwise)
    if args.bottomk_pairwise and Path(args.bottomk_pairwise).exists():
        method_dirs["BK"] = Path(args.bottomk_pairwise)
    if args.fracminhash_pairwise and Path(args.fracminhash_pairwise).exists():
        method_dirs["FMH_ks"] = Path(args.fracminhash_pairwise)
    if args.mg_pairwise and Path(args.mg_pairwise).exists():
        method_dirs["MG"] = Path(args.mg_pairwise)

    if not method_dirs:
        log.error("No method pairwise directories provided / found. "
                  "Provide at least one of --amg-pairwise, --bottomk-pairwise, "
                  "--fracminhash-pairwise, --mg-pairwise.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Load datasets into a shared integer index
    # ------------------------------------------------------------------
    genome_ids = None
    genome_id_to_int = None
    N = None
    method_data = {}   # method_key -> (keys, jac, mc, cqm, cmq)

    for method_key, pdir in method_dirs.items():
        if genome_id_to_int is None:
            genome_ids, genome_id_to_int, N, *rest = load_npz_pairwise(
                pdir, method_key)
            method_data[method_key] = tuple(rest)
        else:
            method_data[method_key] = load_npz_pairwise(
                pdir, method_key, genome_id_to_int, N)

    # ------------------------------------------------------------------
    # Load KMC ground truth
    # ------------------------------------------------------------------
    kmc_keys = kmc_jac = kmc_mc = kmc_cqm = kmc_cmq = None
    if args.kmc_pairwise:
        kmc_dir = Path(args.kmc_pairwise)
        if kmc_dir.exists():
            log.info("Loading KMC ground truth from %s ...", kmc_dir)
            kmc_keys, kmc_jac, kmc_mc, kmc_cqm, kmc_cmq = load_npz_pairwise(
                kmc_dir, "KMC", genome_id_to_int, N)
        else:
            log.warning("KMC pairwise directory not found: %s -- skipping.", kmc_dir)

    # ------------------------------------------------------------------
    # Load Sourmash FracMinHash baseline (optional, off by default)
    # ------------------------------------------------------------------
    smash_keys = smash_jac = smash_mc = None
    if args.include_sourmash and args.sourmash_csv:
        smash_csv = Path(args.sourmash_csv)
        if smash_csv.exists():
            smash_keys, smash_jac, smash_mc, _ = load_sourmash_csv(
                smash_csv, genome_id_to_int, N)
        else:
            log.warning("Sourmash CSV not found: %s", smash_csv)

    # ------------------------------------------------------------------
    # Human-readable method labels
    # ------------------------------------------------------------------
    method_label = {
        "AMG":    "AlphaMaxGeomHash",
        "BK":     "BottomK",
        "FMH_ks": "FracMinHash (kmer-sketch)",
        "MG":     "MaxGeomHash",
    }

    method_keys_list = list(method_data.keys())
    all_stats = []
    n_shared  = {}

    # ---- Each method vs KMC (exact) ---------------------------------------
    if kmc_keys is not None:
        for mk in method_keys_list:
            mkeys, mjac, mmc, _, _ = method_data[mk]
            ia, ib = align_pairs(mkeys, kmc_keys)
            n_shared[(mk, "KMC")] = len(ia)
            if len(ia):
                label = method_label[mk]
                log.info("\n--- %s vs KMC (exact) ---", label)
                all_stats.append(
                    error_stats(kmc_jac[ib], mjac[ia], f"{label} vs KMC -- Jaccard"))
                all_stats.append(
                    error_stats(kmc_mc[ib],  mmc[ia],  f"{label} vs KMC -- max_containment"))
                make_scatter_plot(
                    kmc_jac[ib], mjac[ia],
                    "KMC Jaccard (exact)", f"{label} Jaccard",
                    f"{label} vs KMC (exact) -- Jaccard",
                    str(out_dir / f"scatter_{mk.lower()}_vs_kmc_jaccard.png"))
                make_scatter_plot(
                    kmc_mc[ib], mmc[ia],
                    "KMC max_containment (exact)", f"{label} max_containment",
                    f"{label} vs KMC (exact) -- max_containment",
                    str(out_dir / f"scatter_{mk.lower()}_vs_kmc_maxcont.png"))

    # ---- Each method vs every other method --------------------------------
    for i, mk_a in enumerate(method_keys_list):
        for mk_b in method_keys_list[i+1:]:
            keys_a, jac_a, mc_a, _, _ = method_data[mk_a]
            keys_b, jac_b, mc_b, _, _ = method_data[mk_b]
            ia, ib = align_pairs(keys_a, keys_b)
            n_shared[(mk_a, mk_b)] = len(ia)
            if len(ia):
                label_a = method_label[mk_a]
                label_b = method_label[mk_b]
                log.info("\n--- %s vs %s ---", label_a, label_b)
                all_stats.append(
                    error_stats(jac_b[ib], jac_a[ia],
                                f"{label_a} vs {label_b} -- Jaccard"))
                all_stats.append(
                    error_stats(mc_b[ib],  mc_a[ia],
                                f"{label_a} vs {label_b} -- max_containment"))

    # ---- Sourmash FMH (optional) ------------------------------------------
    if smash_keys is not None:
        for mk in method_keys_list:
            mkeys, mjac, mmc, _, _ = method_data[mk]
            ia, ib = align_pairs(mkeys, smash_keys)
            n_shared[(mk, "Sourmash")] = len(ia)
            if len(ia):
                label = method_label[mk]
                log.info("\n--- %s vs Sourmash FracMinHash ---", label)
                all_stats.append(
                    error_stats(smash_jac[ib], mjac[ia],
                                f"{label} vs Sourmash FMH -- Jaccard"))
                all_stats.append(
                    error_stats(smash_mc[ib],  mmc[ia],
                                f"{label} vs Sourmash FMH -- max_containment"))

        if kmc_keys is not None:
            ia, ib = align_pairs(smash_keys, kmc_keys)
            n_shared[("Sourmash", "KMC")] = len(ia)
            if len(ia):
                log.info("\n--- Sourmash FracMinHash vs KMC (exact) ---")
                all_stats.append(
                    error_stats(kmc_jac[ib], smash_jac[ia],
                                "Sourmash FMH vs KMC -- Jaccard"))
                all_stats.append(
                    error_stats(kmc_mc[ib],  smash_mc[ia],
                                "Sourmash FMH vs KMC -- max_containment"))

    # ------------------------------------------------------------------
    # Save summary
    # ------------------------------------------------------------------
    summary = {
        "n_genomes":          N,
        "methods_loaded":     list(method_dirs.keys()),
        "n_pairs_per_method": {mk: int(len(method_data[mk][0]))
                               for mk in method_data},
        "n_kmc_pairs":        int(len(kmc_keys)) if kmc_keys is not None else 0,
        "n_sourmash_pairs":   int(len(smash_keys)) if smash_keys is not None else 0,
        "n_shared_pairs":     {f"{a}_vs_{b}": int(n)
                               for (a, b), n in n_shared.items()},
        "error_statistics":   all_stats,
    }
    summary_path = str(out_dir / "sanity_check_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # ------------------------------------------------------------------
    # Print summary table
    # ------------------------------------------------------------------
    print()
    print("============================== Sanity Check Summary ==============================")
    for s in all_stats:
        r_str = f"{s['pearson_r']:.4f}" if s["pearson_r"] is not None else "N/A"
        print(f"  {s['metric']:<58}  n={s['n_pairs']:>9,}  "
              f"MAE={s['mae']:.5f}  r={r_str}")
    print()
    print(f"  Full stats written to: {summary_path}")

    # Pass/fail: each method must achieve r >= 0.90 and MAE <= 0.05 vs KMC
    for mk in method_keys_list:
        label = method_label[mk]
        stat  = next(
            (s for s in all_stats
             if f"{label} vs KMC" in s["metric"] and "Jaccard" in s["metric"]),
            None,
        )
        if stat:
            r   = stat.get("pearson_r") or 0
            mae = stat["mae"]
            print()
            if r >= 0.90 and mae <= 0.05:
                print(f"  PASSED  {label}: r={r:.4f}  MAE={mae:.5f}  "
                      "(r >= 0.90 and MAE <= 0.05 for Jaccard vs KMC)")
            else:
                print(f"  WARNING {label}: r={r:.4f}  MAE={mae:.5f}  "
                      "(outside expected range -- investigate before scaling up)")
    print("==================================================================================")


if __name__ == "__main__":
    main()
