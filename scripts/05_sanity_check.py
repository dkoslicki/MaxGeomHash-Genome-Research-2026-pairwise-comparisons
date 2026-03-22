#!/usr/bin/env python3
"""
05_sanity_check.py
==================
Compare AlphaMaxGeomHash pairwise estimates against FracMinHash (sourmash)
estimates and — if available — exact KMC values.

For each shared pair, computes per-metric error statistics:
  - Mean Absolute Error (MAE)
  - Root Mean Square Error (RMSE)
  - Pearson correlation (r)
  - Spearman rank correlation (rho)
  - Median absolute error

Produces a JSON summary and, if matplotlib is available, scatter plots.

Performance design
------------------
The previous version used Python dicts with string-tuple keys, which caused
multi-hour stalls when the AMG dataset grew to 2M pairs (full GTDB run).
This version avoids all string-keyed Python dicts by:

  1. Encoding every genome pair as a single int64 key:
       key = min(genome_idx_A, genome_idx_B) * N_genomes + max(idx_A, idx_B)
     This reduces pair identity to a fast integer comparison.

  2. Loading FracMinHash CSV with pandas (C-level CSV parser) instead of
     csv.DictReader; the entire 2M-row file loads in ~3 seconds.

  3. Aligning datasets with np.intersect1d(return_indices=True) — no Python
     set operations or sorted() on string tuples.

  4. Remapping the KMC genome index to the AMG genome index numerically so
     KMC NPZ arrays can be used directly without any Python loop.

Runtime (full 2M-pair GTDB dataset): ~20 seconds.

Usage:
  python3 05_sanity_check.py \\
      --amg-pairwise  /scratch/.../alphamaxgeom_pairwise \\
      --fracminhash   /scratch/.../gtdb_pairwise_containment.csv \\
      --kmc-pairwise  /scratch/.../kmc_pairwise \\         # optional
      --output        /scratch/.../alphamaxgeom_pairwise/sanity_check
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

def encode_pair_keys(rows: np.ndarray, cols: np.ndarray, N: int) -> np.ndarray:
    """
    Encode (row, col) integer pairs as a single int64 key.
    Pairs are canonicalized so that the smaller index always comes first,
    making keys order-independent (i.e. pair (A,B) == pair (B,A)).

    key = min(row, col) * N + max(row, col)

    N^2 for GTDB: 143614^2 ~ 2.06e10, safely within int64 range.
    """
    rows = rows.astype(np.int64)
    cols = cols.astype(np.int64)
    lo = np.minimum(rows, cols)
    hi = np.maximum(rows, cols)
    return lo * N + hi


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_amg(pairwise_dir: Path):
    """
    Load AMG pairwise NPZ.

    Returns
    -------
    genome_ids : list[str]
        Maps integer genome index -> genome_id string.
    genome_id_to_int : dict[str, int]
        Reverse mapping used by subsequent loaders.
    keys : np.ndarray[int64]
        Sorted canonical pair keys (len = n_pairs).
    jac, mc, cqm, cmq : np.ndarray[float64]
        Jaccard, max_containment, containment_query_in_match,
        containment_match_in_query — aligned with keys.
    """
    npz_path   = pairwise_dir / "pairwise_results.npz"
    index_path = pairwise_dir / "genome_index.json"
    for p in (npz_path, index_path):
        if not p.exists():
            log.error("AMG file not found: %s", p)
            sys.exit(1)

    genome_ids     = json.load(open(index_path))["genomes"]
    genome_id_to_int = {gid: i for i, gid in enumerate(genome_ids)}
    N              = len(genome_ids)

    data = np.load(npz_path)
    keys = encode_pair_keys(data["row"], data["col"], N)
    order = np.argsort(keys, kind="stable")
    keys  = keys[order]

    log.info("AMG: %d pairs over %d genomes (from %s)", len(keys), N, npz_path)
    return (
        genome_ids,
        genome_id_to_int,
        N,
        keys,
        data["jaccard"][order].astype(np.float64),
        data["max_containment"][order].astype(np.float64),
        data["containment_query_in_match"][order].astype(np.float64),
        data["containment_match_in_query"][order].astype(np.float64),
    )


def load_fracminhash(csv_path: Path, genome_id_to_int: dict, N: int):
    """
    Load the FracMinHash pairwise CSV with pandas (fast C-level parser).

    Genome names are mapped to integers via genome_id_to_int so that pair
    keys are int64 rather than string tuples.  Both (A,B) and (B,A) rows
    appear in the CSV; after canonicalization, we keep the row with the
    higher max_containment for each canonical pair.

    Returns: (keys, jac, mc, cont) — all sorted by keys.
    """
    log.info("FracMinHash: reading %s ...", csv_path)
    df = pd.read_csv(
        csv_path,
        usecols=["query_name", "match_name", "jaccard",
                 "max_containment", "containment"],
        dtype={"jaccard": "float32", "max_containment": "float32",
               "containment": "float32"},
    )
    log.info("  %d rows loaded", len(df))

    # Map genome IDs to integers; rows with unknown IDs are dropped.
    q_idx = df["query_name"].map(genome_id_to_int)
    m_idx = df["match_name"].map(genome_id_to_int)
    valid = q_idx.notna() & m_idx.notna() & (q_idx != m_idx)
    df = df[valid].copy()
    q_idx = q_idx[valid].astype(np.int64).values
    m_idx = m_idx[valid].astype(np.int64).values

    lo = np.minimum(q_idx, m_idx)
    hi = np.maximum(q_idx, m_idx)
    raw_keys = lo * N + hi

    jac  = df["jaccard"].values.astype(np.float64)
    mc   = df["max_containment"].values.astype(np.float64)
    cont = df["containment"].values.astype(np.float64)

    # Sort by key, then deduplicate: for each canonical pair keep the row
    # with the highest max_containment (both (A,B) and (B,A) may appear).
    order  = np.argsort(raw_keys, kind="stable")
    raw_keys, jac, mc, cont = raw_keys[order], jac[order], mc[order], cont[order]

    # Find unique keys and, among duplicates, the position with max mc.
    # Strategy: stable-sort by (key ASC, mc DESC) — then np.unique keeps first.
    order2 = np.lexsort((-mc, raw_keys))   # primary key ASC, secondary -mc DESC
    raw_keys, jac, mc, cont = (raw_keys[order2], jac[order2],
                                mc[order2],      cont[order2])
    _, first = np.unique(raw_keys, return_index=True)
    keys = raw_keys[first]
    jac  = jac[first]
    mc   = mc[first]
    cont = cont[first]

    log.info("FracMinHash: %d canonical pairs retained", len(keys))
    return keys, jac, mc, cont


def load_kmc(pairwise_dir: Path, genome_id_to_int: dict, N: int):
    """
    Load KMC exact pairwise NPZ, remapping the KMC genome integer index to
    the AMG genome integer index so pair keys are comparable across datasets.

    Returns: (keys, jac, mc, cqm, cmq) — all sorted by keys.
    Returns empty arrays if KMC files are absent.
    """
    npz_path   = pairwise_dir / "pairwise_results.npz"
    index_path = pairwise_dir / "genome_index.json"

    empty = (np.empty(0, np.int64),) + (np.empty(0, np.float64),) * 4
    if not npz_path.exists() or not index_path.exists():
        log.warning("KMC pairwise files not found in %s — skipping.", pairwise_dir)
        return empty

    kmc_genome_ids = json.load(open(index_path))["genomes"]

    # Build a lookup array: kmc_int -> amg_int  (-1 if genome not in AMG index)
    kmc_to_amg = np.array(
        [genome_id_to_int.get(gid, -1) for gid in kmc_genome_ids],
        dtype=np.int64,
    )

    data      = np.load(npz_path)
    kmc_rows  = data["row"].astype(np.int64)
    kmc_cols  = data["col"].astype(np.int64)

    amg_rows = kmc_to_amg[kmc_rows]
    amg_cols = kmc_to_amg[kmc_cols]

    # Drop pairs where either genome is not in the AMG index
    valid    = (amg_rows >= 0) & (amg_cols >= 0)
    amg_rows = amg_rows[valid]
    amg_cols = amg_cols[valid]

    keys  = encode_pair_keys(amg_rows, amg_cols, N)
    order = np.argsort(keys, kind="stable")
    keys  = keys[order]

    jac = data["jaccard"][valid][order].astype(np.float64)
    mc  = data["max_containment"][valid][order].astype(np.float64)
    cqm = data["containment_query_in_match"][valid][order].astype(np.float64)
    cmq = data["containment_match_in_query"][valid][order].astype(np.float64)

    log.info("KMC: %d pairs (from %d total) from %s", len(keys), len(kmc_rows), npz_path)
    return keys, jac, mc, cqm, cmq


# ---------------------------------------------------------------------------
# Pair alignment
# ---------------------------------------------------------------------------

def align_pairs(keys_a: np.ndarray, keys_b: np.ndarray):
    """
    Return index arrays (ia, ib) such that keys_a[ia] == keys_b[ib].
    Both input arrays must be sorted.
    """
    shared, ia, ib = np.intersect1d(keys_a, keys_b, return_indices=True)
    return ia, ib


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def error_stats(true_vals: np.ndarray, pred_vals: np.ndarray, label: str) -> dict:
    """Compute standard accuracy metrics between reference and predicted arrays."""
    diff = pred_vals - true_vals
    mae  = float(np.mean(np.abs(diff)))
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    med  = float(np.median(np.abs(diff)))

    pearson_r = float(np.corrcoef(true_vals, pred_vals)[0, 1]) \
        if true_vals.std() > 0 and pred_vals.std() > 0 else float("nan")

    # Spearman rank correlation (no scipy dependency)
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
        "  %-40s  n=%7d  MAE=%.5f  RMSE=%.5f  r=%.4f  rho=%.4f",
        label, len(true_vals), mae, rmse,
        pearson_r    if not np.isnan(pearson_r)    else -99,
        spearman_rho if not np.isnan(spearman_rho) else -99,
    )
    return result


# ---------------------------------------------------------------------------
# Plotting (optional — skipped gracefully if matplotlib is absent)
# ---------------------------------------------------------------------------

def make_scatter_plot(x, y, xlabel, ylabel, title, outpath):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.info("  matplotlib not available — skipping: %s", outpath)
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
        description="Sanity check: compare AlphaMaxGeomHash vs FracMinHash and KMC.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--amg-pairwise", required=True,
                   help="Directory with AMG pairwise_results.npz + genome_index.json")
    p.add_argument("--fracminhash", required=True,
                   help="FracMinHash pairwise CSV (gtdb_pairwise_containment.csv)")
    p.add_argument("--kmc-pairwise", default="",
                   help="KMC pairwise directory (optional)")
    p.add_argument("--output", required=True,
                   help="Output directory for JSON summary and scatter plots")
    return p.parse_args()


def main():
    args   = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load datasets — all work with integer pair keys, not string tuples
    # ------------------------------------------------------------------
    (genome_ids, genome_id_to_int, N,
     amg_keys, amg_jac, amg_mc, amg_cqm, amg_cmq) = load_amg(
        Path(args.amg_pairwise))

    fmh_keys, fmh_jac, fmh_mc, fmh_cont = load_fracminhash(
        Path(args.fracminhash), genome_id_to_int, N)

    if args.kmc_pairwise:
        kmc_keys, kmc_jac, kmc_mc, kmc_cqm, kmc_cmq = load_kmc(
            Path(args.kmc_pairwise), genome_id_to_int, N)
    else:
        kmc_keys = np.empty(0, np.int64)
        kmc_jac = kmc_mc = kmc_cqm = kmc_cmq = np.empty(0, np.float64)

    # ------------------------------------------------------------------
    # Align pairs using integer key intersection (fast numpy operation)
    # ------------------------------------------------------------------
    af_ia, af_ib = align_pairs(amg_keys, fmh_keys)
    ak_ia, ak_ib = align_pairs(amg_keys, kmc_keys)
    fk_fa, fk_kb = align_pairs(fmh_keys, kmc_keys)

    log.info("Shared AMG ∩ FMH: %d pairs", len(af_ia))
    if len(kmc_keys):
        log.info("Shared AMG ∩ KMC: %d pairs", len(ak_ia))
        log.info("Shared FMH ∩ KMC: %d pairs", len(fk_fa))

    all_stats = []

    # ---- AMG vs FracMinHash -----------------------------------------------
    if len(af_ia):
        log.info("\n--- AlphaMaxGeomHash vs FracMinHash (sourmash) ---")
        a_jac = amg_jac[af_ia];  f_jac = fmh_jac[af_ib]
        a_mc  = amg_mc[af_ia];   f_mc  = fmh_mc[af_ib]

        all_stats.append(error_stats(f_jac, a_jac, "AMG vs FMH — Jaccard"))
        all_stats.append(error_stats(f_mc,  a_mc,  "AMG vs FMH — max_containment"))

        bad_jac = int(np.sum(np.abs(a_jac - f_jac) > 0.1))
        log.info(
            "  Pairs where |AMG_jaccard − FMH_jaccard| > 0.10:  %d / %d  (%.1f%%)",
            bad_jac, len(af_ia), 100 * bad_jac / max(1, len(af_ia)),
        )

        make_scatter_plot(f_jac, a_jac,
                          "FracMinHash Jaccard", "AlphaMaxGeomHash Jaccard",
                          "AMG vs FMH — Jaccard",
                          str(out_dir / "scatter_amg_vs_fmh_jaccard.png"))
        make_scatter_plot(f_mc, a_mc,
                          "FracMinHash max_containment", "AlphaMaxGeomHash max_containment",
                          "AMG vs FMH — max_containment",
                          str(out_dir / "scatter_amg_vs_fmh_maxcont.png"))

    # ---- AMG vs KMC ----------------------------------------------------------
    if len(ak_ia):
        log.info("\n--- AlphaMaxGeomHash vs KMC (exact) ---")
        a_jac = amg_jac[ak_ia];  k_jac = kmc_jac[ak_ib]
        a_mc  = amg_mc[ak_ia];   k_mc  = kmc_mc[ak_ib]

        all_stats.append(error_stats(k_jac, a_jac, "AMG vs KMC — Jaccard"))
        all_stats.append(error_stats(k_mc,  a_mc,  "AMG vs KMC — max_containment"))

        make_scatter_plot(k_jac, a_jac,
                          "KMC Jaccard (exact)", "AlphaMaxGeomHash Jaccard",
                          "AMG vs KMC (exact) — Jaccard",
                          str(out_dir / "scatter_amg_vs_kmc_jaccard.png"))
        make_scatter_plot(k_mc, a_mc,
                          "KMC max_containment (exact)", "AlphaMaxGeomHash max_containment",
                          "AMG vs KMC (exact) — max_containment",
                          str(out_dir / "scatter_amg_vs_kmc_maxcont.png"))

    # ---- FMH vs KMC (reference baseline) ------------------------------------
    if len(fk_fa):
        log.info("\n--- FracMinHash (sourmash) vs KMC (exact) ---")
        f_jac = fmh_jac[fk_fa];  k_jac = kmc_jac[fk_kb]
        f_mc  = fmh_mc[fk_fa];   k_mc  = kmc_mc[fk_kb]

        all_stats.append(error_stats(k_jac, f_jac, "FMH vs KMC — Jaccard"))
        all_stats.append(error_stats(k_mc,  f_mc,  "FMH vs KMC — max_containment"))

    # ---- Save summary --------------------------------------------------------
    summary = {
        "n_amg_pairs":         int(len(amg_keys)),
        "n_fracminhash_pairs": int(len(fmh_keys)),
        "n_kmc_pairs":         int(len(kmc_keys)),
        "n_amg_vs_fmh_shared": int(len(af_ia)),
        "n_amg_vs_kmc_shared": int(len(ak_ia)),
        "n_fmh_vs_kmc_shared": int(len(fk_fa)),
        "error_statistics":    all_stats,
    }
    summary_path = str(out_dir / "sanity_check_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print()
    print("============================== Sanity Check Summary ==============================")
    for s in all_stats:
        r_str = f"{s['pearson_r']:.4f}" if s["pearson_r"] is not None else "N/A"
        print(f"  {s['metric']:<45}  n={s['n_pairs']:>9,}  "
              f"MAE={s['mae']:.5f}  r={r_str}")
    print()
    print(f"  Full stats written to: {summary_path}")

    fmh_jac_stat = next(
        (s for s in all_stats if "AMG vs FMH" in s["metric"] and "Jaccard" in s["metric"]),
        None,
    )
    if fmh_jac_stat:
        r   = fmh_jac_stat.get("pearson_r") or 0
        mae = fmh_jac_stat["mae"]
        print()
        if r >= 0.90 and mae <= 0.05:
            print("  SANITY CHECK PASSED  (Pearson r >= 0.90 and MAE <= 0.05 for Jaccard)")
        else:
            print(f"  WARNING: AMG–FMH Jaccard correlation r={r:.3f}, MAE={mae:.4f}")
            print("  Results may be outside expected range — investigate before scaling up.")
    print("==================================================================================")


if __name__ == "__main__":
    main()
