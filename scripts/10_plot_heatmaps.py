#!/usr/bin/env python3
"""
10_plot_heatmaps.py
===================
Generate publication-quality pairwise similarity error heatmaps comparing
sketching methods against the KMC exact k-mer ground truth.

Supported methods (all optional; include whichever have been run):
  AlphaMaxGeomHash    --amg-pairwise          (NPZ directory)
  BottomK             --bottomk-pairwise       (NPZ directory)
  FracMinHash         --fracminhash-pairwise   (NPZ directory; kmer-sketch binary)
  Sourmash FracMinHash --sourmash-csv          (CSV file; off by default, use
                                                --include-sourmash to enable)

The heatmap H for a given sketching method S and similarity measure M is
defined per-pair (i,j) as:

    H[i,j] = |S_estimate(i,j) - KMC_exact(i,j)| / KMC_exact(i,j)

i.e., the relative error of the sketch estimate with respect to the ground
truth.  Color encoding: blue (0 = perfect estimate, cold/good) → yellow →
red (high relative error, hot/bad).  Missing pairs (below the 0.01 similarity
threshold) are shown in light gray.

To keep the heatmaps visually comparable:
  - All panels share the same colormap and vmin/vmax.
  - Genomes are clustered by hierarchical clustering (Ward's linkage) on KMC
    Jaccard distances so that similar genomes are adjacent, revealing block
    structure.  The SAME ordering is applied to ALL method panels.

Designing for extensibility
---------------------------
Adding a new method requires only passing its NPZ directory via the
appropriate flag.  The methods list in main() is built dynamically from
whichever flags are provided — no method-specific logic in the plotting loop.

Usage:
    # All three kmer-sketch methods (most common):
    python3 10_plot_heatmaps.py \\
        --amg-pairwise         data/GTDB/alphamaxgeom_pairwise \\
        --bottomk-pairwise     data/GTDB/bottomk_pairwise \\
        --fracminhash-pairwise data/GTDB/fracminhash_pairwise \\
        --kmc-pairwise         data/GTDB/kmc_pairwise \\
        --output               data/GTDB/figures \\
        [--metric jaccard|max_containment] \\
        [--n-genomes 500] \\
        [--show-reference]

    # Also include Sourmash FracMinHash panel:
    python3 10_plot_heatmaps.py ... \\
        --include-sourmash \\
        --sourmash-csv data/GTDB/gtdb_pairwise_containment.csv
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")   # non-interactive backend — safe for headless servers
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Matplotlib publication style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         9,
    "axes.titlesize":    10,
    "axes.labelsize":    9,
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "figure.dpi":        150,   # screen preview; PDF is vector anyway
    "savefig.dpi":       300,
    "savefig.bbox":      "tight",
    "pdf.fonttype":      42,    # embeds TrueType fonts for journal submission
    "ps.fonttype":       42,
})

# ---------------------------------------------------------------------------
# Integer pair-key encoding
# ---------------------------------------------------------------------------

def encode_pair_keys(rows: np.ndarray, cols: np.ndarray, N: int) -> np.ndarray:
    """
    Encode canonical (min, max) genome-index pairs as a single int64 value:
        key = min(row, col) * N + max(row, col)
    Safe for GTDB: N=143614 → N² ≈ 2.06e10, well within int64 range.
    """
    r = rows.astype(np.int64)
    c = cols.astype(np.int64)
    lo = np.minimum(r, c)
    hi = np.maximum(r, c)
    return lo * N + hi


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_npz_sorted(npz_path: Path, index_path: Path) -> dict:
    """
    Load a pairwise NPZ produced by 04_alphamaxgeom_pairwise.py or
    02_kmc_pairwise.py.  Returns a dict with:
        genome_ids       : list[str]  (index → genome_id)
        genome_id_to_int : dict[str, int]
        N                : int
        keys             : int64 array (sorted canonical pair keys)
        jaccard / max_containment / containment_query_in_match /
        containment_match_in_query : float64 arrays aligned with keys
    """
    if not npz_path.exists() or not index_path.exists():
        log.error("File not found: %s or %s", npz_path, index_path)
        sys.exit(1)

    genome_ids       = json.load(open(index_path))["genomes"]
    genome_id_to_int = {gid: i for i, gid in enumerate(genome_ids)}
    N                = len(genome_ids)

    data  = np.load(npz_path)
    keys  = encode_pair_keys(data["row"], data["col"], N)
    order = np.argsort(keys, kind="stable")
    keys  = keys[order]

    # Both AMG and KMC store these four arrays; pick whichever exist.
    def _get(name):
        return data[name][order].astype(np.float64) if name in data else None

    log.info("  Loaded %d pairs from %s", len(keys), npz_path.name)
    return {
        "genome_ids":       genome_ids,
        "genome_id_to_int": genome_id_to_int,
        "N":                N,
        "keys":             keys,
        "jaccard":                    _get("jaccard"),
        "max_containment":            _get("max_containment"),
        "containment_query_in_match": _get("containment_query_in_match"),
        "containment_match_in_query": _get("containment_match_in_query"),
    }


def load_npz_remapped(npz_path: Path, index_path: Path,
                      ref_genome_id_to_int: dict, N: int) -> dict:
    """
    Load a pairwise NPZ (same format as load_npz_sorted) and remap its genome
    integer indices to the reference genome index (from KMC or first-loaded
    dataset).  This is safe even when the NPZ was written with a different
    sorted order than the reference.

    Returns the same dict structure as load_npz_sorted (minus genome_ids /
    genome_id_to_int / N, which come from the reference).
    """
    if not npz_path.exists() or not index_path.exists():
        log.error("File not found: %s or %s", npz_path, index_path)
        sys.exit(1)

    local_genome_ids = json.load(open(index_path))["genomes"]
    # Map each local integer index -> reference integer index (-1 = unknown)
    local_to_ref = np.array(
        [ref_genome_id_to_int.get(gid, -1) for gid in local_genome_ids],
        dtype=np.int64,
    )

    data = np.load(npz_path)
    raw_rows = local_to_ref[data["row"].astype(np.int64)]
    raw_cols = local_to_ref[data["col"].astype(np.int64)]

    valid = (raw_rows >= 0) & (raw_cols >= 0)
    raw_rows = raw_rows[valid]
    raw_cols = raw_cols[valid]

    keys  = encode_pair_keys(raw_rows, raw_cols, N)
    order = np.argsort(keys, kind="stable")
    keys  = keys[order]

    def _get(name):
        return data[name][valid][order].astype(np.float64) if name in data else None

    log.info("  Loaded %d pairs (remapped) from %s", len(keys), npz_path.name)
    return {
        "keys":                       keys,
        "jaccard":                    _get("jaccard"),
        "max_containment":            _get("max_containment"),
        "containment_query_in_match": _get("containment_query_in_match"),
        "containment_match_in_query": _get("containment_match_in_query"),
    }


def load_fracminhash_sorted(csv_path: Path, genome_id_to_int: dict, N: int) -> dict:
    """
    Load the FracMinHash pairwise CSV with pandas.  Maps genome IDs to the
    shared integer index space so pair keys are directly comparable with
    AMG and KMC keys.

    When both (A→B) and (B→A) rows appear in the CSV, the row with the
    higher max_containment is kept for each canonical pair.

    Returns the same dict structure as load_npz_sorted (minus genome_ids /
    genome_id_to_int / N, which are inherited from the AMG/KMC index).
    """
    log.info("  Reading FracMinHash CSV: %s ...", csv_path.name)
    df = pd.read_csv(
        csv_path,
        usecols=["query_name", "match_name", "jaccard",
                 "max_containment", "containment"],
        dtype={"jaccard": "float32", "max_containment": "float32",
               "containment": "float32"},
    )
    log.info("    %d rows loaded", len(df))

    q_idx = df["query_name"].map(genome_id_to_int)
    m_idx = df["match_name"].map(genome_id_to_int)
    valid = q_idx.notna() & m_idx.notna() & (q_idx != m_idx)
    df    = df[valid].copy()
    q_idx = q_idx[valid].astype(np.int64).values
    m_idx = m_idx[valid].astype(np.int64).values

    lo      = np.minimum(q_idx, m_idx)
    hi      = np.maximum(q_idx, m_idx)
    rawkeys = lo * N + hi

    jac  = df["jaccard"].values.astype(np.float64)
    mc   = df["max_containment"].values.astype(np.float64)
    cont = df["containment"].values.astype(np.float64)

    # Sort by (key ASC, max_containment DESC), then unique keeps best row.
    order2  = np.lexsort((-mc, rawkeys))
    rawkeys, jac, mc, cont = rawkeys[order2], jac[order2], mc[order2], cont[order2]
    _, first = np.unique(rawkeys, return_index=True)

    keys = rawkeys[first]
    log.info("    %d canonical pairs retained", len(keys))
    return {
        "keys":                       keys,
        "jaccard":                    jac[first],
        "max_containment":            mc[first],
        "containment_query_in_match": cont[first],   # only one direction stored
        "containment_match_in_query": None,
    }


# ---------------------------------------------------------------------------
# Genome subset selection and submatrix construction
# ---------------------------------------------------------------------------

def select_high_degree_genomes(kmc_keys: np.ndarray, N: int,
                                n_genomes: int) -> np.ndarray:
    """
    Select the n_genomes genomes with the most edges (pairs) in the KMC
    similarity graph.  High-degree genomes produce the densest submatrix,
    minimising the fraction of missing (below-threshold) cells.

    Returns a sorted 1-D array of genome integer indices.
    """
    rows = (kmc_keys // N).astype(np.int64)
    cols = (kmc_keys  % N).astype(np.int64)
    degree = np.zeros(N, dtype=np.int32)
    np.add.at(degree, rows, 1)
    np.add.at(degree, cols, 1)

    top = np.argsort(degree)[-n_genomes:]
    return np.sort(top)   # sorted for reproducibility


def build_submatrix(sorted_keys: np.ndarray, metric_array: np.ndarray,
                    selected_idx: np.ndarray, N: int) -> np.ndarray:
    """
    Build a dense (M × M) similarity matrix for the selected genomes.
    Pairs absent from sorted_keys (below the similarity threshold) are
    filled with NaN.  The diagonal is set to 1.0.

    Uses fully vectorised searchsorted — no Python loop over pairs.
    """
    M   = len(selected_idx)
    mat = np.full((M, M), np.nan)
    np.fill_diagonal(mat, 1.0)

    if metric_array is None:
        return mat

    ii, jj = np.tril_indices(M, k=-1)   # lower-triangle pairs
    g1 = selected_idx[ii].astype(np.int64)
    g2 = selected_idx[jj].astype(np.int64)
    lo = np.minimum(g1, g2)
    hi = np.maximum(g1, g2)
    pair_keys = lo * N + hi

    pos       = np.searchsorted(sorted_keys, pair_keys)
    in_range  = pos < len(sorted_keys)
    # Clip before indexing to avoid out-of-bounds when sorted_keys is small
    # (e.g. test runs with only 50 pairs).  The in_range mask filters these out.
    pos_safe  = np.where(in_range, pos, 0)
    found     = in_range & (sorted_keys[pos_safe] == pair_keys)

    mat[ii[found], jj[found]] = metric_array[pos[found]]
    mat[jj[found], ii[found]] = metric_array[pos[found]]
    return mat


# ---------------------------------------------------------------------------
# Hierarchical clustering
# ---------------------------------------------------------------------------

def cluster_by_similarity(sim_matrix: np.ndarray) -> np.ndarray:
    """
    Cluster genomes by KMC Jaccard similarity using Ward's hierarchical
    linkage.  Missing pairs (NaN) are treated as maximally dissimilar
    (distance = 1.0).  Returns the leaf ordering that groups similar
    genomes together.

    The same ordering must be applied to ALL method panels so that they
    remain directly comparable.
    """
    dist = 1.0 - np.where(np.isnan(sim_matrix), 0.0, sim_matrix)
    np.fill_diagonal(dist, 0.0)

    # squareform converts the symmetric matrix to condensed upper-triangle
    # form required by scipy's linkage.
    condensed = squareform(dist, checks=False)
    Z         = linkage(condensed, method="ward", optimal_ordering=True)
    return leaves_list(Z)


# ---------------------------------------------------------------------------
# Error computation
# ---------------------------------------------------------------------------

def relative_error_matrix(estimate: np.ndarray, ground_truth: np.ndarray,
                           min_gt: float = 1e-6) -> np.ndarray:
    """
    Compute element-wise relative error:
        err[i,j] = |estimate - ground_truth| / ground_truth

    Cells where ground_truth is NaN (pair not in KMC dataset) or below
    min_gt (to avoid division by near-zero) are set to NaN and will
    appear as missing in the heatmap.
    """
    err = np.full_like(estimate, np.nan)
    valid = (~np.isnan(ground_truth)) & (~np.isnan(estimate)) & (ground_truth > min_gt)
    err[valid] = np.abs(estimate[valid] - ground_truth[valid]) / ground_truth[valid]
    return err


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def make_error_colormap():
    """
    Colormap: blue (0 = perfect estimate, cold/good) → yellow → red (high error, hot/bad).
    Uses RdYlBu_r (reversed RdYlBu): low values map to blue, high values to red.
    Perceptually uniform, print-friendly, and accessible to red-green colour-blind readers.
    """
    return plt.cm.RdYlBu_r


def plot_heatmaps(kmc_sim: np.ndarray,
                  methods: list,        # list of {"label": str, "err": ndarray}
                  leaf_order: np.ndarray,
                  metric_label: str,
                  out_dir: Path,
                  show_reference: bool = False):
    """
    Draw the error heatmap figure.

    Parameters
    ----------
    kmc_sim     : M×M KMC ground-truth similarity (reordered by leaf_order)
    methods     : list of dicts {"label": str, "err": M×M relative-error array}
    leaf_order  : clustering index for rows/columns
    metric_label: e.g. "Jaccard"
    out_dir     : directory to write PDF/PNG outputs
    show_reference : if True, add a leftmost panel showing KMC similarity
    """
    n_methods = len(methods)
    n_panels  = n_methods + (1 if show_reference else 0)
    cmap_err  = make_error_colormap()
    cmap_ref  = plt.cm.viridis

    # Determine shared colorscale: 99th percentile of all non-NaN errors
    all_err = np.concatenate([m["err"][~np.isnan(m["err"])].ravel()
                               for m in methods])
    vmax_err = float(np.percentile(all_err, 99)) if len(all_err) else 1.0
    vmax_err = min(vmax_err, 1.0)     # cap at 100 % relative error
    vmax_ref = float(np.nanmax(kmc_sim)) if not np.all(np.isnan(kmc_sim)) else 1.0

    missing_color = "#CCCCCC"   # light gray for pairs below similarity threshold

    panel_w  = 3.8   # inches per panel
    cbar_h   = 0.35  # colorbar height in inches
    fig_w    = panel_w * n_panels + 0.6
    fig_h    = panel_w + cbar_h + 0.7
    fig      = plt.figure(figsize=(fig_w, fig_h))

    # Grid: row 0 = heatmap panels, row 1 = colorbars
    gs = fig.add_gridspec(2, n_panels,
                          height_ratios=[panel_w, cbar_h],
                          hspace=0.08, wspace=0.12,
                          left=0.05, right=0.97,
                          top=0.93, bottom=0.02)

    M = len(leaf_order)

    def _draw_heatmap(ax, matrix, cmap, vmin, vmax, title):
        """Draw one heatmap panel; returns the image object for the colorbar."""
        mat = matrix[np.ix_(leaf_order, leaf_order)].copy()
        # Replace NaN with a sentinel so we can colour them separately.
        nan_mask = np.isnan(mat)
        mat_display = np.where(nan_mask, -1.0, mat)

        # Build a masked-colormap that renders -1 as missing_color.
        cmap_m = mcolors.LinearSegmentedColormap.from_list(
            "with_missing",
            [(0, missing_color), (1e-9, cmap(0.0))] +
            [(i / 100, cmap(i / 100)) for i in range(0, 101)],
        )
        bounds = np.concatenate([[-1.5, -0.5], np.linspace(vmin, vmax, 101)])
        norm   = mcolors.BoundaryNorm(bounds, cmap_m.N)

        # For the error panels we don't want missing_color at vmin=0; use a
        # simple approach: set NaN separately as a masked array.
        img = ax.imshow(
            np.ma.masked_where(nan_mask, np.where(nan_mask, 0, mat)),
            aspect="auto", interpolation="nearest",
            cmap=cmap, vmin=vmin, vmax=vmax,
        )
        # Overlay the NaN cells in the missing colour.
        nan_overlay = np.zeros((*mat.shape, 4), dtype=float)
        c = mcolors.to_rgba(missing_color)
        nan_overlay[nan_mask] = c
        ax.imshow(nan_overlay, aspect="auto", interpolation="nearest")

        ax.set_title(title, fontsize=10, pad=4)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        return img

    panel_idx = 0

    # --- Optional reference panel: KMC exact similarity ---
    if show_reference:
        ax = fig.add_subplot(gs[0, panel_idx])
        img_ref = _draw_heatmap(
            ax, kmc_sim, cmap_ref,
            vmin=0.0, vmax=vmax_ref,
            title=f"KMC exact\n({metric_label})",
        )
        cax = fig.add_subplot(gs[1, panel_idx])
        cb  = fig.colorbar(img_ref, cax=cax, orientation="horizontal")
        cb.set_label("Similarity", fontsize=8)
        cb.ax.tick_params(labelsize=7)
        panel_idx += 1

    # --- Error panels (one per sketching method) ---
    for col, m in enumerate(methods, start=panel_idx):
        ax  = fig.add_subplot(gs[0, col])
        l1_str = f"L1 = {m['l1']:.4f}" if m.get("l1") is not None else ""
        img = _draw_heatmap(
            ax, m["err"], cmap_err,
            vmin=0.0, vmax=vmax_err,
            title=f"{m['label']}\n|error| / KMC ({metric_label})  {l1_str}",
        )

    # Shared colorbar for all error panels (spans all error-panel columns)
    err_start = panel_idx
    err_end   = n_panels - 1
    cax_err   = fig.add_subplot(gs[1, err_start:])
    cb_err    = fig.colorbar(img, cax=cax_err, orientation="horizontal")
    cb_err.set_label(f"Relative error  |estimated − exact| / exact  (0 = blue, ≥{vmax_err:.0%} = red)",
                     fontsize=8)
    cb_err.ax.tick_params(labelsize=7)
    # Add a small legend patch for the missing color
    from matplotlib.patches import Patch
    cax_err.legend(
        handles=[Patch(facecolor=missing_color, edgecolor="gray",
                       label=f"Below threshold / not measured (gray)  |  N = {M} genomes")],
        loc="upper right", fontsize=7, framealpha=0.9,
        bbox_to_anchor=(1.0, -0.6), borderaxespad=0,
    )

    stem = f"heatmap_{metric_label.lower().replace(' ', '_')}_relative_error"
    for suffix in (".pdf", ".png"):
        fpath = out_dir / (stem + suffix)
        fig.savefig(fpath)
        log.info("Saved: %s", fpath)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Argument parsing and main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot pairwise similarity error heatmaps vs KMC ground truth.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # --- method pairwise directories (all optional; include whichever have been run) ---
    p.add_argument("--amg-pairwise", default="",
                   help="AlphaMaxGeomHash pairwise directory "
                        "(pairwise_results.npz + genome_index.json)")
    p.add_argument("--bottomk-pairwise", default="",
                   help="BottomK (kmer-sketch) pairwise directory "
                        "(pairwise_results.npz + genome_index.json)")
    p.add_argument("--fracminhash-pairwise", default="",
                   help="FracMinHash (kmer-sketch) pairwise directory "
                        "(pairwise_results.npz + genome_index.json)")
    # --- Sourmash FracMinHash CSV (kept for backward-compat; off by default) ---
    p.add_argument("--sourmash-csv", default="",
                   help="Sourmash FracMinHash pairwise CSV "
                        "(gtdb_pairwise_containment.csv). "
                        "Included only if --include-sourmash is set.")
    p.add_argument("--include-sourmash", action="store_true",
                   help="Add a Sourmash FracMinHash panel (requires --sourmash-csv). "
                        "Off by default.")
    # backward-compat alias: --fracminhash was the old name for --sourmash-csv
    p.add_argument("--fracminhash", default="",
                   help=argparse.SUPPRESS)   # hidden; use --sourmash-csv instead
    # --- ground truth (required) ---
    p.add_argument("--kmc-pairwise", required=True,
                   help="KMC exact pairwise directory")
    # --- output ---
    p.add_argument("--output", required=True,
                   help="Output directory for figures")
    p.add_argument("--n-genomes", type=int, default=500,
                   help="Number of high-degree genomes to include in the heatmap")
    p.add_argument("--metric",
                   choices=["jaccard", "max_containment",
                            "containment_query_in_match",
                            "containment_match_in_query"],
                   default="jaccard",
                   help="Similarity metric to visualise")
    p.add_argument("--show-reference", action="store_true",
                   help="Add a leftmost panel showing KMC exact similarity")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    metric  = args.metric

    # ------------------------------------------------------------------
    # Load KMC ground truth (defines the reference genome index)
    # ------------------------------------------------------------------
    log.info("Loading KMC exact pairwise ...")
    kmc = load_npz_sorted(
        Path(args.kmc_pairwise) / "pairwise_results.npz",
        Path(args.kmc_pairwise) / "genome_index.json",
    )
    N = kmc["N"]

    # ------------------------------------------------------------------
    # Load each requested sketching method
    # All use load_npz_remapped so they're safely mapped to the KMC index
    # regardless of the order in their own genome_index.json.
    # ------------------------------------------------------------------

    # Helper: load an NPZ directory if the path is provided and exists
    def _load_npz_dir(path_str, label):
        if not path_str:
            return None
        d = Path(path_str)
        if not d.exists():
            log.warning("%s directory not found: %s — skipping.", label, d)
            return None
        log.info("Loading %s pairwise ...", label)
        return load_npz_remapped(
            d / "pairwise_results.npz",
            d / "genome_index.json",
            kmc["genome_id_to_int"],
            N,
        )

    amg  = _load_npz_dir(args.amg_pairwise,          "AlphaMaxGeomHash")
    bk   = _load_npz_dir(args.bottomk_pairwise,       "MinHash")
    fmh_ks = _load_npz_dir(args.fracminhash_pairwise, "FracMinHash (kmer-sketch)")

    # Sourmash FracMinHash CSV (optional; excluded by default)
    sourmash_csv = args.sourmash_csv or args.fracminhash  # backward-compat alias
    fmh_ss = None
    if args.include_sourmash and sourmash_csv:
        csv_path = Path(sourmash_csv)
        if csv_path.exists():
            log.info("Loading Sourmash FracMinHash pairwise ...")
            fmh_ss = load_fracminhash_sorted(csv_path, kmc["genome_id_to_int"], N)
        else:
            log.warning("Sourmash CSV not found: %s — skipping.", csv_path)

    # ------------------------------------------------------------------
    # Select genome subset and build dense submatrices
    # ------------------------------------------------------------------
    log.info("Selecting %d high-degree genomes ...", args.n_genomes)
    selected = select_high_degree_genomes(kmc["keys"], N, args.n_genomes)
    M        = len(selected)

    in_kmc = np.isin(kmc["keys"] // N, selected) & np.isin(kmc["keys"] % N, selected)
    fill   = in_kmc.sum()
    log.info("  KMC pairs within subset: %d / %d possible (%.1f%% fill)",
             fill, M * (M - 1) // 2, 100 * fill / (M * (M - 1) // 2))

    log.info("Building %d×%d submatrices ...", M, M)
    kmc_sim = build_submatrix(kmc["keys"], kmc[metric], selected, N)

    # ------------------------------------------------------------------
    # Build the METHODS list dynamically from whichever datasets loaded
    # Order: BottomK → AlphaMaxGeomHash → FracMinHash (increasing resources)
    # Each entry: {"label": str, "err": M×M relative error array, "l1": float}
    # ------------------------------------------------------------------
    method_specs = [
        (bk,     "MinHash\n(k=1000, kmer=31)"),
        (amg,    "AlphaMaxGeomHash\n(W=64, α=0.45, k=31)"),
        (fmh_ks, "FracMinHash\n(scale=0.01, k=31)"),
        (fmh_ss, "Sourmash FracMinHash\n(scaled=1000, k=31)"),
    ]

    # Precompute lower-triangle indices once for L1 computation
    _ii_tri, _jj_tri = np.tril_indices(M, k=-1)
    kmc_tri = kmc_sim[_ii_tri, _jj_tri]

    methods = []
    for dataset, label in method_specs:
        if dataset is None:
            continue
        sim = build_submatrix(dataset["keys"], dataset[metric], selected, N)
        err = relative_error_matrix(sim, kmc_sim)

        # Raw L1: sum of |estimate - ground_truth| over valid lower-triangle pairs
        est_tri   = sim[_ii_tri, _jj_tri]
        valid_l1  = (~np.isnan(kmc_tri)) & (~np.isnan(est_tri))
        l1_error  = float(np.sum(np.abs(est_tri[valid_l1] - kmc_tri[valid_l1])))

        methods.append({"label": label, "err": err, "l1": l1_error})

        valid = err[~np.isnan(err)]
        if len(valid):
            log.info("  %s  relative error: median=%.4f  99pct=%.4f  max=%.4f  L1=%.4f",
                     label.split("\n")[0], np.median(valid),
                     np.percentile(valid, 99), valid.max(), l1_error)

    if not methods:
        log.error("No method data loaded — provide at least one of "
                  "--amg-pairwise, --bottomk-pairwise, --fracminhash-pairwise.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Cluster using KMC similarity — apply the SAME order to all panels
    # ------------------------------------------------------------------
    log.info("Clustering genomes (Ward hierarchical clustering) ...")
    leaf_order = cluster_by_similarity(kmc_sim)
    log.info("  Clustering complete.")

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    metric_label = {
        "jaccard":                    "Jaccard",
        "max_containment":            "Max containment",
        "containment_query_in_match": "Containment (q in m)",
        "containment_match_in_query": "Containment (m in q)",
    }[metric]

    log.info("Generating heatmap figure (%s, %d panels) ...",
             metric_label, len(methods))
    plot_heatmaps(
        kmc_sim        = kmc_sim,
        methods        = methods,
        leaf_order     = leaf_order,
        metric_label   = metric_label,
        out_dir        = out_dir,
        show_reference = args.show_reference,
    )

    # Save L1 errors so 10_plot_resources.py can build the accuracy bar chart.
    # We load+update rather than overwrite so that jaccard and max_containment
    # runs accumulate into the same file.
    l1_json_path = out_dir / "l1_errors.json"
    existing_l1 = {}
    if l1_json_path.exists():
        with open(l1_json_path) as f:
            existing_l1 = json.load(f)
    existing_l1[metric] = {
        m["label"].split("\n")[0]: m["l1"] for m in methods
    }
    with open(l1_json_path, "w") as f:
        json.dump(existing_l1, f, indent=2)
    log.info("L1 errors written to %s", l1_json_path)

    # Save the genome ordering so other scripts can reuse it
    ordering_path = out_dir / f"genome_subset_order_{metric}.json"
    with open(ordering_path, "w") as f:
        json.dump({
            "metric":         metric,
            "n_genomes":      M,
            "genome_indices": selected.tolist(),
            "leaf_order":     leaf_order.tolist(),
            "genome_ids":     [kmc["genome_ids"][i] for i in selected],
        }, f, indent=2)
    log.info("Genome ordering saved to %s", ordering_path)
    log.info("Done.")


if __name__ == "__main__":
    main()
