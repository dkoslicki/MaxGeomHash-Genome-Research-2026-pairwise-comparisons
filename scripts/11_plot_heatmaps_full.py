#!/usr/bin/env python3
"""
11_plot_heatmaps_full.py
========================
Full-scale pairwise absolute-error heatmaps for all 143,614 GTDB genomes,
rendered with Datashader so that no dense M×M matrix is ever allocated.

Every pair present in the KMC pairwise NPZ is used; nothing is subsampled.
Genomes are ordered on both axes by degree (number of KMC similar pairs),
placing the most-connected genomes at the top-left corner and producing
visible block structure at the resolution Datashader renders (default
2000×2000 pixels per panel, covering ~72×72 genomes per pixel).

Layout and visual style match 10_plot_heatmaps.py:
  - Same RdYlBu_r error colormap, same gray missing-pair colour (#CCCCCC)
  - Same multi-panel layout with shared colorbar
  - Per-panel title includes method label and L1 error over ALL pairs
  - Outputs are named *_full.{pdf,png} so they don't overwrite the 500-genome plots

Usage (activate the sourmash conda env which has datashader):
    conda run -n sourmash python3 scripts/11_plot_heatmaps_full.py \\
        --amg-pairwise         data/GTDB/alphamaxgeom_pairwise \\
        --bottomk-pairwise     data/GTDB/bottomk_pairwise \\
        --fracminhash-pairwise data/GTDB/fracminhash_pairwise \\
        --kmc-pairwise         data/GTDB/kmc_pairwise \\
        --output               data/GTDB/figures

Optional flags:
    --metric jaccard|max_containment   (default: jaccard)
    --ordering degree|rcm|spectral     (default: spectral)
        degree   — sort by KMC pair count; fast, coarse
        rcm      — Reverse Cuthill-McKee bandwidth minimisation; ~seconds
        spectral — Fiedler vector of the normalised Laplacian; ~1-5 min;
                   produces block structure closest to Ward's clustering
    --canvas-size 2000                 (pixels per axis per panel, default 2000)
    --show-reference                   (add a leftmost KMC similarity panel)
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

try:
    # datashader 0.13 uses two APIs removed in Python 3.11+; patch before import.
    import inspect as _inspect
    import warnings as _warnings
    import numpy as _np
    if not hasattr(_np, "warnings"):
        _np.warnings = _warnings                      # removed in NumPy 1.25
    if not hasattr(_inspect, "getargspec"):
        _inspect.getargspec = _inspect.getfullargspec  # removed in Python 3.11
    import datashader as ds
except ImportError:
    sys.exit(
        "datashader is required but not found.\n"
        "Run this script inside the sourmash conda env:\n"
        "  conda run -n sourmash python3 scripts/11_plot_heatmaps_full.py ..."
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Publication style (mirrors 10_plot_heatmaps.py)
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "font.size":       9,
    "axes.titlesize":  10,
    "axes.labelsize":  9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "figure.dpi":      150,
    "savefig.dpi":     300,
    "savefig.bbox":    "tight",
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
})

MISSING_COLOR = "#CCCCCC"

# ---------------------------------------------------------------------------
# Integer pair-key encoding (same convention as 10_plot_heatmaps.py)
# ---------------------------------------------------------------------------

def encode_pair_keys(rows: np.ndarray, cols: np.ndarray, N: int) -> np.ndarray:
    r  = rows.astype(np.int64)
    c  = cols.astype(np.int64)
    lo = np.minimum(r, c)
    hi = np.maximum(r, c)
    return lo * N + hi


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_npz(npz_path: Path, index_path: Path) -> dict:
    """Load a KMC-format NPZ and return sorted canonical keys + metric arrays."""
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

    def _get(name):
        return data[name][order].astype(np.float32) if name in data else None

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
                      ref_id_to_int: dict, N: int) -> dict:
    """Load a sketch-method NPZ and remap genome indices to the KMC reference index."""
    if not npz_path.exists() or not index_path.exists():
        log.error("File not found: %s or %s", npz_path, index_path)
        sys.exit(1)

    local_ids    = json.load(open(index_path))["genomes"]
    local_to_ref = np.array([ref_id_to_int.get(gid, -1) for gid in local_ids],
                             dtype=np.int64)

    data     = np.load(npz_path)
    raw_rows = local_to_ref[data["row"].astype(np.int64)]
    raw_cols = local_to_ref[data["col"].astype(np.int64)]
    valid    = (raw_rows >= 0) & (raw_cols >= 0)
    raw_rows = raw_rows[valid]
    raw_cols = raw_cols[valid]

    keys  = encode_pair_keys(raw_rows, raw_cols, N)
    order = np.argsort(keys, kind="stable")
    keys  = keys[order]

    def _get(name):
        return data[name][valid][order].astype(np.float32) if name in data else None

    log.info("  Loaded %d pairs (remapped) from %s", len(keys), npz_path.name)
    return {
        "keys":                       keys,
        "jaccard":                    _get("jaccard"),
        "max_containment":            _get("max_containment"),
        "containment_query_in_match": _get("containment_query_in_match"),
        "containment_match_in_query": _get("containment_match_in_query"),
    }


# ---------------------------------------------------------------------------
# Genome ordering
# ---------------------------------------------------------------------------

def _build_sparse(kmc_keys: np.ndarray, kmc_vals: np.ndarray, N: int):
    """Build a symmetric scipy CSR similarity matrix from canonical pair keys."""
    from scipy.sparse import csr_matrix
    rows = (kmc_keys // N).astype(np.int32)
    cols = (kmc_keys  % N).astype(np.int32)
    vals = kmc_vals.astype(np.float32)
    r = np.concatenate([rows, cols])
    c = np.concatenate([cols, rows])
    v = np.tile(vals, 2)
    return csr_matrix((v, (r, c)), shape=(N, N), dtype=np.float32)


def degree_order(kmc_keys: np.ndarray, N: int) -> np.ndarray:
    """
    Rank genomes by degree (number of KMC similar pairs), highest first.
    Fast (~seconds) but only captures coarse connectivity structure.
    Returns rank: int32 array of length N.
    """
    rows   = (kmc_keys // N).astype(np.int64)
    cols   = (kmc_keys  % N).astype(np.int64)
    degree = np.zeros(N, dtype=np.int32)
    np.add.at(degree, rows, 1)
    np.add.at(degree, cols, 1)
    order = np.argsort(degree)[::-1].astype(np.int32)
    rank  = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


def rcm_order(kmc_keys: np.ndarray, kmc_vals: np.ndarray, N: int) -> np.ndarray:
    """
    Reverse Cuthill-McKee ordering.  Minimises the bandwidth of the sparse
    similarity matrix (keeps edges close to the diagonal) using a BFS from
    the lowest-degree node.  Runs in ~seconds on 143k nodes.
    Returns rank: int32 array of length N.
    """
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    A     = _build_sparse(kmc_keys, kmc_vals, N)
    order = reverse_cuthill_mckee(A, symmetric_mode=True).astype(np.int32)
    rank  = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


def spectral_order(kmc_keys: np.ndarray, kmc_vals: np.ndarray, N: int) -> np.ndarray:
    """
    Fiedler-vector ordering: sort genomes along the second eigenvector of the
    normalised graph Laplacian.  This is the optimal 1-D relaxation of the
    graph bisection problem and produces block-diagonal structure closely
    mimicking Ward's hierarchical clustering, but operating directly on the
    sparse KMC similarity graph without ever materialising a dense matrix.

    Uses LOBPCG (locally optimal block preconditioned conjugate gradient) with
    a diagonal (Jacobi) preconditioner, which converges robustly for graph
    Laplacians without requiring a sparse factorisation.  Typical runtime on
    143k nodes / ~2M edges: 1–5 minutes.

    Returns rank: int32 array of length N.
    """
    from scipy.sparse import diags
    from scipy.sparse.csgraph import laplacian as sp_laplacian
    from scipy.sparse.linalg import lobpcg

    log.info("  Building sparse Laplacian (%d nodes) ...", N)
    A = _build_sparse(kmc_keys, kmc_vals, N)
    L = sp_laplacian(A, normed=True)          # normalised Laplacian, eigenvalues in [0,2]

    # Jacobi (diagonal) preconditioner: M ≈ diag(L)^{-1}
    diag = np.array(L.diagonal(), dtype=np.float64)
    diag[diag < 1e-10] = 1.0                  # avoid divide-by-zero for isolated nodes
    M = diags(1.0 / diag, format="csr")

    # Initial guess: two random orthonormal vectors
    rng = np.random.default_rng(42)
    X   = rng.standard_normal((N, 2))
    X, _ = np.linalg.qr(X)
    X   = X.astype(np.float64)

    log.info("  Running LOBPCG eigensolver (k=2, tol=1e-4) ...")
    vals_eig, vecs = lobpcg(L.astype(np.float64), X, M=M,
                             largest=False, tol=1e-4, maxiter=400,
                             verbosityLevel=0)

    # Take the eigenvector for the second-smallest eigenvalue (Fiedler vector)
    fiedler = vecs[:, np.argsort(vals_eig)[1]]
    log.info("  Fiedler eigenvalue: %.6f  (0 = disconnected graph)", vals_eig.min())

    order = np.argsort(fiedler).astype(np.int32)
    rank  = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


# ---------------------------------------------------------------------------
# Absolute-error computation → Datashader-ready DataFrame
# ---------------------------------------------------------------------------

def make_error_dataframe(
    kmc_keys:   np.ndarray,
    kmc_vals:   np.ndarray,
    meth_keys:  np.ndarray,
    meth_vals:  np.ndarray,
    rank:       np.ndarray,
    N:          int,
    min_gt:     float = 1e-6,
) -> tuple:
    """
    Intersect the method's pairs with the KMC pairs, compute per-pair
    absolute error |estimate − exact|, and return:

      df    -- pandas DataFrame with columns (x, y, err), float32.
               Both triangles are included so the heatmap is symmetric.
               (x, y) are genome plot ranks; err is the absolute error.
      l1    -- float: sum of absolute errors over the matched pairs
               (each canonical pair counted once, not doubled).

    Datashader's Canvas.points() will aggregate these into a pixel grid
    using the mean of 'err' in each pixel.
    """
    # Find pairs present in both KMC (sorted) and method (sorted)
    pos      = np.searchsorted(kmc_keys, meth_keys)
    in_range = pos < len(kmc_keys)
    pos_safe = np.where(in_range, pos, 0)
    matched  = in_range & (kmc_keys[pos_safe] == meth_keys)

    m_vals = meth_vals[matched]
    k_vals = kmc_vals[pos_safe[matched]]
    keys_m = meth_keys[matched]

    valid = (k_vals > min_gt) & np.isfinite(m_vals) & np.isfinite(k_vals)
    err   = np.abs(m_vals[valid] - k_vals[valid])
    keys_v = keys_m[valid]

    l1   = float(err.sum())
    rows = (keys_v // N).astype(np.int32)   # canonical lo index
    cols = (keys_v  % N).astype(np.int32)   # canonical hi index

    # Both triangles: (rank[rows], rank[cols]) and its transpose
    x = np.concatenate([rank[rows], rank[cols]]).astype(np.float32)
    y = np.concatenate([rank[cols], rank[rows]]).astype(np.float32)
    e = np.tile(err.astype(np.float32), 2)

    return pd.DataFrame({"x": x, "y": y, "err": e}), l1


def make_similarity_dataframe(
    kmc_keys: np.ndarray,
    kmc_vals: np.ndarray,
    rank:     np.ndarray,
    N:        int,
) -> pd.DataFrame:
    """Build a symmetric (x, y, val) DataFrame for the KMC reference panel."""
    rows = (kmc_keys // N).astype(np.int32)
    cols = (kmc_keys  % N).astype(np.int32)
    x = np.concatenate([rank[rows], rank[cols]]).astype(np.float32)
    y = np.concatenate([rank[cols], rank[rows]]).astype(np.float32)
    v = np.tile(kmc_vals.astype(np.float32), 2)
    return pd.DataFrame({"x": x, "y": y, "val": v})


# ---------------------------------------------------------------------------
# Datashader rasterisation → RGBA image
# ---------------------------------------------------------------------------

def rasterise(
    df:         pd.DataFrame,
    col:        str,
    N:          int,
    canvas_size: int,
    vmax:       float,
    lut_rgba:   np.ndarray,   # uint8 (256, 4) RGBA lookup table
) -> np.ndarray:
    """
    Aggregate df[col] into a canvas_size × canvas_size grid using
    Datashader (mean reduction), map through lut_rgba, and paint
    empty pixels with MISSING_COLOR.

    Returns a uint8 RGBA array of shape (canvas_size, canvas_size, 4)
    suitable for matplotlib's imshow.

    Coordinate convention
    ---------------------
    Datashader places y_range[0] at the BOTTOM of the canvas; its output
    array's row 0 corresponds to y_range[0].  matplotlib imshow with the
    default origin="upper" then displays row 0 at the TOP of the plot.
    With y_range=(0, N-1) and rank=0 assigned to the highest-degree genome,
    rank=0 therefore appears at the top-left corner of the image — the same
    convention as the Ward-clustered dense heatmaps in 10_plot_heatmaps.py.
    """
    cvs = ds.Canvas(
        plot_width=canvas_size,
        plot_height=canvas_size,
        x_range=(0, N - 1),
        y_range=(0, N - 1),
    )
    agg  = cvs.points(df, "x", "y", ds.mean(col))
    vals = agg.values.astype(np.float32)    # (canvas_size, canvas_size), NaN = empty

    idx  = np.clip((vals / vmax * 255).astype(np.int32), 0, 255)
    rgba = lut_rgba[idx].copy()             # (h, w, 4)

    missing = (np.array(mcolors.to_rgba(MISSING_COLOR)) * 255).astype(np.uint8)
    rgba[np.isnan(vals)] = missing

    return rgba.astype(np.uint8)


# ---------------------------------------------------------------------------
# Figure composition
# ---------------------------------------------------------------------------

def plot_heatmaps_full(
    kmc_keys:      np.ndarray,
    kmc_sim_vals:  np.ndarray,
    rank:          np.ndarray,
    methods:       list,        # [{"label": str, "df": DataFrame, "l1": float}, ...]
    N:             int,
    canvas_size:   int,
    vmax_err:      float,
    metric_label:  str,
    ordering_label: str,
    out_dir:       Path,
    show_reference: bool = False,
):
    cmap_err = plt.cm.RdYlBu_r
    cmap_ref = plt.cm.viridis

    # Build 8-bit RGBA look-up tables from matplotlib colormaps
    lut_err = (cmap_err(np.linspace(0, 1, 256)) * 255).astype(np.uint8)
    lut_ref = (cmap_ref(np.linspace(0, 1, 256)) * 255).astype(np.uint8)

    n_methods = len(methods)
    n_panels  = n_methods + (1 if show_reference else 0)

    panel_w = 3.8
    cbar_h  = 0.35
    fig_w   = panel_w * n_panels + 0.6
    fig_h   = panel_w + cbar_h + 0.7
    fig     = plt.figure(figsize=(fig_w, fig_h))

    gs = fig.add_gridspec(
        2, n_panels,
        height_ratios=[panel_w, cbar_h],
        hspace=0.08, wspace=0.12,
        left=0.05, right=0.97,
        top=0.93,  bottom=0.02,
    )

    err_col_start = 0    # column index where error panels begin

    # ---- Optional KMC reference panel ----------------------------------------
    if show_reference:
        log.info("  Rasterising KMC reference panel ...")
        vmax_ref = float(np.nanmax(kmc_sim_vals)) or 1.0
        ref_df   = make_similarity_dataframe(kmc_keys, kmc_sim_vals, rank, N)
        rgba_ref = rasterise(ref_df, "val", N, canvas_size, vmax_ref, lut_ref)

        ax = fig.add_subplot(gs[0, 0])
        ax.imshow(rgba_ref, aspect="auto", interpolation="nearest")
        ax.set_title(f"KMC exact\n({metric_label})", fontsize=10, pad=4)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)

        sm_ref = plt.cm.ScalarMappable(
            cmap=cmap_ref, norm=mcolors.Normalize(0, vmax_ref))
        sm_ref.set_array([])
        cax = fig.add_subplot(gs[1, 0])
        cb  = fig.colorbar(sm_ref, cax=cax, orientation="horizontal")
        cb.set_label("Similarity", fontsize=8)
        cb.ax.tick_params(labelsize=7)

        err_col_start = 1

    # ---- Error panels (one per method) ----------------------------------------
    for col_idx, m in enumerate(methods, start=err_col_start):
        log.info("  Rasterising %s ...", m["label"].split("\n")[0])
        rgba = rasterise(m["df"], "err", N, canvas_size, vmax_err, lut_err)

        ax = fig.add_subplot(gs[0, col_idx])
        ax.imshow(rgba, aspect="auto", interpolation="nearest")
        ax.set_title(
            f"{m['label']}\n|error| vs KMC ({metric_label})",
            fontsize=10, pad=4,
        )
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)

    # ---- Shared colorbar for all error panels ---------------------------------
    sm_err = plt.cm.ScalarMappable(
        cmap=cmap_err, norm=mcolors.Normalize(0, vmax_err))
    sm_err.set_array([])
    cax_err = fig.add_subplot(gs[1, err_col_start:])
    cb_err  = fig.colorbar(sm_err, cax=cax_err, orientation="horizontal")
    cb_err.set_label(
        f"Absolute error  |estimated − exact|"
        f"  (0 = blue, ≥{vmax_err:.4f} = red)",
        fontsize=8,
    )
    cb_err.ax.tick_params(labelsize=7)
    cax_err.legend(
        handles=[Patch(
            facecolor=MISSING_COLOR, edgecolor="gray",
            label=f"Below threshold (gray)  |  N = {N:,} genomes",
        )],
        loc="lower left", fontsize=7, framealpha=0.9,
        bbox_to_anchor=(0.0, 1.02), borderaxespad=0,
    )

    metric_stem = metric_label.lower().replace(" ", "_")
    stem = f"heatmap_{metric_stem}_absolute_error_full_{ordering_label}"
    for suffix in (".pdf", ".png"):
        fpath = out_dir / (stem + suffix)
        fig.savefig(fpath)
        log.info("Saved: %s", fpath)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Full-scale pairwise error heatmaps via Datashader.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--amg-pairwise",         default="",
                   help="AlphaMaxGeomHash pairwise directory")
    p.add_argument("--bottomk-pairwise",     default="",
                   help="BottomK/MinHash pairwise directory")
    p.add_argument("--fracminhash-pairwise", default="",
                   help="FracMinHash pairwise directory")
    p.add_argument("--kmc-pairwise",         required=True,
                   help="KMC exact pairwise directory")
    p.add_argument("--output",               required=True,
                   help="Output directory for figures")
    p.add_argument("--metric",
                   choices=["jaccard", "max_containment",
                             "containment_query_in_match",
                             "containment_match_in_query"],
                   default="jaccard",
                   help="Similarity metric to visualise")
    p.add_argument("--ordering",
                   choices=["degree", "rcm", "spectral"],
                   default="spectral",
                   help="Genome axis ordering: degree (fast, coarse), "
                        "rcm (Reverse Cuthill-McKee, ~seconds), or "
                        "spectral (Fiedler vector, ~minutes, best structure)")
    p.add_argument("--canvas-size", type=int, default=2000,
                   help="Pixels per axis per panel")
    p.add_argument("--show-reference", action="store_true",
                   help="Add a leftmost KMC similarity panel")
    return p.parse_args()


def main():
    args    = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    metric  = args.metric

    # ---- Load KMC ground truth ------------------------------------------------
    log.info("Loading KMC pairwise ...")
    kmc = load_npz(
        Path(args.kmc_pairwise) / "pairwise_results.npz",
        Path(args.kmc_pairwise) / "genome_index.json",
    )
    N         = kmc["N"]
    kmc_vals  = kmc[metric]
    if kmc_vals is None:
        log.error("Metric '%s' not found in KMC NPZ.", metric)
        sys.exit(1)

    # ---- Genome ordering ------------------------------------------------------
    log.info("Computing genome ordering: %s (%d genomes) ...", args.ordering, N)
    if args.ordering == "spectral":
        rank = spectral_order(kmc["keys"], kmc_vals, N)
    elif args.ordering == "rcm":
        rank = rcm_order(kmc["keys"], kmc_vals, N)
    else:
        rank = degree_order(kmc["keys"], N)
    log.info("  Ordering done.")

    # ---- Load sketch methods --------------------------------------------------
    def _load_dir(path_str, label):
        if not path_str:
            return None
        d = Path(path_str)
        if not d.exists():
            log.warning("%s directory not found: %s — skipping.", label, d)
            return None
        log.info("Loading %s ...", label)
        return load_npz_remapped(
            d / "pairwise_results.npz",
            d / "genome_index.json",
            kmc["genome_id_to_int"], N,
        )

    datasets = [
        (_load_dir(args.bottomk_pairwise,     "MinHash"),          "MinHash\n(k=1000, kmer=31)"),
        (_load_dir(args.amg_pairwise,          "AlphaMaxGeomHash"), "AlphaMaxGeomHash\n(W=64, α=0.45, k=31)"),
        (_load_dir(args.fracminhash_pairwise,  "FracMinHash"),      "FracMinHash\n(scale=0.01, k=31)"),
    ]

    # ---- Compute absolute errors ----------------------------------------------
    log.info("Computing absolute errors ...")
    methods          = []
    all_err_samples  = []

    for dataset, label in datasets:
        if dataset is None:
            continue
        meth_vals = dataset[metric]
        if meth_vals is None:
            log.warning("%s: metric '%s' not in NPZ — skipping.", label, metric)
            continue

        df, l1 = make_error_dataframe(
            kmc["keys"], kmc_vals,
            dataset["keys"], meth_vals,
            rank, N,
        )
        n_pairs = len(df) // 2
        log.info("  %-30s %10d pairs   L1 = %.4f",
                 label.split("\n")[0], n_pairs, l1)

        methods.append({"label": label, "df": df, "l1": l1})

        # Sample up to 1 M values for vmax estimation (avoid huge alloc)
        err_col = df["err"].values
        step    = max(1, len(err_col) // 1_000_000)
        all_err_samples.append(err_col[::step])

    if not methods:
        log.error("No method data available — provide at least one of "
                  "--amg-pairwise, --bottomk-pairwise, --fracminhash-pairwise.")
        sys.exit(1)

    # Shared colour scale: 99th percentile across all methods, capped at 100 %
    all_err  = np.concatenate(all_err_samples)
    vmax_err = float(np.percentile(all_err, 99))
    vmax_err = min(vmax_err, 1.0)
    log.info("Shared vmax_err (99th percentile): %.4f", vmax_err)

    # ---- Render ---------------------------------------------------------------
    metric_label = {
        "jaccard":                    "Jaccard",
        "max_containment":            "Max containment",
        "containment_query_in_match": "Containment (q in m)",
        "containment_match_in_query": "Containment (m in q)",
    }[metric]

    log.info("Rendering %d panel(s) at %d×%d pixels ...",
             len(methods), args.canvas_size, args.canvas_size)
    plot_heatmaps_full(
        kmc_keys      = kmc["keys"],
        kmc_sim_vals  = kmc_vals,
        rank          = rank,
        methods       = methods,
        N             = N,
        canvas_size   = args.canvas_size,
        vmax_err      = vmax_err,
        metric_label  = metric_label,
        ordering_label = args.ordering,
        out_dir       = out_dir,
        show_reference = args.show_reference,
    )

    # Write mean L1 to l1_errors.json (same file as 10_plot_heatmaps.py)
    l1_json_path = out_dir / "l1_errors.json"
    existing_l1  = {}
    if l1_json_path.exists():
        with open(l1_json_path) as f:
            existing_l1 = json.load(f)
    existing_l1[metric] = {
        m["label"].split("\n")[0]: (
            m["l1"] / (len(m["df"]) // 2) if len(m["df"]) > 0 else 0.0
        )
        for m in methods
    }
    existing_l1[f"{metric}__l1_mode"] = "full"
    with open(l1_json_path, "w") as f:
        json.dump(existing_l1, f, indent=2)
    log.info("Mean L1 errors written to %s", l1_json_path)
    log.info("Done.")


if __name__ == "__main__":
    main()
