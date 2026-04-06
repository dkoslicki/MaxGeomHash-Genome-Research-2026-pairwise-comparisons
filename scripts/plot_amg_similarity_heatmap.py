#!/usr/bin/env python3
"""
plot_amg_similarity_heatmap.py
==============================
One-off: render the raw AlphaMaxGeomHash Jaccard and max_containment similarity
values as full-scale Datashader heatmaps with spectral (Fiedler-vector) genome
ordering.

Unlike the publication heatmaps (11_plot_heatmaps_full.py), this shows the
*estimates themselves* — not errors relative to KMC exact — so the colour scale
runs from 0 to the 90th-percentile observed value and the missing-pair colour
is a neutral grey.

Usage (run inside the sourmash conda env, which ships datashader):

    conda run -n sourmash python3 scripts/plot_amg_similarity_heatmap.py

Optional flags:
    --amg-pairwise   DIR     (default: data/GTDB/alphamaxgeom_pairwise_thr0001)
    --output         DIR     (default: data/GTDB/figures_thr0001)
    --ordering       spectral|rcm|degree   (default: spectral)
    --canvas-size    INT     pixels per axis (default: 2000)
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
    import inspect as _inspect
    import warnings as _warnings
    if not hasattr(np, "warnings"):
        np.warnings = _warnings
    if not hasattr(_inspect, "getargspec"):
        _inspect.getargspec = _inspect.getfullargspec
    import datashader as ds
except ImportError:
    sys.exit(
        "datashader is required but not found.\n"
        "Run inside the sourmash conda env:\n"
        "  conda run -n sourmash python3 scripts/plot_amg_similarity_heatmap.py"
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

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
# Data loading
# ---------------------------------------------------------------------------

def encode_pair_keys(rows: np.ndarray, cols: np.ndarray, N: int) -> np.ndarray:
    r  = rows.astype(np.int64)
    c  = cols.astype(np.int64)
    return np.minimum(r, c) * N + np.maximum(r, c)


def load_npz(npz_path: Path, index_path: Path) -> dict:
    """Load a kmer-sketch pairwise NPZ; return sorted canonical keys + metric arrays."""
    if not npz_path.exists() or not index_path.exists():
        log.error("File not found: %s or %s", npz_path, index_path)
        sys.exit(1)

    genome_ids = json.load(open(index_path))["genomes"]
    N          = len(genome_ids)
    data       = np.load(npz_path)
    keys       = encode_pair_keys(data["row"], data["col"], N)
    order      = np.argsort(keys, kind="stable")
    keys       = keys[order]

    def _get(name):
        return data[name][order].astype(np.float32) if name in data else None

    log.info("  Loaded %d pairs from %s  (%d genomes in index)",
             len(keys), npz_path.name, N)
    return {
        "genome_ids": genome_ids,
        "N":          N,
        "keys":       keys,
        "jaccard":          _get("jaccard"),
        "max_containment":  _get("max_containment"),
    }


# ---------------------------------------------------------------------------
# Genome ordering
# ---------------------------------------------------------------------------

def _build_sparse(keys: np.ndarray, vals: np.ndarray, N: int):
    from scipy.sparse import csr_matrix
    rows = (keys // N).astype(np.int32)
    cols = (keys  % N).astype(np.int32)
    v    = vals.astype(np.float32)
    r    = np.concatenate([rows, cols])
    c    = np.concatenate([cols, rows])
    vv   = np.tile(v, 2)
    return csr_matrix((vv, (r, c)), shape=(N, N), dtype=np.float32)


def spectral_order(keys: np.ndarray, vals: np.ndarray, N: int) -> np.ndarray:
    """Sort genomes by the Fiedler vector of the normalised Laplacian."""
    from scipy.sparse import diags
    from scipy.sparse.csgraph import laplacian as sp_laplacian
    from scipy.sparse.linalg import lobpcg

    log.info("  Building sparse Laplacian (%d nodes) ...", N)
    A = _build_sparse(keys, vals, N)
    L = sp_laplacian(A, normed=True)

    diag = np.array(L.diagonal(), dtype=np.float64)
    diag[diag < 1e-10] = 1.0
    M = diags(1.0 / diag, format="csr")

    rng     = np.random.default_rng(42)
    X, _    = np.linalg.qr(rng.standard_normal((N, 2)))
    X       = X.astype(np.float64)

    log.info("  Running LOBPCG eigensolver (k=2, tol=1e-4) ...")
    vals_eig, vecs = lobpcg(L.astype(np.float64), X, M=M,
                             largest=False, tol=1e-4, maxiter=400,
                             verbosityLevel=0)

    fiedler = vecs[:, np.argsort(vals_eig)[1]]
    log.info("  Fiedler eigenvalue: %.6f", sorted(vals_eig)[1])

    order = np.argsort(fiedler).astype(np.int32)
    rank  = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


def rcm_order(keys: np.ndarray, vals: np.ndarray, N: int) -> np.ndarray:
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    A     = _build_sparse(keys, vals, N)
    order = reverse_cuthill_mckee(A, symmetric_mode=True).astype(np.int32)
    rank  = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


def degree_order(keys: np.ndarray, N: int) -> np.ndarray:
    rows   = (keys // N).astype(np.int64)
    cols   = (keys  % N).astype(np.int64)
    degree = np.zeros(N, dtype=np.int32)
    np.add.at(degree, rows, 1)
    np.add.at(degree, cols, 1)
    order = np.argsort(degree)[::-1].astype(np.int32)
    rank  = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


# ---------------------------------------------------------------------------
# Datashader rasterisation
# ---------------------------------------------------------------------------

def make_similarity_dataframe(
    keys: np.ndarray, vals: np.ndarray, rank: np.ndarray, N: int
) -> pd.DataFrame:
    """Build a symmetric (x, y, val) DataFrame ready for Datashader."""
    rows = (keys // N).astype(np.int32)
    cols = (keys  % N).astype(np.int32)
    x = np.concatenate([rank[rows], rank[cols]]).astype(np.float32)
    y = np.concatenate([rank[cols], rank[rows]]).astype(np.float32)
    v = np.tile(vals.astype(np.float32), 2)
    return pd.DataFrame({"x": x, "y": y, "val": v})


def rasterise(
    df: pd.DataFrame, col: str, N: int, canvas_size: int,
    vmax: float, lut_rgba: np.ndarray,
) -> np.ndarray:
    """Datashader mean-reduction → uint8 RGBA array (canvas_size × canvas_size × 4)."""
    cvs = ds.Canvas(
        plot_width=canvas_size, plot_height=canvas_size,
        x_range=(0, N - 1), y_range=(0, N - 1),
    )
    agg  = cvs.points(df, "x", "y", ds.mean(col))
    vals = agg.values.astype(np.float32)

    idx  = np.clip(np.where(np.isnan(vals), 0, vals / vmax * 255).astype(np.int32), 0, 255)
    rgba = lut_rgba[idx].copy()

    missing_rgba = (np.array(mcolors.to_rgba(MISSING_COLOR)) * 255).astype(np.uint8)
    rgba[np.isnan(vals)] = missing_rgba

    return rgba.astype(np.uint8)


# ---------------------------------------------------------------------------
# Figure composition
# ---------------------------------------------------------------------------

def plot_similarity_heatmaps(
    pairwise:     dict,       # from load_npz()
    rank:         np.ndarray,
    canvas_size:  int,
    ordering_label: str,
    out_dir:      Path,
):
    """
    Two-panel figure: AlphaMaxGeomHash Jaccard (left) | max_containment (right).
    Colour scale: viridis, 0 → 90th-percentile observed value.
    """
    N     = pairwise["N"]
    keys  = pairwise["keys"]
    cmap  = plt.cm.viridis
    lut   = (cmap(np.linspace(0, 1, 256)) * 255).astype(np.uint8)

    panels = []
    for metric, label in [("jaccard", "Jaccard"), ("max_containment", "Max containment")]:
        vals = pairwise[metric]
        if vals is None:
            log.warning("Metric '%s' not found in NPZ — skipping.", metric)
            continue
        vmax = float(np.nanpercentile(vals, 90)) or 1.0
        log.info("Building %s DataFrame (vmax=%.4f) ...", label, vmax)
        df = make_similarity_dataframe(keys, vals, rank, N)
        panels.append({"label": label, "metric": metric, "df": df, "vmax": vmax})

    if not panels:
        log.error("No metrics available; nothing to plot.")
        return

    n_panels = len(panels)
    panel_w  = 4.0
    cbar_h   = 0.35
    fig_w    = panel_w * n_panels + 0.5
    fig_h    = panel_w + cbar_h + 0.7

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = fig.add_gridspec(
        2, n_panels,
        height_ratios=[panel_w, cbar_h],
        hspace=0.08, wspace=0.12,
        left=0.04, right=0.97,
        top=0.93,  bottom=0.02,
    )

    for col_idx, p in enumerate(panels):
        log.info("  Rasterising %s (%d×%d) ...", p["label"], canvas_size, canvas_size)
        rgba = rasterise(p["df"], "val", N, canvas_size, p["vmax"], lut)

        ax = fig.add_subplot(gs[0, col_idx])
        ax.imshow(rgba, aspect="auto", interpolation="nearest")
        ax.set_title(
            f"AlphaMaxGeomHash  —  {p['label']}\n"
            f"spectral ordering  |  {N:,} genomes  |  {len(keys):,} pairs above threshold",
            fontsize=9, pad=4,
        )
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)

        sm  = plt.cm.ScalarMappable(
            cmap=cmap, norm=mcolors.Normalize(0, p["vmax"]))
        sm.set_array([])
        cax = fig.add_subplot(gs[1, col_idx])
        cb  = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(f"{p['label']} (0 → {p['vmax']:.4f} = 90th pctile)", fontsize=8)
        cb.ax.tick_params(labelsize=7)
        cax.legend(
            handles=[Patch(
                facecolor=MISSING_COLOR, edgecolor="gray",
                label="Below threshold / not computed (grey)",
            )],
            loc="lower left", fontsize=7, framealpha=0.9,
            bbox_to_anchor=(0.0, 1.02), borderaxespad=0,
        )

    stem = f"amg_similarity_heatmap_{ordering_label}"
    for suffix in (".pdf", ".png"):
        fpath = out_dir / (stem + suffix)
        fig.savefig(fpath)
        log.info("Saved: %s", fpath)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    base = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(
        description="AMG similarity heatmaps (Jaccard + max_containment) via Datashader.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--amg-pairwise",
                   default=str(base / "data/GTDB/alphamaxgeom_pairwise_thr0001"),
                   help="AlphaMaxGeomHash pairwise directory")
    p.add_argument("--output",
                   default=str(base / "data/GTDB/figures_thr0001"),
                   help="Output directory for figures")
    p.add_argument("--ordering",
                   choices=["spectral", "rcm", "degree"],
                   default="spectral",
                   help="Genome axis ordering")
    p.add_argument("--canvas-size", type=int, default=2000,
                   help="Pixels per axis per panel")
    return p.parse_args()


def main():
    args    = parse_args()
    amg_dir = Path(args.amg_pairwise)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading AlphaMaxGeomHash pairwise from %s ...", amg_dir)
    pw = load_npz(amg_dir / "pairwise_results.npz", amg_dir / "genome_index.json")

    log.info("Computing genome ordering: %s (%d genomes) ...",
             args.ordering, pw["N"])
    if args.ordering == "spectral":
        rank = spectral_order(pw["keys"], pw["jaccard"], pw["N"])
    elif args.ordering == "rcm":
        rank = rcm_order(pw["keys"], pw["jaccard"], pw["N"])
    else:
        rank = degree_order(pw["keys"], pw["N"])
    log.info("  Ordering done.")

    plot_similarity_heatmaps(
        pairwise       = pw,
        rank           = rank,
        canvas_size    = args.canvas_size,
        ordering_label = args.ordering,
        out_dir        = out_dir,
    )
    log.info("Done.")


if __name__ == "__main__":
    main()
