#!/usr/bin/env python3
"""
plot_amg_hierarchical_heatmap.py
=================================
One-off: render the raw AlphaMaxGeomHash Jaccard and max_containment similarity
values as full-scale Datashader heatmaps with hierarchical-clustering genome
ordering.

Five ordering strategies are available (--ordering flag):

  spectral-ward  (default)
      Compute the top-K singular vectors of the symmetrically-normalised
      adjacency matrix D^{-1/2} A D^{-1/2} via randomised SVD (O(N×k×iters)).
      These span the same space as the Fiedler eigenvectors of the normalised
      Laplacian but are robust to graphs with many isolated components (LOBPCG
      on the Laplacian would return only trivial zero-eigenvalue vectors for
      the GTDB graph at threshold=0.001).  Ward agglomerative clustering is
      then applied to the K-D embedding with the sparse AMG similarity graph
      as the connectivity constraint (O(n log n)).

  umap-ward
      Same spectral embedding as above, refined with UMAP (umap-learn) into
      a lower-dimensional non-linear embedding.  UMAP captures manifold
      structure that linear methods miss.  Ward + connectivity is then applied
      to the UMAP features.  Typically the sharpest cluster blocks.
      Requires: umap-learn (pynndescent pulled in as a dependency).

  fastcluster-ward
      Spectral embedding → UMAP → unconstrained Ward linkage via fastcluster's
      C++ nearest-neighbour-chain algorithm.  Unlike the two methods above,
      there is no connectivity constraint, so any two clusters may merge;
      this gives the purest Ward solution in the embedding space but does not
      explicitly enforce that only genomically connected genomes are combined.
      O(n² × k) in time (k = UMAP output dimension), but fastcluster's C++
      implementation is highly optimised and typically fast enough for n ≈ 143k
      with a compact embedding (k ≤ 15).
      Requires: fastcluster, umap-learn.

  hdbscan
      Spectral embedding → UMAP → HDBSCAN (hdbscan package).  HDBSCAN uses
      Borůvka's MST on mutual reachability distances (density-weighted), giving
      a single-linkage-style dendrogram that handles clusters of varying
      density.  Leaf order is extracted from the internal single-linkage tree
      (hdbscan.single_linkage_tree_._linkage → scipy leaves_list).
      O(n log n) with PyNNDescent.
      Requires: hdbscan, umap-learn.

  single
      Minimum spanning tree (MST) on the sparse distance graph (d = 1−Jaccard),
      converted to a single-linkage dendrogram via union-find.  O(E log V) —
      the fastest method — but prone to the chaining effect on sparse graphs.
      No embedding step; works directly from the AMG pairs.

Usage (run inside the sourmash conda env):

    conda run -n sourmash python3 scripts/plot_amg_hierarchical_heatmap.py

Key flags:
    --amg-pairwise   DIR     (default: data/GTDB/alphamaxgeom_pairwise_thr0001)
    --output         DIR     (default: data/GTDB/figures_thr0001)
    --ordering       spectral-ward|umap-ward|fastcluster-ward|hdbscan|single
    --spectral-k     INT     Laplacian eigenvectors (default: 20)
    --umap-components INT    UMAP output dimension (default: 15)
    --umap-neighbors  INT    UMAP n_neighbors (default: 15)
    --hdbscan-min-cluster INT (default: 5)
    --canvas-size    INT     pixels per axis per panel (default: 2000)
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
        "datashader is required.\n"
        "Run inside the sourmash conda env:\n"
        "  conda run -n sourmash python3 scripts/plot_amg_hierarchical_heatmap.py"
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
    r = rows.astype(np.int64)
    c = cols.astype(np.int64)
    return np.minimum(r, c) * N + np.maximum(r, c)


def load_npz(npz_path: Path, index_path: Path) -> dict:
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

    log.info("  Loaded %d pairs (%d genomes in index)", len(keys), N)
    return {
        "genome_ids":      genome_ids,
        "N":               N,
        "keys":            keys,
        "jaccard":         _get("jaccard"),
        "max_containment": _get("max_containment"),
    }


# ---------------------------------------------------------------------------
# Sparse graph helpers
# ---------------------------------------------------------------------------

def build_sparse_similarity(keys, vals, N):
    from scipy.sparse import csr_matrix
    rows = (keys // N).astype(np.int32)
    cols = (keys  % N).astype(np.int32)
    r    = np.concatenate([rows, cols])
    c    = np.concatenate([cols, rows])
    v    = np.tile(vals.astype(np.float32), 2)
    return csr_matrix((v, (r, c)), shape=(N, N), dtype=np.float32)


def build_sparse_distance(keys, vals, N):
    from scipy.sparse import csr_matrix
    rows = (keys // N).astype(np.int32)
    cols = (keys  % N).astype(np.int32)
    d    = np.maximum(1.0 - vals.astype(np.float64), 0.0)
    r    = np.concatenate([rows, cols])
    c    = np.concatenate([cols, rows])
    dv   = np.tile(d, 2)
    return csr_matrix((dv, (r, c)), shape=(N, N), dtype=np.float64)


# ---------------------------------------------------------------------------
# Shared embedding steps
# ---------------------------------------------------------------------------

def spectral_embedding(keys, vals, N, k: int) -> np.ndarray:
    """
    Compact spectral representation of the genome similarity graph.

    Returns an (N, k) float32 feature matrix: the top-k left singular vectors
    of the symmetrically-normalised adjacency matrix A_norm = D^{-1/2} A D^{-1/2},
    scaled by their singular values (X = U[:, :k] * s[:k]).

    Mathematical equivalence with the Laplacian approach
    ----------------------------------------------------
    The eigenvectors of A_norm and those of the normalised Laplacian
    L = I − A_norm are identical; eigenvalues are related by λ_L = 1 − λ_A.
    The top-k singular vectors of A_norm therefore carry the same information
    as the k smallest non-trivial eigenvectors of L.

    Why randomised SVD instead of LOBPCG on the Laplacian
    ------------------------------------------------------
    LOBPCG finds the SMALLEST eigenvalues of L.  For a graph with c connected
    components, c eigenvectors have eigenvalue exactly 0 (trivial).  When
    c ≫ k (which happens at GTDB scale with threshold=0.001 — thousands of
    isolated genomes), LOBPCG fills its entire budget with trivial vectors
    and returns nothing useful.

    Randomised SVD finds the LARGEST singular values of A_norm, which
    correspond to the k most informative non-trivial structures regardless of
    how many disconnected components exist.  Isolated genomes (zero rows in A)
    receive a zero feature vector and are naturally grouped together by Ward.

    Complexity: O(N × k × n_iter) — far faster than LOBPCG for large sparse
    graphs, and immune to trivial eigenvector pollution.
    """
    from sklearn.utils.extmath import randomized_svd
    from scipy.sparse import diags

    A   = build_sparse_similarity(keys, vals, N)
    deg = np.asarray(A.sum(axis=1)).flatten().astype(np.float64)
    deg[deg == 0] = 1.0                          # isolated nodes → identity row
    D_inv_sqrt = diags(1.0 / np.sqrt(deg))
    A_norm     = D_inv_sqrt @ A @ D_inv_sqrt     # symmetric normalised adjacency

    log.info("  Randomised SVD (k=%d, n=%d) ...", k, N)
    U, s, _ = randomized_svd(A_norm, n_components=k, n_iter=5, random_state=42)
    log.info("  Top singular values: %s", np.round(s[:min(5, k)], 4).tolist())
    return (U * s).astype(np.float32)


def umap_embedding(X: np.ndarray, n_components: int, n_neighbors: int) -> np.ndarray:
    """
    Non-linear UMAP refinement of the spectral embedding X.
    Returns an (N, n_components) float32 array.
    """
    import umap as _umap
    import warnings
    log.info("  UMAP: %d → %d dimensions (n_neighbors=%d) ...",
             X.shape[1], n_components, n_neighbors)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reducer = _umap.UMAP(
            n_components  = n_components,
            n_neighbors   = n_neighbors,
            min_dist      = 0.0,   # pack clusters tightly; better for Ward
            random_state  = 42,
            low_memory    = True,
        )
        emb = reducer.fit_transform(X)
    log.info("  UMAP done: output shape %s", emb.shape)
    return emb.astype(np.float32)


# ---------------------------------------------------------------------------
# Ordering method 1: spectral-ward
# ---------------------------------------------------------------------------

def _ward_from_features(X: np.ndarray, connectivity, N: int) -> np.ndarray:
    """
    sklearn Ward agglomerative clustering on feature matrix X with the given
    connectivity constraint.  Returns rank (int32, length N).
    """
    from sklearn.cluster import AgglomerativeClustering

    clf = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=0,
        connectivity=connectivity,
        linkage="ward",
        compute_full_tree=True,
    )
    clf.fit(X)

    # Iterative DFS on the merge tree → leaf order
    children = clf.children_
    root     = 2 * N - 2
    stack    = [root]
    leaves   = []
    while stack:
        node = stack.pop()
        if node < N:
            leaves.append(node)
        else:
            i = node - N
            stack.append(int(children[i, 1]))   # right child (pushed first)
            stack.append(int(children[i, 0]))   # left child (processed first)

    rank = np.empty(N, dtype=np.int32)
    for pos, leaf in enumerate(leaves):
        rank[leaf] = pos
    return rank


def spectral_ward_order(keys, vals, N, k: int) -> np.ndarray:
    """Spectral embedding (K eigenvectors) → Ward + sparse connectivity."""
    log.info("  Computing %d-D spectral embedding ...", k)
    X = spectral_embedding(keys, vals, N, k)
    A = build_sparse_similarity(keys, vals, N)
    log.info("  Fitting Ward (connectivity=sparse, n=%d) ...", N)
    return _ward_from_features(X, A, N)


# ---------------------------------------------------------------------------
# Ordering method 2: umap-ward
# ---------------------------------------------------------------------------

def umap_ward_order(keys, vals, N, k: int,
                    umap_components: int, umap_neighbors: int) -> np.ndarray:
    """Spectral → UMAP (non-linear) → Ward + sparse connectivity."""
    log.info("  Computing %d-D spectral embedding ...", k)
    X_spec = spectral_embedding(keys, vals, N, k)
    X_umap = umap_embedding(X_spec, umap_components, umap_neighbors)
    A      = build_sparse_similarity(keys, vals, N)
    log.info("  Fitting Ward (connectivity=sparse, n=%d) ...", N)
    return _ward_from_features(X_umap, A, N)


# ---------------------------------------------------------------------------
# Ordering method 3: fastcluster-ward (unconstrained)
# ---------------------------------------------------------------------------

def fastcluster_ward_order(keys, vals, N, k: int,
                           umap_components: int, umap_neighbors: int) -> np.ndarray:
    """
    Spectral → UMAP → unconstrained Ward via fastcluster's C++ nn-chain.
    No connectivity constraint: any two clusters may merge.
    O(n² × k) time, O(n × k) space.
    """
    import fastcluster
    from scipy.cluster.hierarchy import leaves_list

    log.info("  Computing %d-D spectral embedding ...", k)
    X_spec = spectral_embedding(keys, vals, N, k)
    X_umap = umap_embedding(X_spec, umap_components, umap_neighbors)

    log.info("  fastcluster Ward linkage (n=%d, k=%d, unconstrained) ...",
             N, X_umap.shape[1])
    Z     = fastcluster.linkage_vector(X_umap, method="ward")
    order = leaves_list(Z)

    rank = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


# ---------------------------------------------------------------------------
# Ordering method 4: hdbscan
# ---------------------------------------------------------------------------

def hdbscan_order(keys, vals, N, k: int,
                  umap_components: int, umap_neighbors: int,
                  min_cluster_size: int) -> np.ndarray:
    """
    Spectral → UMAP → HDBSCAN (Borůvka MST on mutual reachability distances).
    Leaf order from the internal single-linkage tree.
    """
    import hdbscan as _hdbscan
    from scipy.cluster.hierarchy import leaves_list

    log.info("  Computing %d-D spectral embedding ...", k)
    X_spec = spectral_embedding(keys, vals, N, k)
    X_umap = umap_embedding(X_spec, umap_components, umap_neighbors)

    log.info("  HDBSCAN (min_cluster_size=%d, n=%d) ...", min_cluster_size, N)
    clf = _hdbscan.HDBSCAN(
        min_cluster_size = min_cluster_size,
        min_samples      = 1,
        gen_min_span_tree= True,
        core_dist_n_jobs = -1,
    )
    clf.fit(X_umap)

    Z     = clf.single_linkage_tree_._linkage   # scipy-format linkage matrix
    order = leaves_list(Z)

    rank = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


# ---------------------------------------------------------------------------
# Ordering method 5: MST single linkage
# ---------------------------------------------------------------------------

def _mst_to_linkage(mst_csr, n: int) -> np.ndarray:
    """
    Convert a scipy CSR minimum spanning tree to a valid scipy linkage matrix.

    scipy requires that node indices in Z[:,0:2] are proper internal-node IDs:
    leaves are 0..n-1, and each merge i creates a new node with ID n+i.
    We maintain a node_id[] array that maps the current union-find root of each
    component to its linkage-matrix node ID, updated after every merge.

    Disconnected components (graph not fully connected) are joined at
    max_edge_weight + ε so every genome appears in the final dendrogram.
    """
    mst_coo = mst_csr.tocoo()
    idx     = np.argsort(mst_coo.data)
    w_all   = mst_coo.data[idx].tolist()
    i_all   = mst_coo.row[idx].tolist()
    j_all   = mst_coo.col[idx].tolist()

    # Union-Find (path-halving)
    parent  = list(range(n))
    sz      = [1] * n
    # node_id[root] = the scipy linkage node-ID for the cluster rooted there.
    # Initially every genome i is a leaf node with ID i.
    node_id = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    Z       = []
    next_id = n   # next internal node ID to assign (n, n+1, …, 2n-2)

    def merge(ri, rj, w):
        nonlocal next_id
        Z.append([float(node_id[ri]), float(node_id[rj]),
                  float(w), float(sz[ri] + sz[rj])])
        sz[rj]      += sz[ri]
        parent[ri]   = rj
        node_id[rj]  = next_id
        next_id     += 1

    for w, i, j in zip(w_all, i_all, j_all):
        ri, rj = find(i), find(j)
        if ri != rj:
            merge(ri, rj, w)

    # Merge remaining disconnected components
    if len(Z) < n - 1:
        max_w = (Z[-1][2] if Z else 0.0) + 1.0
        roots = []
        seen  = set()
        for k in range(n):
            r = find(k)
            if r not in seen:
                seen.add(r)
                roots.append(r)
        for step in range(len(roots) - 1):
            ri = find(roots[step])
            rj = find(roots[step + 1])
            if ri != rj:
                merge(ri, rj, max_w + step * 1e-9)

    return np.array(Z, dtype=np.float64)


def single_linkage_order(keys, vals, N) -> np.ndarray:
    """MST single linkage on the sparse distance graph.  O(E log V)."""
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.cluster.hierarchy import leaves_list

    log.info("  Building sparse distance graph (%d edges) ...", len(keys))
    D = build_sparse_distance(keys, vals, N)
    log.info("  Computing MST ...")
    mst = minimum_spanning_tree(D)
    log.info("  MST: %d edges", mst.nnz)
    Z     = _mst_to_linkage(mst, N)
    order = leaves_list(Z)
    rank  = np.empty(N, dtype=np.int32)
    rank[order] = np.arange(N, dtype=np.int32)
    return rank


# ---------------------------------------------------------------------------
# Datashader rasterisation
# ---------------------------------------------------------------------------

def make_similarity_dataframe(keys, vals, rank, N) -> pd.DataFrame:
    rows = (keys // N).astype(np.int32)
    cols = (keys  % N).astype(np.int32)
    x    = np.concatenate([rank[rows], rank[cols]]).astype(np.float32)
    y    = np.concatenate([rank[cols], rank[rows]]).astype(np.float32)
    v    = np.tile(vals.astype(np.float32), 2)
    return pd.DataFrame({"x": x, "y": y, "val": v})


def rasterise(df, col, N, canvas_size, vmax, lut_rgba) -> np.ndarray:
    cvs  = ds.Canvas(
        plot_width=canvas_size, plot_height=canvas_size,
        x_range=(0, N - 1), y_range=(0, N - 1),
    )
    agg  = cvs.points(df, "x", "y", ds.mean(col))
    vals = agg.values.astype(np.float32)
    idx  = np.clip(
        np.where(np.isnan(vals), 0, vals / vmax * 255).astype(np.int32), 0, 255)
    rgba = lut_rgba[idx].copy()
    rgba[np.isnan(vals)] = (np.array(mcolors.to_rgba(MISSING_COLOR)) * 255).astype(np.uint8)
    return rgba.astype(np.uint8)


# ---------------------------------------------------------------------------
# Figure composition
# ---------------------------------------------------------------------------

_METHOD_LABELS = {
    "spectral-ward":     "spectral-Ward (connectivity-constrained)",
    "umap-ward":         "UMAP → Ward (connectivity-constrained)",
    "fastcluster-ward":  "UMAP → fastcluster Ward (unconstrained)",
    "hdbscan":           "UMAP → HDBSCAN single-linkage",
    "single":            "MST single-linkage",
}


def plot_heatmaps(pairwise, rank, canvas_size, ordering, out_dir):
    N    = pairwise["N"]
    keys = pairwise["keys"]
    cmap = plt.cm.viridis
    lut  = (cmap(np.linspace(0, 1, 256)) * 255).astype(np.uint8)

    panels = []
    for metric, label in [("jaccard", "Jaccard"), ("max_containment", "Max containment")]:
        v = pairwise[metric]
        if v is None:
            log.warning("Metric '%s' not found — skipping.", metric)
            continue
        vmax = float(np.nanpercentile(v, 90)) or 1.0
        log.info("  Building %s DataFrame (vmax=%.4f) ...", label, vmax)
        df = make_similarity_dataframe(keys, v, rank, N)
        panels.append({"label": label, "df": df, "vmax": vmax})

    if not panels:
        log.error("No metrics to plot.")
        return

    n_panels = len(panels)
    panel_w  = 4.0
    cbar_h   = 0.35
    fig = plt.figure(figsize=(panel_w * n_panels + 0.5, panel_w + cbar_h + 0.7))
    gs  = fig.add_gridspec(
        2, n_panels,
        height_ratios=[panel_w, cbar_h],
        hspace=0.08, wspace=0.12,
        left=0.04, right=0.97,
        top=0.93,  bottom=0.02,
    )

    method_label = _METHOD_LABELS.get(ordering, ordering)

    for col_idx, p in enumerate(panels):
        log.info("  Rasterising %s (%d×%d) ...", p["label"], canvas_size, canvas_size)
        rgba = rasterise(p["df"], "val", N, canvas_size, p["vmax"], lut)

        ax = fig.add_subplot(gs[0, col_idx])
        ax.imshow(rgba, aspect="auto", interpolation="nearest")
        ax.set_title(
            f"AlphaMaxGeomHash  —  {p['label']}\n"
            f"{method_label}  |  "
            f"{N:,} genomes  |  {len(keys):,} pairs above threshold",
            fontsize=9, pad=4,
        )
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)

        sm  = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(0, p["vmax"]))
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

    stem = f"amg_hierarchical_heatmap_{ordering.replace('-', '_')}"
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
        description="AMG hierarchical-clustering heatmaps (Jaccard + max_containment).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--amg-pairwise",
                   default=str(base / "data/GTDB/alphamaxgeom_pairwise_thr0001"),
                   help="AlphaMaxGeomHash pairwise directory")
    p.add_argument("--output",
                   default=str(base / "data/GTDB/figures_thr0001"),
                   help="Output directory")
    p.add_argument("--ordering",
                   choices=["spectral-ward", "umap-ward",
                             "fastcluster-ward", "hdbscan", "single"],
                   default="spectral-ward")
    p.add_argument("--spectral-k",      type=int, default=20,
                   help="Laplacian eigenvectors for spectral embedding")
    p.add_argument("--umap-components", type=int, default=15,
                   help="UMAP output dimensions (umap-ward / fastcluster-ward / hdbscan)")
    p.add_argument("--umap-neighbors",  type=int, default=15,
                   help="UMAP n_neighbors")
    p.add_argument("--hdbscan-min-cluster", type=int, default=5,
                   help="HDBSCAN min_cluster_size")
    p.add_argument("--canvas-size",     type=int, default=2000,
                   help="Pixels per axis per panel")
    return p.parse_args()


def main():
    args    = parse_args()
    amg_dir = Path(args.amg_pairwise)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading AlphaMaxGeomHash pairwise from %s ...", amg_dir)
    pw = load_npz(amg_dir / "pairwise_results.npz", amg_dir / "genome_index.json")

    log.info("Computing genome ordering: %s (n=%d) ...", args.ordering, pw["N"])

    if args.ordering == "spectral-ward":
        rank = spectral_ward_order(pw["keys"], pw["jaccard"], pw["N"],
                                   args.spectral_k)
    elif args.ordering == "umap-ward":
        rank = umap_ward_order(pw["keys"], pw["jaccard"], pw["N"],
                               args.spectral_k, args.umap_components,
                               args.umap_neighbors)
    elif args.ordering == "fastcluster-ward":
        rank = fastcluster_ward_order(pw["keys"], pw["jaccard"], pw["N"],
                                      args.spectral_k, args.umap_components,
                                      args.umap_neighbors)
    elif args.ordering == "hdbscan":
        rank = hdbscan_order(pw["keys"], pw["jaccard"], pw["N"],
                             args.spectral_k, args.umap_components,
                             args.umap_neighbors, args.hdbscan_min_cluster)
    else:  # single
        rank = single_linkage_order(pw["keys"], pw["jaccard"], pw["N"])

    log.info("Ordering done.  Plotting ...")
    plot_heatmaps(
        pairwise    = pw,
        rank        = rank,
        canvas_size = args.canvas_size,
        ordering    = args.ordering,
        out_dir     = out_dir,
    )
    log.info("Done.")


if __name__ == "__main__":
    main()
