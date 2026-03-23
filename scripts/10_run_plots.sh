#!/usr/bin/env bash
# ==============================================================================
# 10_run_plots.sh
#
# Generate all publication figures comparing kmer-sketch methods
# (AlphaMaxGeomHash, MinHash, FracMinHash) against KMC exact ground truth,
# across the 143,614 GTDB representative genomes.
#
# Figures produced (in data/GTDB/figures/):
#   heatmap_jaccard_relative_error.{pdf,png}          — per-pair relative-error
#   heatmap_max_containment_relative_error.{pdf,png}     heatmaps (N-panel)
#   resources_time.{pdf,png}                           — computation time bars
#   resources_disk.{pdf,png}                           — disk-space bars
#   resources_tradeoff.{pdf,png}                       — accuracy vs. resources
#
# Prerequisites:
#   - 03_bottomk_sketch.sh         (BottomK sketches)
#   - 04_fracminhash_sketch.sh     (FracMinHash kmer-sketch sketches)
#   - 05_alphamaxgeom_sketch.sh    (AlphaMaxGeomHash sketches)
#   - 06_bottomk_pairwise.sh       (BottomK pairwise)
#   - 07_fracminhash_pairwise.sh   (FracMinHash pairwise)
#   - 08_alphamaxgeom_pairwise.sh  (AlphaMaxGeomHash pairwise)
#   - 09_sanity_check.py           (sanity_check_summary.json)
#
# Methods that have not yet been run are gracefully skipped — the script
# automatically omits methods whose output directories are absent.
#
# Usage:
#   bash scripts/10_run_plots.sh              # default: 500-genome heatmap
#   N_GENOMES=300 bash scripts/10_run_plots.sh  # smaller heatmap subset
#   SHOW_REF=1    bash scripts/10_run_plots.sh  # add KMC reference panel
#
# To also include Sourmash FracMinHash:
#   INCLUDE_SOURMASH=1 bash scripts/10_run_plots.sh
#
# ==============================================================================

set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --------------- Configurable parameters --------------------------------------
FIGURE_DIR="${BASE}/data/GTDB/figures"
N_GENOMES="${N_GENOMES:-500}"          # genomes to include in heatmaps
SHOW_REF="${SHOW_REF:-0}"             # 1 = add KMC reference panel to heatmap
INCLUDE_SOURMASH="${INCLUDE_SOURMASH:-0}"  # 1 = add Sourmash FracMinHash

KMC_PAIRWISE="${BASE}/data/GTDB/kmc_pairwise"
AMG_PAIRWISE="${BASE}/data/GTDB/alphamaxgeom_pairwise"
BK_PAIRWISE="${BASE}/data/GTDB/bottomk_pairwise"
FMH_KS_PAIRWISE="${BASE}/data/GTDB/fracminhash_pairwise"
SOURMASH_CSV="${BASE}/data/GTDB/gtdb_pairwise_containment.csv"
SANITY_JSON="${AMG_PAIRWISE}/sanity_check/sanity_check_summary.json"
# ------------------------------------------------------------------------------

echo "========================================================"
echo "  10_run_plots.sh — Publication figure generation"
echo "  Start: $(date)"
echo "  Output: ${FIGURE_DIR}"
echo "  Heatmap genomes: ${N_GENOMES}"
echo "========================================================"

mkdir -p "${FIGURE_DIR}"

# Build optional flags
SHOW_REF_FLAG=""
if [[ "${SHOW_REF}" == "1" ]]; then
    SHOW_REF_FLAG="--show-reference"
fi

SOURMASH_FLAGS=""
if [[ "${INCLUDE_SOURMASH}" == "1" ]]; then
    SOURMASH_FLAGS="--include-sourmash --sourmash-csv ${SOURMASH_CSV}"
fi

# Build heatmap method flags (only pass directories that exist)
HEATMAP_METHOD_FLAGS=""
[[ -d "${AMG_PAIRWISE}" ]]    && HEATMAP_METHOD_FLAGS+=" --amg-pairwise ${AMG_PAIRWISE}"
[[ -d "${BK_PAIRWISE}" ]]     && HEATMAP_METHOD_FLAGS+=" --bottomk-pairwise ${BK_PAIRWISE}"
[[ -d "${FMH_KS_PAIRWISE}" ]] && HEATMAP_METHOD_FLAGS+=" --fracminhash-pairwise ${FMH_KS_PAIRWISE}"
if [[ "${INCLUDE_SOURMASH}" == "1" ]]; then
    HEATMAP_METHOD_FLAGS+=" --include-sourmash --sourmash-csv ${SOURMASH_CSV}"
fi

# ---- Step 1: Heatmaps (Jaccard and Max-containment) -------------------------
for METRIC in jaccard max_containment; do
    echo ""
    echo "--- Heatmap: ${METRIC} ---"
    echo "Start: $(date)"
    python3 "${BASE}/scripts/10_plot_heatmaps.py" \
        --kmc-pairwise  "${KMC_PAIRWISE}" \
        --output        "${FIGURE_DIR}" \
        --n-genomes     "${N_GENOMES}" \
        --metric        "${METRIC}" \
        ${SHOW_REF_FLAG} \
        ${HEATMAP_METHOD_FLAGS}
    echo "End: $(date)"
done

# ---- Step 2: Resource comparison figures ------------------------------------
echo ""
echo "--- Resource figures (time, disk, trade-off) ---"
echo "Start: $(date)"
python3 "${BASE}/scripts/10_plot_resources.py" \
    --base-dir    "${BASE}" \
    --sanity-json "${SANITY_JSON}" \
    --output      "${FIGURE_DIR}" \
    ${SOURMASH_FLAGS}
echo "End: $(date)"

echo ""
echo "========================================================"
echo "  Figures written to: ${FIGURE_DIR}"
ls -lh "${FIGURE_DIR}"/*.pdf "${FIGURE_DIR}"/*.png 2>/dev/null || true
echo ""
echo "  Pipeline COMPLETE"
echo "  End: $(date)"
echo "========================================================"
