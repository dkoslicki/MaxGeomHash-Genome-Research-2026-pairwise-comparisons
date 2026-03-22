#!/usr/bin/env bash
# ==============================================================================
# 06_run_plots.sh
#
# Generate all publication figures comparing AlphaMaxGeomHash, FracMinHash,
# and KMC (exact ground truth) across the 143,614 GTDB representative genomes.
#
# Figures produced (in data/GTDB/figures/):
#   heatmap_jaccard_relative_error.{pdf,png}       — per-pair relative-error
#   heatmap_max_containment_relative_error.{pdf,png}  heatmaps (2-panel)
#   resources_time.{pdf,png}                        — computation time bars
#   resources_disk.{pdf,png}                        — disk-space bars
#   resources_tradeoff.{pdf,png}                    — accuracy vs. resources
#
# Prerequisites:
#   - 03_alphamaxgeom_sketch.sh   (AMG sketches)
#   - 04_alphamaxgeom_pairwise.sh (AMG pairwise similarities)
#   - 05_sanity_check.py          (sanity_check_summary.json)
#
# Usage:
#   bash scripts/06_run_plots.sh              # all metrics, 500-genome heatmap
#   N_GENOMES=300 bash scripts/06_run_plots.sh  # smaller heatmap subset
#   SHOW_REF=1    bash scripts/06_run_plots.sh  # include KMC reference panel
#
# ==============================================================================

set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --------------- Configurable parameters --------------------------------------
FIGURE_DIR="${BASE}/data/GTDB/figures"
N_GENOMES="${N_GENOMES:-500}"          # genomes to include in heatmaps
SHOW_REF="${SHOW_REF:-0}"             # 1 = add KMC reference panel to heatmap

AMG_PAIRWISE="${BASE}/data/GTDB/alphamaxgeom_pairwise"
FRACMINHASH="${BASE}/data/GTDB/gtdb_pairwise_containment.csv"
KMC_PAIRWISE="${BASE}/data/GTDB/kmc_pairwise"
SANITY_JSON="${AMG_PAIRWISE}/sanity_check/sanity_check_summary.json"
# ------------------------------------------------------------------------------

echo "========================================================"
echo "  06_run_plots.sh — Publication figure generation"
echo "  Start: $(date)"
echo "  Output: ${FIGURE_DIR}"
echo "  Heatmap genomes: ${N_GENOMES}"
echo "========================================================"

mkdir -p "${FIGURE_DIR}"

# Optional --show-reference flag for the heatmap script
SHOW_REF_FLAG=""
if [[ "${SHOW_REF}" == "1" ]]; then
    SHOW_REF_FLAG="--show-reference"
fi

# ---- Step 1: Heatmaps (Jaccard and Max-containment) -------------------------
for METRIC in jaccard max_containment; do
    echo ""
    echo "--- Heatmap: ${METRIC} ---"
    echo "Start: $(date)"
    python3 "${BASE}/scripts/06_plot_heatmaps.py" \
        --amg-pairwise  "${AMG_PAIRWISE}" \
        --fracminhash   "${FRACMINHASH}" \
        --kmc-pairwise  "${KMC_PAIRWISE}" \
        --output        "${FIGURE_DIR}" \
        --n-genomes     "${N_GENOMES}" \
        --metric        "${METRIC}" \
        ${SHOW_REF_FLAG}
    echo "End: $(date)"
done

# ---- Step 2: Resource comparison figures ------------------------------------
echo ""
echo "--- Resource figures (time, disk, trade-off) ---"
echo "Start: $(date)"
python3 "${BASE}/scripts/06_plot_resources.py" \
    --base-dir    "${BASE}" \
    --sanity-json "${SANITY_JSON}" \
    --output      "${FIGURE_DIR}"
echo "End: $(date)"

echo ""
echo "========================================================"
echo "  Figures written to: ${FIGURE_DIR}"
ls -lh "${FIGURE_DIR}"/*.pdf "${FIGURE_DIR}"/*.png 2>/dev/null || true
echo ""
echo "  Pipeline COMPLETE"
echo "  End: $(date)"
echo "========================================================"
