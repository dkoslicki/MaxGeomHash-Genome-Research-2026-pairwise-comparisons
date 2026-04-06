#!/usr/bin/env bash
# ==============================================================================
# 10_run_plots.sh
#
# Generate all publication figures comparing kmer-sketch methods
# (AlphaMaxGeomHash, MinHash, FracMinHash) against KMC exact ground truth,
# across the 143,614 GTDB representative genomes.
#
# Figures produced (in data/GTDB/figures/):
#   heatmap_jaccard_relative_error_full_{ordering}.{pdf,png}
#   heatmap_max_containment_relative_error_full_{ordering}.{pdf,png}
#   accuracy_l1_jaccard.{pdf,png}
#   accuracy_l1_max_containment.{pdf,png}
#   resources_time.{pdf,png}
#   resources_disk.{pdf,png}
#   resources_tradeoff.{pdf,png}
#
# By default, heatmaps cover ALL ~143k genomes via Datashader (11_plot_heatmaps_full.py).
# Set SUBSET_HEATMAP=1 to additionally produce dense heatmaps for a genome subset
# (10_plot_heatmaps.py).
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
# Usage:
#   bash scripts/10_run_plots.sh                   # full heatmaps (default)
#   SUBSET_HEATMAP=1 bash scripts/10_run_plots.sh  # also produce N-genome subset heatmaps
#   N_GENOMES=300    bash scripts/10_run_plots.sh  # (with SUBSET_HEATMAP=1) subset size
#   ORDERING=rcm     bash scripts/10_run_plots.sh  # genome ordering: degree|rcm|spectral
#   SHOW_REF=1       bash scripts/10_run_plots.sh  # add KMC reference panel
#   INCLUDE_SOURMASH=1 bash scripts/10_run_plots.sh
#
# ==============================================================================

set -euo pipefail

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --------------- Configurable parameters --------------------------------------
FIGURE_DIR="${BASE}/data/GTDB/figures_thr0001"
SUBSET_HEATMAP="${SUBSET_HEATMAP:-0}"      # 1 = also run dense N-genome heatmap
N_GENOMES="${N_GENOMES:-500}"              # genome count for subset heatmap
ORDERING="${ORDERING:-spectral}"           # degree | rcm | spectral
SHOW_REF="${SHOW_REF:-0}"                  # 1 = add KMC reference panel
INCLUDE_SOURMASH="${INCLUDE_SOURMASH:-0}"  # 1 = add Sourmash FracMinHash
CONDA_ENV="${CONDA_ENV:-sourmash}"         # conda env that has datashader

KMC_PAIRWISE="${BASE}/data/GTDB/kmc_pairwise_thr0001"
MG_PAIRWISE="${BASE}/data/GTDB/maxgeom_pairwise_thr0001"
AMG_PAIRWISE="${BASE}/data/GTDB/alphamaxgeom_pairwise_thr0001"
BK_PAIRWISE="${BASE}/data/GTDB/bottomk_pairwise_thr0001"
FMH_KS_PAIRWISE="${BASE}/data/GTDB/fracminhash_pairwise_thr0001"
SOURMASH_CSV="${BASE}/data/GTDB/gtdb_pairwise_containment_thr0001.csv"
SANITY_JSON="${AMG_PAIRWISE}/sanity_check/sanity_check_summary.json"
# ------------------------------------------------------------------------------

echo "========================================================"
echo "  10_run_plots.sh — Publication figure generation"
echo "  Start: $(date)"
echo "  Output: ${FIGURE_DIR}"
echo "  Full heatmap ordering: ${ORDERING}"
if [[ "${SUBSET_HEATMAP}" == "1" ]]; then
    echo "  Subset heatmap: enabled (N=${N_GENOMES})"
fi
echo "========================================================"

mkdir -p "${FIGURE_DIR}"

# Build optional flags
SHOW_REF_FLAG=""
[[ "${SHOW_REF}" == "1" ]] && SHOW_REF_FLAG="--show-reference"

SOURMASH_FLAGS=""
[[ "${INCLUDE_SOURMASH}" == "1" ]] && SOURMASH_FLAGS="--include-sourmash --sourmash-csv ${SOURMASH_CSV}"

# Build method flags (only pass directories that exist)
METHOD_FLAGS=""
[[ -d "${BK_PAIRWISE}" ]]     && METHOD_FLAGS+=" --bottomk-pairwise ${BK_PAIRWISE}"
[[ -d "${MG_PAIRWISE}" ]]     && METHOD_FLAGS+=" --mg-pairwise ${MG_PAIRWISE}"
[[ -d "${AMG_PAIRWISE}" ]]    && METHOD_FLAGS+=" --amg-pairwise ${AMG_PAIRWISE}"
[[ -d "${FMH_KS_PAIRWISE}" ]] && METHOD_FLAGS+=" --fracminhash-pairwise ${FMH_KS_PAIRWISE}"

# ---- Step 1: Full heatmaps via Datashader (all ~143k genomes) ----------------
for METRIC in jaccard max_containment; do
    echo ""
    echo "--- Full heatmap (${METRIC}, ordering=${ORDERING}) ---"
    echo "Start: $(date)"
    conda run -n "${CONDA_ENV}" python3 "${BASE}/scripts/11_plot_heatmaps_full.py" \
        --kmc-pairwise "${KMC_PAIRWISE}" \
        --output       "${FIGURE_DIR}" \
        --metric       "${METRIC}" \
        --ordering     "${ORDERING}" \
        ${SHOW_REF_FLAG} \
        ${METHOD_FLAGS}
    echo "End: $(date)"
done

# ---- Step 2 (optional): Dense subset heatmaps --------------------------------
if [[ "${SUBSET_HEATMAP}" == "1" ]]; then
    for METRIC in jaccard max_containment; do
        echo ""
        echo "--- Subset heatmap (${METRIC}, N=${N_GENOMES}) ---"
        echo "Start: $(date)"
        python3 "${BASE}/scripts/10_plot_heatmaps.py" \
            --kmc-pairwise "${KMC_PAIRWISE}" \
            --output       "${FIGURE_DIR}" \
            --n-genomes    "${N_GENOMES}" \
            --metric       "${METRIC}" \
            --l1-mode      subset \
            ${SHOW_REF_FLAG} \
            ${METHOD_FLAGS}
        echo "End: $(date)"
    done
fi

# ---- Step 3: Resource comparison figures -------------------------------------
echo ""
echo "--- Resource figures (time, disk, accuracy) ---"
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
