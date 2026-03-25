#!/usr/bin/env bash
# ==============================================================================
# run_sketch_and_pairwise_plus_plot.sh
#
# Master pipeline: run ALL steps from sourmash candidate generation through
# kmer-sketch pairwise comparisons to figure generation.
#
# Prerequisite (already done, one-time):
#   bash scripts/00_decompress_genomes.sh   # decompress GTDB genomes
#   bash scripts/01_kmc_count.sh            # KMC k-mer counting
#
# This script runs (in order):
#   1. make_sourmash_sketches.sh            — sourmash FracMinHash, scaled=100
#   2. compute_sourmash_pairwise.sh         — all-vs-all, threshold=0.001 → candidates CSV
#   3. 02_kmc_pairwise.sh                   — KMC exact pairwise on candidates
#   4. 03_bottomk_sketch.sh                 — BottomK sketches (kmer-sketch)
#   5. 04_fracminhash_sketch.sh             — FracMinHash sketches (kmer-sketch)
#   6. 05_alphamaxgeom_sketch.sh            — AlphaMaxGeomHash sketches (kmer-sketch)
#   7. 06_bottomk_pairwise.sh              — BottomK pairwise on candidates
#   8. 07_fracminhash_pairwise.sh          — FracMinHash pairwise on candidates
#   9. 08_alphamaxgeom_pairwise.sh         — AlphaMaxGeomHash pairwise on candidates
#  10. 10_run_plots.sh                     — publication figures (full + subset heatmaps)
#
# Note: 09_sanity_check.py is omitted here because it is slow; run it manually
# when needed.
#
# Steps 3-9 all depend on the candidates CSV from step 2.
# Steps 4-6 (sketching) can run in parallel with each other; steps 7-9
# each depend on their corresponding sketch step.
# Steps 7-9 can all run in parallel once 4-6 are done.
#
# Usage:
#   bash scripts/run_sketch_and_pairwise_plus_plot.sh
#   # Logs are written to scripts/full_pipeline.log
# ==============================================================================

set -euo pipefail

BASE="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
LOG="$BASE/scripts/full_pipeline.log"
GTDB_DIR="$BASE/data/GTDB"

run_step() {
    local label="$1"
    local script="$2"
    echo ""
    echo "=== START ${label} $(date) ===" | tee -a "${LOG}"
    bash "${script}" 2>&1 | tee -a "${LOG}"
    echo "=== END   ${label} $(date) ===" | tee -a "${LOG}"
}

echo "========================================================"  | tee "${LOG}"
echo "  Full Pipeline — scaled=100, threshold=0.001"           | tee -a "${LOG}"
echo "  Start: $(date)"                                         | tee -a "${LOG}"
echo "========================================================"  | tee -a "${LOG}"

# Step 1: Sourmash sketching (scaled=100)
run_step "make_sourmash_sketches" "$BASE/scripts/make_sourmash_sketches.sh"

# Step 2: Sourmash pairwise (threshold=0.001) → gtdb_pairwise_containment_thr0001.csv
run_step "compute_sourmash_pairwise" "$BASE/scripts/compute_sourmash_pairwise.sh"

# Step 3: KMC exact pairwise on new candidate set
run_step "02_kmc_pairwise" "$BASE/scripts/02_kmc_pairwise.sh"

# Steps 4-6: kmer-sketch sketching (can run in parallel, but run serially here
# to avoid saturating I/O; adjust if on a multi-node setup)
run_step "03_bottomk_sketch"     "$BASE/scripts/03_bottomk_sketch.sh"
run_step "04_fracminhash_sketch" "$BASE/scripts/04_fracminhash_sketch.sh"
run_step "05_alphamaxgeom_sketch" "$BASE/scripts/05_alphamaxgeom_sketch.sh"

# Steps 7-9: kmer-sketch pairwise on new candidate set
run_step "06_bottomk_pairwise"     "$BASE/scripts/06_bottomk_pairwise.sh"
run_step "07_fracminhash_pairwise" "$BASE/scripts/07_fracminhash_pairwise.sh"
run_step "08_alphamaxgeom_pairwise" "$BASE/scripts/08_alphamaxgeom_pairwise.sh"

# Step 10: Figures — full Datashader heatmaps (11_plot_heatmaps_full.py) AND
# dense subset heatmaps (10_plot_heatmaps.py, SUBSET_HEATMAP=1)
SUBSET_HEATMAP=1 run_step "10_run_plots" "$BASE/scripts/10_run_plots.sh"

echo "" | tee -a "${LOG}"
echo "========================================================" | tee -a "${LOG}"
echo "  Pipeline COMPLETE"                                       | tee -a "${LOG}"
echo "  End: $(date)"                                           | tee -a "${LOG}"
echo "  Log: ${LOG}"                                            | tee -a "${LOG}"
echo "========================================================" | tee -a "${LOG}"
