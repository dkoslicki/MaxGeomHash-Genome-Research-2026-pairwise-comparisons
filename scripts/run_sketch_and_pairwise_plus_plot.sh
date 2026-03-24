#!/bin/bash
BASE="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
LOG="$BASE/scripts/kmer_sketch_full_run.log"
for s in 05_alphamaxgeom_sketch.sh 03_bottomk_sketch.sh 04_fracminhash_sketch.sh 08_alphamaxgeom_pairwise.sh 06_bottomk_pairwise.sh 07_fracminhash_pairwise.sh; do
    echo ""
    echo "=== START $s $(date) ==="
    bash "$BASE/scripts/$s" && echo "=== END $s $(date) ===" || { echo "FAILED: $s"; exit 1; }
done
bash "$BASE/scripts/10_run_plots.sh"
