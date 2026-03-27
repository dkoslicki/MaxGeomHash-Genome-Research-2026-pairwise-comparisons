#!/usr/bin/env bash
# ==============================================================================
# 02_kmc_pairwise.sh
#
# Compute exact pairwise KMC k-mer similarities for all GTDB genome pairs that
# passed the FracMinHash candidate filter (max_containment >= 0.001).
#
# Prerequisites:
#   bash scripts/01_kmc_count.sh                    (KMC k-mer databases)
#   bash scripts/compute_fracminhash_candidates.sh  (candidates CSV)
# ==============================================================================
set -euo pipefail

BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"

python3 "${BASE_DIR}/scripts/02_kmc_pairwise.py" \
    --metadata   "${GTDB_DIR}/kmc_metadata.csv" \
    --candidates "${GTDB_DIR}/gtdb_pairwise_containment_thr0001.csv" \
    --output     "${GTDB_DIR}/kmc_pairwise_thr0001" \
    --threshold  0.001 \
    --cores      192 \
    --chunk-size 1
