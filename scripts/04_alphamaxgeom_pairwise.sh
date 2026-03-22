#!/usr/bin/env bash
# ==============================================================================
# 04_alphamaxgeom_pairwise.sh
#
# Wrapper: compute AlphaMaxGeomHash pairwise similarities for all genome pairs
# that passed the FracMinHash 0.01 threshold (from gtdb_pairwise_containment.csv).
#
# Only pairs where both genomes have been sketched are evaluated.
# Run 03_alphamaxgeom_sketch.sh first (optionally with TEST_N to sketch a
# subset; this script will then compute only the pairs within that subset).
#
# Full run:    bash 04_alphamaxgeom_pairwise.sh
# After test:  bash 04_alphamaxgeom_pairwise.sh   (same command — adapts automatically)
# ==============================================================================
set -euo pipefail

BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"

SKETCH_DIR="${GTDB_DIR}/alphamaxgeom_sketches"
CANDIDATES="${GTDB_DIR}/gtdb_pairwise_containment.csv"
OUTPUT_DIR="${GTDB_DIR}/alphamaxgeom_pairwise"
CORES=192

echo "Starting AlphaMaxGeomHash pairwise computation ..."
echo "  Sketch dir  : ${SKETCH_DIR}"
echo "  Candidates  : ${CANDIDATES}"
echo "  Output dir  : ${OUTPUT_DIR}"
echo "  Cores       : ${CORES}"
echo ""

python3 "${BASE_DIR}/scripts/04_alphamaxgeom_pairwise.py" \
    --sketch-dir "${SKETCH_DIR}" \
    --candidates "${CANDIDATES}" \
    --output     "${OUTPUT_DIR}" \
    --cores      "${CORES}"
