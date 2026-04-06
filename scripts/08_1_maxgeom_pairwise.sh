#!/usr/bin/env bash
# ==============================================================================
# 08_1_maxgeom_pairwise.sh
#
# Wrapper: compute MaxGeomHash pairwise similarities for all genome pairs
# that passed the FracMinHash all-vs-all candidate filter (max_containment ≥ 0.001,
# from gtdb_pairwise_containment_thr0001.csv; see compute_fracminhash_candidates.sh).
#
# Only pairs where both genomes have been sketched are evaluated.
# Run 05_1_maxgeom_sketch.sh first (optionally with TEST_N to sketch a
# subset; this script will then compute only the pairs within that subset).
#
# Full run:    bash 08_1_maxgeom_pairwise.sh
# After test:  bash 08_1_maxgeom_pairwise.sh   (same command — adapts automatically)
# ==============================================================================
set -euo pipefail

BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"

SKETCH_DIR="${GTDB_DIR}/maxgeom_sketches"
CANDIDATES="${GTDB_DIR}/gtdb_pairwise_containment_thr0001.csv"
OUTPUT_DIR="${GTDB_DIR}/maxgeom_pairwise_thr0001"

CONFIG="${BASE_DIR}/config.json"
CORES=$(python3 -c "import json; print(json.load(open('${CONFIG}'))['pairwise_cores'])")

echo "Starting MaxGeomHash pairwise computation ..."
echo "  Sketch dir  : ${SKETCH_DIR}"
echo "  Candidates  : ${CANDIDATES}"
echo "  Output dir  : ${OUTPUT_DIR}"
echo "  Cores       : ${CORES}"
echo ""

python3 "${BASE_DIR}/scripts/08_1_maxgeom_pairwise.py" \
    --sketch-dir "${SKETCH_DIR}" \
    --candidates "${CANDIDATES}" \
    --output     "${OUTPUT_DIR}" \
    --cores      "${CORES}"
