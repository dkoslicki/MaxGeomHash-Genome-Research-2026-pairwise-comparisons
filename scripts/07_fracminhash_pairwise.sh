#!/usr/bin/env bash
# ==============================================================================
# 07_fracminhash_pairwise.sh
#
# Wrapper: compute FracMinHash (kmer-sketch) pairwise similarities for all
# genome pairs that passed the sourmash 0.001 threshold
# (from gtdb_pairwise_containment_thr0001.csv).
#
# Only pairs where both genomes have been sketched are evaluated.
# Run 04_fracminhash_sketch.sh first (optionally with TEST_N to sketch a
# subset; this script will then compute only the pairs within that subset).
#
# Full run:    bash 07_fracminhash_pairwise.sh
# After test:  bash 07_fracminhash_pairwise.sh   (same command — adapts automatically)
# ==============================================================================
set -euo pipefail

BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"

SKETCH_DIR="${GTDB_DIR}/fracminhash_sketches"
CANDIDATES="${GTDB_DIR}/gtdb_pairwise_containment_thr0001.csv"
OUTPUT_DIR="${GTDB_DIR}/fracminhash_pairwise_thr0001"

CONFIG="${BASE_DIR}/config.json"
CORES=$(python3 -c "import json; print(json.load(open('${CONFIG}'))['pairwise_cores'])")

echo "Starting FracMinHash (kmer-sketch) pairwise computation ..."
echo "  Sketch dir  : ${SKETCH_DIR}"
echo "  Candidates  : ${CANDIDATES}"
echo "  Output dir  : ${OUTPUT_DIR}"
echo "  Cores       : ${CORES}"
echo ""

python3 "${BASE_DIR}/scripts/07_fracminhash_pairwise.py" \
    --sketch-dir "${SKETCH_DIR}" \
    --candidates "${CANDIDATES}" \
    --output     "${OUTPUT_DIR}" \
    --cores      "${CORES}"
