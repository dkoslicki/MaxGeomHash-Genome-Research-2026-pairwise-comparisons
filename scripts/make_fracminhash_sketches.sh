#!/usr/bin/env bash
# ==============================================================================
# make_fracminhash_sketches.sh
#
# Create FracMinHash sketches of all GTDB representative genomes using the
# kmer-sketch binary (replaces make_sourmash_sketches.sh).
#
# Reads scale from config.json (fracminhash.scale = 0.01 → scaled=100).
# Skips sketching entirely if all expected sketches already exist.
#
# Usage:
#   bash scripts/make_fracminhash_sketches.sh
#
# Delegates to 04_fracminhash_sketch.sh; this wrapper just adds the skip-if-done
# guard and the summary reporting expected by the master pipeline.
# ==============================================================================
set -euo pipefail

BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"
SKETCH_DIR="${GTDB_DIR}/fracminhash_sketches"
MANIFEST="${GTDB_DIR}/manysketch_uncompressed.csv"
CONFIG="${BASE_DIR}/config.json"
STATS_FILE="${SKETCH_DIR}/sketch_run_stats.json"

SECONDS=0

# ---------------------------------------------------------------------------
# Derive expected sketch count from manifest
# ---------------------------------------------------------------------------
if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: uncompressed manifest not found: ${MANIFEST}" >&2
    echo "       Run scripts/00_decompress_genomes.sh first." >&2
    exit 1
fi

N_EXPECTED=$(( $(wc -l < "${MANIFEST}") - 1 ))

# ---------------------------------------------------------------------------
# Check if all sketches already exist
# ---------------------------------------------------------------------------
N_EXISTING=0
if [[ -d "${SKETCH_DIR}" ]]; then
    N_EXISTING=$(find "${SKETCH_DIR}" -maxdepth 1 -name "*.fracminhash.sketch" | wc -l)
fi

if [[ "${N_EXISTING}" -ge "${N_EXPECTED}" && -f "${STATS_FILE}" ]]; then
    echo "All ${N_EXISTING} FracMinHash sketches already exist in ${SKETCH_DIR} — skipping."
    echo ""
    echo "Existing run stats:"
    cat "${STATS_FILE}"
    echo ""

    # Still emit timing summary for the master pipeline log
    elapsed=${SECONDS}
    hours=$(( elapsed / 3600 ))
    minutes=$(( (elapsed % 3600) / 60 ))
    seconds_left=$(( elapsed % 60 ))
    sketch_dir_size=$(du -sh "${SKETCH_DIR}" 2>/dev/null | cut -f1 || echo "unavailable")

    echo ""
    echo "============================== Run Summary =============================="
    printf "  Wall-clock time   : %02d:%02d:%02d (hh:mm:ss)  [skipped — already done]\n" \
           "${hours}" "${minutes}" "${seconds_left}"
    printf "  Sketches present  : %d\n" "${N_EXISTING}"
    printf "  Sketch directory  : %s  (%s)\n" "${SKETCH_DIR}" "${sketch_dir_size}"
    echo "========================================================================="
    exit 0
fi

# ---------------------------------------------------------------------------
# Need to (re-)sketch
# ---------------------------------------------------------------------------
echo "Found ${N_EXISTING} / ${N_EXPECTED} sketches. Running 04_fracminhash_sketch.sh ..."
echo ""

bash "${BASE_DIR}/scripts/04_fracminhash_sketch.sh"
