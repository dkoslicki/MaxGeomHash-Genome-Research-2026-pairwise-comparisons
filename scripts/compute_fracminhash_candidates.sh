#!/usr/bin/env bash
# ==============================================================================
# compute_fracminhash_candidates.sh
#
# Wrapper: run all-vs-all FracMinHash pairwise using kmer-sketch sketches to
# produce the candidate pairs CSV (replaces compute_sourmash_pairwise.sh).
#
# Reads configuration from config.json.
# Requires the fracminhash_sketches directory to already be populated
# (run make_fracminhash_sketches.sh / 04_fracminhash_sketch.sh first).
#
# Skip guard: if the output CSV and its stats JSON already exist and record
# the same sketch_dir, threshold, num_shards, and n_genomes as the current
# run, the step is skipped automatically.  Delete the stats JSON or change
# any parameter to force a recompute.
#
# Output:
#   data/GTDB/gtdb_pairwise_containment_thr0001.csv
#   data/GTDB/gtdb_pairwise_containment_thr0001.stats.json
#
# Full run:
#   bash scripts/compute_fracminhash_candidates.sh
# ==============================================================================
set -euo pipefail

BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"
SKETCH_DIR="${GTDB_DIR}/fracminhash_sketches"
OUTPUT_CSV="${GTDB_DIR}/gtdb_pairwise_containment_thr0001.csv"
STATS_JSON="${OUTPUT_CSV%.csv}.stats.json"

CONFIG="${BASE_DIR}/config.json"
cfg() { python3 -c "import json; print(json.load(open('${CONFIG}'))$1)"; }

CORES=$(cfg "['pairwise_cores']")
THRESHOLD=0.001
NUM_SHARDS=128

SECONDS=0

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if [[ ! -d "${SKETCH_DIR}" ]]; then
    echo "ERROR: Sketch directory not found: ${SKETCH_DIR}" >&2
    echo "       Run make_fracminhash_sketches.sh (or 04_fracminhash_sketch.sh) first." >&2
    exit 1
fi

N_SKETCHES=$(find "${SKETCH_DIR}" -maxdepth 1 -name "*.fracminhash.sketch" | wc -l)
if [[ "${N_SKETCHES}" -eq 0 ]]; then
    echo "ERROR: No .fracminhash.sketch files found in ${SKETCH_DIR}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Skip guard: check whether a previous run with matching parameters succeeded
# ---------------------------------------------------------------------------
if [[ -f "${OUTPUT_CSV}" && -f "${STATS_JSON}" ]]; then
    SKIP=$(python3 - <<PYEOF
import json, sys
try:
    with open("${STATS_JSON}") as f:
        s = json.load(f)
    ok = (
        s.get("sketch_dir")  == "${SKETCH_DIR}"   and
        abs(s.get("threshold", -1) - ${THRESHOLD}) < 1e-9 and
        s.get("num_shards")  == ${NUM_SHARDS}      and
        s.get("n_genomes")   == ${N_SKETCHES}
    )
    print("yes" if ok else "no")
except Exception:
    print("no")
PYEOF
)
    if [[ "${SKIP}" == "yes" ]]; then
        echo "All ${N_SKETCHES} FracMinHash candidate pairs already computed with matching parameters — skipping."
        echo ""
        echo "Existing run stats:"
        cat "${STATS_JSON}"
        echo ""

        elapsed=${SECONDS}
        hours=$(( elapsed / 3600 ))
        minutes=$(( (elapsed % 3600) / 60 ))
        secs=$(( elapsed % 60 ))
        output_size=$(du -sh "${OUTPUT_CSV}" 2>/dev/null | cut -f1 || echo "unavailable")
        N_PAIRS=$(python3 -c "import json; print(json.load(open('${STATS_JSON}'))['n_pairs_above_threshold'])")

        echo ""
        echo "============================== Run Summary =============================="
        printf "  Wall-clock time   : %02d:%02d:%02d (hh:mm:ss)  [skipped — already done]\n" \
               "${hours}" "${minutes}" "${secs}"
        printf "  Genomes           : %d\n"  "${N_SKETCHES}"
        printf "  Pairs above thr.  : %s\n" "${N_PAIRS}"
        printf "  Output file       : %s\n" "${OUTPUT_CSV}"
        printf "  Output file size  : %s\n" "${output_size}"
        echo "========================================================================="
        exit 0
    fi
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "FracMinHash all-vs-all candidate detection"
echo "  Sketch dir  : ${SKETCH_DIR}  (${N_SKETCHES} sketches)"
echo "  Output      : ${OUTPUT_CSV}"
echo "  Threshold   : ${THRESHOLD}"
echo "  Num shards  : ${NUM_SHARDS}"
echo "  Cores       : ${CORES}"
echo ""

/usr/bin/time -v \
    python3 "${BASE_DIR}/scripts/compute_fracminhash_candidates.py" \
        --sketch-dir "${SKETCH_DIR}" \
        --output     "${OUTPUT_CSV}" \
        --threshold  "${THRESHOLD}" \
        --cores      "${CORES}" \
        --num-shards "${NUM_SHARDS}" \
    2>&1

echo ""

# ---------------------------------------------------------------------------
# Quick sanity check
# ---------------------------------------------------------------------------
if [[ -f "${OUTPUT_CSV}" ]]; then
    N_ROWS=$(( $(wc -l < "${OUTPUT_CSV}") - 1 ))
    echo "Pairs written (above threshold ${THRESHOLD}): ${N_ROWS}"
else
    echo "WARNING: Output CSV not found — something may have gone wrong." >&2
fi

# ---------------------------------------------------------------------------
# Timing summary
# ---------------------------------------------------------------------------
elapsed=${SECONDS}
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))
output_size=$(du -sh "${OUTPUT_CSV}" 2>/dev/null | cut -f1 || echo "unavailable")

echo ""
echo "============================== Run Summary =============================="
printf "  Wall-clock time   : %02d:%02d:%02d (hh:mm:ss)\n" \
       "${hours}" "${minutes}" "${seconds}"
printf "  Pairs above thr.  : %s\n" "${N_ROWS:-unknown}"
printf "  Output file       : %s\n" "${OUTPUT_CSV}"
printf "  Output file size  : %s\n" "${output_size}"
echo "========================================================================="
