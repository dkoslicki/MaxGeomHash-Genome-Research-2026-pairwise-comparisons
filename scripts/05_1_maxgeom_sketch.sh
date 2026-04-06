#!/usr/bin/env bash
# ==============================================================================
# 05_1_maxgeom_sketch.sh
#
# Create MaxGeomHash sketches of all GTDB representative genomes.
#
# Algorithm:  maxgeom  (W=64, b=90)
# K-mer size: 31, canonical k-mers, seed=42
#
# Reads from the pre-decompressed flat FASTA directory created by
# 00_decompress_genomes.sh.  Run that script once before running this one.
# The sketch binary reads .fna files directly — no per-job decompression,
# no temp files — which keeps all CPUs busy.
#
# Outputs:
#   data/GTDB/maxgeom_sketches/{genome_id}.maxgeom.sketch  (one per genome)
#   data/GTDB/maxgeom_sketches/sketch_run_stats.json            (timing + disk usage)
#
# Prerequisites:
#   bash scripts/00_decompress_genomes.sh   (once, before first sketch run)
#   GNU parallel  (apt: parallel)
#   sketch binary compiled in scripts/kmer-sketch/bin/
#
# Test mode — sketch only the first N genomes:
#   TEST_N=200 bash 05_1_maxgeom_sketch.sh
#
# Full run — sketch all genomes (unset or 0 means "all"):
#   bash 05_1_maxgeom_sketch.sh
# ==============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Paths (absolute throughout — no cd required)
# ---------------------------------------------------------------------------
BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"
MANIFEST="${GTDB_DIR}/manysketch_uncompressed.csv"   # points to plain .fna files
SKETCH_DIR="${GTDB_DIR}/maxgeom_sketches"
SKETCH_BIN="${BASE_DIR}/scripts/kmer-sketch/bin/sketch"

# ---------------------------------------------------------------------------
# Parameters — read from config.json at project root
# ---------------------------------------------------------------------------
CONFIG="${BASE_DIR}/config.json"
if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.json not found at ${CONFIG}" >&2; exit 1
fi
cfg() { python3 -c "import json; print(json.load(open('${CONFIG}'))$1)"; }

ALGO="maxgeom"
KMER=$(cfg "['kmer']")
SEED=$(cfg "['seed']")
W=$(cfg "['maxgeom']['w']")
B=$(cfg "['maxgeom']['b']")
PARALLEL_JOBS=$(cfg "['parallel_jobs']")

# ---------------------------------------------------------------------------
# Test mode
# ---------------------------------------------------------------------------
TEST_N="${TEST_N:-0}"

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if ! command -v parallel &>/dev/null; then
    echo "ERROR: 'parallel' not found in PATH." >&2; exit 1
fi

if [[ ! -f "${SKETCH_BIN}" ]]; then
    echo "ERROR: sketch binary not found: ${SKETCH_BIN}" >&2
    echo "       Run 'make' inside scripts/kmer-sketch/ first." >&2
    exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: uncompressed manifest not found: ${MANIFEST}" >&2
    echo "       Run scripts/00_decompress_genomes.sh first." >&2
    exit 1
fi

mkdir -p "${SKETCH_DIR}"

# ---------------------------------------------------------------------------
# Timing and cleanup
# ---------------------------------------------------------------------------
SECONDS=0
GENOME_LIST=$(mktemp)
GENOME_LIST_TODO=$(mktemp)
trap 'rm -f "${GENOME_LIST}" "${GENOME_LIST_TODO}"' EXIT

# ---------------------------------------------------------------------------
# Extract genome paths from the manifest (column 2: genome_filename)
# ---------------------------------------------------------------------------
if [[ "${TEST_N}" -gt 0 ]]; then
    echo "TEST MODE: restricting to first ${TEST_N} genomes."
    awk -v n="${TEST_N}" -F',' 'NR > 1 && NR <= n+1 {print $2}' "${MANIFEST}" \
        > "${GENOME_LIST}"
else
    awk -F',' 'NR > 1 {print $2}' "${MANIFEST}" > "${GENOME_LIST}"
fi

N_TOTAL=$(wc -l < "${GENOME_LIST}")
echo "Found ${N_TOTAL} genome files to process."

# ---------------------------------------------------------------------------
# Always re-sketch all genomes — no skip logic.
# ---------------------------------------------------------------------------
cp "${GENOME_LIST}" "${GENOME_LIST_TODO}"
N_TODO=${N_TOTAL}
echo "Sketching all ${N_TODO} genomes (existing sketches will be overwritten)."

# ---------------------------------------------------------------------------
# Worker function: sketch a batch of genomes in a single parallel slot.
# ---------------------------------------------------------------------------
export SKETCH_DIR SKETCH_BIN ALGO KMER W B SEED

sketch_genome_batch() {
    for genome_path in "$@"; do
        local genome_name
        genome_name=$(basename "${genome_path}" | sed 's/\.\(fna\|fa\|fasta\)$//')
        local out_sketch="${SKETCH_DIR}/${genome_name}.maxgeom.sketch"

        if ! "${SKETCH_BIN}" \
                --input    "${genome_path}" \
                --kmer     "${KMER}" \
                --algo     "${ALGO}" \
                --w        "${W}" \
                --b        "${B}" \
                --seed     "${SEED}" \
                --canonical \
                --output   "${out_sketch}" 2>/dev/null; then
            echo "ERROR: sketching failed: ${genome_name}" >&2
            rm -f "${out_sketch}"
            return 1
        fi

        echo "DONE ${genome_name}"
    done
}
export -f sketch_genome_batch

BATCH_SIZE=$(( (N_TODO + PARALLEL_JOBS - 1) / PARALLEL_JOBS ))

# ---------------------------------------------------------------------------
# Run in parallel; capture /usr/bin/time -v for peak RAM
# ---------------------------------------------------------------------------
echo ""
echo "Sketching ${N_TODO} genomes with ${PARALLEL_JOBS} parallel jobs ..."
echo "  algo=${ALGO}, k=${KMER}, W=${W}, b=${B}, seed=${SEED}"
echo "  batch size per slot: ${BATCH_SIZE} genomes"
echo ""

TIME_LOG=$(mktemp)
trap 'rm -f "${GENOME_LIST}" "${GENOME_LIST_TODO}" "${TIME_LOG}"' EXIT

/usr/bin/time -v \
    parallel \
        --jobs "${PARALLEL_JOBS}" \
        --line-buffer \
        --halt soon,fail=1 \
        -N "${BATCH_SIZE}" \
        sketch_genome_batch \
        < "${GENOME_LIST_TODO}" \
    2> >(tee "${TIME_LOG}" >&2)

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
elapsed=${SECONDS}
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

N_SKETCHED=$(find "${SKETCH_DIR}" -maxdepth 1 -name "*.maxgeom.sketch" | wc -l)
sketch_dir_size=$(du -sh "${SKETCH_DIR}" 2>/dev/null | cut -f1 || echo "unavailable")

peak_ram_kb=$(grep -m1 "Maximum resident set size" "${TIME_LOG}" \
    | awk '{print $NF}' || echo "0")
if [[ -n "${peak_ram_kb}" && "${peak_ram_kb}" -gt 0 ]]; then
    peak_ram_human=$(awk -v kb="${peak_ram_kb}" 'BEGIN {
        if (kb >= 1024*1024)  printf "%.2f GB", kb/1024/1024
        else if (kb >= 1024)  printf "%.2f MB", kb/1024
        else                  printf "%d KB",   kb
    }')
else
    peak_ram_human="unavailable"
fi

STATS_FILE="${SKETCH_DIR}/sketch_run_stats.json"
python3 - <<PYEOF
import json
stats = {
    "script":                "05_1_maxgeom_sketch.sh",
    "algo":                  "${ALGO}",
    "kmer":                  ${KMER},
    "W":                     ${W},
    "b":                     ${B},
    "seed":                  ${SEED},
    "canonical":             True,
    "parallel_jobs":         ${PARALLEL_JOBS},
    "test_n":                ${TEST_N},
    "n_genomes_in_manifest": ${N_TOTAL},
    "n_genomes_sketched":    ${N_SKETCHED},
    "wall_clock_seconds":    ${elapsed},
    "wall_clock_hms":        "${hours}:${minutes}:${seconds}",
    "peak_ram":              "${peak_ram_human}",
    "sketch_dir":            "${SKETCH_DIR}",
    "sketch_dir_disk_usage": "${sketch_dir_size}",
}
with open("${STATS_FILE}", "w") as f:
    json.dump(stats, f, indent=2)
print(f"Run stats written to ${STATS_FILE}")
PYEOF

echo ""
echo "============================== Run Summary =============================="
printf "  Wall-clock time   : %02d:%02d:%02d (hh:mm:ss)\n" "${hours}" "${minutes}" "${seconds}"
printf "  Peak RAM usage    : %s\n"  "${peak_ram_human}"
printf "  Genomes sketched  : %d total (%d new this run)\n" "${N_SKETCHED}" "${N_TODO}"
printf "  Sketch directory  : %s  (%s)\n" "${SKETCH_DIR}" "${sketch_dir_size}"
printf "  Parameters        : algo=%s  k=%d  W=%d  b=%d  seed=%d\n" \
       "${ALGO}" "${KMER}" "${W}" "${B}" "${SEED}"
echo "========================================================================="
echo ""
echo "Next step:"
echo "  bash ${BASE_DIR}/scripts/08_1_maxgeom_pairwise.sh"
