#!/usr/bin/env bash
# ==============================================================================
# 03_alphamaxgeom_sketch.sh
#
# Create AlphaMaxGeomHash sketches of all GTDB representative genomes.
#
# Algorithm:  alphamaxgeom  (W=64, alpha=0.45)
# K-mer size: 31, canonical k-mers, seed=42
#
# Genomes are stored compressed (.fna.gz).  Each genome is decompressed to
# a per-job temp file, sketched, and the temp file immediately deleted.
# Compressed originals are never modified.
#
# Outputs:
#   data/GTDB/alphamaxgeom_sketches/{genome_id}.alphamaxgeom.sketch  (one per genome)
#   data/GTDB/alphamaxgeom_sketches/sketch_run_stats.json            (timing + disk usage)
#
# Prerequisites:
#   GNU parallel  (apt: parallel)
#   sketch binary compiled in scripts/kmer-sketch/bin/
#
# Test mode — sketch only the first N genomes:
#   TEST_N=200 bash 03_alphamaxgeom_sketch.sh
#
# Full run — sketch all genomes (unset or 0 means "all"):
#   bash 03_alphamaxgeom_sketch.sh
# ==============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Paths (absolute throughout — no cd required)
# ---------------------------------------------------------------------------
BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"
MANIFEST="${GTDB_DIR}/manysketch.csv"           # name,genome_filename,protein_filename
SKETCH_DIR="${GTDB_DIR}/alphamaxgeom_sketches"
SKETCH_BIN="${BASE_DIR}/scripts/kmer-sketch/bin/sketch"

# ---------------------------------------------------------------------------
# AlphaMaxGeomHash parameters — made explicit even where they equal defaults
# ---------------------------------------------------------------------------
ALGO="alphamaxgeom"
KMER=31
W=64          # number of hash buckets
ALPHA=0.45    # geometric sampling parameter
SEED=42

# ---------------------------------------------------------------------------
# Parallelism — one core per job (sketching is CPU-bound, low memory usage)
# ---------------------------------------------------------------------------
PARALLEL_JOBS=192

# ---------------------------------------------------------------------------
# Test mode: set TEST_N to a positive integer to process only the first N
# genomes from the manifest.  Leave at 0 (default) to process all genomes.
# Override at the command line: TEST_N=500 bash 03_alphamaxgeom_sketch.sh
# ---------------------------------------------------------------------------
TEST_N="${TEST_N:-0}"

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
for cmd in parallel gunzip; do
    if ! command -v "${cmd}" &>/dev/null; then
        echo "ERROR: '${cmd}' not found in PATH." >&2
        exit 1
    fi
done

if [[ ! -f "${SKETCH_BIN}" ]]; then
    echo "ERROR: sketch binary not found: ${SKETCH_BIN}" >&2
    echo "       Run 'make' inside scripts/kmer-sketch/ first." >&2
    exit 1
fi

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manysketch manifest not found: ${MANIFEST}" >&2
    exit 1
fi

mkdir -p "${SKETCH_DIR}"

# ---------------------------------------------------------------------------
# Timing and cleanup
# ---------------------------------------------------------------------------
SECONDS=0
TIME_LOG=$(mktemp)
GENOME_LIST=$(mktemp)
GENOME_LIST_TODO=$(mktemp)
trap 'rm -f "${TIME_LOG}" "${GENOME_LIST}" "${GENOME_LIST_TODO}"' EXIT

# ---------------------------------------------------------------------------
# Extract genome paths from the manifest (column 2: genome_filename)
# ---------------------------------------------------------------------------
if [[ "${TEST_N}" -gt 0 ]]; then
    echo "TEST MODE: restricting to first ${TEST_N} genomes."
    # Use awk for line-limiting instead of tail|head|awk — the tail|head combination
    # causes SIGPIPE on tail under set -o pipefail, making the script exit early.
    awk -v n="${TEST_N}" -F',' 'NR > 1 && NR <= n+1 {print $2}' "${MANIFEST}" \
        > "${GENOME_LIST}"
else
    awk -F',' 'NR > 1 {print $2}' "${MANIFEST}" > "${GENOME_LIST}"
fi

N_TOTAL=$(wc -l < "${GENOME_LIST}")
echo "Found ${N_TOTAL} genome files to process."

# ---------------------------------------------------------------------------
# Skip genomes whose sketch already exists (resumable)
# ---------------------------------------------------------------------------
while IFS= read -r genome_path; do
    genome_name=$(basename "${genome_path}" \
        | sed 's/\.\(fna\|fa\|fasta\)\.gz$//; s/\.gz$//')
    out_sketch="${SKETCH_DIR}/${genome_name}.alphamaxgeom.sketch"
    if [[ ! -f "${out_sketch}" ]]; then
        echo "${genome_path}"
    fi
done < "${GENOME_LIST}" > "${GENOME_LIST_TODO}"

N_TODO=$(wc -l < "${GENOME_LIST_TODO}")
echo "${N_TODO} genomes need sketching ($(( N_TOTAL - N_TODO )) already done)."

if [[ "${N_TODO}" -eq 0 ]]; then
    echo "All genomes already sketched — nothing to do."
    exit 0
fi

# ---------------------------------------------------------------------------
# Worker function: decompress one genome and sketch it
# Exported so GNU parallel can call it in a subshell
# ---------------------------------------------------------------------------
export SKETCH_DIR SKETCH_BIN ALGO KMER W ALPHA SEED

sketch_genome() {
    local genome_path="$1"
    local genome_name
    genome_name=$(basename "${genome_path}" \
        | sed 's/\.\(fna\|fa\|fasta\)\.gz$//; s/\.gz$//')
    local out_sketch="${SKETCH_DIR}/${genome_name}.alphamaxgeom.sketch"

    # Per-job temp file: unique name avoids conflicts among parallel workers
    local tmp_fasta
    tmp_fasta=$(mktemp --suffix=".fasta")
    trap "rm -f '${tmp_fasta}'" RETURN   # always delete, even on error

    # Decompress with -c (stdout); never touches the original .gz file
    if ! gunzip -c "${genome_path}" > "${tmp_fasta}" 2>/dev/null; then
        echo "ERROR: decompression failed: ${genome_path}" >&2
        return 1
    fi

    # Create the AlphaMaxGeomHash sketch
    if ! "${SKETCH_BIN}" \
            --input    "${tmp_fasta}" \
            --kmer     "${KMER}" \
            --algo     "${ALGO}" \
            --w        "${W}" \
            --alpha    "${ALPHA}" \
            --seed     "${SEED}" \
            --canonical \
            --output   "${out_sketch}" 2>/dev/null; then
        echo "ERROR: sketching failed: ${genome_name}" >&2
        rm -f "${out_sketch}"   # remove partial output
        return 1
    fi

    echo "DONE ${genome_name}"
}
export -f sketch_genome

# ---------------------------------------------------------------------------
# Run in parallel; capture /usr/bin/time -v for peak RAM
# ---------------------------------------------------------------------------
echo ""
echo "Sketching ${N_TODO} genomes with ${PARALLEL_JOBS} parallel jobs ..."
echo "  algo=${ALGO}, k=${KMER}, W=${W}, alpha=${ALPHA}, seed=${SEED}"
echo ""

/usr/bin/time -v \
    parallel \
        --jobs "${PARALLEL_JOBS}" \
        --line-buffer \
        --halt soon,fail=1 \
        sketch_genome {} \
        < "${GENOME_LIST_TODO}" \
    2> >(tee "${TIME_LOG}" >&2)

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
elapsed=${SECONDS}
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

N_SKETCHED=$(find "${SKETCH_DIR}" -maxdepth 1 -name "*.alphamaxgeom.sketch" | wc -l)
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

# Write JSON stats for downstream comparison (heatmap generation etc.)
STATS_FILE="${SKETCH_DIR}/sketch_run_stats.json"
python3 - <<PYEOF
import json
stats = {
    "script":                "03_alphamaxgeom_sketch.sh",
    "algo":                  "${ALGO}",
    "kmer":                  ${KMER},
    "W":                     ${W},
    "alpha":                 ${ALPHA},
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
printf "  Parameters        : algo=%s  k=%d  W=%d  alpha=%s  seed=%d\n" \
       "${ALGO}" "${KMER}" "${W}" "${ALPHA}" "${SEED}"
echo "========================================================================="
echo ""
echo "Next step:"
echo "  python3 ${BASE_DIR}/scripts/04_alphamaxgeom_pairwise.py \\"
echo "      --sketch-dir ${SKETCH_DIR} \\"
echo "      --candidates ${GTDB_DIR}/gtdb_pairwise_containment.csv \\"
echo "      --output     ${GTDB_DIR}/alphamaxgeom_pairwise \\"
echo "      --cores      192"
echo ""
echo "Or run the wrapper:"
echo "  bash ${BASE_DIR}/scripts/04_alphamaxgeom_pairwise.sh"
