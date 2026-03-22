#!/usr/bin/env bash
# ==============================================================================
# 01_kmc_count.sh
#
# Count k-mers for every GTDB genome using KMC, in parallel.
#
# Produces:
#   - One KMC database (two files: .kmc_pre + .kmc_suf) per genome
#   - kmc_metadata.csv: genome_id, genome_path, db_prefix, n_unique_kmers
#
# Prerequisites: kmc, kmc_tools, GNU parallel
# ==============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration — edit these as needed
# ---------------------------------------------------------------------------
BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GENOME_DIR="${BASE_DIR}/data/GTDB/gtdb_genomes_reps_r226"
KMC_DB_DIR="${BASE_DIR}/data/GTDB/kmc_dbs"      # one subdir per genome
METADATA_CSV="${BASE_DIR}/data/GTDB/kmc_metadata.csv"

KSIZE=31          # k-mer length (must match sourmash sketches for comparability)
KMC_THREADS=4     # threads per KMC job (KMC scales well to ~8; beyond that diminishing returns)
KMC_MEM=20         # RAM per KMC job, in GB
TOTAL_CORES=384   # physical cores available
PARALLEL_JOBS=$(( TOTAL_CORES / KMC_THREADS ))   # = 96 simultaneous KMC jobs

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
for cmd in kmc kmc_tools parallel; do
    if ! command -v "${cmd}" &>/dev/null; then
        echo "ERROR: '${cmd}' not found in PATH." >&2
        exit 1
    fi
done

if [[ ! -d "${GENOME_DIR}" ]]; then
    echo "ERROR: Genome directory not found: ${GENOME_DIR}" >&2
    exit 1
fi

mkdir -p "${KMC_DB_DIR}"

# ---------------------------------------------------------------------------
# Timing
# ---------------------------------------------------------------------------
SECONDS=0
TIME_LOG=$(mktemp)
trap 'rm -f "${TIME_LOG}"' EXIT

# ---------------------------------------------------------------------------
# Build list of genome files to process
# ---------------------------------------------------------------------------
GENOME_LIST=$(mktemp)
find "${GENOME_DIR}" -type f -name "*.gz" | sort > "${GENOME_LIST}"
N_TOTAL=$(wc -l < "${GENOME_LIST}")
echo "Found ${N_TOTAL} genome files to process."

# ---------------------------------------------------------------------------
# Skip genomes whose KMC DB already exists
# ---------------------------------------------------------------------------
GENOME_LIST_TODO=$(mktemp)
while IFS= read -r genome_path; do
    genome_name=$(basename "${genome_path}" | sed 's/\.\(fna\|fa\|fasta\)\.gz$//; s/\.gz$//')
    db_prefix="${KMC_DB_DIR}/${genome_name}"
    if [[ -f "${db_prefix}.kmc_pre" && -f "${db_prefix}.kmc_suf" ]]; then
        : # already done, skip
    else
        echo "${genome_path}"
    fi
done < "${GENOME_LIST}" > "${GENOME_LIST_TODO}"
trap 'rm -f "${TIME_LOG}" "${GENOME_LIST}" "${GENOME_LIST_TODO}"' EXIT

N_TODO=$(wc -l < "${GENOME_LIST_TODO}")
echo "${N_TODO} genomes need k-mer counting ($(( N_TOTAL - N_TODO )) already done)."

# ---------------------------------------------------------------------------
# Worker function: count k-mers for a single genome
# Exported so GNU parallel can call it in a subshell
# ---------------------------------------------------------------------------
export KMC_DB_DIR KSIZE KMC_THREADS KMC_MEM

kmc_count_genome() {
    local genome_path="$1"
    local genome_name
    genome_name=$(basename "${genome_path}" | sed 's/\.\(fna\|fa\|fasta\)\.gz$//; s/\.gz$//')
    local db_prefix="${KMC_DB_DIR}/${genome_name}"

    # Per-job temp dir (avoids conflicts between parallel jobs)
    local tmp_dir
    tmp_dir=$(mktemp -d "${KMC_DB_DIR}/.kmc_tmp_XXXXXX")
    trap "rm -rf '${tmp_dir}'" RETURN

    # Run KMC:
    #   -k${KSIZE}   : k-mer length
    #   -ci1         : include k-mers with count >= 1 (keep singletons — ground truth)
    #   -cs2         : cap counter at 2 (presence/absence; saves disk; sufficient for Jaccard)
    #   -m${KMC_MEM} : RAM limit in GB
    #   -t${KMC_THREADS}: threads
    #   -fm          : multi-fasta format
    #   No -b flag   : canonical k-mers (default; matches sourmash behaviour)
    local kmc_out
    kmc_out=$(kmc \
        -k"${KSIZE}" \
        -ci1 \
        -cs2 \
        -m"${KMC_MEM}" \
        -t"${KMC_THREADS}" \
        -fm \
        "${genome_path}" \
        "${db_prefix}" \
        "${tmp_dir}" 2>&1)

    # Parse unique k-mer count from KMC output
    # KMC prints: "No. of unique k-mers               : 3429876"
    local n_kmers
    n_kmers=$(echo "${kmc_out}" | grep "No. of unique k-mers" | head -1 | awk -F: '{print $NF}' | tr -d ' ')

    if [[ -z "${n_kmers}" || "${n_kmers}" -eq 0 ]]; then
        echo "WARNING: Zero k-mers counted for ${genome_name}. KMC output:" >&2
        echo "${kmc_out}" >&2
        n_kmers=0
    fi

    # Emit one CSV row (tab-separated to allow commas in paths)
    printf '%s\t%s\t%s\t%s\n' "${genome_name}" "${genome_path}" "${db_prefix}" "${n_kmers}"
}
export -f kmc_count_genome

# ---------------------------------------------------------------------------
# Run KMC in parallel via GNU parallel
# --line-buffer: line-buffered output (avoids interleaved partial lines)
# --jobs: how many simultaneous KMC processes
# --halt soon,fail=1: stop all jobs if any single job fails
# ---------------------------------------------------------------------------
echo ""
echo "Running KMC on ${N_TODO} genomes with ${PARALLEL_JOBS} parallel jobs ..."
echo "(${KMC_THREADS} threads each, ${KMC_MEM} GB RAM each)"
echo ""

RESULTS_TMP=$(mktemp)
parallel \
    --jobs "${PARALLEL_JOBS}" \
    --line-buffer \
    --halt soon,fail=1 \
    kmc_count_genome {} \
    < "${GENOME_LIST_TODO}" \
    >> "${RESULTS_TMP}"

# ---------------------------------------------------------------------------
# Merge with any previously completed genomes to produce final metadata CSV
# ---------------------------------------------------------------------------
echo "genome_id,genome_path,db_prefix,n_unique_kmers" > "${METADATA_CSV}"

# Add previously completed genomes (those skipped above)
while IFS= read -r genome_path; do
    genome_name=$(basename "${genome_path}" | sed 's/\.\(fna\|fa\|fasta\)\.gz$//; s/\.gz$//')
    db_prefix="${KMC_DB_DIR}/${genome_name}"
    if [[ -f "${db_prefix}.kmc_pre" ]]; then
        # Try to find existing count in results or leave as 'unknown'
        echo "${genome_name},${genome_path},${db_prefix},unknown"
    fi
done < "${GENOME_LIST}" >> "${METADATA_CSV}"

# Add newly computed genomes (sorted by genome_id for reproducibility)
sort -t$'\t' -k1,1 "${RESULTS_TMP}" | \
    awk -F'\t' '{print $1","$2","$3","$4}' >> "${METADATA_CSV}"

# Deduplicate (prefer the fresh count if both versions exist)
python3 - "${METADATA_CSV}" <<PYEOF
import csv, sys
seen = {}
rows = []
path = sys.argv[1]
with open(path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        gid = row['genome_id']
        # If we have a fresh numeric count, prefer it over 'unknown'
        if gid not in seen or (row['n_unique_kmers'].isdigit()):
            seen[gid] = row
            rows.append(row)
with open(path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['genome_id','genome_path','db_prefix','n_unique_kmers'])
    writer.writeheader()
    writer.writerows(seen.values())
PYEOF

rm -f "${RESULTS_TMP}"

N_WRITTEN=$(( $(wc -l < "${METADATA_CSV}") - 1 ))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
elapsed=${SECONDS}
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))
db_dir_size=$(du -sh "${KMC_DB_DIR}" 2>/dev/null | cut -f1 || echo "unavailable")

echo ""
echo "============================== Run Summary =============================="
printf "  Wall-clock time : %02d:%02d:%02d (hh:mm:ss)\n" "${hours}" "${minutes}" "${seconds}"
printf "  Genomes counted : %d / %d\n" "${N_TODO}" "${N_TOTAL}"
printf "  Metadata CSV    : %s  (%d rows)\n" "${METADATA_CSV}" "${N_WRITTEN}"
printf "  KMC DB dir size : %s  (%s)\n" "${KMC_DB_DIR}" "${db_dir_size}"
echo "========================================================================="
echo ""
echo "Next step:"
echo "  python3 02_kmc_pairwise.py \\"
echo "      --metadata  ${METADATA_CSV} \\"
echo "      --output    ${BASE_DIR}/data/GTDB/kmc_pairwise \\"
echo "      --threshold 0.01 \\"
echo "      --cores     ${TOTAL_CORES}"
echo ""
echo "  # Or, to compute exact values only for candidate pairs from sourmash:"
echo "  python3 02_kmc_pairwise.py \\"
echo "      --metadata   ${METADATA_CSV} \\"
echo "      --candidates ${BASE_DIR}/data/GTDB/gtdb_pairwise_containment.csv \\"
echo "      --output     ${BASE_DIR}/data/GTDB/kmc_pairwise \\"
echo "      --threshold  0.01 \\"
echo "      --cores      ${TOTAL_CORES}"
