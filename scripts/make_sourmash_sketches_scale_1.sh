#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Absolute paths — no relative directory navigation
# ---------------------------------------------------------------------------
BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"
GENOME_DIR="${GTDB_DIR}/gtdb_genomes_reps_r226"
MANIFEST="${GTDB_DIR}/manysketch.csv"
SKETCH_OUT="${GTDB_DIR}/gtdb_genomes_reps_r226_sketches_scale_1.siz.zip"

# ---------------------------------------------------------------------------
# Timing: bash updates $SECONDS automatically; reset it to 0 at script start
# ---------------------------------------------------------------------------
SECONDS=0

# ---------------------------------------------------------------------------
# Temp file to capture /usr/bin/time -v output (goes to stderr by default)
# Cleaned up automatically on exit, including on error
# ---------------------------------------------------------------------------
TIME_LOG=$(mktemp)
trap 'rm -f "${TIME_LOG}"' EXIT

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if [[ ! -d "${GENOME_DIR}" ]]; then
    echo "ERROR: Genome directory not found: ${GENOME_DIR}" >&2
    echo "       Did you decompress gtdb_genomes_reps_r226.tar.gz first?" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Build the manysketch manifest (skip if it already exists)
# ---------------------------------------------------------------------------
if [[ -f "${MANIFEST}" ]]; then
    echo "Manifest already exists at ${MANIFEST} — skipping creation."
else
    echo "Building manifest at ${MANIFEST} ..."
    echo "name,genome_filename,protein_filename" > "${MANIFEST}"

    # Use find so we handle nested subdirectories and filenames with spaces
    find "${GENOME_DIR}" -name "*.gz" -type f | sort | while IFS= read -r filepath; do
        abs_path=$(readlink -f "${filepath}")
        # Strip .fna.gz / .fa.gz / .fasta.gz / .gz for a clean sample name
        name=$(basename "${filepath}" | sed 's/\.\(fna\|fa\|fasta\)\.gz$//; s/\.gz$//')
        printf '%s,%s,\n' "${name}" "${abs_path}"
    done >> "${MANIFEST}"

    n_genomes=$(( $(wc -l < "${MANIFEST}") - 1 ))
    echo "Manifest created with ${n_genomes} genome entries."
fi

# ---------------------------------------------------------------------------
# Sketch all genomes
# Redirect /usr/bin/time -v's stderr to TIME_LOG; tee it back to stderr so
# it still appears in the terminal/log as the job runs
# ---------------------------------------------------------------------------
echo "Starting sourmash manysketch (this will take a while) ..."
/usr/bin/time -v sourmash scripts manysketch "${MANIFEST}" \
    --output "${SKETCH_OUT}" \
    --param-string "dna,k=31,scaled=1,noabund" \
    --cores 700 \
    --force \
    2> >(tee "${TIME_LOG}" >&2)   # capture time -v stderr AND keep it visible

echo "Done. Sketches written to ${SKETCH_OUT}"

# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------
elapsed=${SECONDS}
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

# Peak RAM from /usr/bin/time -v ("Maximum resident set size (kbytes): N")
peak_ram_kb=$(grep -m1 "Maximum resident set size" "${TIME_LOG}" \
              | awk '{print $NF}')
if [[ -n "${peak_ram_kb}" && "${peak_ram_kb}" -gt 0 ]]; then
    peak_ram_human=$(awk -v kb="${peak_ram_kb}" 'BEGIN {
        if (kb >= 1024*1024)      printf "%.2f GB", kb/1024/1024
        else if (kb >= 1024)      printf "%.2f MB", kb/1024
        else                      printf "%d KB",   kb
    }')
else
    peak_ram_human="unavailable"
fi

# Disk usage: output sketch file and input genome directory
sketch_size=$(du -sh "${SKETCH_OUT}" 2>/dev/null | cut -f1 || echo "unavailable")
genome_dir_size=$(du -sh "${GENOME_DIR}" 2>/dev/null | cut -f1 || echo "unavailable")

echo ""
echo "============================== Run Summary =============================="
printf "  Wall-clock time : %02d:%02d:%02d (hh:mm:ss)\n" \
       "${hours}" "${minutes}" "${seconds}"
printf "  Peak RAM usage  : %s\n"  "${peak_ram_human}"
printf "  Sketch output   : %s  (%s)\n" "${SKETCH_OUT}" "${sketch_size}"
printf "  Genome dir size : %s  (%s)\n" "${GENOME_DIR}" "${genome_dir_size}"
echo "========================================================================="
