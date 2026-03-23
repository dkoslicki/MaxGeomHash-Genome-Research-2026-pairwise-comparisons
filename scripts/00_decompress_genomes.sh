#!/usr/bin/env bash
# ==============================================================================
# 00_decompress_genomes.sh
#
# ONE-TIME SETUP: Decompress all GTDB representative genome .fna.gz files
# into a flat directory of plain .fna files, then generate a new sourmash-
# compatible manifest (manysketch_uncompressed.csv) pointing to them.
#
# This step exists because the kmer-sketch binary does not support stdin /
# piped input, so the sketch scripts previously had to decompress each genome
# to a per-job temp file.  With 143,614 small genomes that creates massive
# process-spawn overhead and leaves CPUs mostly idle.  Storing uncompressed
# files lets the sketch binary read them directly with zero decompression
# overhead per sketch job.
#
# Run this ONCE before running any sketch script.  Subsequent runs skip files
# that already exist, making it safe to re-run after partial failures.
#
# Outputs:
#   data/GTDB/gtdb_genomes_reps_r226_uncompressed/  (one .fna per genome)
#   data/GTDB/manysketch_uncompressed.csv            (updated manifest)
#
# Prerequisites:
#   GNU parallel  (apt: parallel)
#   gunzip
#
# Usage:
#   bash scripts/00_decompress_genomes.sh
# ==============================================================================
set -euo pipefail

BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"
MANIFEST="${GTDB_DIR}/manysketch.csv"
UNCOMP_DIR="${GTDB_DIR}/gtdb_genomes_reps_r226_uncompressed"
UNCOMP_MANIFEST="${GTDB_DIR}/manysketch_uncompressed.csv"

CONFIG="${BASE_DIR}/config.json"
if [[ ! -f "${CONFIG}" ]]; then
    echo "ERROR: config.json not found at ${CONFIG}" >&2; exit 1
fi
PARALLEL_JOBS=$(python3 -c "import json; print(json.load(open('${CONFIG}'))['parallel_jobs'])")

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
for cmd in parallel gunzip; do
    if ! command -v "${cmd}" &>/dev/null; then
        echo "ERROR: '${cmd}' not found in PATH." >&2; exit 1
    fi
done

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: manifest not found: ${MANIFEST}" >&2; exit 1
fi

mkdir -p "${UNCOMP_DIR}"

# ---------------------------------------------------------------------------
# Count work to do
# ---------------------------------------------------------------------------
N_TOTAL=$(awk -F',' 'NR > 1' "${MANIFEST}" | wc -l)
N_EXISTING=$(find "${UNCOMP_DIR}" -maxdepth 1 -name "*.fna" | wc -l)
N_TODO=$(( N_TOTAL - N_EXISTING ))

echo "========================================================"
echo "  00_decompress_genomes.sh — One-time genome decompression"
echo "  Total genomes  : ${N_TOTAL}"
echo "  Already done   : ${N_EXISTING}"
echo "  To decompress  : ${N_TODO}"
echo "  Output dir     : ${UNCOMP_DIR}"
echo "  Parallel jobs  : ${PARALLEL_JOBS}"
echo "========================================================"

if [[ "${N_TODO}" -eq 0 ]]; then
    echo "All genomes already decompressed — skipping decompression."
else
    # -----------------------------------------------------------------------
    # Worker function: decompress one .fna.gz to the flat output directory.
    # Skips files that already exist (safe to re-run after partial failure).
    # -----------------------------------------------------------------------
    export UNCOMP_DIR

    decompress_genome() {
        local gz_path="$1"
        local genome_name
        genome_name=$(basename "${gz_path}" | sed 's/\.gz$//')   # e.g. GCA_xxx.fna
        local out_fasta="${UNCOMP_DIR}/${genome_name}"

        if [[ -f "${out_fasta}" ]]; then
            echo "SKIP ${genome_name}"
            return 0
        fi

        if ! gunzip -c "${gz_path}" > "${out_fasta}" 2>/dev/null; then
            echo "ERROR: decompression failed: ${gz_path}" >&2
            rm -f "${out_fasta}"
            return 1
        fi
        echo "DONE ${genome_name}"
    }
    export -f decompress_genome

    SECONDS=0
    awk -F',' 'NR > 1 {print $2}' "${MANIFEST}" \
        | parallel \
            --jobs "${PARALLEL_JOBS}" \
            --line-buffer \
            --halt soon,fail=1 \
            decompress_genome {}

    elapsed=${SECONDS}
    hours=$(( elapsed / 3600 ))
    minutes=$(( (elapsed % 3600) / 60 ))
    seconds=$(( elapsed % 60 ))
    echo ""
    printf "  Decompression wall-clock time : %02d:%02d:%02d (hh:mm:ss)\n" \
           "${hours}" "${minutes}" "${seconds}"
fi

# ---------------------------------------------------------------------------
# Generate (or refresh) manysketch_uncompressed.csv
# ---------------------------------------------------------------------------
echo ""
echo "Writing uncompressed manifest: ${UNCOMP_MANIFEST}"
python3 - <<PYEOF
import csv, pathlib

manifest_in  = "${MANIFEST}"
manifest_out = "${UNCOMP_MANIFEST}"
uncomp_dir   = pathlib.Path("${UNCOMP_DIR}")

with open(manifest_in, newline="") as fin, \
     open(manifest_out, "w", newline="") as fout:

    reader = csv.reader(fin)
    writer = csv.writer(fout)

    header = next(reader)
    writer.writerow(header)

    for row in reader:
        if len(row) < 2:
            continue
        name     = row[0]
        gz_path  = row[1]
        fna_name = pathlib.Path(gz_path).name.rstrip(".gz").rstrip("")
        # Strip .gz from the filename
        fna_name = pathlib.Path(gz_path).stem   # e.g. GCA_xxx.fna
        # pathlib.stem only strips one suffix, .fna.gz needs two passes
        # Use string replacement instead:
        fna_name = pathlib.Path(gz_path).name
        if fna_name.endswith(".gz"):
            fna_name = fna_name[:-3]   # strip .gz → leaves .fna

        fna_path = str(uncomp_dir / fna_name)
        protein  = row[2] if len(row) > 2 else ""
        writer.writerow([name, fna_path, protein])

print(f"Manifest written: {manifest_out}")
PYEOF

# ---------------------------------------------------------------------------
# Final count
# ---------------------------------------------------------------------------
N_DECOMPRESSED=$(find "${UNCOMP_DIR}" -maxdepth 1 -name "*.fna" | wc -l)
uncomp_size=$(du -sh "${UNCOMP_DIR}" 2>/dev/null | cut -f1 || echo "unavailable")

echo ""
echo "========================================================"
echo "  Genomes decompressed : ${N_DECOMPRESSED} / ${N_TOTAL}"
echo "  Uncompressed size    : ${uncomp_size}"
echo "  Manifest             : ${UNCOMP_MANIFEST}"
echo ""
echo "  Next step: run sketch scripts (03, 04, 05)"
echo "  They will read from ${UNCOMP_DIR}"
echo "========================================================"
