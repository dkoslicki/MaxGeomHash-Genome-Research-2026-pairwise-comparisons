#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Absolute paths
# ---------------------------------------------------------------------------
BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
GTDB_DIR="${BASE_DIR}/data/GTDB"
SKETCH_IN="${GTDB_DIR}/gtdb_genomes_reps_r226_sketches.siz.zip"
PAIRWISE_OUT="${GTDB_DIR}/gtdb_pairwise_containment.csv"

# ---------------------------------------------------------------------------
# Parameters — adjust as needed
# ---------------------------------------------------------------------------
CORES=500       # physical cores available on the H200 node
KSIZE=31
# Only emit pairs whose max-containment is at or above this threshold.
# At 143k genomes, every pair above 0 would be ~10 billion rows.
# 0.01 = 1% shared k-mers; tune down if you want more distant relationships.
THRESHOLD=0.01

# ---------------------------------------------------------------------------
# Timing + resource tracking
# ---------------------------------------------------------------------------
SECONDS=0
TIME_LOG=$(mktemp)
trap 'rm -f "${TIME_LOG}"' EXIT

# ---------------------------------------------------------------------------
# Sanity check
# ---------------------------------------------------------------------------
if [[ ! -f "${SKETCH_IN}" ]]; then
    echo "ERROR: Sketch file not found: ${SKETCH_IN}" >&2
    echo "       Run make_sourmash_sketches.sh first." >&2
    exit 1
fi

# Quick inventory so we know what we're dealing with before committing
echo "Summarising sketch collection ..."
sourmash sig summarize "${SKETCH_IN}"
echo ""

# ---------------------------------------------------------------------------
# Rough scale warning
# N*(N-1)/2 pairs — just informational
# ---------------------------------------------------------------------------
N_SKETCHES=$(sourmash sig summarize "${SKETCH_IN}" 2>/dev/null \
             | grep -oP '(?<=loaded )\d+(?= sketches)' || echo "unknown")
if [[ "${N_SKETCHES}" =~ ^[0-9]+$ ]]; then
    N_PAIRS=$(( N_SKETCHES * (N_SKETCHES - 1) / 2 ))
    echo "INFO: ${N_SKETCHES} sketches → ${N_PAIRS} unique pairs to evaluate."
    echo "INFO: Only pairs with max-containment >= ${THRESHOLD} will be written."
    echo ""
fi

# ---------------------------------------------------------------------------
# Check that the branchwater pairwise subcommand is available
# ---------------------------------------------------------------------------
if ! sourmash scripts pairwise --help &>/dev/null; then
    echo "ERROR: 'sourmash scripts pairwise' not found." >&2
    echo "       Install the branchwater plugin:" >&2
    echo "         conda install -c conda-forge sourmash-plugin-branchwater" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Run all-vs-all pairwise containment
#
# Key flags:
#   --cores          multithreaded Rust backend
#   -k               k-mer size (must match what was used in sketching)
#   --threshold      min max-containment to emit; keeps output manageable
#   --write-all      (omitted intentionally) — would emit self-comparisons too
#
# Output CSV columns (branchwater pairwise):
#   query_name, query_md5, match_name, match_md5,
#   containment, max_containment, jaccard, intersect_hashes, ksize, scaled,
#   query_containment_ani, match_containment_ani, average_containment_ani,
#   max_containment_ani
# ---------------------------------------------------------------------------
echo "Starting pairwise containment computation ..."
echo "Output: ${PAIRWISE_OUT}"
echo ""

/usr/bin/time -v \
    sourmash scripts pairwise \
        "${SKETCH_IN}" \
        --output "${PAIRWISE_OUT}" \
        --ksize "${KSIZE}" \
        --threshold "${THRESHOLD}" \
        --cores "${CORES}" \
    2> >(tee "${TIME_LOG}" >&2)

echo ""
echo "Pairwise computation complete."

# ---------------------------------------------------------------------------
# Quick sanity check on the output
# ---------------------------------------------------------------------------
if [[ -f "${PAIRWISE_OUT}" ]]; then
    N_ROWS=$(( $(wc -l < "${PAIRWISE_OUT}") - 1 ))
    echo "Pairs written (above threshold ${THRESHOLD}): ${N_ROWS}"
else
    echo "WARNING: Output file not found — something may have gone wrong." >&2
fi

# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------
elapsed=${SECONDS}
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

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

output_size=$(du -sh "${PAIRWISE_OUT}" 2>/dev/null | cut -f1 || echo "unavailable")

echo ""
echo "============================== Run Summary =============================="
printf "  Wall-clock time   : %02d:%02d:%02d (hh:mm:ss)\n" \
       "${hours}" "${minutes}" "${seconds}"
printf "  Peak RAM usage    : %s\n"    "${peak_ram_human}"
printf "  Pairs above thr.  : %s\n"   "${N_ROWS:-unknown}"
printf "  Output file       : %s\n"   "${PAIRWISE_OUT}"
printf "  Output file size  : %s\n"   "${output_size}"
echo "========================================================================="
echo ""
echo "Downstream suggestions:"
echo "  # Convert to numpy matrix for sourmash plot:"
echo "  sourmash scripts pairwise_to_matrix ${PAIRWISE_OUT} -o gtdb_matrix.numpy"
echo ""
echo "  # MDS or t-SNE plot (requires sourmash-plugin-betterplot):"
echo "  sourmash scripts mds2 ${PAIRWISE_OUT} -o gtdb_mds.png"
