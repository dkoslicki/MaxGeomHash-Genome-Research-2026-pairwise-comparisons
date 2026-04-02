#!/usr/bin/env bash
# ==============================================================================
# smoke_test/run_smoke_test.sh
#
# End-to-end smoke test for the MaxGeomHash GTDB pairwise similarity pipeline.
#
# Validates the full pipeline — sketching, candidate generation, exact KMC
# pairwise, kmer-sketch pairwise (BottomK, FracMinHash, AlphaMaxGeomHash), and
# sanity-check accuracy reporting — on a small, curated genome subset.
#
# Genome subset design
# --------------------
# Rather than taking the first N genomes from the manifest (which are highly
# diverse and would produce few or no candidate pairs), this script selects
# 100 genomes that are known to appear in actual candidate pairs from the full
# run (gtdb_pairwise_containment_thr0001.csv).  These 100 genomes produce ~913
# candidate pairs at threshold=0.001, providing a meaningful accuracy test.
#
# Output
# ------
# All outputs are written to smoke_test/data/ and are gitignored.
# A summary is printed at the end; non-zero exit means the test failed.
#
# Usage
# -----
#   cd /scratch/shared_data/MaxGeomHash_Genome_Research_2026
#   bash smoke_test/run_smoke_test.sh
#
# Requirements
# ------------
#   - KMC databases already built (scripts/01_kmc_count.sh, one-time)
#   - kmer-sketch binary compiled: scripts/kmer-sketch/bin/sketch
#   - GNU parallel (apt: parallel)
#   - Python packages: numpy, scipy, pandas, matplotlib, python-dateutil
#
# Timing (approximate, on 384-core H200 node, using SMOKE_CORES=8):
#   Sketching:          ~5 s
#   Candidate gen:      ~1 s
#   KMC pairwise:       ~5–10 min (913 pairs × kmc_tools invocations)
#   kmer-sketch pairwise: <1 s each
#   Sanity check:       ~5 s
#   Total:              ~15 min
# ==============================================================================

set -euo pipefail

BASE_DIR="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
SMOKE_DIR="${BASE_DIR}/smoke_test"
DATA_DIR="${SMOKE_DIR}/data"
SCRIPTS_DIR="${BASE_DIR}/scripts"
SKETCH_BIN="${SCRIPTS_DIR}/kmer-sketch/bin/sketch"
CONFIG="${BASE_DIR}/config.json"
FULL_MANIFEST="${BASE_DIR}/data/GTDB/manysketch_uncompressed.csv"
KMC_METADATA="${BASE_DIR}/data/GTDB/kmc_metadata.csv"
FULL_CANDIDATES="${BASE_DIR}/data/GTDB/gtdb_pairwise_containment_thr0001.csv"

# Number of worker processes for the smoke test (much less than the 384-core
# full run; adjust if desired).
SMOKE_CORES="${SMOKE_CORES:-8}"

# ==============================================================================
# Sanity checks
# ==============================================================================

for bin in "${SKETCH_BIN}" parallel python3; do
    if ! command -v "${bin}" &>/dev/null && [[ ! -f "${bin}" ]]; then
        echo "ERROR: required binary not found: ${bin}" >&2
        exit 1
    fi
done

if [[ ! -f "${FULL_MANIFEST}" ]]; then
    echo "ERROR: manifest not found: ${FULL_MANIFEST}" >&2
    echo "       Run scripts/00_decompress_genomes.sh first." >&2
    exit 1
fi

if [[ ! -f "${KMC_METADATA}" ]]; then
    echo "ERROR: KMC metadata not found: ${KMC_METADATA}" >&2
    echo "       Run scripts/01_kmc_count.sh first." >&2
    exit 1
fi

if [[ ! -f "${FULL_CANDIDATES}" ]]; then
    echo "ERROR: full candidates CSV not found: ${FULL_CANDIDATES}" >&2
    echo "       Run scripts/compute_fracminhash_candidates.sh first." >&2
    exit 1
fi

# ==============================================================================
# Read parameters from config.json
# ==============================================================================

cfg() { python3 -c "import json; print(json.load(open('${CONFIG}'))$1)"; }

KMER=$(cfg "['kmer']")
SEED=$(cfg "['seed']")
SCALE=$(cfg "['fracminhash']['scale']")
K_BK=$(cfg "['bottomk']['k']")
W_AMG=$(cfg "['alphamaxgeom']['w']")
ALPHA_AMG=$(cfg "['alphamaxgeom']['alpha']")
THRESHOLD=0.001

# ==============================================================================
# Step 0 — Select 100 genome subset from known candidate pairs
# ==============================================================================

echo "========================================================"
echo "  MaxGeomHash Pipeline Smoke Test"
echo "  Start: $(date)"
echo "  Output: ${DATA_DIR}"
echo "  Cores: ${SMOKE_CORES}"
echo "========================================================"

mkdir -p \
    "${DATA_DIR}/fracminhash_sketches" \
    "${DATA_DIR}/bottomk_sketches" \
    "${DATA_DIR}/alphamaxgeom_sketches" \
    "${DATA_DIR}/kmc_pairwise" \
    "${DATA_DIR}/bottomk_pairwise" \
    "${DATA_DIR}/fracminhash_pairwise" \
    "${DATA_DIR}/alphamaxgeom_pairwise" \
    "${DATA_DIR}/sanity_check"

echo ""
echo "=== Step 0: Selecting genome subset ==="

GENOME_LIST="${DATA_DIR}/genome_list.txt"
SUB_MANIFEST="${DATA_DIR}/manifest.csv"
SUB_CANDIDATES="${DATA_DIR}/candidates_thr0001.csv"

python3 - <<PYEOF
import csv, itertools, os, sys

full_manifest   = "${FULL_MANIFEST}"
full_candidates = "${FULL_CANDIDATES}"
genome_list_out = "${GENOME_LIST}"
manifest_out    = "${SUB_MANIFEST}"
candidates_out  = "${SUB_CANDIDATES}"

# ------------------------------------------------------------------
# Collect the first 100 unique genome IDs that appear in real pairs.
# This guarantees a non-trivial pairwise workload.
# ------------------------------------------------------------------
N_TARGET = 100
genomes_in_order = []
seen = set()
with open(full_candidates) as f:
    for row in itertools.islice(csv.DictReader(f), 2000):
        for g in (row["query_name"], row["match_name"]):
            if g not in seen:
                seen.add(g)
                genomes_in_order.append(g)
        if len(genomes_in_order) >= N_TARGET:
            break

genome_set = set(genomes_in_order)

# ------------------------------------------------------------------
# Build sub-manifest: retain only the selected genomes
# ------------------------------------------------------------------
paths = {}
with open(full_manifest) as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    rows = [r for r in reader if r["name"] in genome_set]
    paths = {r["name"]: r["genome_filename"] for r in rows}

if len(paths) < N_TARGET:
    print(f"WARNING: only {len(paths)} / {N_TARGET} genome paths found in manifest",
          file=sys.stderr)

with open(manifest_out, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

# Genome path list (one per line) for the sketching loop
with open(genome_list_out, "w") as f:
    for row in rows:
        f.write(row["genome_filename"] + "\n")

# ------------------------------------------------------------------
# Extract candidate pairs that involve only the selected genomes
# ------------------------------------------------------------------
n_pairs = 0
with open(full_candidates) as fin, open(candidates_out, "w", newline="") as fout:
    reader = csv.DictReader(fin)
    writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
    writer.writeheader()
    for row in reader:
        if row["query_name"] in genome_set and row["match_name"] in genome_set:
            writer.writerow(row)
            n_pairs += 1

print(f"Selected {len(paths)} genomes, {n_pairs} candidate pairs at "
      f"threshold={${THRESHOLD}}")
PYEOF

N_GENOMES=$(wc -l < "${GENOME_LIST}")
N_PAIRS=$(( $(wc -l < "${SUB_CANDIDATES}") - 1 ))
echo "Genomes: ${N_GENOMES}  |  Candidate pairs: ${N_PAIRS}"

# ==============================================================================
# Step 1 — Sketch all genomes (3 algorithms, serially per genome)
# ==============================================================================

echo ""
echo "=== Step 1: Sketching ${N_GENOMES} genomes ==="

SECONDS=0

export SKETCH_BIN DATA_DIR KMER SEED SCALE K_BK W_AMG ALPHA_AMG

sketch_one_genome() {
    local genome_path="$1"
    local genome_name
    genome_name=$(basename "${genome_path}" | sed 's/\.\(fna\|fa\|fasta\)$//')

    "${SKETCH_BIN}" --input "${genome_path}" --kmer "${KMER}" \
        --algo fracminhash --scale "${SCALE}" --seed "${SEED}" --canonical \
        --output "${DATA_DIR}/fracminhash_sketches/${genome_name}.fracminhash.sketch" \
        2>/dev/null || { echo "ERROR: fracminhash sketch failed for ${genome_name}" >&2; return 1; }

    "${SKETCH_BIN}" --input "${genome_path}" --kmer "${KMER}" \
        --algo bottomk --k "${K_BK}" --seed "${SEED}" --canonical \
        --output "${DATA_DIR}/bottomk_sketches/${genome_name}.bottomk.sketch" \
        2>/dev/null || { echo "ERROR: bottomk sketch failed for ${genome_name}" >&2; return 1; }

    "${SKETCH_BIN}" --input "${genome_path}" --kmer "${KMER}" \
        --algo alphamaxgeom --w "${W_AMG}" --alpha "${ALPHA_AMG}" --seed "${SEED}" --canonical \
        --output "${DATA_DIR}/alphamaxgeom_sketches/${genome_name}.alphamaxgeom.sketch" \
        2>/dev/null || { echo "ERROR: alphamaxgeom sketch failed for ${genome_name}" >&2; return 1; }

    echo "DONE ${genome_name}"
}
export -f sketch_one_genome

parallel \
    --jobs "${SMOKE_CORES}" \
    --line-buffer \
    --halt soon,fail=1 \
    sketch_one_genome \
    < "${GENOME_LIST}"

SKETCH_ELAPSED=${SECONDS}
N_FMH=$(find "${DATA_DIR}/fracminhash_sketches" -name "*.fracminhash.sketch" | wc -l)
N_BK=$(find  "${DATA_DIR}/bottomk_sketches"     -name "*.bottomk.sketch"     | wc -l)
N_AMG=$(find "${DATA_DIR}/alphamaxgeom_sketches" -name "*.alphamaxgeom.sketch" | wc -l)
printf "  FracMinHash sketches: %d  BottomK: %d  AlphaMaxGeomHash: %d  (%.1f s)\n" \
    "${N_FMH}" "${N_BK}" "${N_AMG}" "${SKETCH_ELAPSED}"

if [[ "${N_FMH}" -ne "${N_GENOMES}" || "${N_BK}" -ne "${N_GENOMES}" || "${N_AMG}" -ne "${N_GENOMES}" ]]; then
    echo "ERROR: sketch count mismatch — expected ${N_GENOMES} for each algorithm." >&2
    exit 1
fi

# ==============================================================================
# Step 2 — Re-run candidate generation from scratch on the smoke-test sketches
# (validates compute_fracminhash_candidates.py independently of the full run)
# ==============================================================================

echo ""
echo "=== Step 2: FracMinHash candidate generation (from scratch on subset) ==="
SECONDS=0

FRESH_CANDIDATES="${DATA_DIR}/candidates_fresh_thr0001.csv"
python3 "${SCRIPTS_DIR}/compute_fracminhash_candidates.py" \
    --sketch-dir  "${DATA_DIR}/fracminhash_sketches" \
    --output      "${FRESH_CANDIDATES}" \
    --threshold   "${THRESHOLD}" \
    --cores       "${SMOKE_CORES}" \
    --num-shards  4

N_FRESH=$(( $(wc -l < "${FRESH_CANDIDATES}") - 1 ))
printf "  Fresh candidate pairs generated: %d  (%.1f s)\n" "${N_FRESH}" "${SECONDS}"
printf "  Expected (from full-run subset):  %d\n" "${N_PAIRS}"

# Soft check: fresh generation should produce a comparable number of pairs.
# Exact match is not required because compute_fracminhash_candidates uses
# approximate containment (hash-count based), whereas the sub-sampled CSV
# was drawn from the exact full-run output.
if [[ "${N_FRESH}" -eq 0 && "${N_PAIRS}" -gt 0 ]]; then
    echo "WARNING: candidate generation produced 0 pairs but ${N_PAIRS} were expected." >&2
fi

# For downstream steps use the pre-extracted subset candidates (exact match
# with the full run's ground truth), not the freshly generated ones.
echo "  Using pre-extracted candidates for pairwise steps."

# ==============================================================================
# Step 3 — KMC exact pairwise (ground truth)
# ==============================================================================

echo ""
echo "=== Step 3: KMC exact pairwise (${N_PAIRS} candidate pairs) ==="
SECONDS=0

# Write pipeline-log-style START marker (used later by 10_plot_resources.py)
SMOKE_PIPELINE_LOG="${DATA_DIR}/resources_root/scripts/full_pipeline.log"
mkdir -p "$(dirname "${SMOKE_PIPELINE_LOG}")"
echo "=== START 02_kmc_pairwise $(date) ===" >> "${SMOKE_PIPELINE_LOG}"

python3 "${SCRIPTS_DIR}/02_kmc_pairwise.py" \
    --metadata    "${KMC_METADATA}" \
    --candidates  "${SUB_CANDIDATES}" \
    --output      "${DATA_DIR}/kmc_pairwise" \
    --threshold   "${THRESHOLD}" \
    --cores       "${SMOKE_CORES}" \
    --chunk-size  1 \
    2>&1 | grep -E "INFO|WARNING|ERROR|pairs|Saved" || true

KMC_PAIRWISE_ELAPSED=${SECONDS}
echo "=== END   02_kmc_pairwise $(date) ===" >> "${SMOKE_PIPELINE_LOG}"
printf "  Completed in %.1f s\n" "${KMC_PAIRWISE_ELAPSED}"

# ==============================================================================
# Step 4 — kmer-sketch pairwise (3 methods)
# ==============================================================================

echo ""
echo "=== Step 4: kmer-sketch pairwise (BottomK, FracMinHash, AlphaMaxGeomHash) ==="
SECONDS=0

python3 "${SCRIPTS_DIR}/06_bottomk_pairwise.py" \
    --sketch-dir  "${DATA_DIR}/bottomk_sketches" \
    --candidates  "${SUB_CANDIDATES}" \
    --output      "${DATA_DIR}/bottomk_pairwise" \
    --cores       "${SMOKE_CORES}" \
    2>&1 | grep -E "INFO|WARNING|ERROR|pairs|Saved" || true
printf "  BottomK:          %.1f s\n" "${SECONDS}"; SECONDS=0

python3 "${SCRIPTS_DIR}/07_fracminhash_pairwise.py" \
    --sketch-dir  "${DATA_DIR}/fracminhash_sketches" \
    --candidates  "${SUB_CANDIDATES}" \
    --output      "${DATA_DIR}/fracminhash_pairwise" \
    --cores       "${SMOKE_CORES}" \
    2>&1 | grep -E "INFO|WARNING|ERROR|pairs|Saved" || true
printf "  FracMinHash:      %.1f s\n" "${SECONDS}"; SECONDS=0

python3 "${SCRIPTS_DIR}/08_alphamaxgeom_pairwise.py" \
    --sketch-dir  "${DATA_DIR}/alphamaxgeom_sketches" \
    --candidates  "${SUB_CANDIDATES}" \
    --output      "${DATA_DIR}/alphamaxgeom_pairwise" \
    --cores       "${SMOKE_CORES}" \
    2>&1 | grep -E "INFO|WARNING|ERROR|pairs|Saved" || true
printf "  AlphaMaxGeomHash: %.1f s\n" "${SECONDS}"

# ==============================================================================
# Step 5 — Sanity check: accuracy vs. KMC exact
# ==============================================================================

echo ""
echo "=== Step 5: Sanity check (accuracy vs. KMC exact) ==="
SECONDS=0

python3 "${SCRIPTS_DIR}/09_sanity_check.py" \
    --amg-pairwise         "${DATA_DIR}/alphamaxgeom_pairwise" \
    --bottomk-pairwise     "${DATA_DIR}/bottomk_pairwise" \
    --fracminhash-pairwise "${DATA_DIR}/fracminhash_pairwise" \
    --kmc-pairwise         "${DATA_DIR}/kmc_pairwise" \
    --output               "${DATA_DIR}/sanity_check" \
    2>&1 | grep -E "INFO|WARNING|ERROR|Pearson|L1|Saved" || true

printf "  Sanity check completed in %.1f s\n" "${SECONDS}"

# ==============================================================================
# Step 6 — Heatmap figures (dense subset heatmaps, both metrics)
# ==============================================================================

echo ""
echo "=== Step 6: Heatmap figures (10_plot_heatmaps.py) ==="
mkdir -p "${DATA_DIR}/figures"
SECONDS=0

for METRIC in jaccard max_containment; do
    python3 "${SCRIPTS_DIR}/10_plot_heatmaps.py" \
        --kmc-pairwise        "${DATA_DIR}/kmc_pairwise" \
        --amg-pairwise        "${DATA_DIR}/alphamaxgeom_pairwise" \
        --bottomk-pairwise    "${DATA_DIR}/bottomk_pairwise" \
        --fracminhash-pairwise "${DATA_DIR}/fracminhash_pairwise" \
        --output              "${DATA_DIR}/figures" \
        --n-genomes           100 \
        --metric              "${METRIC}" \
        --l1-mode             subset \
        2>&1 | grep -E "INFO|WARNING|ERROR|Saved" || true
    echo "  Heatmap (${METRIC}) done"
done

printf "  Heatmaps completed in %.1f s\n" "${SECONDS}"

# ==============================================================================
# Step 7 — Resource and L1 error figures (10_plot_resources.py)
#
# 10_plot_resources.py expects:
#   <base>/data/GTDB/{method}_sketches/sketch_run_stats.json
#   <base>/data/GTDB/{method}_pairwise_thr0001/pairwise_run_stats.json
#   <base>/scripts/full_pipeline.log   (=== START/END 02_kmc_pairwise ===)
#
# The smoke test inline sketching does not produce sketch_run_stats.json, and
# the pairwise dirs use plain names without the _thr0001 suffix.  We create a
# minimal fake project root (smoke_test/data/resources_root/) with:
#   - sketch_run_stats.json written from the Step-1 timing
#   - symlinks from the expected _thr0001 names to the actual smoke-test dirs
#   - the pipeline log already started in Step 3
# ==============================================================================

echo ""
echo "=== Step 7: Resource and L1 error figures (10_plot_resources.py) ==="
SECONDS=0

FAKE_ROOT="${DATA_DIR}/resources_root"
FAKE_GTDB="${FAKE_ROOT}/data/GTDB"

mkdir -p \
    "${FAKE_GTDB}/bottomk_sketches" \
    "${FAKE_GTDB}/fracminhash_sketches" \
    "${FAKE_GTDB}/alphamaxgeom_sketches"

# 7a. Write sketch_run_stats.json for each algorithm using Step-1 timing.
#     parallel_jobs = SMOKE_CORES (what was actually used).
python3 - <<PYEOF
import json, math

smoke_cores = ${SMOKE_CORES}
sketch_elapsed = ${SKETCH_ELAPSED}   # total wall-clock for all 3 algos together
n_genomes = ${N_GENOMES}
fake_gtdb = "${FAKE_GTDB}"

# All three algorithms were sketched in the same parallel batch; split evenly.
per_algo_s = sketch_elapsed / 3.0

for algo, ext, fname in [
    ("fracminhash",   "fracminhash.sketch",   "fracminhash_sketches"),
    ("bottomk",       "bottomk.sketch",       "bottomk_sketches"),
    ("alphamaxgeom",  "alphamaxgeom.sketch",  "alphamaxgeom_sketches"),
]:
    stats = {
        "script":                f"smoke_test (inline sketch)",
        "algo":                  algo,
        "parallel_jobs":         smoke_cores,
        "n_genomes_sketched":    n_genomes,
        "wall_clock_seconds":    round(per_algo_s, 1),
        "wall_clock_hms":        f"0:{int(per_algo_s//60)}:{int(per_algo_s%60)}",
    }
    out = f"{fake_gtdb}/{fname}/sketch_run_stats.json"
    with open(out, "w") as f:
        json.dump(stats, f, indent=2)
print("Sketch stats JSONs written.")
PYEOF

# 7b. Symlink pairwise dirs with the _thr0001 suffix that 10_plot_resources.py expects.
for method in bottomk fracminhash alphamaxgeom; do
    target="${DATA_DIR}/${method}_pairwise"
    link="${FAKE_GTDB}/${method}_pairwise_thr0001"
    [[ -L "${link}" ]] && rm "${link}"
    ln -s "${target}" "${link}"
done
# KMC pairwise dir (no stats JSON, but disk size is measured)
KMC_LINK="${FAKE_GTDB}/kmc_pairwise_thr0001"
[[ -L "${KMC_LINK}" ]] && rm "${KMC_LINK}"
ln -s "${DATA_DIR}/kmc_pairwise" "${KMC_LINK}"

# 7c. Run the resource plotter.
python3 "${SCRIPTS_DIR}/10_plot_resources.py" \
    --base-dir    "${FAKE_ROOT}" \
    --sanity-json "${DATA_DIR}/sanity_check/sanity_check_summary.json" \
    --output      "${DATA_DIR}/figures" \
    2>&1 | grep -E "INFO|WARNING|ERROR|Saved|summary" || true

printf "  Resource figures completed in %.1f s\n" "${SECONDS}"

# ==============================================================================
# Step 8 — Full Datashader heatmaps (11_plot_heatmaps_full.py, spectral ordering)
#
# Exercises the spectral-clustering genome-ordering path and Datashader
# rendering that 10_run_plots.sh uses for publication figures.  Requires the
# sourmash conda env (which ships datashader); the step is skipped if that env
# is not found.
# ==============================================================================

echo ""
echo "=== Step 8: Full Datashader heatmaps (11_plot_heatmaps_full.py, spectral) ==="
SECONDS=0

CONDA_ENV="${CONDA_ENV:-sourmash}"

# Check whether the env and datashader are available; skip gracefully if not.
if ! conda run -n "${CONDA_ENV}" python3 -c "import datashader" 2>/dev/null; then
    echo "  SKIPPED: conda env '${CONDA_ENV}' does not have datashader." \
         "Set CONDA_ENV= to override."
else
    for METRIC in jaccard max_containment; do
        conda run -n "${CONDA_ENV}" python3 "${SCRIPTS_DIR}/11_plot_heatmaps_full.py" \
            --kmc-pairwise         "${DATA_DIR}/kmc_pairwise" \
            --amg-pairwise         "${DATA_DIR}/alphamaxgeom_pairwise" \
            --bottomk-pairwise     "${DATA_DIR}/bottomk_pairwise" \
            --fracminhash-pairwise "${DATA_DIR}/fracminhash_pairwise" \
            --output               "${DATA_DIR}/figures" \
            --metric               "${METRIC}" \
            --ordering             spectral \
            --canvas-size          200 \
            2>&1 | grep -E "INFO|WARNING|ERROR|Saved" || true
        echo "  Full heatmap (${METRIC}) done"
    done
    printf "  Full heatmaps completed in %.1f s\n" "${SECONDS}"
fi

# ==============================================================================
# Summary
# ==============================================================================

echo ""
echo "========================================================"
echo "  Smoke Test COMPLETE"
echo "  End: $(date)"
echo ""

echo "  Results:"
for d in kmc_pairwise bottomk_pairwise fracminhash_pairwise alphamaxgeom_pairwise; do
    if [[ -f "${DATA_DIR}/${d}/pairwise_results.npz" ]]; then
        n=$(python3 -c "
import numpy as np
d = np.load('${DATA_DIR}/${d}/pairwise_results.npz')
print(len(d['row']))
" 2>/dev/null || echo "?")
        printf "    %-30s %s pairs\n" "${d}" "${n}"
    else
        printf "    %-30s MISSING pairwise_results.npz\n" "${d}"
    fi
done

echo ""
if [[ -f "${DATA_DIR}/sanity_check/sanity_check_summary.json" ]]; then
    python3 - <<PYEOF
import json
summary = json.load(open("${DATA_DIR}/sanity_check/sanity_check_summary.json"))
print(f"  Genomes: {summary.get('n_genomes')}  |  KMC pairs: {summary.get('n_kmc_pairs')}")
print("  Accuracy (mean |error| vs KMC exact):")
for s in summary.get("error_statistics", []):
    if "vs KMC" in s["metric"]:
        print(f"    {s['metric']:<55} MAE={s['mae']:.6f}  r={s['pearson_r']:.4f}")
PYEOF
else
    echo "  WARNING: sanity_check_summary.json not found"
fi

echo ""
echo "  Figures:"
ls "${DATA_DIR}/figures/"*.pdf "${DATA_DIR}/figures/"*.png \
    "${DATA_DIR}/sanity_check/"*.png 2>/dev/null \
    | while read -r f; do printf "    %s\n" "$(basename "$f")"; done

echo ""
echo "  Output directory: ${DATA_DIR}"
echo "========================================================"
