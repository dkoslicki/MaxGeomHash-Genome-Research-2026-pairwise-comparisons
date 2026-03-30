# MaxGeomHash Genome Research 2026

Pairwise genome similarity study comparing three sketching algorithms against
exact k-mer counting (KMC) across all 143,614 GTDB r226 representative genomes.

**Methods compared**

| Method | Algorithm | Tool |
|--------|-----------|------|
| KMC (exact) | Exact k-mer counting | kmc + custom pairwise |
| BottomK | BottomK, k=1000, kmer=31 | kmer-sketch binary |
| FracMinHash | FracMinHash, scale=0.01 (= scaled=100), kmer=31 | kmer-sketch binary |
| AlphaMaxGeomHash | AlphaMaxGeomHash, W=64, α=0.45, kmer=31 | kmer-sketch binary |
| Sourmash FracMinHash | FracMinHash, scaled=1000, kmer=31 | sourmash | *(baseline; code preserved, excluded from default figures)* |

---

## Pipeline overview

```
00_decompress_genomes.sh       ← ONE-TIME SETUP: decompress all .fna.gz → flat .fna dir
                                  also writes manysketch_uncompressed.csv

01_kmc_count.sh                ← KMC k-mer counting (ground truth)

make_fracminhash_sketches.sh   ← FracMinHash sketching for candidate generation
compute_fracminhash_candidates.sh ← all-vs-all; writes gtdb_pairwise_containment_thr0001.csv
02_kmc_pairwise.sh             ← KMC exact pairwise on candidates → kmc_pairwise_thr0001/

03_bottomk_sketch.sh           ← BottomK sketching        (reads uncompressed .fna)
04_fracminhash_sketch.sh       ← FracMinHash sketching     (reads uncompressed .fna)
05_alphamaxgeom_sketch.sh      ← AlphaMaxGeomHash sketching(reads uncompressed .fna)

06_bottomk_pairwise.sh         ← BottomK pairwise on candidates
07_fracminhash_pairwise.sh     ← FracMinHash pairwise on candidates
08_alphamaxgeom_pairwise.sh    ← AlphaMaxGeomHash pairwise on candidates

09_sanity_check.py             ← Compare all methods vs KMC; produce accuracy stats
10_run_plots.sh                ← Generate all publication figures
  └── 11_plot_heatmaps_full.py ← Full-dataset Datashader heatmaps (all ~143k genomes)
  └── 10_plot_heatmaps.py      ← Dense subset heatmaps (SUBSET_HEATMAP=1)
  └── 10_plot_resources.py     ← Time + disk bar charts + accuracy vs. resources scatter
```

Steps 00–02 are prerequisites for accuracy comparisons.
**The recommended way to (re-)run everything from step 2 onward (if on the cbdmk02 GPU server)** is:

```bash
bash scripts/run_sketch_and_pairwise_plus_plot.sh
# Logs → scripts/full_pipeline.log
```

Steps 03–10 can also be run independently per method.

---

## Prerequisites

```bash
# Compile the kmer-sketch binary (if not already done):
cd scripts/kmer-sketch && make && cd ../..

# Python dependencies:
pip install -r scripts/requirements.txt

# ONE-TIME: decompress all 143,614 genomes to a flat .fna directory.
# Takes ~10-20 minutes; safe to re-run after partial failures (skips existing).
bash scripts/00_decompress_genomes.sh
```

The decompressed genomes are stored in:
`data/GTDB/gtdb_genomes_reps_r226_uncompressed/`  (~400–650 GB)

A new manifest pointing to these files is written to:
`data/GTDB/manysketch_uncompressed.csv`

All sketch scripts (03, 04, 05) read from this manifest automatically.

---

## Running a new method from scratch

### Test mode (first 200 genomes)

```bash
# BottomK
TEST_N=200 bash scripts/03_bottomk_sketch.sh
bash scripts/06_bottomk_pairwise.sh
python3 scripts/09_sanity_check.py \
    --bottomk-pairwise data/GTDB/bottomk_pairwise_thr0001 \
    --kmc-pairwise     data/GTDB/kmc_pairwise_thr0001 \
    --output           data/GTDB/bottomk_pairwise_thr0001/sanity_check

# FracMinHash
TEST_N=200 bash scripts/04_fracminhash_sketch.sh
bash scripts/07_fracminhash_pairwise.sh

# AlphaMaxGeomHash
TEST_N=200 bash scripts/05_alphamaxgeom_sketch.sh
bash scripts/08_alphamaxgeom_pairwise.sh
```

### Full run (all 143,614 genomes)

Identical commands — just omit `TEST_N`.  The pairwise scripts automatically
use however many sketches are present.

---

## Running all three kmer-sketch methods in serial

```bash
cat > /tmp/run_pipeline.sh << 'EOF'
BASE="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"
LOG="$BASE/scripts/kmer_sketch_full_run.log"
for s in 05_alphamaxgeom_sketch.sh 03_bottomk_sketch.sh 04_fracminhash_sketch.sh \
          08_alphamaxgeom_pairwise.sh 06_bottomk_pairwise.sh 07_fracminhash_pairwise.sh; do
    echo ""
    echo "=== START $s $(date) ==="
    bash "$BASE/scripts/$s" && echo "=== END $s $(date) ===" || { echo "FAILED: $s"; exit 1; }
done
EOF
nohup bash /tmp/run_pipeline.sh \
    > /scratch/shared_data/MaxGeomHash_Genome_Research_2026/scripts/kmer_sketch_full_run.log 2>&1 \
    & disown $!; echo "PID: $!"
```

---

## Sanity check (after all methods are run)

```bash
python3 scripts/09_sanity_check.py \
    --amg-pairwise         data/GTDB/alphamaxgeom_pairwise_thr0001 \
    --bottomk-pairwise     data/GTDB/bottomk_pairwise_thr0001 \
    --fracminhash-pairwise data/GTDB/fracminhash_pairwise_thr0001 \
    --kmc-pairwise         data/GTDB/kmc_pairwise_thr0001 \
    --output               data/GTDB/alphamaxgeom_pairwise_thr0001/sanity_check
```

To also include the Sourmash FracMinHash baseline:

```bash
python3 scripts/09_sanity_check.py ... \
    --include-sourmash \
    --sourmash-csv data/GTDB/gtdb_pairwise_containment_thr0001.csv
```

---

## Generating publication figures

```bash
bash scripts/10_run_plots.sh
```

This produces all figures in `data/GTDB/figures_thr0001/`.  Environment variables:

| Variable | Default | Effect |
|----------|---------|--------|
| `SUBSET_HEATMAP` | 0 | Set to 1 to also produce dense N-genome subset heatmaps |
| `N_GENOMES` | 500 | Number of genomes in subset heatmap panels |
| `ORDERING` | spectral | Genome ordering for full heatmaps: `degree`, `rcm`, or `spectral` |
| `SHOW_REF` | 0 | Set to 1 to add a KMC reference panel to heatmaps |
| `INCLUDE_SOURMASH` | 0 | Set to 1 to add Sourmash FracMinHash to figures |

The script automatically skips methods whose output directories do not yet
exist, so it works correctly even if only some methods have been run.

---

## Script reference

### One-time setup

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `00_decompress_genomes.sh` | `gtdb_genomes_reps_r226/*.fna.gz` | `gtdb_genomes_reps_r226_uncompressed/*.fna` + `manysketch_uncompressed.csv` | **Once**, before any sketch script |

### KMC (ground truth) (already run on cbdmk02; takes a long time)

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `01_kmc_count.sh` | GTDB genome FASTA files | `data/GTDB/kmc_dbs/` | Once; already done on cbdmk02 as it takes a while |
| `make_fracminhash_sketches.sh` | `manysketch_uncompressed.csv` | `data/GTDB/fracminhash_sketches/` | Once (skips if done) |
| `compute_fracminhash_candidates.sh` | fracminhash_sketches | `data/GTDB/gtdb_pairwise_containment_thr0001.csv` | Once|
| `02_kmc_pairwise.sh` | kmc_dbs + candidates CSV | `data/GTDB/kmc_pairwise_thr0001/` | Once|

### BottomK (kmer-sketch)

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `03_bottomk_sketch.sh` | `gtdb_genomes_reps_r226_uncompressed/*.fna` | `data/GTDB/bottomk_sketches/` | After `00_decompress_genomes.sh` |
| `06_bottomk_pairwise.py/.sh` | bottomk_sketches + candidates CSV | `data/GTDB/bottomk_pairwise_thr0001/` | After `03_bottomk_sketch.sh` |

Parameters: `--algo bottomk --k 1000 --kmer 31 --canonical --seed 42`

### FracMinHash (kmer-sketch)

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `04_fracminhash_sketch.sh` | `gtdb_genomes_reps_r226_uncompressed/*.fna` | `data/GTDB/fracminhash_sketches/` | After `00_decompress_genomes.sh` |
| `07_fracminhash_pairwise.py/.sh` | fracminhash_sketches + candidates CSV | `data/GTDB/fracminhash_pairwise_thr0001/` | After `04_fracminhash_sketch.sh` |

Parameters: `--algo fracminhash --scale 0.01 --kmer 31 --canonical --seed 42`

### AlphaMaxGeomHash

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `05_alphamaxgeom_sketch.sh` | `gtdb_genomes_reps_r226_uncompressed/*.fna` | `data/GTDB/alphamaxgeom_sketches/` | After `00_decompress_genomes.sh` |
| `08_alphamaxgeom_pairwise.py/.sh` | alphamaxgeom_sketches + candidates CSV | `data/GTDB/alphamaxgeom_pairwise_thr0001/` | After `05_alphamaxgeom_sketch.sh` |

Parameters: `--algo alphamaxgeom --w 64 --alpha 0.45 --kmer 31 --canonical --seed 42`

### Analysis and visualization

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `09_sanity_check.py` | Any combination of pairwise dirs + KMC | `sanity_check_summary.json` + scatter PNGs | After pairwise step(s) |
| `10_plot_heatmaps.py` | Pairwise dirs + KMC pairwise | Heatmap PDF/PNG (subset) | After `09_sanity_check.py` |
| `11_plot_heatmaps_full.py` | Pairwise dirs + KMC pairwise | Full-dataset Datashader heatmaps | After pairwise step(s) |
| `10_plot_resources.py` | Project root (reads JSON stats + logs) | Resource bar chart PDF/PNG | After pairwise step(s) |
| `10_run_plots.sh` | — | All figures in `data/GTDB/figures_thr0001/` | Wrapper for all `10_plot_*.py` + `11_plot_heatmaps_full.py` |

---

## Output data layout

```
data/GTDB/
├── gtdb_genomes_reps_r226/               ← original compressed genomes (.fna.gz)
├── gtdb_genomes_reps_r226_uncompressed/  ← flat dir of plain .fna files (~400-650 GB)
│                                            created by 00_decompress_genomes.sh
├── manysketch_uncompressed.csv           ← manifest pointing to .fna files
├── kmc_dbs/                              ← KMC k-mer databases (3.9 TB)
├── kmc_pairwise_thr0001/                 ← KMC exact pairwise (NPZ)
├── gtdb_pairwise_containment_thr0001.csv ← FracMinHash candidate pairs (threshold=0.001)
├── bottomk_sketches/                     ← BottomK sketches + sketch_run_stats.json
├── bottomk_pairwise_thr0001/             ← BottomK pairwise NPZ + run stats
│   └── sanity_check/
├── fracminhash_sketches/                 ← FracMinHash (kmer-sketch) sketches
├── fracminhash_pairwise_thr0001/         ← FracMinHash pairwise NPZ + run stats
├── alphamaxgeom_sketches/                ← AlphaMaxGeomHash sketches
├── alphamaxgeom_pairwise_thr0001/        ← AlphaMaxGeomHash pairwise NPZ + run stats
│   └── sanity_check/
└── figures_thr0001/                      ← all publication figures (PDF + PNG)
    ├── heatmap_jaccard_relative_error_full_{ordering}.{pdf,png}
    ├── heatmap_max_containment_relative_error_full_{ordering}.{pdf,png}
    ├── resources_time.{pdf,png}
    ├── resources_disk.{pdf,png}
    └── resources_tradeoff.{pdf,png}
```

---

## Notes on script numbering and archive

The pipeline was extended after the initial AMG-only run.  MinHash was later
replaced by BottomK.  Superseded scripts were moved to `scripts/archive/`
(still tracked by git for full reproducibility):

| Archived script | Superseded by | Key difference |
|-----------------|---------------|----------------|
| `archive/03_alphamaxgeom_sketch.sh` | `05_alphamaxgeom_sketch.sh` | Same logic; renumbered |
| `archive/04_alphamaxgeom_pairwise.*` | `08_alphamaxgeom_pairwise.*` | Same logic; renumbered |
| `archive/05_sanity_check.py` | `09_sanity_check.py` | AMG+Sourmash only → multi-method |
| `archive/06_plot_heatmaps.py` | `10_plot_heatmaps.py` | Hardcoded → any kmer-sketch method |
| `archive/06_plot_resources.py` | `10_plot_resources.py` | Hardcoded → all four methods |
| `archive/06_run_plots.sh` | `10_run_plots.sh` | Calls archived scripts → current scripts |

Use the **active scripts** (00–10 in `scripts/`) for all current and future work.

---

## Notes on the Sourmash FracMinHash baseline

The Sourmash FracMinHash pipeline was run as an independent baseline before
the kmer-sketch implementation.  Its pairwise CSV is preserved for reference.

Candidate pair generation now uses `compute_fracminhash_candidates.py` (the
kmer-sketch FracMinHash all-vs-all), which replaced sourmash manysketch after
sourmash was found to corrupt `.sig.zip` files at scaled=100.  The candidates
CSV (`gtdb_pairwise_containment_thr0001.csv`, threshold=0.001) is used as the
input filter for all kmer-sketch pairwise scripts.

The Sourmash method is **excluded from default figures**.  To include it,
pass `--include-sourmash` to `09_sanity_check.py`, `10_plot_heatmaps.py`,
`10_plot_resources.py`, or set `INCLUDE_SOURMASH=1` when running
`10_run_plots.sh`.
