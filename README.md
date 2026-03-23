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
00_decompress_genomes.sh  ← ONE-TIME SETUP: decompress all .fna.gz → flat .fna dir
                            also writes manysketch_uncompressed.csv

01_kmc_count.sh           ← KMC k-mer counting (ground truth)
02_kmc_pairwise.sh        ← KMC exact pairwise similarity

03_bottomk_sketch.sh      ← BottomK sketching        (reads uncompressed .fna)
04_fracminhash_sketch.sh  ← FracMinHash sketching     (reads uncompressed .fna)
05_alphamaxgeom_sketch.sh ← AlphaMaxGeomHash sketching(reads uncompressed .fna)

06_bottomk_pairwise.sh    ← BottomK pairwise (calls 06_bottomk_pairwise.py)
07_fracminhash_pairwise.sh← FracMinHash pairwise (calls 07_fracminhash_pairwise.py)
08_alphamaxgeom_pairwise.sh← AlphaMaxGeomHash pairwise (calls 08_alphamaxgeom_pairwise.py)

09_sanity_check.py        ← Compare all methods vs KMC; produce accuracy stats
10_run_plots.sh           ← Generate all publication figures
  └── 10_plot_heatmaps.py ← Error heatmaps (relative error vs KMC)
  └── 10_plot_resources.py← Time + disk bar charts + accuracy vs. resources scatter
```

Steps 01–02 (KMC) are prerequisites for accuracy comparisons and have already
been run.  **Step 00 must be run once** before any sketch script.
Steps 03–10 can then be run independently for each method.

---

## Prerequisites

```bash
# Compile the kmer-sketch binary (if not already done):
cd scripts/kmer-sketch && make && cd ../..

# Python dependencies:
pip install -r scripts/requirements.txt

# ONE-TIME: decompress all 143,614 genomes to a flat .fna directory.
# This eliminates per-job gunzip overhead in the sketch scripts,
# allowing the sketch binary to read files directly and fully saturate CPUs.
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
    --bottomk-pairwise data/GTDB/bottomk_pairwise \
    --kmc-pairwise     data/GTDB/kmc_pairwise \
    --output           data/GTDB/bottomk_pairwise/sanity_check

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
    --amg-pairwise         data/GTDB/alphamaxgeom_pairwise \
    --bottomk-pairwise     data/GTDB/bottomk_pairwise \
    --fracminhash-pairwise data/GTDB/fracminhash_pairwise \
    --kmc-pairwise         data/GTDB/kmc_pairwise \
    --output               data/GTDB/sanity_check
```

To also include the Sourmash FracMinHash baseline:

```bash
python3 scripts/09_sanity_check.py ... \
    --include-sourmash \
    --sourmash-csv data/GTDB/gtdb_pairwise_containment.csv
```

---

## Generating publication figures

```bash
bash scripts/10_run_plots.sh
```

This produces all figures in `data/GTDB/figures/`.  Environment variables:

| Variable | Default | Effect |
|----------|---------|--------|
| `N_GENOMES` | 500 | Number of genomes in heatmap panels |
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

### KMC (ground truth) — already run

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `01_kmc_count.sh` | GTDB genome FASTA files | `data/GTDB/kmc_dbs/` | Once; already done |
| `02_kmc_pairwise.sh` | kmc_dbs | `data/GTDB/kmc_pairwise/` | Once; already done |

### BottomK (kmer-sketch)

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `03_bottomk_sketch.sh` | `gtdb_genomes_reps_r226_uncompressed/*.fna` | `data/GTDB/bottomk_sketches/` | After `00_decompress_genomes.sh` |
| `06_bottomk_pairwise.py/.sh` | bottomk_sketches + candidates CSV | `data/GTDB/bottomk_pairwise/` | After `03_bottomk_sketch.sh` |

Parameters: `--algo bottomk --k 1000 --kmer 31 --canonical --seed 42`

### FracMinHash (kmer-sketch)

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `04_fracminhash_sketch.sh` | `gtdb_genomes_reps_r226_uncompressed/*.fna` | `data/GTDB/fracminhash_sketches/` | After `00_decompress_genomes.sh` |
| `07_fracminhash_pairwise.py/.sh` | fracminhash_sketches + candidates CSV | `data/GTDB/fracminhash_pairwise/` | After `04_fracminhash_sketch.sh` |

Parameters: `--algo fracminhash --scale 0.01 --kmer 31 --canonical --seed 42`

### AlphaMaxGeomHash

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `05_alphamaxgeom_sketch.sh` | `gtdb_genomes_reps_r226_uncompressed/*.fna` | `data/GTDB/alphamaxgeom_sketches/` | After `00_decompress_genomes.sh` |
| `08_alphamaxgeom_pairwise.py/.sh` | alphamaxgeom_sketches + candidates CSV | `data/GTDB/alphamaxgeom_pairwise/` | After `05_alphamaxgeom_sketch.sh` |

Parameters: `--algo alphamaxgeom --w 64 --alpha 0.45 --kmer 31 --canonical --seed 42`

### Analysis and visualization

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `09_sanity_check.py` | Any combination of pairwise dirs + KMC | `sanity_check_summary.json` + scatter PNGs | After pairwise step(s) |
| `10_plot_heatmaps.py` | Pairwise dirs + KMC pairwise | Heatmap PDF/PNG | After `09_sanity_check.py` |
| `10_plot_resources.py` | Project root (reads JSON stats + logs) | Resource bar chart PDF/PNG | After pairwise step(s) |
| `10_run_plots.sh` | — | All figures in `data/GTDB/figures/` | Wrapper for both `10_plot_*.py` |

---

## Output data layout

```
data/GTDB/
├── gtdb_genomes_reps_r226/               ← original compressed genomes (.fna.gz)
├── gtdb_genomes_reps_r226_uncompressed/  ← flat dir of plain .fna files (~400-650 GB)
│                                            created by 00_decompress_genomes.sh
├── manysketch_uncompressed.csv           ← manifest pointing to .fna files
├── kmc_dbs/                              ← KMC k-mer databases (3.9 TB)
├── kmc_pairwise/                         ← KMC exact pairwise (41 MB NPZ)
├── bottomk_sketches/                     ← BottomK sketches + sketch_run_stats.json
├── bottomk_pairwise/                     ← BottomK pairwise NPZ + run stats
│   └── sanity_check/
├── fracminhash_sketches/                 ← FracMinHash (kmer-sketch) sketches
├── fracminhash_pairwise/                 ← FracMinHash pairwise NPZ + run stats
├── alphamaxgeom_sketches/                ← AlphaMaxGeomHash sketches
├── alphamaxgeom_pairwise/                ← AlphaMaxGeomHash pairwise NPZ + run stats
│   └── sanity_check/
├── gtdb_pairwise_containment.csv         ← Sourmash FracMinHash pairwise (candidate pairs)
├── sanity_check/                         ← combined sanity check (all methods)
└── figures/                              ← all publication figures (PDF + PNG)
    ├── heatmap_jaccard_relative_error.{pdf,png}
    ├── heatmap_max_containment_relative_error.{pdf,png}
    ├── resources_time.{pdf,png}
    ├── resources_disk.{pdf,png}
    ├── resources_tradeoff.{pdf,png}
    └── resource_summary.json
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
the kmer-sketch implementation.  Its pairwise CSV
(`gtdb_pairwise_containment.csv`) is still used as the **candidate pair
filter** for all kmer-sketch pairwise scripts — only pairs with
`max_containment > 0.01` in the Sourmash run are evaluated.

The Sourmash method is **excluded from default figures**.  To include it,
pass `--include-sourmash` to `09_sanity_check.py`, `10_plot_heatmaps.py`,
`10_plot_resources.py`, or set `INCLUDE_SOURMASH=1` when running
`10_run_plots.sh`.
