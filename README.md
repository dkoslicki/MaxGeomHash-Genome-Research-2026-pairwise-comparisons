# MaxGeomHash Genome Research 2026

Pairwise genome similarity study comparing three sketching algorithms against
exact k-mer counting (KMC) across all 143,614 GTDB r226 representative genomes.

**Methods compared**

| Method | Algorithm | Tool |
|--------|-----------|------|
| KMC (exact) | Exact k-mer counting | kmc + custom pairwise |
| MinHash | MinHash, num-perm=1000, k=31 | kmer-sketch binary |
| FracMinHash | FracMinHash, scale=0.001 (= scaled=1000), k=31 | kmer-sketch binary |
| AlphaMaxGeomHash | AlphaMaxGeomHash, W=64, α=0.45, k=31 | kmer-sketch binary |
| Sourmash FracMinHash | FracMinHash, scaled=1000, k=31 | sourmash | *(baseline; code preserved, excluded from default figures)* |

---

## Pipeline overview

```
01_kmc_count.sh           ← KMC k-mer counting (ground truth)
02_kmc_pairwise.sh        ← KMC exact pairwise similarity

03_minhash_sketch.sh      ← MinHash sketching
04_fracminhash_sketch.sh  ← FracMinHash (kmer-sketch) sketching
05_alphamaxgeom_sketch.sh ← AlphaMaxGeomHash sketching

06_minhash_pairwise.sh    ← MinHash pairwise (calls 06_minhash_pairwise.py)
07_fracminhash_pairwise.sh← FracMinHash pairwise (calls 07_fracminhash_pairwise.py)
08_alphamaxgeom_pairwise.sh← AlphaMaxGeomHash pairwise (calls 08_alphamaxgeom_pairwise.py)

09_sanity_check.py        ← Compare all methods vs KMC; produce accuracy stats
10_run_plots.sh           ← Generate all publication figures
  └── 10_plot_heatmaps.py ← Error heatmaps (relative error vs KMC)
  └── 10_plot_resources.py← Time + disk bar charts + accuracy vs. resources scatter
```

Steps 01–02 (KMC) are prerequisites for accuracy comparisons and have already
been run.  Steps 03–10 can be run independently for each method.

---

## Prerequisites

```bash
# Compile the kmer-sketch binary (if not already done):
cd scripts/kmer-sketch && make && cd ../..

# Python dependencies:
pip install -r scripts/requirements.txt
```

---

## Running a new method from scratch

All three kmer-sketch method pipelines follow the same pattern.
Replace `METHOD` with `minhash`, `fracminhash`, or `alphamaxgeom`:

### Test mode (first 200 genomes)

```bash
# Sketch
TEST_N=200 bash scripts/03_minhash_sketch.sh

# Pairwise
bash scripts/06_minhash_pairwise.sh

# Sanity check (compare vs KMC)
python3 scripts/09_sanity_check.py \
    --minhash-pairwise data/GTDB/minhash_pairwise \
    --kmc-pairwise     data/GTDB/kmc_pairwise \
    --output           data/GTDB/minhash_pairwise/sanity_check
```

### Full run (all 143,614 genomes)

Identical commands — just omit `TEST_N`:

```bash
bash scripts/03_minhash_sketch.sh
bash scripts/06_minhash_pairwise.sh
```

The pairwise scripts automatically use however many sketches are present, so
test mode and full-run mode use the exact same command.

---

## Running all three kmer-sketch methods in parallel (full run)

```bash
nohup bash -c '
  bash scripts/03_minhash_sketch.sh      > scripts/minhash_sketch.log 2>&1
  bash scripts/06_minhash_pairwise.sh    >> scripts/minhash_sketch.log 2>&1
' &

nohup bash -c '
  bash scripts/04_fracminhash_sketch.sh      > scripts/fracminhash_sketch.log 2>&1
  bash scripts/07_fracminhash_pairwise.sh    >> scripts/fracminhash_sketch.log 2>&1
' &

# (AlphaMaxGeomHash has already been run; re-run only if needed)
# nohup bash -c '
#   bash scripts/05_alphamaxgeom_sketch.sh    > scripts/alphamaxgeom_full.log 2>&1
#   bash scripts/08_alphamaxgeom_pairwise.sh >> scripts/alphamaxgeom_full.log 2>&1
# ' &
```

---

## Sanity check (after all methods are run)

```bash
python3 scripts/09_sanity_check.py \
    --amg-pairwise         data/GTDB/alphamaxgeom_pairwise \
    --minhash-pairwise     data/GTDB/minhash_pairwise \
    --fracminhash-pairwise data/GTDB/fracminhash_pairwise \
    --kmc-pairwise         data/GTDB/kmc_pairwise \
    --output               data/GTDB/sanity_check
```

To also include the Sourmash FracMinHash baseline in the output:

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

### KMC (ground truth) — already run

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `01_kmc_count.sh` | GTDB genome FASTA files | `data/GTDB/kmc_dbs/` | Once; already done |
| `02_kmc_pairwise.sh` | kmc_dbs | `data/GTDB/kmc_pairwise/` | Once; already done |

### MinHash (kmer-sketch)

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `03_minhash_sketch.sh` | GTDB genome FASTA.gz files | `data/GTDB/minhash_sketches/` | Before `06_minhash_pairwise.sh` |
| `06_minhash_pairwise.py/.sh` | minhash_sketches + candidates CSV | `data/GTDB/minhash_pairwise/` | After `03_minhash_sketch.sh` |

Parameters: `--algo minhash --num-perm 1000 --kmer 31 --canonical --seed 42`

### FracMinHash (kmer-sketch)

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `04_fracminhash_sketch.sh` | GTDB genome FASTA.gz files | `data/GTDB/fracminhash_sketches/` | Before `07_fracminhash_pairwise.sh` |
| `07_fracminhash_pairwise.py/.sh` | fracminhash_sketches + candidates CSV | `data/GTDB/fracminhash_pairwise/` | After `04_fracminhash_sketch.sh` |

Parameters: `--algo fracminhash --scale 0.001 --kmer 31 --canonical --seed 42`
(scale=0.001 is equivalent to sourmash's scaled=1000)

### AlphaMaxGeomHash — already run

| Script | Input | Output | When to run |
|--------|-------|--------|-------------|
| `05_alphamaxgeom_sketch.sh` | GTDB genome FASTA.gz files | `data/GTDB/alphamaxgeom_sketches/` | Already run (2026-03-21) |
| `08_alphamaxgeom_pairwise.py/.sh` | alphamaxgeom_sketches + candidates CSV | `data/GTDB/alphamaxgeom_pairwise/` | Already run (2026-03-21) |

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
├── kmc_dbs/                          ← KMC k-mer databases (3.9 TB)
├── kmc_pairwise/                     ← KMC exact pairwise (41 MB NPZ)
├── minhash_sketches/                 ← MinHash sketches + sketch_run_stats.json
├── minhash_pairwise/                 ← MinHash pairwise NPZ + run stats
│   └── sanity_check/                 ← sanity_check_summary.json + scatter PNGs
├── fracminhash_sketches/             ← FracMinHash (kmer-sketch) sketches
├── fracminhash_pairwise/             ← FracMinHash pairwise NPZ + run stats
├── alphamaxgeom_sketches/            ← AlphaMaxGeomHash sketches (14.8 GB)
├── alphamaxgeom_pairwise/            ← AlphaMaxGeomHash pairwise (37 MB NPZ)
│   └── sanity_check/                 ← per-method sanity check output
├── gtdb_pairwise_containment.csv     ← Sourmash FracMinHash pairwise (candidate pairs)
├── sanity_check/                     ← combined sanity check (all methods)
└── figures/                          ← all publication figures (PDF + PNG)
    ├── heatmap_jaccard_relative_error.{pdf,png}
    ├── heatmap_max_containment_relative_error.{pdf,png}
    ├── resources_time.{pdf,png}
    ├── resources_disk.{pdf,png}
    ├── resources_tradeoff.{pdf,png}
    └── resource_summary.json
```

---

## Full run results (AlphaMaxGeomHash, 2026-03-21)

| Metric | Value |
|--------|-------|
| AMG sketching time (192 cores) | 15m 33s |
| AMG pairwise time (192 cores) | 2m 02s |
| AMG sketches disk | 14.8 GB |
| AMG vs KMC Pearson r (Jaccard) | 0.9980 |
| AMG vs KMC MAE (Jaccard) | 0.00234 |
| Pairs compared (AMG vs KMC) | 1,911,078 |

---

## Notes on script numbering and archive

The pipeline was extended after the initial AMG-only run to also support
MinHash and FracMinHash (kmer-sketch).  The original scripts were renumbered
to make room, and the superseded originals were moved to `scripts/archive/`
(still tracked by git for full reproducibility):

| Archived script | Superseded by | Key difference |
|-----------------|---------------|----------------|
| `archive/03_alphamaxgeom_sketch.sh` | `05_alphamaxgeom_sketch.sh` | Same logic; new number |
| `archive/04_alphamaxgeom_pairwise.*` | `08_alphamaxgeom_pairwise.*` | Same logic; new number |
| `archive/05_sanity_check.py` | `09_sanity_check.py` | AMG+Sourmash only → multi-method |
| `archive/06_plot_heatmaps.py` | `10_plot_heatmaps.py` | AMG+Sourmash hardcoded → any kmer-sketch method; Sourmash optional |
| `archive/06_plot_resources.py` | `10_plot_resources.py` | KMC+Sourmash+AMG hardcoded → all four methods; Sourmash optional |
| `archive/06_run_plots.sh` | `10_run_plots.sh` | Calls archived scripts → calls current scripts |

Use the **active scripts** (03–10 in `scripts/`) for all current and future work.

---

## Notes on the Sourmash FracMinHash baseline

The Sourmash FracMinHash pipeline (scripts `make_sourmash_sketches.sh`,
`compute_sourmash_pairwise.sh`) was run as an independent baseline before
the kmer-sketch implementation.  Its pairwise CSV
(`gtdb_pairwise_containment.csv`) is still used as the **candidate pair
filter** for all kmer-sketch pairwise scripts — only pairs with
`max_containment > 0.01` in the Sourmash run are evaluated.

The Sourmash method is **excluded from default figures** because we are now
comparing the three kmer-sketch implementations directly.  To include it,
pass `--include-sourmash` to `09_sanity_check.py`, `10_plot_heatmaps.py`,
`10_plot_resources.py`, or set `INCLUDE_SOURMASH=1` when running
`10_run_plots.sh`.
