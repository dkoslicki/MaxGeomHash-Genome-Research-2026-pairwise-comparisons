# Refactor: sourmash scaled=100, threshold=0.001

**Date:** 2026-03-24
**Status:** Code updated; re-run in progress or pending.
**NOT committed to git** — this file is for recovery reference only.

---

## What changed and why

### Previous state (before this refactor)

| Parameter | Old value | File |
|-----------|-----------|------|
| Sourmash sketch scale | `scaled=1000` | `make_sourmash_sketches.sh` |
| Sourmash sketch file | `gtdb_genomes_reps_r226_sketches.siz.zip` (typo: `.siz.` not `.sig.`) | `make_sourmash_sketches.sh`, `compute_sourmash_pairwise.sh` |
| Pairwise candidate threshold | `0.01` | `compute_sourmash_pairwise.sh` |
| Candidate pairs CSV | `gtdb_pairwise_containment.csv` | all pairwise scripts |
| KMC/BottomK/FMH/AMG output dirs | `kmc_pairwise/`, `bottomk_pairwise/`, `fracminhash_pairwise/`, `alphamaxgeom_pairwise/` | bash wrappers + `10_run_plots.sh` |
| Figures dir | `data/GTDB/figures/` | `10_run_plots.sh` |

The old sketches at `scaled=1000` do not have sufficient resolution to reliably
detect `max_containment` as low as 0.001.  Each k-mer is retained with
probability 1/1000, so very low-similarity pairs produce zero intersection even
when true containment is ~0.001.  Moving to `scaled=100` retains 10× more
hashes, making the 0.001 threshold detectable.

The `.siz.zip` filename was a typo; the correct sourmash extension is `.sig.zip`.

### New state (after this refactor)

| Parameter | New value | File |
|-----------|-----------|------|
| Sourmash sketch scale | `scaled=100` | `make_sourmash_sketches.sh` |
| Sourmash sketch file | `gtdb_genomes_reps_r226_sketches_scale_100.sig.zip` | `make_sourmash_sketches.sh`, `compute_sourmash_pairwise.sh` |
| Pairwise candidate threshold | `0.001` | `compute_sourmash_pairwise.sh` |
| Candidate pairs CSV | `gtdb_pairwise_containment_thr0001.csv` | all pairwise scripts |
| KMC output dir | `kmc_pairwise_thr0001/` | `02_kmc_pairwise.sh` |
| BottomK output dir | `bottomk_pairwise_thr0001/` | `06_bottomk_pairwise.sh` |
| FracMinHash output dir | `fracminhash_pairwise_thr0001/` | `07_fracminhash_pairwise.sh` |
| AlphaMaxGeomHash output dir | `alphamaxgeom_pairwise_thr0001/` | `08_alphamaxgeom_pairwise.sh` |
| Figures dir | `data/GTDB/figures_thr0001/` | `10_run_plots.sh` |

---

## Old results preserved (do not delete)

All old results remain untouched at their original locations:

```
data/GTDB/gtdb_genomes_reps_r226_sketches.siz.zip   ← old scaled=1000 sketches
data/GTDB/gtdb_pairwise_containment.csv              ← old threshold=0.01 candidates
data/GTDB/kmc_pairwise/                              ← old KMC results
data/GTDB/bottomk_pairwise/                          ← old BottomK results
data/GTDB/fracminhash_pairwise/                      ← old FracMinHash results
data/GTDB/alphamaxgeom_pairwise/                     ← old AlphaMaxGeomHash results
data/GTDB/figures/                                   ← old figures
```

---

## Files modified in this refactor

### Scripts directory

| File | What changed |
|------|-------------|
| `scripts/make_sourmash_sketches.sh` | Output filename → `*_scale_100.sig.zip`; `scaled=1000` → `scaled=100` |
| `scripts/compute_sourmash_pairwise.sh` | Sketch input → `*_scale_100.sig.zip`; CSV output → `*_thr0001.csv`; `THRESHOLD=0.01` → `0.001`; updated comment |
| `scripts/02_kmc_pairwise.sh` | `--candidates` → `*_thr0001.csv`; `--output` → `kmc_pairwise_thr0001`; `--threshold 0.01` → `0.001` |
| `scripts/06_bottomk_pairwise.sh` | `CANDIDATES` → `*_thr0001.csv`; `OUTPUT_DIR` → `bottomk_pairwise_thr0001`; comment updated |
| `scripts/07_fracminhash_pairwise.sh` | Same pattern as above |
| `scripts/08_alphamaxgeom_pairwise.sh` | Same pattern as above |
| `scripts/06_bottomk_pairwise.py` | Docstring: "0.01" → "0.001" |
| `scripts/07_fracminhash_pairwise.py` | Docstring: "0.01" → "0.001" |
| `scripts/08_alphamaxgeom_pairwise.py` | Docstring: "0.01" → "0.001" |
| `scripts/10_run_plots.sh` | All pairwise dir vars → `*_thr0001`; `SOURMASH_CSV` → `*_thr0001.csv`; `FIGURE_DIR` → `figures_thr0001` |
| `scripts/run_sketch_and_pairwise_plus_plot.sh` | Full rewrite: fixed stale script names (`03_alphamaxgeom_sketch.sh` → `05_*`, etc.); added sourmash steps 1–2 as explicit steps; added `09_sanity_check.py` call; all paths use `_thr0001` directories |

---

## How to re-run everything

```bash
nohup bash /scratch/shared_data/MaxGeomHash_Genome_Research_2026/scripts/run_sketch_and_pairwise_plus_plot.sh \
    > /scratch/shared_data/MaxGeomHash_Genome_Research_2026/scripts/full_pipeline.log 2>&1 &
echo "PID: $!"
# Monitor with:
tail -f /scratch/shared_data/MaxGeomHash_Genome_Research_2026/scripts/full_pipeline.log
```

Or run individual steps manually in order:

```bash
BASE="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"

# Step 1: Sourmash sketches (scaled=100)  ~requires re-run from scratch
bash $BASE/scripts/make_sourmash_sketches.sh

# Step 2: Sourmash pairwise (threshold=0.001)
bash $BASE/scripts/compute_sourmash_pairwise.sh

# Step 3: KMC exact pairwise on new candidates
bash $BASE/scripts/02_kmc_pairwise.sh

# Steps 4-6: kmer-sketch sketching (can run in parallel)
bash $BASE/scripts/03_bottomk_sketch.sh
bash $BASE/scripts/04_fracminhash_sketch.sh
bash $BASE/scripts/05_alphamaxgeom_sketch.sh

# Steps 7-9: kmer-sketch pairwise (can run in parallel)
bash $BASE/scripts/06_bottomk_pairwise.sh
bash $BASE/scripts/07_fracminhash_pairwise.sh
bash $BASE/scripts/08_alphamaxgeom_pairwise.sh

# Step 10: Accuracy comparison
GTDB="$BASE/data/GTDB"
python3 $BASE/scripts/09_sanity_check.py \
    --amg-pairwise         $GTDB/alphamaxgeom_pairwise_thr0001 \
    --bottomk-pairwise     $GTDB/bottomk_pairwise_thr0001 \
    --fracminhash-pairwise $GTDB/fracminhash_pairwise_thr0001 \
    --kmc-pairwise         $GTDB/kmc_pairwise_thr0001 \
    --output               $GTDB/alphamaxgeom_pairwise_thr0001/sanity_check

# Step 11: Figures
bash $BASE/scripts/10_run_plots.sh
```

---

## How to recover the old pipeline state

If you need to go back to `scaled=1000` / threshold=`0.01`, reverse these diffs:

```bash
cd /scratch/shared_data/MaxGeomHash_Genome_Research_2026
git diff scripts/
```

The old data files are still present and untouched — only the scripts changed.
To revert the scripts:

```bash
git checkout -- scripts/make_sourmash_sketches.sh
git checkout -- scripts/compute_sourmash_pairwise.sh
git checkout -- scripts/02_kmc_pairwise.sh
git checkout -- scripts/06_bottomk_pairwise.sh
git checkout -- scripts/07_fracminhash_pairwise.sh
git checkout -- scripts/08_alphamaxgeom_pairwise.sh
git checkout -- scripts/06_bottomk_pairwise.py
git checkout -- scripts/07_fracminhash_pairwise.py
git checkout -- scripts/08_alphamaxgeom_pairwise.py
git checkout -- scripts/10_run_plots.sh
git checkout -- scripts/run_sketch_and_pairwise_plus_plot.sh
```

---

## Expected new output locations

```
data/GTDB/
├── gtdb_genomes_reps_r226_sketches_scale_100.sig.zip   ← NEW sourmash sketches
├── gtdb_pairwise_containment_thr0001.csv               ← NEW candidate pairs
├── kmc_pairwise_thr0001/                               ← NEW KMC results
├── bottomk_pairwise_thr0001/                           ← NEW BottomK results
├── fracminhash_pairwise_thr0001/                       ← NEW FracMinHash results
│   └── sanity_check/                                   ← accuracy summary JSON
├── alphamaxgeom_pairwise_thr0001/                      ← NEW AlphaMaxGeomHash results
│   └── sanity_check/
└── figures_thr0001/                                    ← NEW publication figures
```
