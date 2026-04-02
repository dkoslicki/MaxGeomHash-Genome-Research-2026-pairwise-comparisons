# deprecated/

This directory preserves scripts, logs, and documentation from earlier
iterations of the MaxGeomHash pipeline that have been superseded but are
retained for reproducibility.  **Do not use these for new work.**

---

## Contents

### `scripts/make_sourmash_sketches.sh`

**Superseded by:** `scripts/make_fracminhash_sketches.sh`

Produced `data/GTDB/gtdb_genomes_reps_r226_sketches_scale_100.sig.zip` using
`sourmash scripts manysketch` at scaled=100 (scale=0.01).  Abandoned because
`manysketch` corrupted the `.sig.zip` at scaled=100 (FileNotFoundError on the
internal `.sig.gz` entry) and parallelised extremely poorly (1–2 cores used
despite `--cores 500`).  The kmer-sketch binary is used instead.

### `scripts/compute_sourmash_pairwise.sh`

**Superseded by:** `scripts/compute_fracminhash_candidates.sh`

Ran `sourmash scripts pairwise` (branchwater plugin) to produce the all-vs-all
candidate pairs CSV.  Abandoned due to the corrupted `.sig.zip` produced by
`make_sourmash_sketches.sh`.  The kmer-sketch FracMinHash all-vs-all Python
script (`compute_fracminhash_candidates.py`) is used instead.

### `scripts/full_maxgeomhash_pipeline.sh`

**Superseded by:** `scripts/run_sketch_and_pairwise_plus_plot.sh`

Original AMG-only pipeline that called the old-numbered scripts
(`03_alphamaxgeom_sketch.sh`, `04_alphamaxgeom_pairwise.sh`,
`05_sanity_check.py`) — all of which have been archived.  The current master
pipeline adds BottomK and FracMinHash, updates all script numbers, and uses
kmer-sketch throughout.

### `scripts/archive/`

Superseded scripts from pipeline iterations prior to the full three-method
comparison.  All were renumbered or replaced when BottomK was added and when
MinHash was replaced by BottomK.

| File | Superseded by | Key change |
|------|---------------|------------|
| `03_alphamaxgeom_sketch.sh` | `05_alphamaxgeom_sketch.sh` | Renumbered; logic identical |
| `04_alphamaxgeom_pairwise.py` | `08_alphamaxgeom_pairwise.py` | Renumbered; logic identical |
| `04_alphamaxgeom_pairwise.sh` | `08_alphamaxgeom_pairwise.sh` | Renumbered; logic identical |
| `05_sanity_check.py` | `09_sanity_check.py` | AMG+Sourmash only → all four methods |
| `06_plot_heatmaps.py` | `10_plot_heatmaps.py` | AMG-hardcoded → any kmer-sketch method |
| `06_plot_resources.py` | `10_plot_resources.py` | AMG-hardcoded → all four methods |
| `06_run_plots.sh` | `10_run_plots.sh` | Calls archived scripts → current scripts |

### `scripts/logs/`

Run logs from earlier pipeline iterations (scaled=1000 / threshold=0.01 era
and the initial AMG-only runs).  The canonical log for the current full run
(scaled=100, threshold=0.001) is `scripts/full_pipeline.log`.

### `REFACTOR_SCALE100_THR0001.md`

Migration notes written during the 2026-03-24 refactor from
`scaled=1000, threshold=0.01` to `scaled=100, threshold=0.001`.  Retained for
context; the active codebase has been fully updated.

---

## Data preserved alongside deprecated scripts

The following data files in `data/GTDB/` correspond to these deprecated runs
and are retained for reference:

| File / Directory | Era | Notes |
|------------------|-----|-------|
| `gtdb_pairwise_containment.csv` | threshold=0.01 | Old candidate pairs CSV |
| `kmc_pairwise/` | threshold=0.01 | Old KMC exact pairwise results |
| `bottomk_pairwise/` | threshold=0.01 | Old BottomK pairwise |
| `fracminhash_pairwise/` | threshold=0.01 | Old FracMinHash pairwise |
| `alphamaxgeom_pairwise/` | threshold=0.01 | Old AlphaMaxGeomHash pairwise |
| `minhash_sketches/` | MinHash era | MinHash replaced by BottomK |
| `alphamaxgeom_sketches_zipped/` | AMG-only era | Superseded by uncompressed dir |
| `gtdb_genomes_reps_r226_sketches.siz.zip` | sourmash scaled=1000 | Corrupted; do not use |
| `gtdb_genomes_reps_r226_sketches_scale_100.sig.zip` | sourmash scaled=100 | Corrupt at scaled=100 |
| `gtdb_pairwise_ani.csv` | sourmash ANI | Reference only |
