# smoke_test/

End-to-end validation for the MaxGeomHash GTDB pairwise similarity pipeline.

---

## Purpose

Runs the complete pipeline — FracMinHash sketching, candidate pair generation,
KMC exact pairwise, kmer-sketch pairwise (BottomK, FracMinHash,
AlphaMaxGeomHash), sanity-check accuracy reporting, and all publication figure
scripts — on a small, curated 100-genome subset.  Verifies correctness before
or after any change to the pipeline scripts.

## Genome subset design

Rather than using the first N genomes from the manifest (which are highly
diverse and produce few or no candidate pairs), the smoke test selects 100
genomes that appear in actual candidate pairs from the full run
(`data/GTDB/gtdb_pairwise_containment_thr0001.csv`).  These 100 genomes
produce 913 candidate pairs at threshold=0.001, providing meaningful pairwise
and accuracy work in every step.

## Running the smoke test

```bash
cd /scratch/shared_data/MaxGeomHash_Genome_Research_2026
bash smoke_test/run_smoke_test.sh
```

All output is written to `smoke_test/data/` (gitignored).  To adjust the
core count (default 8):

```bash
SMOKE_CORES=16 bash smoke_test/run_smoke_test.sh
```

To override the conda environment used for the Datashader step (default
`sourmash`):

```bash
CONDA_ENV=myenv bash smoke_test/run_smoke_test.sh
```

## Prerequisites

The following one-time setup steps must have been completed before running
the smoke test:

| Prerequisite | Script |
|---|---|
| KMC k-mer databases | `scripts/01_kmc_count.sh` |
| kmer-sketch binary | `cd scripts/kmer-sketch && make` |
| Decompressed FASTA files | `scripts/00_decompress_genomes.sh` |
| Full candidates CSV | `scripts/compute_fracminhash_candidates.sh` |

Step 8 (full Datashader heatmaps) additionally requires the `sourmash` conda
environment (or whichever env is specified by `CONDA_ENV`) to have `datashader`
installed.  If it is absent the step is skipped with a warning rather than
failing.

## Expected output

The smoke test passes when all of the following are true:

- **Sketching**: 100 FracMinHash, BottomK, and AlphaMaxGeomHash sketches created.
- **Candidate generation**: 913 candidate pairs at threshold=0.001 (this is
  validated against the full-run CSV; exact match is expected since the same
  100 genomes produce the same hash-based pairs).
- **KMC pairwise**: 681 pairs above threshold in exact counting (some of the
  913 candidates fall below threshold in the exact computation).
- **kmer-sketch pairwise**: 913 pairs for each of the three methods.
- **Accuracy** (approximate values from 2026-03-27 reference run):

| Comparison | Metric | MAE | Pearson r |
|---|---|---|---|
| AlphaMaxGeomHash vs KMC | Jaccard | ~0.000371 | ~0.9982 |
| AlphaMaxGeomHash vs KMC | max_containment | ~0.001249 | ~0.9972 |
| BottomK vs KMC | Jaccard | ~0.000518 | ~0.9949 |
| BottomK vs KMC | max_containment | ~0.001560 | ~0.9949 |
| FracMinHash vs KMC | Jaccard | ~0.000086 | ~0.9999 |
| FracMinHash vs KMC | max_containment | ~0.000268 | ~0.9998 |

## Approximate timing (8 cores)

| Step | Time |
|---|---|
| Sketching (100 genomes × 3 algorithms) | ~22 s |
| Candidate generation | ~1 s |
| KMC pairwise (913 pairs) | ~30 s |
| kmer-sketch pairwise (3 × 913 pairs) | ~7 s |
| Sanity check | ~2 s |
| Heatmap figures — dense subset (2 metrics) | ~4 s |
| Resource figures (time, disk, accuracy L1) | ~2 s |
| Full Datashader heatmaps — spectral (2 metrics) | ~10 s |
| **Total** | **~1–2 min** |

## Output layout

```
smoke_test/data/
├── manifest.csv                  ← 100-genome sub-manifest
├── genome_list.txt               ← paths to the 100 FASTA files
├── candidates_thr0001.csv        ← candidate pairs extracted from full run
├── candidates_fresh_thr0001.csv  ← candidates regenerated from scratch (Step 2)
├── fracminhash_sketches/         ← 100 FracMinHash sketches
├── bottomk_sketches/             ← 100 BottomK sketches
├── alphamaxgeom_sketches/        ← 100 AlphaMaxGeomHash sketches
├── kmc_pairwise/                 ← KMC exact results (681 pairs)
│   └── pairwise_results.npz
├── bottomk_pairwise/             ← BottomK results (913 pairs)
│   └── pairwise_results.npz
├── fracminhash_pairwise/         ← FracMinHash results (913 pairs)
│   └── pairwise_results.npz
├── alphamaxgeom_pairwise/        ← AlphaMaxGeomHash results (913 pairs)
│   └── pairwise_results.npz
├── resources_root/               ← fake project root used by 10_plot_resources.py
│   └── data/GTDB/                   (symlinks + per-algorithm stats JSONs)
├── sanity_check/
│   ├── sanity_check_summary.json           ← accuracy metrics (JSON)
│   ├── scatter_amg_vs_kmc_jaccard.png
│   ├── scatter_amg_vs_kmc_maxcont.png
│   ├── scatter_bk_vs_kmc_jaccard.png
│   ├── scatter_bk_vs_kmc_maxcont.png
│   ├── scatter_fmh_ks_vs_kmc_jaccard.png
│   └── scatter_fmh_ks_vs_kmc_maxcont.png
└── figures/
    ├── heatmap_jaccard_relative_error.{pdf,png}         ← dense 100-genome heatmaps
    ├── heatmap_max_containment_relative_error.{pdf,png}
    ├── heatmap_jaccard_relative_error_full_spectral.{pdf,png}   ← Datashader heatmaps
    ├── heatmap_max_containment_relative_error_full_spectral.{pdf,png}
    ├── resources_time.{pdf,png}
    ├── resources_disk.{pdf,png}
    ├── resources_tradeoff.{pdf,png}
    ├── accuracy_l1_jaccard.{pdf,png}
    └── accuracy_l1_max_containment.{pdf,png}
```
