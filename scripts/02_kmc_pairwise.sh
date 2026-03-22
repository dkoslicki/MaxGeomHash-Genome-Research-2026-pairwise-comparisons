#!/bin/bash
python3 02_kmc_pairwise.py \
    --metadata   /scratch/shared_data/MaxGeomHash_Genome_Research_2026/data/GTDB/kmc_metadata.csv \
    --candidates /scratch/shared_data/MaxGeomHash_Genome_Research_2026/data/GTDB/gtdb_pairwise_containment.csv \
    --output     /scratch/shared_data/MaxGeomHash_Genome_Research_2026/data/GTDB/kmc_pairwise \
    --threshold  0.01 \
    --cores      192 \
    --chunk-size 1
