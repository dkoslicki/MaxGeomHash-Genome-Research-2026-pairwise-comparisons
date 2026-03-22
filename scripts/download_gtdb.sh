#!/bin/bash
wget -P ../data/GTDB https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/genomic_files_reps/gtdb_genomes_reps_r226.tar.gz
cd ../data/GTDB
tar -xzvf gtdb_genomes_reps_r226.tar.gz
