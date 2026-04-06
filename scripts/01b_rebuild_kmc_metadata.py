#!/usr/bin/env python3
"""
rebuild_kmc_metadata.py
=======================
One-shot recovery tool. Scans KMC_DB_DIR for completed .kmc_pre/.kmc_suf
pairs, looks up the original genome path from the manysketch manifest,
reads the unique k-mer count from each database header via the kmc_db_info
C++ helper (O(1) per database, no k-mer iteration), and writes kmc_metadata.csv.

Run from anywhere:
    python3 rebuild_kmc_metadata.py

All paths are hard-coded to match 01_kmc_count.sh.

PREREQUISITES:
    Compile kmc_db_info once before running this script:
        cd /scratch/shared_data/MaxGeomHash_Genome_Research_2026/scripts
        git clone --depth 1 --filter=blob:none --sparse \\
            https://github.com/refresh-bio/KMC.git kmc_api_src
        cd kmc_api_src && git sparse-checkout set kmc_api && cd ..
        cp kmc_db_info.cpp kmc_api_src/kmc_api/
        g++ -O2 -std=c++17 -o kmc_db_info \\
            kmc_api_src/kmc_api/kmc_db_info.cpp \\
            kmc_api_src/kmc_api/kmc_file.cpp \\
            kmc_api_src/kmc_api/kmer_api.cpp \\
            kmc_api_src/kmc_api/mmer.cpp \\
            -I kmc_api_src/kmc_api
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths — must match 01_kmc_count.sh
# ---------------------------------------------------------------------------
BASE_DIR     = Path("/scratch/shared_data/MaxGeomHash_Genome_Research_2026")
GTDB_DIR     = BASE_DIR / "data" / "GTDB"
KMC_DB_DIR   = GTDB_DIR / "kmc_dbs"
MANIFEST     = GTDB_DIR / "manysketch.csv"   # name,genome_filename,protein_filename
METADATA_CSV = GTDB_DIR / "kmc_metadata.csv"

# kmc_db_info binary — expected alongside this script
SCRIPT_DIR      = Path(__file__).resolve().parent
KMC_DB_INFO_BIN = SCRIPT_DIR / "kmc_db_info"


# ---------------------------------------------------------------------------
# Check for the C++ helper binary
# ---------------------------------------------------------------------------
def check_kmc_db_info() -> bool:
    if not KMC_DB_INFO_BIN.exists():
        print(f"ERROR: kmc_db_info binary not found at {KMC_DB_INFO_BIN}")
        print()
        print("Compile it once with:")
        print(f"  cd {SCRIPT_DIR}")
        print("  git clone --depth 1 --filter=blob:none --sparse \\")
        print("      https://github.com/refresh-bio/KMC.git kmc_api_src")
        print("  cd kmc_api_src && git sparse-checkout set kmc_api && cd ..")
        print("  cp kmc_db_info.cpp kmc_api_src/kmc_api/")
        print("  g++ -O2 -std=c++17 -o kmc_db_info \\")
        print("      kmc_api_src/kmc_api/kmc_db_info.cpp \\")
        print("      kmc_api_src/kmc_api/kmc_file.cpp \\")
        print("      kmc_api_src/kmc_api/kmer_api.cpp \\")
        print("      kmc_api_src/kmc_api/mmer.cpp \\")
        print("      -I kmc_api_src/kmc_api")
        return False
    return True


# ---------------------------------------------------------------------------
# Batch k-mer count read via kmc_db_info
# We chunk the input so we can report progress every CHUNK_SIZE databases.
# ---------------------------------------------------------------------------
CHUNK_SIZE = 1_000

def get_kmer_counts_batch(db_prefixes: list[str]) -> dict[str, int]:
    """
    Pipe prefixes to kmc_db_info in chunks of CHUNK_SIZE.
    Returns a dict: prefix -> total_kmers.
    Reads only the .kmc_pre header per database — works for -cs1 and -cs2.
    """
    total   = len(db_prefixes)
    counts: dict[str, int] = {}

    for chunk_start in range(0, total, CHUNK_SIZE):
        chunk = db_prefixes[chunk_start : chunk_start + CHUNK_SIZE]
        chunk_end = min(chunk_start + CHUNK_SIZE, total)

        result = subprocess.run(
            [str(KMC_DB_INFO_BIN)],
            input="\n".join(chunk).encode(),
            capture_output=True,
        )

        if result.stderr:
            sys.stderr.write(result.stderr.decode())

        for line in result.stdout.decode().splitlines():
            parts = line.strip().split("\t")
            if len(parts) == 2:
                prefix, count_str = parts
                counts[prefix] = int(count_str) if count_str.isdigit() else 0

        print(f"  Read headers: {chunk_end:,} / {total:,} ({100*chunk_end/total:.1f}%)",
              flush=True)

    return counts


# ---------------------------------------------------------------------------
# Load genome_name -> genome_path from manysketch manifest
# ---------------------------------------------------------------------------
print(f"Loading genome paths from {MANIFEST} ...")
name_to_path: dict[str, str] = {}
with open(MANIFEST, newline="") as f:
    for row in csv.DictReader(f):
        name_to_path[row["name"]] = row["genome_filename"]
print(f"  {len(name_to_path):,} entries loaded.")

# ---------------------------------------------------------------------------
# Check for the helper binary before doing any real work
# ---------------------------------------------------------------------------
if not check_kmc_db_info():
    sys.exit(1)

# ---------------------------------------------------------------------------
# Scan KMC_DB_DIR for completed databases
# ---------------------------------------------------------------------------
print(f"\nScanning {KMC_DB_DIR} for completed KMC databases ...")
pre_files = sorted(KMC_DB_DIR.glob("*.kmc_pre"))
print(f"  Found {len(pre_files):,} .kmc_pre files.")

complete_prefixes: list[str] = []
for i, p in enumerate(pre_files):
    # p.stem on "GCF_001234.kmc_pre" -> "GCF_001234"
    bare = KMC_DB_DIR / p.stem
    if (KMC_DB_DIR / (p.stem + ".kmc_suf")).exists():
        complete_prefixes.append(str(bare))
    if (i + 1) % CHUNK_SIZE == 0 or (i + 1) == len(pre_files):
        print(f"  Scanned: {i+1:,} / {len(pre_files):,} .kmc_pre files ...",
              flush=True)

print(f"  {len(complete_prefixes):,} databases have both .kmc_pre and .kmc_suf.")

# ---------------------------------------------------------------------------
# Batch-read all k-mer counts (one subprocess call total)
# ---------------------------------------------------------------------------
print(f"\nReading k-mer counts from database headers via kmc_db_info ...")
print("  (Reads only .kmc_pre header per database — no k-mer iteration.)")

counts_map = get_kmer_counts_batch(complete_prefixes)

n_zero = sum(1 for v in counts_map.values() if v == 0)
if n_zero:
    print(f"  WARNING: {n_zero} databases returned 0 k-mers.")

print(f"  Done. Read counts for {len(counts_map):,} databases.")

# ---------------------------------------------------------------------------
# Build metadata rows
# ---------------------------------------------------------------------------
rows: list[dict] = []
missing_path: list[str] = []
total = len(complete_prefixes)

for i, db_prefix in enumerate(sorted(complete_prefixes)):
    genome_name = Path(db_prefix).name

    if genome_name not in name_to_path:
        missing_path.append(genome_name)
        genome_path = "UNKNOWN"
    else:
        genome_path = name_to_path[genome_name]

    rows.append({
        "genome_id":      genome_name,
        "genome_path":    genome_path,
        "db_prefix":      db_prefix,
        "n_unique_kmers": counts_map.get(db_prefix, 0),
    })

    if (i + 1) % CHUNK_SIZE == 0 or (i + 1) == total:
        print(f"  Built rows: {i+1:,} / {total:,} ({100*(i+1)/total:.1f}%)",
              flush=True)

rows.sort(key=lambda r: r["genome_id"])

# ---------------------------------------------------------------------------
# Write CSV
# ---------------------------------------------------------------------------
print(f"\nWriting {METADATA_CSV} ...")
with open(METADATA_CSV, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["genome_id", "genome_path", "db_prefix", "n_unique_kmers"],
    )
    writer.writeheader()
    writer.writerows(rows)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print(f"\n{'='*60}")
print(f"  Rows written            : {len(rows):,}")
print(f"  Paths not in manifest   : {len(missing_path):,}")
print(f"  Databases with 0 k-mers : {n_zero:,}")
print(f"  Output                  : {METADATA_CSV}")
print(f"{'='*60}")

if missing_path:
    print(f"\nWARNING: {len(missing_path)} names not in manysketch.csv (genome_path='UNKNOWN').")
    print("  First 5:", missing_path[:5])

print()
print("Next step:")
print(f"  python3 02_kmc_pairwise.py \\")
print(f"      --metadata   {METADATA_CSV} \\")
print(f"      --candidates {GTDB_DIR}/gtdb_pairwise_containment.csv \\")
print(f"      --output     {GTDB_DIR}/kmc_pairwise \\")
print(f"      --threshold  0.01 \\")
print(f"      --cores      192")
