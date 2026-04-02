  nohup bash -c '
  BASE="/scratch/shared_data/MaxGeomHash_Genome_Research_2026"

  echo "========================================================"
  echo "  AlphaMaxGeomHash Full Pipeline"
  echo "  Start: $(date)"
  echo "========================================================"

  echo ""
  echo "--- STEP 1/3: Sketching all 143,614 genomes ---"
  echo "Start: $(date)"
  bash "${BASE}/scripts/03_alphamaxgeom_sketch.sh"
  echo "End:   $(date)"

  echo ""
  echo "--- STEP 2/3: Pairwise similarity (2,084,021 candidate pairs) ---"
  echo "Start: $(date)"
  bash "${BASE}/scripts/04_alphamaxgeom_pairwise.sh"
  echo "End:   $(date)"

  echo ""
  echo "--- STEP 3/3: Sanity check vs FracMinHash and KMC ---"
  echo "Start: $(date)"
  python3 "${BASE}/scripts/05_sanity_check.py" \
      --amg-pairwise  "${BASE}/data/GTDB/alphamaxgeom_pairwise" \
      --fracminhash   "${BASE}/data/GTDB/gtdb_pairwise_containment.csv" \
      --kmc-pairwise  "${BASE}/data/GTDB/kmc_pairwise" \
      --output        "${BASE}/data/GTDB/alphamaxgeom_pairwise/sanity_check"
  echo "End:   $(date)"

  echo ""
  echo "========================================================"
  echo "  Pipeline COMPLETE"
  echo "  End: $(date)"
  echo "========================================================"
  ' > /scratch/shared_data/MaxGeomHash_Genome_Research_2026/scripts/alphamaxgeom_full_run.log 2>&1 &
  disown $!
  echo "Pipeline PID: $!  —  monitor with:"
  echo "  tail -f /scratch/shared_data/MaxGeomHash_Genome_Research_2026/scripts/alphamaxgeom_full_run.log"
