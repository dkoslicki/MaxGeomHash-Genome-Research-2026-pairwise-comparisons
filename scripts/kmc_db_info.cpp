// kmc_db_info.cpp
// ================
// Reads total_kmers from KMC database headers in batch.
// Uses ONLY the kmc_api source — no kmc_core library needed.
//
// SETUP (one-time):
//   # Download the 4 kmc_api source files alongside this file:
//   BASE=https://raw.githubusercontent.com/refresh-bio/KMC/master/kmc_api
//   wget ${BASE}/kmc_file.h ${BASE}/kmc_file.cpp \
//        ${BASE}/kmer_api.h ${BASE}/kmer_defs.h
//
//   # Compile:
//   g++ -O2 -std=c++17 -o kmc_db_info kmc_db_info.cpp kmc_file.cpp
//
// USAGE:
//   # One prefix per line on stdin; outputs "prefix<TAB>total_kmers" per line.
//   cat list_of_prefixes.txt | ./kmc_db_info
//   echo "/path/to/kmc_dbs/GCF_001234" | ./kmc_db_info
//
// Works for both -cs1 databases (existing ones) and -cs2 databases (new ones).
// total_kmers is always stored in the .kmc_pre header regardless of counter format.

#include "kmc_file.h"
#include <iostream>
#include <string>
#include <cstdint>

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::string prefix;
    while (std::getline(std::cin, prefix))
    {
        if (prefix.empty()) continue;

        CKMCFile db;

        // OpenForListing reads the .kmc_pre header and LUT (needed for the API),
        // but does NOT load the entire .kmc_suf into RAM — only buffers it.
        // For header-only reads this is the correct lightweight open mode.
        if (!db.OpenForListing(prefix))
        {
            std::cerr << "Warning: cannot open KMC database: " << prefix << "\n";
            std::cout << prefix << "\t0\n";
            continue;
        }

        // Info() reads directly from the in-memory header populated by OpenForListing.
        // total_kmers is the exact count of unique k-mers stored, regardless of
        // whether -cs1 (no counters) or -cs2 (capped counters) was used during counting.
        uint32_t kmer_length, mode, counter_size, lut_prefix_length, signature_len;
        uint32_t min_count;
        uint64_t max_count, total_kmers;
        db.Info(kmer_length, mode, counter_size, lut_prefix_length,
                signature_len, min_count, max_count, total_kmers);

        db.Close();

        std::cout << prefix << "\t" << total_kmers << "\n";
    }

    return 0;
}
