# Clone only the kmc_api subdirectory (no full repo history, no large binaries)
git clone --depth 1 --filter=blob:none --sparse \
    https://github.com/refresh-bio/KMC.git kmc_api_src
cd kmc_api_src
git sparse-checkout set kmc_api
cd ..

# Copy our helper into the source tree so the relative includes resolve
cp kmc_db_info.cpp kmc_api_src/kmc_api/

# Compile from within the directory where all relative paths are correct
g++ -O2 -std=c++17 \
    -o kmc_db_info \
    kmc_api_src/kmc_api/kmc_db_info.cpp \
    kmc_api_src/kmc_api/kmc_file.cpp \
    kmc_api_src/kmc_api/kmer_api.cpp \
    kmc_api_src/kmc_api/mmer.cpp \
    -I kmc_api_src/kmc_api

# Move the binary back to the scripts directory
mv kmc_db_info kmc_db_info  # already here since -o is relative

# Optional: keep the source around or remove it
# rm -rf kmc_api_src
