#!/bin/bash

# Set file type (pdb or xyz)
file_type=${1:-pdb}

# Move to project root
cd ..

# Define data files
if [ "$file_type" = "xyz" ]; then
    data_file="./src/tests-data/xyz/helium.xyz"
elif [ "$file_type" = "pdb" ]; then
    data_file="./src/tests-data/pdb/water.pdb"
else
    echo "Unknown file type: $file_type"
    exit 1
fi

# Check if data file exists
if [ ! -f "$data_file" ]; then
    echo "Error: Data file $data_file not found in the current directory"
    exit 1
fi

# Check if cargo-flamegraph is installed
if ! command -v cargo-flamegraph &> /dev/null; then
    echo "cargo-flamegraph is not installed. Install it with:"
    echo "cargo install flamegraph"
    exit 1
fi

# Set environment variables to override the release profile settings
# This enables debug symbols and disables optimizations for better profiling
export CARGO_PROFILE_RELEASE_DEBUG=true
export CARGO_PROFILE_RELEASE_OPT_LEVEL=1
export CARGO_PROFILE_RELEASE_LTO=false
export CARGO_PROFILE_RELEASE_STRIP=false
export CARGO_PROFILE_RELEASE_CODEGEN_UNITS=16

echo "Set environment variables for optimal profiling"

# Build and run flamegraph with RUSTFLAGS for better profiling
echo "Running flamegraph..."
RUSTFLAGS="-C force-frame-pointers=yes" cargo flamegraph --bin molio -- --profile "$data_file"

# Move the flamegraph.svg file to the benchmark directory
if [ -f "flamegraph.svg" ]; then
    mv flamegraph.svg "benchmark_chem/flamegraph_${file_type}.svg"
    echo "Flamegraph saved to benchmark_chem/flamegraph_${file_type}.svg"
else
    echo "Error: flamegraph.svg was not generated"
fi 
