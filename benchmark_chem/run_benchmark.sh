#!/bin/bash

# Set file type (pdb or xyz)
file_type=${1:-pdb}

# Define data files and make sure they exist
if [ "$file_type" = "xyz" ]; then
    data_file="helium.xyz"
elif [ "$file_type" = "pdb" ]; then
    data_file="water.pdb"
else
    echo "Unknown file type: $file_type"
    exit 1
fi

# Check if data file exists
if [ ! -f "$data_file" ]; then
    echo "Error: Data file $data_file not found in the current directory"
    exit 1
fi

# Compile C++ program
(cd cpp && ./compile.sh)

# Build Rust program
(cd .. && cargo build --release)

# Run benchmarks
hyperfine \
    "../target/release/molio $data_file" \
    "./cpp/read_chem $data_file" \
    --export-markdown ./output.md \
    --warmup 10
