# Chemical File Format Benchmarks

This directory contains benchmarks comparing `molio` (Rust) with [`chemfiles`](https://github.com/chemfiles/chemfiles/) (C++) for reading chemical file formats.

## Requirements

- Rust (for `molio`)
- C++ compiler and `chemfiles` library (for the C++ benchmark)
- [hyperfine](https://github.com/sharkdp/hyperfine) for running hyperfine benchmarks
- (Optional) [cargo-flamegraph](https://github.com/flamegraph-rs/flamegraph) for profiling

## Running Benchmarks

You need to install `chemfiles` first and then point the `compile.sh` script to the location of where you've installed `chemfiles`.

### Hyperfine Benchmarks (comparing with C++)

To run the benchmarks using [`hyperfine`](https://github.com/sharkdp/hyperfine) (comparing with C++):

```bash
# For PDB benchmarks (default)
./run_benchmark.sh pdb

# For XYZ benchmarks
./run_benchmark.sh xyz
```

The benchmark results will be saved to `output.md` in Markdown format.

### Criterion Benchmarks (Rust only)

For more detailed Rust-only benchmarks using Criterion:

```bash
# From the project root
cargo bench
```

Criterion will generate HTML reports in `target/criterion/`.

### Profiling with Flamegraph

To profile the code and generate a [`flamegraph`](https://www.brendangregg.com/flamegraphs.html):

```bash
# For PDB benchmarks (default)
./run_flamegraph.sh pdb

# For XYZ benchmarks
./run_flamegraph.sh xyz
```

The flamegraph will be saved as `flamegraph_pdb.svg` or `flamegraph_xyz.svg` in the benchmark directory.

## How It Works

The benchmarks compare:
1. A Rust implementation using `molio`
2. A C++ implementation using `chemfiles`

Both implementations read the same file (either `water.pdb` or `helium.xyz`) and perform the same operations (iterate through all frames and count atoms).

## Directory Structure

- `cpp/` - C++ benchmark code using chemfiles
  - `read_chem.cpp` - C++ implementation for all format types
  - `compile.sh` - Script to compile the C++ benchmark
- `water.pdb` - Sample PDB file for benchmarking
- `helium.xyz` - Sample XYZ file for benchmarking
- `run_benchmark.sh` - Script to run `hyperfine` benchmarks
- `run_flamegraph.sh` - Script to run profiling with `flamegraph`
- `output.md` - Benchmark results 
