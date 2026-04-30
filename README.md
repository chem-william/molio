[![Codecov](https://codecov.io/github/chem-william/molio/coverage.svg?branch=main)](https://codecov.io/gh/chem-william/molio)
[![dependency status](https://deps.rs/repo/github/chem-william/molio/status.svg)](https://deps.rs/repo/github/chem-william/molio)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19705648.svg)](https://doi.org/10.5281/zenodo.19705648)

# molio

`molio` is a Rust library for reading and writing chemistry files used in computational chemistry simulations. It provides a unified interface to access atomic (positions, velocities, atomic symbols) and trajectory (frames, topology) information and other data across various chemical file formats.

Currently supports

| Format              | Extension |
| ------              | --------- |
| (Extended) XYZ      | .xyz      |
| PDB                 | .pdb      |
| SDF                 | .sdf      |
| SMILES              | .smi      |
| AMBER Trajectories  | .nc       |
| AMBER Traj          | .ncrst    |

This project is heavily inspired by [`chemfiles`](https://github.com/chemfiles/chemfiles/), a modern C++ library with the same purpose.

## Usage

### Reading Trajectories

```rust
use std::path::Path;
use molio::trajectory::Trajectory;

// Open a trajectory file
let path = Path::new("./src/tests-data/xyz/extended.xyz");
let mut trajectory = Trajectory::open(path).unwrap();

// Read a specific frame
let frame = trajectory.read_at(0).unwrap().unwrap();

// Access frame properties
let energy = frame.properties["ENERGY"].expect_double();

// Access atomic properties
let atom = &frame[0];
let unit_cell = frame.unit_cell;
```

### Writing Trajectories

```rust
use std::path::Path;
use molio::{trajectory::Trajectory, frame::Frame, atom::Atom};

// Create a trajectory file for writing
let path = Path::new("output.pdb");
let mut writer = Trajectory::create(path).unwrap();

// Populate Frame with a single C atom with position [1.0, 2.0, 3.0]
let mut frame = Frame::new();
frame.add_atom(Atom::new("C".to_string()), [1.0, 2.0, 3.0]);
// ..add other properties here

// Write frames
writer.write(&frame).unwrap();

// Optionally, explicitly finalize the file (for example, adds END record for PDB)
// If not called explicitly, finalization happens when writer is dropped
writer.finish().unwrap();
```

### Format specific notes

- AMBER: we use [netcdf3](https://docs.rs/netcdf3/latest/netcdf3/), a pure Rust implementation for reading and writing NetCDF-3 files, to read AMBER trajectories. However, due to the way that lib performs, writes to disk won't happen until we hit `finish()`. There are no partial writes.

This also means that we can only read AMBER trajectories based on NetCDF-3 files.

When reading the topology, we assume that the topology file has the name as the `.nc` file. For example, if the trajectory file is titled `trajectory.nc`, the topology must exist in `trajectory.parm7`. Otherwise, topology data won't be read.

## Contributions
[Contributions](https://github.com/chem-william/molio/edit/main/CONTRIBUTING.md) are very welcome! Please open an [issue](https://github.com/chem-william/molio/issues/new) to discuss bugs or new features.

## Citation
Consider citing the original [`chemfiles`](https://github.com/chemfiles/chemfiles/) library.
