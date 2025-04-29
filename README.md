[![Codecov](https://codecov.io/github/chem-william/molio/coverage.svg?branch=main)](https://codecov.io/gh/chem-william/molio)

# molio

`molio` is a Rust library for reading and writing chemistry files used in computational chemistry simulations. It provides a unified interface to access atomic (positions, velocities, atomic symbols) and trajectory (frames, topology) information and other data across various chemical file formats. Currently supports

| Format         | Extension |
| ------         | --------- |
| (Extended) XYZ | .xyz      |
| PDB            | .pdb      |

This project is a Rust port of [`chemfiles`](https://github.com/chemfiles/chemfiles/), a modern C++ library with the same purpose.

## Usage

### Reading Trajectories

```rust
use std::path::Path;
use molio::Trajectory;

// Open a trajectory file
let path = Path::new("trajectory.xyz");
let mut trajectory = Trajectory::open(path)?;

// Read a specific frame
let frame = trajectory.read_at(0)?;

// Access frame properties
let energy = frame.properties["Energy"].expect_double();

// Access atomic properties
let atom = &frame[0];
let charge = atom.properties["charge"].expect_double();
```

### Writing Trajectories

```rust
use std::path::Path;
use molio::{Trajectory, Frame};

// Create a trajectory file for writing
let path = Path::new("output.pdb");
let mut writer = Trajectory::create(path)?;

// Write frames
writer.write(&frame1)?;
writer.write(&frame2)?;

// Optionally, explicitly finalize the file (adds END record for PDB)
// If not called explicitly, finalization happens when writer is dropped
writer.finish()?;
```

## Contributions
[Contributions](https://github.com/chem-william/molio/edit/main/CONTRIBUTING.md) are very welcome! Please open an [issue](https://github.com/chem-william/molio/issues/new) to discuss bugs or new features.

## Citation
If you use this library in a scientific publication, please cite it using the DOI and BibTeX entry below
```bibtex
TODO
```

Also, consider citing the original [`chemfiles`](https://github.com/chemfiles/chemfiles/) library.