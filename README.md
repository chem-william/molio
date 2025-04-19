# molio

A Rust library for parsing molecular trajectory files.

## Usage

```rust
use std::path::Path;
use your_crate::Trajectory;

// Open a trajectory file
let path = Path::new("trajectory.xyz");
let mut trajectory = Trajectory::new(path)?;

// Read a specific frame
let frame = trajectory.read_at(0)?;

// Access frame properties
let energy = frame.properties["Energy"].expect_double();

// Access atomic properties
let atom = &frame[0];
let charge = atom.properties["charge"].expect_double();
```
