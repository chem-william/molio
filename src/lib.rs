pub mod angle;
pub mod atom;
pub mod bond;
pub mod connectivity;
pub mod dihedral;
pub mod error;
pub mod extendedxyzparser;
pub mod format;
pub mod formats;
pub mod frame;
pub mod improper;
pub mod property;
pub mod residue;
pub mod topology;
pub mod trajectory;
pub mod unit_cell;

use std::{hint::black_box, path::Path};
use trajectory::Trajectory;

/// Read a trajectory file and return the total number of atoms processed
pub fn read_trajectory(path: &Path) -> usize {
    let mut trajectory = Trajectory::open(path).unwrap();
    let mut total_atoms = 0;
    while let Some(next_frame) = trajectory.read().unwrap() {
        total_atoms += next_frame.size();
    }
    black_box(total_atoms)
}
