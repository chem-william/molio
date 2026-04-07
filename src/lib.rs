// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
//
// See LICENSE at the project root for full text.

#![doc = include_str!("../README.md")]
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

use std::{fs::File, io::BufReader};
use std::{hint::black_box, path::Path};

use format::Format;

/// Read a trajectory file and return the total number of atoms processed
pub fn read_trajectory(path: &Path) -> usize {
    let mut format = Format::new(path).unwrap();
    let file = File::open(path).unwrap();
    let mut reader = BufReader::new(file);
    let mut total_atoms = 0;
    while let Some(next_frame) = format.read(&mut reader).unwrap() {
        total_atoms += next_frame.size();
    }
    black_box(total_atoms)
}
