// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.
//

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Seek},
};

use crate::{error::CError, format::FileFormat, frame::Frame};

pub struct SDFFormat;

impl FileFormat for SDFFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        todo!();
    }

    fn read(&self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        todo!();
    }

    fn write_next(&self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        todo!();
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        let position = reader.stream_position()?;
        let mut line = String::new();
        let mut read_bytes = 0;

        // Ignore header lines: molecule name, metadata, and general comment line
        for _ in 0..3 {
            read_bytes = reader.read_line(&mut line)?;
            line.clear();
        }

        if read_bytes == 0 && line.trim().is_empty() {
            return Ok(None);
        }

        read_bytes = reader.read_line(&mut line)?;
        if read_bytes < 10 {
            return Err(CError::GenericError(format!("counts line must have at least 10 characters in SDF file. It has {read_bytes}: '{line}'")));
        }

        let natoms = line[..3].trim().parse::<usize>().map_err(|_| {
            CError::GenericError(format!("could not parse counts line in SDF file: '{line}'"))
        })?;
        let nbonds = line[3..6].trim().parse::<usize>().map_err(|_| {
            CError::GenericError(format!("could not parse counts line in SDF file: '{line}'"))
        })?;

        for _ in 0..(natoms + nbonds) {
            read_bytes = reader.read_line(&mut line)?;
            if read_bytes == 0 && line.is_empty() {
                return Err(CError::GenericError(
                    "not enough lines in for SDF format".to_string(),
                ));
            }
            line.clear()
        }

        // Search for ending character, updating the cursor in the file for the next
        // call to forward
        while reader.read_line(&mut line)? > 0 {
            if line.trim() == "$$$$" {
                break;
            }
            line.clear();
        }

        // We have enough data to parse an entire molecule.
        // So, even if the file may not have an ending string,
        // return the start of this step
        Ok(Some(position))
    }

    fn finalize(&self, writer: &mut BufWriter<File>) -> Result<(), CError> {
        todo!();
    }
}

mod tests {
    use std::path::Path;

    use crate::trajectory::Trajectory;

    #[test]
    fn check_nsteps_aspirin() {
        let path = Path::new("./src/tests-data/sdf/aspirin.sdf");
        let trajectory = Trajectory::open(path).unwrap();

        assert_eq!(trajectory.size, 1);
    }

    #[test]
    fn check_nsteps_kinases() {
        let path = Path::new("./src/tests-data/sdf/kinases.sdf");
        let trajectory = Trajectory::open(path).unwrap();

        assert_eq!(trajectory.size, 6);
    }
}
