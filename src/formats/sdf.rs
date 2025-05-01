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
        // TODO: replace with has_data_left when stabilized
        if reader.fill_buf().map(|b| !b.is_empty()).unwrap() {
            Ok(Some(self.read_next(reader).unwrap()))
        } else {
            Ok(None)
        }
    }

    fn write_next(&self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        todo!();
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        let position = reader.stream_position()?;
        let mut buf = Vec::with_capacity(256); // Pre-allocate a reasonable size

        // Ignore header lines: molecule name, metadata, and general comment line
        for _ in 0..3 {
            buf.clear();
            let bytes_read = reader.read_until(b'\n', &mut buf)?;
            if bytes_read == 0 {
                // EOF reached
                return Ok(None);
            }
        }

        // Read counts line
        buf.clear();
        let bytes_read = reader.read_until(b'\n', &mut buf)?;
        if bytes_read == 0 {
            return Ok(None);
        }
        if bytes_read < 10 {
            return Err(CError::GenericError(format!(
                "counts line must have at least 10 characters in SDF file. It has {bytes_read} bytes: '{}'", String::from_utf8_lossy(&buf)
            )));
        }

        // Parse atom and bond counts
        let counts_line = std::str::from_utf8(&buf)
            .map_err(|_| CError::GenericError("invalid UTF-8 in counts line".to_string()))?;

        let natoms = counts_line[..3].trim().parse::<usize>().map_err(|_| {
            CError::GenericError(format!(
                "could not parse atom count in SDF file: '{}'",
                &counts_line[..3]
            ))
        })?;
        let nbonds = counts_line[3..6].trim().parse::<usize>().map_err(|_| {
            CError::GenericError(format!(
                "could not parse bond count in SDF file: '{}'",
                &counts_line[3..6]
            ))
        })?;

        for _ in 0..(natoms + nbonds) {
            buf.clear();
            if reader.read_until(b'\n', &mut buf)? == 0 {
                return Err(CError::GenericError(
                    "unexpected EOF in SDF format".to_string(),
                ));
            }
        }

        // Search for ending character, updating the cursor in the file for the next
        // call to forward
        loop {
            buf.clear();
            let bytes_read = reader.read_until(b'\n', &mut buf)?;
            if bytes_read == 0 {
                break; // EOF reached
            }

            // Check if the line starts with "$$$$"
            if buf.starts_with(b"$$$$") {
                break;
            }
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
