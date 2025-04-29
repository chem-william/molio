// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Seek},
};

use log::warn;

use crate::{error::CError, format::FileFormat, frame::Frame};
pub struct SMIFormat;

impl SMIFormat {}

impl FileFormat for SMIFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        todo!()
    }

    fn read(&self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        todo!()
    }

    fn write_next(
        &self,
        writer: &mut BufWriter<File>,
        frame: &crate::frame::Frame,
    ) -> Result<(), crate::error::CError> {
        todo!()
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        let pos: u64;

        let mut line = String::new();
        loop {
            match reader.read_line(&mut line)? {
                0 => return Ok(None),                    // EOF
                _ if line.trim().is_empty() => continue, // skip blank
                _ => {
                    pos = reader.stream_position()?; // position after this line
                    break;
                }
            }
        }

        Ok(Some(pos))
    }
    // optional<uint64_t> SMIFormat::forward() {
    //     auto position = file_.tellpos();

    //     auto line = file_.readline();
    //     while (trim(line).empty()) {
    //         if (file_.eof()) {
    //             return nullopt;
    //         }
    //         line = file_.readline();
    //     }

    //     return position;
    // }

    fn finalize(
        &self,
        writer: &mut std::io::BufWriter<std::fs::File>,
    ) -> Result<(), crate::error::CError> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::trajectory::Trajectory;
    use std::path::Path;

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/smi/test.smi");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 8);
    }

    #[test]
    fn check_nsteps_with_newlines() {
        let path = Path::new("./src/tests-data/smi/spaces.smi");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 8);
    }
}
