// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.
//

use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter},
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
        todo!();
    }

    fn finalize(&self, writer: &mut BufWriter<File>) -> Result<(), CError> {
        todo!();
    }
}

mod tests {
    use std::path::Path;

    use crate::trajectory::Trajectory;

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/sdf/aspirin.sdf");
        let trajectory = Trajectory::open(path).unwrap();
    }
}
