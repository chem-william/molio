use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
};

use crate::{error::CError, format::FileFormat, frame::Frame};

pub struct XTCFormat;

impl XTCFormat {}

impl FileFormat for XTCFormat {
    fn read_next(&mut self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        todo!();
    }

    fn read(&mut self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        todo!();
    }

    fn write_next(&mut self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        todo!();
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        todo!();
    }

    fn finalize(&self, writer: &mut BufWriter<File>) -> Result<(), CError> {
        todo!();
    }
}
