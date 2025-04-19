use crate::error::CError;
use crate::formats::pdb::PDBFormat;
use crate::formats::xyz::XYZFormat;
use crate::frame::Frame;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub enum TextFormat {
    XYZ,
    PDB,
    Guess,
}

pub enum Format {
    XYZ(XYZFormat),
    PDB(PDBFormat),
}

impl Format {
    pub fn new(path: &Path) -> Result<Self, CError> {
        let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");

        match ext.to_lowercase().as_str() {
            "xyz" => Ok(Format::XYZ(XYZFormat)),
            "pdb" => Ok(Format::PDB(PDBFormat)),
            _ => Err(CError::GenericError("unknown file format".to_string())),
        }
    }
    pub fn new_from_format(fmt: TextFormat, path: &Path) -> Result<Self, CError> {
        match fmt {
            TextFormat::XYZ => Ok(Format::XYZ(XYZFormat)),
            TextFormat::PDB => Ok(Format::PDB(PDBFormat)),
            TextFormat::Guess => Self::new(path),
        }
    }
}

pub trait FileFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError>;
    fn read(&self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError>;
    // fn read_at(&mut self, index: usize) -> Result<Frame, CError>;
    fn write(&self, path: &Path, frame: &Frame) -> Result<(), CError>;
    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError>;
}

impl FileFormat for Format {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        match self {
            Format::XYZ(format) => format.read_next(reader),
            Format::PDB(format) => format.read_next(reader),
        }
    }
    fn read(&self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        match self {
            Format::XYZ(format) => format.read(reader),
            Format::PDB(format) => format.read(reader),
        }
    }

    fn write(&self, path: &Path, frame: &Frame) -> Result<(), CError> {
        match self {
            Format::XYZ(format) => format.write(path, frame),
            Format::PDB(format) => format.write(path, frame),
        }
    }
    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        match self {
            Format::XYZ(format) => format.forward(reader),
            Format::PDB(format) => format.forward(reader),
        }
    }
}
