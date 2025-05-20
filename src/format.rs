// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::error::CError;
use crate::formats::pdb::PDBFormat;
use crate::formats::sdf::SDFFormat;
use crate::formats::smi::SMIFormat;
use crate::formats::xyz::XYZFormat;
use crate::frame::Frame;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;

/// Supported text-based trajectory file formats for reading and writing.
///
/// - `XYZ`: plain-text XYZ coordinate format.
/// - `PDB`: Protein Data Bank format.
/// - `SMI`: SMILES string format.
/// - `SDF`: Structure Data File format.
/// - `Guess`: autodetect format from file extension.
#[derive(Clone, Copy)]
pub enum TextFormat {
    /// XYZ file format.
    XYZ,
    /// PDB file format.
    PDB,
    /// SMILES file format.
    SMI,
    /// SDF file format.
    SDF,
    /// Automatically detect format from file extension.
    Guess,
}

/// Concrete file format strategy for reading and writing trajectory data.
pub enum Format<'a> {
    /// Handler for the XYZ format.
    XYZ(XYZFormat),
    /// Handler for the PDB format.
    PDB(PDBFormat<'a>),
    /// Handler for the SMILES format.
    SMI(SMIFormat),
    /// Handler for the SDF format.
    SDF(SDFFormat),
}

impl Format<'_> {
    /// Creates a new [`Format`] by inferring the format from the provided file `path`.
    ///
    /// # Errors
    ///
    /// Returns an error if the file extension is unrecognized.
    pub fn new(path: &Path) -> Result<Self, CError> {
        let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");

        match ext.to_lowercase().as_str() {
            "xyz" => Ok(Format::XYZ(XYZFormat)),
            "pdb" => Ok(Format::PDB(PDBFormat::new())),
            "smi" => Ok(Format::SMI(SMIFormat::default())),
            "sdf" => Ok(Format::SDF(SDFFormat)),
            _ => Err(CError::GenericError("unknown file format".to_string())),
        }
    }
    /// Creates a new `Format` using the specified `TextFormat` and file `path`.
    ///
    /// `TextFormat::Guess` delegates to [`Format::new`].
    ///
    /// # Errors
    ///
    /// Returns an error if format initialization fails or guessing cannot detect the format.
    pub fn new_from_format(fmt: &TextFormat, path: &Path) -> Result<Self, CError> {
        match fmt {
            TextFormat::XYZ => Ok(Format::XYZ(XYZFormat)),
            TextFormat::PDB => Ok(Format::PDB(PDBFormat::new())),
            TextFormat::SMI => Ok(Format::SMI(SMIFormat::default())),
            TextFormat::SDF => Ok(Format::SDF(SDFFormat)),
            TextFormat::Guess => Self::new(path),
        }
    }
}

/// Common interface for reading and writing trajectory file formats.
pub trait FileFormat {
    /// Reads the next [`Frame`] from `reader`.
    ///
    /// # Errors
    /// Returns an error if reading or parsing the frame fails.
    fn read_next(&mut self, reader: &mut BufReader<File>) -> Result<Frame, CError>;

    /// Reads a single [`Frame`], returning `None` at end-of-file.
    ///
    /// # Errors
    ///
    /// Returns an error if an I/O or parsing error occurs.
    fn read(&mut self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError>;
    // fn read_at(&mut self, index: usize) -> Result<Frame, CError>;
    // fn write(&self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError>;

    /// Writes the next [`Frame`] to `writer`.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails.
    fn write_next(&mut self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError>;

    /// Advances to the next frame in `reader`, returning its byte offset.
    ///
    /// # Errors
    ///
    /// Returns an error if an I/O or parsing error occurs.
    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError>;

    /// Finalizes the file output, performing any necessary cleanup (e.g., writing end records).
    /// This should be called when done writing to a file to ensure it's properly closed
    ///
    /// # Errors
    ///
    /// Returns an error if finalization fails.
    fn finalize(&self, writer: &mut BufWriter<File>) -> Result<(), CError>;
}

impl FileFormat for Format<'_> {
    fn read_next(&mut self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        match self {
            Format::XYZ(format) => format.read_next(reader),
            Format::PDB(format) => format.read_next(reader),
            Format::SMI(format) => format.read_next(reader),
            Format::SDF(format) => format.read_next(reader),
        }
    }

    fn read(&mut self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        // TODO: replace with has_data_left when stabilized
        if reader.fill_buf().map(|b| !b.is_empty()).unwrap() {
            match self {
                Format::XYZ(format) => format.read(reader),
                Format::PDB(format) => format.read(reader),
                Format::SMI(format) => format.read(reader),
                Format::SDF(format) => format.read(reader),
            }
        } else {
            Ok(None)
        }
    }

    fn write_next(&mut self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        match self {
            Format::XYZ(format) => format.write_next(writer, frame),
            Format::PDB(format) => format.write_next(writer, frame),
            Format::SMI(format) => format.write_next(writer, frame),
            Format::SDF(format) => format.write_next(writer, frame),
        }
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        match self {
            Format::XYZ(format) => format.forward(reader),
            Format::PDB(format) => format.forward(reader),
            Format::SMI(format) => format.forward(reader),
            Format::SDF(format) => format.forward(reader),
        }
    }

    fn finalize(&self, writer: &mut BufWriter<File>) -> Result<(), CError> {
        match self {
            Format::XYZ(format) => format.finalize(writer),
            Format::PDB(format) => format.finalize(writer),
            Format::SMI(format) => format.finalize(writer),
            Format::SDF(format) => format.finalize(writer),
        }
    }
}
