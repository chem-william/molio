// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::error::CError;
use crate::formats::amber::AMBERTrajFormat;
use crate::formats::pdb::PDBFormat;
use crate::formats::sdf::SDFFormat;
use crate::formats::smi::SMIFormat;
use crate::formats::xyz::XYZFormat;
use crate::frame::Frame;
use std::fs::File;
use std::io::{BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

/// Supported trajectory file formats.
///
/// `FormatKind` is only responsible for format selection and guessing. The
/// actual format machinery lives in [`FormatReader`] and [`FormatWriter`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FormatKind {
    /// XYZ file format.
    XYZ,
    /// PDB file format.
    PDB,
    /// SMILES file format.
    SMI,
    /// SDF file format.
    SDF,
    /// AMBER NetCDF binary format.
    AMBER,
    /// Automatically detect format from file extension.
    Guess,
}

impl FormatKind {
    /// Guess a format from a file extension.
    ///
    /// # Errors
    ///
    /// Returns an error if the extension is unknown.
    pub fn from_path(path: &Path) -> Result<Self, CError> {
        let ext = path
            .extension()
            .and_then(|value| value.to_str())
            .ok_or_else(|| CError::UnknownFormat(path.display().to_string()))?;

        match ext.to_ascii_lowercase().as_str() {
            "xyz" => Ok(Self::XYZ),
            "pdb" => Ok(Self::PDB),
            "smi" => Ok(Self::SMI),
            "sdf" => Ok(Self::SDF),
            "ncrst" | "nc" => Ok(Self::AMBER),
            _ => Err(CError::UnknownFormat(ext.to_string())),
        }
    }

    pub(crate) fn resolve(self, path: &Path) -> Result<Self, CError> {
        match self {
            Self::Guess => Self::from_path(path),
            kind => Ok(kind),
        }
    }
}

/// Low-level codec trait for text-based formats.
///
/// Text formats implement this trait with the minimal hooks needed to
/// parse/write individual frames. [`TextReader`] and [`TextWriter`] handle
/// frame indexing, seeking, and buffered I/O on top of this.
pub trait Codec {
    /// Reads the next [`Frame`] from `reader`.
    ///
    /// # Errors
    /// Returns an error if reading or parsing the frame fails.
    fn read_next(&mut self, reader: &mut BufReader<File>) -> Result<Frame, CError>;

    /// Advances past the current frame, returning the byte offset of the next
    /// frame start. Returns `None` at end-of-file.
    ///
    /// # Errors
    ///
    /// Returns an error if an I/O or parsing error occurs.
    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError>;

    /// Writes the next [`Frame`] to `writer`.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails.
    fn write_next(&mut self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError>;

    /// Finalizes the file output, performing any necessary cleanup.
    ///
    /// # Errors
    ///
    /// Returns an error if finalization fails.
    fn finalize(&self, writer: &mut BufWriter<File>) -> Result<(), CError>;
}

/// Text-format reader with frame indexing via byte offset scanning.
/// Generic over the codec type — each text format implements [`Codec`].
pub(crate) struct TextReader<C: Codec> {
    codec: C,
    reader: BufReader<File>,
    frame_positions: Vec<u64>,
}

impl<C: Codec> TextReader<C> {
    pub(crate) fn open(path: &Path, codec: C) -> Result<Self, CError> {
        let mut reader = BufReader::with_capacity(64 * 1024, File::open(path)?);

        let mut frame_positions = vec![0];
        while let Some(pos) = codec.forward(&mut reader)? {
            frame_positions.push(pos);
        }
        reader.rewind()?;

        Ok(Self {
            codec,
            reader,
            frame_positions,
        })
    }

    pub(crate) fn read(&mut self) -> Result<Frame, CError> {
        self.codec.read_next(&mut self.reader)
    }

    pub(crate) fn read_at(&mut self, index: usize) -> Result<Frame, CError> {
        let position = *self
            .frame_positions
            .get(index)
            .ok_or_else(|| CError::GenericError("frame index out of bounds".to_string()))?;
        self.reader.seek(SeekFrom::Start(position))?;
        self.codec.read_next(&mut self.reader)
    }

    pub(crate) fn len(&self) -> Result<usize, CError> {
        Ok(self.frame_positions.len().saturating_sub(1))
    }
}

/// Text-format writer. Generic over the codec type.
pub(crate) struct TextWriter<C: Codec> {
    codec: C,
    writer: BufWriter<File>,
}

impl<C: Codec> TextWriter<C> {
    pub(crate) fn create(path: &Path, codec: C) -> Result<Self, CError> {
        Ok(Self {
            codec,
            writer: BufWriter::new(File::create(path)?),
        })
    }

    pub(crate) fn write(&mut self, frame: &Frame) -> Result<(), CError> {
        self.codec.write_next(&mut self.writer, frame)?;
        self.writer.flush()?;
        Ok(())
    }

    pub(crate) fn finish(&mut self) -> Result<(), CError> {
        self.codec.finalize(&mut self.writer)?;
        self.writer.flush()?;
        Ok(())
    }
}

// Dispatches a method call across all variants of a format enum.
macro_rules! format_dispatch {
    ($self:expr, $method:ident $(, $arg:expr)*) => {
        match $self {
            Self::Xyz(inner) => inner.$method($($arg),*),
            Self::Pdb(inner) => inner.$method($($arg),*),
            Self::Smi(inner) => inner.$method($($arg),*),
            Self::Sdf(inner) => inner.$method($($arg),*),
            Self::Amber(inner) => inner.$method($($arg),*),
        }
    };
}

/// Format reader — one variant per supported format.
/// Text formats are wrapped in [`TextReader`], binary formats appear directly.
#[allow(clippy::large_enum_variant)]
pub(crate) enum FormatReader {
    Xyz(TextReader<XYZFormat>),
    Pdb(TextReader<PDBFormat>),
    Smi(TextReader<SMIFormat>),
    Sdf(TextReader<SDFFormat>),
    Amber(AMBERTrajFormat),
}

impl FormatReader {
    pub(crate) fn open(path: &Path, kind: FormatKind) -> Result<Self, CError> {
        match kind {
            FormatKind::XYZ => Ok(Self::Xyz(TextReader::open(path, XYZFormat)?)),
            FormatKind::PDB => Ok(Self::Pdb(TextReader::open(path, PDBFormat::new())?)),
            FormatKind::SMI => Ok(Self::Smi(TextReader::open(path, SMIFormat::default())?)),
            FormatKind::SDF => Ok(Self::Sdf(TextReader::open(path, SDFFormat)?)),
            FormatKind::AMBER => Ok(Self::Amber(AMBERTrajFormat::open(path)?)),
            FormatKind::Guess => unreachable!("Guess should be resolved before reaching FormatReader"),
        }
    }

    pub(crate) fn read(&mut self) -> Result<Frame, CError> {
        format_dispatch!(self, read)
    }

    pub(crate) fn read_at(&mut self, index: usize) -> Result<Frame, CError> {
        format_dispatch!(self, read_at, index)
    }

    pub(crate) fn len(&self) -> Result<usize, CError> {
        format_dispatch!(self, len)
    }
}

/// Format writer — one variant per supported format.
/// Text formats are wrapped in [`TextWriter`], binary formats appear directly.
pub(crate) enum FormatWriter {
    Xyz(TextWriter<XYZFormat>),
    Pdb(TextWriter<PDBFormat>),
    Smi(TextWriter<SMIFormat>),
    Sdf(TextWriter<SDFFormat>),
    Amber(AMBERTrajFormat),
}

impl FormatWriter {
    pub(crate) fn create(path: &Path, kind: FormatKind) -> Result<Self, CError> {
        match kind {
            FormatKind::XYZ => Ok(Self::Xyz(TextWriter::create(path, XYZFormat)?)),
            FormatKind::PDB => Ok(Self::Pdb(TextWriter::create(path, PDBFormat::new())?)),
            FormatKind::SMI => Ok(Self::Smi(TextWriter::create(path, SMIFormat::default())?)),
            FormatKind::SDF => Ok(Self::Sdf(TextWriter::create(path, SDFFormat)?)),
            FormatKind::AMBER => Ok(Self::Amber(AMBERTrajFormat::create(path)?)),
            FormatKind::Guess => unreachable!("Guess should be resolved before reaching FormatWriter"),
        }
    }

    pub(crate) fn write(&mut self, frame: &Frame) -> Result<(), CError> {
        format_dispatch!(self, write, frame)
    }

    pub(crate) fn finish(&mut self) -> Result<(), CError> {
        format_dispatch!(self, finish)
    }
}
