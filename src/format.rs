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
/// actual text-format machinery lives in [`TextReader`] and [`TextWriter`].
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
/// parse/write individual frames. The trajectory reader/writer handles frame
/// indexing, seeking, and buffered I/O.
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

// Dispatches a method call across all variants of a codec enum.
macro_rules! dispatch {
    ($self:expr, $method:ident $(, $arg:expr)*) => {
        match $self {
            Self::Xyz(inner) => inner.$method($($arg),*),
            Self::Pdb(inner) => inner.$method($($arg),*),
            Self::Smi(inner) => inner.$method($($arg),*),
            Self::Sdf(inner) => inner.$method($($arg),*),
        }
    };
}

/// Enum wrapping text-based codecs. The I/O state (reader, frame positions)
/// lives in [`TextReader`] / [`TextWriter`], not here.
pub(crate) enum CodecImpl {
    Xyz(XYZFormat),
    Pdb(PDBFormat),
    Smi(SMIFormat),
    Sdf(SDFFormat),
}

impl CodecImpl {
    pub(crate) fn from_kind(kind: FormatKind) -> Result<Self, CError> {
        match kind {
            FormatKind::XYZ => Ok(Self::Xyz(XYZFormat)),
            FormatKind::PDB => Ok(Self::Pdb(PDBFormat::new())),
            FormatKind::SMI => Ok(Self::Smi(SMIFormat::default())),
            FormatKind::SDF => Ok(Self::Sdf(SDFFormat)),
            FormatKind::AMBER => {
                unreachable!("AMBER uses its own binary reader, not CodecImpl")
            }
            FormatKind::Guess => {
                unreachable!("Guess should be handled before reaching CodecImpl")
            }
        }
    }

    pub(crate) fn read_next(&mut self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        dispatch!(self, read_next, reader)
    }

    pub(crate) fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        dispatch!(self, forward, reader)
    }

    pub(crate) fn write_next(
        &mut self,
        writer: &mut BufWriter<File>,
        frame: &Frame,
    ) -> Result<(), CError> {
        dispatch!(self, write_next, writer, frame)
    }

    pub(crate) fn finalize(&self, writer: &mut BufWriter<File>) -> Result<(), CError> {
        dispatch!(self, finalize, writer)
    }
}

/// Text-format reader with frame indexing via byte offset scanning.
/// Owns the codec, the buffered reader, and the frame position index.
pub(crate) struct TextReader {
    codec: CodecImpl,
    reader: BufReader<File>,
    frame_positions: Vec<u64>,
}

impl TextReader {
    pub(crate) fn open(path: &Path, kind: FormatKind) -> Result<Self, CError> {
        let codec = CodecImpl::from_kind(kind)?;
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

    pub(crate) fn len(&self) -> usize {
        self.frame_positions.len().saturating_sub(1)
    }
}

/// Text-format writer.
pub(crate) struct TextWriter {
    codec: CodecImpl,
    writer: BufWriter<File>,
}

impl TextWriter {
    pub(crate) fn create(path: &Path, kind: FormatKind) -> Result<Self, CError> {
        let codec = CodecImpl::from_kind(kind)?;
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

/// Strategy enum dispatching reads to either text or binary format backends.
#[allow(clippy::large_enum_variant)]
pub(crate) enum ReaderStrategy {
    Text(TextReader),
    Binary(AMBERTrajFormat),
}

impl ReaderStrategy {
    pub(crate) fn read(&mut self) -> Result<Frame, CError> {
        match self {
            Self::Text(t) => t.read(),
            Self::Binary(b) => b.read(),
        }
    }

    pub(crate) fn read_at(&mut self, index: usize) -> Result<Frame, CError> {
        match self {
            Self::Text(t) => t.read_at(index),
            Self::Binary(b) => b.read_at(index),
        }
    }

    pub(crate) fn len(&self) -> Result<usize, CError> {
        match self {
            Self::Text(t) => Ok(t.len()),
            Self::Binary(b) => b.len(),
        }
    }
}

/// Strategy enum dispatching writes to either text or binary format backends.
pub(crate) enum WriterStrategy {
    Text(TextWriter),
    Binary(AMBERTrajFormat),
}

impl WriterStrategy {
    pub(crate) fn write(&mut self, frame: &Frame) -> Result<(), CError> {
        match self {
            Self::Text(t) => t.write(frame),
            Self::Binary(b) => b.write(frame),
        }
    }

    pub(crate) fn finish(&mut self) -> Result<(), CError> {
        match self {
            Self::Text(t) => t.finish(),
            Self::Binary(b) => b.finish(),
        }
    }
}
