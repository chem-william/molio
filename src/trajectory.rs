// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::error::CError;
use crate::format::{FormatKind, FormatReader, FormatWriter};
use crate::frame::Frame;
use log::error;
use std::path::Path;

/// A handle to a trajectory file for reading.
pub struct TrajectoryReader {
    /// Number of frames in the file
    size: usize,

    strategy: FormatReader,
    /// Next index that will be read.
    index: usize,
}

/// A handle to a trajectory file for writing.
pub struct TrajectoryWriter {
    strategy: FormatWriter,
    frame_count: usize,
}

/// Index type that guarantees the frame exists in a trajectory.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct FrameIndex(usize);

impl FrameIndex {
    /// Creates a new frame index if it is within the valid range.
    #[must_use]
    pub fn new(index: usize, max: usize) -> Option<Self> {
        if index < max { Some(Self(index)) } else { None }
    }

    /// Get the underlying index value.
    #[must_use]
    pub fn value(self) -> usize {
        self.0
    }
}

/// Factory functions for creating trajectory readers and writers.
pub struct Trajectory;

impl Trajectory {
    /// Opens a trajectory file for reading, using format autodetection.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or the format is unknown.
    pub fn open(path: impl AsRef<Path>) -> Result<TrajectoryReader, CError> {
        Self::open_with_format(path, FormatKind::Guess)
    }

    /// Opens a trajectory file for reading with an explicit format.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or the format fails to
    /// initialize.
    pub fn open_with_format(
        path: impl AsRef<Path>,
        format: FormatKind,
    ) -> Result<TrajectoryReader, CError> {
        let kind = format.resolve(path.as_ref())?;

        let strategy = FormatReader::open(path.as_ref(), kind)?;
        let size = strategy.len()?;

        Ok(TrajectoryReader {
            size,
            strategy,
            index: 0,
        })
    }

    pub fn append(path: impl AsRef<Path>) -> Result<TrajectoryWriter, CError> {
        Self::append_with_format(path, FormatKind::Guess)
    }

    pub fn append_with_format(
        path: impl AsRef<Path>,
        format: FormatKind,
    ) -> Result<TrajectoryWriter, CError> {
        let kind = format.resolve(path.as_ref())?;

        let strategy = FormatWriter::open(path.as_ref(), kind)?;

        Ok(TrajectoryWriter {
            strategy,
            frame_count: 0,
        })
    }

    /// Creates a new trajectory file for writing, using format autodetection.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created or the format is
    /// unknown.
    pub fn create(path: impl AsRef<Path>) -> Result<TrajectoryWriter, CError> {
        Self::create_with_format(path, FormatKind::Guess)
    }

    /// Creates a new trajectory file for writing with an explicit format.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created or the format
    /// fails to initialize.
    pub fn create_with_format(
        path: impl AsRef<Path>,
        format: FormatKind,
    ) -> Result<TrajectoryWriter, CError> {
        let kind = format.resolve(path.as_ref())?;

        let strategy = FormatWriter::create(path.as_ref(), kind)?;

        Ok(TrajectoryWriter {
            strategy,
            frame_count: 0,
        })
    }
}

impl TrajectoryReader {
    /// Reads the next frame from the trajectory.
    ///
    /// # Errors
    ///
    /// Returns an error if the format fails while reading.
    pub fn read(&mut self) -> Result<Option<Frame>, CError> {
        if self.index >= self.size {
            return Ok(None);
        }

        let mut frame = self.strategy.read()?;
        frame.set_frame_index(self.index);
        self.index += 1;
        Ok(Some(frame))
    }

    /// Reads a specific frame from the trajectory.
    ///
    /// # Errors
    ///
    /// Returns an error if the frame cannot be read.
    pub fn read_frame(&mut self, index: FrameIndex) -> Result<Frame, CError> {
        let index = index.value();
        let mut frame = self.strategy.read_at(index)?;
        frame.set_frame_index(index);
        self.index = index + 1;
        Ok(frame)
    }

    /// Reads a specific frame by index.
    ///
    /// Returns `Ok(None)` when the index is out of bounds.
    ///
    /// # Errors
    ///
    /// Returns an error if the frame cannot be read.
    pub fn read_at(&mut self, index: usize) -> Result<Option<Frame>, CError> {
        match self.frame_index(index) {
            Some(index) => self.read_frame(index).map(Some),
            None => Ok(None),
        }
    }

    /// Gets a frame index if it is valid for this trajectory.
    #[must_use]
    pub fn frame_index(&self, index: usize) -> Option<FrameIndex> {
        FrameIndex::new(index, self.size)
    }

    /// Returns an iterator over all frames in the trajectory.
    pub fn frames(&mut self) -> impl Iterator<Item = Result<Frame, CError>> + '_ {
        let indices = 0..self.size;
        indices.filter_map(move |i| self.frame_index(i).map(|index| self.read_frame(index)))
    }

    /// Returns the number of frames in [`Self`].
    pub fn len(&self) -> usize {
        self.size
    }

    /// Returns whether [`Self`] is empty.
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }
}

impl TrajectoryWriter {
    /// Writes a frame to the trajectory.
    ///
    /// # Errors
    ///
    /// Returns an error if writing fails.
    pub fn write(&mut self, frame: &Frame) -> Result<(), CError> {
        self.strategy.write(frame)?;
        self.frame_count += 1;
        Ok(())
    }

    /// Finalizes the trajectory file.
    ///
    /// # Errors
    ///
    /// Returns an error if finalization fails.
    pub fn finish(&mut self) -> Result<(), CError> {
        self.strategy.finish()
    }

    /// Returns the number of frames written so far.
    #[must_use]
    pub fn frame_count(&self) -> usize {
        self.frame_count
    }
}

impl Drop for TrajectoryWriter {
    fn drop(&mut self) {
        if let Err(error) = self.finish() {
            error!("warning: Failed to finalize trajectory file: {error}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    const TEST_XYZ_PATH: &str = "./src/tests-data/xyz/extended.xyz";

    #[test]
    fn test_trajectory_open() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let reader = Trajectory::open(&path).unwrap();
        assert_eq!(reader.size, 3);
    }

    #[test]
    fn test_trajectory_with_format() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let reader = Trajectory::open_with_format(&path, FormatKind::XYZ).unwrap();
        assert_eq!(reader.size, 3);
    }

    #[test]
    fn test_read_frame() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let mut reader = Trajectory::open(&path).unwrap();

        // Read first frame
        let idx0 = reader.frame_index(0).unwrap();
        let frame0 = reader.read_frame(idx0).unwrap();
        assert_eq!(frame0.size(), 192);

        // Read second frame
        let idx1 = reader.frame_index(1).unwrap();
        let frame1 = reader.read_frame(idx1).unwrap();
        assert_eq!(frame1.size(), 62);
    }

    #[test]
    fn test_invalid_frame_index() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let reader = Trajectory::open(&path).unwrap();
        assert!(reader.frame_index(4).is_none());
    }

    #[test]
    fn test_read_at() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let mut reader = Trajectory::open(&path).unwrap();

        // Read first frame
        let frame0 = reader.read_at(0).unwrap().unwrap();
        assert_eq!(frame0.size(), 192);

        // Read second frame
        let frame1 = reader.read_at(1).unwrap().unwrap();
        assert_eq!(frame1.size(), 62);

        // Try to read frame at invalid index
        let invalid_frame = reader.read_at(4).unwrap();
        assert!(invalid_frame.is_none());
    }
}
