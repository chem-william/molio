use crate::error::CError;
use crate::format::FileFormat;
use crate::format::Format;
use crate::format::TextFormat;
use crate::frame::Frame;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::io::{BufReader, Seek};
use std::path::Path;

/// A handle to a trajectory file for reading
pub struct TrajectoryReader<'a> {
    pub size: usize,
    strategy: Format<'a>,
    reader: BufReader<File>,
    frame_positions: Vec<u64>,
}

/// A handle to a trajectory file for writing
pub struct TrajectoryWriter<'a> {
    strategy: Format<'a>,
    writer: BufWriter<File>,
    frame_count: usize,
}

/// Index type that guarantees the frame exists in a trajectory
#[derive(Debug, Clone, Copy)]
pub struct FrameIndex(usize);

impl FrameIndex {
    /// Creates a new frame index if it's within the valid range for the trajectory
    pub fn new(index: usize, max: usize) -> Option<Self> {
        if index < max {
            Some(FrameIndex(index))
        } else {
            None
        }
    }

    /// Get the underlying index value
    pub fn value(&self) -> usize {
        self.0
    }
}

/// Factory functions for creating trajectory readers and writers
pub struct Trajectory;

impl Trajectory {
    /// Opens a trajectory file for reading, using format autodetection
    pub fn open(path: &Path) -> Result<TrajectoryReader, CError> {
        Self::open_with_format(path, TextFormat::Guess)
    }

    /// Opens a trajectory file for reading with an explicit format
    pub fn open_with_format(path: &Path, fmt: TextFormat) -> Result<TrajectoryReader, CError> {
        let strategy = Format::new_from_format(&fmt, path)?;
        let file = File::open(path).map_err(CError::IoError)?;
        let mut reader = BufReader::new(file);

        let frame_positions = TrajectoryReader::scan_all(&mut reader, &strategy);
        let size = frame_positions.len() - 1;

        Ok(TrajectoryReader {
            size,
            strategy,
            reader,
            frame_positions,
        })
    }

    /// Creates a new trajectory file for writing, using format autodetection
    pub fn create(path: &Path) -> Result<TrajectoryWriter, CError> {
        Self::create_with_format(path, TextFormat::Guess)
    }

    /// Creates a new trajectory file for writing with an explicit format
    pub fn create_with_format(path: &Path, fmt: TextFormat) -> Result<TrajectoryWriter, CError> {
        let strategy = Format::new_from_format(&fmt, path)?;
        let file = File::create(path).map_err(CError::IoError)?;
        let writer = BufWriter::new(file);

        Ok(TrajectoryWriter {
            strategy,
            writer,
            frame_count: 0,
        })
    }
}

impl<'a> TrajectoryReader<'a> {
    /// Reads the next frame from the trajectory
    pub fn read(&mut self) -> Result<Option<Frame>, CError> {
        self.strategy.read(&mut self.reader)
    }

    /// Reads a specific frame from the trajectory
    pub fn read_frame(&mut self, index: FrameIndex) -> Result<Frame, CError> {
        let position = self.frame_positions[index.value()];
        self.reader.seek(std::io::SeekFrom::Start(position))?;
        self.strategy.read_next(&mut self.reader)
    }

    /// Reads a specific frame from the trajectory by index
    /// Returns None if the index is out of bounds
    pub fn read_at(&mut self, index: usize) -> Result<Option<Frame>, CError> {
        match self.frame_index(index) {
            Some(idx) => self.read_frame(idx).map(Some),
            None => Ok(None),
        }
    }

    /// Gets a frame index if it's valid for this trajectory
    pub fn frame_index(&self, index: usize) -> Option<FrameIndex> {
        FrameIndex::new(index, self.size)
    }

    /// Returns an iterator over all frames in the trajectory
    pub fn frames(&'a mut self) -> impl Iterator<Item = Result<Frame, CError>> + 'a {
        let indices = 0..self.size;
        indices.filter_map(move |i| self.frame_index(i).map(|idx| self.read_frame(idx)))
    }

    /// Scans the file to find the position of each frame
    fn scan_all(reader: &mut BufReader<File>, strategy: &Format) -> Vec<u64> {
        let mut frame_positions = vec![0];
        while let Some(pos) = strategy.forward(reader).unwrap() {
            frame_positions.push(pos);
        }
        reader.rewind().unwrap();
        frame_positions
    }
}

impl TrajectoryWriter<'_> {
    /// Writes a frame to the trajectory
    pub fn write(&mut self, frame: &Frame) -> Result<(), CError> {
        self.strategy.write_next(&mut self.writer, frame)?;
        self.writer.flush()?;
        self.frame_count += 1;
        Ok(())
    }

    /// Finalizes the trajectory file by writing any necessary closing records
    ///
    /// This should be called when done writing to ensure the file is properly finalized.
    /// If not called explicitly, the file will be finalized when the writer is dropped.
    pub fn finish(&mut self) -> Result<(), CError> {
        self.strategy.finalize(&mut self.writer)?;
        self.writer.flush()?;
        Ok(())
    }

    /// Returns the number of frames written so far
    pub fn frame_count(&self) -> usize {
        self.frame_count
    }
}

impl Drop for TrajectoryWriter<'_> {
    fn drop(&mut self) {
        // Attempt to finalize the file when the writer is dropped.
        // Ignore errors since we can't return them from drop().
        if let Err(e) = self.finish() {
            eprintln!("Warning: Failed to finalize trajectory file: {}", e);
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
        let reader = Trajectory::open_with_format(&path, TextFormat::XYZ).unwrap();
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

        // Try to get a frame index that's too large
        let invalid_idx = reader.frame_index(4);
        assert!(invalid_idx.is_none());
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
