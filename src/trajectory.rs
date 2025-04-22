use crate::error::CError;
use crate::format::FileFormat;
use crate::format::Format;
use crate::format::TextFormat;
use crate::frame::Frame;
use std::fs::File;
use std::io::{BufReader, Seek};
use std::path::Path;

pub struct Trajectory<'a> {
    pub frames: Vec<Frame>,
    pub size: usize,

    path: &'a Path,
    strategy: Format,
    reader: BufReader<File>,
    frame_positions: Vec<u64>,
}

impl<'a> Trajectory<'a> {
    /// Constructs a `Trajectory` that will guess the format.
    pub fn new(path: &'a Path) -> Result<Self, CError> {
        Self::with_format(path, TextFormat::Guess)
    }

    /// Constructs a `Trajectory` with an explicitly provided format.
    pub fn with_format(path: &'a Path, fmt: TextFormat) -> Result<Self, CError> {
        let strategy = Format::new_from_format(fmt, path)?;

        let file = File::open(path).map_err(CError::IoError)?;
        let mut reader = BufReader::new(file);
        let frame_positions = Trajectory::scan_all(&mut reader, &strategy);
        let size = frame_positions.len() - 1;

        Ok(Trajectory {
            frames: vec![],
            path,
            strategy,
            frame_positions,
            reader,
            size,
        })
    }

    pub fn read(&mut self) -> Result<Option<Frame>, CError> {
        self.strategy.read(&mut self.reader)
    }

    pub fn read_at(&mut self, index: usize) -> Result<Frame, CError> {
        if index >= self.frame_positions.len() {
            if self.frame_positions.is_empty() {
                return Err(CError::GenericError("frame_positions is empty".to_string()));
            }
            return Err(CError::GenericError(
                "index is bigger than frame_positions".to_string(),
            ));
        }
        self.reader
            .seek(std::io::SeekFrom::Start(self.frame_positions[index]))
            .unwrap();
        self.strategy.read_next(&mut self.reader)
    }

    pub fn write(&self, frame: &Frame) -> Result<(), CError> {
        self.strategy.write(self.path, frame)
    }

    fn scan_all(reader: &mut BufReader<File>, strategy: &Format) -> Vec<u64> {
        let mut frame_positions = vec![0];
        while let Some(pos) = strategy.forward(reader).unwrap() {
            frame_positions.push(pos)
        }
        reader.rewind().unwrap();
        frame_positions
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    const TEST_XYZ_PATH: &str = "./src/tests-data/xyz/extended.xyz";

    #[test]
    fn test_trajectory_new() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let trajectory = Trajectory::new(&path).unwrap();
        assert_eq!(trajectory.size, 3);
    }

    #[test]
    fn test_trajectory_with_format() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let trajectory = Trajectory::with_format(&path, TextFormat::XYZ).unwrap();
        assert_eq!(trajectory.size, 3);
    }

    #[test]
    fn test_read_at() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let mut trajectory = Trajectory::new(&path).unwrap();

        // Read first frame
        let frame0 = trajectory.read_at(0).unwrap();
        assert_eq!(frame0.size(), 192);

        // Read second frame
        let frame1 = trajectory.read_at(1).unwrap();
        assert_eq!(frame1.size(), 62);
    }

    #[test]
    fn test_read_at_invalid_index() {
        let path = PathBuf::from(TEST_XYZ_PATH);
        let mut trajectory = Trajectory::new(&path).unwrap();

        // Try to read frame at invalid index
        let result = trajectory.read_at(4);
        assert!(result.is_err());
    }
}
