use thiserror::Error;

#[derive(Error, Debug)]
pub enum CError {
    #[error("Unsupported file format: `{0}`")]
    UnsupportedFileFormat(String),
    #[error("{0}")]
    IoError(#[from] std::io::Error),
    #[error("generic error: {0}")]
    GenericError(String),
    #[error("{format} format: not enough lines at step {step} (expected {expected}, got {got})")]
    UnexpectedEof {
        format: String,
        step: usize,
        expected: usize,
        got: usize,
    },
    #[error("unknown format: {0}")]
    UnknownFormat(String),
    #[error("")]
    UnexpectedSymbol,
    #[error("Failed to parse float: {0}")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("Missing token")]
    MissingToken,
}
