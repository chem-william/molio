// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
//
// See LICENSE at the project root for full text.

use thiserror::Error;

#[derive(Error, Debug)]
pub enum CError {
    #[error("Unsupported file format: `{0}`")]
    UnsupportedFileFormat(String),
    #[error("{0}")]
    IoError(#[from] std::io::Error),
    #[error("generic error: {0}")]
    GenericError(String),
    #[error("{format} format: not enough lines (expected {expected}, got {got})")]
    UnexpectedEof {
        format: String,
        expected: usize,
        got: usize,
    },
    #[error("unknown format: {0}")]
    UnknownFormat(String),
    #[error("")]
    UnexpectedSymbol,
    #[error("Failed to parse float: {0}")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("Failed to parse integer: {0}")]
    ParseIntegerError(#[from] std::num::ParseIntError),
    #[error("Missing token")]
    MissingToken,
    #[error("{expected_record_type} record '{actual_record_type}': {reason}")]
    InvalidRecord {
        expected_record_type: String,
        actual_record_type: String,
        reason: String,
    },
}
