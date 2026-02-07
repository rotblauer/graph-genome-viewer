//! File I/O utilities for loading graph genomes and alignment data

pub mod gfa;
pub mod gaf;
pub mod alignment;

pub use gfa::GfaLoader;
pub use alignment::{AlignmentData, ReadAlignment};

use thiserror::Error;

/// Errors that can occur during file I/O
#[derive(Error, Debug)]
pub enum IoError {
    #[error("File not found: {0}")]
    FileNotFound(String),

    #[error("Invalid file format: {0}")]
    InvalidFormat(String),

    #[error("Parse error at line {line}: {message}")]
    ParseError { line: usize, message: String },

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Unsupported GFA version: {0}")]
    UnsupportedVersion(String),
}

/// Result type for I/O operations
pub type IoResult<T> = Result<T, IoError>;

