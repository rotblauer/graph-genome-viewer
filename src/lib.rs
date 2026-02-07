//! Graph Genome Viewer Library
//!
//! This library provides core functionality for parsing and manipulating
//! graph genomes in GFA format.
//!
//! # ⚠️ AI-Generated Prototype
//!
//! This code was generated with AI assistance and is an **experimental prototype**.
//! **NOT FOR PRODUCTION USE** - Not validated for scientific accuracy.
//! Use at your own risk for research/educational purposes only.

pub mod graph;
pub mod io;
pub mod layout;
pub mod analysis;

// Re-export commonly used types
pub use graph::{GraphGenome, Segment, Link, Path, PathSegment, Orientation, TagValue, SegmentSequence};
pub use io::gfa::GfaLoader;


