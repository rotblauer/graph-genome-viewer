//! GAF (Graph Alignment Format) file parser

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::graph::Orientation;
use super::{IoError, IoResult};

/// A single alignment from a GAF file
#[derive(Clone, Debug)]
pub struct GafAlignment {
    /// Query sequence name
    pub query_name: String,
    /// Query sequence length
    pub query_length: usize,
    /// Query start position (0-based)
    pub query_start: usize,
    /// Query end position
    pub query_end: usize,
    /// Strand relative to path (+ or -)
    pub strand: Orientation,
    /// Path string (>s1>s2<s3)
    pub path: String,
    /// Path length
    pub path_length: usize,
    /// Path start position
    pub path_start: usize,
    /// Path end position
    pub path_end: usize,
    /// Number of matching bases
    pub matches: usize,
    /// Alignment block length
    pub block_length: usize,
    /// Mapping quality
    pub mapq: u8,
    /// Optional tags
    pub tags: Vec<(String, String)>,
}

impl GafAlignment {
    /// Calculate alignment identity (matches / block_length)
    pub fn identity(&self) -> f64 {
        if self.block_length == 0 {
            0.0
        } else {
            self.matches as f64 / self.block_length as f64
        }
    }

    /// Parse the path string into segment references
    pub fn path_segments(&self) -> Vec<(String, Orientation)> {
        parse_path_string(&self.path)
    }
}

/// Parse a path string (>s1>s2<s3) into segment names and orientations
fn parse_path_string(path: &str) -> Vec<(String, Orientation)> {
    let mut segments = Vec::new();
    let mut current_name = String::new();
    let mut current_orient = Orientation::Forward;

    for c in path.chars() {
        match c {
            '>' => {
                if !current_name.is_empty() {
                    segments.push((current_name.clone(), current_orient));
                    current_name.clear();
                }
                current_orient = Orientation::Forward;
            }
            '<' => {
                if !current_name.is_empty() {
                    segments.push((current_name.clone(), current_orient));
                    current_name.clear();
                }
                current_orient = Orientation::Reverse;
            }
            _ => {
                current_name.push(c);
            }
        }
    }

    if !current_name.is_empty() {
        segments.push((current_name, current_orient));
    }

    segments
}

/// GAF file loader
pub struct GafLoader;

impl GafLoader {
    /// Load a GAF file
    pub fn load(path: &Path) -> IoResult<Vec<GafAlignment>> {
        if !path.exists() {
            return Err(IoError::FileNotFound(path.display().to_string()));
        }

        let file = File::open(path)?;
        let reader = BufReader::new(file);

        Self::parse(reader)
    }

    /// Parse GAF from a reader
    pub fn parse<R: BufRead>(reader: R) -> IoResult<Vec<GafAlignment>> {
        let mut alignments = Vec::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result?;
            let line = line.trim();

            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let alignment = Self::parse_line(line, line_num + 1)?;
            alignments.push(alignment);
        }

        Ok(alignments)
    }

    /// Parse a single GAF line
    fn parse_line(line: &str, line_num: usize) -> IoResult<GafAlignment> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 12 {
            return Err(IoError::ParseError {
                line: line_num,
                message: format!("GAF line requires at least 12 fields, got {}", fields.len()),
            });
        }

        let query_name = fields[0].to_string();
        let query_length = fields[1].parse().unwrap_or(0);
        let query_start = fields[2].parse().unwrap_or(0);
        let query_end = fields[3].parse().unwrap_or(0);

        let strand = match fields[4] {
            "+" => Orientation::Forward,
            "-" => Orientation::Reverse,
            _ => Orientation::Forward,
        };

        let path = fields[5].to_string();
        let path_length = fields[6].parse().unwrap_or(0);
        let path_start = fields[7].parse().unwrap_or(0);
        let path_end = fields[8].parse().unwrap_or(0);
        let matches = fields[9].parse().unwrap_or(0);
        let block_length = fields[10].parse().unwrap_or(0);
        let mapq = fields[11].parse().unwrap_or(0);

        // Parse optional tags
        let tags: Vec<(String, String)> = fields
            .iter()
            .skip(12)
            .filter_map(|f| {
                let parts: Vec<&str> = f.splitn(3, ':').collect();
                if parts.len() >= 3 {
                    Some((parts[0].to_string(), parts[2].to_string()))
                } else {
                    None
                }
            })
            .collect();

        Ok(GafAlignment {
            query_name,
            query_length,
            query_start,
            query_end,
            strand,
            path,
            path_length,
            path_start,
            path_end,
            matches,
            block_length,
            mapq,
            tags,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_path_string() {
        let segments = parse_path_string(">s1>s2<s3>s4");
        assert_eq!(segments.len(), 4);
        assert_eq!(segments[0], ("s1".to_string(), Orientation::Forward));
        assert_eq!(segments[2], ("s3".to_string(), Orientation::Reverse));
    }

    #[test]
    fn test_identity_calculation() {
        let alignment = GafAlignment {
            query_name: "read1".to_string(),
            query_length: 100,
            query_start: 0,
            query_end: 100,
            strand: Orientation::Forward,
            path: ">s1".to_string(),
            path_length: 100,
            path_start: 0,
            path_end: 100,
            matches: 95,
            block_length: 100,
            mapq: 60,
            tags: vec![],
        };

        assert!((alignment.identity() - 0.95).abs() < 0.001);
    }
}

