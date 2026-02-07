//! GFA file parser
//!
//! Supports both GFA1 and GFA2 formats for graph genome files.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::graph::{GraphGenome, GraphHeader, Segment, Link, Orientation, Path as GraphPath, PathSegment, TagValue};
use super::{IoError, IoResult};

/// GFA file loader
pub struct GfaLoader;

impl GfaLoader {
    /// Load a GFA file from path
    pub fn load(path: &Path) -> IoResult<GraphGenome> {
        if !path.exists() {
            return Err(IoError::FileNotFound(path.display().to_string()));
        }

        let file = File::open(path)?;
        let reader = BufReader::new(file);

        Self::parse(reader)
    }

    /// Parse GFA from a reader
    pub fn parse<R: BufRead>(reader: R) -> IoResult<GraphGenome> {
        let mut graph = GraphGenome::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line = line_result?;
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.is_empty() {
                continue;
            }

            match fields[0] {
                "H" => {
                    // Header line
                    graph.header = Some(Self::parse_header(&fields)?);
                }
                "S" => {
                    // Segment line
                    let segment = Self::parse_segment(&fields, line_num + 1)?;
                    graph.add_segment(segment);
                }
                "L" => {
                    // Link line
                    let link = Self::parse_link(&fields, line_num + 1)?;
                    graph.add_link(link);
                }
                "P" => {
                    // Path line
                    let path = Self::parse_path(&fields, line_num + 1)?;
                    graph.add_path(path);
                }
                "W" => {
                    // Walk line (GFA 1.1+)
                    let path = Self::parse_walk(&fields, line_num + 1)?;
                    graph.add_path(path);
                }
                "E" => {
                    // Edge line (GFA2)
                    let link = Self::parse_edge(&fields, line_num + 1)?;
                    graph.add_link(link);
                }
                "C" | "F" | "G" | "U" | "O" => {
                    // GFA2 specific lines - log and skip for now
                    log::debug!("Skipping GFA2 {} line at {}", fields[0], line_num + 1);
                }
                _ => {
                    // Unknown line type
                    log::warn!("Unknown line type '{}' at line {}", fields[0], line_num + 1);
                }
            }
        }

        Ok(graph)
    }

    /// Parse header line
    fn parse_header(fields: &[&str]) -> IoResult<GraphHeader> {
        let mut header = GraphHeader::default();

        for field in fields.iter().skip(1) {
            if let Some(version) = field.strip_prefix("VN:Z:") {
                header.version = version.to_string();
            } else if let Some(ts) = field.strip_prefix("TS:i:") {
                header.trace_spacing = ts.parse().ok();
            }
        }

        Ok(header)
    }

    /// Parse segment line: S <name> <sequence> [tags...]
    fn parse_segment(fields: &[&str], line_num: usize) -> IoResult<Segment> {
        if fields.len() < 3 {
            return Err(IoError::ParseError {
                line: line_num,
                message: "Segment line requires at least 3 fields".to_string(),
            });
        }

        let name = fields[1];
        let sequence = fields[2];

        let mut segment = Segment::new(name, sequence);

        // Parse optional tags
        for field in fields.iter().skip(3) {
            if let Some((key, value)) = Self::parse_tag(field) {
                segment.add_tag(key, value);
            }
        }

        // Check for LN tag (length) if sequence is missing
        if sequence == "*" {
            if let Some(TagValue::Integer(len)) = segment.get_tag("LN") {
                segment.length = Some(*len as usize);
            }
        }

        Ok(segment)
    }

    /// Parse link line: L <from> <from_orient> <to> <to_orient> <overlap>
    fn parse_link(fields: &[&str], line_num: usize) -> IoResult<Link> {
        if fields.len() < 6 {
            return Err(IoError::ParseError {
                line: line_num,
                message: "Link line requires at least 6 fields".to_string(),
            });
        }

        let from_segment = fields[1];
        let from_orient = Orientation::from_char(fields[2].chars().next().unwrap_or('+'))
            .ok_or_else(|| IoError::ParseError {
                line: line_num,
                message: format!("Invalid orientation: {}", fields[2]),
            })?;
        let to_segment = fields[3];
        let to_orient = Orientation::from_char(fields[4].chars().next().unwrap_or('+'))
            .ok_or_else(|| IoError::ParseError {
                line: line_num,
                message: format!("Invalid orientation: {}", fields[4]),
            })?;
        let overlap = fields[5];

        Ok(Link::new(from_segment, from_orient, to_segment, to_orient, overlap))
    }

    /// Parse path line: P <name> <segment_names> <overlaps>
    fn parse_path(fields: &[&str], line_num: usize) -> IoResult<GraphPath> {
        if fields.len() < 3 {
            return Err(IoError::ParseError {
                line: line_num,
                message: "Path line requires at least 3 fields".to_string(),
            });
        }

        let name = fields[1];
        let segment_str = fields[2];

        let segments: Vec<PathSegment> = segment_str
            .split(',')
            .filter_map(|s| PathSegment::parse(s.trim()))
            .collect();

        let mut path = GraphPath::with_segments(name, segments);

        // Parse overlaps if present
        if fields.len() > 3 && fields[3] != "*" {
            path.overlaps = fields[3].split(',').map(|s| s.to_string()).collect();
        }

        Ok(path)
    }

    /// Parse walk line (GFA 1.1): W <sample> <haplotype> <seqname> <start> <end> <walk>
    fn parse_walk(fields: &[&str], line_num: usize) -> IoResult<GraphPath> {
        if fields.len() < 7 {
            return Err(IoError::ParseError {
                line: line_num,
                message: "Walk line requires at least 7 fields".to_string(),
            });
        }

        let sample = fields[1];
        let haplotype: u32 = fields[2].parse().unwrap_or(0);
        let seq_name = fields[3];
        let start: usize = fields[4].parse().unwrap_or(0);
        let end: usize = fields[5].parse().unwrap_or(0);
        let walk_str = fields[6];

        // Parse walk string (>seg1<seg2>seg3 format)
        let segments = Self::parse_walk_string(walk_str);

        let name = format!("{}#{}#{}", sample, haplotype, seq_name);
        let path = GraphPath::with_segments(name, segments)
            .set_walk_metadata(sample, haplotype, seq_name, start, end);

        Ok(path)
    }

    /// Parse walk string (>seg1<seg2>seg3)
    fn parse_walk_string(walk: &str) -> Vec<PathSegment> {
        let mut segments = Vec::new();
        let mut current_name = String::new();
        let mut current_orient = Orientation::Forward;

        for c in walk.chars() {
            match c {
                '>' => {
                    if !current_name.is_empty() {
                        segments.push(PathSegment::new(&current_name, current_orient));
                        current_name.clear();
                    }
                    current_orient = Orientation::Forward;
                }
                '<' => {
                    if !current_name.is_empty() {
                        segments.push(PathSegment::new(&current_name, current_orient));
                        current_name.clear();
                    }
                    current_orient = Orientation::Reverse;
                }
                _ => {
                    current_name.push(c);
                }
            }
        }

        // Don't forget the last segment
        if !current_name.is_empty() {
            segments.push(PathSegment::new(&current_name, current_orient));
        }

        segments
    }

    /// Parse GFA2 edge line: E <id> <from><orient> <to><orient> <from_start> <from_end> <to_start> <to_end> <alignment>
    fn parse_edge(fields: &[&str], line_num: usize) -> IoResult<Link> {
        if fields.len() < 9 {
            return Err(IoError::ParseError {
                line: line_num,
                message: "Edge line requires at least 9 fields".to_string(),
            });
        }

        let edge_id = fields[1];

        // Parse from reference (segment+orientation)
        let (from_segment, from_orient) = Self::parse_ref(fields[2])?;
        let (to_segment, to_orient) = Self::parse_ref(fields[3])?;

        let alignment = if fields.len() > 8 { fields[8] } else { "*" };

        Ok(Link::new(from_segment, from_orient, to_segment, to_orient, alignment)
            .with_edge_id(edge_id))
    }

    /// Parse a GFA2 reference (segment+orient, e.g., "seg1+")
    fn parse_ref(s: &str) -> IoResult<(String, Orientation)> {
        if s.is_empty() {
            return Err(IoError::InvalidFormat("Empty reference".to_string()));
        }

        let last = s.chars().last().unwrap();
        let orient = Orientation::from_char(last).unwrap_or(Orientation::Forward);
        let name = if last == '+' || last == '-' {
            &s[..s.len() - 1]
        } else {
            s
        };

        Ok((name.to_string(), orient))
    }

    /// Parse a GFA tag (e.g., "LN:i:100", "RC:i:5", "DP:f:1.5")
    fn parse_tag(field: &str) -> Option<(String, TagValue)> {
        let parts: Vec<&str> = field.splitn(3, ':').collect();
        if parts.len() != 3 {
            return None;
        }

        let key = parts[0];
        let type_char = parts[1];
        let value_str = parts[2];

        let value = match type_char {
            "A" => value_str.chars().next().map(TagValue::Char)?,
            "i" => value_str.parse().ok().map(TagValue::Integer)?,
            "f" => value_str.parse().ok().map(TagValue::Float)?,
            "Z" | "H" => TagValue::String(value_str.to_string()),
            _ => TagValue::String(value_str.to_string()),
        };

        Some((key.to_string(), value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_simple_gfa() {
        let gfa = r#"H	VN:Z:1.0
S	s1	ACGT
S	s2	TGCA
L	s1	+	s2	+	0M
P	path1	s1+,s2+	*
"#;

        let reader = BufReader::new(Cursor::new(gfa));
        let graph = GfaLoader::parse(reader).unwrap();

        assert_eq!(graph.segments.len(), 2);
        assert_eq!(graph.links.len(), 1);
        assert_eq!(graph.paths.len(), 1);
    }

    #[test]
    fn test_parse_walk_string() {
        let segments = GfaLoader::parse_walk_string(">s1>s2<s3>s4");
        assert_eq!(segments.len(), 4);
        assert_eq!(segments[0].name, "s1");
        assert_eq!(segments[0].orientation, Orientation::Forward);
        assert_eq!(segments[2].name, "s3");
        assert_eq!(segments[2].orientation, Orientation::Reverse);
    }
}

