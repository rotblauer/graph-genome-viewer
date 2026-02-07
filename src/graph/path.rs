//! Graph path representation

use super::Orientation;

/// A path through the graph (corresponds to GFA P-line or W-line)
#[derive(Clone, Debug)]
pub struct Path {
    /// Path name/identifier
    pub name: String,

    /// Ordered list of segment references with orientations
    pub segments: Vec<PathSegment>,

    /// Optional overlaps between consecutive segments
    pub overlaps: Vec<String>,

    /// Is this a circular path?
    pub circular: bool,

    /// Optional sample name (for W-lines)
    pub sample_name: Option<String>,

    /// Optional haplotype index (for W-lines)
    pub haplotype: Option<u32>,

    /// Optional sequence name (for W-lines)
    pub sequence_name: Option<String>,

    /// Optional start position (for W-lines)
    pub start: Option<usize>,

    /// Optional end position (for W-lines)
    pub end: Option<usize>,
}

/// A reference to a segment within a path
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathSegment {
    /// Segment name
    pub name: String,
    /// Orientation in the path
    pub orientation: Orientation,
}

impl PathSegment {
    /// Create a new path segment reference
    pub fn new(name: impl Into<String>, orientation: Orientation) -> Self {
        Self {
            name: name.into(),
            orientation,
        }
    }

    /// Parse from GFA format (e.g., "segment+" or "segment-")
    pub fn parse(s: &str) -> Option<Self> {
        if s.is_empty() {
            return None;
        }

        let last_char = s.chars().last()?;
        let orientation = Orientation::from_char(last_char)?;
        let name = &s[..s.len() - 1];

        Some(Self {
            name: name.to_string(),
            orientation,
        })
    }
}

impl std::fmt::Display for PathSegment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.name, self.orientation)
    }
}

impl Path {
    /// Create a new path
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            segments: Vec::new(),
            overlaps: Vec::new(),
            circular: false,
            sample_name: None,
            haplotype: None,
            sequence_name: None,
            start: None,
            end: None,
        }
    }

    /// Create a path with segments
    pub fn with_segments(name: impl Into<String>, segments: Vec<PathSegment>) -> Self {
        Self {
            name: name.into(),
            segments,
            overlaps: Vec::new(),
            circular: false,
            sample_name: None,
            haplotype: None,
            sequence_name: None,
            start: None,
            end: None,
        }
    }

    /// Add a segment to the path
    pub fn add_segment(&mut self, segment: PathSegment) {
        self.segments.push(segment);
    }

    /// Add a segment by name and orientation
    pub fn add(&mut self, name: impl Into<String>, orientation: Orientation) {
        self.segments.push(PathSegment::new(name, orientation));
    }

    /// Get the number of segments in the path
    pub fn len(&self) -> usize {
        self.segments.len()
    }

    /// Check if path is empty
    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    /// Check if a segment is in this path
    pub fn contains_segment(&self, segment_name: &str) -> bool {
        self.segments.iter().any(|s| s.name == segment_name)
    }

    /// Get the index of a segment in this path (first occurrence)
    pub fn segment_index(&self, segment_name: &str) -> Option<usize> {
        self.segments.iter().position(|s| s.name == segment_name)
    }

    /// Get all unique segment names in this path
    pub fn unique_segments(&self) -> Vec<&str> {
        let mut seen = std::collections::HashSet::new();
        self.segments
            .iter()
            .filter(|s| seen.insert(s.name.as_str()))
            .map(|s| s.name.as_str())
            .collect()
    }

    /// Set as circular path
    pub fn set_circular(mut self, circular: bool) -> Self {
        self.circular = circular;
        self
    }

    /// Set walk metadata (GFA W-line)
    pub fn set_walk_metadata(
        mut self,
        sample_name: impl Into<String>,
        haplotype: u32,
        sequence_name: impl Into<String>,
        start: usize,
        end: usize,
    ) -> Self {
        self.sample_name = Some(sample_name.into());
        self.haplotype = Some(haplotype);
        self.sequence_name = Some(sequence_name.into());
        self.start = Some(start);
        self.end = Some(end);
        self
    }

    /// Get a human-readable description
    pub fn description(&self) -> String {
        let mut desc = self.name.clone();
        if let Some(ref sample) = self.sample_name {
            desc = format!("{} ({})", desc, sample);
        }
        if self.circular {
            desc.push_str(" [circular]");
        }
        desc
    }
}

impl std::fmt::Display for Path {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "P\t{}\t", self.name)?;

        let segments: Vec<String> = self.segments.iter().map(|s| s.to_string()).collect();
        write!(f, "{}", segments.join(","))?;

        if !self.overlaps.is_empty() {
            write!(f, "\t{}", self.overlaps.join(","))?;
        } else {
            write!(f, "\t*")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_segment_parse() {
        let ps = PathSegment::parse("segment1+").unwrap();
        assert_eq!(ps.name, "segment1");
        assert_eq!(ps.orientation, Orientation::Forward);

        let ps2 = PathSegment::parse("seg-").unwrap();
        assert_eq!(ps2.name, "seg");
        assert_eq!(ps2.orientation, Orientation::Reverse);
    }

    #[test]
    fn test_path_creation() {
        let mut path = Path::new("path1");
        path.add("s1", Orientation::Forward);
        path.add("s2", Orientation::Reverse);
        path.add("s3", Orientation::Forward);

        assert_eq!(path.len(), 3);
        assert!(path.contains_segment("s2"));
        assert!(!path.contains_segment("s4"));
    }
}

