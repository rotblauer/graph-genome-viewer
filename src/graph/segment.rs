//! Graph segment (node) representation

use std::collections::HashMap;

/// A segment in the graph genome (corresponds to GFA S-line)
#[derive(Clone, Debug)]
pub struct Segment {
    /// Unique segment identifier
    pub name: String,

    /// DNA sequence (may be empty or "*" for placeholder)
    pub sequence: SegmentSequence,

    /// Optional segment length (for when sequence is not stored)
    pub length: Option<usize>,

    /// Optional tags (key-value pairs)
    pub tags: HashMap<String, TagValue>,

    /// Computed GC content (cached)
    gc_content: Option<f32>,
}

/// Segment sequence representation
#[derive(Clone, Debug)]
pub enum SegmentSequence {
    /// Full sequence stored
    Sequence(String),
    /// Sequence not stored (placeholder)
    Missing,
}

impl SegmentSequence {
    /// Check if sequence is present
    pub fn is_present(&self) -> bool {
        matches!(self, Self::Sequence(_))
    }

    /// Get sequence string if present
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Self::Sequence(s) => Some(s),
            Self::Missing => None,
        }
    }

    /// Get sequence length
    pub fn len(&self) -> usize {
        match self {
            Self::Sequence(s) => s.len(),
            Self::Missing => 0,
        }
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// GFA tag value types
#[derive(Clone, Debug)]
pub enum TagValue {
    /// Character
    Char(char),
    /// Signed integer
    Integer(i64),
    /// Float
    Float(f32),
    /// String/printable string
    String(String),
    /// Hex byte array
    Hex(Vec<u8>),
    /// Array of values
    Array(Vec<TagValue>),
}

impl Segment {
    /// Create a new segment with sequence
    pub fn new(name: impl Into<String>, sequence: impl Into<String>) -> Self {
        let seq_string = sequence.into();
        let sequence = if seq_string == "*" {
            SegmentSequence::Missing
        } else {
            SegmentSequence::Sequence(seq_string)
        };

        Self {
            name: name.into(),
            sequence,
            length: None,
            tags: HashMap::new(),
            gc_content: None,
        }
    }

    /// Create a segment with only length (no sequence)
    pub fn with_length(name: impl Into<String>, length: usize) -> Self {
        Self {
            name: name.into(),
            sequence: SegmentSequence::Missing,
            length: Some(length),
            tags: HashMap::new(),
            gc_content: None,
        }
    }

    /// Get the sequence length
    pub fn sequence_length(&self) -> usize {
        self.length.unwrap_or_else(|| self.sequence.len())
    }

    /// Calculate GC content (percentage of G and C bases)
    pub fn gc_content(&mut self) -> f32 {
        if let Some(gc) = self.gc_content {
            return gc;
        }

        let gc = match &self.sequence {
            SegmentSequence::Sequence(seq) => {
                if seq.is_empty() {
                    0.0
                } else {
                    let gc_count = seq
                        .chars()
                        .filter(|&c| c == 'G' || c == 'g' || c == 'C' || c == 'c')
                        .count();
                    (gc_count as f32) / (seq.len() as f32) * 100.0
                }
            }
            SegmentSequence::Missing => 0.0,
        };

        self.gc_content = Some(gc);
        gc
    }

    /// Get a subsequence
    pub fn subsequence(&self, start: usize, end: usize) -> Option<&str> {
        match &self.sequence {
            SegmentSequence::Sequence(seq) => seq.get(start..end),
            SegmentSequence::Missing => None,
        }
    }

    /// Add a tag
    pub fn add_tag(&mut self, key: impl Into<String>, value: TagValue) {
        self.tags.insert(key.into(), value);
    }

    /// Get a tag value
    pub fn get_tag(&self, key: &str) -> Option<&TagValue> {
        self.tags.get(key)
    }

    /// Get the depth/coverage tag (DP or KC) if present
    pub fn coverage(&self) -> Option<f32> {
        // Try DP (depth) tag first
        if let Some(TagValue::Integer(dp)) = self.tags.get("DP") {
            return Some(*dp as f32);
        }
        if let Some(TagValue::Float(dp)) = self.tags.get("DP") {
            return Some(*dp);
        }
        // Try KC (k-mer count) tag
        if let Some(TagValue::Integer(kc)) = self.tags.get("KC") {
            return Some(*kc as f32);
        }
        None
    }
}

impl std::fmt::Display for Segment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "S\t{}\t", self.name)?;
        match &self.sequence {
            SegmentSequence::Sequence(seq) => write!(f, "{}", seq)?,
            SegmentSequence::Missing => write!(f, "*")?,
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_creation() {
        let seg = Segment::new("s1", "ACGTACGT");
        assert_eq!(seg.name, "s1");
        assert_eq!(seg.sequence_length(), 8);
    }

    #[test]
    fn test_gc_content() {
        let mut seg = Segment::new("s1", "GCGCGC");
        assert!((seg.gc_content() - 100.0).abs() < 0.01);

        let mut seg2 = Segment::new("s2", "ATATAT");
        assert!((seg2.gc_content() - 0.0).abs() < 0.01);

        let mut seg3 = Segment::new("s3", "ACGT");
        assert!((seg3.gc_content() - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_missing_sequence() {
        let seg = Segment::new("s1", "*");
        assert!(!seg.sequence.is_present());
    }
}

