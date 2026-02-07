//! Graph link (edge) representation

use super::Orientation;

/// A link between two segments in the graph (corresponds to GFA L-line)
#[derive(Clone, Debug)]
pub struct Link {
    /// Source segment name
    pub from_segment: String,

    /// Source segment orientation
    pub from_orient: Orientation,

    /// Target segment name
    pub to_segment: String,

    /// Target segment orientation
    pub to_orient: Orientation,

    /// Overlap CIGAR string (e.g., "0M", "5M", "*")
    pub overlap: String,

    /// Optional edge identifier (for GFA2)
    pub edge_id: Option<String>,

    /// Optional alignment positions (for GFA2 E-lines)
    pub alignment: Option<EdgeAlignment>,
}

/// Edge alignment details (GFA2 format)
#[derive(Clone, Debug)]
pub struct EdgeAlignment {
    /// Start position on source segment
    pub from_start: usize,
    /// End position on source segment (exclusive)
    pub from_end: usize,
    /// Start position on target segment
    pub to_start: usize,
    /// End position on target segment (exclusive)
    pub to_end: usize,
}

impl Link {
    /// Create a new link
    pub fn new(
        from_segment: impl Into<String>,
        from_orient: Orientation,
        to_segment: impl Into<String>,
        to_orient: Orientation,
        overlap: impl Into<String>,
    ) -> Self {
        Self {
            from_segment: from_segment.into(),
            from_orient,
            to_segment: to_segment.into(),
            to_orient,
            overlap: overlap.into(),
            edge_id: None,
            alignment: None,
        }
    }

    /// Create a link with edge ID (GFA2)
    pub fn with_edge_id(mut self, edge_id: impl Into<String>) -> Self {
        self.edge_id = Some(edge_id.into());
        self
    }

    /// Create a link with alignment positions (GFA2)
    pub fn with_alignment(mut self, alignment: EdgeAlignment) -> Self {
        self.alignment = Some(alignment);
        self
    }

    /// Get the overlap length from CIGAR
    pub fn overlap_length(&self) -> usize {
        parse_cigar_length(&self.overlap)
    }

    /// Check if this link connects in the forward direction
    pub fn is_forward(&self) -> bool {
        self.from_orient == Orientation::Forward && self.to_orient == Orientation::Forward
    }

    /// Get the reverse of this link
    pub fn reverse(&self) -> Self {
        Self {
            from_segment: self.to_segment.clone(),
            from_orient: self.to_orient.flip(),
            to_segment: self.from_segment.clone(),
            to_orient: self.from_orient.flip(),
            overlap: self.overlap.clone(),
            edge_id: self.edge_id.clone(),
            alignment: self.alignment.as_ref().map(|a| EdgeAlignment {
                from_start: a.to_start,
                from_end: a.to_end,
                to_start: a.from_start,
                to_end: a.from_end,
            }),
        }
    }

    /// Check if this link involves a specific segment
    pub fn involves_segment(&self, segment: &str) -> bool {
        self.from_segment == segment || self.to_segment == segment
    }

    /// Get the other segment in this link
    pub fn other_segment(&self, segment: &str) -> Option<(&str, Orientation)> {
        if self.from_segment == segment {
            Some((&self.to_segment, self.to_orient))
        } else if self.to_segment == segment {
            Some((&self.from_segment, self.from_orient))
        } else {
            None
        }
    }
}

impl std::fmt::Display for Link {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "L\t{}\t{}\t{}\t{}\t{}",
            self.from_segment,
            self.from_orient,
            self.to_segment,
            self.to_orient,
            self.overlap
        )
    }
}

/// Parse a simple CIGAR string to get the total length
fn parse_cigar_length(cigar: &str) -> usize {
    if cigar == "*" {
        return 0;
    }

    let mut length = 0;
    let mut num_str = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num_str.push(c);
        } else if !num_str.is_empty() {
            if let Ok(num) = num_str.parse::<usize>() {
                // M, =, X contribute to length
                if c == 'M' || c == '=' || c == 'X' {
                    length += num;
                }
            }
            num_str.clear();
        }
    }

    length
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_link_creation() {
        let link = Link::new("s1", Orientation::Forward, "s2", Orientation::Forward, "5M");
        assert_eq!(link.from_segment, "s1");
        assert_eq!(link.to_segment, "s2");
        assert_eq!(link.overlap_length(), 5);
    }

    #[test]
    fn test_link_reverse() {
        let link = Link::new("s1", Orientation::Forward, "s2", Orientation::Reverse, "0M");
        let rev = link.reverse();
        assert_eq!(rev.from_segment, "s2");
        assert_eq!(rev.from_orient, Orientation::Forward);
        assert_eq!(rev.to_segment, "s1");
        assert_eq!(rev.to_orient, Orientation::Reverse);
    }

    #[test]
    fn test_cigar_parsing() {
        assert_eq!(parse_cigar_length("10M"), 10);
        assert_eq!(parse_cigar_length("5M2I3M"), 8);
        assert_eq!(parse_cigar_length("*"), 0);
        assert_eq!(parse_cigar_length("0M"), 0);
    }
}

