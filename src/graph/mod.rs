//! Core graph genome data structures
//!
//! This module defines the fundamental types for representing graph genomes,
//! following the GFA (Graphical Fragment Assembly) specification.

pub mod segment;
pub mod link;
pub mod path;

pub use segment::{Segment, TagValue};
pub use link::Link;
pub use path::{Path, PathSegment};

use std::collections::HashMap;

/// A complete graph genome representation
#[derive(Clone, Debug, Default)]
pub struct GraphGenome {
    /// Graph segments (nodes) - keyed by segment name
    pub segments: HashMap<String, Segment>,

    /// Links (edges) between segments
    pub links: Vec<Link>,

    /// Paths through the graph (haplotypes, references, etc.)
    pub paths: Vec<Path>,

    /// Optional header information
    pub header: Option<GraphHeader>,

    /// Adjacency list for efficient traversal (segment_id -> connected links)
    adjacency: HashMap<String, Vec<usize>>,
}

/// GFA header information
#[derive(Clone, Debug, Default)]
pub struct GraphHeader {
    /// GFA version (1.0, 1.1, 2.0)
    pub version: String,
    /// Optional trace spacing
    pub trace_spacing: Option<u32>,
}

impl GraphGenome {
    /// Create a new empty graph genome
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a segment to the graph
    pub fn add_segment(&mut self, segment: Segment) {
        self.segments.insert(segment.name.clone(), segment);
    }

    /// Add a link to the graph
    pub fn add_link(&mut self, link: Link) {
        let link_idx = self.links.len();

        // Update adjacency list
        self.adjacency
            .entry(link.from_segment.clone())
            .or_default()
            .push(link_idx);
        self.adjacency
            .entry(link.to_segment.clone())
            .or_default()
            .push(link_idx);

        self.links.push(link);
    }

    /// Add a path to the graph
    pub fn add_path(&mut self, path: Path) {
        self.paths.push(path);
    }

    /// Get a segment by name
    pub fn get_segment(&self, name: &str) -> Option<&Segment> {
        self.segments.get(name)
    }

    /// Get all links connected to a segment
    pub fn get_links_for_segment(&self, segment_name: &str) -> Vec<&Link> {
        self.adjacency
            .get(segment_name)
            .map(|indices| indices.iter().map(|&i| &self.links[i]).collect())
            .unwrap_or_default()
    }

    /// Get outgoing links from a segment (in a specific orientation)
    pub fn get_outgoing_links(&self, segment_name: &str, orientation: Orientation) -> Vec<&Link> {
        self.get_links_for_segment(segment_name)
            .into_iter()
            .filter(|link| {
                link.from_segment == segment_name && link.from_orient == orientation
            })
            .collect()
    }

    /// Get incoming links to a segment (in a specific orientation)
    pub fn get_incoming_links(&self, segment_name: &str, orientation: Orientation) -> Vec<&Link> {
        self.get_links_for_segment(segment_name)
            .into_iter()
            .filter(|link| {
                link.to_segment == segment_name && link.to_orient == orientation
            })
            .collect()
    }

    /// Calculate total sequence length
    pub fn total_sequence_length(&self) -> usize {
        self.segments.values().map(|s| s.sequence_length()).sum()
    }

    /// Get the number of nodes (segments)
    pub fn node_count(&self) -> usize {
        self.segments.len()
    }

    /// Get the number of edges (links)
    pub fn edge_count(&self) -> usize {
        self.links.len()
    }

    /// Find all source segments (no incoming edges)
    pub fn find_sources(&self) -> Vec<&str> {
        self.segments
            .keys()
            .filter(|name| {
                self.get_incoming_links(name, Orientation::Forward).is_empty()
                    && self.get_incoming_links(name, Orientation::Reverse).is_empty()
            })
            .map(|s| s.as_str())
            .collect()
    }

    /// Find all sink segments (no outgoing edges)
    pub fn find_sinks(&self) -> Vec<&str> {
        self.segments
            .keys()
            .filter(|name| {
                self.get_outgoing_links(name, Orientation::Forward).is_empty()
                    && self.get_outgoing_links(name, Orientation::Reverse).is_empty()
            })
            .map(|s| s.as_str())
            .collect()
    }

    /// Check if the graph is a DAG (directed acyclic graph)
    pub fn is_dag(&self) -> bool {
        // Simple cycle detection using DFS
        let mut visited = HashMap::new();
        let mut rec_stack = HashMap::new();

        for segment_name in self.segments.keys() {
            if self.has_cycle_dfs(segment_name, &mut visited, &mut rec_stack) {
                return false;
            }
        }
        true
    }

    fn has_cycle_dfs(
        &self,
        node: &str,
        visited: &mut HashMap<String, bool>,
        rec_stack: &mut HashMap<String, bool>,
    ) -> bool {
        if rec_stack.get(node).copied().unwrap_or(false) {
            return true;
        }
        if visited.get(node).copied().unwrap_or(false) {
            return false;
        }

        visited.insert(node.to_string(), true);
        rec_stack.insert(node.to_string(), true);

        for link in self.get_outgoing_links(node, Orientation::Forward) {
            if self.has_cycle_dfs(&link.to_segment, visited, rec_stack) {
                return true;
            }
        }

        rec_stack.insert(node.to_string(), false);
        false
    }
}

/// Segment orientation
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Default)]
pub enum Orientation {
    #[default]
    Forward,
    Reverse,
}

impl Orientation {
    /// Parse from GFA character
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            '+' => Some(Self::Forward),
            '-' => Some(Self::Reverse),
            _ => None,
        }
    }

    /// Convert to GFA character
    pub fn to_char(self) -> char {
        match self {
            Self::Forward => '+',
            Self::Reverse => '-',
        }
    }

    /// Get the opposite orientation
    pub fn flip(self) -> Self {
        match self {
            Self::Forward => Self::Reverse,
            Self::Reverse => Self::Forward,
        }
    }
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orientation_parsing() {
        assert_eq!(Orientation::from_char('+'), Some(Orientation::Forward));
        assert_eq!(Orientation::from_char('-'), Some(Orientation::Reverse));
        assert_eq!(Orientation::from_char('x'), None);
    }

    #[test]
    fn test_graph_basic_operations() {
        let mut graph = GraphGenome::new();

        graph.add_segment(Segment::new("s1", "ACGT"));
        graph.add_segment(Segment::new("s2", "TGCA"));
        graph.add_link(Link::new("s1", Orientation::Forward, "s2", Orientation::Forward, "0M"));

        assert_eq!(graph.node_count(), 2);
        assert_eq!(graph.edge_count(), 1);
        assert_eq!(graph.total_sequence_length(), 8);
    }
}

