//! Alignment data structures for long-read sequencing data

use crate::graph::Orientation;
use crate::io::gaf::GafAlignment;

/// Container for alignment data
#[derive(Clone, Debug, Default)]
pub struct AlignmentData {
    /// Individual read alignments
    pub alignments: Vec<ReadAlignment>,

    /// Per-segment coverage statistics
    pub segment_coverage: std::collections::HashMap<String, SegmentCoverage>,
}

/// A single read alignment
#[derive(Clone, Debug)]
pub struct ReadAlignment {
    /// Read/query name
    pub name: String,
    /// Read length
    pub length: usize,
    /// Path through the graph (segment names with orientations)
    pub path: Vec<(String, Orientation)>,
    /// Alignment identity (0.0-1.0)
    pub identity: f64,
    /// Mapping quality
    pub mapq: u8,
    /// Start position on read
    pub read_start: usize,
    /// End position on read
    pub read_end: usize,
    /// Start position on path
    pub path_start: usize,
    /// End position on path
    pub path_end: usize,
}

/// Coverage statistics for a segment
#[derive(Clone, Debug, Default)]
pub struct SegmentCoverage {
    /// Number of reads covering this segment
    pub read_count: usize,
    /// Total bases covering this segment
    pub total_bases: usize,
    /// Average coverage depth
    pub depth: f32,
    /// Minimum coverage at any position
    pub min_coverage: u32,
    /// Maximum coverage at any position
    pub max_coverage: u32,
}

impl AlignmentData {
    /// Create new empty alignment data
    pub fn new() -> Self {
        Self::default()
    }

    /// Load from GAF alignments
    pub fn from_gaf(gaf_alignments: Vec<GafAlignment>) -> Self {
        let mut data = Self::new();

        for gaf in gaf_alignments {
            let path = gaf.path_segments();
            let alignment = ReadAlignment {
                name: gaf.query_name.clone(),
                length: gaf.query_length,
                path,
                identity: gaf.identity(),
                mapq: gaf.mapq,
                read_start: gaf.query_start,
                read_end: gaf.query_end,
                path_start: gaf.path_start,
                path_end: gaf.path_end,
            };

            // Update segment coverage
            for (segment_name, _) in &alignment.path {
                let coverage = data.segment_coverage
                    .entry(segment_name.clone())
                    .or_default();
                coverage.read_count += 1;
            }

            data.alignments.push(alignment);
        }

        data
    }

    /// Get alignments that pass through a specific segment
    pub fn alignments_through_segment(&self, segment: &str) -> Vec<&ReadAlignment> {
        self.alignments
            .iter()
            .filter(|a| a.path.iter().any(|(s, _)| s == segment))
            .collect()
    }

    /// Get coverage for a segment
    pub fn get_coverage(&self, segment: &str) -> Option<&SegmentCoverage> {
        self.segment_coverage.get(segment)
    }

    /// Calculate overall alignment statistics
    pub fn overall_stats(&self) -> AlignmentSummary {
        if self.alignments.is_empty() {
            return AlignmentSummary::default();
        }

        let total_reads = self.alignments.len();
        let total_bases: usize = self.alignments.iter().map(|a| a.length).sum();
        let mean_identity: f64 = self.alignments.iter().map(|a| a.identity).sum::<f64>() / total_reads as f64;
        let mean_mapq: f32 = self.alignments.iter().map(|a| a.mapq as f32).sum::<f32>() / total_reads as f32;

        AlignmentSummary {
            total_reads,
            total_bases,
            mean_identity,
            mean_mapq,
            segments_covered: self.segment_coverage.len(),
        }
    }
}

/// Summary statistics for alignments
#[derive(Clone, Debug, Default)]
pub struct AlignmentSummary {
    pub total_reads: usize,
    pub total_bases: usize,
    pub mean_identity: f64,
    pub mean_mapq: f32,
    pub segments_covered: usize,
}

