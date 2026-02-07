//! Goodness-of-fit metrics for evaluating alignments against the graph

use crate::graph::GraphGenome;
use crate::io::AlignmentData;
use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
pub struct FitMetrics {
    pub overall_score: f64,
    pub segment_scores: HashMap<String, SegmentFit>,
    pub problematic_paths: Vec<String>,
    pub uncovered_segments: Vec<String>,
    pub high_coverage_segments: Vec<String>,
}

#[derive(Clone, Debug, Default)]
pub struct SegmentFit {
    pub expected_coverage: f32,
    pub observed_coverage: f32,
    pub coverage_ratio: f32,
    pub mean_identity: f64,
    pub score: f64,
}

impl FitMetrics {
    pub fn compute(graph: &GraphGenome, alignments: &AlignmentData) -> Self {
        let mut metrics = Self::default();

        if graph.segments.is_empty() {
            return metrics;
        }

        let total_reads = alignments.alignments.len() as f32;
        let num_segments = graph.segments.len() as f32;
        let expected_per_segment = (total_reads / num_segments).max(1.0);

        let mut total_score = 0.0;
        let mut scored_segments = 0;

        for (name, _segment) in &graph.segments {
            let observed = alignments.segment_coverage
                .get(name)
                .map(|c| c.read_count as f32)
                .unwrap_or(0.0);

            let coverage_ratio = observed / expected_per_segment;

            let coverage_score = if coverage_ratio < 0.1 {
                coverage_ratio * 5.0
            } else if coverage_ratio > 5.0 {
                1.0 / coverage_ratio
            } else {
                1.0 - (1.0 - coverage_ratio).abs().min(1.0) * 0.5
            };

            let segment_alignments: Vec<_> = alignments.alignments.iter()
                .filter(|a| a.path.iter().any(|(s, _)| s == name))
                .collect();

            let mean_identity = if segment_alignments.is_empty() {
                0.0
            } else {
                segment_alignments.iter().map(|a| a.identity).sum::<f64>()
                    / segment_alignments.len() as f64
            };

            let segment_score = (coverage_score as f64 * 0.4 + mean_identity * 0.6).clamp(0.0, 1.0);

            let segment_fit = SegmentFit {
                expected_coverage: expected_per_segment,
                observed_coverage: observed,
                coverage_ratio,
                mean_identity,
                score: segment_score,
            };

            metrics.segment_scores.insert(name.clone(), segment_fit);
            total_score += segment_score;
            scored_segments += 1;

            if observed < 0.1 {
                metrics.uncovered_segments.push(name.clone());
            } else if coverage_ratio > 5.0 {
                metrics.high_coverage_segments.push(name.clone());
            }
        }

        metrics.overall_score = if scored_segments > 0 {
            total_score / scored_segments as f64
        } else {
            0.0
        };

        metrics
    }

    pub fn summary(&self) -> String {
        format!(
            "Overall Fit: {:.1}%\nUncovered: {}\nHigh Coverage: {}",
            self.overall_score * 100.0,
            self.uncovered_segments.len(),
            self.high_coverage_segments.len()
        )
    }
}

