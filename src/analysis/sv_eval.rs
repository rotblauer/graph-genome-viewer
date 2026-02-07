//! Structural variant representation evaluation
//!
//! Evaluates how well structural variants are represented in a graph genome,
//! including bubble detection, path analysis, and coverage metrics.

use crate::graph::{GraphGenome, Orientation, Segment, TagValue};
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Region of interest for SV evaluation
#[derive(Clone, Debug)]
pub struct GenomicRegion {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl GenomicRegion {
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
        }
    }

    /// Parse region from string like "chr17:10984564-10993960"
    pub fn parse(region_str: &str) -> Option<Self> {
        let parts: Vec<&str> = region_str.split(':').collect();
        if parts.len() != 2 {
            return None;
        }
        let chrom = parts[0].to_string();
        let coords: Vec<&str> = parts[1].split('-').collect();
        if coords.len() != 2 {
            return None;
        }
        let start = coords[0].parse().ok()?;
        let end = coords[1].parse().ok()?;
        Some(Self { chrom, start, end })
    }

    pub fn length(&self) -> u64 {
        self.end - self.start
    }
}

impl fmt::Display for GenomicRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

/// Detected structural variant in the graph
#[derive(Clone, Debug)]
pub struct DetectedSV {
    pub sv_type: SVType,
    pub ref_segments: Vec<String>,
    pub alt_segments: Vec<String>,
    pub ref_length: usize,
    pub alt_length: usize,
    pub size_diff: i64,
    pub samples_with_ref: Vec<String>,
    pub samples_with_alt: Vec<String>,
}

/// Types of structural variants
#[derive(Clone, Debug, PartialEq)]
pub enum SVType {
    Deletion,
    Insertion,
    Inversion,
    Complex,
    SNV,
    Unknown,
}

impl fmt::Display for SVType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SVType::Deletion => write!(f, "DEL"),
            SVType::Insertion => write!(f, "INS"),
            SVType::Inversion => write!(f, "INV"),
            SVType::Complex => write!(f, "COMPLEX"),
            SVType::SNV => write!(f, "SNV"),
            SVType::Unknown => write!(f, "UNKNOWN"),
        }
    }
}

/// A bubble structure in the graph (represents an SV)
#[derive(Clone, Debug)]
pub struct Bubble {
    pub start_segment: String,
    pub end_segment: String,
    pub paths: Vec<BubblePath>,
}

#[derive(Clone, Debug)]
pub struct BubblePath {
    pub segments: Vec<(String, Orientation)>,
    pub total_length: usize,
    pub is_reference: bool,
}

/// Complete SV evaluation report
#[derive(Clone, Debug, Default)]
pub struct SVEvalReport {
    pub region: Option<GenomicRegion>,
    pub total_segments: usize,
    pub total_links: usize,
    pub total_paths: usize,
    pub reference_path_length: usize,
    pub detected_svs: Vec<DetectedSV>,
    pub bubbles: Vec<Bubble>,
    pub samples: Vec<SampleSummary>,
    pub graph_complexity: GraphComplexity,
    pub warnings: Vec<String>,
}

#[derive(Clone, Debug, Default)]
pub struct SampleSummary {
    pub name: String,
    pub haplotype: u32,
    pub path_length: usize,
    pub segment_count: usize,
    pub sv_count: usize,
}

#[derive(Clone, Debug, Default)]
pub struct GraphComplexity {
    pub max_degree: usize,
    pub avg_degree: f64,
    pub branching_segments: usize,
    pub is_dag: bool,
}

/// SV representation evaluator
pub struct SVEvaluator<'a> {
    graph: &'a GraphGenome,
    region: Option<GenomicRegion>,
}

impl<'a> SVEvaluator<'a> {
    pub fn new(graph: &'a GraphGenome) -> Self {
        Self {
            graph,
            region: None,
        }
    }

    pub fn with_region(mut self, region: GenomicRegion) -> Self {
        self.region = Some(region);
        self
    }

    /// Run full evaluation and produce a report
    pub fn evaluate(&self) -> SVEvalReport {
        let mut report = SVEvalReport::default();
        report.region = self.region.clone();

        // Basic graph stats
        report.total_segments = self.graph.segments.len();
        report.total_links = self.graph.links.len();
        report.total_paths = self.graph.paths.len();

        // Detect SVs from segment tags (VT tags)
        report.detected_svs = self.detect_svs_from_tags();

        // Find bubbles in the graph
        report.bubbles = self.find_bubbles();

        // Analyze paths/samples
        report.samples = self.analyze_samples();

        // Compute reference path length
        report.reference_path_length = self.compute_reference_length();

        // Compute graph complexity
        report.graph_complexity = self.compute_complexity();

        // Check for issues
        report.warnings = self.check_for_issues();

        report
    }

    /// Detect SVs by looking at segment VT (variant type) tags
    fn detect_svs_from_tags(&self) -> Vec<DetectedSV> {
        let mut svs = Vec::new();

        for (name, segment) in &self.graph.segments {
            if let Some(TagValue::String(vt)) = segment.get_tag("VT") {
                let sv_type = match vt.as_str() {
                    "DEL" => SVType::Deletion,
                    "INS" => SVType::Insertion,
                    "INV" => SVType::Inversion,
                    "COMPLEX" => SVType::Complex,
                    "SNV" => SVType::SNV,
                    _ => SVType::Unknown,
                };

                let alt_length = segment.sequence_length();

                // Find corresponding reference segments by looking at paths
                let (ref_segs, ref_len) = self.find_reference_equivalent(name);
                let (samples_ref, samples_alt) = self.find_samples_for_variant(name, &ref_segs);

                svs.push(DetectedSV {
                    sv_type,
                    ref_segments: ref_segs,
                    alt_segments: vec![name.clone()],
                    ref_length: ref_len,
                    alt_length,
                    size_diff: alt_length as i64 - ref_len as i64,
                    samples_with_ref: samples_ref,
                    samples_with_alt: samples_alt,
                });
            }
        }

        svs
    }

    /// Find reference segments that correspond to a variant segment
    fn find_reference_equivalent(&self, _variant_seg: &str) -> (Vec<String>, usize) {
        // Look for segments with "ref" in name or on the GRCh38 path
        let mut ref_segs = Vec::new();
        let mut ref_len = 0;

        for path in &self.graph.paths {
            if path.name.contains("GRCh38") || path.name.contains("grch38") {
                for seg in &path.segments {
                    if seg.name.contains("var_ref") {
                        ref_segs.push(seg.name.clone());
                        if let Some(segment) = self.graph.get_segment(&seg.name) {
                            ref_len += segment.sequence_length();
                        }
                    }
                }
                break;
            }
        }

        (ref_segs, ref_len)
    }

    /// Find which samples carry reference vs alternate alleles
    fn find_samples_for_variant(&self, variant_seg: &str, ref_segs: &[String]) -> (Vec<String>, Vec<String>) {
        let mut samples_ref = Vec::new();
        let mut samples_alt = Vec::new();

        for path in &self.graph.paths {
            let path_segments: HashSet<_> = path.segments.iter().map(|s| &s.name).collect();

            let has_variant = path_segments.contains(&variant_seg.to_string());
            let has_all_refs = !ref_segs.is_empty() && ref_segs.iter().all(|r| path_segments.contains(r));

            let sample_name = format!("{}_{}", path.name, path.haplotype.unwrap_or(0));

            if has_variant {
                samples_alt.push(sample_name);
            } else if has_all_refs {
                samples_ref.push(sample_name);
            }
        }

        (samples_ref, samples_alt)
    }

    /// Find bubble structures (branching points where paths diverge and reconverge)
    fn find_bubbles(&self) -> Vec<Bubble> {
        let mut bubbles = Vec::new();

        // Build in-degree and out-degree maps
        let mut out_degree: HashMap<String, usize> = HashMap::new();
        let mut in_degree: HashMap<String, usize> = HashMap::new();

        for link in &self.graph.links {
            *out_degree.entry(link.from_segment.clone()).or_default() += 1;
            *in_degree.entry(link.to_segment.clone()).or_default() += 1;
        }

        // Find segments with out-degree > 1 (bubble starts)
        let bubble_starts: Vec<_> = out_degree
            .iter()
            .filter(|(_, &deg)| deg > 1)
            .map(|(name, _)| name.clone())
            .collect();

        for start in bubble_starts {
            // Find where paths from this start converge
            if let Some(bubble) = self.trace_bubble(&start, &in_degree) {
                bubbles.push(bubble);
            }
        }

        bubbles
    }

    /// Trace a bubble from a start segment to find where paths reconverge
    fn trace_bubble(&self, start: &str, in_degree: &HashMap<String, usize>) -> Option<Bubble> {
        let outgoing = self.graph.get_outgoing_links(start, Orientation::Forward);
        if outgoing.len() < 2 {
            return None;
        }

        let mut paths = Vec::new();
        let mut end_candidates: HashMap<String, usize> = HashMap::new();

        // Trace each path from the bubble start
        for link in outgoing {
            let mut path_segs = vec![(link.to_segment.clone(), link.to_orient)];
            let mut current = &link.to_segment;
            let mut current_orient = link.to_orient;

            // Follow the path until we find a convergence point
            for _ in 0..100 {
                // Safety limit
                let next_links = self.graph.get_outgoing_links(current, current_orient);
                if next_links.is_empty() {
                    break;
                }
                if next_links.len() > 1 {
                    // Another branching point - stop here
                    break;
                }
                let next = &next_links[0];
                path_segs.push((next.to_segment.clone(), next.to_orient));
                current = &next.to_segment;
                current_orient = next.to_orient;

                // Check if this segment has in-degree > 1 (convergence point)
                if in_degree.get(current).copied().unwrap_or(0) > 1 {
                    *end_candidates.entry(current.clone()).or_default() += 1;
                    break;
                }
            }

            let total_length: usize = path_segs
                .iter()
                .filter_map(|(name, _)| self.graph.get_segment(name))
                .map(|s| s.sequence_length())
                .sum();

            let is_reference = path_segs.iter().any(|(name, _)| name.contains("ref"));

            paths.push(BubblePath {
                segments: path_segs,
                total_length,
                is_reference,
            });
        }

        // Find the common end point
        let end = end_candidates
            .into_iter()
            .filter(|(_, count)| *count >= 2)
            .max_by_key(|(_, count)| *count)
            .map(|(name, _)| name)?;

        Some(Bubble {
            start_segment: start.to_string(),
            end_segment: end,
            paths,
        })
    }

    /// Analyze sample/haplotype paths
    fn analyze_samples(&self) -> Vec<SampleSummary> {
        self.graph
            .paths
            .iter()
            .map(|path| {
                let path_length: usize = path
                    .segments
                    .iter()
                    .filter_map(|s| self.graph.get_segment(&s.name))
                    .map(|seg| seg.sequence_length())
                    .sum();

                let sv_count = path
                    .segments
                    .iter()
                    .filter(|s| {
                        self.graph
                            .get_segment(&s.name)
                            .and_then(|seg| seg.get_tag("VT"))
                            .is_some()
                    })
                    .count();

                SampleSummary {
                    name: path.name.clone(),
                    haplotype: path.haplotype.unwrap_or(0),
                    path_length,
                    segment_count: path.segments.len(),
                    sv_count,
                }
            })
            .collect()
    }

    /// Compute total reference path length
    fn compute_reference_length(&self) -> usize {
        for path in &self.graph.paths {
            if path.name.contains("GRCh38") || path.name.contains("grch38") {
                return path
                    .segments
                    .iter()
                    .filter_map(|s| self.graph.get_segment(&s.name))
                    .map(|seg| seg.sequence_length())
                    .sum();
            }
        }
        0
    }

    /// Compute graph complexity metrics
    fn compute_complexity(&self) -> GraphComplexity {
        let mut degree_map: HashMap<String, usize> = HashMap::new();

        for link in &self.graph.links {
            *degree_map.entry(link.from_segment.clone()).or_default() += 1;
            *degree_map.entry(link.to_segment.clone()).or_default() += 1;
        }

        let degrees: Vec<usize> = degree_map.values().copied().collect();
        let max_degree = degrees.iter().copied().max().unwrap_or(0);
        let avg_degree = if degrees.is_empty() {
            0.0
        } else {
            degrees.iter().sum::<usize>() as f64 / degrees.len() as f64
        };

        let branching_segments = degree_map.values().filter(|&&d| d > 2).count();

        GraphComplexity {
            max_degree,
            avg_degree,
            branching_segments,
            is_dag: self.is_dag(),
        }
    }

    /// Check if the graph is a DAG (directed acyclic graph)
    fn is_dag(&self) -> bool {
        // Simple cycle detection using DFS
        let mut visited = HashSet::new();
        let mut rec_stack = HashSet::new();

        for name in self.graph.segments.keys() {
            if self.has_cycle_dfs(name, &mut visited, &mut rec_stack) {
                return false;
            }
        }
        true
    }

    fn has_cycle_dfs(
        &self,
        node: &str,
        visited: &mut HashSet<String>,
        rec_stack: &mut HashSet<String>,
    ) -> bool {
        if rec_stack.contains(node) {
            return true;
        }
        if visited.contains(node) {
            return false;
        }

        visited.insert(node.to_string());
        rec_stack.insert(node.to_string());

        for link in self.graph.get_outgoing_links(node, Orientation::Forward) {
            if self.has_cycle_dfs(&link.to_segment, visited, rec_stack) {
                return true;
            }
        }

        rec_stack.remove(node);
        false
    }

    /// Check for common issues in the graph
    fn check_for_issues(&self) -> Vec<String> {
        let mut warnings = Vec::new();

        // Check for disconnected segments
        let linked_segments: HashSet<_> = self
            .graph
            .links
            .iter()
            .flat_map(|l| vec![l.from_segment.clone(), l.to_segment.clone()])
            .collect();

        for name in self.graph.segments.keys() {
            if !linked_segments.contains(name) {
                warnings.push(format!("Segment '{}' is not connected to any other segment", name));
            }
        }

        // Check for segments in paths but not in segments list
        for path in &self.graph.paths {
            for seg in &path.segments {
                if !self.graph.segments.contains_key(&seg.name) {
                    warnings.push(format!(
                        "Path '{}' references unknown segment '{}'",
                        path.name, seg.name
                    ));
                }
            }
        }

        warnings
    }
}

impl SVEvalReport {
    /// Generate a human-readable summary
    pub fn summary(&self) -> String {
        let mut lines = Vec::new();

        lines.push("=== SV Representation Evaluation Report ===".to_string());
        lines.push(String::new());

        if let Some(ref region) = self.region {
            lines.push(format!("Region: {}", region));
        }

        lines.push(String::new());
        lines.push("--- Graph Statistics ---".to_string());
        lines.push(format!("Segments: {}", self.total_segments));
        lines.push(format!("Links: {}", self.total_links));
        lines.push(format!("Paths/Haplotypes: {}", self.total_paths));
        lines.push(format!("Reference path length: {} bp", self.reference_path_length));

        lines.push(String::new());
        lines.push("--- Graph Complexity ---".to_string());
        lines.push(format!("Max node degree: {}", self.graph_complexity.max_degree));
        lines.push(format!("Avg node degree: {:.2}", self.graph_complexity.avg_degree));
        lines.push(format!("Branching nodes: {}", self.graph_complexity.branching_segments));
        lines.push(format!("Is DAG: {}", self.graph_complexity.is_dag));

        lines.push(String::new());
        lines.push("--- Detected Structural Variants ---".to_string());
        if self.detected_svs.is_empty() {
            lines.push("No SVs detected from segment tags".to_string());
        } else {
            for (i, sv) in self.detected_svs.iter().enumerate() {
                lines.push(format!(
                    "  {}. {} - Alt segments: {:?}",
                    i + 1,
                    sv.sv_type,
                    sv.alt_segments
                ));
                lines.push(format!(
                    "     Ref length: {} bp, Alt length: {} bp, Size diff: {:+} bp",
                    sv.ref_length, sv.alt_length, sv.size_diff
                ));
                if !sv.samples_with_alt.is_empty() {
                    lines.push(format!("     Samples with ALT: {:?}", sv.samples_with_alt));
                }
            }
        }

        lines.push(String::new());
        lines.push("--- Bubble Structures ---".to_string());
        if self.bubbles.is_empty() {
            lines.push("No bubbles detected".to_string());
        } else {
            for (i, bubble) in self.bubbles.iter().enumerate() {
                lines.push(format!(
                    "  {}. {} -> {} ({} alternative paths)",
                    i + 1,
                    bubble.start_segment,
                    bubble.end_segment,
                    bubble.paths.len()
                ));
                for (j, path) in bubble.paths.iter().enumerate() {
                    let path_type = if path.is_reference { "REF" } else { "ALT" };
                    lines.push(format!(
                        "     Path {}: {} ({} bp, {} segments)",
                        j + 1,
                        path_type,
                        path.total_length,
                        path.segments.len()
                    ));
                }
            }
        }

        lines.push(String::new());
        lines.push("--- Sample/Haplotype Paths ---".to_string());
        for sample in &self.samples {
            lines.push(format!(
                "  {} (hap {}): {} bp, {} segments, {} SV segments",
                sample.name, sample.haplotype, sample.path_length, sample.segment_count, sample.sv_count
            ));
        }

        if !self.warnings.is_empty() {
            lines.push(String::new());
            lines.push("--- Warnings ---".to_string());
            for warning in &self.warnings {
                lines.push(format!("  ⚠ {}", warning));
            }
        }

        lines.join("\n")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genomic_region_parse() {
        let region = GenomicRegion::parse("chr17:10984564-10993960").unwrap();
        assert_eq!(region.chrom, "chr17");
        assert_eq!(region.start, 10984564);
        assert_eq!(region.end, 10993960);
        assert_eq!(region.length(), 9396);
    }
}

