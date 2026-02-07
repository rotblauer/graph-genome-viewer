//! Alignment statistics calculations

use crate::io::AlignmentData;
use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
pub struct AlignmentStats {
    pub total_reads: usize,
    pub total_aligned_bases: usize,
    pub mean_identity: f64,
    pub median_identity: f64,
    pub mean_mapq: f32,
    pub coverage_by_segment: HashMap<String, f32>,
    pub identity_distribution: Vec<(f64, usize)>,
    pub mapq_distribution: Vec<(u8, usize)>,
}

impl AlignmentStats {
    pub fn compute(data: &AlignmentData) -> Self {
        if data.alignments.is_empty() {
            return Self::default();
        }

        let total_reads = data.alignments.len();
        let total_aligned_bases: usize = data.alignments.iter()
            .map(|a| a.read_end - a.read_start).sum();

        let mut identities: Vec<f64> = data.alignments.iter()
            .map(|a| a.identity).collect();
        identities.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mean_identity = identities.iter().sum::<f64>() / total_reads as f64;
        let median_identity = identities[total_reads / 2];

        let mean_mapq = data.alignments.iter()
            .map(|a| a.mapq as f32).sum::<f32>() / total_reads as f32;

        let mut coverage_by_segment = HashMap::new();
        for (name, cov) in &data.segment_coverage {
            coverage_by_segment.insert(name.clone(), cov.read_count as f32);
        }

        let mut identity_bins = vec![0usize; 20];
        for &id in &identities {
            let bin = ((id * 20.0) as usize).min(19);
            identity_bins[bin] += 1;
        }
        let identity_distribution: Vec<(f64, usize)> = identity_bins.iter()
            .enumerate()
            .map(|(i, &count)| (i as f64 * 0.05, count))
            .collect();

        let mut mapq_bins = HashMap::new();
        for a in &data.alignments {
            *mapq_bins.entry(a.mapq).or_insert(0usize) += 1;
        }
        let mut mapq_distribution: Vec<(u8, usize)> = mapq_bins.into_iter().collect();
        mapq_distribution.sort_by_key(|&(q, _)| q);

        Self {
            total_reads,
            total_aligned_bases,
            mean_identity,
            median_identity,
            mean_mapq,
            coverage_by_segment,
            identity_distribution,
            mapq_distribution,
        }
    }
}

