//! Analysis module for alignment statistics and fit metrics

mod stats;
mod fit;
pub mod sv_eval;

pub use stats::AlignmentStats;
pub use fit::FitMetrics;
pub use sv_eval::{SVEvaluator, SVEvalReport, GenomicRegion, DetectedSV, SVType, Bubble};

