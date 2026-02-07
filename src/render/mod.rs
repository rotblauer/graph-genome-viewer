//! Rendering module for graph visualization

mod graph_renderer;
pub mod colors;

pub use graph_renderer::GraphRenderer;
pub use colors::{ColorPalette, CoverageGradient, VariantColors};

