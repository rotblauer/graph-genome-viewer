//! Color schemes and gradients for visualization

use egui::Color32;

/// Color palette for graph visualization
pub struct ColorPalette;

impl ColorPalette {
    // Segment colors - professional and subtle
    pub const SEGMENT_FILL: Color32 = Color32::from_rgb(70, 130, 180);
    pub const SEGMENT_STROKE: Color32 = Color32::from_rgb(100, 160, 210);
    pub const SEGMENT_SELECTED: Color32 = Color32::from_rgb(255, 180, 100);
    pub const SEGMENT_HOVER: Color32 = Color32::from_rgb(90, 150, 200);

    // Link colors
    pub const LINK_DEFAULT: Color32 = Color32::from_rgb(120, 120, 140);
    pub const LINK_HOVER: Color32 = Color32::from_rgb(180, 180, 200);
    pub const LINK_SELECTED: Color32 = Color32::from_rgb(255, 200, 100);

    // Path colors - distinct, colorblind-friendly palette
    pub const PATH_COLORS: [Color32; 8] = [
        Color32::from_rgb(86, 180, 233),   // Sky blue
        Color32::from_rgb(230, 159, 0),    // Orange
        Color32::from_rgb(0, 158, 115),    // Bluish green
        Color32::from_rgb(240, 228, 66),   // Yellow
        Color32::from_rgb(0, 114, 178),    // Blue
        Color32::from_rgb(213, 94, 0),     // Vermilion
        Color32::from_rgb(204, 121, 167),  // Reddish purple
        Color32::from_rgb(170, 170, 170),  // Gray
    ];

    // Background colors
    pub const BACKGROUND: Color32 = Color32::from_rgb(30, 30, 35);
    pub const GRID: Color32 = Color32::from_rgb(50, 50, 55);

    // Text colors
    pub const TEXT_PRIMARY: Color32 = Color32::from_rgb(230, 230, 240);
    pub const TEXT_SECONDARY: Color32 = Color32::from_rgb(160, 160, 170);

    /// Get a path color by index (cycles through available colors)
    pub fn path_color(index: usize) -> Color32 {
        Self::PATH_COLORS[index % Self::PATH_COLORS.len()]
    }
}

/// Gradient for coverage visualization
pub struct CoverageGradient;

impl CoverageGradient {
    /// Get color for coverage value (normalized 0-1)
    pub fn color(normalized_coverage: f32) -> Color32 {
        let t = normalized_coverage.clamp(0.0, 1.0);

        // Blue (low) -> Green (medium) -> Yellow -> Red (high)
        if t < 0.33 {
            let tt = t / 0.33;
            Color32::from_rgb(
                (30.0 + 20.0 * tt) as u8,
                (60.0 + 100.0 * tt) as u8,
                (180.0 - 80.0 * tt) as u8,
            )
        } else if t < 0.66 {
            let tt = (t - 0.33) / 0.33;
            Color32::from_rgb(
                (50.0 + 150.0 * tt) as u8,
                (160.0 + 40.0 * tt) as u8,
                (100.0 - 50.0 * tt) as u8,
            )
        } else {
            let tt = (t - 0.66) / 0.34;
            Color32::from_rgb(
                (200.0 + 55.0 * tt) as u8,
                (200.0 - 100.0 * tt) as u8,
                (50.0 - 50.0 * tt) as u8,
            )
        }
    }
}

/// Variant type colors
pub struct VariantColors;

impl VariantColors {
    pub const SNV: Color32 = Color32::from_rgb(100, 149, 237);
    pub const INSERTION: Color32 = Color32::from_rgb(50, 205, 50);
    pub const DELETION: Color32 = Color32::from_rgb(220, 20, 60);
    pub const COMPLEX: Color32 = Color32::from_rgb(255, 165, 0);
    pub const INVERSION: Color32 = Color32::from_rgb(148, 0, 211);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_color_cycles() {
        let c0 = ColorPalette::path_color(0);
        let c8 = ColorPalette::path_color(8);
        assert_eq!(c0, c8);
    }

    #[test]
    fn test_coverage_gradient() {
        let low = CoverageGradient::color(0.0);
        let high = CoverageGradient::color(1.0);
        assert_ne!(low, high);
    }
}

