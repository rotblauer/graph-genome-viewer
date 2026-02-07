//! Color schemes and gradients for visualization

use egui::Color32;

/// Color palette for graph visualization
pub struct ColorPalette;

impl ColorPalette {
    pub const SEGMENT_FILL: Color32 = Color32::from_rgb(70, 130, 180);
    pub const SEGMENT_STROKE: Color32 = Color32::from_rgb(100, 160, 210);
    pub const SEGMENT_SELECTED: Color32 = Color32::from_rgb(255, 180, 100);
    pub const SEGMENT_HOVER: Color32 = Color32::from_rgb(90, 150, 200);

    pub const LINK_DEFAULT: Color32 = Color32::from_rgb(120, 120, 140);
    pub const LINK_HOVER: Color32 = Color32::from_rgb(180, 180, 200);
    pub const LINK_SELECTED: Color32 = Color32::from_rgb(255, 200, 100);

    pub const PATH_COLORS: [Color32; 8] = [
        Color32::from_rgb(86, 180, 233),
        Color32::from_rgb(230, 159, 0),
        Color32::from_rgb(0, 158, 115),
        Color32::from_rgb(240, 228, 66),
        Color32::from_rgb(0, 114, 178),
        Color32::from_rgb(213, 94, 0),
        Color32::from_rgb(204, 121, 167),
        Color32::from_rgb(150, 150, 150),
    ];

    pub const TEXT_PRIMARY: Color32 = Color32::from_rgb(220, 220, 230);
    pub const TEXT_SECONDARY: Color32 = Color32::from_rgb(160, 160, 170);
    pub const TEXT_DIM: Color32 = Color32::from_rgb(100, 100, 110);

    pub const BG_DARK: Color32 = Color32::from_rgb(20, 22, 28);
    pub const BG_PANEL: Color32 = Color32::from_rgb(28, 30, 36);

    pub fn path_color(index: usize) -> Color32 {
        Self::PATH_COLORS[index % Self::PATH_COLORS.len()]
    }
}

/// Gradient for coverage heatmap
pub struct CoverageGradient;

impl CoverageGradient {
    pub fn color(value: f32) -> Color32 {
        let value = value.clamp(0.0, 1.0);
        let (r, g, b) = if value < 0.25 {
            let t = value * 4.0;
            (0, (t * 200.0) as u8, 200)
        } else if value < 0.5 {
            let t = (value - 0.25) * 4.0;
            (0, 200, (200.0 - t * 200.0) as u8)
        } else if value < 0.75 {
            let t = (value - 0.5) * 4.0;
            ((t * 255.0) as u8, 200, 0)
        } else {
            let t = (value - 0.75) * 4.0;
            (255, (200.0 - t * 200.0) as u8, 0)
        };
        Color32::from_rgb(r, g, b)
    }

    pub fn identity_color(identity: f64) -> Color32 {
        let identity = identity.clamp(0.0, 1.0) as f32;
        if identity < 0.9 {
            let t = identity / 0.9;
            Color32::from_rgb(200, (t * 200.0) as u8, 0)
        } else {
            let t = (identity - 0.9) / 0.1;
            Color32::from_rgb((200.0 - t * 200.0) as u8, 200, 0)
        }
    }
}

