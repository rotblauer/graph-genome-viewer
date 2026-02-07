//! Graph rendering engine

use egui::{Color32, Painter, Pos2, Rect, Stroke, Vec2};
use crate::graph::GraphGenome;
use crate::layout::{Layout, SegmentLayout, LinkLayout};
use crate::app::{Viewport, UiState, ColorScheme};
use super::colors::{ColorPalette, CoverageGradient};

pub struct GraphRenderer {
    stroke_width: f32,
    link_width: f32,
}

impl Default for GraphRenderer {
    fn default() -> Self {
        Self { stroke_width: 2.0, link_width: 2.0 }
    }
}

impl GraphRenderer {
    pub fn new() -> Self { Self::default() }

    pub fn render(
        &self, painter: &Painter, rect: Rect,
        graph: &GraphGenome, layout: &Layout,
        viewport: &Viewport, ui_state: &UiState,
    ) {
        let transform = ViewTransform::new(rect, &layout.bounds, viewport);

        for link_layout in &layout.links {
            self.render_link(painter, link_layout, &transform, ui_state);
        }

        for (name, seg_layout) in &layout.segments {
            let segment = graph.get_segment(name);
            self.render_segment(painter, seg_layout, segment, &transform, ui_state);
        }

        for (idx, path) in graph.paths.iter().enumerate() {
            self.render_path_overlay(painter, path, layout, &transform, idx);
        }
    }

    fn render_segment(
        &self, painter: &Painter, seg_layout: &SegmentLayout,
        segment: Option<&crate::graph::Segment>,
        transform: &ViewTransform, ui_state: &UiState,
    ) {
        let screen_rect = transform.segment_to_screen(seg_layout);
        if !transform.clip_rect.intersects(screen_rect) { return; }

        let (fill_color, stroke_color) = self.segment_colors(&seg_layout.name, segment, ui_state);
        painter.rect(screen_rect, 6.0, fill_color, Stroke::new(self.stroke_width, stroke_color));

        if screen_rect.width() > 30.0 {
            let label = if screen_rect.width() > 80.0 {
                if let Some(seg) = segment {
                    format!("{} ({}bp)", seg_layout.name, seg.sequence_length())
                } else { seg_layout.name.clone() }
            } else { seg_layout.name.clone() };

            painter.text(screen_rect.center(), egui::Align2::CENTER_CENTER,
                &label, egui::FontId::proportional(12.0), ColorPalette::TEXT_PRIMARY);
        }
    }

    fn segment_colors(&self, name: &str, segment: Option<&crate::graph::Segment>,
        ui_state: &UiState) -> (Color32, Color32) {
        if ui_state.selected_segment.as_ref() == Some(&name.to_string()) {
            return (ColorPalette::SEGMENT_SELECTED, ColorPalette::SEGMENT_SELECTED);
        }

        match ui_state.color_scheme {
            ColorScheme::Coverage => {
                if let Some(seg) = segment {
                    if let Some(cov) = seg.coverage() {
                        let color = CoverageGradient::color((cov / 100.0).clamp(0.0, 1.0));
                        return (color, color);
                    }
                }
            }
            ColorScheme::GcContent => {
                if let Some(seg) = segment {
                    let mut seg = seg.clone();
                    let color = CoverageGradient::color(seg.gc_content() / 100.0);
                    return (color, color);
                }
            }
            _ => {}
        }
        (ColorPalette::SEGMENT_FILL, ColorPalette::SEGMENT_STROKE)
    }

    fn render_link(&self, painter: &Painter, link: &LinkLayout,
        transform: &ViewTransform, ui_state: &UiState) {
        if link.control_points.len() < 2 { return; }

        let screen_points: Vec<Pos2> = link.control_points.iter()
            .map(|&(x, y)| transform.to_screen(x, y)).collect();

        let color = if ui_state.selected_segment.as_ref() == Some(&link.from)
            || ui_state.selected_segment.as_ref() == Some(&link.to)
        { ColorPalette::LINK_SELECTED } else { ColorPalette::LINK_DEFAULT };

        if screen_points.len() == 4 {
            self.draw_bezier(painter, &screen_points, color);
        } else {
            painter.line_segment([screen_points[0], *screen_points.last().unwrap()],
                Stroke::new(self.link_width, color));
        }
    }

    fn draw_bezier(&self, painter: &Painter, points: &[Pos2], color: Color32) {
        if points.len() != 4 { return; }
        let steps = 20;
        let mut prev = points[0];
        for i in 1..=steps {
            let t = i as f32 / steps as f32;
            let mt = 1.0 - t;
            let x = mt.powi(3)*points[0].x + 3.0*mt.powi(2)*t*points[1].x
                + 3.0*mt*t.powi(2)*points[2].x + t.powi(3)*points[3].x;
            let y = mt.powi(3)*points[0].y + 3.0*mt.powi(2)*t*points[1].y
                + 3.0*mt*t.powi(2)*points[2].y + t.powi(3)*points[3].y;
            let current = Pos2::new(x, y);
            painter.line_segment([prev, current], Stroke::new(self.link_width, color));
            prev = current;
        }
    }

    fn render_path_overlay(&self, painter: &Painter, path: &crate::graph::Path,
        layout: &Layout, transform: &ViewTransform, path_index: usize) {
        let color = ColorPalette::path_color(path_index);
        let offset = (path_index as f32 - 1.5) * 3.0;
        for window in path.segments.windows(2) {
            if let (Some(from), Some(to)) = (layout.get_segment(&window[0].name),
                layout.get_segment(&window[1].name)) {
                let start = transform.to_screen(from.right(), from.y + offset);
                let end = transform.to_screen(to.left(), to.y + offset);
                painter.line_segment([start, end], Stroke::new(2.0, color));
            }
        }
    }
}

struct ViewTransform {
    clip_rect: Rect,
    graph_center: (f32, f32),
    scale: f32,
    screen_center: Pos2,
    offset: Vec2,
}

impl ViewTransform {
    fn new(rect: Rect, bounds: &crate::layout::LayoutBounds, viewport: &Viewport) -> Self {
        let graph_width = bounds.width().max(1.0);
        let graph_height = bounds.height().max(1.0);
        let scale_x = rect.width() / graph_width;
        let scale_y = rect.height() / graph_height;
        let base_scale = scale_x.min(scale_y) * 0.9;
        Self {
            clip_rect: rect,
            graph_center: bounds.center(),
            scale: base_scale * viewport.zoom,
            screen_center: rect.center(),
            offset: viewport.center * viewport.zoom,
        }
    }

    fn to_screen(&self, x: f32, y: f32) -> Pos2 {
        let dx = (x - self.graph_center.0) * self.scale - self.offset.x;
        let dy = (y - self.graph_center.1) * self.scale - self.offset.y;
        Pos2::new(self.screen_center.x + dx, self.screen_center.y + dy)
    }

    fn segment_to_screen(&self, seg: &SegmentLayout) -> Rect {
        let center = self.to_screen(seg.x, seg.y);
        let half_width = seg.width * self.scale / 2.0;
        let half_height = seg.height * self.scale / 2.0;
        Rect::from_center_size(center, Vec2::new(half_width * 2.0, half_height * 2.0))
    }
}

