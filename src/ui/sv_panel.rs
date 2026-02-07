//! SV (Structural Variant) visualization panel
//!
//! Provides clear visualization of structural variant events in the graph

use egui::{Context, Ui, Color32, Stroke, Rect, Pos2, Vec2, Painter};
use crate::app::GraphGenomeApp;
use crate::analysis::sv_eval::{SVEvaluator, SVEvalReport, SVType, DetectedSV};
use crate::graph::{GraphGenome, Orientation};
use crate::render::colors::VariantColors;
use std::collections::HashMap;

/// SV Visualization Panel
pub struct SVPanel;

impl SVPanel {
    /// Show the SV panel as a bottom panel
    pub fn show(ctx: &Context, app: &mut GraphGenomeApp) {
        egui::TopBottomPanel::bottom("sv_panel")
            .default_height(250.0)
            .resizable(true)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.heading("🧬 Structural Variant View");
                    ui.separator();
                    if ui.button("Refresh Analysis").clicked() {
                        if let Some(ref graph) = app.graph {
                            app.sv_report = Some(SVEvaluator::new(graph).evaluate());
                        }
                    }
                });
                ui.separator();

                if let Some(ref graph) = app.graph {
                    // Run analysis if not already done
                    if app.sv_report.is_none() {
                        app.sv_report = Some(SVEvaluator::new(graph).evaluate());
                    }

                    egui::ScrollArea::horizontal().show(ui, |ui| {
                        Self::render_sv_diagram(ui, graph, app.sv_report.as_ref());
                    });
                } else {
                    ui.centered_and_justified(|ui| {
                        ui.label("Load a GFA file to view structural variants");
                    });
                }
            });
    }

    /// Render the main SV diagram showing all paths through the graph
    fn render_sv_diagram(ui: &mut Ui, graph: &GraphGenome, report: Option<&SVEvalReport>) {
        let available_width = ui.available_width().max(1200.0);
        let height = 180.0;

        let (response, painter) = ui.allocate_painter(
            Vec2::new(available_width, height),
            egui::Sense::hover()
        );
        let rect = response.rect;

        // Dark background
        painter.rect_filled(rect, 4.0, Color32::from_rgb(25, 27, 32));

        // Get segment ordering from reference path or topological sort
        let segment_order = Self::get_segment_order(graph);
        if segment_order.is_empty() {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "No paths found in graph",
                egui::FontId::proportional(14.0),
                Color32::GRAY
            );
            return;
        }

        // Calculate segment positions
        let margin = 60.0;
        let usable_width = available_width - margin * 2.0;

        // Get segment lengths for proportional spacing
        let total_len: usize = segment_order.iter()
            .filter_map(|name| graph.get_segment(name))
            .map(|s| s.sequence_length())
            .sum();

        let mut segment_positions: HashMap<String, (f32, f32)> = HashMap::new();
        let mut x = margin;

        for name in &segment_order {
            if let Some(seg) = graph.get_segment(name) {
                let width = (seg.sequence_length() as f32 / total_len as f32 * usable_width)
                    .max(40.0)
                    .min(150.0);
                segment_positions.insert(name.clone(), (x, width));
                x += width + 8.0;
            }
        }

        // Draw each path as a horizontal track
        let path_spacing = 28.0;
        let start_y = rect.min.y + 35.0;

        // Draw header with segment names
        Self::draw_segment_headers(&painter, rect, &segment_order, &segment_positions, graph);

        // Draw each sample/path
        for (idx, path) in graph.paths.iter().enumerate() {
            let y = start_y + idx as f32 * path_spacing;
            if y > rect.max.y - 20.0 {
                break;
            }
            Self::draw_path_track(&painter, rect, path, &segment_positions, graph, y, idx);
        }

        // Draw legend
        Self::draw_legend(&painter, rect);
    }

    /// Get segment order from reference path or by topological order
    fn get_segment_order(graph: &GraphGenome) -> Vec<String> {
        // Try to find reference path first
        for path in &graph.paths {
            if path.name.contains("GRCh38") || path.name.contains("grch38") || path.name.contains("ref") {
                return path.segments.iter().map(|s| s.name.clone()).collect();
            }
        }
        // Fall back to first path
        if let Some(path) = graph.paths.first() {
            return path.segments.iter().map(|s| s.name.clone()).collect();
        }
        // Fall back to all segments
        graph.segments.keys().cloned().collect()
    }

    /// Draw segment header labels
    fn draw_segment_headers(
        painter: &Painter,
        rect: Rect,
        segment_order: &[String],
        positions: &HashMap<String, (f32, f32)>,
        graph: &GraphGenome
    ) {
        let y = rect.min.y + 12.0;

        for name in segment_order {
            if let Some(&(x, width)) = positions.get(name) {
                let center_x = x + width / 2.0;

                // Determine segment type for coloring
                let color = if let Some(seg) = graph.get_segment(name) {
                    Self::get_segment_color(seg)
                } else {
                    Color32::from_rgb(150, 150, 160)
                };

                // Shortened name for display
                let display_name = if name.len() > 12 {
                    format!("{}...", &name[..10])
                } else {
                    name.clone()
                };

                painter.text(
                    Pos2::new(center_x, y),
                    egui::Align2::CENTER_CENTER,
                    &display_name,
                    egui::FontId::proportional(9.0),
                    color
                );
            }
        }
    }

    /// Get color for a segment based on its variant type
    fn get_segment_color(segment: &crate::graph::Segment) -> Color32 {
        if let Some(vt) = segment.get_tag("VT") {
            match format!("{:?}", vt).as_str() {
                s if s.contains("DEL") => VariantColors::DELETION,
                s if s.contains("INS") => VariantColors::INSERTION,
                s if s.contains("INV") => VariantColors::INVERSION,
                s if s.contains("COMPLEX") => VariantColors::COMPLEX,
                s if s.contains("SNV") => VariantColors::SNV,
                _ => Color32::from_rgb(150, 150, 160),
            }
        } else if segment.name.contains("ref") {
            Color32::from_rgb(100, 149, 237)  // Reference blue
        } else if segment.name.contains("flank") {
            Color32::from_rgb(128, 128, 140)  // Flanking gray
        } else {
            Color32::from_rgb(150, 150, 160)
        }
    }

    /// Draw a single path track
    fn draw_path_track(
        painter: &Painter,
        rect: Rect,
        path: &crate::graph::Path,
        positions: &HashMap<String, (f32, f32)>,
        graph: &GraphGenome,
        y: f32,
        path_idx: usize
    ) {
        let track_height = 18.0;
        let is_reference = path.name.contains("GRCh38") || path.name.contains("grch38");

        // Draw path label
        let label = if let Some(hap) = path.haplotype {
            format!("{} h{}", Self::short_name(&path.name), hap)
        } else {
            Self::short_name(&path.name)
        };

        let label_color = if is_reference {
            Color32::from_rgb(100, 200, 255)
        } else {
            crate::render::ColorPalette::path_color(path_idx)
        };

        painter.text(
            Pos2::new(rect.min.x + 5.0, y + track_height / 2.0),
            egui::Align2::LEFT_CENTER,
            &label,
            egui::FontId::proportional(10.0),
            label_color
        );

        // Draw segments in this path
        let mut prev_end: Option<f32> = None;

        for seg_ref in &path.segments {
            if let Some(&(seg_x, seg_width)) = positions.get(&seg_ref.name) {
                let segment = graph.get_segment(&seg_ref.name);
                let color = segment.map(|s| Self::get_segment_color(s))
                    .unwrap_or(Color32::GRAY);

                // Draw connector from previous segment
                if let Some(prev_x) = prev_end {
                    let connector_color = Color32::from_rgb(80, 80, 90);
                    painter.line_segment(
                        [Pos2::new(prev_x, y + track_height / 2.0),
                         Pos2::new(seg_x, y + track_height / 2.0)],
                        Stroke::new(1.5, connector_color)
                    );
                }

                // Draw segment box
                let seg_rect = Rect::from_min_size(
                    Pos2::new(seg_x, y),
                    Vec2::new(seg_width, track_height)
                );

                // Fill with segment color
                let fill_color = Color32::from_rgba_unmultiplied(
                    color.r(), color.g(), color.b(), 180
                );
                painter.rect_filled(seg_rect, 3.0, fill_color);

                // Border
                painter.rect_stroke(seg_rect, 3.0, Stroke::new(1.0, color));

                // Inversion indicator (arrow pointing left)
                if seg_ref.orientation == Orientation::Reverse {
                    let arrow_y = y + track_height / 2.0;
                    let arrow_x = seg_x + seg_width / 2.0;

                    // Draw inversion arrow
                    painter.line_segment(
                        [Pos2::new(arrow_x + 6.0, arrow_y),
                         Pos2::new(arrow_x - 6.0, arrow_y)],
                        Stroke::new(2.0, Color32::WHITE)
                    );
                    painter.line_segment(
                        [Pos2::new(arrow_x - 6.0, arrow_y),
                         Pos2::new(arrow_x - 2.0, arrow_y - 4.0)],
                        Stroke::new(2.0, Color32::WHITE)
                    );
                    painter.line_segment(
                        [Pos2::new(arrow_x - 6.0, arrow_y),
                         Pos2::new(arrow_x - 2.0, arrow_y + 4.0)],
                        Stroke::new(2.0, Color32::WHITE)
                    );
                }

                // Show segment size if room
                if seg_width > 35.0 {
                    if let Some(seg) = segment {
                        let size_text = Self::format_bp(seg.sequence_length());
                        painter.text(
                            seg_rect.center(),
                            egui::Align2::CENTER_CENTER,
                            &size_text,
                            egui::FontId::proportional(8.0),
                            Color32::WHITE
                        );
                    }
                }

                prev_end = Some(seg_x + seg_width);
            }
        }
    }

    /// Draw the legend explaining colors
    fn draw_legend(painter: &Painter, rect: Rect) {
        let legend_items = [
            ("REF", Color32::from_rgb(100, 149, 237)),
            ("DEL", VariantColors::DELETION),
            ("INS", VariantColors::INSERTION),
            ("INV", VariantColors::INVERSION),
            ("COMPLEX", VariantColors::COMPLEX),
            ("SNV", VariantColors::SNV),
        ];

        let mut x = rect.max.x - 350.0;
        let y = rect.max.y - 18.0;

        for (label, color) in legend_items {
            // Color box
            let box_rect = Rect::from_min_size(
                Pos2::new(x, y),
                Vec2::new(12.0, 12.0)
            );
            painter.rect_filled(box_rect, 2.0, color);

            // Label
            painter.text(
                Pos2::new(x + 16.0, y + 6.0),
                egui::Align2::LEFT_CENTER,
                label,
                egui::FontId::proportional(9.0),
                Color32::from_rgb(180, 180, 190)
            );

            x += 55.0;
        }
    }

    /// Shorten sample names for display
    fn short_name(name: &str) -> String {
        // Handle W-line format like "NA19240#1#chr17"
        let parts: Vec<&str> = name.split('#').collect();
        if parts.len() >= 1 {
            parts[0].to_string()
        } else {
            name.to_string()
        }
    }

    /// Format base pairs for display
    fn format_bp(bp: usize) -> String {
        if bp >= 1000 {
            format!("{:.1}kb", bp as f32 / 1000.0)
        } else {
            format!("{}bp", bp)
        }
    }
}

/// SV Summary panel showing statistics
pub struct SVSummaryPanel;

impl SVSummaryPanel {
    pub fn show(ui: &mut Ui, report: &SVEvalReport) {
        ui.collapsing("SV Summary", |ui| {
            ui.horizontal(|ui| {
                ui.label("Total SVs detected:");
                ui.strong(format!("{}", report.detected_svs.len()));
            });

            if !report.detected_svs.is_empty() {
                ui.separator();

                for sv in &report.detected_svs {
                    ui.horizontal(|ui| {
                        let (icon, color) = match sv.sv_type {
                            SVType::Deletion => ("⊖", VariantColors::DELETION),
                            SVType::Insertion => ("⊕", VariantColors::INSERTION),
                            SVType::Inversion => ("↺", VariantColors::INVERSION),
                            SVType::Complex => ("◈", VariantColors::COMPLEX),
                            SVType::SNV => ("•", VariantColors::SNV),
                            SVType::Unknown => ("?", Color32::GRAY),
                        };

                        ui.colored_label(color, icon);
                        ui.label(format!("{}", sv.sv_type));
                        ui.label("-");
                        ui.label(format!("{:?}", sv.alt_segments));

                        if sv.size_diff != 0 {
                            let size_str = if sv.size_diff > 0 {
                                format!("+{} bp", sv.size_diff)
                            } else {
                                format!("{} bp", sv.size_diff)
                            };
                            ui.small(size_str);
                        }
                    });
                }
            }

            if !report.bubbles.is_empty() {
                ui.separator();
                ui.label(format!("Bubble structures: {}", report.bubbles.len()));
            }
        });
    }
}

