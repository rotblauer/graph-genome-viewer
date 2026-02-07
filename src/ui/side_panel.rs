//! Side panel with graph information and controls

use egui::{Context, Ui};
use crate::app::{GraphGenomeApp, LayoutAlgorithmChoice, ColorScheme};

pub struct SidePanel;

impl SidePanel {
    pub fn show(ctx: &Context, app: &mut GraphGenomeApp) {
        egui::SidePanel::left("side_panel")
            .default_width(280.0)
            .resizable(true)
            .show(ctx, |ui| {
                ui.heading("Graph Genome Viewer");
                ui.separator();

                Self::graph_info_section(ui, app);
                ui.separator();

                Self::layout_section(ui, app);
                ui.separator();

                Self::display_section(ui, app);
                ui.separator();

                Self::paths_section(ui, app);
            });
    }

    fn graph_info_section(ui: &mut Ui, app: &GraphGenomeApp) {
        ui.collapsing("Graph Info", |ui| {
            if let Some(ref graph) = app.graph {
                ui.horizontal(|ui| {
                    ui.label("Segments:");
                    ui.label(format!("{}", graph.node_count()));
                });
                ui.horizontal(|ui| {
                    ui.label("Links:");
                    ui.label(format!("{}", graph.edge_count()));
                });
                ui.horizontal(|ui| {
                    ui.label("Paths:");
                    ui.label(format!("{}", graph.paths.len()));
                });
                ui.horizontal(|ui| {
                    ui.label("Total bp:");
                    ui.label(format!("{}", graph.total_sequence_length()));
                });
                if graph.is_dag() {
                    ui.label("✓ Acyclic (DAG)");
                } else {
                    ui.label("○ Contains cycles");
                }
            } else {
                ui.label("No graph loaded");
            }
        });
    }

    fn layout_section(ui: &mut Ui, app: &mut GraphGenomeApp) {
        ui.collapsing("Layout", |ui| {
            ui.horizontal(|ui| {
                ui.label("Algorithm:");
                egui::ComboBox::from_id_salt("layout_algo")
                    .selected_text(app.layout_algorithm.name())
                    .show_ui(ui, |ui| {
                        let changed = ui.selectable_value(
                            &mut app.layout_algorithm,
                            LayoutAlgorithmChoice::TubeMap,
                            "Tube Map (Linear)"
                        ).changed();
                        let changed = changed || ui.selectable_value(
                            &mut app.layout_algorithm,
                            LayoutAlgorithmChoice::ForceDirected,
                            "Force-Directed"
                        ).changed();
                        if changed {
                            app.recompute_layout();
                        }
                    });
            });

            if ui.button("Recompute Layout").clicked() {
                app.recompute_layout();
            }
        });
    }

    fn display_section(ui: &mut Ui, app: &mut GraphGenomeApp) {
        ui.collapsing("Display Options", |ui| {
            ui.horizontal(|ui| {
                ui.label("Color by:");
                egui::ComboBox::from_id_salt("color_scheme")
                    .selected_text(app.ui_state.color_scheme.name())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut app.ui_state.color_scheme, ColorScheme::Default, "Default");
                        ui.selectable_value(&mut app.ui_state.color_scheme, ColorScheme::Coverage, "Coverage");
                        ui.selectable_value(&mut app.ui_state.color_scheme, ColorScheme::GcContent, "GC Content");
                        ui.selectable_value(&mut app.ui_state.color_scheme, ColorScheme::Identity, "Identity");
                    });
            });

            ui.checkbox(&mut app.ui_state.show_alignments, "Show Alignments");
            ui.checkbox(&mut app.ui_state.show_coverage, "Show Coverage Heatmap");
        });
    }

    fn paths_section(ui: &mut Ui, app: &mut GraphGenomeApp) {
        ui.collapsing("Paths", |ui| {
            if let Some(ref graph) = app.graph {
                if graph.paths.is_empty() {
                    ui.label("No paths in graph");
                } else {
                    egui::ScrollArea::vertical().max_height(200.0).show(ui, |ui| {
                        for (idx, path) in graph.paths.iter().enumerate() {
                            let color = crate::render::ColorPalette::path_color(idx);
                            ui.horizontal(|ui| {
                                let rect = ui.available_rect_before_wrap();
                                let size = egui::vec2(12.0, 12.0);
                                let pos = egui::pos2(rect.min.x, rect.min.y + 2.0);
                                ui.painter().rect_filled(
                                    egui::Rect::from_min_size(pos, size),
                                    2.0, color
                                );
                                ui.add_space(16.0);

                                let selected = app.ui_state.selected_path.as_ref() == Some(&path.name);
                                if ui.selectable_label(selected, &path.name).clicked() {
                                    app.ui_state.selected_path = if selected {
                                        None
                                    } else {
                                        Some(path.name.clone())
                                    };
                                }
                            });
                        }
                    });
                }
            } else {
                ui.label("Load a graph first");
            }
        });
    }
}

