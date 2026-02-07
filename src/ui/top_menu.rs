//! Top menu bar

use egui::Context;
use crate::app::GraphGenomeApp;

pub struct TopMenu;

impl TopMenu {
    pub fn show(ctx: &Context, app: &mut GraphGenomeApp) {
        egui::TopBottomPanel::top("top_menu").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open GFA...").clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("GFA files", &["gfa", "gfa2"])
                            .pick_file()
                        {
                            app.load_gfa(path);
                        }
                        ui.close_menu();
                    }

                    if ui.button("Open Alignments (GAF)...").clicked() {
                        if let Some(path) = rfd::FileDialog::new()
                            .add_filter("GAF files", &["gaf"])
                            .pick_file()
                        {
                            log::info!("Would load alignments from {:?}", path);
                        }
                        ui.close_menu();
                    }

                    ui.separator();

                    if ui.button("Export PNG...").clicked() {
                        ui.close_menu();
                    }

                    ui.separator();

                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });

                ui.menu_button("View", |ui| {
                    ui.checkbox(&mut app.ui_state.show_side_panel, "Side Panel");
                    ui.checkbox(&mut app.ui_state.show_metrics_panel, "Metrics Panel");

                    ui.separator();

                    if ui.button("Reset View").clicked() {
                        app.viewport = crate::app::Viewport::default();
                        ui.close_menu();
                    }

                    if ui.button("Fit to Window").clicked() {
                        app.viewport.zoom = 1.0;
                        app.viewport.center = egui::Vec2::ZERO;
                        ui.close_menu();
                    }
                });

                ui.menu_button("Layout", |ui| {
                    if ui.button("Tube Map").clicked() {
                        app.layout_algorithm = crate::app::LayoutAlgorithmChoice::TubeMap;
                        app.recompute_layout();
                        ui.close_menu();
                    }
                    if ui.button("Force-Directed").clicked() {
                        app.layout_algorithm = crate::app::LayoutAlgorithmChoice::ForceDirected;
                        app.recompute_layout();
                        ui.close_menu();
                    }
                });

                ui.menu_button("Help", |ui| {
                    if ui.button("About").clicked() {
                        ui.close_menu();
                    }
                    if ui.button("Keyboard Shortcuts").clicked() {
                        ui.close_menu();
                    }
                });
            });
        });
    }
}

