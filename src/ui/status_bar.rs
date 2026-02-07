//! Bottom status bar

use egui::Context;
use crate::app::{GraphGenomeApp, LoadingState};

pub struct StatusBar;

impl StatusBar {
    pub fn show(ctx: &Context, app: &GraphGenomeApp) {
        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                match &app.loading_state {
                    LoadingState::Idle => {
                        if app.graph.is_some() {
                            ui.label("Ready");
                        } else {
                            ui.label("No graph loaded");
                        }
                    }
                    LoadingState::Loading(msg) => {
                        ui.spinner();
                        ui.label(msg);
                    }
                    LoadingState::Error(msg) => {
                        ui.colored_label(egui::Color32::RED, format!("Error: {}", msg));
                    }
                }

                ui.separator();
                ui.label(format!("Zoom: {:.0}%", app.viewport.zoom * 100.0));

                if let Some(ref graph) = app.graph {
                    ui.separator();
                    ui.label(format!("{} segments", graph.node_count()));
                    ui.label(format!("{} links", graph.edge_count()));
                }

                if let Some(ref seg) = app.ui_state.selected_segment {
                    ui.separator();
                    ui.label(format!("Selected: {}", seg));
                }
            });
        });
    }
}

