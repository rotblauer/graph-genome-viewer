//! Graph Genome Viewer - A professional pangenome visualization tool
//!
//! This application provides intuitive visualization of graph genomes (GFA format)
//! with support for overlaying long-read sequencing alignments and evaluating
//! alignment quality metrics.

mod app;
mod graph;
mod io;
mod layout;
mod render;
mod analysis;
mod ui;

use anyhow::Result;
use eframe::egui;

fn main() -> Result<()> {
    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    log::info!("Starting Graph Genome Viewer");

    // Configure native window options
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_min_inner_size([800.0, 600.0])
            .with_title("Graph Genome Viewer"),
        ..Default::default()
    };

    // Run the application
    eframe::run_native(
        "Graph Genome Viewer",
        native_options,
        Box::new(|cc| Ok(Box::new(app::GraphGenomeApp::new(cc)))),
    )
    .map_err(|e| anyhow::anyhow!("Failed to run application: {}", e))
}

