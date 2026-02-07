//! Main application state and eframe App implementation

use eframe::egui;
use std::path::PathBuf;

use crate::graph::GraphGenome;
use crate::io::gfa::GfaLoader;
use crate::layout::{Layout, LayoutAlgorithm, TubeMapLayout, ForceDirectedLayout};
use crate::render::GraphRenderer;
use crate::analysis::{AlignmentStats, SVEvalReport};
use crate::ui::{SidePanel, TopMenu, StatusBar, SVPanel};

/// Main application state
pub struct GraphGenomeApp {
    /// Currently loaded graph genome
    pub graph: Option<GraphGenome>,

    /// Computed layout for visualization
    pub layout: Option<Layout>,

    /// Current layout algorithm selection
    pub layout_algorithm: LayoutAlgorithmChoice,

    /// Graph renderer
    pub renderer: GraphRenderer,

    /// Alignment statistics (if alignments loaded)
    pub alignment_stats: Option<AlignmentStats>,

    /// SV evaluation report
    pub sv_report: Option<SVEvalReport>,

    /// Current viewport transform (pan/zoom)
    pub viewport: Viewport,

    /// UI state
    pub ui_state: UiState,

    /// File loading state
    pub loading_state: LoadingState,
}

/// Viewport for pan/zoom navigation
#[derive(Clone, Debug)]
pub struct Viewport {
    /// Center position in graph coordinates
    pub center: egui::Vec2,
    /// Zoom level (1.0 = 100%)
    pub zoom: f32,
}

impl Default for Viewport {
    fn default() -> Self {
        Self {
            center: egui::Vec2::ZERO,
            zoom: 1.0,
        }
    }
}

/// Layout algorithm choices
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum LayoutAlgorithmChoice {
    #[default]
    TubeMap,
    ForceDirected,
    Hierarchical,
}

impl LayoutAlgorithmChoice {
    pub fn name(&self) -> &'static str {
        match self {
            Self::TubeMap => "Tube Map (Linear)",
            Self::ForceDirected => "Force-Directed",
            Self::Hierarchical => "Hierarchical",
        }
    }
}

/// UI panel visibility state
#[derive(Clone, Debug)]
pub struct UiState {
    /// Show left side panel
    pub show_side_panel: bool,
    /// Show bottom metrics panel
    pub show_metrics_panel: bool,
    /// Show SV visualization panel
    pub show_sv_panel: bool,
    /// Currently selected node/segment ID
    pub selected_segment: Option<String>,
    /// Currently selected path
    pub selected_path: Option<String>,
    /// Color scheme for rendering
    pub color_scheme: ColorScheme,
    /// Show alignment overlay
    pub show_alignments: bool,
    /// Show coverage heatmap
    pub show_coverage: bool,
}

impl Default for UiState {
    fn default() -> Self {
        Self {
            show_side_panel: true,
            show_metrics_panel: true,
            show_sv_panel: true,
            selected_segment: None,
            selected_path: None,
            color_scheme: ColorScheme::default(),
            show_alignments: true,
            show_coverage: false,
        }
    }
}

/// Color scheme options
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ColorScheme {
    #[default]
    Default,
    Coverage,
    Identity,
    GcContent,
    PathHighlight,
}

impl ColorScheme {
    pub fn name(&self) -> &'static str {
        match self {
            Self::Default => "Default",
            Self::Coverage => "Coverage Heatmap",
            Self::Identity => "Alignment Identity",
            Self::GcContent => "GC Content",
            Self::PathHighlight => "Path Highlight",
        }
    }
}

/// File loading state for async operations
#[derive(Clone, Debug, Default)]
pub enum LoadingState {
    #[default]
    Idle,
    Loading(String),
    Error(String),
}

impl GraphGenomeApp {
    /// Create a new application instance
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        // Configure custom fonts and styles
        configure_styles(&cc.egui_ctx);

        Self {
            graph: None,
            layout: None,
            layout_algorithm: LayoutAlgorithmChoice::default(),
            renderer: GraphRenderer::new(),
            alignment_stats: None,
            sv_report: None,
            viewport: Viewport::default(),
            ui_state: UiState::default(),
            loading_state: LoadingState::Idle,
        }
    }

    /// Load a GFA file
    pub fn load_gfa(&mut self, path: PathBuf) {
        self.loading_state = LoadingState::Loading(format!("Loading {:?}...", path.file_name()));

        match GfaLoader::load(&path) {
            Ok(graph) => {
                log::info!("Loaded graph with {} segments and {} links",
                    graph.segments.len(), graph.links.len());

                // Compute initial layout
                let layout = self.compute_layout(&graph);
                self.layout = Some(layout);
                self.graph = Some(graph);
                self.sv_report = None;  // Reset SV report for new graph
                self.viewport = Viewport::default();
                self.loading_state = LoadingState::Idle;
            }
            Err(e) => {
                log::error!("Failed to load GFA: {}", e);
                self.loading_state = LoadingState::Error(format!("Failed to load: {}", e));
            }
        }
    }

    /// Compute layout based on current algorithm choice
    fn compute_layout(&self, graph: &GraphGenome) -> Layout {
        match self.layout_algorithm {
            LayoutAlgorithmChoice::TubeMap => {
                TubeMapLayout::new().compute(graph)
            }
            LayoutAlgorithmChoice::ForceDirected => {
                ForceDirectedLayout::new().compute(graph)
            }
            LayoutAlgorithmChoice::Hierarchical => {
                // TODO: Implement hierarchical layout
                TubeMapLayout::new().compute(graph)
            }
        }
    }

    /// Recompute layout with current algorithm
    pub fn recompute_layout(&mut self) {
        if let Some(ref graph) = self.graph {
            self.layout = Some(self.compute_layout(graph));
        }
    }

    /// Handle viewport pan
    fn handle_pan(&mut self, delta: egui::Vec2) {
        self.viewport.center -= delta / self.viewport.zoom;
    }

    /// Handle viewport zoom
    fn handle_zoom(&mut self, delta: f32, mouse_pos: egui::Vec2) {
        let old_zoom = self.viewport.zoom;
        self.viewport.zoom = (self.viewport.zoom * (1.0 + delta * 0.1)).clamp(0.1, 10.0);

        // Zoom towards mouse position
        let zoom_factor = self.viewport.zoom / old_zoom;
        self.viewport.center = mouse_pos + (self.viewport.center - mouse_pos) * zoom_factor;
    }
}

impl eframe::App for GraphGenomeApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Top menu bar
        TopMenu::show(ctx, self);

        // Side panel (graph info, paths, settings)
        if self.ui_state.show_side_panel {
            SidePanel::show(ctx, self);
        }

        // SV visualization panel (bottom)
        if self.ui_state.show_sv_panel {
            SVPanel::show(ctx, self);
        }

        // Bottom status bar
        StatusBar::show(ctx, self);

        // Central graph view
        egui::CentralPanel::default().show(ctx, |ui| {
            self.show_graph_view(ui);
        });

        // Handle file drop
        self.handle_dropped_files(ctx);
    }
}

impl GraphGenomeApp {
    /// Show the main graph visualization area
    fn show_graph_view(&mut self, ui: &mut egui::Ui) {
        let available_size = ui.available_size();

        // Allocate space for the graph view
        let (response, painter) = ui.allocate_painter(available_size, egui::Sense::click_and_drag());
        let rect = response.rect;

        // Draw background
        painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(20, 22, 28));

        // Handle input
        if response.dragged() {
            self.handle_pan(response.drag_delta());
        }

        if let Some(hover_pos) = response.hover_pos() {
            let scroll_delta = ui.input(|i| i.raw_scroll_delta.y);
            if scroll_delta != 0.0 {
                let mouse_vec = hover_pos - rect.min;
                self.handle_zoom(scroll_delta * 0.01, mouse_vec);
            }
        }

        // Render graph
        if let (Some(ref graph), Some(ref layout)) = (&self.graph, &self.layout) {
            self.renderer.render(
                &painter,
                rect,
                graph,
                layout,
                &self.viewport,
                &self.ui_state,
            );
        } else {
            // Show welcome message
            self.show_welcome_message(&painter, rect);
        }
    }

    /// Show welcome message when no graph is loaded
    fn show_welcome_message(&self, painter: &egui::Painter, rect: egui::Rect) {
        let center = rect.center();

        // Title
        painter.text(
            center - egui::vec2(0.0, 60.0),
            egui::Align2::CENTER_CENTER,
            "Graph Genome Viewer",
            egui::FontId::proportional(32.0),
            egui::Color32::from_rgb(200, 200, 210),
        );

        // Prototype warning
        painter.text(
            center - egui::vec2(0.0, 30.0),
            egui::Align2::CENTER_CENTER,
            "⚠️ AI-Generated Prototype - Not for Production Use",
            egui::FontId::proportional(12.0),
            egui::Color32::from_rgb(255, 180, 100),
        );

        // Subtitle
        painter.text(
            center + egui::vec2(0.0, 10.0),
            egui::Align2::CENTER_CENTER,
            "Drop a GFA file here or use File → Open",
            egui::FontId::proportional(16.0),
            egui::Color32::from_rgb(128, 128, 140),
        );

        // Instructions
        let instructions = "Supported formats: GFA, GFA2";
        painter.text(
            center + egui::vec2(0.0, 40.0),
            egui::Align2::CENTER_CENTER,
            instructions,
            egui::FontId::proportional(12.0),
            egui::Color32::from_rgb(100, 100, 110),
        );
    }

    /// Handle files dropped onto the window
    fn handle_dropped_files(&mut self, ctx: &egui::Context) {
        ctx.input(|i| {
            for file in &i.raw.dropped_files {
                if let Some(path) = &file.path {
                    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
                    if ext == "gfa" || ext == "gfa2" {
                        self.load_gfa(path.clone());
                    }
                }
            }
        });
    }
}

/// Configure application styles and fonts
fn configure_styles(ctx: &egui::Context) {
    let mut style = (*ctx.style()).clone();

    // Use a dark theme with subtle colors
    style.visuals = egui::Visuals::dark();
    style.visuals.panel_fill = egui::Color32::from_rgb(28, 30, 36);
    style.visuals.window_fill = egui::Color32::from_rgb(32, 34, 40);
    style.visuals.extreme_bg_color = egui::Color32::from_rgb(20, 22, 28);

    // Selection colors
    style.visuals.selection.bg_fill = egui::Color32::from_rgb(60, 100, 180);
    style.visuals.selection.stroke = egui::Stroke::new(1.0, egui::Color32::from_rgb(100, 140, 220));

    // Widget styling
    style.visuals.widgets.inactive.bg_fill = egui::Color32::from_rgb(45, 48, 55);
    style.visuals.widgets.hovered.bg_fill = egui::Color32::from_rgb(55, 58, 68);
    style.visuals.widgets.active.bg_fill = egui::Color32::from_rgb(65, 68, 78);

    ctx.set_style(style);
}

