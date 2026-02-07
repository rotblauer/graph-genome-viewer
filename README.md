# Graph Genome Viewer

A professional pangenome visualization tool built in Rust with egui for GPU-accelerated rendering.

## Features

- **Graph Visualization**: Stunning, professional visualizations of graph genomes (GFA/GFA2 format)
- **Multiple Layout Algorithms**:
  - Tube Map (linear layout inspired by transit maps)
  - Force-Directed (Fruchterman-Reingold algorithm)
- **Long-Read Alignment Analysis**: Load and visualize alignment data (GAF format)
- **Interactive Navigation**: Pan, zoom, and select segments
- **Color Schemes**: Coverage heatmaps, GC content, alignment identity
- **Cross-Platform**: Compiles to a static binary for Linux, macOS, and Windows

## Building

### Prerequisites

- Rust 1.70+ (install via [rustup](https://rustup.rs))

### Build Release Binary

```bash
cargo build --release
```

The binary will be at `target/release/graph-genome-viewer`.

### Build Static Binary (Linux)

```bash
rustup target add x86_64-unknown-linux-musl
cargo build --release --target x86_64-unknown-linux-musl
```

## Usage

### Running the Application

```bash
cargo run --release
```

Or run the binary directly:

```bash
./target/release/graph-genome-viewer
```

### Loading Files

- **File → Open GFA**: Load a graph genome in GFA or GFA2 format
- **File → Open Alignments (GAF)**: Load alignment data
- **Drag & Drop**: Drop GFA files directly onto the window

### Navigation

- **Pan**: Click and drag
- **Zoom**: Mouse wheel
- **Select**: Click on segments

### Keyboard Shortcuts

- `Ctrl+O`: Open file
- `Ctrl+Q`: Quit
- `R`: Reset view
- `F`: Fit to window

## Supported Formats

### Input
- GFA 1.0/1.1/2.0 (graph genomes)
- GAF (graph alignments)

### Export
- PNG (planned)
- SVG (planned)

## Project Structure

```
src/
├── main.rs           # Entry point
├── app.rs            # Application state and eframe App impl
├── graph/            # Core graph data structures
│   ├── mod.rs        # GraphGenome, Orientation
│   ├── segment.rs    # Segment (node) type
│   ├── link.rs       # Link (edge) type
│   └── path.rs       # Path type
├── io/               # File parsing
│   ├── gfa.rs        # GFA parser
│   ├── gaf.rs        # GAF parser
│   └── alignment.rs  # Alignment data structures
├── layout/           # Graph layout algorithms
│   ├── tubemap.rs    # Linear tube map layout
│   └── force.rs      # Force-directed layout
├── render/           # Visualization rendering
│   ├── graph_renderer.rs
│   └── colors.rs     # Color palettes
├── analysis/         # Alignment statistics
│   ├── stats.rs      # Alignment statistics
│   └── fit.rs        # Goodness-of-fit metrics
└── ui/               # UI components
    ├── side_panel.rs
    ├── top_menu.rs
    └── status_bar.rs
```

## License

MIT
