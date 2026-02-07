# Graph Genome Viewer

> ⚠️ **AI-GENERATED PROTOTYPE** ⚠️
> 
> This code was generated with AI assistance (GitHub Copilot/Claude) as an **experimental prototype**.
> 
> **NOT FOR PRODUCTION USE** - This is a proof-of-concept for exploring pangenome visualization approaches.
> 
> - Code may contain bugs, inefficiencies, or incorrect implementations
> - Not validated for scientific accuracy
> - No warranty or support provided
> - Use at your own risk for research/educational purposes only

---

A pangenome visualization tool built in Rust with egui for GPU-accelerated rendering.

## Features

- **Graph Visualization**: Visualizations of graph genomes (GFA/GFA2 format)
- **Structural Variant View**: Clear visualization of SVs (deletions, insertions, inversions, complex events)
- **Multiple Layout Algorithms**:
  - Tube Map (linear layout inspired by transit maps)
  - Force-Directed (Fruchterman-Reingold algorithm)
- **Long-Read Alignment Analysis**: Load and visualize alignment data (GAF format)
- **Interactive Navigation**: Pan, zoom, and select segments
- **Color Schemes**: Coverage heatmaps, GC content, alignment identity
- **Cross-Platform**: Compiles to a static binary for Linux, macOS, and Windows

## Quick Start: NA19240 SV Visualization

To visualize the structural variant representation at chr17:10984564-10993960:

```bash
# Run the GUI and load the example file
cargo run --release

# Then drag-and-drop data/na19240_chr17_sv.gfa or use File > Open
```

Or run the command-line evaluation:

```bash
cargo run --example eval_na19240_sv
```

This will show you the known SVs in the NA19240 sample with an ASCII diagram:

```
Legend: [REF]=Reference  [DEL]=Deletion  [INS]=Insertion  [INV]=Inversion  [CPX]=Complex

  GRCh38 hap0     │ [===]─[===]─[REF]─[REF]─[REF]─[---]─[===]─[===]  (5520 bp)
  NA19240 hap1    │ [===]─[===]─[DEL]─[---]─[===]─[===]              (3760 bp)
  NA19240 hap2    │ [===]─[===]─[INS]─[REF]─[REF]─[---]─[===]─[===]  (7263 bp)
  HG005 hap1      │ [===]─[===]─[CPX]─<[CPX]>─[REF]─[---]─[===]─[===] (4960 bp)
```

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

## Testing

### Running Tests

```bash
# Run all tests
cargo test

# Run specific test module
cargo test --test gfa_tests

# Run with verbose output
cargo test --verbose
```

### Test Data

Sample GFA files are provided in the `data/` directory:

- `simple_bubble.gfa`: A simple bubble structure for basic testing
- `na19240_chr17_sv.gfa`: Simulated chr17 structural variant region for NA19240

## Data Processing Scripts

Python scripts in `scripts/` help process real-world pangenome data.

### Setup

```bash
cd scripts
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### Download Real NA19240 Data

The `download_na19240_data.py` script downloads real pangenome data from the
Human Pangenome Reference Consortium (HPRC) and 1000 Genomes Project.

```bash
# Show information about available data sources
python download_na19240_data.py --info

# Download small test files (quick, for testing)
python download_na19240_data.py --output ../data/ --test

# Download HPRC chr17 pangenome data (~500MB-1GB)
python download_na19240_data.py --output ../data/ --chr17

# Download 1000 Genomes structural variant VCF
python download_na19240_data.py --output ../data/ --sv-vcf
```

### Process Pangenome Data

The `process_pangenome.py` script provides additional processing capabilities:

```bash
# Generate sample test data
python process_pangenome.py generate-sample --output ../data/sample.gfa

# Generate realistic chr17 structural variant data
python process_pangenome.py generate-realistic --output ../data/chr17_sv.gfa

# Display GFA statistics
python process_pangenome.py stats --input ../data/na19240_chr17_sv.gfa

# Validate a GFA file
python process_pangenome.py validate --input ../data/sample.gfa

# Extract a genomic region from a larger GFA
python process_pangenome.py extract-region --input full.gfa --output region.gfa \
    --chrom chr17 --start 10984564 --end 10993960
```

### Target Region

The scripts are configured to process the chr17:10984564-10993960 (hg38) region,
which contains a known structural variant of interest in NA19240.

## Continuous Integration

The project uses GitHub Actions for CI/CD:

- **Test Suite**: Runs on every push/PR (Linux + macOS, stable + beta Rust)
- **Release Builds**: Creates binaries for Linux and macOS (x86_64 + ARM64)
- **Data Processing**: Validates Python scripts and sample data generation

