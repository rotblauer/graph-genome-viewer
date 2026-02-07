//! NA19240 SV Evaluation Example
//!
//! Evaluates the graph representation of known SVs near chr17:10984564-10993960
//!
//! Run with: cargo run --example eval_na19240_sv
//!
//! For visualization, run the main app and load data/na19240_chr17_sv.gfa:
//!   cargo run --release
//!   Then drag-and-drop the GFA file or use File > Open

use graph_genome_viewer::analysis::{SVEvaluator, GenomicRegion};
use graph_genome_viewer::io::GfaLoader;
use std::path::Path;

fn main() {
    // Load the NA19240 SV graph
    let gfa_path = Path::new("data/na19240_chr17_sv.gfa");

    println!("╔══════════════════════════════════════════════════════════════════╗");
    println!("║         NA19240 Structural Variant Graph Evaluation              ║");
    println!("║              chr17:10984564-10993960                             ║");
    println!("╚══════════════════════════════════════════════════════════════════╝\n");

    println!("Loading GFA file: {}\n", gfa_path.display());

    let graph = match GfaLoader::load(gfa_path) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("Failed to load GFA file: {:?}", e);
            eprintln!("\nMake sure you're running from the project root directory.");
            return;
        }
    };

    // Define the region of interest
    let region = GenomicRegion::new("chr17", 10984564, 10993960);

    // Run the SV evaluation
    let evaluator = SVEvaluator::new(&graph).with_region(region);
    let report = evaluator.evaluate();

    // Print the summary report
    println!("{}", report.summary());

    // Visual representation (ASCII)
    println!("\n\n╔══════════════════════════════════════════════════════════════════╗");
    println!("║                    Graph Structure (ASCII)                       ║");
    println!("╚══════════════════════════════════════════════════════════════════╝\n");

    print_ascii_graph(&graph);

    println!("\n┌──────────────────────────────────────────────────────────────────┐");
    println!("│  For interactive visualization, run the GUI:                     │");
    println!("│    cargo run --release                                           │");
    println!("│  Then load: data/na19240_chr17_sv.gfa                            │");
    println!("└──────────────────────────────────────────────────────────────────┘");
}

fn print_ascii_graph(graph: &graph_genome_viewer::GraphGenome) {
    println!("Legend: [REF]=Reference  [DEL]=Deletion  [INS]=Insertion  [INV]=Inversion  [CPX]=Complex  [SNV]=SNV\n");

    // Print each path
    for path in &graph.paths {
        let name = format_path_name(&path.name, path.haplotype);
        let is_ref = path.name.contains("GRCh38");

        if is_ref {
            print!("  {:15} │ ", name);
        } else {
            print!("  {:15} │ ", name);
        }

        for (i, seg) in path.segments.iter().enumerate() {
            if i > 0 {
                print!("─");
            }

            let seg_info = graph.get_segment(&seg.name);
            let vt = seg_info.and_then(|s| s.get_tag("VT"));

            let (prefix, suffix) = match seg.orientation {
                graph_genome_viewer::Orientation::Forward => ("", ""),
                graph_genome_viewer::Orientation::Reverse => ("<", ">"),
            };

            let display = if let Some(vt) = vt {
                let vt_str = format!("{:?}", vt);
                if vt_str.contains("DEL") {
                    format!("{}[DEL]{}", prefix, suffix)
                } else if vt_str.contains("INS") {
                    format!("{}[INS]{}", prefix, suffix)
                } else if vt_str.contains("COMPLEX") {
                    format!("{}[CPX]{}", prefix, suffix)
                } else if vt_str.contains("SNV") {
                    format!("{}[SNV]{}", prefix, suffix)
                } else {
                    format!("{}[???]{}", prefix, suffix)
                }
            } else if seg.name.contains("flank") {
                format!("{}[===]{}", prefix, suffix)
            } else if seg.name.contains("ref") {
                format!("{}[REF]{}", prefix, suffix)
            } else {
                format!("{}[---]{}", prefix, suffix)
            };

            print!("{}", display);
        }

        // Calculate total length
        let total_len: usize = path.segments.iter()
            .filter_map(|s| graph.get_segment(&s.name))
            .map(|seg| seg.sequence_length())
            .sum();

        println!("  ({} bp)", total_len);
    }

    println!("\n  Flanking regions shown as [===], Reference segments as [REF]");
}

fn format_path_name(name: &str, haplotype: Option<u32>) -> String {
    let parts: Vec<&str> = name.split('#').collect();
    let short_name = parts.first().unwrap_or(&name);

    if let Some(hap) = haplotype {
        format!("{} hap{}", short_name, hap)
    } else {
        short_name.to_string()
    }
}

