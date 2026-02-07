//! Force-directed layout using Fruchterman-Reingold algorithm

use crate::graph::GraphGenome;
use super::{Layout, LayoutAlgorithm, SegmentLayout, LinkLayout};
use rand::Rng;
use std::collections::HashMap;

pub struct ForceDirectedLayout {
    iterations: usize,
    ideal_length: f32,
    cooling: f32,
    temperature: f32,
    segment_height: f32,
}

impl Default for ForceDirectedLayout {
    fn default() -> Self {
        Self {
            iterations: 100,
            ideal_length: 150.0,
            cooling: 0.95,
            temperature: 200.0,
            segment_height: 30.0,
        }
    }
}

impl ForceDirectedLayout {
    pub fn new() -> Self { Self::default() }
}

#[derive(Clone)]
struct NodeState { x: f32, y: f32, dx: f32, dy: f32 }

impl LayoutAlgorithm for ForceDirectedLayout {
    fn compute(&self, graph: &GraphGenome) -> Layout {
        let mut layout = Layout::new();
        if graph.segments.is_empty() { return layout; }

        let mut rng = rand::thread_rng();
        let segment_names: Vec<String> = graph.segments.keys().cloned().collect();
        let n = segment_names.len();

        let mut states: HashMap<String, NodeState> = segment_names.iter()
            .map(|name| (name.clone(), NodeState {
                x: rng.gen_range(-200.0..200.0),
                y: rng.gen_range(-200.0..200.0),
                dx: 0.0, dy: 0.0,
            }))
            .collect();

        let area = (n as f32 * self.ideal_length * self.ideal_length).sqrt();
        let k = (area / n as f32).sqrt();
        let mut temp = self.temperature;

        for _ in 0..self.iterations {
            for state in states.values_mut() { state.dx = 0.0; state.dy = 0.0; }

            // Repulsive forces
            for i in 0..n {
                for j in (i + 1)..n {
                    let (dx, dy) = {
                        let si = &states[&segment_names[i]];
                        let sj = &states[&segment_names[j]];
                        (sj.x - si.x, sj.y - si.y)
                    };
                    let dist = (dx * dx + dy * dy).sqrt().max(0.01);
                    let force = k * k / dist;
                    let fx = dx / dist * force;
                    let fy = dy / dist * force;

                    if let Some(si) = states.get_mut(&segment_names[i]) { si.dx -= fx; si.dy -= fy; }
                    if let Some(sj) = states.get_mut(&segment_names[j]) { sj.dx += fx; sj.dy += fy; }
                }
            }

            // Attractive forces
            for link in &graph.links {
                if let (Some(si), Some(sj)) = (states.get(&link.from_segment).cloned(), states.get(&link.to_segment).cloned()) {
                    let dx = sj.x - si.x;
                    let dy = sj.y - si.y;
                    let dist = (dx * dx + dy * dy).sqrt().max(0.01);
                    let force = dist * dist / k;
                    let fx = dx / dist * force;
                    let fy = dy / dist * force;

                    if let Some(s) = states.get_mut(&link.from_segment) { s.dx += fx; s.dy += fy; }
                    if let Some(s) = states.get_mut(&link.to_segment) { s.dx -= fx; s.dy -= fy; }
                }
            }

            // Apply displacements
            for state in states.values_mut() {
                let dist = (state.dx * state.dx + state.dy * state.dy).sqrt().max(0.01);
                let capped = dist.min(temp);
                state.x += state.dx / dist * capped;
                state.y += state.dy / dist * capped;
            }
            temp *= self.cooling;
        }

        // Create segment layouts
        for (rank, name) in segment_names.iter().enumerate() {
            if let (Some(state), Some(segment)) = (states.get(name), graph.get_segment(name)) {
                let width = (segment.sequence_length() as f32 * 0.1).max(40.0).min(200.0);
                layout.segments.insert(name.clone(), SegmentLayout {
                    name: name.clone(), x: state.x, y: state.y,
                    width, height: self.segment_height, rank,
                });
            }
        }

        // Create link layouts
        for link in &graph.links {
            if let (Some(from_seg), Some(to_seg)) = (layout.segments.get(&link.from_segment), layout.segments.get(&link.to_segment)) {
                layout.links.push(LinkLayout {
                    from: link.from_segment.clone(),
                    to: link.to_segment.clone(),
                    control_points: vec![from_seg.center(), to_seg.center()],
                });
            }
        }

        layout.recalculate_bounds();
        layout
    }

    fn name(&self) -> &'static str { "Force-Directed" }
}

