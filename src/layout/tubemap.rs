//! Tube Map layout - linear layout inspired by transit maps
//! 
//! This layout arranges segments horizontally based on topological order,
//! with parallel paths shown as parallel lines (like a transit map).

use crate::graph::{GraphGenome, Orientation};
use super::{Layout, LayoutAlgorithm, SegmentLayout, LinkLayout};
use std::collections::{HashMap, HashSet, VecDeque};

pub struct TubeMapLayout {
    segment_height: f32,
    path_spacing: f32,
    base_spacing: f32,
    min_segment_width: f32,
}

impl Default for TubeMapLayout {
    fn default() -> Self {
        Self {
            segment_height: 30.0,
            path_spacing: 40.0,
            base_spacing: 0.1,
            min_segment_width: 20.0,
        }
    }
}

impl TubeMapLayout {
    pub fn new() -> Self { Self::default() }

    fn topological_order(&self, graph: &GraphGenome) -> Vec<String> {
        let mut in_degree: HashMap<&str, usize> = HashMap::new();
        let mut order = Vec::new();
        let mut queue = VecDeque::new();

        for name in graph.segments.keys() {
            in_degree.insert(name.as_str(), 0);
        }

        for link in &graph.links {
            if link.from_orient == Orientation::Forward && link.to_orient == Orientation::Forward {
                *in_degree.entry(&link.to_segment).or_default() += 1;
            }
        }

        for (name, &degree) in &in_degree {
            if degree == 0 {
                queue.push_back(*name);
            }
        }

        while let Some(node) = queue.pop_front() {
            order.push(node.to_string());
            for link in graph.get_outgoing_links(node, Orientation::Forward) {
                if let Some(deg) = in_degree.get_mut(link.to_segment.as_str()) {
                    *deg = deg.saturating_sub(1);
                    if *deg == 0 {
                        queue.push_back(&link.to_segment);
                    }
                }
            }
        }

        for name in graph.segments.keys() {
            if !order.contains(name) {
                order.push(name.clone());
            }
        }

        order
    }

    fn assign_lanes(&self, graph: &GraphGenome, order: &[String]) -> HashMap<String, usize> {
        let mut lanes: HashMap<String, usize> = HashMap::new();
        let mut used_lanes: HashSet<usize> = HashSet::new();

        for (path_idx, path) in graph.paths.iter().enumerate() {
            for seg in &path.segments {
                if !lanes.contains_key(&seg.name) {
                    lanes.insert(seg.name.clone(), path_idx);
                    used_lanes.insert(path_idx);
                }
            }
        }

        let mut next_lane = used_lanes.len();
        for name in order {
            if !lanes.contains_key(name) {
                lanes.insert(name.clone(), next_lane);
                next_lane += 1;
            }
        }

        lanes
    }
}

impl LayoutAlgorithm for TubeMapLayout {
    fn compute(&self, graph: &GraphGenome) -> Layout {
        let mut layout = Layout::new();

        if graph.segments.is_empty() { return layout; }

        let order = self.topological_order(graph);
        let lanes = self.assign_lanes(graph, &order);

        let mut x_position = 0.0;

        for (rank, name) in order.iter().enumerate() {
            if let Some(segment) = graph.get_segment(name) {
                let lane = *lanes.get(name).unwrap_or(&0);
                let width = (segment.sequence_length() as f32 * self.base_spacing)
                    .max(self.min_segment_width);

                let seg_layout = SegmentLayout {
                    name: name.clone(),
                    x: x_position + width / 2.0,
                    y: lane as f32 * self.path_spacing,
                    width,
                    height: self.segment_height,
                    rank,
                };

                layout.segments.insert(name.clone(), seg_layout);
                x_position += width + 20.0;
            }
        }

        for link in &graph.links {
            if let (Some(from_seg), Some(to_seg)) = (
                layout.segments.get(&link.from_segment),
                layout.segments.get(&link.to_segment),
            ) {
                let start = (from_seg.right(), from_seg.y);
                let end = (to_seg.left(), to_seg.y);
                let control_offset = (end.0 - start.0).abs() * 0.3;
                let control_points = vec![
                    start,
                    (start.0 + control_offset, start.1),
                    (end.0 - control_offset, end.1),
                    end,
                ];

                layout.links.push(LinkLayout {
                    from: link.from_segment.clone(),
                    to: link.to_segment.clone(),
                    control_points,
                });
            }
        }

        layout.recalculate_bounds();
        layout
    }

    fn name(&self) -> &'static str { "Tube Map" }
}
