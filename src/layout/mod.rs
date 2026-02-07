//! Graph layout algorithms for visualization
mod tubemap;
mod force;
pub use tubemap::TubeMapLayout;
pub use force::ForceDirectedLayout;
use crate::graph::GraphGenome;
use std::collections::HashMap;
#[derive(Clone, Debug, Default)]
pub struct Layout {
    pub segments: HashMap<String, SegmentLayout>,
    pub links: Vec<LinkLayout>,
    pub bounds: LayoutBounds,
}
#[derive(Clone, Debug)]
pub struct SegmentLayout {
    pub name: String,
    pub x: f32,
    pub y: f32,
    pub width: f32,
    pub height: f32,
    pub rank: usize,
}
impl SegmentLayout {
    pub fn left(&self) -> f32 { self.x - self.width / 2.0 }
    pub fn right(&self) -> f32 { self.x + self.width / 2.0 }
    pub fn top(&self) -> f32 { self.y - self.height / 2.0 }
    pub fn bottom(&self) -> f32 { self.y + self.height / 2.0 }
    pub fn center(&self) -> (f32, f32) { (self.x, self.y) }
}
#[derive(Clone, Debug)]
pub struct LinkLayout {
    pub from: String,
    pub to: String,
    pub control_points: Vec<(f32, f32)>,
}
#[derive(Clone, Debug, Default)]
pub struct LayoutBounds {
    pub min_x: f32,
    pub min_y: f32,
    pub max_x: f32,
    pub max_y: f32,
}
impl LayoutBounds {
    pub fn width(&self) -> f32 { self.max_x - self.min_x }
    pub fn height(&self) -> f32 { self.max_y - self.min_y }
    pub fn center(&self) -> (f32, f32) {
        ((self.min_x + self.max_x) / 2.0, (self.min_y + self.max_y) / 2.0)
    }
    pub fn include(&mut self, x: f32, y: f32) {
        self.min_x = self.min_x.min(x);
        self.min_y = self.min_y.min(y);
        self.max_x = self.max_x.max(x);
        self.max_y = self.max_y.max(y);
    }
}
pub trait LayoutAlgorithm {
    fn compute(&self, graph: &GraphGenome) -> Layout;
    fn name(&self) -> &'static str;
}
impl Layout {
    pub fn new() -> Self { Self::default() }
    pub fn recalculate_bounds(&mut self) {
        self.bounds = LayoutBounds {
            min_x: f32::MAX, min_y: f32::MAX,
            max_x: f32::MIN, max_y: f32::MIN,
        };
        for segment in self.segments.values() {
            self.bounds.include(segment.left(), segment.top());
            self.bounds.include(segment.right(), segment.bottom());
        }
        self.bounds.min_x -= 50.0;
        self.bounds.min_y -= 50.0;
        self.bounds.max_x += 50.0;
        self.bounds.max_y += 50.0;
    }
    pub fn get_segment(&self, name: &str) -> Option<&SegmentLayout> {
        self.segments.get(name)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GraphGenome, Segment, Link, Orientation};

    fn create_test_graph() -> GraphGenome {
        let mut graph = GraphGenome::new();
        graph.add_segment(Segment::new("s1", "ACGTACGT"));
        graph.add_segment(Segment::new("s2", "TGCATGCA"));
        graph.add_segment(Segment::new("s3", "AAAAAAAA"));
        graph.add_link(Link::new("s1", Orientation::Forward, "s2", Orientation::Forward, "0M"));
        graph.add_link(Link::new("s2", Orientation::Forward, "s3", Orientation::Forward, "0M"));
        graph
    }

    #[test]
    fn test_tubemap_layout() {
        let graph = create_test_graph();
        let layout = TubeMapLayout::new().compute(&graph);

        assert_eq!(layout.segments.len(), 3);
        assert_eq!(layout.links.len(), 2);

        // All segments should be laid out
        assert!(layout.get_segment("s1").is_some());
        assert!(layout.get_segment("s2").is_some());
        assert!(layout.get_segment("s3").is_some());
    }

    #[test]
    fn test_force_layout() {
        let graph = create_test_graph();
        let layout = ForceDirectedLayout::new().compute(&graph);

        assert_eq!(layout.segments.len(), 3);
        assert_eq!(layout.links.len(), 2);
    }

    #[test]
    fn test_layout_bounds() {
        let mut layout = Layout::new();
        layout.segments.insert("test".to_string(), SegmentLayout {
            name: "test".to_string(),
            x: 100.0,
            y: 50.0,
            width: 40.0,
            height: 20.0,
            rank: 0,
        });
        layout.recalculate_bounds();

        assert!(layout.bounds.width() > 0.0);
        assert!(layout.bounds.height() > 0.0);
    }

    #[test]
    fn test_segment_layout_edges() {
        let seg = SegmentLayout {
            name: "test".to_string(),
            x: 100.0,
            y: 50.0,
            width: 40.0,
            height: 20.0,
            rank: 0,
        };

        assert_eq!(seg.left(), 80.0);
        assert_eq!(seg.right(), 120.0);
        assert_eq!(seg.top(), 40.0);
        assert_eq!(seg.bottom(), 60.0);
        assert_eq!(seg.center(), (100.0, 50.0));
    }

    #[test]
    fn test_empty_graph_layout() {
        let graph = GraphGenome::new();
        let layout = TubeMapLayout::new().compute(&graph);

        assert_eq!(layout.segments.len(), 0);
        assert_eq!(layout.links.len(), 0);
    }
}
