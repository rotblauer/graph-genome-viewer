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
