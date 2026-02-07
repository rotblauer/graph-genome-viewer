//! Integration tests for GFA parsing and graph operations

use graph_genome_viewer::graph::{GraphGenome, Segment, Link, Path, PathSegment, TagValue, SegmentSequence, Orientation};
use graph_genome_viewer::io::gfa::GfaLoader;
use std::io::Cursor;

// ============================================================================
// Test Data Constants
// ============================================================================

const SIMPLE_GFA: &str = r#"H	VN:Z:1.0
S	s1	ACGTACGT
S	s2	TGCATGCA
S	s3	AAAAAAAA
L	s1	+	s2	+	0M
L	s1	+	s3	+	0M
L	s2	+	s3	+	0M
P	path1	s1+,s2+,s3+	*
P	path2	s1+,s3+	*
"#;

const GFA_WITH_TAGS: &str = r#"H	VN:Z:1.0
S	seg1	ACGTACGTACGT	LN:i:12	DP:i:25	RC:i:100
S	seg2	TGCA	LN:i:4	DP:i:10
L	seg1	+	seg2	+	4M
"#;

const GFA_WITH_WALKS: &str = r#"H	VN:Z:1.1
S	s1	ACGT
S	s2	TGCA
S	s3	AAAA
L	s1	+	s2	+	0M
L	s2	+	s3	+	0M
W	sample1	0	chr1	0	100	>s1>s2>s3
W	sample1	1	chr1	0	100	>s1<s2>s3
"#;

const BUBBLE_GFA: &str = r#"H	VN:Z:1.0
S	ref1	ACGTACGTACGT
S	alt1a	TGCA
S	alt1b	AAAA
S	ref2	CCCCCCCC
L	ref1	+	alt1a	+	0M
L	ref1	+	alt1b	+	0M
L	alt1a	+	ref2	+	0M
L	alt1b	+	ref2	+	0M
P	reference	ref1+,alt1a+,ref2+	*
P	alternate	ref1+,alt1b+,ref2+	*
"#;

const COMPLEX_SV_GFA: &str = r#"H	VN:Z:1.1
S	flank5	ACGTACGTACGTACGT
S	ref_var	TGCATGCATGCA
S	alt_del	TG
S	alt_ins	TGCATGCATGCAAAAAAAAAAATGCATGCA
S	flank3	CCCCCCCCCCCCCCCC
L	flank5	+	ref_var	+	0M
L	flank5	+	alt_del	+	0M
L	flank5	+	alt_ins	+	0M
L	ref_var	+	flank3	+	0M
L	alt_del	+	flank3	+	0M
L	alt_ins	+	flank3	+	0M
W	GRCh38	0	chr1	0	44	>flank5>ref_var>flank3
W	NA19240	1	chr1	0	34	>flank5>alt_del>flank3
W	NA19240	2	chr1	0	62	>flank5>alt_ins>flank3
"#;

// ============================================================================
// GFA Parsing Tests
// ============================================================================

#[cfg(test)]
mod gfa_parsing_tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_parse_simple_gfa() {
        let reader = BufReader::new(Cursor::new(SIMPLE_GFA));
        let graph = GfaLoader::parse(reader).expect("Failed to parse simple GFA");

        assert_eq!(graph.node_count(), 3, "Should have 3 segments");
        assert_eq!(graph.edge_count(), 3, "Should have 3 links");
        assert_eq!(graph.paths.len(), 2, "Should have 2 paths");

        // Verify segment sequences
        let s1 = graph.get_segment("s1").expect("s1 should exist");
        assert_eq!(s1.sequence_length(), 8);
    }

    #[test]
    fn test_parse_gfa_with_tags() {
        let reader = BufReader::new(Cursor::new(GFA_WITH_TAGS));
        let graph = GfaLoader::parse(reader).expect("Failed to parse GFA with tags");

        let seg1 = graph.get_segment("seg1").expect("seg1 should exist");

        // Check tag parsing
        if let Some(TagValue::Integer(dp)) = seg1.get_tag("DP") {
            assert_eq!(*dp, 25, "DP tag should be 25");
        } else {
            panic!("DP tag not found or wrong type");
        }

        if let Some(TagValue::Integer(rc)) = seg1.get_tag("RC") {
            assert_eq!(*rc, 100, "RC tag should be 100");
        } else {
            panic!("RC tag not found or wrong type");
        }
    }

    #[test]
    fn test_parse_walks() {
        let reader = BufReader::new(Cursor::new(GFA_WITH_WALKS));
        let graph = GfaLoader::parse(reader).expect("Failed to parse GFA with walks");

        assert_eq!(graph.paths.len(), 2, "Should have 2 walks converted to paths");

        // Check walk metadata
        let path = &graph.paths[0];
        assert!(path.sample_name.is_some());
        assert_eq!(path.sample_name.as_deref(), Some("sample1"));
    }

    #[test]
    fn test_parse_bubble_structure() {
        let reader = BufReader::new(Cursor::new(BUBBLE_GFA));
        let graph = GfaLoader::parse(reader).expect("Failed to parse bubble GFA");

        assert_eq!(graph.node_count(), 4, "Bubble should have 4 segments");
        assert_eq!(graph.edge_count(), 4, "Bubble should have 4 links");

        // Verify bubble structure: ref1 should have 2 outgoing links
        let outgoing = graph.get_outgoing_links("ref1", Orientation::Forward);
        assert_eq!(outgoing.len(), 2, "ref1 should have 2 outgoing links");

        // ref2 should have 2 incoming links
        let incoming = graph.get_incoming_links("ref2", Orientation::Forward);
        assert_eq!(incoming.len(), 2, "ref2 should have 2 incoming links");
    }

    #[test]
    fn test_parse_complex_sv() {
        let reader = BufReader::new(Cursor::new(COMPLEX_SV_GFA));
        let graph = GfaLoader::parse(reader).expect("Failed to parse complex SV GFA");

        assert_eq!(graph.node_count(), 5);
        assert_eq!(graph.edge_count(), 6);
        assert_eq!(graph.paths.len(), 3);

        // Verify different allele lengths
        let ref_var = graph.get_segment("ref_var").unwrap();
        let alt_del = graph.get_segment("alt_del").unwrap();
        let alt_ins = graph.get_segment("alt_ins").unwrap();

        assert!(alt_del.sequence_length() < ref_var.sequence_length(), "Deletion should be shorter");
        assert!(alt_ins.sequence_length() > ref_var.sequence_length(), "Insertion should be longer");
    }

    #[test]
    fn test_empty_gfa() {
        let empty = "H\tVN:Z:1.0\n";
        let reader = BufReader::new(Cursor::new(empty));
        let graph = GfaLoader::parse(reader).expect("Failed to parse empty GFA");

        assert_eq!(graph.node_count(), 0);
        assert_eq!(graph.edge_count(), 0);
        assert_eq!(graph.paths.len(), 0);
    }

    #[test]
    fn test_gfa_with_comments() {
        let gfa_with_comments = r#"# This is a comment
H	VN:Z:1.0
# Another comment
S	s1	ACGT
"#;
        let reader = BufReader::new(Cursor::new(gfa_with_comments));
        let graph = GfaLoader::parse(reader).expect("Failed to parse GFA with comments");

        assert_eq!(graph.node_count(), 1);
    }
}

// ============================================================================
// Graph Data Structure Tests
// ============================================================================

#[cfg(test)]
mod graph_structure_tests {
    use super::*;

    #[test]
    fn test_create_empty_graph() {
        let graph = GraphGenome::new();
        assert_eq!(graph.node_count(), 0);
        assert_eq!(graph.edge_count(), 0);
    }

    #[test]
    fn test_add_segment() {
        let mut graph = GraphGenome::new();
        let segment = Segment::new("seg1", "ACGTACGT");
        graph.add_segment(segment);

        assert_eq!(graph.node_count(), 1);
        assert!(graph.get_segment("seg1").is_some());
        assert!(graph.get_segment("nonexistent").is_none());
    }

    #[test]
    fn test_add_link() {
        let mut graph = GraphGenome::new();
        graph.add_segment(Segment::new("s1", "ACGT"));
        graph.add_segment(Segment::new("s2", "TGCA"));

        let link = Link::new("s1", Orientation::Forward, "s2", Orientation::Forward, "0M");
        graph.add_link(link);

        assert_eq!(graph.edge_count(), 1);

        let links = graph.get_links_for_segment("s1");
        assert_eq!(links.len(), 1);
    }

    #[test]
    fn test_add_path() {
        let mut graph = GraphGenome::new();
        graph.add_segment(Segment::new("s1", "ACGT"));
        graph.add_segment(Segment::new("s2", "TGCA"));

        let path = Path::with_segments("path1", vec![
            PathSegment::new("s1", Orientation::Forward),
            PathSegment::new("s2", Orientation::Forward),
        ]);
        graph.add_path(path);

        assert_eq!(graph.paths.len(), 1);
        assert_eq!(graph.paths[0].name, "path1");
        assert_eq!(graph.paths[0].segments.len(), 2);
    }

    #[test]
    fn test_total_sequence_length() {
        let mut graph = GraphGenome::new();
        graph.add_segment(Segment::new("s1", "ACGT"));      // 4bp
        graph.add_segment(Segment::new("s2", "TGCATGCA")); // 8bp
        graph.add_segment(Segment::new("s3", "AA"));       // 2bp

        assert_eq!(graph.total_sequence_length(), 14);
    }

    #[test]
    fn test_outgoing_incoming_links() {
        let mut graph = GraphGenome::new();
        graph.add_segment(Segment::new("a", "ACGT"));
        graph.add_segment(Segment::new("b", "TGCA"));
        graph.add_segment(Segment::new("c", "AAAA"));

        // a -> b, a -> c
        graph.add_link(Link::new("a", Orientation::Forward, "b", Orientation::Forward, "0M"));
        graph.add_link(Link::new("a", Orientation::Forward, "c", Orientation::Forward, "0M"));

        let outgoing_a = graph.get_outgoing_links("a", Orientation::Forward);
        assert_eq!(outgoing_a.len(), 2);

        let incoming_b = graph.get_incoming_links("b", Orientation::Forward);
        assert_eq!(incoming_b.len(), 1);

        let incoming_a = graph.get_incoming_links("a", Orientation::Forward);
        assert_eq!(incoming_a.len(), 0);
    }

    #[test]
    fn test_find_sources_sinks() {
        let mut graph = GraphGenome::new();
        graph.add_segment(Segment::new("source", "ACGT"));
        graph.add_segment(Segment::new("middle", "TGCA"));
        graph.add_segment(Segment::new("sink", "AAAA"));

        graph.add_link(Link::new("source", Orientation::Forward, "middle", Orientation::Forward, "0M"));
        graph.add_link(Link::new("middle", Orientation::Forward, "sink", Orientation::Forward, "0M"));

        let sources = graph.find_sources();
        let sinks = graph.find_sinks();

        assert!(sources.contains(&"source"));
        assert!(!sources.contains(&"middle"));
        assert!(!sources.contains(&"sink"));

        assert!(sinks.contains(&"sink"));
        assert!(!sinks.contains(&"middle"));
        assert!(!sinks.contains(&"source"));
    }
}

// ============================================================================
// Segment Tests
// ============================================================================

#[cfg(test)]
mod segment_tests {
    use super::*;

    #[test]
    fn test_segment_creation() {
        let seg = Segment::new("test_seg", "ACGTACGT");
        assert_eq!(seg.name, "test_seg");
        assert_eq!(seg.sequence_length(), 8);
    }

    #[test]
    fn test_segment_gc_content() {
        // 50% GC (4 GC out of 8) - gc_content returns percentage
        let mut seg = Segment::new("gc50", "ACGTACGT");
        let gc = seg.gc_content();
        assert!((gc - 50.0).abs() < 0.01);

        // 100% GC
        let mut seg_all_gc = Segment::new("gc100", "GGGGCCCC");
        assert!((seg_all_gc.gc_content() - 100.0).abs() < 0.01);

        // 0% GC
        let mut seg_no_gc = Segment::new("gc0", "AAAATTTT");
        assert!(seg_no_gc.gc_content().abs() < 0.01);
    }

    #[test]
    fn test_segment_with_tags() {
        let mut seg = Segment::new("tagged", "ACGT");
        seg.add_tag("DP", TagValue::Integer(42));
        seg.add_tag("RC", TagValue::Integer(100));
        seg.add_tag("SN", TagValue::String("chr1".to_string()));

        if let Some(TagValue::Integer(dp)) = seg.get_tag("DP") {
            assert_eq!(*dp, 42);
        } else {
            panic!("DP tag not found");
        }

        if let Some(TagValue::String(sn)) = seg.get_tag("SN") {
            assert_eq!(sn, "chr1");
        } else {
            panic!("SN tag not found");
        }
    }

    #[test]
    fn test_segment_missing_sequence() {
        let seg = Segment::with_length("no_seq", 1000);
        assert_eq!(seg.sequence_length(), 1000);
        // Check that sequence is missing using the sequence field directly
        assert!(!seg.sequence.is_present());
    }
}

// ============================================================================
// Link Tests
// ============================================================================

#[cfg(test)]
mod link_tests {
    use super::*;

    #[test]
    fn test_link_creation() {
        let link = Link::new("s1", Orientation::Forward, "s2", Orientation::Forward, "4M");

        assert_eq!(link.from_segment, "s1");
        assert_eq!(link.to_segment, "s2");
        assert_eq!(link.from_orient, Orientation::Forward);
        assert_eq!(link.to_orient, Orientation::Forward);
        assert_eq!(link.overlap, "4M");
    }

    #[test]
    fn test_link_reverse() {
        let link = Link::new("s1", Orientation::Forward, "s2", Orientation::Reverse, "0M");
        let reversed = link.reverse();

        assert_eq!(reversed.from_segment, "s2");
        assert_eq!(reversed.to_segment, "s1");
        assert_eq!(reversed.from_orient, Orientation::Forward); // Reverse flipped
        assert_eq!(reversed.to_orient, Orientation::Reverse);   // Forward flipped
    }

    #[test]
    fn test_overlap_length() {
        let link_4m = Link::new("s1", Orientation::Forward, "s2", Orientation::Forward, "4M");
        assert_eq!(link_4m.overlap_length(), 4);

        let link_0m = Link::new("s1", Orientation::Forward, "s2", Orientation::Forward, "0M");
        assert_eq!(link_0m.overlap_length(), 0);

        let link_star = Link::new("s1", Orientation::Forward, "s2", Orientation::Forward, "*");
        assert_eq!(link_star.overlap_length(), 0);
    }
}

// ============================================================================
// Path Tests
// ============================================================================

#[cfg(test)]
mod path_tests {
    use super::*;

    #[test]
    fn test_path_creation() {
        let path = Path::new("test_path");
        assert_eq!(path.name, "test_path");
        assert!(path.segments.is_empty());
    }

    #[test]
    fn test_path_with_segments() {
        let path = Path::with_segments("path1", vec![
            PathSegment::new("s1", Orientation::Forward),
            PathSegment::new("s2", Orientation::Forward),
            PathSegment::new("s3", Orientation::Reverse),
        ]);

        assert_eq!(path.segments.len(), 3);
        assert_eq!(path.segments[2].orientation, Orientation::Reverse);
    }

    #[test]
    fn test_path_segment_parse() {
        let seg_plus = PathSegment::parse("segment1+").unwrap();
        assert_eq!(seg_plus.name, "segment1");
        assert_eq!(seg_plus.orientation, Orientation::Forward);

        let seg_minus = PathSegment::parse("segment2-").unwrap();
        assert_eq!(seg_minus.name, "segment2");
        assert_eq!(seg_minus.orientation, Orientation::Reverse);

        assert!(PathSegment::parse("").is_none());
    }

    #[test]
    fn test_walk_metadata() {
        let mut path = Path::new("NA19240#1#chr17");
        path.sample_name = Some("NA19240".to_string());
        path.haplotype = Some(1);
        path.sequence_name = Some("chr17".to_string());
        path.start = Some(10984564);
        path.end = Some(10993960);

        assert_eq!(path.sample_name.as_deref(), Some("NA19240"));
        assert_eq!(path.haplotype, Some(1));
        assert_eq!(path.sequence_name.as_deref(), Some("chr17"));
    }
}

// ============================================================================
// Orientation Tests
// ============================================================================

#[cfg(test)]
mod orientation_tests {
    use super::*;

    #[test]
    fn test_orientation_from_char() {
        assert_eq!(Orientation::from_char('+'), Some(Orientation::Forward));
        assert_eq!(Orientation::from_char('-'), Some(Orientation::Reverse));
        assert_eq!(Orientation::from_char('x'), None);
    }

    #[test]
    fn test_orientation_flip() {
        assert_eq!(Orientation::Forward.flip(), Orientation::Reverse);
        assert_eq!(Orientation::Reverse.flip(), Orientation::Forward);
    }

    #[test]
    fn test_orientation_display() {
        assert_eq!(format!("{}", Orientation::Forward), "+");
        assert_eq!(format!("{}", Orientation::Reverse), "-");
    }
}

