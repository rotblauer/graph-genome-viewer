#!/usr/bin/env python3
"""
Process pangenome data for visualization in graph-genome-viewer.

This script downloads and processes NA19240 assembly data from the Human Pangenome
Reference Consortium (HPRC) and extracts regions of interest for visualization.

Target region: chr17:10984564-10993960 (hg38) - a known structural variant region
"""

import argparse
import gzip
import hashlib
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from urllib.parse import urlparse

import requests
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# Data Sources and Configuration
# ============================================================================

@dataclass
class DataSource:
    """Configuration for a data source."""
    name: str
    url: str
    description: str
    file_type: str  # 'gfa', 'fasta', 'paf', etc.
    md5: Optional[str] = None


# HPRC year 1 pangenome data sources
HPRC_SOURCES = {
    # Minigraph-Cactus pangenome GFA (full graph)
    "hprc_mc_grch38": DataSource(
        name="hprc_mc_grch38",
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gfa.gz",
        description="HPRC Year 1 Minigraph-Cactus pangenome (GRCh38)",
        file_type="gfa"
    ),
    # Chromosome-specific GFAs (smaller, more manageable)
    "hprc_mc_chr17": DataSource(
        name="hprc_mc_chr17",
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.chroms/chr17.vg",
        description="HPRC Year 1 Minigraph-Cactus chr17 (vg format)",
        file_type="vg"
    ),
}

# Alternative: Direct sample assembly links
NA19240_ASSEMBLIES = {
    "maternal": DataSource(
        name="NA19240_maternal",
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG01978/assemblies/year1_f1_assembly_v2_genbank/HG01978.maternal.f1_assembly_v2_genbank.fa.gz",
        description="NA19240 maternal assembly (placeholder URL - update with actual)",
        file_type="fasta"
    ),
    "paternal": DataSource(
        name="NA19240_paternal",
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG01978/assemblies/year1_f1_assembly_v2_genbank/HG01978.paternal.f1_assembly_v2_genbank.fa.gz",
        description="NA19240 paternal assembly (placeholder URL - update with actual)",
        file_type="fasta"
    ),
}

# Default target region
DEFAULT_REGION = {
    "chrom": "chr17",
    "start": 10984564,
    "end": 10993960,
    "padding": 5000,  # Add flanking sequence for context
}


# ============================================================================
# Utility Functions
# ============================================================================

def download_file(url: str, output_path: Path, show_progress: bool = True) -> Path:
    """Download a file with progress bar."""
    logger.info(f"Downloading: {url}")

    response = requests.get(url, stream=True)
    response.raise_for_status()

    total_size = int(response.headers.get('content-length', 0))

    with open(output_path, 'wb') as f:
        if show_progress and total_size > 0:
            with tqdm(total=total_size, unit='iB', unit_scale=True) as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    size = f.write(chunk)
                    pbar.update(size)
        else:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

    logger.info(f"Downloaded to: {output_path}")
    return output_path


def decompress_gzip(input_path: Path, output_path: Optional[Path] = None) -> Path:
    """Decompress a gzip file."""
    if output_path is None:
        output_path = input_path.with_suffix('')

    logger.info(f"Decompressing: {input_path}")

    with gzip.open(input_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return output_path


def check_tool_available(tool: str) -> bool:
    """Check if a command-line tool is available."""
    return shutil.which(tool) is not None


def run_command(cmd: List[str], capture_output: bool = True, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command."""
    logger.debug(f"Running: {' '.join(cmd)}")
    return subprocess.run(cmd, capture_output=capture_output, check=check, text=True)


# ============================================================================
# GFA Processing Classes
# ============================================================================

@dataclass
class GFASegment:
    """Represents a GFA segment (S line)."""
    name: str
    sequence: str
    tags: Dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        if self.sequence == '*':
            return int(self.tags.get('LN', 0))
        return len(self.sequence)

    def to_gfa(self) -> str:
        """Convert to GFA format string."""
        parts = ['S', self.name, self.sequence]
        for key, value in self.tags.items():
            parts.append(f"{key}:{value}")
        return '\t'.join(parts)


@dataclass
class GFALink:
    """Represents a GFA link (L line)."""
    from_segment: str
    from_orient: str
    to_segment: str
    to_orient: str
    overlap: str
    tags: Dict[str, str] = field(default_factory=dict)

    def to_gfa(self) -> str:
        """Convert to GFA format string."""
        parts = ['L', self.from_segment, self.from_orient,
                 self.to_segment, self.to_orient, self.overlap]
        for key, value in self.tags.items():
            parts.append(f"{key}:{value}")
        return '\t'.join(parts)


@dataclass
class GFAPath:
    """Represents a GFA path (P line) or walk (W line)."""
    name: str
    segments: List[Tuple[str, str]]  # List of (segment_name, orientation)
    overlaps: Optional[str] = None
    is_walk: bool = False
    sample: Optional[str] = None
    haplotype: Optional[int] = None
    seq_name: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None

    def to_gfa(self) -> str:
        """Convert to GFA format string."""
        if self.is_walk:
            walk_str = ''.join(f"{'>' if o == '+' else '<'}{s}" for s, o in self.segments)
            return f"W\t{self.sample}\t{self.haplotype}\t{self.seq_name}\t{self.start}\t{self.end}\t{walk_str}"
        else:
            seg_str = ','.join(f"{s}{o}" for s, o in self.segments)
            overlap_str = self.overlaps if self.overlaps else '*'
            return f"P\t{self.name}\t{seg_str}\t{overlap_str}"


@dataclass
class GFAGraph:
    """In-memory representation of a GFA graph."""
    header: Optional[str] = None
    segments: Dict[str, GFASegment] = field(default_factory=dict)
    links: List[GFALink] = field(default_factory=list)
    paths: List[GFAPath] = field(default_factory=list)

    def add_segment(self, segment: GFASegment):
        self.segments[segment.name] = segment

    def add_link(self, link: GFALink):
        self.links.append(link)

    def add_path(self, path: GFAPath):
        self.paths.append(path)

    def get_connected_segments(self, segment_names: Set[str]) -> Set[str]:
        """Get all segments connected to the given segments."""
        connected = set(segment_names)

        for link in self.links:
            if link.from_segment in connected:
                connected.add(link.to_segment)
            if link.to_segment in connected:
                connected.add(link.from_segment)

        return connected

    def extract_subgraph(self, segment_names: Set[str], include_paths: bool = True) -> 'GFAGraph':
        """Extract a subgraph containing only the specified segments."""
        subgraph = GFAGraph()
        subgraph.header = self.header

        # Add segments
        for name in segment_names:
            if name in self.segments:
                subgraph.add_segment(self.segments[name])

        # Add links between included segments
        for link in self.links:
            if link.from_segment in segment_names and link.to_segment in segment_names:
                subgraph.add_link(link)

        # Add paths that traverse the included segments
        if include_paths:
            for path in self.paths:
                path_segments = {s for s, _ in path.segments}
                if path_segments & segment_names:  # If any overlap
                    # Filter path to only included segments
                    filtered_segments = [(s, o) for s, o in path.segments if s in segment_names]
                    if filtered_segments:
                        new_path = GFAPath(
                            name=path.name,
                            segments=filtered_segments,
                            is_walk=path.is_walk,
                            sample=path.sample,
                            haplotype=path.haplotype,
                            seq_name=path.seq_name,
                            start=path.start,
                            end=path.end
                        )
                        subgraph.add_path(new_path)

        return subgraph

    def to_gfa(self) -> str:
        """Convert graph to GFA format string."""
        lines = []

        # Header
        if self.header:
            lines.append(self.header)
        else:
            lines.append("H\tVN:Z:1.0")

        # Segments (sorted for reproducibility)
        for name in sorted(self.segments.keys()):
            lines.append(self.segments[name].to_gfa())

        # Links
        for link in self.links:
            lines.append(link.to_gfa())

        # Paths
        for path in self.paths:
            lines.append(path.to_gfa())

        return '\n'.join(lines) + '\n'

    def write(self, output_path: Path):
        """Write graph to file."""
        with open(output_path, 'w') as f:
            f.write(self.to_gfa())
        logger.info(f"Wrote GFA to: {output_path}")

    def stats(self) -> Dict[str, int]:
        """Return basic statistics about the graph."""
        total_seq_len = sum(s.length for s in self.segments.values())
        return {
            'segments': len(self.segments),
            'links': len(self.links),
            'paths': len(self.paths),
            'total_sequence_length': total_seq_len,
        }


# ============================================================================
# GFA Parsing
# ============================================================================

def parse_gfa_line(line: str) -> Optional[tuple]:
    """Parse a single GFA line."""
    line = line.strip()
    if not line or line.startswith('#'):
        return None

    fields = line.split('\t')
    if not fields:
        return None

    record_type = fields[0]

    if record_type == 'H':
        return ('H', line)

    elif record_type == 'S':
        # Segment: S <name> <sequence> [tags...]
        if len(fields) < 3:
            return None
        name = fields[1]
        sequence = fields[2]
        tags = {}
        for field in fields[3:]:
            if ':' in field:
                parts = field.split(':', 2)
                if len(parts) >= 2:
                    tags[parts[0]] = ':'.join(parts[1:])
        return ('S', GFASegment(name=name, sequence=sequence, tags=tags))

    elif record_type == 'L':
        # Link: L <from> <from_orient> <to> <to_orient> <overlap>
        if len(fields) < 6:
            return None
        tags = {}
        for field in fields[6:]:
            if ':' in field:
                parts = field.split(':', 2)
                if len(parts) >= 2:
                    tags[parts[0]] = ':'.join(parts[1:])
        return ('L', GFALink(
            from_segment=fields[1],
            from_orient=fields[2],
            to_segment=fields[3],
            to_orient=fields[4],
            overlap=fields[5],
            tags=tags
        ))

    elif record_type == 'P':
        # Path: P <name> <segment_names> <overlaps>
        if len(fields) < 3:
            return None
        name = fields[1]
        seg_str = fields[2]
        overlaps = fields[3] if len(fields) > 3 else None

        segments = []
        for seg in seg_str.split(','):
            seg = seg.strip()
            if seg:
                if seg.endswith('+') or seg.endswith('-'):
                    segments.append((seg[:-1], seg[-1]))
                else:
                    segments.append((seg, '+'))

        return ('P', GFAPath(name=name, segments=segments, overlaps=overlaps))

    elif record_type == 'W':
        # Walk: W <sample> <haplotype> <seqname> <start> <end> <walk>
        if len(fields) < 7:
            return None
        sample = fields[1]
        haplotype = int(fields[2])
        seq_name = fields[3]
        start = int(fields[4])
        end = int(fields[5])
        walk_str = fields[6]

        # Parse walk string: >seg1>seg2<seg3...
        segments = []
        current_seg = ""
        current_orient = "+"
        for char in walk_str:
            if char == '>':
                if current_seg:
                    segments.append((current_seg, current_orient))
                current_seg = ""
                current_orient = "+"
            elif char == '<':
                if current_seg:
                    segments.append((current_seg, current_orient))
                current_seg = ""
                current_orient = "-"
            else:
                current_seg += char
        if current_seg:
            segments.append((current_seg, current_orient))

        path_name = f"{sample}#{haplotype}#{seq_name}"
        return ('W', GFAPath(
            name=path_name,
            segments=segments,
            is_walk=True,
            sample=sample,
            haplotype=haplotype,
            seq_name=seq_name,
            start=start,
            end=end
        ))

    return None


def parse_gfa_file(input_path: Path) -> GFAGraph:
    """Parse a GFA file into a GFAGraph object."""
    logger.info(f"Parsing GFA file: {input_path}")

    graph = GFAGraph()

    # Handle gzipped files
    if str(input_path).endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

    with opener(input_path, mode) as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 100000 == 0:
                logger.debug(f"Processed {line_num} lines...")

            result = parse_gfa_line(line)
            if result is None:
                continue

            record_type, data = result

            if record_type == 'H':
                graph.header = data
            elif record_type == 'S':
                graph.add_segment(data)
            elif record_type == 'L':
                graph.add_link(data)
            elif record_type in ('P', 'W'):
                graph.add_path(data)

    stats = graph.stats()
    logger.info(f"Parsed GFA: {stats['segments']} segments, {stats['links']} links, {stats['paths']} paths")

    return graph


# ============================================================================
# Region Extraction
# ============================================================================

def find_segments_in_region(
    graph: GFAGraph,
    chrom: str,
    start: int,
    end: int,
    reference_path: Optional[str] = None
) -> Set[str]:
    """
    Find all segments that overlap with a genomic region.

    This uses the path/walk annotations to map coordinates to segments.
    """
    logger.info(f"Finding segments in region: {chrom}:{start}-{end}")

    matching_segments = set()

    # Find paths that contain the chromosome
    for path in graph.paths:
        # Check if this path is for our chromosome
        path_chrom = None
        if path.seq_name:
            path_chrom = path.seq_name
        elif path.name:
            # Try to extract chromosome from path name
            if chrom in path.name:
                path_chrom = chrom

        if path_chrom and chrom in path_chrom:
            # This path is relevant
            # For walks with coordinates, check overlap
            if path.start is not None and path.end is not None:
                if path.end < start or path.start > end:
                    continue

            # Add all segments from this path
            for seg_name, _ in path.segments:
                matching_segments.add(seg_name)

    # If no paths matched, try reference path name patterns
    if not matching_segments and reference_path:
        for path in graph.paths:
            if reference_path in path.name:
                for seg_name, _ in path.segments:
                    matching_segments.add(seg_name)

    logger.info(f"Found {len(matching_segments)} segments in region")
    return matching_segments


def extract_region_subgraph(
    graph: GFAGraph,
    chrom: str,
    start: int,
    end: int,
    expand_hops: int = 1
) -> GFAGraph:
    """
    Extract a subgraph for a specific genomic region.

    Args:
        graph: The full GFA graph
        chrom: Chromosome name (e.g., 'chr17')
        start: Start position
        end: End position
        expand_hops: Number of hops to expand the subgraph by (for context)

    Returns:
        A new GFAGraph containing only the relevant segments
    """
    # Find initial segments
    segment_names = find_segments_in_region(graph, chrom, start, end)

    # Expand by connected segments
    for _ in range(expand_hops):
        segment_names = graph.get_connected_segments(segment_names)

    logger.info(f"Extracting subgraph with {len(segment_names)} segments (after {expand_hops} expansion hops)")

    return graph.extract_subgraph(segment_names)


# ============================================================================
# Sample Data Generation
# ============================================================================

def create_sample_sv_graph(output_path: Path) -> GFAGraph:
    """
    Create a sample GFA graph representing a structural variant.

    This creates a synthetic graph similar to what would be found at
    chr17:10984564-10993960 with a complex structural variant.
    """
    logger.info("Creating sample structural variant graph")

    graph = GFAGraph()
    graph.header = "H\tVN:Z:1.0\tbn:Z:NA19240_chr17_sv_sample"

    # Create a graph structure with:
    # - Reference backbone segments
    # - An alternate path with a deletion
    # - An alternate path with an insertion
    # - Complex nested structure

    # Reference-like sequences (simplified, using placeholder sequences)
    ref_segments = [
        GFASegment("s1_ref_flank5", "A" * 500, {"SN:Z": "chr17", "SO:i": "10984564"}),
        GFASegment("s2_ref_pre", "C" * 300, {"SN:Z": "chr17", "SO:i": "10985064"}),
        GFASegment("s3_ref_var", "G" * 200, {"SN:Z": "chr17", "SO:i": "10985364"}),  # Variant region
        GFASegment("s4_ref_post", "T" * 300, {"SN:Z": "chr17", "SO:i": "10985564"}),
        GFASegment("s5_ref_flank3", "A" * 500, {"SN:Z": "chr17", "SO:i": "10985864"}),
    ]

    # Alternate allele segments
    alt_segments = [
        GFASegment("s3_alt1_del", "G" * 50, {"VT:Z": "DEL"}),  # Deletion variant
        GFASegment("s3_alt2_ins", "G" * 200 + "ACGTACGTACGT" * 10 + "G" * 200, {"VT:Z": "INS"}),  # Insertion variant
        GFASegment("s3_alt3_snv", "A" * 200, {"VT:Z": "SNV"}),  # SNV variant
    ]

    # Complex nested variant
    nested_segments = [
        GFASegment("s3_nest_a", "CAGT" * 50, {"VT:Z": "COMPLEX"}),
        GFASegment("s3_nest_b", "TCGA" * 25, {"VT:Z": "COMPLEX"}),
        GFASegment("s3_nest_c", "ATAT" * 25, {"VT:Z": "COMPLEX"}),
    ]

    # Add all segments
    for seg in ref_segments + alt_segments + nested_segments:
        graph.add_segment(seg)

    # Reference backbone links
    graph.add_link(GFALink("s1_ref_flank5", "+", "s2_ref_pre", "+", "0M"))
    graph.add_link(GFALink("s2_ref_pre", "+", "s3_ref_var", "+", "0M"))
    graph.add_link(GFALink("s3_ref_var", "+", "s4_ref_post", "+", "0M"))
    graph.add_link(GFALink("s4_ref_post", "+", "s5_ref_flank3", "+", "0M"))

    # Alternate allele links (bubble structure)
    # Deletion path
    graph.add_link(GFALink("s2_ref_pre", "+", "s3_alt1_del", "+", "0M"))
    graph.add_link(GFALink("s3_alt1_del", "+", "s4_ref_post", "+", "0M"))

    # Insertion path
    graph.add_link(GFALink("s2_ref_pre", "+", "s3_alt2_ins", "+", "0M"))
    graph.add_link(GFALink("s3_alt2_ins", "+", "s4_ref_post", "+", "0M"))

    # SNV path
    graph.add_link(GFALink("s2_ref_pre", "+", "s3_alt3_snv", "+", "0M"))
    graph.add_link(GFALink("s3_alt3_snv", "+", "s4_ref_post", "+", "0M"))

    # Complex nested variant (within the variant region)
    graph.add_link(GFALink("s2_ref_pre", "+", "s3_nest_a", "+", "0M"))
    graph.add_link(GFALink("s3_nest_a", "+", "s3_nest_b", "+", "0M"))
    graph.add_link(GFALink("s3_nest_a", "+", "s3_nest_c", "+", "0M"))
    graph.add_link(GFALink("s3_nest_b", "+", "s4_ref_post", "+", "0M"))
    graph.add_link(GFALink("s3_nest_c", "+", "s4_ref_post", "+", "0M"))

    # Add paths representing different haplotypes/samples
    paths = [
        GFAPath("GRCh38#0#chr17", [
            ("s1_ref_flank5", "+"), ("s2_ref_pre", "+"), ("s3_ref_var", "+"),
            ("s4_ref_post", "+"), ("s5_ref_flank3", "+")
        ], is_walk=True, sample="GRCh38", haplotype=0, seq_name="chr17", start=10984564, end=10986364),

        GFAPath("NA19240#1#chr17", [
            ("s1_ref_flank5", "+"), ("s2_ref_pre", "+"), ("s3_alt1_del", "+"),
            ("s4_ref_post", "+"), ("s5_ref_flank3", "+")
        ], is_walk=True, sample="NA19240", haplotype=1, seq_name="chr17", start=10984564, end=10986214),

        GFAPath("NA19240#2#chr17", [
            ("s1_ref_flank5", "+"), ("s2_ref_pre", "+"), ("s3_alt2_ins", "+"),
            ("s4_ref_post", "+"), ("s5_ref_flank3", "+")
        ], is_walk=True, sample="NA19240", haplotype=2, seq_name="chr17", start=10984564, end=10986484),

        GFAPath("HG002#1#chr17", [
            ("s1_ref_flank5", "+"), ("s2_ref_pre", "+"), ("s3_alt3_snv", "+"),
            ("s4_ref_post", "+"), ("s5_ref_flank3", "+")
        ], is_walk=True, sample="HG002", haplotype=1, seq_name="chr17", start=10984564, end=10986364),

        GFAPath("HG005#1#chr17", [
            ("s1_ref_flank5", "+"), ("s2_ref_pre", "+"), ("s3_nest_a", "+"),
            ("s3_nest_b", "+"), ("s4_ref_post", "+"), ("s5_ref_flank3", "+")
        ], is_walk=True, sample="HG005", haplotype=1, seq_name="chr17", start=10984564, end=10986364),
    ]

    for path in paths:
        graph.add_path(path)

    # Write to file
    graph.write(output_path)

    return graph


def create_realistic_chr17_graph(output_path: Path) -> GFAGraph:
    """
    Create a more realistic chr17 region graph with actual sequence patterns.

    This simulates a structural variant region with more realistic features.
    """
    import random
    random.seed(42)  # For reproducibility

    logger.info("Creating realistic chr17 structural variant graph")

    graph = GFAGraph()
    graph.header = "H\tVN:Z:1.1\tbn:Z:NA19240_chr17_10984564_10993960"

    def random_seq(length: int) -> str:
        """Generate random DNA sequence."""
        return ''.join(random.choices('ACGT', k=length))

    # Create a more complex graph with multiple samples
    # Approximate the structure at chr17:10984564-10993960

    # Backbone segments (representing ~9.4kb region)
    segments = []

    # 5' flanking region (shared by all)
    segments.append(GFASegment("flank5_1", random_seq(1000), {"SN:Z": "chr17", "SO:i": "10984564", "RC:i": "50"}))
    segments.append(GFASegment("flank5_2", random_seq(500), {"SN:Z": "chr17", "SO:i": "10985564", "RC:i": "48"}))

    # Variable region - this is where the structural variants occur
    # Reference path
    segments.append(GFASegment("var_ref_1", random_seq(800), {"SN:Z": "chr17", "SO:i": "10986064", "RC:i": "25"}))
    segments.append(GFASegment("var_ref_2", random_seq(600), {"SN:Z": "chr17", "SO:i": "10986864", "RC:i": "24"}))
    segments.append(GFASegment("var_ref_3", random_seq(700), {"SN:Z": "chr17", "SO:i": "10987464", "RC:i": "26"}))

    # Alternate allele 1: Large deletion (~1.5kb deleted)
    segments.append(GFASegment("var_alt1_del", random_seq(100), {"VT:Z": "DEL", "RC:i": "12"}))

    # Alternate allele 2: Large insertion (~2kb inserted)
    segments.append(GFASegment("var_alt2_ins", random_seq(2000), {"VT:Z": "INS", "RC:i": "8"}))

    # Alternate allele 3: Complex rearrangement with inversion
    segments.append(GFASegment("var_alt3_a", random_seq(300), {"VT:Z": "COMPLEX", "RC:i": "5"}))
    segments.append(GFASegment("var_alt3_b", random_seq(400), {"VT:Z": "COMPLEX", "RC:i": "5"}))
    segments.append(GFASegment("var_alt3_c", random_seq(250), {"VT:Z": "COMPLEX", "RC:i": "5"}))

    # Nested bubble within var_ref_2
    segments.append(GFASegment("var_ref_2_alt", random_seq(580), {"VT:Z": "SNV", "RC:i": "10"}))

    # 3' flanking region (shared by all)
    segments.append(GFASegment("var_post", random_seq(500), {"SN:Z": "chr17", "SO:i": "10988164", "RC:i": "50"}))
    segments.append(GFASegment("flank3_1", random_seq(1000), {"SN:Z": "chr17", "SO:i": "10988664", "RC:i": "52"}))
    segments.append(GFASegment("flank3_2", random_seq(800), {"SN:Z": "chr17", "SO:i": "10989664", "RC:i": "49"}))

    for seg in segments:
        graph.add_segment(seg)

    # Build links
    links = [
        # Reference backbone
        ("flank5_1", "+", "flank5_2", "+"),
        ("flank5_2", "+", "var_ref_1", "+"),
        ("var_ref_1", "+", "var_ref_2", "+"),
        ("var_ref_2", "+", "var_ref_3", "+"),
        ("var_ref_3", "+", "var_post", "+"),
        ("var_post", "+", "flank3_1", "+"),
        ("flank3_1", "+", "flank3_2", "+"),

        # Alternate allele 1: Deletion (skips var_ref_1, var_ref_2, var_ref_3)
        ("flank5_2", "+", "var_alt1_del", "+"),
        ("var_alt1_del", "+", "var_post", "+"),

        # Alternate allele 2: Insertion
        ("flank5_2", "+", "var_alt2_ins", "+"),
        ("var_alt2_ins", "+", "var_ref_2", "+"),

        # Alternate allele 3: Complex rearrangement
        ("flank5_2", "+", "var_alt3_a", "+"),
        ("var_alt3_a", "+", "var_alt3_b", "-"),  # Inversion
        ("var_alt3_a", "+", "var_alt3_c", "+"),  # Branch
        ("var_alt3_b", "-", "var_ref_3", "+"),
        ("var_alt3_c", "+", "var_ref_3", "+"),

        # Nested bubble (SNV within var_ref_2)
        ("var_ref_1", "+", "var_ref_2_alt", "+"),
        ("var_ref_2_alt", "+", "var_ref_3", "+"),
    ]

    for from_seg, from_orient, to_seg, to_orient in links:
        graph.add_link(GFALink(from_seg, from_orient, to_seg, to_orient, "0M"))

    # Add sample paths
    sample_paths = [
        # Reference (GRCh38)
        ("GRCh38", 0, "chr17", ["flank5_1+", "flank5_2+", "var_ref_1+", "var_ref_2+", "var_ref_3+", "var_post+", "flank3_1+", "flank3_2+"]),

        # NA19240 haplotype 1 - has the deletion
        ("NA19240", 1, "chr17", ["flank5_1+", "flank5_2+", "var_alt1_del+", "var_post+", "flank3_1+", "flank3_2+"]),

        # NA19240 haplotype 2 - has the insertion
        ("NA19240", 2, "chr17", ["flank5_1+", "flank5_2+", "var_alt2_ins+", "var_ref_2+", "var_ref_3+", "var_post+", "flank3_1+", "flank3_2+"]),

        # HG002 haplotype 1 - reference-like with SNV
        ("HG002", 1, "chr17", ["flank5_1+", "flank5_2+", "var_ref_1+", "var_ref_2_alt+", "var_ref_3+", "var_post+", "flank3_1+", "flank3_2+"]),

        # HG002 haplotype 2 - reference
        ("HG002", 2, "chr17", ["flank5_1+", "flank5_2+", "var_ref_1+", "var_ref_2+", "var_ref_3+", "var_post+", "flank3_1+", "flank3_2+"]),

        # HG005 - complex rearrangement
        ("HG005", 1, "chr17", ["flank5_1+", "flank5_2+", "var_alt3_a+", "var_alt3_b-", "var_ref_3+", "var_post+", "flank3_1+", "flank3_2+"]),

        # HG005 haplotype 2 - other branch of complex
        ("HG005", 2, "chr17", ["flank5_1+", "flank5_2+", "var_alt3_a+", "var_alt3_c+", "var_ref_3+", "var_post+", "flank3_1+", "flank3_2+"]),
    ]

    for sample, hap, seq_name, path_segs in sample_paths:
        segments_parsed = []
        for seg in path_segs:
            if seg.endswith('+'):
                segments_parsed.append((seg[:-1], '+'))
            elif seg.endswith('-'):
                segments_parsed.append((seg[:-1], '-'))
            else:
                segments_parsed.append((seg, '+'))

        path = GFAPath(
            name=f"{sample}#{hap}#{seq_name}",
            segments=segments_parsed,
            is_walk=True,
            sample=sample,
            haplotype=hap,
            seq_name=seq_name,
            start=10984564,
            end=10993960
        )
        graph.add_path(path)

    graph.write(output_path)

    stats = graph.stats()
    logger.info(f"Created graph: {stats}")

    return graph


# ============================================================================
# HPRC Data Processing
# ============================================================================

def download_hprc_chr17(output_dir: Path) -> Optional[Path]:
    """
    Download chr17 data from HPRC.

    Note: The full pangenome files are very large. This function provides
    options for different data sources.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Try to download chr17-specific GFA if available
    # HPRC provides chromosome-specific files in some releases

    chr17_urls = [
        # Year 1 release minigraph-cactus by chromosome
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.chroms/chr17.gfa.gz",
        # Alternative: pggb output
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr17.hprc-v1.0-pggb.gfa.gz",
    ]

    for url in chr17_urls:
        filename = Path(urlparse(url).path).name
        output_path = output_dir / filename

        try:
            logger.info(f"Attempting to download: {url}")
            download_file(url, output_path)
            return output_path
        except Exception as e:
            logger.warning(f"Failed to download {url}: {e}")
            continue

    logger.error("Could not download chr17 data from any source")
    return None


def extract_region_from_hprc(
    gfa_path: Path,
    output_path: Path,
    chrom: str,
    start: int,
    end: int
) -> Optional[Path]:
    """
    Extract a specific region from an HPRC GFA file.
    """
    logger.info(f"Extracting region {chrom}:{start}-{end} from {gfa_path}")

    # Parse the GFA
    graph = parse_gfa_file(gfa_path)

    # Extract subgraph for region
    subgraph = extract_region_subgraph(graph, chrom, start, end, expand_hops=2)

    if not subgraph.segments:
        logger.error("No segments found in the specified region")
        return None

    # Write output
    subgraph.write(output_path)

    return output_path


# ============================================================================
# Main CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Process pangenome data for graph-genome-viewer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate sample data for testing
  %(prog)s generate-sample --output data/sample_sv.gfa

  # Generate realistic chr17 region data
  %(prog)s generate-realistic --output data/na19240_chr17.gfa

  # Download and process HPRC chr17 data (requires internet)
  %(prog)s download-hprc --output-dir data/hprc

  # Extract a region from an existing GFA file
  %(prog)s extract-region --input full.gfa --output region.gfa \\
      --chrom chr17 --start 10984564 --end 10993960

  # Parse and display stats for a GFA file
  %(prog)s stats --input data/sample.gfa
        """
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Generate sample data
    gen_sample = subparsers.add_parser("generate-sample", help="Generate sample SV graph")
    gen_sample.add_argument("--output", "-o", type=Path, required=True,
                           help="Output GFA file path")

    # Generate realistic data
    gen_real = subparsers.add_parser("generate-realistic", help="Generate realistic chr17 graph")
    gen_real.add_argument("--output", "-o", type=Path, required=True,
                         help="Output GFA file path")

    # Download HPRC data
    dl_hprc = subparsers.add_parser("download-hprc", help="Download HPRC chr17 data")
    dl_hprc.add_argument("--output-dir", "-o", type=Path, required=True,
                        help="Output directory for downloaded files")

    # Extract region
    extract = subparsers.add_parser("extract-region", help="Extract a genomic region from GFA")
    extract.add_argument("--input", "-i", type=Path, required=True,
                        help="Input GFA file")
    extract.add_argument("--output", "-o", type=Path, required=True,
                        help="Output GFA file")
    extract.add_argument("--chrom", "-c", type=str, default="chr17",
                        help="Chromosome name (default: chr17)")
    extract.add_argument("--start", "-s", type=int, default=10984564,
                        help="Start position (default: 10984564)")
    extract.add_argument("--end", "-e", type=int, default=10993960,
                        help="End position (default: 10993960)")
    extract.add_argument("--expand", type=int, default=2,
                        help="Number of hops to expand the subgraph (default: 2)")

    # Stats command
    stats_cmd = subparsers.add_parser("stats", help="Display GFA file statistics")
    stats_cmd.add_argument("--input", "-i", type=Path, required=True,
                          help="Input GFA file")

    # Validate command
    validate = subparsers.add_parser("validate", help="Validate a GFA file")
    validate.add_argument("--input", "-i", type=Path, required=True,
                         help="Input GFA file to validate")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 0

    try:
        if args.command == "generate-sample":
            args.output.parent.mkdir(parents=True, exist_ok=True)
            create_sample_sv_graph(args.output)
            print(f"Generated sample graph: {args.output}")

        elif args.command == "generate-realistic":
            args.output.parent.mkdir(parents=True, exist_ok=True)
            create_realistic_chr17_graph(args.output)
            print(f"Generated realistic graph: {args.output}")

        elif args.command == "download-hprc":
            result = download_hprc_chr17(args.output_dir)
            if result:
                print(f"Downloaded HPRC data to: {result}")
            else:
                print("Failed to download HPRC data")
                return 1

        elif args.command == "extract-region":
            result = extract_region_from_hprc(
                args.input, args.output,
                args.chrom, args.start, args.end
            )
            if result:
                print(f"Extracted region to: {result}")
            else:
                print("Failed to extract region")
                return 1

        elif args.command == "stats":
            graph = parse_gfa_file(args.input)
            stats = graph.stats()
            print(f"\nGFA Statistics for: {args.input}")
            print("-" * 40)
            print(f"Segments:               {stats['segments']:,}")
            print(f"Links:                  {stats['links']:,}")
            print(f"Paths:                  {stats['paths']:,}")
            print(f"Total sequence length:  {stats['total_sequence_length']:,} bp")

        elif args.command == "validate":
            graph = parse_gfa_file(args.input)
            stats = graph.stats()

            errors = []
            warnings = []

            # Check for orphan segments (no links)
            linked_segments = set()
            for link in graph.links:
                linked_segments.add(link.from_segment)
                linked_segments.add(link.to_segment)

            orphans = set(graph.segments.keys()) - linked_segments
            if orphans and len(graph.segments) > 1:
                warnings.append(f"Found {len(orphans)} orphan segments (no links): {list(orphans)[:5]}...")

            # Check for links to non-existent segments
            for link in graph.links:
                if link.from_segment not in graph.segments:
                    errors.append(f"Link references non-existent segment: {link.from_segment}")
                if link.to_segment not in graph.segments:
                    errors.append(f"Link references non-existent segment: {link.to_segment}")

            # Check paths
            for path in graph.paths:
                for seg_name, _ in path.segments:
                    if seg_name not in graph.segments:
                        errors.append(f"Path '{path.name}' references non-existent segment: {seg_name}")

            print(f"\nValidation Results for: {args.input}")
            print("-" * 40)
            print(f"Segments: {stats['segments']}, Links: {stats['links']}, Paths: {stats['paths']}")

            if errors:
                print(f"\nErrors ({len(errors)}):")
                for err in errors[:10]:
                    print(f"  ❌ {err}")
                if len(errors) > 10:
                    print(f"  ... and {len(errors) - 10} more errors")

            if warnings:
                print(f"\nWarnings ({len(warnings)}):")
                for warn in warnings[:10]:
                    print(f"  ⚠️  {warn}")

            if not errors and not warnings:
                print("\n✅ GFA file is valid!")

            return 1 if errors else 0

    except Exception as e:
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

