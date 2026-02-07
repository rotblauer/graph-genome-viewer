#!/usr/bin/env python3
"""
Download real NA19240 pangenome data from public repositories.

This script downloads actual NA19240 data from the Human Pangenome Reference
Consortium (HPRC) and related sources.

Usage:
    cd scripts
    source .venv/bin/activate  # or however your venv is set up
    python download_na19240_data.py --output ../data/

Data Sources:
- HPRC Year 1 Pangenome release (Minigraph-Cactus)
- 1000 Genomes Project NA19240 variants
- HPRC sample-specific assemblies
"""

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict
from urllib.parse import urlparse

try:
    import requests
    from tqdm import tqdm
except ImportError:
    print("Missing required packages. Install with:")
    print("  pip install requests tqdm")
    sys.exit(1)


# ============================================================================
# Configuration: Real HPRC and 1000 Genomes Data URLs
# ============================================================================

@dataclass
class DataFile:
    """Information about a downloadable data file."""
    name: str
    url: str
    description: str
    size_mb: Optional[float] = None
    md5: Optional[str] = None


# HPRC Year 1 Minigraph-Cactus pangenome data
# Source: https://github.com/human-pangenomics/HPP_Year1_Assemblies
HPRC_DATA = {
    # Full pangenome GFA (very large, ~150GB uncompressed)
    "hprc_full": DataFile(
        name="hprc-v1.0-mc-grch38.gfa.gz",
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gfa.gz",
        description="HPRC Year 1 Full Minigraph-Cactus pangenome (GRCh38 coordinates)",
        size_mb=15000,  # ~15GB compressed
    ),

    # Chromosome-specific VG files (smaller, can be converted to GFA)
    "chr17_vg": DataFile(
        name="chr17.vg",
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.chroms/chr17.vg",
        description="HPRC chr17 in VG format (requires vg tools to convert to GFA)",
        size_mb=800,
    ),
}

# HPRC PGGB (Progressive Graphical Pangenome Builder) data
# These are often available as GFA directly
PGGB_DATA = {
    "chr17_pggb": DataFile(
        name="chr17.hprc-v1.0-pggb.gfa.gz",
        url="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr17.hprc-v1.0-pggb.gfa.gz",
        description="HPRC chr17 PGGB pangenome graph (GFA format)",
        size_mb=500,
    ),
}

# 1000 Genomes Project structural variants
# NA19240 is a Yoruba (YRI) individual from Ibadan, Nigeria
THOUSAND_GENOMES_DATA = {
    "sv_vcf": DataFile(
        name="ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
        url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz",
        description="1000 Genomes Phase 3 structural variants",
        size_mb=50,
    ),
}

# Smaller test datasets from GFA examples
TEST_DATA = {
    # Small GFA test file from vg project
    "vg_tiny": DataFile(
        name="tiny.gfa",
        url="https://raw.githubusercontent.com/vgteam/vg/refs/heads/master/test/tiny/tiny.gfa",
        description="Tiny test GFA from vg project",
        size_mb=0.001,
    ),
}


# ============================================================================
# Download Functions
# ============================================================================

def download_file(url: str, output_path: Path, show_progress: bool = True) -> bool:
    """Download a file with progress bar."""
    print(f"Downloading: {output_path.name}")
    print(f"  From: {url}")

    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))

        with open(output_path, 'wb') as f:
            if show_progress and total_size > 0:
                with tqdm(total=total_size, unit='iB', unit_scale=True, unit_divisor=1024) as pbar:
                    for chunk in response.iter_content(chunk_size=8192):
                        size = f.write(chunk)
                        pbar.update(size)
            else:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

        print(f"  ✓ Downloaded: {output_path} ({output_path.stat().st_size / 1024 / 1024:.1f} MB)")
        return True

    except requests.exceptions.RequestException as e:
        print(f"  ✗ Failed: {e}")
        if output_path.exists():
            output_path.unlink()
        return False


def decompress_gzip(input_path: Path, output_path: Optional[Path] = None) -> Optional[Path]:
    """Decompress a gzip file."""
    if output_path is None:
        output_path = input_path.with_suffix('')

    print(f"Decompressing: {input_path.name} -> {output_path.name}")

    try:
        with gzip.open(input_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"  ✓ Decompressed: {output_path}")
        return output_path
    except Exception as e:
        print(f"  ✗ Failed to decompress: {e}")
        return None


def check_vg_available() -> bool:
    """Check if vg tools are available for format conversion."""
    return shutil.which('vg') is not None


def convert_vg_to_gfa(vg_path: Path, gfa_path: Path) -> bool:
    """Convert VG format to GFA using vg tools."""
    if not check_vg_available():
        print("  ⚠ vg tools not found - cannot convert VG to GFA")
        print("    Install from: https://github.com/vgteam/vg")
        return False

    print(f"Converting: {vg_path.name} -> {gfa_path.name}")
    try:
        result = subprocess.run(
            ['vg', 'view', '-g', str(vg_path)],
            capture_output=True,
            check=True
        )
        with open(gfa_path, 'wb') as f:
            f.write(result.stdout)
        print(f"  ✓ Converted: {gfa_path}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Conversion failed: {e}")
        return False


# ============================================================================
# Main Download Functions
# ============================================================================

def download_test_data(output_dir: Path) -> List[Path]:
    """Download small test datasets for quick testing."""
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    print("\n=== Downloading Test Data ===\n")

    for key, data in TEST_DATA.items():
        output_path = output_dir / data.name
        if download_file(data.url, output_path):
            downloaded.append(output_path)

    return downloaded


def download_chr17_data(output_dir: Path, prefer_gfa: bool = True) -> Optional[Path]:
    """
    Download chromosome 17 pangenome data.

    This is the most useful data for NA19240 chr17 SV visualization.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n=== Downloading HPRC chr17 Pangenome Data ===\n")

    # Try PGGB GFA first (direct GFA format)
    if prefer_gfa:
        data = PGGB_DATA.get("chr17_pggb")
        if data:
            output_path = output_dir / data.name
            if download_file(data.url, output_path):
                # Decompress
                gfa_path = decompress_gzip(output_path)
                if gfa_path:
                    return gfa_path

    # Fall back to VG format (requires conversion)
    data = HPRC_DATA.get("chr17_vg")
    if data:
        output_path = output_dir / data.name
        if download_file(data.url, output_path):
            if check_vg_available():
                gfa_path = output_path.with_suffix('.gfa')
                if convert_vg_to_gfa(output_path, gfa_path):
                    return gfa_path
            else:
                print("\n  Note: Downloaded VG file. To convert to GFA:")
                print(f"    vg view -g {output_path} > {output_path.with_suffix('.gfa')}")
                return output_path

    return None


def download_sv_vcf(output_dir: Path) -> Optional[Path]:
    """Download 1000 Genomes structural variant VCF."""
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n=== Downloading 1000 Genomes SV Data ===\n")

    data = THOUSAND_GENOMES_DATA.get("sv_vcf")
    if data:
        output_path = output_dir / data.name
        if download_file(data.url, output_path):
            return output_path

    return None


def create_na19240_sample_from_vcf(vcf_path: Path, output_gfa: Path,
                                    chrom: str = "chr17",
                                    start: int = 10984564,
                                    end: int = 10993960) -> bool:
    """
    Extract NA19240-specific variants from VCF and create a sample GFA.

    Note: This creates a simplified representation. For full pangenome
    visualization, use the HPRC pangenome GFA files.
    """
    print(f"\nExtracting NA19240 variants for {chrom}:{start}-{end}")

    # This would require bcftools or similar - simplified version here
    print("  Note: Full VCF processing requires bcftools")
    print("  Creating a sample file with region metadata instead")

    # Create a minimal GFA with region info
    gfa_content = f"""H\tVN:Z:1.0
# NA19240 region: {chrom}:{start}-{end}
# Source: 1000 Genomes Project
# Sample: NA19240 (Yoruba, YRI population)
S\tref_segment\t{'N' * min(1000, end - start)}\tLN:i:{end - start}
"""

    with open(output_gfa, 'w') as f:
        f.write(gfa_content)

    print(f"  ✓ Created sample GFA: {output_gfa}")
    return True


def show_data_info():
    """Display information about available data sources."""
    print("""
╔═══════════════════════════════════════════════════════════════════════════╗
║                    NA19240 Pangenome Data Sources                         ║
╠═══════════════════════════════════════════════════════════════════════════╣
║                                                                           ║
║  NA19240 is a Yoruba (YRI) individual from Ibadan, Nigeria, included in   ║
║  both the 1000 Genomes Project and the Human Pangenome Reference          ║
║  Consortium (HPRC) Year 1 release.                                        ║
║                                                                           ║
╠═══════════════════════════════════════════════════════════════════════════╣
║  RECOMMENDED DATA:                                                        ║
╠═══════════════════════════════════════════════════════════════════════════╣
║                                                                           ║
║  1. HPRC PGGB chr17 (~500MB compressed)                                   ║
║     - Best for: Graph genome visualization                                ║
║     - Format: GFA (ready to use)                                          ║
║     - Contains: All 47 HPRC samples including NA19240                     ║
║                                                                           ║
║  2. HPRC Minigraph-Cactus chr17 (~800MB)                                  ║
║     - Best for: Higher-resolution graph                                   ║
║     - Format: VG (requires vg tools to convert to GFA)                    ║
║                                                                           ║
║  3. 1000 Genomes SV VCF (~50MB)                                           ║
║     - Best for: Variant annotation                                        ║
║     - Contains: Structural variants for NA19240                           ║
║                                                                           ║
╠═══════════════════════════════════════════════════════════════════════════╣
║  QUICK START:                                                             ║
║    python download_na19240_data.py --output ../data/ --test               ║
║                                                                           ║
║  FULL chr17 DATA:                                                         ║
║    python download_na19240_data.py --output ../data/ --chr17              ║
╚═══════════════════════════════════════════════════════════════════════════╝
""")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Download real NA19240 pangenome data from public sources",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show information about available data
  python download_na19240_data.py --info

  # Download small test files only (quick)
  python download_na19240_data.py --output ../data/ --test

  # Download chr17 pangenome data (~500MB-1GB)
  python download_na19240_data.py --output ../data/ --chr17

  # Download all available data
  python download_na19240_data.py --output ../data/ --all

Target Region of Interest:
  chr17:10,984,564-10,993,960 (hg38)
  This region contains a known structural variant in NA19240.
"""
    )

    parser.add_argument("--output", "-o", type=Path, default=Path("../data"),
                        help="Output directory for downloaded files (default: ../data)")
    parser.add_argument("--info", action="store_true",
                        help="Show information about available data sources")
    parser.add_argument("--test", action="store_true",
                        help="Download small test files only (quick)")
    parser.add_argument("--chr17", action="store_true",
                        help="Download HPRC chr17 pangenome data")
    parser.add_argument("--sv-vcf", action="store_true",
                        help="Download 1000 Genomes SV VCF")
    parser.add_argument("--all", action="store_true",
                        help="Download all available data")

    args = parser.parse_args()

    if args.info:
        show_data_info()
        return 0

    if not any([args.test, args.chr17, args.sv_vcf, args.all]):
        parser.print_help()
        print("\n⚠  No data selected. Use --test for quick start or --info for details.")
        return 1

    output_dir = args.output.resolve()
    print(f"\nOutput directory: {output_dir}")

    downloaded_files = []

    if args.test or args.all:
        files = download_test_data(output_dir)
        downloaded_files.extend(files)

    if args.chr17 or args.all:
        gfa_file = download_chr17_data(output_dir)
        if gfa_file:
            downloaded_files.append(gfa_file)

    if args.sv_vcf or args.all:
        vcf_file = download_sv_vcf(output_dir)
        if vcf_file:
            downloaded_files.append(vcf_file)

    # Summary
    print("\n" + "=" * 60)
    print("Download Summary")
    print("=" * 60)

    if downloaded_files:
        print(f"\n✓ Downloaded {len(downloaded_files)} file(s):\n")
        for f in downloaded_files:
            size_mb = f.stat().st_size / 1024 / 1024
            print(f"  • {f.name} ({size_mb:.1f} MB)")

        print(f"\nFiles are in: {output_dir}")
        print("\nTo load in graph-genome-viewer:")
        print("  cargo run --release")
        print("  File → Open GFA → select a .gfa file")
    else:
        print("\n✗ No files were downloaded successfully")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())




