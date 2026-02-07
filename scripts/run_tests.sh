#!/bin/bash
# Test script for graph-genome-viewer
set -e

echo "=== Running Rust Tests ==="
cd "$(dirname "$0")/.."

echo "1. Checking compilation..."
cargo check --lib 2>&1

echo "2. Running unit tests..."
cargo test --lib 2>&1

echo "3. Running integration tests..."
cargo test --test gfa_tests 2>&1

echo "=== Running Python Scripts ==="
cd scripts

echo "4. Testing Python script..."
python3 process_pangenome.py --help

echo "5. Generating sample data..."
python3 process_pangenome.py generate-sample --output ../data/test_sample.gfa

echo "6. Validating sample data..."
python3 process_pangenome.py validate --input ../data/test_sample.gfa

echo "=== All tests passed! ==="

