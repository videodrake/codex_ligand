#!/bin/bash
# Phase 3 Diversity-Aware Docking — Auto-generated run script
# Generated: 2026-03-11T05:50:42.407913
# Total active jobs: 18
# Rounds: 3
# Max workers: 16

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Activate environment
# conda activate vina  # Uncomment if needed

OUTPUT_DIR="/home/user/codex_ligand/output/phase3_docking"

echo '=== Round 1: 6 jobs ==='

mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/173940"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/97806"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/VAX-C12_0"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/173940"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/97806"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/VAX-C12_0"

# Execute round 1 via Phase 3 runner
python -m egfr_pipeline.phase3.run_diverse_docking \
    --execute \
    --round 1 \
    --max_workers 16 \
    --output_dir "$OUTPUT_DIR"

# Evaluate saturation after round
python -m egfr_pipeline.phase3.run_diverse_docking \
    --evaluate-round 1 \
    --output_dir "$OUTPUT_DIR"

echo '=== Round 2: 6 jobs ==='

mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/173940"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/97806"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/VAX-C12_0"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/173940"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/97806"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/VAX-C12_0"

# Execute round 2 via Phase 3 runner
python -m egfr_pipeline.phase3.run_diverse_docking \
    --execute \
    --round 2 \
    --max_workers 16 \
    --output_dir "$OUTPUT_DIR"

# Evaluate saturation after round
python -m egfr_pipeline.phase3.run_diverse_docking \
    --evaluate-round 2 \
    --output_dir "$OUTPUT_DIR"

echo '=== Round 3: 6 jobs ==='

mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/173940"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/97806"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT01/VAX-C12_0"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/173940"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/97806"
mkdir -p "/home/user/codex_ligand/output/phase3_docking/runs/3GT8_raw/3GT8_raw_PKT02/VAX-C12_0"

# Execute round 3 via Phase 3 runner
python -m egfr_pipeline.phase3.run_diverse_docking \
    --execute \
    --round 3 \
    --max_workers 16 \
    --output_dir "$OUTPUT_DIR"

# Evaluate saturation after round
python -m egfr_pipeline.phase3.run_diverse_docking \
    --evaluate-round 3 \
    --output_dir "$OUTPUT_DIR"

echo '=== Phase 3 Docking Complete ==='
echo "Results in: $OUTPUT_DIR"