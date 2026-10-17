#!/usr/bin/env bash
# Phase 2 (fragments). Requires phase 1 outputs (run_phase1.sh) for structures_main.parquet
# and subpocket_map.csv. The ligand pull is rate-limited and takes ~1 h cold.
set -euo pipefail
cd "$(dirname "$0")"
PY="conda run --no-capture-output -n pmml-klifs python -W ignore"
$PY scripts/08_pull_ligands.py       # data pull: ligand mol2 for the main set (rate-limited)
$PY scripts/09_fragment_assign.py    # processing: BRICS fragments, region assignment, categories
