#!/usr/bin/env bash
# Phase 1 end to end. Each stage reads the previous stage's files, so any stage can be rerun alone.
# API responses are cached under data/raw/; the pocket pull is rate-limited and takes ~1.5 h cold.
set -euo pipefail
cd "$(dirname "$0")"
PY="conda run --no-capture-output -n pmml-klifs python -W ignore"
$PY scripts/01_pull_metadata.py        # data pull: kinases, structures, ligands, IFPs
$PY scripts/02_select_structures.py    # processing: filters + (kinase, ligand) dedup, exclusions logged
$PY scripts/03_pull_pockets.py         # data pull: pocket coordinates for the main set (rate-limited)
$PY scripts/04_assign_subpockets.py    # processing: KinFragLib subpocket lookup table, near-ties
$PY scripts/05_count.py                # processing: region x interaction counts, kinase table
$PY scripts/06_model.py                # modelling: Dirichlet posteriors, pairwise diffs, prior sensitivity
$PY scripts/07_figures.py              # figures: fig1_posteriors
