# pmml-generative-chem

conda env: pmml-klifs

Build a data pipeline and Bayesian analysis for a course assignment on
Box's loop (Probabilistic Models & ML).

RESEARCH QUESTION
Do different subpockets of the kinase ATP binding site have distinct
chemical preferences? This is a pilot for a larger project on whether
pocket-conditioned generative models generalize compositionally.

DATA
KLIFS (https://klifs.net), which provides a sequence-aligned 85-residue
kinase binding site with consistent numbering across all kinases.
Use the `opencadd.databases.klifs` Python package rather than raw API
calls. Relevant endpoints: kinase_names, structures_list,
interactions_get_IFP (85 residues x 7 interaction types = 595 bits),
interactions_match_residues, structure_get_pocket, structure_get_ligand.

BUILD IN TWO PHASES. Phase 1 must be complete and working before
starting phase 2.

Phase 1 - IFP-based analysis:
1. Pull all human kinase structures with their IFPs.
2. Group the 85 KLIFS positions into structural regions (hinge,
   G-loop, alphaC, DFG, front pocket, back pocket). Check whether KLIFS
   exposes these labels; if not, build a lookup table from the KLIFS
   numbering scheme and cite the source.
3. For each region, count interaction types observed across structures.
4. Model: per region, interaction counts ~ Categorical(pi_region),
   with a Dirichlet(alpha) prior. Conjugate, so posteriors are closed
   form: Dirichlet(alpha_v + n_v).
5. Report posterior distributions per region AND credible intervals on
   pairwise DIFFERENCES between regions. No frequentist tests, no
   p-values.
6. Run with alpha = 1 (uniform) and alpha = 0.5 (Jeffreys); report
   whether conclusions change.

Phase 2 - fragment-based analysis (only after phase 1 works):
1. Fetch ligand structures, fragment with RDKit BRICS
   (rdkit.Chem.BRICS.BRICSDecompose).
2. Assign each fragment to a region by minimum atom-to-residue
   distance. LOG ambiguous assignments and near-ties explicitly --
   these matter for the writeup.
3. Bin fragments into 4-5 chemical categories (keep the number small;
   counts per region will be modest and the Dirichlet prior dominates
   when counts are low).
4. Same Dirichlet-Categorical treatment as phase 1.

REQUIREMENTS
- Cache all API responses to disk. Do not re-fetch on every run.
- Separate scripts: data pull, processing, modeling, figures.
- Save intermediate data as parquet or CSV so each stage runs alone.
- One primary figure: posterior distributions by region, suitable for
  a 2-page writeup.
- Record counts per region/category in a table. Flag any cell with
  fewer than 10 observations as prior-dominated.
- Log every data-exclusion decision with its reason; these go in the
  writeup's criticism section.

DO NOT
- Do not use scipy.stats hypothesis tests. This is a Bayesian analysis.
- Do not silently drop structures. Log exclusions.
- Do not expand the category count to improve apparent resolution.

Start by exploring the KLIFS API and confirming what region annotations
are available, then show me the data pull plan before writing the full
pipeline.
---

## Pipeline (phase 1)

```
conda env create -f environment.yml   # or reuse the existing pmml-klifs env
./run_phase1.sh
```

| Stage | Script | Outputs |
|---|---|---|
| data pull | `01_pull_metadata.py` | `data/raw/{kinases,structures,ligands,ifps}_*.parquet` |
| processing | `02_select_structures.py` | `data/processed/structures_{full,main}.parquet`, `exclusions_{select,dedup}.csv` |
| data pull | `03_pull_pockets.py` | `data/raw/pockets/<klifs_id>.parquet` (rate-limited) |
| processing | `04_assign_subpockets.py` | `subpocket_map.csv`, `residue_assignments.parquet`, `near_ties.csv`, `results/subpocket_map.md` |
| processing | `05_count.py` | `counts_long.csv`, `results/counts_table.csv`, `results/kinase_counts.csv` |
| modelling | `06_model.py` | `results/posterior_summary.csv`, `pairwise_diffs.csv`, `pairwise_tv.csv`, `prior_sensitivity*.csv` |
| figures | `07_figures.py` | `figures/fig1_posteriors.{pdf,png}` |

Regions are KinFragLib's subpockets (AP, FP, SE, GA, B1, B2; Sydow et al., JCIM 2020) plus
`other` for positions >= 8 A from every subpocket centre. `other` is a heterogeneous mixture, so
it is a reference distribution, not a seventh subpocket. All exclusions are in
`data/processed/exclusions_*.csv`.
