# Writeup plan: one anchor number per loop step

Rule: each step of the loop is carried by ONE number. Anything that doesn't anchor a step is
provenance: it lives in the repo, is cited by file at most, and gets no sentence of its own.
Step names are inferred (step 3 = "can the data answer the question"); rename to the course's.

| Step | Anchor number | Source | Status |
|---|---|---|---|
| 1. Question | 49 estimands: 7 regions (6 KinFragLib subpockets + `other` reference) x 7 IFP interaction types | design | fixed |
| 2. Data | 3,635 of 13,325 human KLIFS structures (27%) survive the scope filters + (kinase, ligand) dedup | `exclusions_*.csv` | TBD: minus assignment drops |
| 3. Can the data answer it? | share of per-structure residue assignments within 1 A of the runner-up centre | `near_ties.csv` | TBD (~25% in test run) |
| 4. Model | prior share of posterior mass, 7 alpha / (7 alpha + N_r), for the least-populated subpocket | `counts_table.csv` | TBD |
| 5. Inference / result | smallest posterior TV distance between two well-populated subpockets, with 95% CrI | `pairwise_tv.csv` | TBD |
| 6. Criticise / revise | number of kinases holding 50% of main-set observations (imbalance that dedup doesn't fix -> hierarchical model) | `kinase_counts.csv` | TBD |

## What each anchor carries

1. Defines the estimand, and states that `other` is a reference mixture, not a seventh subpocket.
2. One funnel number. Rules are named in a single sentence; per-rule counts are provenance.
3. The step-3 limitation: the partition is not well defined at residue level, so phase 1 tests
   "residues nearest different cavity centres" (the weaker question). Phase 2 asks the real one.
4. Says where alpha matters (small N_r) and where it can't (everywhere else), in one number,
   replacing a vacuous global sensitivity table.
5. If even the closest pair of well-populated subpockets has a CrI clear of zero, all are
   distinct; if not, that pair is the result. TV is biased upward at small N, which is
   negligible for well-populated regions and is one reason they're restricted to those.
6. Motivates the next model: kinase-level profiles partially pooled toward a global subpocket
   profile.

## Figure

Fig. 1 (two panels) serves steps 5 and 1: (a) the posterior TV matrix between region profiles;
(b) interaction-type profiles with 95% CrI for the four subpockets with the most main-set
observations (rule fixed before the full numbers), against the `other` reference band.

## Provenance only (repo, not prose)

- Exclusion counts per rule; the 27 kept-with-flag structures
- 85-position lookup table, modal shares, unstable positions (`results/subpocket_map.md`)
- Closest-atom vs centroid choice (10% vs 20% of bits in `other`): one clause in step 3, no number
- Per-cell prior sensitivity, including the ~2x zero-count ratio and the B1 (alpha+10)/(7alpha+10) example
- Pairwise per-type differences (`pairwise_diffs.csv`)
- Sensitivity analyses: per-structure labels, full (no dedup) set, top-3 families, alpha = 0.5
  (a single clause in step 6 says whether any changes the step-5 anchor)
- DFG-in filter and B1 sparsity, covalent gap, substrate rule (bullets in criticism notes)
