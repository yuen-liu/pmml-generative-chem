# Writeup plan: one anchor number per step, under Box's loop

Structure follows the course text's three activities (model design, posterior inference, model
criticism) with revision closing the criticism section. Six anchors; anything that does not
anchor one is provenance (repo only, cited by file at most). Remap if the course's Figure 1.1
has more boxes.

| Activity | # | Anchor number | Source | Status |
|---|---|---|---|---|
| Model design | 1 | 49 estimands: 7 regions (6 KinFragLib subpockets + `other` reference) x 7 IFP interaction types | design | fixed |
| Model design | 2 | 3,635 of 13,325 human KLIFS structures (27%) survive the scope filters + (kinase, ligand) dedup | `exclusions_*.csv` | TBD: minus assignment drops |
| Posterior inference | 3 | smallest posterior TV distance between two well-populated subpockets, with 95% CrI | `pairwise_tv.csv` | TBD |
| Model criticism | 4 | share of per-structure residue assignments within 1 A of the runner-up centre | `near_ties.csv` | TBD (~25% in test run) |
| Model criticism | 5 | the n >= 10 cell whose posterior mean moves most between alpha = 1 and 0.5 | `prior_sensitivity_by_cell.csv` | TBD (see definition) |
| Model criticism | 6 | posterior predictive overdispersion ratio R across kinases, pooled over the six subpockets (main, alpha = 1), with 95% interval -> revision | `ppc_overdispersion.csv` | TBD |

## Anchor 5: definition fixed before the full run

The test-run example (B1, N_r = 10, all apolar, posterior mean 0.647 -> 0.778) came from 142
AGC-only structures and will probably not survive the full run. Defined now:
- Anchor 5 = the `main`-analysis cell with n >= 10 and the largest |posterior-mean shift| between priors.
- If no main cell shifts by more than 0.01, take it from the family subsets (smaller N_r), and
  say it is a subset.
- If no cell shifts anywhere, anchor 5 becomes "the count threshold happens to suffice here", and
  the paragraph shrinks to the prior-share statistic plus the zero-count limit.

Paragraph order: (a) a cell that passes the n >= 10 threshold yet moves with the prior;
(b) why: with V = 7 the prior adds 7 alpha of mass however the data concentrate,
posterior mean (alpha + n)/(7 alpha + N_r) when all mass is in one type; (c) the statistic that
captures it, prior share 7 alpha/(7 alpha + N_r), with the n < 10 flag secondary; (d) the
zero-count cell as the limit (uniform / Jeffreys ratio -> 2; 1.98 observed); (e) why the
main result holds: anchor 3 rests only on well-populated regions, whose prior share is
negligible. Do NOT claim every reported cell clears the threshold: the zero-count charge
cells appear in the figure and tables.

## What the other anchors carry

1. Defines the estimand, and states that `other` is a reference mixture, not a seventh subpocket.
2. One funnel number. Scope filters are design assumptions; rules are named in one sentence.
3. If even the closest pair of well-populated subpockets has a CrI clear of zero, all are
   distinct; if not, that pair is the result. TV is biased upward at small N, which is
   negligible for well-populated regions.
4. The partition is not well defined at residue level, so phase 1 answers the weaker question
   ("residues nearest different cavity centres"). Phase 2 asks the real one.
6. The model check. R = observed / replicated multinomial deviance of per-kinase counts,
   per posterior draw (realized discrepancy; deviance, not Pearson, because of near-zero
   cells). Calibrated on simulated data: 0.95 [0.88, 1.04] when the model holds; 3.7 and 18
   under mild and strong kinase-level variation. Reading: the step-3 intervals are too narrow by
   roughly sqrt(R) (quasi-likelihood approximation; say so). That is the cost of conjugate
   inference showing up in criticism, and it motivates the revision: kinase-level profiles
   partially pooled toward a global subpocket profile. Kinases-for-50% of observations
   (`kinase_counts.csv`) is supporting detail here. No tail probability is reported.

## Figure

Fig. 1 (two panels) serves anchors 1 and 3: (a) the posterior TV matrix between region
profiles; (b) interaction-type profiles with 95% CrI for the four subpockets with the most
main-set observations (rule fixed before the full numbers), against the `other` reference band.

## Provenance only (repo, not prose)

- Exclusion counts per rule; the 27 kept-with-flag structures
- 85-position lookup table, modal shares, unstable positions (`results/subpocket_map.md`)
- Closest-atom vs centroid choice (10% vs 20% of bits in `other`): one clause under anchor 4
- Pairwise per-type differences (`pairwise_diffs.csv`)
- Sensitivity analyses: per-structure labels, full (no dedup) set, top-3 families, alpha = 0.5
  (one clause under anchor 6 says whether any changes anchor 3)
- DFG-in filter and B1 sparsity, covalent gap, substrate rule (bullets in criticism notes)
