# Writeup plan: one anchor number per step, under Box's loop

Structure follows the course text's three activities (model design, posterior inference, model
criticism) with revision closing the criticism section. Six anchors; anything that does not
anchor one is provenance (repo only, cited by file at most). Remap if the course's Figure 1.1
has more boxes.

| Activity | # | Anchor number | Source | Status |
|---|---|---|---|---|
| Model design | 1 | 49 estimands: 7 regions (6 KinFragLib subpockets + `other` reference) x 7 IFP interaction types | design | fixed |
| Model design | 2 | 3,635 of 13,325 human KLIFS structures (27%) survive the scope filters + (kinase, ligand) dedup | `exclusions_*.csv` | **3,619 of 13,325 (27.2%)** after 16 assignment drops |
| Posterior inference | 3 | pairwise subpocket distinguishability vs a same-profile null at N/R (all six subpockets: all have prior share < 1%) | `pairwise_tv_null.csv` | **14 of 15 subpocket pairs distinguishable; FP-SE is not: TV 0.031, widened lower bound 0.021 vs N/R null 0.023** (same verdict at alpha = 0.5) |
| Model criticism | 4 | share of per-structure residue assignments within 1 A of the runner-up centre | `near_ties.csv` | **25.4%** of subpocket-labelled assignments (51,804 / 203,858); 21 unstable positions |
| Model criticism | 5 | the n >= 10 cell whose posterior mean moves most between alpha = 1 and 0.5 | `prior_sensitivity_by_cell.csv` | main max shift 0.001 (< 0.01), so taken from a family subset: **JakA, B1 x apolar, n = 22, N_r = 24: 0.742 -> 0.818**, prior share 0.23 |
| Model criticism | 6 | posterior predictive overdispersion ratio R across kinases, pooled over the six subpockets (main, alpha = 1), with 95% interval -> revision | `ppc_overdispersion.csv` | **R = 2.73 [2.62, 2.85]**, CrIs ~1.65x too narrow; full (no dedup) set R = 5.57 |

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

## Full-run notes (2026-10-17)

- Kinases for 50% of main observations: 17 of 227 (CDK2 alone 9.8%). Supporting detail for anchor 6.
- Dedup DID work, but the evidence is the PPC, not B1's share: R = 5.57 (full) vs 2.73 (main).
  B1's share is 3.1% in both main and full, so there is no B1 drop.
- Zero-count cells in main: 6; uniform / Jeffreys mean ratio 1.997-1.9998 (the limit, 2).
- Anchor 3 issue: TV > 0 by construction, so "CrI clear of zero" tests nothing. Null baseline
  (two regions sharing one pooled profile, same N): FP-SE observed 0.031 [0.025, 0.037] vs null
  97.5% 0.013 (naive) / 0.023 (at N/R). AP-B2 0.073 vs 0.016 / 0.027.
- The figure rule's top four are AP, GA, SE, B2. FP (6,448) missed B2 (6,749) by 301
  observations, and FP-SE is the closest subpocket pair overall.

## Anchor 3 (decided 2026-10-17)

Population: all six subpockets. The prior-share rule decides it (every subpocket's prior share is
< 1%). The figure's "top four" is a display choice, not the anchor population.

Rule (in `06b_posterior_predictive.py`, fixed before pairwise results were inspected): a pair is
distinguishable if the lower end of its 95% CrI, widened about the mean by sqrt(R), lies above
the 97.5% point of the null TV for two regions sharing the pair's pooled profile, simulated at
N/R. Say in one clause that N/R as an effective sample size is an approximation (quasi-likelihood).

Result: 14 of 15 pairs distinguishable under both priors. FP-SE is not: TV 0.031 [0.025, 0.037],
widened lower bound 0.021 vs null 0.023 (alpha = 1) / 0.024 (alpha = 0.5); 1,000 null replicates.
The margin is narrow. At 200 replicates the Jeffreys run passed, which was Monte Carlo noise.
Undeduplicated set: FP-SE passes (0.025 vs 0.021); it has more observations, even after dividing by the larger R.

Families (state as separate observations, NOT corroboration): with 10-30x fewer observations,
FP-SE is not distinguishable in MAPK or JakA but is in CDK. Those subsets fail several other pairs
too (MAPK: FP-SE, FP-GA, FP-B1, SE-GA; JakA: 10 of 15), so they are underpowered and don't single
FP-SE out. Drop the earlier "within MAPK, SE-GA is similarly marginal" line: under the formal rule
MAPK fails FP-SE itself, plus three others.

External check (qualitative): in KinFragLib's own fragment statistics, FP and SE have the two
lowest median fragment logP of the six subpockets (v1.0.0: FP 0.73, SE 0.99; next lowest B1 1.16;
`notebooks/2_2_fragment_analysis_statistics.ipynb`, cell 126; current master: 0.78, 1.01).
Caveats to state: (i) a different quantity on different objects (fragment lipophilicity vs residue
IFP types), so it supports SIMILARITY of FP and SE only, not the direction of any difference;
(ii) medians with SD ~1.2-1.3, no intervals; (iii) looked up after the FP-SE result, so it is a
post-hoc plausibility check. Our profiles agree on similarity: apolar share 0.895 (FP) vs 0.893 (SE).
They are not identical: SE has protein-anion : ligand-cation contacts at 0.017 vs 0.001 in FP
(consistent with basic solubilising groups at the solvent front).

Headline: all subpocket pairs are distinguishable except FP-SE, which fails narrowly once
overdispersion is accounted for, and FP and SE are also the two subpockets whose fragments
KinFragLib found least lipophilic.
