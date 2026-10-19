# Writeup plan: eight anchors under Box's loop

Structure follows the course text's three activities: **in model design we bake our goals and
assumptions into the model; in posterior inference we solve the algorithmic problem of
calculating the posterior; in model criticism we ask how the model has succeeded and where it
has fallen short.** Revision closes criticism.

**Hard cap: eight anchors.** Two pages cannot report both phases symmetrically — that yields two
thin analyses instead of one argument. Anything that does not serve the flip (anchor 1) or the
criticism is provenance: repo only, cited by file at most.

**Posterior inference carries no numeric anchor, and that is the point.** Conjugacy trades
inference difficulty for modelling fidelity, and the bill comes due in criticism (anchors 6-8).

The phase-1-only plan, with the full decision record for the old anchors 3 and 5, is preserved at
`writeup_plan_phase1.md`.

| Activity | # | Anchor | Source |
|---|---|---|---|
| Criticism (lead) | 1 | **The flip.** FP-SE is the least-separated subpocket pair under both measurements, and only the ligand-side one resolves it | `pairwise_tv_null*.csv` |
| Model design | 2 | One model, fit twice: 7 regions x 7 IFP types = 49 estimands; 7 regions x 4 chemical categories = 28 | design |
| Model design | 3 | 3,619 of 13,325 human KLIFS structures (27.2%) -> 3,424 ligands -> 10,141 BRICS fragments | `exclusions_*.csv` |
| Criticism | 4 | **AP is 95.3% aromatic** [0.946, 0.960] - external validation | `posterior_summary_fragments.csv` |
| Criticism | 5 | Metric dependence: 12.0% (centroid) vs 88.6% (min-atom) near-ties, the two rules agreeing on 54.1% of fragments | `near_ties*.csv` |
| Criticism | 6 | Overdispersion R = 2.44 [2.31, 2.58] (fragments) and 2.73 [2.62, 2.85] (IFP) -> revision | `ppc_overdispersion*.csv` |
| Criticism | 7 | Sparse-region power: all five of the fifteen phase-2 pair failures involve B1 (192 fragments) or B2 (299) | `pairwise_tv_null_fragments.csv` |
| Criticism | 8 | Prior dominance: the n >= 10 cell whose posterior mean moves most between alpha = 1 and 0.5 | `prior_sensitivity_by_cell*.csv` |

## Anchor 1: the flip (the spine)

FP and SE are the hardest pair to separate under both measurements. Only the ligand-side
measurement separates them. This is the only finding neither phase produces alone.

**Numbers.**
- Phase 1 (IFP): TV 0.031 [0.025, 0.037], widened lower bound 0.021, null 97.5% 0.023 -> **not
  distinguishable**. Rank 1 of 15 (smallest TV; median 0.137).
- Phase 2 (fragments): TV 0.123 [0.098, 0.148], widened lower bound 0.084, null 97.5% 0.060 ->
  **distinguishable**. Rank 2 of 15 (median 0.221).
- Same verdicts at alpha = 0.5.

**What separates them (settled before the sentence was written).** No single category drives it;
two do, in opposition, carrying 75.5% of the TV between them:

| Category | FP | SE | diff (95% CrI) | share of TV |
|---|---|---|---|---|
| Aromatic | 0.459 | 0.556 | -0.097 [-0.126, -0.068] | 39.3% |
| Polar, non-aromatic | 0.235 | 0.146 | +0.089 [+0.066, +0.111] | 36.2% |
| Aliphatic | 0.080 | 0.046 | +0.034 [+0.020, +0.048] | 13.8% |
| Charged | 0.226 | 0.252 | -0.026 [-0.051, -0.001] | 10.7% |

State it as an aromatic-for-polar **substitution** at matched charged content: FP trades aromatic
rings for polar non-aromatic groups relative to SE. Do **not** attribute the difference to charge
- the charged shares are nearly equal (0.226 vs 0.252) and that CrI's upper bound is -0.0013.
Charged separates FP and SE *together* from AP (0.017), which is a different claim.

**The defence: the flip runs against the power gradient.** Phase 1 had 3.3x more observations on
this exact pair (FP 6,448 + SE 8,659 = 15,107) than phase 2 (2,640 + 1,981 = 4,621), and still
found no difference, while phase 2 found one with less. A power artefact would run the other way.
Use the per-pair ratio, not the global one (68,493 vs 10,141 observations); the per-pair figure is
the one the argument needs. Normalised against each phase's own null - the only scale-free
comparison available, since the two vocabularies put TV on different scales (phase 2's median
pairwise TV is 1.6x phase 1's) - TV/null is 1.32 in phase 1 against 2.05 in phase 2.

**Caution, one clause, mandatory.** The two phases use different observation units and different
nulls, so the dissociation is **between measurements, not a proven fact about the pocket**.

**Also disclose, one clause.** Phase 1's undeduplicated set *does* find FP-SE distinguishable
(0.025 widened vs 0.021 null). Its effect size is essentially unchanged (0.034 vs 0.031); it
passes only because the null shrinks with N. Pseudo-replication passing a test the main set
fails is consistent with the spine, but it is in `pairwise_tv_null.csv` and should be disclosed
rather than left to be found.

**What the headline is NOT.** Not "subpockets have distinct chemical preferences" - that is the
background result (14 of 15 pairs in phase 1, 10 of 15 in phase 2) and it goes in one sentence.

## Anchor 2: one model, fit twice

The Dirichlet-Categorical model is stated once and fit to two observation units: protein-side
interaction types (phase 1) and ligand-side fragment chemistry (phase 2). Describing it twice
costs almost nothing, which is the economy that makes anchor 1 affordable.

    observation in region r ~ Categorical(pi_r),  pi_r ~ Dirichlet(alpha * 1)
    =>  pi_r | data ~ Dirichlet(alpha + n_r)

Assumptions carry over unchanged: regions independent, observations exchangeable within a region
(violated across kinases - anchor 6), `other` a reference mixture rather than a seventh subpocket.

**Disclose the category-definition revision as a specification fix, not post-hoc tuning.** The
first ordering (charged -> hbond -> aromatic -> aliphatic) sent every N-heterocycle to `hbond`:
pyridine, pyrimidine, imidazole, aminopyrimidine - the core hinge-binding chemistry of the whole
dataset. Reordered to charged -> aromatic -> hbond -> aliphatic. The test applied: would the
change be made if the shares were even? Yes - it is a validity error, independent of outcome.
Category count was never expanded.

## Anchor 3: scope

One funnel sentence. 13,325 human KLIFS structures -> 3,619 after the scope filters (non-human,
atypical kinase, DFG not in, no orthosteric ligand, missing/empty IFP, substrate-like ligand) and
(kinase, ligand) deduplication -> 3,424 parsable ligands -> 10,141 BRICS fragments after the
small-fragment merge. Scope filters are design assumptions, named in one sentence. Every drop is
logged with its rule in `exclusions_*.csv`; nothing is silently dropped.

Phase 2 exclusions: 139 oversized BRICS fragments (>22 heavy atoms, KinFragLib's rule), 56
unparsable ligands, 14 missing anchor residues, 2 unavailable.

## Anchor 4: AP is 95.3% aromatic (external validation)

Posterior mean 0.953 [0.946, 0.960] over 3,322 fragments. The adenine pocket is where
hinge-binding heteroaromatics sit - this is what medicinal chemistry says should be there, and it
was not fit to. It is external evidence that the pipeline measures something real, which is worth
more than most numbers competing for this space. Contrast: AP's charged share is 0.017 against
0.226 (FP) and 0.252 (SE), the solubilising groups at the solvent front.

## Anchor 5: metric dependence

The subpocket partition is partly an artefact of the distance rule used to define it.

- Phase 2, centroid-to-centre: **12.0%** of fragment assignments within 1 A of the runner-up centre.
- Phase 2, the spec's minimum atom-to-residue rule: **88.6%**.
- The two rules agree on only **54.1%** of fragments.
- Phase 1, residue level: 25.4% near-ties (51,804 / 203,858), 21 unstable positions.

Stronger than phase 1's 25.4% on its own, and it is the argument for the larger project: if the
partition is this sensitive to an arbitrary geometric choice, subpockets are better treated as
latent variables to be learned than as a fixed labelling to condition on.

## Anchor 6: overdispersion -> revision

R = observed / replicated multinomial deviance of per-kinase counts, per posterior draw (realized
discrepancy, Gelman-Meng-Stern; deviance rather than Pearson because of near-zero cells).

- Fragments: **R = 2.44 [2.31, 2.58]** -> CrIs roughly 1.56x too narrow.
- IFP: **R = 2.73 [2.62, 2.85]** -> roughly 1.65x too narrow. Undeduplicated set 5.57.
- Calibration on simulated data: 0.95 [0.88, 1.04] when the model holds; 3.7 and 18 under mild
  and strong kinase-level variation.

Say explicitly that sqrt(R) as an interval inflation factor, and N/R as an effective sample size,
are **quasi-likelihood approximations, not derivations**. Both phases exceed the
model-holds range in the same direction and by a similar amount, so the misfit is a property of
the model rather than of either observation unit. Supporting detail: 17 of 227 kinases hold 50% of
IFP observations; 17 of 214 hold 50% of fragments.

**Revision (unchanged):** kinase-level profiles partially pooled toward a global subpocket
profile. This is where conjugacy's bill comes due - the revision has no closed form and needs MCMC
or variational inference, which is exactly the trade the inference section describes. No tail
probability is reported.

## Anchor 7: sparse-region power

Phase 2 distinguishes 10 of 15 subpocket pairs against 14 of 15 in phase 1. **All five failures
involve B1 (192 fragments) or B2 (299)**: GA-B1, SE-B2, FP-B2, GA-B2, B1-B2. Every pair involving
AP or FP separates. So the weaker phase-2 result is a power limit in two sparse regions, not an
absence of signal - and the sparsity is traceable to a design decision, the DFG-in filter, which
excludes the structures whose ligands reach the back pocket. Identical verdicts at alpha = 0.5.

This is the honest reading of the 10-vs-14 drop and it belongs here rather than being left to
look like a contradiction of anchor 1.

## Anchor 8: prior dominance

A cell can clear the n >= 10 flag and still move with the prior if its region's total N_r is
small. The statistic that captures this is the prior's share of posterior concentration,
V*alpha / (V*alpha + N_r); the n < 10 count flag is secondary.

Paragraph order: (a) a cell that passes the n >= 10 threshold yet moves with the prior; (b) why -
with V categories the prior adds V*alpha of mass however the data concentrate; (c) the prior-share
statistic; (d) the zero-count cell as the limit (uniform/Jeffreys posterior-mean ratio -> 2, 1.998
observed); (e) why anchor 1 holds regardless - FP and SE both have prior share < 1%.

Do **not** claim every reported cell clears the threshold: zero-count cells appear in the figures
and tables.

## Figures

**Fig. 1 - the flip. Two panels, one pair, two measurements.** Panel (a) phase 1 (IFP
interactions), panel (b) phase 2 (fragment chemistry). Each panel, on a TV axis:

1. The FP-SE posterior TV: point at the mean, solid bar for the 95% CrI, lighter extension to the
   sqrt(R)-widened interval, so the overdispersion correction is visible rather than asserted.
2. The null: simulated same-profile TV at N/R as a filled density, with its 97.5% point as a
   vertical rule spanning the panel. The verdict reads off as "does the widened bar clear the
   rule": 0.084 vs 0.060 clears; 0.021 vs 0.023 does not.
3. A context strip along the top: the other 14 pair TVs as grey ticks, so the figure is honest
   that FP-SE is rank 1 and rank 2 of 15.

Independent x-axes per panel (phase 1 ~[0, 0.25], phase 2 ~[0, 0.5]).
**Caption must state that the two panels use different vocabularies and that the axes are not
comparable** - the panels otherwise invite direct comparison. The same note belongs on the
secondary figures, whose panel (b) is log-scaled for IFP and linear for fragments.

Optional panel (c): the four per-category FP-SE differences as dots with CrIs - the anchor 1 table
plotted. Include if Fig. 1 gets a half page; otherwise it moves to the secondary figure.

**Figs. S1/S2 (secondary, demoted from primary):** the existing two-panel per-phase figures - TV
matrix across all 7 regions, and interaction/category profiles with CrIs for the four
best-populated subpockets against the `other` reference band. These serve anchors 2, 4 and 7 and
are no longer the lead.

## Provenance only (repo, not prose)

- Exclusion counts per rule; the 27 kept-with-flag structures
- 85-position lookup table, modal shares, 21 unstable positions (`results/subpocket_map.md`)
- Closest-atom vs centroid residue rule (10% vs 20% of bits in `other`): one clause under anchor 5
- Small-fragment merge accounting: 46.8% of BRICS fragments tiny, 8.0% relabelled, holding 2.5%
  of heavy atoms - the merge moved very little chemistry. One clause if space allows.
- Pairwise per-type differences for pairs other than FP-SE (`pairwise_diffs*.csv`)
- Sensitivity analyses: per-structure labels, full (no dedup) set, top-3 families, alpha = 0.5
- Family subsets: with 10-30x fewer observations they fail many pairs and corroborate nothing;
  state as separate observations, not support
- KinFragLib's lowest-median-logP check on FP/SE: a post-hoc qualitative plausibility note on
  similarity only, not direction. Now largely superseded by anchor 1's own decomposition.
- DFG-in filter and B1/B2 sparsity (anchor 7 carries the one-sentence version), covalent gap,
  substrate rule
