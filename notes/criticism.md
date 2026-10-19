# Criticism-section notes (phase 1)

Draft points for the writeup. Numbers marked TBD are filled from the full run; numbers from
the partial (AGC-only) test run are not to be quoted.

## 1. Can the data answer the question? Residue-level partition is not well defined

KinFragLib's subpocket centres were built to partition *ligand* space: fragments sit inside the
cavity, so nearest-centre assignment suits them. Residues line the cavity from outside, and one
residue can face two subpockets. That is why the anchor sets overlap (51 anchors AP, SE and FP;
81 anchors FP, GA and B1). Forcing a partition onto residues imposes structure the geometry
lacks.

Evidence: share of per-structure assignments within 1 A of the runner-up centre (TBD; ~25% in
the test run); positions whose modal label holds in < 75% of structures (TBD; 16 in the test run).
This instability says the partition is not well defined at residue level; it is not noise.

Consequence for the claim: phase 1 does not test "do subpockets have distinct chemical
preferences". It tests the adjacent, weaker question "do residues nearest different cavity
centres have different IFP interaction-type profiles". Phase 2 asks the question properly,
because fragments are what the centres were built for.

Also: closest-heavy-atom distance (not centroid) was chosen after seeing the share of IFP bits
in `other` (10% vs 20%), before looking at any subpocket profile. State this as a post-hoc
design choice.

## 2. Prior sensitivity is cell-conditional, not global

With order 10^3 observations per region, alpha = 1 vs 0.5 changes nothing where counts are
large. Report it conditional on counts (`results/prior_sensitivity_by_count.csv`,
`prior_sensitivity_by_cell.csv`):
- n_rv = 0 cells: posterior mean = alpha / (7 alpha + N_r), i.e. pure prior; uniform / Jeffreys
  -> 2 as N_r grows. Example: AP x protein-cation/ligand-anion has 0 observations even in the
  full set (TBD final); its posterior is the prior's choice.
- PRIMARY diagnostic: the prior's share of posterior mass, 7 alpha / (7 alpha + N_r). The prior
  contributes 7 alpha regardless of how the data concentrate, so V matters as much as n. The
  n < 10 cell flag is secondary. A cell can clear n >= 10 and still move if N_r is small:
  with all N_r = 10 observations in one type the posterior mean is (alpha + 10)/(7 alpha + 10),
  0.647 at alpha = 1 and 0.778 at alpha = 0.5 (seen for B1 in the test run).

## 3. The DFG-in filter partly causes back-pocket sparsity

The KinFragLib scope filter keeps DFG-in only (1,821 non-DFG-in human structures dropped).
Type-II inhibitors bind DFG-out and reach deepest into the back pocket, so B1 sparsity is partly
an artefact of the structural filter, not a fact about kinase chemistry. In the test run B1 holds
0.3% (main) / 3.0% (full) of observations; B2 holds 8-11% and is not sparse (TBD final).
B1/B2 stay unmerged in phase 1. Revisit in phase 2, where fragment counts are what matters; if
merged there, report the unmerged result alongside.

Do NOT quantify the filter by applying the lookup table to DFG-out structures: DFG-out flips the
motif and moves the back-pocket geometry (83 is xDFG; 81 anchors FP, GA and B1), so a
DFG-in-derived table there measures its own misapplication, not the filter's effect. Doing it
properly means recomputing centres per DFG-out structure (out of scope for this writeup).

Do not compare these shares to KinFragLib's B1/B2 fragment shares: theirs are shares of
fragments, ours of IFP interactions at residue positions. Both can hold at once, because residues
can line a subpocket that few fragment centroids land in (the point of section 1 again).

Main vs full B1 share: the test-run 0.3% vs 3.0% is NOT a dedup effect. Main then covered
only the cached AGC kinases. Compare only on the full run, and note that even then main vs full
mixes dedup with a shift in kinase mix (dedup removes most from the most-crystallised kinases).

## 4. Other exclusions and gaps

- KinFragLib's covalent-ligand filter not applied (it needs RCSB downloads).
- The substrate rule (phosphate/ribose SMARTS) also removes nucleoside-mimetic inhibitors.
- Kinase imbalance survives dedup (`results/kinase_counts.csv`; TBD: kinases needed for 50% of
  observations). The family-level posteriors check whether the pooled result is the
  most-crystallised families. The proposed next model is hierarchical: kinase-level profiles
  partially pooled toward a global subpocket profile.
- `other` is a mixture of heterogeneous positions: a reference distribution, not a seventh
  subpocket.

---

# Phase 2 (fragments) — full run, 2026-10-18

Scope: 10,141 BRICS fragments from 3,424 ligands, 214 kinases. Exclusions: 139 ligands with a
BRICS fragment > 22 heavy atoms (KinFragLib rule), 56 unparsable mol2, 14 missing anchor
residues, 2 un-fetched (KLIFS 2776, 2777).

## 1. The decomposition is partly a property of the metric

Centroid-to-subpocket-centre: 12.0% of fragments within 1 A of a runner-up.
Minimum atom-to-residue (the original spec's rule): 88.6% near-ties, and it agrees with the
centroid rule on only 54.1% of fragments. A fragment inside the cavity is at contact distance
from residues lining several regions at once, so that rule barely discriminates.
With phase 1's 25.4% residue-level near-ties, this says the subpocket decomposition is partly a
property of the metric, not only of the pocket — an argument FOR learning subpockets as latents.

## 2. Merge accounting

46.8% of BRICS fragments were tiny (< 3 heavy atoms), but only 8.0% were relabelled and those
held just 2.5% of heavy atoms. The merge moved little chemistry. Median final fragment: 9 heavy
atoms; 0.3% still below 3.

## 3. Result: 10 of 15 subpocket pairs distinguishable (vs 14/15 in phase 1)

Same rule as phase 1 (sqrt(R)-widened CrI lower bound vs the 97.5% point of the same-profile
null at N/R; R = 2.44 [2.31, 2.58], so intervals are ~1.56x too narrow). Identical verdict at
alpha = 0.5. Every AP and FP pair separates; the five failures all involve B1 (192 fragments) or
B2 (299), whose nulls are large because counts are small. Not a null result.

FP-SE flips between phases: the one failure in phase 1 (IFP), clearly distinguishable in phase 2
(TV 0.123, widened lower bound 0.084 vs null 0.060). The two subpockets' protein-side interaction
patterns are alike while their fragment chemistry differs, which is worth a sentence.

## 4. Category distribution and what drives the result

aromatic 68.4%, charged 14.2%, hbond 12.4%, aliphatic 5.0%. AP is 95.3% aromatic (the
hinge-binding heteroaromatics); FP and SE carry the charged fragments (0.226 / 0.252 vs 0.017 in
AP) — solubilising amines at the solvent front. Chemically recognisable, which is the point of
the revised priority ordering.

## 5. Overdispersion is lower than phase 1 but still real

R = 2.44 (fragments) vs 2.73 (IFP, deduplicated). Same conclusion: the independence assumption
fails and the hierarchical revision is motivated. 17 of 214 kinases hold 50% of fragments.

## 6. B1/B2 remain the sparse regions

B1: 192 fragments, prior share 2.0%; B2: 299. The DFG-in filter is still the likely cause
(type-II ligands reach the back pocket). Phase 2 has ~10x less data per region than phase 1,
so this bites harder here.
