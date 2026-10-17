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
- The prior's total weight is 7 alpha / (7 alpha + N_r), so a cell can clear n >= 10 and still
  move if its region's total N_r is small. The test run showed this for B1 (N_r = 10, all apolar;
  posterior mean 0.65 vs 0.78). The n < 10 cell flag alone misses it; report the region-level
  prior weight too.

## 3. The DFG-in filter partly causes back-pocket sparsity

The KinFragLib scope filter keeps DFG-in only (1,821 non-DFG-in human structures dropped).
Type-II inhibitors bind DFG-out and reach deepest into the back pocket, so B1 sparsity is partly
an artefact of the structural filter, not a fact about kinase chemistry. In the test run B1 holds
0.3% (main) / 3.0% (full) of observations; B2 holds 8-11% and is not sparse (TBD final).
Collapsing B1+B2 would be a post-hoc merge driven by B1's counts; if done (likely in phase 2,
where fragment counts are smaller), say so and keep the unmerged result alongside.

## 4. Other exclusions and gaps

- KinFragLib's covalent-ligand filter not applied (it needs RCSB downloads).
- The substrate rule (phosphate/ribose SMARTS) also removes nucleoside-mimetic inhibitors.
- Kinase imbalance survives dedup (`results/kinase_counts.csv`; TBD: kinases needed for 50% of
  observations). The family-level posteriors check whether the pooled result is the
  most-crystallised families. The proposed next model is hierarchical: kinase-level profiles
  partially pooled toward a global subpocket profile.
- `other` is a mixture of heterogeneous positions: a reference distribution, not a seventh
  subpocket.
