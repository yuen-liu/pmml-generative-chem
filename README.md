# Do kinase ATP-site subpockets have distinct chemical preferences?

A Dirichlet–Categorical analysis of the kinase ATP binding site, built around Box's loop for a
Probabilistic Models & ML assignment. Pilot for a larger project on whether pocket-conditioned
generative models generalise compositionally.

The same model is fit twice, to opposite sides of the interface:

- **Phase 1 (protein side)** — KLIFS interaction-fingerprint bits, typed into 7 interaction types.
- **Phase 2 (ligand side)** — BRICS fragments of co-crystallised ligands, typed into 4 chemical
  categories.

**Result.** Subpocket profiles are distinguishable — 14 of 15 pairs on interaction patterns and 10
of 15 on fragment chemistry, with every failure in a sparse region — but the decomposition is
partly an artefact of the assignment metric, and the credible intervals are overconfident by
roughly 1.6×.

The writeup is `notes/writeup.md`; the primary figure is `figures/fig1_main.pdf`.

## Running it

```sh
conda env create -f environment.yml    # creates the pmml-klifs env
conda activate pmml-klifs
./run_phase1.sh
./run_phase2.sh                        # phase 1 must have run first
```

Every KLIFS response is cached under `data/raw/` and never re-fetched; `data/` is
gitignored, so a fresh clone re-pulls (several hours, rate-limited). Each stage reads its
predecessor's files from disk and can be run alone. The coordinate pulls (stages 3 and 8) are
throttled to one request every 1.5 s by default (`--interval`), write atomically, and verify
cached files on startup, so they are safe to interrupt and resume. klifs.net is a public academic
server; please leave the throttle on.

### Phase 1 — interaction fingerprints

| Stage | Script | Outputs |
|---|---|---|
| data pull | `01_pull_metadata.py` | `data/raw/{kinases,structures,ligands,ifps}_*.parquet` |
| processing | `02_select_structures.py` | `structures_{full,main}.parquet`, `exclusions_{select,dedup}.csv` |
| data pull | `03_pull_pockets.py` | `data/raw/pockets/<klifs_id>.parquet` (rate-limited) |
| processing | `04_assign_subpockets.py` | `subpocket_map.csv`, `residue_assignments.parquet`, `near_ties.csv`, `results/subpocket_map.md` |
| processing | `05_count.py` | `counts_long.csv`, `results/counts_table.csv`, `results/kinase_counts.csv` |
| modelling | `06_model.py` | `results/posterior_summary.csv`, `pairwise_diffs.csv`, `pairwise_tv.csv`, `prior_sensitivity*.csv` |
| criticism | `06b_posterior_predictive.py` | `results/ppc_overdispersion.csv`, `pairwise_tv_null.csv` |
| figures | `07_figures.py` | `figures/fig1_posteriors.{pdf,png}`, `figS1_posteriors_by_type.{pdf,png}` |

### Phase 2 — ligand fragments

| Stage | Script | Outputs |
|---|---|---|
| data pull | `08_pull_ligands.py` | `data/raw/ligands/<klifs_id>.mol2` (rate-limited) |
| processing | `09_fragment_assign.py` | `fragments.parquet`, `exclusions_fragment.csv`, `fragment_flags.csv` (near-ties, both assignment rules), `merge_accounting.csv` |
| processing | `10_count_fragments.py` | `counts_fragments_long.csv`, `results/fragment_counts_table.csv` |
| modelling | `06_model.py --phase fragment` | `results/*_fragments.csv` |
| criticism | `06b_posterior_predictive.py --phase fragment` | `results/ppc_overdispersion_fragments.csv`, `pairwise_tv_null_fragments.csv` |
| figures | `07_figures.py --phase fragment` | `figures/fig1_posteriors_fragments.{pdf,png}` |

### Cross-phase

| Stage | Script | Outputs |
|---|---|---|
| figures | `11_fig_main.py` | `figures/fig1_main.{pdf,png}` — the primary figure |

`11_fig_main.py` reads both phases' results, so it runs at the end of `run_phase2.sh`.

Stages 6, 6b and 7 are shared: `--phase {ifp,fragment}` swaps the observation unit and category
vocabulary and nothing else.

## Regions

Regions are KinFragLib's six subpockets — AP (adenine), FP (front), SE (solvent-exposed), GA (gate
area), B1 and B2 (back) — each defined as the mean Cα position of a set of anchor residues in
KLIFS's 85-position alignment, recomputed per structure (Sydow, Schmiel, Mortier & Volkamer,
*JCIM* 2020, 60, 6081–6094; anchor sets taken from their `fragmentation/cli.py`, not from memory).
A seventh region, `other`, holds everything ≥ 8 Å from every centre. **`other` is a heterogeneous
mixture and therefore a reference distribution, not a seventh subpocket.**

## Model

For each region *r*, independently:

```
observation ~ Categorical(pi_r),   pi_r ~ Dirichlet(alpha * 1)
=>  pi_r | data ~ Dirichlet(alpha + n_r)        (conjugate, closed form)
```

Monte Carlo is used only for functionals with no closed form: pairwise differences, total-variation
distance between region profiles, and the simulated nulls. Everything is run at alpha = 1 and
alpha = 0.5; no reported conclusion differs between them. No frequentist tests anywhere.

A pair of regions is called **resolved** when the lower end of its 95% CrI on TV distance, widened
about the mean by √R, lies above the 97.5% point of the TV distribution simulated for two regions
sharing the pair's pooled profile at effective sample size N/R. The rule was fixed before any
pairwise result was inspected. TV is positive by construction, so "the CrI excludes zero" would
test nothing — hence the explicit null.

## Known limitations

These are the writeup's criticism section; the numbers are in `notes/criticism.md`.

- **Metric dependence.** 12.0% of fragments lie within 1 Å of a runner-up subpocket centre under
  the centroid rule, and 88.6% under the minimum atom-to-residue rule; the two rules disagree on
  46% of fragments. At residue level, 25.4% near-ties and 21 of 85 positions with no stable label.
  The partition is imposed by a distance rule, not by the data.
- **Overdispersion.** Posterior predictive check across kinases gives R = 2.73 (interactions) and
  2.44 (fragments) against 0.95 when the model holds, so exchangeability across kinases fails and
  the intervals are ~1.6× too narrow. √R and N/R are quasi-likelihood approximations, not
  derivations.
- **Sparse back pockets.** All five unresolved phase-2 pairs involve B1 (192 fragments) or B2
  (299). The DFG-in filter excludes the structures whose ligands reach there, so this is a power
  limit the design imposed on itself.
- **Prior dominance.** The n ≥ 10 flag is secondary; the primary diagnostic is the prior's share of
  posterior concentration, V·alpha / (V·alpha + N_r). Reported per cell, not globally.

## Notes

Each was a deliberate decision, not an oversight:

1. **Regions are KinFragLib subpockets, not the sequence segments the brief listed** (hinge,
   G-loop, αC, DFG, …). KLIFS does not expose subpocket labels, and the sequence-segment framing
   does not correspond to cavity geometry. Adopted on request after checking what KLIFS provides.
2. **Fragmentation uses `FindBRICSBonds` plus graph splitting, not `BRICSDecompose`.**
   `BRICSDecompose` returns SMILES and discards the 3D coordinates that region assignment needs.
   No dummy atoms are ever created; the pipeline works on atom-index sets.
3. **Fragments are assigned by centroid-to-subpocket-centre, not minimum atom-to-residue
   distance.** The brief's rule barely discriminates — 88.6% of its assignments are near-ties, at a
   median margin of 0.41 Å — because a fragment inside the cavity is at contact distance from
   residues lining several subpockets. The brief's rule is still computed and reported as a logged
   sensitivity analysis, and the disagreement between the two rules is itself a result.
4. **Category definitions were revised once**, from charged → H-bond → aromatic → aliphatic to
   charged → aromatic → H-bond → aliphatic. The original ordering sent every N-heterocycle
   (pyridine, pyrimidine, imidazole, aminopyrimidine) to the H-bond bin — a validity error in the
   specification, independent of outcome. The category count was never expanded.
5. **KinFragLib's "discard ligands with no AP fragment" rule was not adopted**; their >22-heavy-atom
   discard was. Both are logged.

## Repository layout

```
scripts/          numbered pipeline stages; common.py holds paths, throttle, region definitions
data/raw/         cached KLIFS responses (gitignored, never edited)
data/processed/   stage outputs and exclusion logs
results/          posterior summaries, pairwise tables, diagnostics
figures/          fig1_main is primary; fig1_posteriors* and figS1* are supplementary
notes/            writeup.md, writeup_plan.md (eight anchors), criticism.md
```

Environment pins that matter: `pandas < 3` (under pandas 3, opencadd's in-place cleanups fail
silently — an unknown resolution becomes 0 rather than NaN), and opencadd installed from GitHub
with `--no-deps`, because the released builds pin `biopython <= 1.77`, which has no osx-arm64
build.
