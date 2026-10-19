# Do kinase ATP-site subpockets have distinct chemical preferences?

A Dirichlet–Categorical model, fit twice, around Box's loop.

**Result.** Subpocket profiles are distinguishable — 14 of 15 pairs on protein-side interaction
patterns and 10 of 15 on ligand-side fragment chemistry, with every failure in a sparse region —
but the decomposition is partly an artefact of the assignment metric, and the intervals are
overconfident by roughly 1.6×.

---

## Model design

*In model design we bake our goals and assumptions into the model.*

The kinase ATP site is conventionally divided into six subpockets — adenine (AP), front (FP),
solvent-exposed (SE), gate area (GA) and two back pockets (B1, B2) — defined by Sydow et al.
(KinFragLib, *JCIM* 2020) as points in space: the mean Cα position of a small set of anchor
residues in KLIFS's 85-position sequence alignment, recomputed for every structure. Everything
beyond 8 Å of the nearest centre falls into a seventh region, `other`, which is kept as a
reference mixture rather than dropped, and is not a subpocket.

The question is whether the composition of what occupies each region differs between regions, so
the estimand is one composition per region. For each region *r*:

$$ \text{observation} \sim \mathrm{Categorical}(\pi_r), \qquad \pi_r \sim \mathrm{Dirichlet}(\alpha \mathbf{1}) $$

The model is stated once and fit to two observation units. **Phase 1** (protein side): each
observation is one KLIFS interaction-fingerprint bit set at one residue, typed into 7 interaction
types — 7 regions × 7 types = 49 estimands. **Phase 2** (ligand side): each observation is one
BRICS fragment of a co-crystallised ligand, assigned to the region of its nearest centre and
typed into 4 chemical categories — 28 estimands. Describing the model twice costs almost nothing,
which is what makes the two-sided comparison affordable.

The assumptions are deliberately flat and deliberately checkable: regions are independent (no
pooling), observations are exchangeable within a region, and the prior is symmetric. Exchangeability
across kinases is the one we expect to fail, and §criticism shows it does.

**Scope** (every drop logged with its rule in `exclusions_*.csv`; nothing is silently excluded).
Of 13,325 human KLIFS structures, 3,619 (27.2%) survive filters on species, kinase group, DFG-in
conformation, presence of an orthosteric non-substrate-like ligand, and a non-empty fingerprint,
followed by deduplication to one structure per (kinase, ligand) pair. Those yield 68,493 typed
interactions, and 3,424 parsable ligands yielding 10,141 fragments after KinFragLib's
small-fragment merge.

**One specification was revised, and the revision is disclosed as a fix rather than a result.**
The first category ordering (charged → H-bond-capable → aromatic → aliphatic) sent every
N-heterocycle to the H-bond bin: pyridine, pyrimidine, imidazole, aminopyrimidine — the core
hinge-binding chemistry of the entire dataset. It was reordered to charged → aromatic → H-bond →
aliphatic. The test applied was whether the change was outcome-driven: it would have been made had
the shares been perfectly even, so it is a validity error, not tuning. The category count was
never expanded.

## Posterior inference

*In posterior inference we solve the algorithmic problem of calculating the posterior.*

There is no algorithmic problem. The Dirichlet prior is conjugate to the Categorical likelihood,
so the posterior is $\pi_r \mid \text{data} \sim \mathrm{Dirichlet}(\alpha + n_r)$ in closed form,
and every marginal is Beta. Monte Carlo (40,000 draws) is used only to push that exact posterior
through functionals that have no closed form: pairwise differences $\pi_{a,v} - \pi_{b,v}$, the
total-variation distance $\tfrac12\sum_v |\pi_{a,v} - \pi_{b,v}|$ between two region profiles, and
the simulated null distributions below. Everything is reported at $\alpha = 1$ and repeated at
$\alpha = 0.5$; no conclusion in this report differs between them.

**This section has no numeric anchor, and that is the point.** Conjugacy buys exact inference by
spending modelling fidelity, and the bill comes due in criticism.

## Model criticism

*In model criticism we ask how the model has succeeded, where it has fallen short.*

### Where it succeeded

**The profiles recover known medicinal chemistry.** AP — the adenine pocket, where hinge-binding
heteroaromatics sit — is 95.3% aromatic [0.946, 0.960]. Nothing in the pipeline was fit to that
expectation, so it is external evidence that the measurement tracks something real. The
complementary contrast is equally recognisable: AP's charged share is 0.017, against 0.226 (FP)
and 0.252 (SE), the solubilising amines at the solvent front.

**Most pairs separate, judged against an explicit null.** TV distance is strictly positive by
construction, so "the credible interval excludes zero" tests nothing. The rule used instead, fixed
before any pairwise result was inspected: a pair is *resolved* if the lower end of its 95% CrI,
widened about the mean by $\sqrt{R}$ (below), lies above the 97.5% point of the TV distribution
simulated for two regions that share the pair's pooled profile at effective sample size $N/R$.
Fourteen of 15 subpocket pairs clear it in phase 1, 10 of 15 in phase 2 (Fig. 1a,b), identically at
$\alpha = 0.5$. Prior dominance does not drive this: every subpocket's prior share of posterior
concentration is below 1% in phase 1 and below 2.1% in phase 2.

**FP–SE, in one sentence.** It is the most similar pair under both measurements — rank 1 of 15 in
phase 1, rank 2 of 15 in phase 2 — and only fragment chemistry resolves it (TV 0.123, widened lower
bound 0.084 against a null of 0.060) where interaction patterns do not (0.031, 0.021 against 0.023).
This is not a power artefact: phase 1 had 3.3× more observations on this exact pair (15,107 against
4,621) and still did not resolve it, so the difference runs against the power gradient. It is best
read as a difference in resolution between two measurements, not a proven fact about the pocket:
the phases use different observation units, different vocabularies, and different nulls. Fig. 1c
shows no single category carries it — aromatic (39% of the TV) and polar non-aromatic (36%) do,
in opposition, while charged, the category one might expect at the solvent front, contributes
least (11%) and is near-equal between the two. Phase 1's undeduplicated set does resolve FP–SE,
at essentially unchanged effect size (0.034 against 0.031), passing only because pseudo-replication
shrinks the null; the deduplicated main analysis is the one reported.

### Where it fell short

**1. The decomposition is partly an artefact of the assignment metric.** Under the
centroid-to-centre rule, 12.0% of fragments lie within 1 Å of a runner-up subpocket centre. Under
the minimum atom-to-residue rule specified at the outset, 88.6% do — a fragment inside a cavity is
at contact distance from residues lining several subpockets — and the two rules disagree on the
region of 46% of all fragments. At residue level in phase 1, 25.4% of assignments are near-ties and
21 of the 85 positions have no stable label. The regions are therefore not a partition the data
impose; they are a labelling a distance rule imposes, and a different defensible rule relabels
nearly half the observations. This is the central limitation, and it is the argument for treating
subpockets as latent variables to be learned rather than as a fixed conditioning set.

**2. The intervals are too narrow, because exchangeability across kinases fails.** A posterior
predictive check with a realized discrepancy (Gelman, Meng & Stern 1996) compares the multinomial
deviance of the observed per-kinase counts to that of counts replicated from the posterior, per
draw. Deviance rather than Pearson $\chi^2$, because near-zero cells destabilise the latter. The
ratio is $R = 2.73$ [2.62, 2.85] for interactions and $2.44$ [2.31, 2.58] for fragments; on
simulated data from the fitted model it is $0.95$ [0.88, 1.04], rising to 3.7 and 18 under mild and
strong kinase-level variation. Both phases miss in the same direction by a similar amount, so the
misfit belongs to the model, not to either observation unit. Reading it as a quasi-likelihood
overdispersion — and $\sqrt{R}$ as an interval inflation factor, $N/R$ as an effective sample size,
are approximations, not derivations — the credible intervals are about 1.6× too narrow. The cause
is visible in the data: 17 of 227 kinases supply half the interactions, 17 of 214 half the
fragments. The undeduplicated set gives $R = 5.57$, which is what deduplication bought.

**3. Two regions are too sparse to compare.** All five unresolved phase-2 pairs involve B1 (192
fragments) or B2 (299); every pair involving AP or FP resolves. The sparsity is a consequence of a
design decision — the DFG-in filter excludes precisely the structures whose ligands reach the back
pocket — so the 10-of-15 result is a power limit this design imposed on itself, not evidence that
the back pockets are chemically undifferentiated.

**4. The prior still speaks where counts are thin.** The $n \ge 10$ flag is not sufficient: what
matters is the prior's share of posterior concentration, $V\alpha / (V\alpha + N_r)$, because the
prior contributes $V\alpha$ however the data concentrate. The JakA subset's B1 × apolar cell has
$n = 22$ and clears the flag, yet its posterior mean moves 0.742 → 0.818 between $\alpha = 1$ and
$0.5$, because $N_r = 24$ and the prior holds 23% of the concentration. In the limit, a zero-count
cell has posterior mean $\alpha/(V\alpha + N_r)$ — pure prior, with a uniform-to-Jeffreys ratio
tending to 2 (1.998 observed). Six main-analysis cells are of this kind and appear in the figures.
The headline results are unaffected, since they rest only on regions whose prior share is under 2%.

### Revision

The model's one structural failure points at one revision: replace independent per-region
compositions with kinase-level profiles partially pooled toward a global subpocket profile — a
hierarchical Dirichlet. That absorbs the overdispersion rather than correcting for it after the
fact, and it widens exactly the intervals that are currently too narrow. It also has no conjugate
closed form, so the next turn of the loop pays in inference what this one saved. The deeper
revision, which the metric-dependence result argues for and this model cannot reach, is to stop
conditioning on a fixed subpocket labelling at all.

---

**Figure 1.** *(a, b)* Posterior total-variation distance between region profiles for all 15
subpocket pairs, phase 1 (IFP interaction types) and phase 2 (BRICS fragment chemistry); main
analysis, $\alpha = 1$. Grey span: the 97.5% point of the same-profile null simulated at $N/R$.
Thick line: 95% CrI; thin line: the same interval widened by $\sqrt{R}$; dot: posterior mean,
hollow where the widened lower bound does not clear the null. Rows share the phase 1 ordering; the
marginal pair FP–SE is highlighted. **The two x-axes are not comparable** — the phases use
different vocabularies (7 interaction types against 4 chemical categories), which places TV on
different scales. *(c)* Per-category posterior differences $\pi_{FP} - \pi_{SE}$ (phase 2), mean
and 95% CrI, with each category's share of the pair's TV.

*Supplementary figures S1/S2 give the full region × type posterior profiles for each phase; note
that their phase 1 panels use a log scale (IFP shares span 10⁻⁴ to 0.9) and their phase 2 panels a
linear one.*
