"""Stage 6 (modelling): Dirichlet-Categorical posteriors per region, and pairwise region differences.

Model, independently for each region r:
    interaction type of each observation ~ Categorical(pi_r),   pi_r ~ Dirichlet(alpha * 1)
    =>  pi_r | data ~ Dirichlet(alpha + n_r)                         (conjugate, closed form)
Marginals are Beta(alpha + n_rv, sum(alpha + n_r) - alpha - n_rv); intervals and all functionals
(differences, total-variation distance) are computed from Monte Carlo draws of the posterior.
Run for alpha = 1 (uniform) and alpha = 0.5 (Jeffreys). No hypothesis tests.

Outputs (results/)
  posterior_summary.csv   per analysis x prior x region x interaction: mean, 95% CrI
  pairwise_diffs.csv      pi_a,v - pi_b,v: mean, 95% CrI, P(diff > 0), whether the CrI excludes 0
  pairwise_tv.csv         total-variation distance 0.5 * sum_v |pi_a,v - pi_b,v|: mean, 95% CrI
  prior_sensitivity.csv   per analysis: decision flips and max shift between the two priors
"""
from itertools import combinations

import numpy as np
import pandas as pd

from common import ALPHAS, INTERACTION_TYPES, PRIOR_DOMINATED_N, PROC, REGIONS, RESULTS, log

N_DRAWS = 40_000
SEED = 20261017
CRI = (0.025, 0.975)


def posterior_draws(counts, alpha, rng):
    """counts: region -> length-V array. Returns region -> (N_DRAWS, V) draws from Dirichlet(alpha + n)."""
    return {r: rng.dirichlet(alpha + np.asarray(n, float), size=N_DRAWS) for r, n in counts.items()}


def main():
    long = pd.read_csv(PROC / "counts_long.csv")
    rng = np.random.default_rng(SEED)
    summ, diffs, tvs = [], [], []
    for analysis, g in long.groupby("analysis", sort=False):
        wide = g.pivot(index="region", columns="interaction", values="n").reindex(
            index=REGIONS, columns=INTERACTION_TYPES).fillna(0).astype(int)
        counts = {r: wide.loc[r].to_numpy() for r in REGIONS}
        for prior, alpha in ALPHAS.items():
            draws = posterior_draws(counts, alpha, rng)
            for r in REGIONS:
                post = alpha + counts[r]
                mean_exact = post / post.sum()
                lo, hi = np.quantile(draws[r], CRI, axis=0)
                for v, it in enumerate(INTERACTION_TYPES):
                    summ.append(dict(analysis=analysis, prior=prior, alpha=alpha, region=r, interaction=it,
                                     n=int(counts[r][v]), n_region=int(counts[r].sum()),
                                     post_mean=mean_exact[v], cri_lo=lo[v], cri_hi=hi[v],
                                     prior_dominated=bool(counts[r][v] < PRIOR_DOMINATED_N)))
            for a, b in combinations(REGIONS, 2):
                d = draws[a] - draws[b]
                lo, hi = np.quantile(d, CRI, axis=0)
                for v, it in enumerate(INTERACTION_TYPES):
                    diffs.append(dict(analysis=analysis, prior=prior, region_a=a, region_b=b, interaction=it,
                                      diff_mean=d[:, v].mean(), cri_lo=lo[v], cri_hi=hi[v],
                                      p_a_gt_b=(d[:, v] > 0).mean(),
                                      cri_excludes_0=bool(lo[v] > 0 or hi[v] < 0),
                                      either_prior_dominated=bool(min(counts[a][v], counts[b][v]) < PRIOR_DOMINATED_N)))
                tv = 0.5 * np.abs(d).sum(axis=1)
                tlo, thi = np.quantile(tv, CRI)
                tvs.append(dict(analysis=analysis, prior=prior, region_a=a, region_b=b,
                                n_a=int(counts[a].sum()), n_b=int(counts[b].sum()),
                                tv_mean=tv.mean(), cri_lo=tlo, cri_hi=thi))
        log.info("modelled %s", analysis)

    summ, diffs, tvs = pd.DataFrame(summ), pd.DataFrame(diffs), pd.DataFrame(tvs)
    summ.to_csv(RESULTS / "posterior_summary.csv", index=False)
    diffs.to_csv(RESULTS / "pairwise_diffs.csv", index=False)
    tvs.to_csv(RESULTS / "pairwise_tv.csv", index=False)

    # Prior sensitivity: do any "CrI excludes 0" conclusions change between alpha = 1 and 0.5?
    key = ["analysis", "region_a", "region_b", "interaction"]
    w = diffs.pivot_table(index=key, columns="prior", values="cri_excludes_0").reset_index()
    flips = w[w["uniform"] != w["jeffreys"]]
    flips.merge(diffs[diffs.prior == "uniform"][key + ["either_prior_dominated"]], on=key).to_csv(
        RESULTS / "prior_sensitivity_flips.csv", index=False)
    pm = summ.pivot_table(index=["analysis", "region", "interaction"], columns="prior", values="post_mean")
    shift = (pm["uniform"] - pm["jeffreys"]).abs().groupby(level="analysis").max()
    sens = pd.DataFrame({
        "n_pairwise_comparisons": w.groupby("analysis").size(),
        "n_conclusion_flips": flips.groupby("analysis").size(),
        "max_abs_posterior_mean_shift": shift,
    }).fillna(0).astype({"n_conclusion_flips": int})
    sens.to_csv(RESULTS / "prior_sensitivity.csv")
    log.info("prior sensitivity:\n%s", sens.to_string())


if __name__ == "__main__":
    main()
