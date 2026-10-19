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
import argparse
from itertools import combinations

import numpy as np
import pandas as pd

from common import ALPHAS, PHASES, PRIOR_DOMINATED_N, PROC, REGIONS, RESULTS, log

N_DRAWS = 40_000
SEED = 20261017
CRI = (0.025, 0.975)


def posterior_draws(counts, alpha, rng):
    """counts: region -> length-V array. Returns region -> (N_DRAWS, V) draws from Dirichlet(alpha + n)."""
    return {r: rng.dirichlet(alpha + np.asarray(n, float), size=N_DRAWS) for r, n in counts.items()}


def main(phase):
    cfg = PHASES[phase]
    INTERACTION_TYPES, sfx = cfg["vocab"], cfg["suffix"]
    long = pd.read_csv(PROC / cfg["counts"])
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
    summ.to_csv(RESULTS / f"posterior_summary{sfx}.csv", index=False)
    diffs.to_csv(RESULTS / f"pairwise_diffs{sfx}.csv", index=False)
    tvs.to_csv(RESULTS / f"pairwise_tv{sfx}.csv", index=False)

    # Prior sensitivity: do any "CrI excludes 0" conclusions change between alpha = 1 and 0.5?
    key = ["analysis", "region_a", "region_b", "interaction"]
    w = diffs.pivot_table(index=key, columns="prior", values="cri_excludes_0").reset_index()
    flips = w[w["uniform"] != w["jeffreys"]]
    flips.merge(diffs[diffs.prior == "uniform"][key + ["either_prior_dominated"]], on=key).to_csv(
        RESULTS / f"prior_sensitivity_flips{sfx}.csv", index=False)
    pm = summ.pivot_table(index=["analysis", "region", "interaction"], columns="prior", values="post_mean")
    shift = (pm["uniform"] - pm["jeffreys"]).abs().groupby(level="analysis").max()
    sens = pd.DataFrame({
        "n_pairwise_comparisons": w.groupby("analysis").size(),
        "n_conclusion_flips": flips.groupby("analysis").size(),
        "max_abs_posterior_mean_shift": shift,
    }).fillna(0).astype({"n_conclusion_flips": int})
    sens.to_csv(RESULTS / f"prior_sensitivity{sfx}.csv")
    log.info("prior sensitivity (global):\n%s", sens.to_string())

    # Per-cell sensitivity, conditional on the cell count. For n_rv = 0 the posterior mean is
    # alpha / (V*alpha + N_r): pure prior, and the two priors differ by ~2x when N_r >> V.
    idx = ["analysis", "region", "interaction", "n", "n_region"]
    cell = summ.pivot_table(index=idx, columns="prior", values=["post_mean", "cri_hi"])
    cell.columns = [f"{v}_{p}" for v, p in cell.columns]
    cell = cell.reset_index()
    cell["mean_ratio_uniform_over_jeffreys"] = cell["post_mean_uniform"] / cell["post_mean_jeffreys"]
    cell["mean_abs_diff"] = (cell["post_mean_uniform"] - cell["post_mean_jeffreys"]).abs()
    cell["cri_hi_ratio_uniform_over_jeffreys"] = cell["cri_hi_uniform"] / cell["cri_hi_jeffreys"]
    # Share of the posterior's total concentration contributed by the prior. A cell can clear
    # n >= 10 and still move with the prior if its region's total N_r is small.
    V = len(INTERACTION_TYPES)
    for prior, alpha in ALPHAS.items():
        cell[f"prior_weight_{prior}"] = V * alpha / (V * alpha + cell["n_region"])
    cell["n_bin"] = pd.cut(cell["n"], [-1, 0, 9, 99, 999, np.inf], labels=["0", "1-9", "10-99", "100-999", ">=1000"])
    cell.sort_values(["analysis", "n"]).to_csv(RESULTS / f"prior_sensitivity_by_cell{sfx}.csv", index=False)
    by_bin = cell.groupby(["analysis", "n_bin"], observed=True).agg(
        n_cells=("n", "size"),
        median_mean_ratio=("mean_ratio_uniform_over_jeffreys", "median"),
        max_mean_ratio=("mean_ratio_uniform_over_jeffreys", "max"),
        max_mean_abs_diff=("mean_abs_diff", "max"),
        median_cri_hi_ratio=("cri_hi_ratio_uniform_over_jeffreys", "median"))
    by_bin.to_csv(RESULTS / f"prior_sensitivity_by_count{sfx}.csv")
    log.info("prior sensitivity by cell count (main):\n%s",
             by_bin.loc["main"].to_string() if "main" in by_bin.index.get_level_values(0) else by_bin.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", default="ifp", choices=list(PHASES))
    main(ap.parse_args().phase)
