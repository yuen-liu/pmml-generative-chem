"""Stage 6b (model criticism): posterior predictive check for overdispersion across kinases.

The model treats every IFP observation in region r as an independent draw from pi_r. If
kinases have different profiles, per-kinase counts vary more than a multinomial allows, and
the Dirichlet credible intervals are too narrow.

Discrepancy (realized, parameter-dependent; Gelman, Meng & Stern 1996): the multinomial
deviance of per-kinase counts around pi,
    D(X; pi) = 2 * sum_k sum_v X_kv * log(X_kv / (N_k * pi_v)),   0 * log 0 = 0.
The deviance is used instead of Pearson chi-square, which is unstable for near-zero pi_v
(the zero-count charge cells). For each posterior draw pi^(s) ~ Dirichlet(alpha + n_r):
    D_obs^(s) = D(X_obs; pi^(s)),  X_rep^(s)_k ~ Multinomial(N_k, pi^(s)),  D_rep^(s) = D(X_rep^(s); pi^(s))
and the overdispersion ratio R^(s) = D_obs^(s) / D_rep^(s). R near 1: kinase-level variation is
what the model expects. R >> 1: the independence assumption fails; credible intervals are too
narrow by roughly sqrt(R) (quasi-likelihood approximation). The distribution of R is
reported; no tail probability is computed.

Anchor: R pooled over the six subpockets (sum of D_obs / sum of D_rep per draw), main analysis,
alpha = 1. Per-region rows, the full set and Jeffreys are provenance.

  results/ppc_overdispersion.csv
"""
import numpy as np
import pandas as pd

from common import ALPHAS, INTERACTION_TYPES, PROC, REGIONS, RESULTS, SUBPOCKETS, log

N_DRAWS = 2_000
SEED = 1996
CRI = (0.025, 0.975)
POOLED = "pooled_subpockets"


def deviance(X, E):
    """X, E: (..., K, V). Multinomial deviance summed over kinases and types."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(X > 0, X * np.log(X / E), 0.0)
    return 2.0 * t.sum(axis=(-2, -1))


def region_ratios(X, n_pooled, alpha, rng):
    """X: (K, V) observed per-kinase counts. Returns (D_obs, D_rep) arrays over draws, or None."""
    X = X[X.sum(axis=1) > 0]
    if len(X) < 2:
        return None
    Nk = X.sum(axis=1)
    pis = rng.dirichlet(alpha + n_pooled.astype(float), size=N_DRAWS)        # S x V
    d_obs = np.empty(N_DRAWS)
    d_rep = np.empty(N_DRAWS)
    for s, pi in enumerate(pis):
        E = Nk[:, None] * pi[None, :]
        d_obs[s] = deviance(X, E)
        d_rep[s] = deviance(rng.multinomial(Nk, pi).astype(float), E)
    return d_obs, d_rep, len(X), int(Nk.sum())


def summarize(analysis, prior, region, d_obs, d_rep, n_kin, n_obs):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = d_obs / d_rep
    r = r[np.isfinite(r)]
    lo, hi = np.quantile(r, CRI)
    med = float(np.median(r))
    return dict(analysis=analysis, prior=prior, region=region, n_kinases=n_kin, n_obs=n_obs,
                d_obs_median=float(np.median(d_obs)), d_rep_median=float(np.median(d_rep)),
                ratio_median=med, ratio_cri_lo=lo, ratio_cri_hi=hi,
                approx_cri_width_inflation=float(np.sqrt(med)), n_valid_draws=len(r))


def main():
    pooled = pd.read_csv(PROC / "counts_long.csv")
    by_kin = pd.read_csv(PROC / "counts_by_kinase_long.csv")
    rng = np.random.default_rng(SEED)
    rows = []
    for analysis in ("main", "full"):
        for prior, alpha in ALPHAS.items():
            tot_obs = np.zeros(N_DRAWS)
            tot_rep = np.zeros(N_DRAWS)
            pooled_kin, pooled_n = set(), 0
            for region in REGIONS:
                n_r = (pooled[(pooled.analysis == analysis) & (pooled.region == region)]
                       .set_index("interaction")["n"].reindex(INTERACTION_TYPES).fillna(0).to_numpy())
                X = (by_kin[(by_kin.analysis == analysis) & (by_kin.region == region)]
                     .pivot(index="kinase", columns="interaction", values="n")
                     .reindex(columns=INTERACTION_TYPES).fillna(0).to_numpy())
                res = region_ratios(X, n_r, alpha, rng)
                if res is None:
                    log.warning("%s/%s/%s: fewer than 2 kinases with observations; skipped", analysis, prior, region)
                    continue
                d_obs, d_rep, n_kin, n_obs = res
                rows.append(summarize(analysis, prior, region, d_obs, d_rep, n_kin, n_obs))
                if region in SUBPOCKETS:
                    tot_obs += d_obs
                    tot_rep += d_rep
                    pooled_n += n_obs
                    pooled_kin.add(n_kin)
            rows.append(summarize(analysis, prior, POOLED, tot_obs, tot_rep, max(pooled_kin, default=0), pooled_n))
            log.info("PPC %s/%s done", analysis, prior)
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "ppc_overdispersion.csv", index=False)
    a = out[(out.analysis == "main") & (out.prior == "uniform")]
    log.info("overdispersion ratio, main, alpha = 1:\n%s",
             a[["region", "n_kinases", "n_obs", "ratio_median", "ratio_cri_lo", "ratio_cri_hi",
                "approx_cri_width_inflation"]].to_string(index=False))


if __name__ == "__main__":
    main()
