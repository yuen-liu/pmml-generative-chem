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
alpha = 1. Per-region rows, the full set, families and Jeffreys are provenance.

Pairwise null baseline (decides anchor 3). TV distance is positive by construction, so a CrI
"clear of zero" tests nothing. For each pair, the null is the posterior-mean TV of two regions
that share one profile (the pair's pooled counts), simulated at the observed sizes N_a, N_b and
at N/R, with R this analysis's pooled overdispersion ratio. N/R as an effective sample size is
an approximation (quasi-likelihood), not a derivation. Decision rule, fixed before
looking at pairwise results: a pair is distinguishable if the lower end of its 95% CrI, widened
about the mean by sqrt(R), lies above the 97.5% point of the N/R null.

  results/ppc_overdispersion.csv
  results/pairwise_tv_null.csv
"""
import argparse

import numpy as np
import pandas as pd

from common import ALPHAS, PHASES, PROC, REGIONS, RESULTS, SUBPOCKETS, log

N_DRAWS = 2_000
SEED = 1996
CRI = (0.025, 0.975)
POOLED = "pooled_subpockets"
NULL_REPS, NULL_DRAWS = 1000, 400


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


def null_tv(n_a, n_b, n_eff_a, n_eff_b, alpha, rng):
    """TV posterior means for two regions that truly share one profile (the pair's pooled counts),
    observed with n_eff_a / n_eff_b observations. Returns the null distribution (NULL_REPS values)."""
    pi = (n_a + n_b + alpha) / (n_a + n_b + alpha).sum()
    out = np.empty(NULL_REPS)
    for i in range(NULL_REPS):
        xa = rng.multinomial(max(int(n_eff_a), 1), pi)
        xb = rng.multinomial(max(int(n_eff_b), 1), pi)
        da = rng.dirichlet(alpha + xa, size=NULL_DRAWS)
        db = rng.dirichlet(alpha + xb, size=NULL_DRAWS)
        out[i] = 0.5 * np.abs(da - db).sum(axis=1).mean()
    return out


def pairwise_null(pooled, ppc, rng, vocab, sfx):
    """Null baseline and decision rule for every region pair (see module docstring)."""
    INTERACTION_TYPES = vocab
    tv = pd.read_csv(RESULTS / f"pairwise_tv{sfx}.csv")
    rows = []
    for (analysis, prior), t in tv.groupby(["analysis", "prior"], sort=False):
        r_row = ppc[(ppc.analysis == analysis) & (ppc.prior == prior) & (ppc.region == POOLED)]
        if r_row.empty:
            continue
        R = float(r_row["ratio_median"].iloc[0])
        alpha = ALPHAS[prior]
        cnt = (pooled[pooled.analysis == analysis].pivot(index="region", columns="interaction", values="n")
               .reindex(index=REGIONS, columns=INTERACTION_TYPES).fillna(0).to_numpy())
        for _, p in t.iterrows():
            n_a, n_b = cnt[REGIONS.index(p.region_a)], cnt[REGIONS.index(p.region_b)]
            naive = null_tv(n_a, n_b, n_a.sum(), n_b.sum(), alpha, rng)
            eff = null_tv(n_a, n_b, n_a.sum() / R, n_b.sum() / R, alpha, rng)
            s = np.sqrt(R)
            lo_infl = p.tv_mean - s * (p.tv_mean - p.cri_lo)
            hi_infl = p.tv_mean + s * (p.cri_hi - p.tv_mean)
            null_q = float(np.quantile(eff, 0.975))
            rows.append(dict(
                analysis=analysis, prior=prior, region_a=p.region_a, region_b=p.region_b,
                both_subpockets=p.region_a in SUBPOCKETS and p.region_b in SUBPOCKETS,
                n_a=int(n_a.sum()), n_b=int(n_b.sum()), R=R,
                tv_mean=p.tv_mean, cri_lo=p.cri_lo, cri_hi=p.cri_hi,
                cri_lo_inflated=lo_infl, cri_hi_inflated=hi_infl,
                null_median=float(np.median(naive)), null_q975=float(np.quantile(naive, 0.975)),
                null_eff_median=float(np.median(eff)), null_eff_q975=null_q,
                distinguishable=bool(lo_infl > null_q)))
        log.info("null TV %s/%s done (R = %.2f)", analysis, prior, R)
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / f"pairwise_tv_null{sfx}.csv", index=False)
    m = out[(out.analysis == "main") & (out.prior == "uniform") & out.both_subpockets].sort_values("tv_mean")
    log.info("subpocket pairs, main, alpha = 1 (%d of %d distinguishable):\n%s",
             m.distinguishable.sum(), len(m),
             m[["region_a", "region_b", "tv_mean", "cri_lo_inflated", "null_eff_q975", "distinguishable"]]
             .round(3).to_string(index=False))


def main(phase):
    cfg = PHASES[phase]
    INTERACTION_TYPES, sfx = cfg["vocab"], cfg["suffix"]
    pooled = pd.read_csv(PROC / cfg["counts"])
    by_kin = pd.read_csv(PROC / cfg["by_kinase"])
    fam_of = (pd.read_parquet(PROC / "structures_main.parquet")
              .drop_duplicates("kinase.klifs_name").set_index("kinase.klifs_name")["kinase.family"])
    rng = np.random.default_rng(SEED)
    rows = []
    families = [a for a in pooled["analysis"].unique() if a.startswith("family:")]
    base = [a for a in ("main", "full") if a in set(pooled["analysis"])]
    for analysis in base + families:
        if analysis.startswith("family:"):
            fam = analysis.split(":", 1)[1]
            kin_rows = by_kin[(by_kin.analysis == "main") & by_kin.kinase.map(fam_of).eq(fam)]
        else:
            kin_rows = by_kin[by_kin.analysis == analysis]
        for prior, alpha in ALPHAS.items():
            tot_obs = np.zeros(N_DRAWS)
            tot_rep = np.zeros(N_DRAWS)
            pooled_kin, pooled_n = set(), 0
            for region in REGIONS:
                n_r = (pooled[(pooled.analysis == analysis) & (pooled.region == region)]
                       .set_index("interaction")["n"].reindex(INTERACTION_TYPES).fillna(0).to_numpy())
                X = (kin_rows[kin_rows.region == region]
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
    out.to_csv(RESULTS / f"ppc_overdispersion{sfx}.csv", index=False)
    a = out[(out.analysis == "main") & (out.prior == "uniform")]
    log.info("overdispersion ratio, main, alpha = 1:\n%s",
             a[["region", "n_kinases", "n_obs", "ratio_median", "ratio_cri_lo", "ratio_cri_hi",
                "approx_cri_width_inflation"]].to_string(index=False))
    pairwise_null(pooled, out, rng, INTERACTION_TYPES, sfx)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", default="ifp", choices=list(PHASES))
    main(ap.parse_args().phase)
