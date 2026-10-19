"""Phase 2, stage 10 (processing): count fragments per region x chemical category.

An observation is one BRICS fragment (after the KinFragLib small-fragment merge), assigned to one
region and one of the four mutually exclusive categories. Mirrors stage 5's output schema so the
phase 1 modelling stages (06, 06b, 07) can run on it with --phase fragment.

Analyses
  main            all fragments, centroid-to-centre assignment
  main_min_atom   same fragments, the spec's minimum atom-to-residue assignment (sensitivity)
  family:<name>   the 3 families with the most fragments

Outputs
  data/processed/counts_fragments_long.csv            analysis, region, interaction, n
  data/processed/counts_fragments_by_kinase_long.csv  same, per kinase (for the PPC)
  results/fragment_counts_table.csv                   wide table, prior share + n < 10 flags
  results/fragment_kinase_counts.csv                  per-kinase fragment counts (imbalance)
"""
import pandas as pd

from common import ALPHAS, FRAGMENT_CATEGORIES, PRIOR_DOMINATED_N, PROC, REGIONS, RESULTS, log

N_FAMILIES = 3


def wide(df, region_col="region"):
    return (pd.crosstab(df[region_col], df["category"])
            .reindex(index=REGIONS, columns=FRAGMENT_CATEGORIES).fillna(0).astype(int))


def to_long(w, analysis):
    out = w.reset_index().melt(id_vars="region", var_name="interaction", value_name="n")
    out.insert(0, "analysis", analysis)
    return out


def main():
    f = pd.read_parquet(PROC / "fragments.parquet")
    log.info("%d fragments from %d ligands, %d kinases", len(f),
             f["structure.klifs_id"].nunique(), f["kinase.klifs_name"].nunique())

    parts = [to_long(wide(f), "main")]
    m = wide(f.rename(columns={"region": "_c", "region_min_atom": "region"}))
    parts.append(to_long(m, "main_min_atom"))
    fam_sizes = f["kinase.family"].value_counts()
    log.info("top families (fragments): %s", fam_sizes.head(N_FAMILIES).to_dict())
    for fam in fam_sizes.index[:N_FAMILIES]:
        parts.append(to_long(wide(f[f["kinase.family"] == fam]), f"family:{fam}"))
    long = pd.concat(parts, ignore_index=True)
    long.to_csv(PROC / "counts_fragments_long.csv", index=False)

    by_kin = []
    for k, g in f.groupby("kinase.klifs_name"):
        by_kin.append(to_long(wide(g), "main").assign(kinase=k))
    pd.concat(by_kin, ignore_index=True).to_csv(PROC / "counts_fragments_by_kinase_long.csv", index=False)

    w = long.pivot_table(index=["analysis", "region"], columns="interaction", values="n", aggfunc="sum")
    w = w[FRAGMENT_CATEGORIES]
    w["total"] = w.sum(axis=1)
    flagged = w[FRAGMENT_CATEGORIES].apply(lambda c: c.map(lambda n: f"{n}*" if n < PRIOR_DOMINATED_N else str(n)))
    flagged["total"] = w["total"]
    V = len(FRAGMENT_CATEGORIES)
    for prior, alpha in ALPHAS.items():
        flagged[f"prior_share_{prior}"] = (V * alpha / (V * alpha + w["total"])).round(5)
    flagged["n_cells_lt10"] = (w[FRAGMENT_CATEGORIES] < PRIOR_DOMINATED_N).sum(axis=1)
    flagged.to_csv(RESULTS / "fragment_counts_table.csv")
    log.info("main counts ('*' = n < %d):\n%s", PRIOR_DOMINATED_N, flagged.loc["main"].to_string())

    kc = (f.groupby(["kinase.klifs_name", "kinase.family"]).size().rename("n_fragments")
          .reset_index().sort_values("n_fragments", ascending=False).reset_index(drop=True))
    kc["share"] = (kc["n_fragments"] / kc["n_fragments"].sum()).round(4)
    kc["cum_share"] = kc["share"].cumsum().round(4)
    kc.to_csv(RESULTS / "fragment_kinase_counts.csv", index=False)
    log.info("kinases for 50%% of fragments: %d of %d", int((kc["cum_share"] < 0.5).sum()) + 1, len(kc))


if __name__ == "__main__":
    main()
