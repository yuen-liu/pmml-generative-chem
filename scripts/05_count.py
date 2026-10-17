"""Stage 5 (processing): count IFP interaction types per region, for each analysis.

An observation is one set IFP bit: interaction type v between the ligand and some residue in
region r, in one structure. (KLIFS IFPs are binary per residue x type, so a residue making two
H-bonds of the same type counts once.)

Analyses (all use the same 7 regions: AP FP SE GA B1 B2 other)
  main               one structure per (kinase, ligand); fixed 85-position lookup table
  main_perstructure  same structures; per-structure residue labels (mapping sensitivity)
  full               all structures passing the filters, no dedup; fixed table (dedup sensitivity)
  family:<name>      main set restricted to each of the 3 best-represented kinase families

Outputs
  data/processed/counts_long.csv        analysis, region, interaction, n
  results/counts_table.csv              wide table per analysis; cells < 10 flagged prior-dominated
  results/kinase_counts.csv             per-kinase structure and observation counts (imbalance)
  data/processed/exclusions_count.csv
"""
import argparse
import json

import numpy as np
import pandas as pd

from common import (INTERACTION_TYPES, N_RESIDUES, OTHER, PRIOR_DOMINATED_N, PROC, REGIONS, RESULTS, log,
                    write_exclusions)

N_FAMILIES = 3


def decode(fps):
    """Series of 595-char bit strings -> (n_structures, 85, 7) uint8 array."""
    raw = np.frombuffer("".join(fps).encode(), dtype=np.uint8) - ord("0")
    return raw.reshape(len(fps), N_RESIDUES, len(INTERACTION_TYPES))


def count_fixed(bits, pos_region):
    """bits: (n, 85, 7); pos_region: length-85 array of region labels -> DataFrame region x type."""
    per_pos = bits.sum(axis=0)                                         # 85 x 7
    df = pd.DataFrame(per_pos, columns=INTERACTION_TYPES)
    df["region"] = pos_region
    return df.groupby("region")[INTERACTION_TYPES].sum().reindex(REGIONS, fill_value=0)


def to_long(wide, analysis):
    out = wide.reset_index().melt(id_vars="region", var_name="interaction", value_name="n")
    out.insert(0, "analysis", analysis)
    return out


def main(allow_partial):
    meta_f = PROC / "subpocket_map_meta.json"
    if not meta_f.exists():
        raise SystemExit("no subpocket_map_meta.json; run 04_assign_subpockets.py first")
    meta = json.loads(meta_f.read_text())
    if meta["partial"] and not allow_partial:
        raise SystemExit(f"subpocket_map.csv was built from a partial pocket pull ({meta}); rerun "
                         "04_assign_subpockets.py after 03_pull_pockets.py finishes")
    main_set = pd.read_parquet(PROC / "structures_main.parquet")
    full = pd.read_parquet(PROC / "structures_full.parquet")
    smap = pd.read_csv(PROC / "subpocket_map.csv").sort_values("klifs_pos")
    pos_region = smap["subpocket"].to_numpy()
    ra = pd.read_parquet(PROC / "residue_assignments.parquet")

    # Main analysis uses only structures whose subpocket centres could be computed, so that the
    # fixed-table and per-structure analyses cover identical structures.
    assigned = set(ra["structure.klifs_id"].unique())
    dropped = main_set[~main_set["structure.klifs_id"].isin(assigned)]
    write_exclusions([[r["structure.klifs_id"], r["structure.pdb_id"], r["kinase.klifs_name"],
                       "no_subpocket_assignment",
                       "No per-structure subpocket assignment (see exclusions_assign.csv); dropped from main "
                       "so main and main_perstructure cover the same structures"]
                      for _, r in dropped.iterrows()], "count")
    main_set = main_set[main_set["structure.klifs_id"].isin(assigned)].reset_index(drop=True)

    parts = []
    bits_main = decode(main_set["interaction.fingerprint"].tolist())
    parts.append(to_long(count_fixed(bits_main, pos_region), "main"))
    parts.append(to_long(count_fixed(decode(full["interaction.fingerprint"].tolist()), pos_region), "full"))

    # Per-structure labels: label matrix (n, 85); positions absent from a structure have no label.
    lab = ra.pivot(index="structure.klifs_id", columns="klifs_pos", values="label")
    lab = lab.reindex(index=main_set["structure.klifs_id"], columns=range(1, N_RESIDUES + 1))
    unlabelled_bits = int(bits_main[lab.isna().to_numpy()].sum())
    if unlabelled_bits:
        log.warning("main_perstructure: %d set bits on positions with no coordinates; counted as 'other'",
                    unlabelled_bits)
    lab = lab.fillna(OTHER).to_numpy()
    rows = {r: np.zeros(len(INTERACTION_TYPES), dtype=int) for r in REGIONS}
    for r in REGIONS:
        rows[r] = bits_main[lab == r].sum(axis=0)
    per = pd.DataFrame.from_dict(rows, orient="index", columns=INTERACTION_TYPES)
    per.index.name = "region"
    parts.append(to_long(per, "main_perstructure"))

    # Family subsets (best-represented by number of main-set structures).
    fam_sizes = main_set["kinase.family"].value_counts()
    top_fams = fam_sizes.index[:N_FAMILIES].tolist()
    log.info("top families (main-set structures): %s", fam_sizes.head(N_FAMILIES).to_dict())
    for fam in top_fams:
        m = (main_set["kinase.family"] == fam).to_numpy()
        parts.append(to_long(count_fixed(bits_main[m], pos_region), f"family:{fam}"))

    long = pd.concat(parts, ignore_index=True)
    long.to_csv(PROC / "counts_long.csv", index=False)

    # Human-readable counts table with prior-dominated flags.
    wide = long.pivot_table(index=["analysis", "region"], columns="interaction", values="n", aggfunc="sum")
    wide = wide[INTERACTION_TYPES]
    wide["total"] = wide.sum(axis=1)
    wide["n_cells_lt10"] = (wide[INTERACTION_TYPES] < PRIOR_DOMINATED_N).sum(axis=1)
    flagged = wide[INTERACTION_TYPES].apply(
        lambda col: col.map(lambda n: f"{n}*" if n < PRIOR_DOMINATED_N else str(n)))
    flagged["total"] = wide["total"]
    flagged["n_cells_lt10"] = wide["n_cells_lt10"]
    flagged.to_csv(RESULTS / "counts_table.csv")
    log.info("counts written; '*' marks cells with n < %d (prior-dominated)", PRIOR_DOMINATED_N)

    # Kinase imbalance table.
    n_bits = bits_main.sum(axis=(1, 2))
    kc = (main_set.assign(n_obs=n_bits)
          .groupby(["kinase.klifs_name", "kinase.family", "kinase.group"])
          .agg(n_structures_main=("structure.klifs_id", "size"), n_obs_main=("n_obs", "sum"))
          .reset_index())
    kc = kc.merge(full.groupby("kinase.klifs_name").size().rename("n_structures_full").reset_index(),
                  on="kinase.klifs_name", how="left")
    kc = kc.sort_values("n_obs_main", ascending=False).reset_index(drop=True)
    kc["share_obs_main"] = (kc["n_obs_main"] / kc["n_obs_main"].sum()).round(4)
    kc["cum_share_obs_main"] = kc["share_obs_main"].cumsum().round(4)
    kc.to_csv(RESULTS / "kinase_counts.csv", index=False)
    log.info("kinases in main: %d; top 5 by observations:\n%s", len(kc), kc.head(5).to_string())
    log.info("kinases needed for 50%% of main observations: %d", int((kc["cum_share_obs_main"] < 0.5).sum()) + 1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--allow-partial", action="store_true", help="accept a partial lookup table (test only)")
    main(ap.parse_args().allow_partial)
