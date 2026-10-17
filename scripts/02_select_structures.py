"""Stage 2 (processing): apply scope filters and deduplication; log every excluded structure.

Filters follow KinFragLib's structure selection (KinaseFocusedFragmentLibrary @88b61d8,
preprocessing/cli.py lines 32-41 and preprocessing/discard.py), because the subpocket scheme
is only defined on that domain. Each structure is logged under the FIRST rule it fails.

Outputs
  data/processed/structures_full.parquet  all structures passing the filters (sensitivity set)
  data/processed/structures_main.parquet  one structure per (kinase, ligand) (main analysis)
  data/processed/exclusions_select.csv, exclusions_dedup.csv
  data/processed/flags_select.csv          kept structures with an unverifiable check
"""
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from common import N_RESIDUES, PROC, RAW, INTERACTION_TYPES, log, write_exclusions

RDLogger.DisableLog("rdApp.*")

# Verbatim from KinFragLib preprocessing/discard.py (substrate filter).
PHOSPHATES = [Chem.MolFromSmiles(s) for s in ("O[P](O)(O)O", "O[P](=O)(O)O", "O[P+](O)(O)O")]
RIBOSE = Chem.MolFromSmiles("OC1COCC1O")


def is_substrate(smiles):
    """True/False per KinFragLib's phosphate/ribose rule; None if the SMILES is missing or unparsable."""
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return any(mol.HasSubstructMatch(p) for p in PHOSPHATES) or mol.HasSubstructMatch(RIBOSE)


def main():
    st = pd.read_parquet(RAW / "structures_all.parquet")
    kin = pd.read_parquet(RAW / "kinases_detail.parquet")[["kinase.klifs_id", "kinase.group", "kinase.family"]]
    lig = pd.read_parquet(RAW / "ligands_all.parquet")[["ligand.klifs_id", "ligand.smiles"]]
    ifp = pd.read_parquet(RAW / "ifps_human.parquet")

    st = st.drop(columns=["kinase.group", "kinase.family", "interaction.fingerprint"])  # empty in this endpoint
    st = st.merge(kin, on="kinase.klifs_id", how="left")
    st = st.merge(lig, on="ligand.klifs_id", how="left")
    st = st.merge(ifp, on="structure.klifs_id", how="left")
    fp = st["interaction.fingerprint"]
    fp_ok = fp.notna() & (fp.astype(str).str.len() == N_RESIDUES * len(INTERACTION_TYPES))
    has_ligand = ~st["ligand.expo_id"].astype(str).isin(["-", "0", "", "None", "<NA>"])
    substrate = st["ligand.smiles"].map(is_substrate)

    rules = [
        ("non_human", st["species.klifs"] != "Human",
         "Not a human kinase (KinFragLib scope; the question concerns the human kinome)"),
        ("atypical_kinase", st["kinase.group"] == "Atypical",
         "Atypical kinase group; excluded by KinFragLib (preprocessing/cli.py:35)"),
        ("dfg_not_in", st["structure.dfg"] != "in",
         "DFG conformation is not 'in' (out / out-like / na); KinFragLib subpocket centres are "
         "defined on DFG-in structures only (preprocessing/cli.py:41)"),
        ("no_orthosteric_ligand", ~has_ligand,
         "No orthosteric ligand in the ATP site, so KLIFS provides no IFP"),
        ("ifp_missing", has_ligand & ~fp_ok,
         "Orthosteric ligand present but KLIFS returned no valid 595-bit IFP"),
        ("ifp_all_zero", fp_ok & (fp.astype(str).str.count("1") == 0),
         "IFP has no interactions set (all 595 bits zero)"),
        ("substrate_like_ligand", substrate.eq(True),
         "Ligand contains a phosphate or ribose (substrate/cofactor); KinFragLib discard.py rule"),
    ]

    remaining = pd.Series(True, index=st.index)
    excl = []
    for rule, mask, reason in rules:
        hit = remaining & mask.fillna(False).astype(bool)
        for _, r in st[hit].iterrows():
            excl.append([r["structure.klifs_id"], r["structure.pdb_id"], r["kinase.klifs_name"], rule, reason])
        remaining &= ~hit
    write_exclusions(excl, "select")

    full = st[remaining].copy()
    flags = []
    for _, r in full[substrate[remaining].isna()].iterrows():
        flags.append([r["structure.klifs_id"], r["structure.pdb_id"], r["kinase.klifs_name"],
                      "substrate_check_unverifiable", "No parsable ligand SMILES; kept without substrate check"])
    for _, r in full[full["kinase.group"].isna()].iterrows():
        flags.append([r["structure.klifs_id"], r["structure.pdb_id"], r["kinase.klifs_name"],
                      "group_unknown", "Kinase group missing; Atypical check not applied"])
    pd.DataFrame(flags, columns=["structure.klifs_id", "structure.pdb_id", "kinase.klifs_name", "flag", "reason"]
                 ).to_csv(PROC / "flags_select.csv", index=False)
    log.info("kept-with-flag: %d", len(flags))

    # Main set: one structure per (kinase, ligand). Best KLIFS quality score, then best resolution,
    # then fewest missing residues, then lowest KLIFS ID (deterministic).
    full["_res"] = full["structure.resolution"].fillna(np.inf)
    order = full.sort_values(
        ["kinase.klifs_id", "ligand.expo_id", "structure.qualityscore", "_res", "structure.missing_residues",
         "structure.klifs_id"],
        ascending=[True, True, False, True, True, True])
    keep = ~order.duplicated(["kinase.klifs_id", "ligand.expo_id"], keep="first")
    main_set = order[keep]
    kept_for = main_set.set_index(["kinase.klifs_id", "ligand.expo_id"])["structure.klifs_id"]
    dedup = []
    for _, r in order[~keep].iterrows():
        k = kept_for[(r["kinase.klifs_id"], r["ligand.expo_id"])]
        dedup.append([r["structure.klifs_id"], r["structure.pdb_id"], r["kinase.klifs_name"],
                      "duplicate_kinase_ligand",
                      f"Same (kinase, ligand {r['ligand.expo_id']}) as kept structure {k}; "
                      "lower quality score / resolution (or alternate model / chain)"])
    write_exclusions(dedup, "dedup")

    cols = [c for c in full.columns if c != "_res"]
    full[cols].to_parquet(PROC / "structures_full.parquet", index=False)
    main_set[cols].to_parquet(PROC / "structures_main.parquet", index=False)
    for name, d in (("full", full), ("main", main_set)):
        log.info("%s set: %d structures, %d kinases, %d families, %d ligands", name, len(d),
                 d["kinase.klifs_id"].nunique(), d["kinase.family"].nunique(), d["ligand.expo_id"].nunique())


if __name__ == "__main__":
    main()
