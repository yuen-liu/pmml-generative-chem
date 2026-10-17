"""Phase 2, stage 9 (processing): BRICS-fragment each ligand and assign fragments to regions.

Fragmentation. `BRICS.BRICSDecompose` returns SMILES and discards 3D coordinates, which the
spatial assignment needs, so fragments are cut with `BRICS.FindBRICSBonds` +
`Chem.FragmentOnBonds` instead. These are the same BRICS bonds `BRICSDecompose` would cut
(it calls FindBRICSBonds internally), but the conformer is kept. One difference: BRICSDecompose
deduplicates and returns unique fragment SMILES, whereas splitting on bonds keeps every piece
of every ligand, which is what a per-region count needs. KinFragLib fragments the same way.

Assignment. For each fragment, the distance to a region is the minimum distance between any
fragment heavy atom and any heavy atom of any pocket residue labelled with that region (the
fixed lookup table from stage 4). The fragment goes to the nearest region. This is the spec's
"minimum atom-to-residue distance", and matches the residue-level rule used in phase 1.

Logged (these matter for the writeup):
  near-tie        nearest and second-nearest regions within NEAR_TIE_MARGIN_A
  tiny fragment   < MIN_HEAVY_ATOMS heavy atoms (linkers; assignment is noisy)
  unfragmented    ligand has no BRICS bond, so the whole ligand is one fragment
  exclusions      ligand file missing, unparsable, or pocket coordinates missing

Outputs
  data/processed/fragments.parquet        one row per fragment, with its region and distances
  data/processed/fragment_flags.csv       near-ties, tiny fragments, unfragmented ligands
  data/processed/exclusions_fragment.csv
"""
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS

from common import NEAR_TIE_MARGIN_A, PROC, RAW, REGIONS, log, write_exclusions

RDLogger.DisableLog("rdApp.*")

LIGAND_DIR = RAW / "ligands"
POCKET_DIR = RAW / "pockets"
MIN_HEAVY_ATOMS = 3
XYZ = ["atom.x", "atom.y", "atom.z"]


def brics_pieces(mol):
    """Split on BRICS bonds, keeping coordinates. Returns [(piece_mol, heavy_atom_indices), ...]."""
    bonds = list(BRICS.FindBRICSBonds(mol))
    if not bonds:
        return [(mol, list(range(mol.GetNumAtoms())))], 0
    bidx = [mol.GetBondBetweenAtoms(a, b).GetIdx() for (a, b), _ in bonds]
    split = Chem.FragmentOnBonds(mol, bidx, addDummies=True)
    frags = Chem.GetMolFrags(split, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=[])
    mapping = []
    Chem.GetMolFrags(split, asMols=False)
    idx_groups = Chem.GetMolFrags(split)
    for piece, idxs in zip(frags, idx_groups):
        # Keep only indices that exist in the original molecule (dummies are appended at the end).
        orig = [i for i in idxs if i < mol.GetNumAtoms()]
        mapping.append((piece, orig))
    return mapping, len(bonds)


def main():
    main_set = pd.read_parquet(PROC / "structures_main.parquet")
    smap = pd.read_csv(PROC / "subpocket_map.csv").sort_values("klifs_pos")
    pos2region = dict(zip(smap["klifs_pos"], smap["subpocket"]))

    rows, flags, excl = [], [], []
    for _, st in main_set.iterrows():
        sid = int(st["structure.klifs_id"])
        meta = [sid, st["structure.pdb_id"], st["kinase.klifs_name"]]
        lig_f, poc_f = LIGAND_DIR / f"{sid}.mol2", POCKET_DIR / f"{sid}.parquet"
        if not lig_f.exists():
            excl.append(meta + ["ligand_unavailable", "Ligand mol2 not fetched (see data/raw/ligands/_failures.csv)"])
            continue
        if not poc_f.exists():
            excl.append(meta + ["pocket_unavailable", "Pocket coordinates not fetched"])
            continue
        mol = Chem.MolFromMol2File(str(lig_f), removeHs=True)
        if mol is None or mol.GetNumHeavyAtoms() == 0:
            excl.append(meta + ["ligand_unparsable", "RDKit could not parse the ligand mol2"])
            continue

        # Pocket heavy atoms grouped by region (fixed lookup table).
        poc = pd.read_parquet(poc_f)
        poc = poc[poc["residue.klifs_id"].notna()].copy()
        poc = poc[~poc["atom.name"].str.match(r"^\d*H")]
        poc["region"] = poc["residue.klifs_id"].astype(int).map(pos2region)
        region_xyz = {r: g[XYZ].to_numpy(float) for r, g in poc.groupby("region") if len(g)}
        present = [r for r in REGIONS if r in region_xyz]
        if not present:
            excl.append(meta + ["no_pocket_regions", "No labelled pocket residues with coordinates"])
            continue

        conf = mol.GetConformer()
        pieces, n_bonds = brics_pieces(mol)
        if n_bonds == 0:
            flags.append(meta + ["unfragmented_ligand", f"No BRICS bond; whole ligand ({mol.GetNumHeavyAtoms()} "
                                                        "heavy atoms) treated as one fragment"])
        for k, (piece, orig_idx) in enumerate(pieces):
            heavy = [i for i in orig_idx if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
            if not heavy:
                continue
            P = np.array([list(conf.GetAtomPosition(i)) for i in heavy])
            d = np.array([np.linalg.norm(P[:, None, :] - region_xyz[r][None], axis=2).min() for r in present])
            order = np.argsort(d)
            nearest, d1 = present[order[0]], d[order[0]]
            second = present[order[1]] if len(present) > 1 else None
            d2 = d[order[1]] if len(present) > 1 else np.nan
            smiles = Chem.MolToSmiles(Chem.RemoveHs(piece, sanitize=False))
            rows.append(dict(**{"structure.klifs_id": sid, "kinase.klifs_name": st["kinase.klifs_name"],
                                "kinase.family": st["kinase.family"], "ligand.expo_id": st["ligand.expo_id"]},
                             fragment_ix=k, n_heavy=len(heavy), smiles=smiles,
                             region=nearest, d_nearest=d1, second=second, d_second=d2,
                             margin=d2 - d1, unfragmented=n_bonds == 0))
            if second is not None and (d2 - d1) < NEAR_TIE_MARGIN_A:
                flags.append(meta + ["near_tie", f"fragment {k} ({smiles}): {nearest} at {d1:.2f} A vs "
                                                 f"{second} at {d2:.2f} A (margin {d2 - d1:.2f} A)"])
            if len(heavy) < MIN_HEAVY_ATOMS:
                flags.append(meta + ["tiny_fragment", f"fragment {k} ({smiles}) has {len(heavy)} heavy atoms; "
                                                      "assignment is noisy (linker-sized)"])

    write_exclusions(excl, "fragment")
    pd.DataFrame(flags, columns=["structure.klifs_id", "structure.pdb_id", "kinase.klifs_name", "flag", "detail"]
                 ).to_csv(PROC / "fragment_flags.csv", index=False)
    frags = pd.DataFrame(rows)
    frags.to_parquet(PROC / "fragments.parquet", index=False)

    log.info("fragments: %d from %d ligands (%.1f per ligand)", len(frags),
             frags["structure.klifs_id"].nunique(), len(frags) / max(frags["structure.klifs_id"].nunique(), 1))
    log.info("heavy-atom size: median %.0f, %% with < %d: %.1f%%", frags["n_heavy"].median(), MIN_HEAVY_ATOMS,
             100 * (frags["n_heavy"] < MIN_HEAVY_ATOMS).mean())
    log.info("near-ties: %.1f%% of fragments (< %.1f A margin)",
             100 * (frags["margin"] < NEAR_TIE_MARGIN_A).mean(), NEAR_TIE_MARGIN_A)
    log.info("fragments per region:\n%s", frags["region"].value_counts().reindex(REGIONS).fillna(0).astype(int).to_string())
    log.info("unfragmented ligands: %d", frags.loc[frags["unfragmented"], "structure.klifs_id"].nunique())


if __name__ == "__main__":
    main()
