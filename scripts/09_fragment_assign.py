"""Phase 2, stage 9 (processing): BRICS-fragment each ligand, assign fragments to regions, categorise.

FRAGMENTATION. `BRICS.BRICSDecompose` returns SMILES and discards the 3D conformer the spatial
assignment needs, so BRICS bonds are found with `BRICS.FindBRICSBonds` (which `BRICSDecompose`
calls internally) and the molecule is split by removing those bonds from the atom graph. Working
on atom-index sets rather than `FragmentOnBonds` means no dummy atoms are ever created, so
centroids are over real atoms only (including hydrogens, as in KinFragLib).

ASSIGNMENT (KinFragLib, Sydow et al. 2020; code @88b61d8 fragmentation/cli.py):
  1. Each BRICS fragment's centroid (all real atoms, H included) goes to the nearest subpocket
     centre; >= OTHER_CUTOFF_A from every centre -> 'other' (their X pool).
  2. `fix_small_fragments`: a BRICS fragment with < MIN_HEAVY_ATOMS heavy atoms is RELABELLED to a
     neighbour's region (never cut away on its own). Ported from their
     fragmentation/pocketIdentification.py.
  3. Only bonds between fragments in DIFFERENT regions are cleaved, so same-region neighbours stay
     joined. Final fragments are maximal same-region connected groups. This is the merge: it moves
     chemistry between categories, so the accounting is logged (see below).
  4. Ligands with any BRICS fragment > MAX_BRICS_HEAVY heavy atoms are dropped (their rule, e.g.
     staurosporine). Their "discard ligands without an AP fragment" rule is NOT applied: it serves
     their recombination use case, not a per-region count. Logged as a deviation.

SENSITIVITY. Each final fragment is also assigned by minimum atom-to-residue distance (the
original spec's rule) and the two are compared. That rule leaves ~84% of fragments within 1 A of
a runner-up, because a fragment in the cavity is at contact distance from residues lining several
regions at once; it is reported as a sensitivity analysis, not used.

CATEGORIES (4, mutually exclusive by a stated priority rule; first match wins):
    charged    any atom with a non-zero formal charge (KLIFS mol2 carries protonation)
    hbond      any N or O atom (donor or acceptor capable)
    aromatic   any aromatic atom (so: aromatic carbocycles, since N/O-containing aromatics
               match `hbond` first)
    aliphatic  everything else
Chemical categories are not naturally exclusive (pyridine is both aromatic and an H-bond
acceptor); the Dirichlet-Categorical needs exactly one per fragment. The priority rule is a
design assumption, reported as such alongside the partition assumption.

Outputs
  data/processed/fragments.parquet         one row per final fragment
  data/processed/fragment_flags.csv        near-ties, merge ties, unfragmented ligands
  data/processed/merge_accounting.csv      how much chemistry the merge moved
  data/processed/exclusions_fragment.csv
"""
from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS

from common import (NEAR_TIE_MARGIN_A, OTHER, OTHER_CUTOFF_A, PROC, RAW, REGIONS, SUBPOCKETS, log,
                    write_exclusions)

RDLogger.DisableLog("rdApp.*")

LIGAND_DIR = RAW / "ligands"
POCKET_DIR = RAW / "pockets"
MIN_HEAVY_ATOMS = 3      # KinFragLib fix_small_fragments min_size
MAX_BRICS_HEAVY = 22     # KinFragLib large-fragment discard
CATEGORIES = ["charged", "hbond", "aromatic", "aliphatic"]
XYZ = ["atom.x", "atom.y", "atom.z"]


def components(n_atoms, bonds):
    """Connected components over `bonds` (list of (i, j)). Returns a label per atom index."""
    adj = defaultdict(list)
    for i, j in bonds:
        adj[i].append(j)
        adj[j].append(i)
    label = [-1] * n_atoms
    c = 0
    for s in range(n_atoms):
        if label[s] != -1:
            continue
        stack, label[s] = [s], c
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if label[v] == -1:
                    label[v] = c
                    stack.append(v)
        c += 1
    return label, c


def categorise(mol, atom_idx):
    """First match wins: charged -> hbond -> aromatic -> aliphatic."""
    atoms = [mol.GetAtomWithIdx(i) for i in atom_idx]
    if any(a.GetFormalCharge() != 0 for a in atoms):
        return "charged"
    if any(a.GetSymbol() in ("N", "O") for a in atoms):
        return "hbond"
    if any(a.GetIsAromatic() for a in atoms):
        return "aromatic"
    return "aliphatic"


def fix_small_fragments(regions, heavy, adj, flags, meta):
    """Relabel BRICS fragments with < MIN_HEAVY_ATOMS heavy atoms (KinFragLib, in place)."""
    for _ in range(20):
        changed = False
        for f in range(len(regions)):
            if heavy[f] >= MIN_HEAVY_ATOMS or not adj[f]:
                continue
            same = [n for n in adj[f] if regions[n] == regions[f]]
            other = [n for n in adj[f] if regions[n] != regions[f]]
            if not other:
                continue
            if same:
                group, stack = {f}, [f]          # grow the same-region group transitively
                while stack:
                    u = stack.pop()
                    for v in adj[u]:
                        if v not in group and regions[v] == regions[f]:
                            group.add(v)
                            stack.append(v)
                if sum(heavy[g] for g in group) >= MIN_HEAVY_ATOMS:
                    continue
                pool = other
            else:
                pool = adj[f]
            best = max(heavy[n] for n in pool)
            tied = [n for n in pool if heavy[n] == best]
            if len(tied) > 1:
                flags.append(meta + ["merge_tie", f"BRICS fragment {f} ({heavy[f]} heavy atoms) has "
                                                  f"{len(tied)} largest neighbours of {best} heavy atoms "
                                                  f"in regions {sorted({regions[n] for n in tied})}; "
                                                  "took the lowest fragment index"])
            regions[f] = regions[min(tied)]
            changed = True
        if not changed:
            break
    return regions


def main():
    main_set = pd.read_parquet(PROC / "structures_main.parquet")
    smap = pd.read_csv(PROC / "subpocket_map.csv").sort_values("klifs_pos")
    pos2region = dict(zip(smap["klifs_pos"], smap["subpocket"]))
    import importlib
    centres_of = importlib.import_module("04_assign_subpockets").subpocket_centres

    rows, flags, excl, merge_rows = [], [], [], []
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
        mol = Chem.MolFromMol2File(str(lig_f), removeHs=False)
        if mol is None or mol.GetNumHeavyAtoms() == 0:
            excl.append(meta + ["ligand_unparsable", "RDKit could not parse the ligand mol2"])
            continue

        poc_all = pd.read_parquet(poc_f)
        poc_all = poc_all[poc_all["residue.klifs_id"].notna()].copy()
        poc_all["res"] = poc_all["residue.klifs_id"].astype(int)
        ca = {r: g[XYZ].to_numpy(float)[0] for r, g in poc_all[poc_all["atom.name"] == "CA"].groupby("res")}
        centres, err = centres_of(ca)
        if err:
            excl.append(meta + ["anchor_residue_missing", f"Subpocket centre undefined: {err}"])
            continue
        C = np.stack([centres[s] for s in SUBPOCKETS])

        pos = mol.GetConformer().GetPositions()
        is_heavy = np.array([a.GetAtomicNum() > 1 for a in mol.GetAtoms()])
        brics = [(b[0][0], b[0][1]) for b in BRICS.FindBRICSBonds(mol)]
        brics_set = {frozenset(b) for b in brics}
        keep = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
                if frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) not in brics_set]
        label, n_frag = components(mol.GetNumAtoms(), keep)
        label = np.array(label)
        if not brics:
            flags.append(meta + ["unfragmented_ligand", f"No BRICS bond; whole ligand "
                                                        f"({int(is_heavy.sum())} heavy atoms) is one fragment"])

        heavy_counts = [int(is_heavy[label == f].sum()) for f in range(n_frag)]
        if max(heavy_counts) > MAX_BRICS_HEAVY:
            excl.append(meta + ["large_brics_fragment", f"A BRICS fragment has {max(heavy_counts)} heavy atoms "
                                                        f"(> {MAX_BRICS_HEAVY}); KinFragLib discards these ligands"])
            continue

        # 1. Assign each BRICS fragment by centroid (real atoms, H included).
        regions = []
        for f in range(n_frag):
            cen = pos[label == f].mean(axis=0)
            d = np.linalg.norm(cen - C, axis=1)
            regions.append(SUBPOCKETS[int(d.argmin())] if d.min() < OTHER_CUTOFF_A else OTHER)
        before = list(regions)

        # 2. Relabel small fragments.
        adj = defaultdict(list)
        for i, j in brics:
            a, b = label[i], label[j]
            if a != b:
                adj[a].append(b)
                adj[b].append(a)
        regions = fix_small_fragments(regions, heavy_counts, adj, flags, meta)
        moved = [f for f in range(n_frag) if regions[f] != before[f]]
        tot_heavy = int(is_heavy.sum())
        merge_rows.append({"structure.klifs_id": sid, "n_brics": n_frag, "n_relabelled": len(moved),
                           "heavy_total": tot_heavy, "heavy_relabelled": sum(heavy_counts[f] for f in moved),
                           "n_tiny_brics": sum(h < MIN_HEAVY_ATOMS for h in heavy_counts)})

        # 3. Cleave only between regions: final fragments are same-region connected groups.
        same_region_bonds = keep + [(i, j) for i, j in brics if regions[label[i]] == regions[label[j]]]
        final, n_final = components(mol.GetNumAtoms(), same_region_bonds)
        final = np.array(final)

        region_xyz = {}
        poc_h = poc_all[~poc_all["atom.name"].str.match(r"^\d*H")].copy()
        poc_h["region"] = poc_h["res"].map(pos2region)
        for r, g in poc_h.groupby("region"):
            region_xyz[r] = g[XYZ].to_numpy(float)
        present = [r for r in REGIONS if r in region_xyz]

        for k in range(n_final):
            idx = np.flatnonzero(final == k)
            hidx = [int(i) for i in idx if is_heavy[i]]
            if not hidx:
                continue
            region = regions[label[idx[0]]]
            cen = pos[idx].mean(axis=0)
            d = np.linalg.norm(cen - C, axis=1)
            order = np.argsort(d)
            margin = float(d[order[1]] - d[order[0]]) if len(d) > 1 else np.nan
            # Sensitivity: the spec's minimum atom-to-residue rule on the same fragment.
            P = pos[hidx]
            dmin = np.array([np.linalg.norm(P[:, None, :] - region_xyz[r][None], axis=2).min() for r in present])
            mo = np.argsort(dmin)
            smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=hidx)
            rows.append({"structure.klifs_id": sid, "kinase.klifs_name": st["kinase.klifs_name"],
                         "kinase.family": st["kinase.family"], "ligand.expo_id": st["ligand.expo_id"],
                         "fragment_ix": k, "n_heavy": len(hidx), "smiles": smiles,
                         "category": categorise(mol, hidx), "region": region,
                         "d_centre": float(d[order[0]]), "second": SUBPOCKETS[int(order[1])],
                         "margin": margin,
                         "region_min_atom": present[mo[0]], "d_min_atom": float(dmin[mo[0]]),
                         "margin_min_atom": float(dmin[mo[1]] - dmin[mo[0]]) if len(dmin) > 1 else np.nan,
                         "unfragmented": not brics})
            if margin is not np.nan and margin < NEAR_TIE_MARGIN_A and region != OTHER:
                flags.append(meta + ["near_tie", f"fragment {k} ({smiles}): {region} at {d[order[0]]:.2f} A vs "
                                                 f"{SUBPOCKETS[int(order[1])]} at {d[order[1]]:.2f} A "
                                                 f"(margin {margin:.2f} A)"])

    write_exclusions(excl, "fragment")
    pd.DataFrame(flags, columns=["structure.klifs_id", "structure.pdb_id", "kinase.klifs_name", "flag", "detail"]
                 ).to_csv(PROC / "fragment_flags.csv", index=False)
    mr = pd.DataFrame(merge_rows)
    mr.to_csv(PROC / "merge_accounting.csv", index=False)
    frags = pd.DataFrame(rows)
    frags.to_parquet(PROC / "fragments.parquet", index=False)

    n_lig = frags["structure.klifs_id"].nunique()
    log.info("fragments: %d from %d ligands (%.2f per ligand)", len(frags), n_lig, len(frags) / max(n_lig, 1))
    log.info("MERGE ACCOUNTING: %d of %d BRICS fragments relabelled (%.1f%%); "
             "%.1f%% of heavy atoms sat in a relabelled fragment; %.1f%% of BRICS fragments were tiny (< %d heavy)",
             mr.n_relabelled.sum(), mr.n_brics.sum(), 100 * mr.n_relabelled.sum() / max(mr.n_brics.sum(), 1),
             100 * mr.heavy_relabelled.sum() / max(mr.heavy_total.sum(), 1),
             100 * mr.n_tiny_brics.sum() / max(mr.n_brics.sum(), 1), MIN_HEAVY_ATOMS)
    log.info("near-ties (centroid rule): %.1f%% of fragments within %.1f A",
             100 * (frags["margin"] < NEAR_TIE_MARGIN_A).mean(), NEAR_TIE_MARGIN_A)
    log.info("SENSITIVITY: min-atom rule agrees on %.1f%% of fragments; its near-tie rate is %.1f%%",
             100 * (frags["region"] == frags["region_min_atom"]).mean(),
             100 * (frags["margin_min_atom"] < NEAR_TIE_MARGIN_A).mean())
    log.info("fragment size: median %.0f heavy atoms; %% below %d: %.1f%%",
             frags["n_heavy"].median(), MIN_HEAVY_ATOMS, 100 * (frags["n_heavy"] < MIN_HEAVY_ATOMS).mean())
    log.info("counts by region x category:\n%s",
             pd.crosstab(frags["region"], frags["category"]).reindex(index=REGIONS, columns=CATEGORIES)
             .fillna(0).astype(int).to_string())


if __name__ == "__main__":
    main()
