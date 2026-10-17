"""Stage 4 (processing): assign each of the 85 KLIFS positions to a KinFragLib subpocket or 'other'.

Per structure:
  1. Subpocket centre = mean C-alpha position of its anchor residues (common.SUBPOCKET_ANCHORS).
     Missing anchor C-alpha: use the mean of the C-alphas of residues i-1 and i+1, or the one that
     is present; if both are missing the structure is skipped. This is KinFragLib's own fallback
     (fragmentation/pocketIdentification.py, calc_subpocket_center).
  2. Residue-to-centre distance = distance from the centre to the residue's CLOSEST heavy atom
     (the same minimum atom-to-residue geometry as phase 2's fragment assignment). Each residue goes
     to its nearest centre; residues >= 8 A from every centre become 'other' (KinFragLib's X-pool
     cutoff, applied to residues). Using heavy-atom centroids instead puts 45 positions / 20% of IFP
     bits in 'other' versus 27 / 10% (checked on the first 120 structures before the profiles were looked at).
Across structures:
  3. The fixed lookup table used in the main analysis is the modal label per position.

Outputs
  data/processed/residue_assignments.parquet  per-structure labels (sensitivity analysis)
  data/processed/subpocket_map.csv             fixed 85-position lookup table + stability
  data/processed/near_ties.csv                 near-tie and cutoff-borderline assignments
  data/processed/exclusions_assign.csv
  results/subpocket_map.md                     the lookup table with citation (for the writeup)
"""
import argparse
import json

import numpy as np
import pandas as pd

from common import (NEAR_TIE_MARGIN_A, OTHER, OTHER_CUTOFF_A, PROC, RAW, RESULTS, SUBPOCKET_ANCHORS,
                    SUBPOCKETS, log, write_exclusions)

POCKET_DIR = RAW / "pockets"
CUTOFF_BORDER_A = 0.5  # |d_nearest - 8 A| below this is logged as borderline for 'other'

# KLIFS sequence segments, for display only (KLIFS pocket endpoint, residue.klifs_region;
# van Linden et al., J. Med. Chem. 2014, 57, 249-277). Not used for the analysis.
KLIFS_SEGMENTS = [("I", 1, 3), ("g.l", 4, 9), ("II", 10, 13), ("III", 14, 19), ("αC", 20, 30),
                  ("b.l", 31, 37), ("IV", 38, 41), ("V", 42, 44), ("GK", 45, 45), ("hinge", 46, 48),
                  ("linker", 49, 52), ("αD", 53, 59), ("αE", 60, 64), ("VI", 65, 67), ("c.l", 68, 75),
                  ("VII", 76, 78), ("VIII", 79, 79), ("xDFG", 80, 83), ("a.l", 84, 85)]


def segment(res):
    return next(name for name, lo, hi in KLIFS_SEGMENTS if lo <= res <= hi)


def subpocket_centres(ca):
    """ca: {klifs_id: xyz}. Returns ({subpocket: xyz}, None) or (None, reason)."""
    centres = {}
    for sp, anchors in SUBPOCKET_ANCHORS.items():
        pts = []
        for res in anchors:
            if res in ca:
                pts.append(ca[res])
                continue
            nei = [ca[r] for r in (res - 1, res + 1) if r in ca]
            if not nei:
                return None, f"anchor residue {res} ({sp}) and both neighbours missing"
            pts.append(np.mean(nei, axis=0))
        centres[sp] = np.mean(pts, axis=0)
    return centres, None


def assign_structure(pocket):
    pocket = pocket[pocket["residue.klifs_id"].notna()].copy()
    pocket["res"] = pocket["residue.klifs_id"].astype(int)
    xyz_cols = ["atom.x", "atom.y", "atom.z"]
    ca = {r: g[xyz_cols].to_numpy(float)[0] for r, g in pocket[pocket["atom.name"] == "CA"].groupby("res")}
    centres, err = subpocket_centres(ca)
    if err:
        return None, err
    heavy = pocket[~pocket["atom.name"].str.match(r"^\d*H")]
    C = np.stack([centres[sp] for sp in SUBPOCKETS])                                  # 6 x 3
    atom_d = np.linalg.norm(heavy[xyz_cols].to_numpy(float)[:, None, :] - C[None], axis=2)
    dmin = pd.DataFrame(atom_d, columns=SUBPOCKETS).groupby(heavy["res"].to_numpy()).min()
    D = dmin.to_numpy()                                                               # n_res x 6
    order = np.argsort(D, axis=1)
    d1 = D[np.arange(len(D)), order[:, 0]]
    d2 = D[np.arange(len(D)), order[:, 1]]
    out = pd.DataFrame({
        "klifs_pos": dmin.index.to_numpy(),
        "nearest": np.array(SUBPOCKETS)[order[:, 0]],
        "d_nearest": d1,
        "second": np.array(SUBPOCKETS)[order[:, 1]],
        "d_second": d2,
    })
    out["margin"] = out["d_second"] - out["d_nearest"]
    out["label"] = np.where(out["d_nearest"] < OTHER_CUTOFF_A, out["nearest"], OTHER)
    return out, None


def main(allow_partial):
    st = pd.read_parquet(PROC / "structures_main.parquet")
    # The lookup table must be built from the whole main set. Refuse to run while the pocket pull
    # is incomplete: structures neither cached nor recorded as failed fetches mean it is still going.
    failed_f = POCKET_DIR / "_failures.csv"
    failed = set(pd.read_csv(failed_f)["structure.klifs_id"]) if failed_f.exists() else set()
    pending = [s for s in st["structure.klifs_id"].astype(int)
               if not (POCKET_DIR / f"{s}.parquet").exists() and s not in failed]
    if pending and not allow_partial:
        raise SystemExit(f"{len(pending)} main-set pockets not yet fetched; finish 03_pull_pockets.py "
                         "(or pass --allow-partial for a test run, whose outputs are marked partial)")
    excl, parts = [], []
    for _, r in st.iterrows():
        sid = int(r["structure.klifs_id"])
        f = POCKET_DIR / f"{sid}.parquet"
        meta = [sid, r["structure.pdb_id"], r["kinase.klifs_name"]]
        if not f.exists():
            excl.append(meta + ["pocket_unavailable", "Pocket coordinates not fetched (see data/raw/pockets/_failures.csv)"])
            continue
        res, err = assign_structure(pd.read_parquet(f))
        if err:
            excl.append(meta + ["anchor_residue_missing", f"Subpocket centre undefined: {err} (KinFragLib skips these)"])
            continue
        res.insert(0, "structure.klifs_id", sid)
        parts.append(res)
    write_exclusions(excl, "assign")
    ra = pd.concat(parts, ignore_index=True)
    ra.to_parquet(PROC / "residue_assignments.parquet", index=False)
    n_struct = ra["structure.klifs_id"].nunique()
    log.info("assigned %d structures", n_struct)

    # Near-ties: nearest and second-nearest centres within NEAR_TIE_MARGIN_A, for residues not in 'other'.
    tie = ra[(ra["label"] != OTHER) & (ra["margin"] < NEAR_TIE_MARGIN_A)].assign(
        issue=f"near_tie (<{NEAR_TIE_MARGIN_A} A between nearest and second centre)")
    border = ra[(ra["d_nearest"] - OTHER_CUTOFF_A).abs() < CUTOFF_BORDER_A].assign(
        issue=f"cutoff_borderline (nearest centre within {CUTOFF_BORDER_A} A of the {OTHER_CUTOFF_A} A cutoff)")
    ties = pd.concat([tie, border], ignore_index=True)
    ties.to_csv(PROC / "near_ties.csv", index=False)
    log.info("near-ties: %d residue assignments; cutoff-borderline: %d", len(tie), len(border))

    # Fixed lookup table: modal label per position across structures.
    freq = ra.groupby(["klifs_pos", "label"]).size().unstack(fill_value=0)
    share = freq.div(freq.sum(axis=1), axis=0)
    top2 = share.apply(lambda row: row.sort_values(ascending=False).iloc[:2], axis=1)
    rows = []
    for pos in range(1, 86):
        if pos not in share.index:
            rows.append({"klifs_pos": pos, "klifs_segment": segment(pos), "subpocket": OTHER,
                         "modal_share": np.nan, "runner_up": None, "runner_up_share": np.nan,
                         "n_structures": 0, "median_d_nearest_A": np.nan,
                         "note": "never observed; defaulted to other"})
            continue
        s = share.loc[pos].sort_values(ascending=False)
        sub = ra[ra["klifs_pos"] == pos]
        rows.append({
            "klifs_pos": pos, "klifs_segment": segment(pos), "subpocket": s.index[0],
            "modal_share": round(s.iloc[0], 3),
            "runner_up": s.index[1] if len(s) > 1 and s.iloc[1] > 0 else None,
            "runner_up_share": round(s.iloc[1], 3) if len(s) > 1 else 0.0,
            "n_structures": int(freq.loc[pos].sum()),
            "median_d_nearest_A": round(sub["d_nearest"].median(), 2),
            "note": "unstable (modal share < 0.75)" if s.iloc[0] < 0.75 else "",
        })
    smap = pd.DataFrame(rows)
    smap.to_csv(PROC / "subpocket_map.csv", index=False)
    meta = {"n_main_set": len(st), "n_assigned": int(n_struct), "n_pending_fetch": len(pending),
            "partial": bool(pending)}
    (PROC / "subpocket_map_meta.json").write_text(json.dumps(meta, indent=2))
    if pending:
        log.warning("PARTIAL lookup table: built from %d structures with %d pockets still unfetched",
                    n_struct, len(pending))
    log.info("lookup table sizes:\n%s", smap["subpocket"].value_counts().to_string())
    log.info("unstable positions (modal share < 0.75): %d", (smap["modal_share"] < 0.75).sum())

    with open(RESULTS / "subpocket_map.md", "w") as fh:
        fh.write("# KLIFS position -> KinFragLib subpocket lookup table\n\n")
        fh.write(
            "Subpocket centres: mean C-alpha of KinFragLib anchor residues "
            f"({', '.join(f'{k}: {v}' for k, v in SUBPOCKET_ANCHORS.items())}), from Sydow et al., "
            "*J. Chem. Inf. Model.* 2020, 60, 6081-6094 (doi:10.1021/acs.jcim.0c00839), code at "
            "github.com/volkamerlab/KinaseFocusedFragmentLibrary @88b61d8 (`fragmentation/cli.py` L28-33). "
            f"Each residue is assigned to the centre nearest to any of its heavy atoms, or to `other` if its "
            f"closest heavy atom is >= {OTHER_CUTOFF_A} A from all centres (KinFragLib's X-pool cutoff). The table gives the modal "
            f"label over {n_struct} main-set structures; `modal_share` is that label's frequency. "
            "Segment names are KLIFS's (van Linden et al., *J. Med. Chem.* 2014) and are shown for reference only.\n\n")
        fh.write(smap.to_markdown(index=False))
        fh.write("\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--allow-partial", action="store_true", help="run before the pocket pull finishes (test only)")
    main(ap.parse_args().allow_partial)
