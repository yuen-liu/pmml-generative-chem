"""Phase 2, stage 8 (data pull): orthosteric ligand coordinates (mol2) for every main-set structure.

Raw mol2 text is cached per structure in data/raw/ligands/<klifs_id>.mol2 and never re-fetched.
Rate-limited (default 1.0 s between calls; klifs.net is a public academic server). Failures go to
data/raw/ligands/_failures.csv and are retried on the next run.
"""
import argparse

import pandas as pd

from common import PROC, RAW, Throttle, klifs_session, log

LIGAND_DIR = RAW / "ligands"


def complete_mol2(text):
    """A mol2 record is complete when the ATOM block holds the atom count the header declares."""
    try:
        lines = text.splitlines()
        i = lines.index("@<TRIPOS>MOLECULE")
        n_atoms = int(lines[i + 2].split()[0])
        j = lines.index("@<TRIPOS>ATOM")
    except (ValueError, IndexError):
        return False
    seen = 0
    for line in lines[j + 1:]:
        if line.startswith("@<TRIPOS>"):
            break
        if line.strip():
            seen += 1
    return seen == n_atoms and n_atoms > 0


def verify_cache():
    """Delete cached files that are truncated or unreadable so the next run re-fetches them."""
    bad = []
    for f in sorted(LIGAND_DIR.glob("*.mol2")):
        try:
            ok = complete_mol2(f.read_text())
        except OSError:
            ok = False
        if not ok:
            bad.append(f.name)
            f.unlink()
    for tmp in LIGAND_DIR.glob(".*.mol2.tmp"):
        tmp.unlink()
    if bad:
        log.warning("removed %d incomplete cached ligand(s), will re-fetch: %s", len(bad), ", ".join(bad))
    else:
        log.info("cache verified: every cached ligand file is complete")
    return bad


def main(interval, limit):
    LIGAND_DIR.mkdir(parents=True, exist_ok=True)
    verify_cache()
    ids = pd.read_parquet(PROC / "structures_main.parquet")["structure.klifs_id"].astype(int).tolist()
    todo = [i for i in ids if not (LIGAND_DIR / f"{i}.mol2").exists()]
    n_cached = len(ids) - len(todo)
    if limit:
        todo = todo[:limit]
    log.info("ligands: %d in main set, %d cached, %d to fetch (~%.0f min at %.1fs/call)",
             len(ids), n_cached, len(todo), len(todo) * interval / 60, interval)
    if not todo:
        return
    s = klifs_session()
    throttle = Throttle(min_interval_s=interval)
    failures = []
    for n, sid in enumerate(todo, 1):
        try:
            text = throttle.call(s.coordinates.to_text, sid, entity="ligand", extension="mol2")
            if not complete_mol2(text):
                raise ValueError("response is not a complete mol2 record")
            # Write atomically: a kill mid-write can never leave a partial file that the next
            # run would skip as already cached.
            tmp = LIGAND_DIR / f".{sid}.mol2.tmp"
            tmp.write_text(text)
            tmp.replace(LIGAND_DIR / f"{sid}.mol2")
        except Exception as exc:
            failures.append({"structure.klifs_id": sid, "error": f"{type(exc).__name__}: {exc}"[:300]})
            log.warning("ligand %d failed: %s", sid, exc)
        if n % 100 == 0 or n == len(todo):
            log.info("ligands: %d/%d fetched this run (%d failures)", n, len(todo), len(failures))
    pd.DataFrame(failures, columns=["structure.klifs_id", "error"]).to_csv(LIGAND_DIR / "_failures.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--interval", type=float, default=1.0, help="seconds between calls")
    ap.add_argument("--limit", type=int, default=0, help="fetch at most N (0 = no limit)")
    a = ap.parse_args()
    main(a.interval, a.limit)
