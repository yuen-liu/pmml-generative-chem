"""Phase 2, stage 8 (data pull): orthosteric ligand coordinates (mol2) for every main-set structure.

Raw mol2 text is cached per structure in data/raw/ligands/<klifs_id>.mol2 and never re-fetched.
Rate-limited (default 1.0 s between calls; klifs.net is a public academic server). Failures go to
data/raw/ligands/_failures.csv and are retried on the next run.
"""
import argparse

import pandas as pd

from common import PROC, RAW, Throttle, klifs_session, log

LIGAND_DIR = RAW / "ligands"


def main(interval, limit):
    LIGAND_DIR.mkdir(parents=True, exist_ok=True)
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
            (LIGAND_DIR / f"{sid}.mol2").write_text(text)
        except Exception as exc:
            failures.append({"structure.klifs_id": sid, "error": f"{type(exc).__name__}: {exc}"[:300]})
            log.warning("ligand %d failed: %s", sid, exc)
        if n % 100 == 0 or n == len(todo):
            log.info("ligands: %d/%d fetched this run (%d failures)", n, len(todo), len(failures))
    pd.DataFrame(failures, columns=["structure.klifs_id", "error"]).to_csv(LIGAND_DIR / "_failures.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--interval", type=float, default=1.0, help="seconds between calls")
    ap.add_argument("--limit", type=int, default=0, help="fetch at most N (for testing)")
    a = ap.parse_args()
    main(a.interval, a.limit)
