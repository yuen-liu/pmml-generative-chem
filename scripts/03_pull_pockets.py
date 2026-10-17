"""Stage 3 (data pull): pocket coordinates (mol2) for every structure in the main set.

Needed to place KinFragLib subpocket centres in each structure. One file per structure in
data/raw/pockets/<klifs_id>.parquet; cached files are never re-fetched. Rate-limited: each
opencadd call makes two HTTP requests (mol2 + residue-ID table), so calls are spaced
--interval seconds apart (default 1.5 s, i.e. about 1.3 requests/s). Failures go to
data/raw/pockets/_failures.csv and are retried on the next run.
"""
import argparse

import pandas as pd

from common import PROC, RAW, Throttle, klifs_session, log

POCKET_DIR = RAW / "pockets"
KEEP_COLS = ["atom.name", "atom.x", "atom.y", "atom.z", "residue.id", "residue.name", "residue.klifs_id"]


def main(interval, limit):
    POCKET_DIR.mkdir(parents=True, exist_ok=True)
    ids = pd.read_parquet(PROC / "structures_main.parquet")["structure.klifs_id"].astype(int).tolist()
    todo = [i for i in ids if not (POCKET_DIR / f"{i}.parquet").exists()]
    n_cached = len(ids) - len(todo)
    if limit:
        todo = todo[:limit]
    log.info("pockets: %d in main set, %d cached, %d to fetch (~%.0f min at %.1fs/call)",
             len(ids), n_cached, len(todo),
             len(todo) * interval / 60, interval)
    if not todo:
        return
    s = klifs_session()
    throttle = Throttle(min_interval_s=interval)
    failures = []
    for n, sid in enumerate(todo, 1):
        try:
            df = throttle.call(s.coordinates.to_dataframe, sid, entity="pocket", extension="mol2")
            df[KEEP_COLS].to_parquet(POCKET_DIR / f"{sid}.parquet", index=False)
        except Exception as exc:
            failures.append({"structure.klifs_id": sid, "error": f"{type(exc).__name__}: {exc}"[:300]})
            log.warning("pocket %d failed: %s", sid, exc)
        if n % 100 == 0 or n == len(todo):
            log.info("pockets: %d/%d fetched this run (%d failures)", n, len(todo), len(failures))
    pd.DataFrame(failures, columns=["structure.klifs_id", "error"]).to_csv(POCKET_DIR / "_failures.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--interval", type=float, default=1.5, help="seconds between opencadd calls")
    ap.add_argument("--limit", type=int, default=0, help="fetch at most N (for testing)")
    a = ap.parse_args()
    main(a.interval, a.limit)
