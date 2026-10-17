"""Stage 1 (data pull): kinase, structure and ligand metadata plus IFPs for all human structures.

Everything is cached under data/raw/; reruns read from disk. --refresh re-fetches.
About 30 requests total (IFPs in batches of 600 IDs; opencadd's all_interactions() fails with HTTP 414).
"""
import argparse

import pandas as pd

from common import RAW, Throttle, cached, klifs_session, log

IFP_BATCH = 600


def main(refresh):
    s = klifs_session()
    throttle = Throttle(min_interval_s=1.0)

    kin = cached(RAW / "kinases_all.parquet", lambda: throttle.call(s.kinases.all_kinases), refresh)
    # all_kinases() lacks group/family; by_kinase_klifs_id() has them.
    ids = kin["kinase.klifs_id"].astype(int).tolist()
    kin_detail = cached(
        RAW / "kinases_detail.parquet",
        lambda: pd.concat([throttle.call(s.kinases.by_kinase_klifs_id, ids[i:i + 200])
                           for i in range(0, len(ids), 200)], ignore_index=True),
        refresh,
    )
    st = cached(RAW / "structures_all.parquet", lambda: throttle.call(s.structures.all_structures), refresh)
    cached(RAW / "ligands_all.parquet", lambda: throttle.call(s.ligands.all_ligands), refresh)

    human_ids = st.loc[st["species.klifs"] == "Human", "structure.klifs_id"].astype(int).tolist()

    def fetch_ifps():
        parts = []
        for i in range(0, len(human_ids), IFP_BATCH):
            parts.append(throttle.call(s.interactions.by_structure_klifs_id, human_ids[i:i + IFP_BATCH]))
            log.info("IFPs: %d/%d structure IDs requested", min(i + IFP_BATCH, len(human_ids)), len(human_ids))
        return pd.concat(parts, ignore_index=True)

    ifp = cached(RAW / "ifps_human.parquet", fetch_ifps, refresh)
    log.info("kinases %d (detail %d), structures %d, human structures %d, IFPs returned %d",
             len(kin), len(kin_detail), len(st), len(human_ids), len(ifp))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true", help="ignore the cache and re-fetch")
    main(ap.parse_args().refresh)
