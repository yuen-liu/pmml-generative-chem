"""Shared paths, KLIFS caching/rate limiting, subpocket definitions and exclusion logging."""
import logging
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"            # cached API responses (never edited)
PROC = ROOT / "data" / "processed"     # stage outputs
RESULTS = ROOT / "results"             # model outputs and tables
FIGS = ROOT / "figures"
for _d in (RAW, PROC, RESULTS, FIGS):
    _d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("pmml")

N_RESIDUES = 85
# KLIFS interaction.id 1..7, in IFP bit order within each residue (residue-major: bit = 7*(res-1) + (type-1))
INTERACTION_TYPES = [
    "apolar",
    "aromatic_face_to_face",
    "aromatic_edge_to_face",
    "hbond_donor_protein",
    "hbond_acceptor_protein",
    "cation_protein_anion_ligand",
    "anion_protein_cation_ligand",
]
INTERACTION_LABELS = {
    "apolar": "Apolar contact",
    "aromatic_face_to_face": "Aromatic face-to-face",
    "aromatic_edge_to_face": "Aromatic edge-to-face",
    "hbond_donor_protein": "H-bond, protein donor",
    "hbond_acceptor_protein": "H-bond, protein acceptor",
    "cation_protein_anion_ligand": "Protein cation : ligand anion",
    "anion_protein_cation_ligand": "Protein anion : ligand cation",
}

# KinFragLib subpockets. Each subpocket is a point: the mean C-alpha position of its anchor residues
# (KLIFS numbering), recomputed for every structure.
# Source: Sydow, Schmiel, Mortier & Volkamer, J. Chem. Inf. Model. 2020, 60, 6081-6094,
#   doi:10.1021/acs.jcim.0c00839; code: github.com/volkamerlab/KinaseFocusedFragmentLibrary
#   @88b61d8, kinase_focused_fragment_library/fragmentation/cli.py lines 28-33.
SUBPOCKET_ANCHORS = {
    "AP": [46, 51, 75, 15],   # adenine pocket
    "FP": [72, 51, 10, 81],   # front pocket
    "SE": [51],               # solvent-exposed
    "GA": [45, 17, 81],       # gate area
    "B1": [81, 28, 43, 38],   # back pocket 1
    "B2": [18, 24, 70, 83],   # back pocket 2
}
SUBPOCKETS = list(SUBPOCKET_ANCHORS)
OTHER = "other"
REGIONS = SUBPOCKETS + [OTHER]
# KinFragLib sends a fragment to the "X" bin pool when it is >= 8 A from the nearest subpocket
# centre (fragmentation/cli.py line 158). We use the same cutoff for residues -> "other".
OTHER_CUTOFF_A = 8.0
# Nearest-vs-second-nearest centre margin below which an assignment is logged as a near-tie.
NEAR_TIE_MARGIN_A = 1.0

ALPHAS = {"uniform": 1.0, "jeffreys": 0.5}
PRIOR_DOMINATED_N = 10

# Phase 2 fragment categories: mutually exclusive by priority, first match wins.
# See 09_fragment_assign.py for the definitions and the ordering's rationale.
FRAGMENT_CATEGORIES = ["charged", "aromatic", "hbond", "aliphatic"]
FRAGMENT_LABELS = {
    "charged": "Charged",
    "aromatic": "Aromatic",
    "hbond": "Polar, non-aromatic",
    "aliphatic": "Aliphatic",
}

# The two analyses share the Dirichlet-Categorical machinery (06, 06b, 07) and differ only in
# their observation unit and category vocabulary. `--phase` selects one.
PHASES = {
    "ifp": dict(vocab=INTERACTION_TYPES, labels=INTERACTION_LABELS,
                counts="counts_long.csv", by_kinase="counts_by_kinase_long.csv",
                suffix="", unit="IFP interaction"),
    "fragment": dict(vocab=FRAGMENT_CATEGORIES, labels=FRAGMENT_LABELS,
                     counts="counts_fragments_long.csv", by_kinase="counts_fragments_by_kinase_long.csv",
                     suffix="_fragments", unit="BRICS fragment"),
}


class Throttle:
    """Minimum wall-clock spacing between calls to klifs.net (a public academic server)."""

    def __init__(self, min_interval_s):
        self.min_interval_s = min_interval_s
        self._last = 0.0

    def call(self, fn, *args, retries=4, **kwargs):
        for attempt in range(retries):
            wait = self.min_interval_s - (time.monotonic() - self._last)
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()
            try:
                return fn(*args, **kwargs)
            except Exception as exc:  # network / server errors: back off and retry
                if attempt == retries - 1:
                    raise
                backoff = self.min_interval_s * 4 ** (attempt + 1)
                log.warning("KLIFS call failed (%s: %s); retrying in %.0fs", type(exc).__name__, exc, backoff)
                time.sleep(backoff)


def klifs_session():
    from opencadd.databases.klifs import setup_remote
    return setup_remote()


def cached(path, fetch, refresh=False):
    """Return the parquet at `path`, fetching and writing it first if absent (or if refresh)."""
    path = Path(path)
    if path.exists() and not refresh:
        return pd.read_parquet(path)
    df = fetch()
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    log.info("cached %s (%d rows)", path.relative_to(ROOT), len(df))
    return df


def write_exclusions(rows, stage):
    """Write one stage's exclusion log. Every dropped structure gets a row with its rule and reason."""
    cols = ["stage", "structure.klifs_id", "structure.pdb_id", "kinase.klifs_name", "rule", "reason"]
    df = pd.DataFrame(rows, columns=cols[1:])
    df.insert(0, "stage", stage)
    out = PROC / f"exclusions_{stage}.csv"
    df.to_csv(out, index=False)
    log.info("%s: logged %d exclusions -> %s", stage, len(df), out.relative_to(ROOT))
    if len(df):
        log.info("\n%s", df.groupby("rule").size().to_string())
    return df
