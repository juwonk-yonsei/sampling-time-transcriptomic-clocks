"""tAge clocks (rodent and multispecies; composite and module) as linear weights on
scaled-diff input in mouse Entrez space."""
from __future__ import annotations

from functools import lru_cache

import pandas as pd
from scmodage.clocks import Clock, effective_weights

from .paths import TAGE, TAGE_META

OUTCOME_FILE = {"mortality": "Mortality", "chronological": "Chronoage"}
COMPOSITE_DIR = {"rodent": TAGE / "tage_zenodo", "multispecies": TAGE / "tage_small"}
MODULE_DIR = {"rodent": TAGE / "tage_modules" / "rodent" / "Rodent module clocks",
              "multispecies": TAGE / "tage_modules" / "multispecies" / "Multispecies module clocks"}
SPECIES_LABEL = {"rodent": ("Rodents", "Rodents"), "multispecies": ("Multispecies", "Multispecies")}


def composite(clock_set: str, outcome: str) -> Clock:
    fname = f"EN_{OUTCOME_FILE[outcome]}_{SPECIES_LABEL[clock_set][1]}_Multitissue_scaleddiff.pkl"
    path = COMPOSITE_DIR[clock_set] / fname
    if not path.exists():
        path = TAGE / "tage_zenodo" / fname
    return Clock(name=f"{clock_set}_{outcome}_composite", outcome=outcome, module="composite",
                 function="Composite clock", path=path)


def modules(clock_set: str, outcome: str) -> list[Clock]:
    meta = pd.read_csv(TAGE_META / "module_clocks_metadata.csv")
    meta = meta[(meta["species"] == SPECIES_LABEL[clock_set][0])
                & (meta["outcome"] == ("Mortality" if outcome == "mortality" else "Chronological"))
                & (meta["color"] != "allmodulegenes")]
    return [Clock(name=f"{clock_set}_{outcome}_{r.color}", outcome=outcome, module=r.color, function=r.function,
                  path=MODULE_DIR[clock_set] / r.filename) for r in meta.itertuples()]


@lru_cache(maxsize=None)
def weights(clock: Clock) -> pd.Series:
    """Selected genes (mouse Entrez, str) -> linear weight on scaled-diff input."""
    w = effective_weights(clock)
    w = w[w["selected"] & (w["weight"] != 0)]
    return w.set_index("entrez_mouse")["weight"]


TAGE_OWN_SERIES = {"GSE292884", "GSE292885"}


def training_sources() -> set[str]:
    """GEO series in the tAge training annotation; source 'Current' is tAge's own ITP liver data."""
    from .paths import TAGE_TRAINING
    return set(pd.read_excel(TAGE_TRAINING)["Source"].astype(str)) | TAGE_OWN_SERIES
