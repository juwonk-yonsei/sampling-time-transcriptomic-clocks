"""Loading and applying the published tAge clocks (composite and module clocks).

All clocks are relative (scaled-diff) elastic-net pipelines in mouse Entrez space.
Mortality clocks output log10(hazard ratio); chronological clocks output age as a
fraction of species maximum lifespan.
"""
from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

from .paths import RETRAINED_CLOCKS, TAGE_META, TAGE_MODULES, TAGE_ZENODO

OUTCOMES = {"mortality": "Mortality", "chronological": "Chronoage"}
UNITS = {"mortality": "log10 hazard ratio", "chronological": "fraction of maximum lifespan"}


@dataclass(frozen=True)
class Clock:
    name: str
    outcome: str
    module: str
    function: str
    path: Path


def _patch(model: Pipeline) -> Pipeline:
    for _, step in model.steps:
        if isinstance(step, SimpleImputer) and not hasattr(step, "_fill_dtype"):
            step._fill_dtype = step.statistics_.dtype
    return model


@lru_cache(maxsize=None)
def load_model(path: str) -> Pipeline:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _patch(joblib.load(path))


def module_registry(outcome: str = "mortality", include_all: bool = False, clock_set: str | None = None) -> list[Clock]:
    """clock_set: "published" (tAge, default) or "retrained" (manuscript specification,
    scripts/23_retrain_clocks.py); defaults to the SCMODAGE_CLOCKS environment variable."""
    clock_set = clock_set or os.environ.get("SCMODAGE_CLOCKS", "published")
    meta = pd.read_csv(TAGE_META / "module_clocks_metadata.csv")
    meta = meta[(meta["species"] == "Rodents") & (meta["outcome"] == ("Mortality" if outcome == "mortality" else "Chronological"))]
    if not include_all:
        meta = meta[meta["color"] != "allmodulegenes"]
    root = TAGE_MODULES if clock_set == "published" else RETRAINED_CLOCKS
    return [
        Clock(name=r.color, outcome=outcome, module=r.color, function=r.function, path=root / r.filename)
        for r in meta.itertuples()
    ]


def composite_clock(outcome: str = "mortality") -> Clock:
    fname = f"EN_{OUTCOMES[outcome]}_Rodents_Multitissue_scaleddiff.pkl"
    return Clock(name="composite", outcome=outcome, module="composite", function="Composite clock", path=TAGE_ZENODO / fname)


def module_functions() -> pd.Series:
    m = pd.read_csv(TAGE_META / "modules" / "Module_to_function_map_rodent.csv", index_col=0)
    return m["Function"]


def short_label(module: str) -> str:
    f = module_functions().get(module, module)
    return f if f != "-" else module


def features(clock: Clock) -> list[str]:
    return [str(g) for g in load_model(str(clock.path)).feature_names_in_]


def predict(clock: Clock, scaled_diff: pd.DataFrame) -> pd.Series:
    """`scaled_diff`: samples x mouse-Entrez genes (NaN = undetected; imputed by the clock)."""
    model = load_model(str(clock.path))
    feats = [str(g) for g in model.feature_names_in_]
    X = scaled_diff.reindex(columns=feats)
    X.columns = model.feature_names_in_
    if not np.any(model.steps[-1][1].coef_):
        return pd.Series(np.nan, index=scaled_diff.index, name=clock.name)
    return pd.Series(model.predict(X), index=scaled_diff.index, name=clock.name)


def predict_all(clocks: list[Clock], scaled_diff: pd.DataFrame) -> pd.DataFrame:
    return pd.concat([predict(c, scaled_diff) for c in clocks], axis=1)


def effective_weights(clock: Clock) -> pd.DataFrame:
    """Per-gene linear weight of the clock on scaled-diff input (after centring/selection)."""
    model = load_model(str(clock.path))
    feats = np.array([str(g) for g in model.feature_names_in_])
    steps = dict(model.steps)
    imputer = next(s for s in steps.values() if isinstance(s, SimpleImputer))
    scaler = steps.get("scaler")
    selector = steps.get("featureselection")
    est = model.steps[-1][1]
    keep = selector.get_support() if selector is not None else np.ones(len(feats), bool)
    scale = getattr(scaler, "scale_", None)
    scale = np.ones(len(feats)) if scale is None else scale
    coef = np.zeros(len(feats))
    coef[keep] = est.coef_
    coef = coef / scale
    return pd.DataFrame(
        {"entrez_mouse": feats, "weight": coef, "training_median": imputer.statistics_, "selected": keep}
    )


def detection_fraction(clock: Clock, scaled_diff: pd.DataFrame) -> pd.Series:
    feats = features(clock)
    X = scaled_diff.reindex(columns=feats)
    return X.notna().mean(axis=1)


def mask_low_detection(scores: pd.DataFrame, min_detection: float, by: str = "celltype") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Set a clock's scores to NaN in every stratum where its median feature detection is
    below `min_detection` (near-constant predictions carry no information). Returns the
    masked scores and the per-stratum median detection table."""
    scores = scores.copy()
    score_cols = [c for c in scores.columns if c.startswith(("mortality:", "chronological:"))]
    det_cols = {c: f"detected_{c}" for c in score_cols if f"detected_{c}" in scores.columns}
    med = scores.groupby(by)[list(det_cols.values())].median()
    for col, dcol in det_cols.items():
        low = med.index[med[dcol] < min_detection]
        scores.loc[scores[by].isin(low), col] = np.nan
    return scores, med
