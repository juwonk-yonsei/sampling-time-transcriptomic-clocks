#!/usr/bin/env python3
"""Recompute run-5 family-P tests outside the confirmation code."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
OUT = ROOT / "results" / "stage7c" / "independent_check.tsv"


def lrt_row(frame_path: Path, r_stat: float) -> dict:
    frame = pd.read_csv(frame_path, sep="\t")
    full = smf.mixedlm("score ~ cos_t + sin_t", frame, groups=frame["subject_id"]).fit(reml=False, method="bfgs")
    null = smf.mixedlm("score ~ 1", frame, groups=frame["subject_id"]).fit(reml=False, method="bfgs")
    stat = 2 * (full.llf - null.llf)
    rel = abs(stat - r_stat) / abs(r_stat) if r_stat != 0 else math.nan
    return {"r_stat": r_stat, "py_stat": stat, "diff": rel, "criterion": "relative", "pass": bool(rel < 1e-2)}


def main() -> None:
    family = pd.read_csv(ROOT / "results/confirm/family_P.tsv", sep="\t")
    rows = []
    for clock in ["TAms_chrono", "TAms_mort", "RA_default"]:
        r_stat = float(family.loc[(family["dataset"] == "L1") & (family["clock"] == clock), "lrt_stat"].iloc[0])
        out = lrt_row(ROOT / "results/confirm/L1" / f"frame_P_{clock}.tsv", r_stat)
        out.update({"dataset": "L1", "clock": clock, "quantity": "lrt_stat", "optimizer": "bfgs"})
        rows.append(out)
        print("L1", clock, out["diff"])
    pairs = pd.read_csv(ROOT / "results/confirm/L3/paired_deltas.tsv", sep="\t")
    for clock in ["TAms_chrono", "TAms_mort", "RA_default"]:
        block = pairs.loc[pairs["clock"] == clock]
        py = stats.ttest_rel(block["morning"], block["evening"])
        r_stat = float(family.loc[(family["dataset"] == "L3") & (family["clock"] == clock), "t_stat"].iloc[0])
        diff = abs(float(py.statistic) - r_stat)
        rows.append({
            "dataset": "L3",
            "clock": clock,
            "quantity": "t_stat",
            "r_stat": r_stat,
            "py_stat": float(py.statistic),
            "diff": diff,
            "criterion": "absolute",
            "pass": bool(diff < 1e-8),
            "optimizer": "",
        })
        print("L3", clock, diff)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)


if __name__ == "__main__":
    main()
