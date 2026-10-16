#!/usr/bin/env python3
"""Recompute R1 model-M LRT and R3 paired t outside the confirmation code."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
OUT = ROOT / "results" / "stage6" / "rehearsal" / "independent"
OUT.mkdir(parents=True, exist_ok=True)


def lrt_row(frame_path: Path, r_stat: float) -> dict:
    frame = pd.read_csv(frame_path, sep="\t")
    full = smf.mixedlm("score ~ cos_t + sin_t", frame, groups=frame["subject_id"]).fit(reml=False, method="bfgs")
    null = smf.mixedlm("score ~ 1", frame, groups=frame["subject_id"]).fit(reml=False, method="bfgs")
    stat = 2 * (full.llf - null.llf)
    rel = abs(stat - r_stat) / abs(r_stat) if r_stat != 0 else math.nan
    return {"r_lrt_stat": r_stat, "py_lrt_stat": stat, "relative_diff": rel, "optimizer": "bfgs"}


def main() -> None:
    rows = []
    point = pd.read_csv(ROOT / "results/stage6/rehearsal/R1/P_point.tsv", sep="\t")
    for _, row in point.iterrows():
        clock = row["clock"]
        out = lrt_row(ROOT / "results/stage6/rehearsal/R1" / f"frame_P_{clock}.tsv", float(row["lrt_stat"]))
        out.update({"rehearsal": "R1", "clock": clock, "quantity": "lrt_stat"})
        rows.append(out)
        print(clock, out["relative_diff"])
    for label in ("gene_column", "transcript_column"):
        pairs = pd.read_csv(ROOT / f"results/stage6/rehearsal/R3/{label}/paired_deltas.tsv", sep="\t")
        r_tab = pd.read_csv(ROOT / f"results/stage6/rehearsal/R3/{label}/paired_tests.tsv", sep="\t")
        for clock in ["TAms_chrono", "TAms_mort", "RA_default"]:
            block = pairs.loc[pairs["clock"] == clock]
            py = stats.ttest_rel(block["morning"], block["evening"])
            r_stat = float(r_tab.loc[r_tab["clock"] == clock, "t_stat"].iloc[0])
            diff = abs(float(py.statistic) - r_stat)
            rows.append({
                "rehearsal": f"R3_{label}",
                "clock": clock,
                "quantity": "t_stat",
                "r_lrt_stat": r_stat,
                "py_lrt_stat": float(py.statistic),
                "relative_diff": diff,
            })
            print(label, clock, diff)
    pd.DataFrame(rows).to_csv(OUT / "independent_check.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
