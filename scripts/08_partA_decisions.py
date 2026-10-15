#!/usr/bin/env python3
"""Part A. Rejudge the positive control from stage-1 tables. Do not refit.
Check GSE48113 exclusion samples from the metadata table only.
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage2"
OUT.mkdir(parents=True, exist_ok=True)

pc = pd.read_csv(ROOT / "results" / "stage1" / "positive_control.tsv", sep="\t")
decision = pd.read_csv(
    ROOT / "results" / "stage1" / "positive_control_decision.tsv", sep="\t"
)

rows = []
for dataset, block in pc.groupby("dataset", sort=False):
    peaks = {r.gene: r.peak_h for r in block.itertuples(index=False)}
    dist = float(
        decision.loc[decision.dataset == dataset, "per3_arntl_circ_dist_h"].iloc[0]
    )
    a_prime = (4 <= peaks["PER1"] <= 12) and (4 <= peaks["PER3"] <= 12)
    b = dist >= 6
    rows.append(
        {
            "dataset": dataset,
            "PER1_peak_h": peaks["PER1"],
            "PER3_peak_h": peaks["PER3"],
            "criterion_a_prime": a_prime,
            "criterion_b": b,
            "per3_arntl_circ_dist_h": dist,
            "decision": "pass" if a_prime and b else "review",
            "refit": False,
            "source": "results/stage1/positive_control.tsv",
        }
    )
rejudge = pd.DataFrame(rows)
rejudge.to_csv(OUT / "partA_positive_control_rejudge.tsv", sep="\t", index=False)

meta = pd.read_csv(ROOT / "data" / "processed" / "GSE48113_samples.tsv", sep="\t")
in_phase = meta.loc[meta.condition == "control"].copy()
has_clock = in_phase.clock_time_h.notna()
excluded = meta.loc[meta.gsm.isin(["GSM1168837", "GSM1168735"])].copy()
check = pd.DataFrame(
    [
        {
            "n_samples": int(len(meta)),
            "n_in_phase": int((meta.condition == "control").sum()),
            "n_in_phase_with_clock": int(has_clock.sum()),
            "reference_n": 147,
            "matches_reference": bool(int(has_clock.sum()) == 147),
            "n_clock_na_in_phase": int((~has_clock).sum()),
            "gsm1168837_in_table": bool((meta.gsm == "GSM1168837").any()),
            "gsm1168735_in_table": bool((meta.gsm == "GSM1168735").any()),
        }
    ]
)
check.to_csv(OUT / "gse48113_inphase_clock_count.tsv", sep="\t", index=False)
excluded.to_csv(OUT / "gse48113_excluded_samples.tsv", sep="\t", index=False)

gse = pd.read_csv(ROOT / "data" / "processed" / "GSE122725_samples.tsv", sep="\t")


def field(text, key):
    for part in str(text).split("||"):
        part = part.strip()
        prefix = key + ":"
        if part.startswith(prefix):
            return part[len(prefix) :].strip()
    return ""


wake = gse.raw_characteristics.map(lambda s: field(s, "relclocktime"))
wake_num = pd.to_numeric(wake, errors="coerce")
out = pd.DataFrame(
    {
        "gsm": gse.gsm,
        "subject_id": gse.subject_id,
        "relclocktime": wake,
        "wake_ref_h": wake_num,
        "wake_ref_source": "given",
        "clock_time_h": pd.NA,
        "hours_awake": pd.NA,
    }
)
out.to_csv(OUT / "GSE122725_wake_ref.tsv", sep="\t", index=False)
print(rejudge.to_string(index=False))
print(check.to_string(index=False))
print("wake_ref non-null", int(out.wake_ref_h.notna().sum()), "of", len(out))
print("wake_ref min/median/max", out.wake_ref_h.min(), out.wake_ref_h.median(), out.wake_ref_h.max())
