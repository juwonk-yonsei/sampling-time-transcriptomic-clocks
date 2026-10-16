#!/usr/bin/env python3
"""Part A. Record instruction-4 decisions, rejudge the gate, define GSE113883 times.

Does not refit models and does not open locked expression.
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = ROOT / "config" / "decisions.tsv"
OUT = ROOT / "results" / "stage4"
OUT.mkdir(parents=True, exist_ok=True)
DATE = "2026-10-16"

ROWS = [
    (
        "4.1.1",
        "Three-platform coverage is approved. The four-platform sentence in instruction 3 section 1.3 was an error. The three coverages in decision 3.3 stand.",
        "Instruction 4 section 1.1 item 1. Report 3 section 11.",
        "prior",
    ),
    (
        "4.1.2",
        "GSE113883 condition=control and period=time_course for every sample is approved.",
        "Instruction 4 section 1.1 item 2. Report 3 section 11.",
        "prior",
    ),
    (
        "4.1.3",
        "The GSE196728 sea-level control label is not applicable. The dataset is excluded by decision 4.3.",
        "Instruction 4 section 1.1 item 3.",
        "prior",
    ),
    (
        "4.1.4",
        "The species-call thresholds are not applicable. The record is kept. GSE196728 is excluded.",
        "Instruction 4 section 1.1 item 4.",
        "prior",
    ),
    (
        "4.1.5",
        "The GSE196728 count-column case match is not applicable. The dataset is excluded.",
        "Instruction 4 section 1.1 item 5.",
        "prior",
    ),
    (
        "4.1.6",
        "Dropping 583 ENSEMBL keys that mapped to more than one SYMBOL is approved.",
        "Instruction 4 section 1.1 item 6. Report 3 section 11.",
        "prior",
    ),
    (
        "4.1.7",
        "Using circadianphase as an hour-unit axis is approved. The stored values are 3-hour bin labels. State that on every result that uses them.",
        "Instruction 4 section 1.1 item 7. Report 3 section 11.",
        "prior",
    ),
    (
        "4.1.8",
        "The unresolved zero point of salivary_dlmo_circadian_time is handled by instruction 4 section 1.4. Stage-3 T_dlmo rows are not used.",
        "Instruction 4 section 1.1 item 8.",
        "prior",
    ),
    (
        "4.1.9",
        "Keeping expression matrices out of git is approved.",
        "Instruction 4 section 1.1 item 9. Report 3 section 11.",
        "prior",
    ),
    (
        "4.1.10",
        "The module-clock path was not found in stage 3. Part B repeats the search with Zenodo checksums.",
        "Instruction 4 section 1.1 item 10.",
        "prior",
    ),
    (
        "4.1.11",
        "GSE113883 sample-median MAD of 0 is not used as the QC flag. Recompute sample medians on genes with mean TPM at least 1.",
        "Instruction 4 section 1.1 item 11 and section 3.1.",
        "prior",
    ),
    (
        "4.1.12",
        "Leaving sections 3.2-3.5 unrun after the stage-3 positive-control failure is approved. After decision 4.2, Part C runs those sections for GSE113883.",
        "Instruction 4 section 1.1 item 12. Report 3 section 12 is also approved.",
        "prior",
    ),
    (
        "4.2",
        "Development positive-control criteria are (a') PER1 and PER3 peaks inside 04-12 inclusive, and PER1 and PER3 LRT p both < 0.01. Criterion (b) is removed. Locked-data gates are set at preregistration and do not use criterion (b). Rejudge uses stored peaks and does not refit. GSE39445 pass, GSE56931 pass, GSE113883 pass.",
        "Instruction 4 section 1.2. Post-hoc correction. GSE113883 clock results had not been calculated.",
        "사후",
    ),
    (
        "4.3",
        "GSE196728 is excluded from every later analysis. Prior outputs are not deleted.",
        "Instruction 4 section 1.3. Sea-level positive control was weak (PER1 LRT p 0.187843896897575, PER3 LRT p 0.0573355762517646, 8 subjects). Processing text maps reads to mm10. ARNTL FPKM is missing for every sample because gene length was unavailable.",
        "사후",
    ),
]


def append_decisions():
    existing = DECISIONS.read_text()
    if not existing.endswith("\n"):
        existing += "\n"
        DECISIONS.write_text(existing)
    present = set()
    for line in existing.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            present.add(parts[1])
    new_lines = []
    for item_id, content, reason, timing in ROWS:
        if item_id in present:
            continue
        new_lines.append(f"{DATE}\t{item_id}\t{content}\t{reason}\t{timing}")
    if new_lines:
        with DECISIONS.open("a") as handle:
            handle.write("\n".join(new_lines) + "\n")
    out = OUT / "partA_decisions.tsv"
    with out.open("w") as handle:
        handle.write("date\tid\tcontent\treason\ttiming\n")
        for item_id, content, reason, timing in ROWS:
            handle.write(f"{DATE}\t{item_id}\t{content}\t{reason}\t{timing}\n")
    print(f"appended {len(new_lines)} decision rows")


def rejudge():
    fits = pd.read_csv(ROOT / "results" / "stage3" / "gates" / "g1_gene_fits.tsv", sep="\t")
    rows = []
    for dataset in ["GSE39445", "GSE56931", "GSE113883"]:
        per1 = fits[(fits.dataset == dataset) & (fits.gene == "PER1")].iloc[0]
        per3 = fits[(fits.dataset == dataset) & (fits.gene == "PER3")].iloc[0]
        a_prime = (4 <= per1.peak_h <= 12) and (4 <= per3.peak_h <= 12)
        p_ok = (per1.lrt_p < 0.01) and (per3.lrt_p < 0.01)
        rows.append(
            {
                "dataset": dataset,
                "per1_peak_h": per1.peak_h,
                "per1_lrt_p": per1.lrt_p,
                "per3_peak_h": per3.peak_h,
                "per3_lrt_p": per3.lrt_p,
                "criterion_a_prime": bool(a_prime),
                "both_lrt_p_lt_0.01": bool(p_ok),
                "criterion_b": "removed",
                "refit": False,
                "pass": bool(a_prime and p_ok),
                "source": "results/stage3/gates/g1_gene_fits.tsv",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "gate_rejudge.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    if not out["pass"].all():
        raise SystemExit("stored positive-control rejudge did not pass")


def parse_clock(text):
    hour, minute = str(text).split(":")
    return int(hour) + int(minute) / 60.0


def dlmo_and_wake():
    samples = pd.read_csv(ROOT / "results" / "stage3" / "GSE113883_samples.tsv", sep="\t")
    dlmo_rows = []
    for subject, group in samples.groupby("subject_id", sort=True):
        values = group["salivary_dlmo_circadian_time"].dropna().unique()
        if len(values) == 0:
            continue
        if len(values) != 1:
            raise SystemExit(f"more than one DLMO value for {subject}")
        value = float(values[0])
        in_range = 18 <= value <= 26
        dlmo_rows.append(
            {
                "subject_id": subject,
                "salivary_dlmo_circadian_time": value,
                "n_samples": int(group["salivary_dlmo_circadian_time"].notna().sum()),
                "n_unique": 1,
                "hour_decimal": value,
                "in_18_26": in_range,
            }
        )
    dlmo = pd.DataFrame(dlmo_rows)
    use_t_dlmo = bool(len(dlmo) == 7 and dlmo["in_18_26"].all())
    dlmo["use_as_dlmo_clock_h"] = use_t_dlmo
    dlmo.to_csv(OUT / "gse113883_dlmo_values.tsv", sep="\t", index=False)
    decision = pd.DataFrame(
        [
            {
                "n_subjects_with_value": int(len(dlmo)),
                "all_in_18_26": bool(dlmo["in_18_26"].all()) if len(dlmo) else False,
                "use_t_dlmo": use_t_dlmo,
                "rule": "numeric hour in [18, 26]; values at or above 24 are next-day morning labels",
                "t_dlmo_definition": "(clock_time_h - dlmo_clock_h) mod 24" if use_t_dlmo else "",
            }
        ]
    )
    decision.to_csv(OUT / "gse113883_dlmo_decision.tsv", sep="\t", index=False)
    print(dlmo.to_string(index=False))
    print(decision.to_string(index=False))

    wake_rows = []
    for subject, group in samples.groupby("subject_id", sort=True):
        group = group.sort_values(["hour", "gsm"])
        times = [parse_clock(x) for x in group["local_time"]]
        elapsed = [0.0]
        for i in range(1, len(times)):
            delta = times[i] - times[i - 1]
            if delta < -1e-9:
                delta += 24.0
            if delta < -1e-9:
                raise SystemExit(f"negative step after wrap for {subject}")
            elapsed.append(elapsed[-1] + delta)
        for (_, row), step in zip(group.iterrows(), elapsed):
            wake = 2.0 + step
            wake_rows.append(
                {
                    "gsm": row.gsm,
                    "subject_id": subject,
                    "local_time": row.local_time,
                    "clock_time_h": row.clock_time_h,
                    "hour": row.hour,
                    "elapsed_from_first_draw_h": step,
                    "hours_since_wake": wake,
                    "hour_minus_hours_since_wake": float(row.hour) - wake,
                }
            )
    wake = pd.DataFrame(wake_rows)
    wake.to_csv(OUT / "gse113883_hours_since_wake.tsv", sep="\t", index=False)
    diff = wake["hour_minus_hours_since_wake"]
    summary = pd.DataFrame(
        [
            {
                "n_samples": int(len(wake)),
                "n_subjects": int(wake.subject_id.nunique()),
                "diff_min": float(diff.min()),
                "diff_median": float(diff.median()),
                "diff_max": float(diff.max()),
                "n_unique_diff": int(diff.round(10).nunique()),
                "definition": "hours_since_wake = 2 + cumulative local_time steps within subject, ordered by hour, with a negative step wrapped by +24",
            }
        ]
    )
    summary.to_csv(OUT / "gse113883_hour_minus_wake_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    append_decisions()
    rejudge()
    dlmo_and_wake()
