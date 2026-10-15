#!/usr/bin/env python3
"""Build sample metadata tables from GEO header files only.

Header files were produced by 01_extract_geo_headers.py, which stops at
!series_matrix_table_begin. This script does not open series matrices.
"""

from __future__ import annotations

import csv
import math
import re
import statistics
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER_DIR = ROOT / "data" / "processed" / "geo_headers"
OUT_DIR = ROOT / "data" / "processed"
RESULT_DIR = ROOT / "results" / "stage1"

COLUMNS = [
    "dataset",
    "gsm",
    "subject_id_raw",
    "subject_id",
    "condition_raw",
    "condition",
    "period",
    "light_group",
    "sample_type",
    "platform",
    "channel_count",
    "clock_time_h",
    "elapsed_h",
    "hours_awake",
    "circadian_ref_h",
    "age",
    "sex",
    "batch_fields",
    "time_source",
    "raw_characteristics",
]


def parse_header(path: Path) -> dict[str, list[list[str]]]:
    rows: dict[str, list[list[str]]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("!"):
            continue
        if line.startswith("!series_matrix_table_begin"):
            break
        parts = line.split("\t")
        key = parts[0][1:]
        vals = [v.strip().strip('"') for v in parts[1:]]
        rows.setdefault(key, []).append(vals)
    return rows


def characteristics_by_sample(rows: dict, n: int) -> list[dict[str, str]]:
    char_rows = rows.get("Sample_characteristics_ch1", [])
    out = []
    for i in range(n):
        acc: dict[str, list[str]] = {}
        order: list[str] = []
        for row in char_rows:
            if i >= len(row):
                continue
            val = row[i].strip()
            if not val:
                continue
            if ":" in val:
                key, value = val.split(":", 1)
                key = key.strip()
                value = value.strip()
            else:
                key, value = val, ""
            acc.setdefault(key, []).append(value)
            order.append(f"{key}: {value}" if value else key)
        flat = {}
        for key, values in acc.items():
            uniq = list(dict.fromkeys(values))
            flat[key] = uniq[0] if len(uniq) == 1 else " | ".join(uniq)
        flat["__raw__"] = " || ".join(order)
        out.append(flat)
    return out


def parse_clock(text: str) -> float | None:
    text = text.strip()
    if not text or text.upper() in {"NA", "N/A", "NULL", "NONE"}:
        return None
    match = re.fullmatch(r"(\d{1,2}):(\d{2})", text)
    if not match:
        return None
    hour = int(match.group(1))
    minute = int(match.group(2))
    if hour > 24 or minute > 59 or (hour == 24 and minute != 0):
        return None
    return hour + minute / 60.0


def is_blank_token(text: str) -> bool:
    return text.strip().upper() in {"", "NA", "N/A", "NULL", "NONE"}


def to_float(text: str) -> float | None:
    text = text.strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def fmt(value) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.10g}"
    return str(value)


def subject_from_39445_title(title: str) -> str | None:
    match = re.search(r"(?:Subject_)?(AF\d+)", title)
    return match.group(1) if match else None


def assign_elapsed_from_numeric(records: list[dict], group_keys: list[str], value_key: str) -> None:
    groups: dict[tuple, list[dict]] = {}
    for rec in records:
        if rec[value_key] is None:
            rec["elapsed_h"] = None
            continue
        key = tuple(rec[k] for k in group_keys)
        groups.setdefault(key, []).append(rec)
    for group in groups.values():
        origin = min(rec[value_key] for rec in group)
        for rec in group:
            rec["elapsed_h"] = rec[value_key] - origin


def assign_elapsed_unwrapped(records: list[dict], group_keys: list[str], order_key: str) -> None:
    groups: dict[tuple, list[dict]] = {}
    for rec in records:
        key = tuple(rec[k] for k in group_keys)
        groups.setdefault(key, []).append(rec)
    for group in groups.values():
        ordered = sorted(group, key=lambda rec: (rec[order_key] is None, rec[order_key] if rec[order_key] is not None else 0))
        elapsed = 0.0
        prev_clock = None
        for rec in ordered:
            clock = rec["clock_time_h"]
            if clock is None or rec[order_key] is None:
                rec["elapsed_h"] = None
                prev_clock = None
                continue
            if prev_clock is None:
                rec["elapsed_h"] = 0.0
            else:
                elapsed += (clock - prev_clock) % 24.0
                rec["elapsed_h"] = elapsed
            prev_clock = clock


def build_gse39445(rows: dict) -> tuple[list[dict], dict]:
    n = len(rows["Sample_geo_accession"][0])
    chars = characteristics_by_sample(rows, n)
    titles = rows["Sample_title"][0]
    notes = {"subject_from_title": 0, "title_subject_mismatch": 0, "clock_parse_fail": 0}
    records = []
    for i in range(n):
        ch = chars[i]
        title = titles[i]
        raw_subject = ch.get("subject", "")
        if raw_subject:
            subject_raw = raw_subject
        else:
            subject_raw = subject_from_39445_title(title) or ""
            if subject_raw:
                notes["subject_from_title"] += 1
        title_subject = subject_from_39445_title(title)
        if raw_subject and title_subject and raw_subject != title_subject:
            notes["title_subject_mismatch"] += 1
        protocol = ch.get("sleepprotocol", "")
        if protocol == "Sleep Extension":
            condition = "control"
        elif protocol == "Sleep Restriction":
            condition = "sleep_restriction"
        else:
            condition = "other"
        clock = parse_clock(ch.get("timesampletaken", ""))
        if not is_blank_token(ch.get("timesampletaken", "")) and clock is None:
            notes["clock_parse_fail"] += 1
        hours = to_float(ch.get("hoursawake", ""))
        phase = to_float(ch.get("circadianphase", ""))
        records.append(
            {
                "dataset": "GSE39445",
                "gsm": rows["Sample_geo_accession"][0][i],
                "subject_id_raw": subject_raw,
                "subject_id": subject_raw,
                "condition_raw": protocol,
                "condition": condition,
                "period": "extended_wakefulness",
                "light_group": "NA",
                "sample_type": ch.get("tissue", ""),
                "platform": rows["Sample_platform_id"][0][i],
                "channel_count": rows["Sample_channel_count"][0][i],
                "clock_time_h": clock,
                "elapsed_h": None,
                "hours_awake": hours,
                "circadian_ref_h": phase,
                "age": "NA",
                "sex": "NA",
                "batch_fields": "NA",
                "time_source": "given" if clock is not None else "NA",
                "raw_characteristics": ch.get("__raw__", ""),
                "_order": hours,
            }
        )
    assign_elapsed_from_numeric(records, ["subject_id", "condition"], "hours_awake")
    return records, notes


def build_gse56931(rows: dict) -> tuple[list[dict], dict]:
    n = len(rows["Sample_geo_accession"][0])
    chars = characteristics_by_sample(rows, n)
    notes = {"clock_parse_fail": 0}
    records = []
    for i in range(n):
        ch = chars[i]
        group = ch.get("biological group", "")
        if group == "baseline":
            condition, period = "control", "baseline"
        elif group == "sleep deprivation":
            condition, period = "sleep_deprivation", "deprivation"
        elif group == "recovery":
            condition, period = "recovery", "recovery"
        else:
            condition, period = "other", "other"
        hour = to_float(ch.get("hour", ""))
        day = to_float(ch.get("day", ""))
        if not is_blank_token(ch.get("hour", "")) and hour is None:
            notes["clock_parse_fail"] += 1
        absolute = None if hour is None or day is None else (day - 1.0) * 24.0 + hour
        records.append(
            {
                "dataset": "GSE56931",
                "gsm": rows["Sample_geo_accession"][0][i],
                "subject_id_raw": ch.get("patient", ""),
                "subject_id": ch.get("patient", ""),
                "condition_raw": group,
                "condition": condition,
                "period": period,
                "light_group": "NA",
                "sample_type": "blood",
                "platform": rows["Sample_platform_id"][0][i],
                "channel_count": rows["Sample_channel_count"][0][i],
                "clock_time_h": hour,
                "elapsed_h": None,
                "hours_awake": None,
                "circadian_ref_h": None,
                "age": "NA",
                "sex": ch.get("Sex", "NA") or "NA",
                "batch_fields": "NA",
                "time_source": "given" if hour is not None else "NA",
                "raw_characteristics": ch.get("__raw__", ""),
                "_absolute_h": absolute,
            }
        )
    assign_elapsed_from_numeric(records, ["subject_id", "period"], "_absolute_h")
    return records, notes


def build_gse48113(rows: dict) -> tuple[list[dict], dict]:
    n = len(rows["Sample_geo_accession"][0])
    chars = characteristics_by_sample(rows, n)
    notes = {"clock_parse_fail": 0}
    records = []
    for i in range(n):
        ch = chars[i]
        timing = ch.get("sleep timing", "")
        if timing == "In phase with respect to melatonin":
            condition = "control"
        elif timing == "Out of phase with respect to melatonin":
            condition = "out_of_phase"
        else:
            condition = "other"
        day = ch.get("day sample was taken in protocol", "")
        if day == "1":
            period = "D1"
        elif day == "4":
            period = "D4"
        else:
            period = "other"
        clock_text = ch.get("time sample taken", "")
        clock = parse_clock(clock_text)
        if not is_blank_token(clock_text) and clock is None:
            notes["clock_parse_fail"] += 1
        point = to_float(ch.get("time point", ""))
        records.append(
            {
                "dataset": "GSE48113",
                "gsm": rows["Sample_geo_accession"][0][i],
                "subject_id_raw": ch.get("subject", ""),
                "subject_id": ch.get("subject", ""),
                "condition_raw": timing,
                "condition": condition,
                "period": period,
                "light_group": "NA",
                "sample_type": ch.get("tissue", ""),
                "platform": rows["Sample_platform_id"][0][i],
                "channel_count": rows["Sample_channel_count"][0][i],
                "clock_time_h": clock,
                "elapsed_h": None,
                "hours_awake": None,
                "circadian_ref_h": None,
                "age": "NA",
                "sex": "NA",
                "batch_fields": "NA",
                "time_source": "given" if clock is not None else "NA",
                "raw_characteristics": ch.get("__raw__", ""),
                "_time_point": point,
            }
        )
    assign_elapsed_unwrapped(records, ["subject_id", "period"], "_time_point")
    return records, notes


def build_gse122725(rows: dict) -> tuple[list[dict], dict]:
    n = len(rows["Sample_geo_accession"][0])
    chars = characteristics_by_sample(rows, n)
    notes = {}
    records = []
    for i in range(n):
        ch = chars[i]
        cond_raw = ch.get("condition", "")
        if cond_raw == "DayOrientedSchedule":
            condition, period = "control", "period1"
        elif cond_raw == "NightOrientedSchedule":
            condition, period = "night_oriented", "period2"
        else:
            condition, period = "other", "other"
        group = ch.get("experimentalgroup", "")
        if group == "BrightLight":
            light = "bright"
        elif group == "Control":
            light = "dim"
        else:
            light = "NA"
        rel = to_float(ch.get("relclocktime", ""))
        records.append(
            {
                "dataset": "GSE122725",
                "gsm": rows["Sample_geo_accession"][0][i],
                "subject_id_raw": ch.get("subject", ""),
                "subject_id": ch.get("subject", ""),
                "condition_raw": cond_raw,
                "condition": condition,
                "period": period,
                "light_group": light,
                "sample_type": "PBMC",
                "platform": rows["Sample_platform_id"][0][i],
                "channel_count": rows["Sample_channel_count"][0][i],
                "clock_time_h": None,
                "elapsed_h": None,
                "hours_awake": None,
                "circadian_ref_h": None,
                "age": "NA",
                "sex": "NA",
                "batch_fields": "NA",
                "time_source": "NA",
                "raw_characteristics": ch.get("__raw__", ""),
                "_rel": rel,
            }
        )
    assign_elapsed_from_numeric(records, ["subject_id", "period"], "_rel")
    return records, notes


def write_table(path: Path, records: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, extrasaction="ignore", delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for rec in records:
            writer.writerow({col: fmt(rec.get(col)) for col in COLUMNS})


def consecutive_intervals(records: list[dict]) -> list[float]:
    groups: dict[tuple, list[dict]] = {}
    for rec in records:
        if rec.get("elapsed_h") is None or not rec.get("subject_id"):
            continue
        key = (rec["subject_id"], rec["period"], rec["condition"])
        groups.setdefault(key, []).append(rec)
    diffs = []
    for group in groups.values():
        ordered = sorted(group, key=lambda rec: rec["elapsed_h"])
        for left, right in zip(ordered, ordered[1:]):
            diffs.append(right["elapsed_h"] - left["elapsed_h"])
    return diffs


def dist_summary(values: list[float]) -> str:
    if not values:
        return "n=0"
    ordered = sorted(values)
    def q(p):
        if len(ordered) == 1:
            return ordered[0]
        pos = (len(ordered) - 1) * p
        lo = math.floor(pos)
        hi = math.ceil(pos)
        if lo == hi:
            return ordered[lo]
        return ordered[lo] * (hi - pos) + ordered[hi] * (pos - lo)
    return (
        f"n={len(ordered)} min={ordered[0]:.6g} q25={q(0.25):.6g} "
        f"median={q(0.5):.6g} q75={q(0.75):.6g} max={ordered[-1]:.6g}"
    )


def summarize(records: list[dict], notes: dict) -> list[str]:
    dataset = records[0]["dataset"]
    lines = [f"dataset\t{dataset}", f"n_samples\t{len(records)}"]
    subjects = [rec["subject_id"] for rec in records if rec["subject_id"]]
    lines.append(f"n_subjects\t{len(set(subjects))}")
    missing_subject = sum(1 for rec in records if not rec["subject_id"])
    lines.append(f"n_missing_subject\t{missing_subject}")
    per_subject = Counter(subjects)
    counts = sorted(per_subject.values())
    lines.append(
        "samples_per_subject\t"
        + (f"min={counts[0]} median={statistics.median(counts):.6g} max={counts[-1]}" if counts else "NA")
    )
    combo = Counter((rec["condition"], rec["period"], rec["light_group"]) for rec in records)
    lines.append("condition_period_light_counts")
    for key in sorted(combo):
        lines.append(f"\t{key[0]}\t{key[1]}\t{key[2]}\t{combo[key]}")
    clocks = [rec["clock_time_h"] for rec in records if rec["clock_time_h"] is not None]
    if clocks:
        lines.append(f"clock_time_h\tn={len(clocks)} min={min(clocks):.6g} max={max(clocks):.6g}")
    else:
        lines.append("clock_time_h\tNA")
    lines.append("consecutive_elapsed_interval_h\t" + dist_summary(consecutive_intervals(records)))
    for key, value in notes.items():
        lines.append(f"note_{key}\t{value}")
    lines.append("")
    return lines


def main() -> None:
    builders = {
        "GSE39445": build_gse39445,
        "GSE56931": build_gse56931,
        "GSE48113": build_gse48113,
        "GSE122725": build_gse122725,
    }
    summary_lines = []
    for gse, builder in builders.items():
        rows = parse_header(HEADER_DIR / f"{gse}_header.txt")
        records, notes = builder(rows)
        write_table(OUT_DIR / f"{gse}_samples.tsv", records)
        summary_lines.extend(summarize(records, notes))
        print(gse, len(records), "subjects", len({r['subject_id'] for r in records if r['subject_id']}))
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    (RESULT_DIR / "metadata_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
