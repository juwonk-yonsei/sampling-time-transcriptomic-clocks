#!/usr/bin/env python3
"""Part B metadata, circadian-phase check, species check, checksums.

Does not open GSE295758 expression. Stops at !series_matrix_table_begin.
Does not download module-clock files.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage3"
OUT.mkdir(parents=True, exist_ok=True)
RAW = ROOT / "data" / "raw"
LOCKED = ROOT / "data" / "locked"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def append_checksums(rows: list[tuple[str, Path]]) -> None:
    path = ROOT / "config" / "raw_checksums.tsv"
    text = path.read_text()
    present = set()
    for line in text.splitlines()[1:]:
        parts = line.split("\t")
        if parts:
            present.add(parts[0])
    new = []
    for rel, file_path in rows:
        if rel in present:
            continue
        new.append(f"{rel}\t{file_path.stat().st_size}\t{sha256(file_path)}")
    if new:
        if not text.endswith("\n"):
            text += "\n"
        path.write_text(text + "\n".join(new) + "\n")
    print(f"checksums appended {len(new)}")


def read_bang_lines(path: Path) -> list[str]:
    lines = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                if line.strip() == "":
                    continue
                break
            lines.append(line.rstrip("\n"))
    return lines


def split_geo(line: str) -> tuple[str, list[str]]:
    key, _, rest = line.partition("\t")
    vals = [piece.strip().strip('"') for piece in rest.split("\t")] if rest else []
    return key, vals


def sample_frame(lines: list[str]) -> tuple[pd.DataFrame, list[str]]:
    rows: dict[str, list[list[str]]] = {}
    series_lines = []
    n = None
    for line in lines:
        key, vals = split_geo(line)
        if key.startswith("!Series_"):
            series_lines.append(line)
            continue
        if not key.startswith("!Sample_"):
            continue
        if n is None:
            n = len(vals)
        if key == "!Sample_characteristics_ch1":
            rows.setdefault(key, []).append(vals)
        else:
            rows.setdefault(key, []).append(vals)
    if n is None:
        raise SystemExit("no sample rows")
    base = {key[1:]: vals[0] for key, vals in rows.items() if key != "!Sample_characteristics_ch1"}
    char_rows = rows.get("!Sample_characteristics_ch1", [])
    records = []
    for i in range(n):
        char = {}
        raw_parts = []
        for vals in char_rows:
            text = vals[i] if i < len(vals) else ""
            raw_parts.append(text)
            if ":" in text:
                name, value = text.split(":", 1)
                char[name.strip()] = value.strip()
        rec = {key: (vals[i] if i < len(vals) else "") for key, vals in base.items()}
        rec["characteristics"] = char
        rec["raw_characteristics"] = " || ".join(raw_parts)
        records.append(rec)
    return pd.DataFrame(records), series_lines


def clock_from_hhmm(text: str) -> float:
    text = str(text).strip()
    if text.upper() in {"", "NA", "NAN", "NONE"}:
        return np.nan
    if ":" not in text:
        return float(text)
    hour, minute = text.split(":", 1)
    return float(hour) + float(minute) / 60.0


def to_num(text: str) -> float:
    text = str(text).strip()
    if text.upper() in {"", "NA", "NAN", "NONE"}:
        return np.nan
    return float(text)


def circ_mean_h(hours: np.ndarray) -> float:
    hours = hours[np.isfinite(hours)]
    if len(hours) == 0:
        return np.nan
    ang = hours / 24.0 * 2.0 * np.pi
    return float((np.arctan2(np.mean(np.sin(ang)), np.mean(np.cos(ang))) % (2.0 * np.pi)) * 24.0 / (2.0 * np.pi))


def circ_sd_h(hours: np.ndarray) -> float:
    hours = hours[np.isfinite(hours)]
    if len(hours) < 2:
        return np.nan
    ang = hours / 24.0 * 2.0 * np.pi
    r_bar = float(np.hypot(np.mean(np.cos(ang)), np.mean(np.sin(ang))))
    if r_bar <= 0:
        return np.nan
    if r_bar >= 1:
        return 0.0
    return float(np.sqrt(-2.0 * np.log(r_bar)) * 24.0 / (2.0 * np.pi))


def mod24(values: np.ndarray) -> np.ndarray:
    return np.mod(values, 24.0)


def circadianphase() -> None:
    meta = pd.read_csv(ROOT / "data" / "processed" / "GSE39445_samples.tsv", sep="\t")
    meta["clock_time_h"] = pd.to_numeric(meta["clock_time_h"], errors="coerce")
    meta["hours_awake"] = pd.to_numeric(meta["hours_awake"], errors="coerce")
    meta["circadian_ref_h"] = pd.to_numeric(meta["circadian_ref_h"], errors="coerce")
    rows = []
    for (condition, subject), block in meta.groupby(["condition", "subject_id"], sort=False):
        clock_off = mod24(block["clock_time_h"].to_numpy() - block["circadian_ref_h"].to_numpy())
        wake_off = mod24(block["hours_awake"].to_numpy() - block["circadian_ref_h"].to_numpy())
        raw_wake = block["hours_awake"].to_numpy() - block["circadian_ref_h"].to_numpy()
        rows.append(
            {
                "condition": condition,
                "subject_id": subject,
                "n_samples": int(len(block)),
                "circ_sd_clock_minus_phase_h": circ_sd_h(clock_off),
                "circ_mean_clock_minus_phase_h": circ_mean_h(clock_off),
                "circ_sd_awake_minus_phase_h": circ_sd_h(wake_off),
                "circ_mean_awake_minus_phase_h": circ_mean_h(wake_off),
                "linear_sd_awake_minus_phase_h": float(np.nanstd(raw_wake, ddof=1)) if len(block) > 1 else np.nan,
                "phase_step_unique": ",".join(str(x) for x in sorted(block["circadian_ref_h"].dropna().unique())),
            }
        )
    per = pd.DataFrame(rows)
    per.to_csv(OUT / "gse39445_circadianphase_offset_by_subject.tsv", sep="\t", index=False, float_format="%.10g")
    summary = []
    for condition, block in per.groupby("condition", sort=False):
        means = block["circ_mean_clock_minus_phase_h"].to_numpy()
        wake_means = block["circ_mean_awake_minus_phase_h"].to_numpy()
        summary.append(
            {
                "condition": condition,
                "n_subjects": int(block["subject_id"].nunique()),
                "clock_offset_circ_sd_min": float(np.nanmin(block["circ_sd_clock_minus_phase_h"])),
                "clock_offset_circ_sd_median": float(np.nanmedian(block["circ_sd_clock_minus_phase_h"])),
                "clock_offset_circ_sd_max": float(np.nanmax(block["circ_sd_clock_minus_phase_h"])),
                "clock_offset_circ_mean_min": float(np.nanmin(means)),
                "clock_offset_circ_mean_median": float(np.nanmedian(means)),
                "clock_offset_circ_mean_max": float(np.nanmax(means)),
                "awake_offset_circ_sd_min": float(np.nanmin(block["circ_sd_awake_minus_phase_h"])),
                "awake_offset_circ_sd_median": float(np.nanmedian(block["circ_sd_awake_minus_phase_h"])),
                "awake_offset_circ_sd_max": float(np.nanmax(block["circ_sd_awake_minus_phase_h"])),
                "awake_offset_circ_mean_min": float(np.nanmin(wake_means)),
                "awake_offset_circ_mean_median": float(np.nanmedian(wake_means)),
                "awake_offset_circ_mean_max": float(np.nanmax(wake_means)),
                "awake_linear_sd_median": float(np.nanmedian(block["linear_sd_awake_minus_phase_h"])),
                "awake_linear_sd_max": float(np.nanmax(block["linear_sd_awake_minus_phase_h"])),
            }
        )
    pd.DataFrame(summary).to_csv(
        OUT / "gse39445_circadianphase_offset_summary.tsv", sep="\t", index=False, float_format="%.10g"
    )
    values = meta["circadian_ref_h"].dropna().to_numpy()
    step = np.unique(np.round(np.diff(np.sort(np.unique(values))), 10))
    multiple_of_3 = bool(np.all(np.isclose(values / 3.0, np.round(values / 3.0))))
    awake_sd = per["linear_sd_awake_minus_phase_h"].to_numpy()
    n_nonzero_awake_sd = int(np.sum(awake_sd > 1e-6))
    # Hour scale: values are multiples of 3, and within-subject awake-minus-phase
    # SD is 0 for every subject except those already missing a regular step.
    data_matches_hours = multiple_of_3 and float(np.nanmedian(awake_sd)) == 0.0
    pd.DataFrame(
        [
            {
                "n_values": int(len(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "unique_sorted": ",".join(str(x) for x in np.sort(np.unique(values))),
                "all_multiples_of_3": multiple_of_3,
                "unique_sorted_differences": ",".join(str(x) for x in step),
                "max_within_subject_linear_sd_awake_minus_phase": float(np.nanmax(awake_sd)),
                "n_subjects_awake_linear_sd_gt_0": n_nonzero_awake_sd,
                "data_matches_hour_scale": data_matches_hours,
                "use_as_time_axis_in_section_5_1": data_matches_hours,
                "zero_point_from_si": "melatonin peak; tp = 24*t - 24*p with t and p in fractional days; tp is hours away from the melatonin peak",
                "binning_from_si": "3-h-wide circadian bins, also described as 45 degree bins",
            }
        ]
    ).to_csv(OUT / "gse39445_circadianphase_decision.tsv", sep="\t", index=False)
    print("circadianphase data_matches_hours", data_matches_hours, "n_nonzero", n_nonzero_awake_sd)


def build_gse113883() -> None:
    lines = read_bang_lines(RAW / "GSE113883" / "GSE113883_series_matrix.txt.gz")
    (OUT / "GSE113883_bang_lines.txt").write_text("\n".join(lines) + "\n")
    frame, _ = sample_frame(lines)
    char = frame["characteristics"]
    hour = char.map(lambda d: to_num(d.get("hour", "")))
    local = char.map(lambda d: d.get("local_time", ""))
    clock = local.map(clock_from_hhmm)
    dlmo = char.map(lambda d: to_num(d.get("salivary_dlmo_circadian_time", "")))
    subject = char.map(lambda d: d.get("subjectid", ""))
    elapsed = hour.groupby(subject).transform(lambda s: s - s.min())
    out = pd.DataFrame(
        {
            "dataset": "GSE113883",
            "gsm": frame["Sample_geo_accession"],
            "sample_title": frame["Sample_title"],
            "subject_id_raw": subject,
            "subject_id": subject,
            "condition_raw": "",
            "condition": "control",
            "period": "time_course",
            "light_group": "NA",
            "sample_type": char.map(lambda d: d.get("tissue", "")),
            "platform": frame["Sample_platform_id"],
            "channel_count": frame["Sample_channel_count"],
            "clock_time_h": clock,
            "elapsed_h": elapsed,
            "hours_awake": np.nan,
            "circadian_ref_h": dlmo,
            "age": np.nan,
            "sex": "NA",
            "batch_fields": "NA",
            "time_source": "given",
            "local_time": local,
            "hour": hour,
            "salivary_dlmo_circadian_time": char.map(lambda d: d.get("salivary_dlmo_circadian_time", "")),
            "salivary_dlmo_zt": char.map(lambda d: d.get("salivary_dlmo_zt", "")),
            "serum_melatonin": char.map(lambda d: d.get("serum_melatonin", "")),
            "raw_characteristics": frame["raw_characteristics"],
        }
    )
    out.to_csv(OUT / "GSE113883_samples.tsv", sep="\t", index=False, float_format="%.10g")
    protocols = []
    for line in lines:
        if line.startswith("!Sample_extract_protocol_ch1") or line.startswith("!Sample_data_processing") or line.startswith("!Series_overall_design") or line.startswith("!Series_summary"):
            _, vals = split_geo(line)
            text = vals[0] if vals else ""
            if text not in protocols:
                protocols.append(line.split("\t", 1)[0] + "\t" + text)
    (OUT / "GSE113883_protocol_lines.txt").write_text("\n".join(protocols) + "\n")
    print("GSE113883 samples", len(out), "dlmo non-null", int(out["circadian_ref_h"].notna().sum()))


def classify_gene_id(gene_id: str) -> str:
    text = gene_id.strip()
    if text.startswith("ENSMUSG"):
        return "mouse_ensembl"
    if text.startswith("ENSG"):
        return "human_ensembl"
    letters = [ch for ch in text if ch.isalpha()]
    if not letters:
        return "other"
    if letters[0].isupper() and all(ch.islower() for ch in letters[1:]):
        return "mouse_symbol_style"
    if all(ch.isupper() for ch in letters):
        return "human_symbol_style"
    return "mixed_symbol"


def build_gse196728() -> None:
    lines = read_bang_lines(RAW / "GSE196728" / "GSE196728_series_matrix.txt.gz")
    (OUT / "GSE196728_bang_lines.txt").write_text("\n".join(lines) + "\n")
    frame, _ = sample_frame(lines)
    char = frame["characteristics"]
    subject = char.map(lambda d: d.get("subject", ""))
    altitude = char.map(lambda d: d.get("altitude", ""))
    clock = char.map(lambda d: to_num(d.get("time", "")))
    out = pd.DataFrame(
        {
            "dataset": "GSE196728",
            "gsm": frame["Sample_geo_accession"],
            "sample_title": frame["Sample_title"],
            "subject_id_raw": subject,
            "subject_id": subject,
            "condition_raw": altitude,
            "condition": np.where(altitude.eq("SL"), "control", altitude),
            "period": altitude,
            "light_group": "NA",
            "sample_type": "Whole blood",
            "platform": frame["Sample_platform_id"],
            "channel_count": frame["Sample_channel_count"],
            "clock_time_h": clock,
            "elapsed_h": np.nan,
            "hours_awake": np.nan,
            "circadian_ref_h": np.nan,
            "age": np.nan,
            "sex": "NA",
            "batch_fields": "NA",
            "time_source": "given",
            "altitude": altitude,
            "time": char.map(lambda d: d.get("time", "")),
            "raw_characteristics": frame["raw_characteristics"],
        }
    )
    out.to_csv(OUT / "GSE196728_samples.tsv", sep="\t", index=False, float_format="%.10g")
    protocols = []
    for line in lines:
        if line.startswith(("!Sample_data_processing", "!Sample_extract_protocol_ch1", "!Series_overall_design", "!Series_summary", "!Sample_organism_ch1", "!Sample_source_name_ch1")):
            _, vals = split_geo(line)
            text = vals[0] if vals else ""
            if text not in [p.split("\t", 1)[-1] for p in protocols]:
                protocols.append(line.split("\t", 1)[0] + "\t" + text)
    (OUT / "GSE196728_protocol_lines.txt").write_text("\n".join(protocols) + "\n")

    xlsx = RAW / "GSE196728" / "GSE196728_Peru_raw.xlsx"
    book = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    sheet = book[book.sheetnames[0]]
    ids = []
    header = None
    for i, row in enumerate(sheet.iter_rows(values_only=True)):
        if i == 0:
            header = ["" if cell is None else str(cell) for cell in row]
            continue
        gene_id = "" if row[0] is None else str(row[0])
        if gene_id:
            ids.append(gene_id)
    book.close()
    classes = pd.Series([classify_gene_id(gene_id) for gene_id in ids])
    counts = classes.value_counts()
    n = int(len(ids))
    n_mouse_ens = int(counts.get("mouse_ensembl", 0))
    n_mouse_sym = int(counts.get("mouse_symbol_style", 0))
    n_human = int(counts.get("human_ensembl", 0) + counts.get("human_symbol_style", 0))
    n_mixed = int(counts.get("mixed_symbol", 0) + counts.get("other", 0))
    # Human symbols are mostly all capitals. A minority (C#orf#) have a lowercase tail
    # and match the letter rule used for mouse symbols. Require a human-style majority.
    use = n > 0 and n_mouse_ens == 0 and (n_human / n) >= 0.8 and (n_mouse_sym / n) < 0.2
    decision = pd.DataFrame(
        [
            {
                "n_gene_ids": n,
                "n_human_ensembl": int(counts.get("human_ensembl", 0)),
                "n_human_symbol_style": int(counts.get("human_symbol_style", 0)),
                "n_mouse_ensembl": n_mouse_ens,
                "n_mouse_symbol_style": n_mouse_sym,
                "n_mixed_or_other": n_mixed,
                "fraction_human_style": n_human / n if n else np.nan,
                "fraction_mouse_symbol_style": n_mouse_sym / n if n else np.nan,
                "example_10": "|".join(ids[:10]),
                "use_for_analysis": use,
                "rule": "use only if ENSMUSG count is 0, all-cap or ENSG fraction is at least 0.8, and first-cap-rest-lower fraction is below 0.2",
            }
        ]
    )
    decision.to_csv(OUT / "gse196728_species_decision.tsv", sep="\t", index=False, float_format="%.10g")
    print(decision.to_string(index=False))
    if not use:
        print("GSE196728 not used")
        return
    counts_df = pd.read_excel(xlsx, sheet_name=0, index_col=0)
    counts_df.index = counts_df.index.astype(str)
    counts_df.columns = [str(col) for col in counts_df.columns]
    titles = list(out["sample_title"])
    stripped = {}
    for col in counts_df.columns:
        key = col[:-4] if col.endswith(".raw") else col
        stripped.setdefault(key.casefold(), []).append(col)
    matched = []
    unmatched_titles = []
    for title in titles:
        hits = stripped.get(title.casefold(), [])
        if len(hits) == 1:
            matched.append(title)
        else:
            unmatched_titles.append(title)
    used = {stripped[title.casefold()][0] for title in matched}
    unmatched_cols = [col for cols in stripped.values() for col in cols if col not in used]
    pd.DataFrame(
        [
            {
                "n_matrix_columns": int(counts_df.shape[1]),
                "n_samples": int(len(out)),
                "n_matched_titles": int(len(matched)),
                "n_unmatched_titles": int(len(unmatched_titles)),
                "n_unmatched_columns": int(len(unmatched_cols)),
                "unmatched_title_examples": "|".join(unmatched_titles[:10]),
                "unmatched_column_examples": "|".join(unmatched_cols[:10]),
            }
        ]
    ).to_csv(OUT / "gse196728_column_match.tsv", sep="\t", index=False)
    if unmatched_titles or unmatched_cols:
        print("GSE196728 column match incomplete; counts not written")
        return
    ordered = [stripped[title.casefold()][0] for title in out["sample_title"]]
    slim = counts_df.loc[:, ordered]
    slim.columns = out["gsm"].tolist()
    expr_dir = OUT / "expr"
    expr_dir.mkdir(parents=True, exist_ok=True)
    slim.to_csv(expr_dir / "GSE196728_counts.tsv.gz", sep="\t", compression="gzip")
    print("wrote counts", slim.shape)


def build_gse295758() -> None:
    lines = read_bang_lines(LOCKED / "GSE295758" / "GSE295758_series_matrix.txt.gz")
    (OUT / "GSE295758_bang_lines.txt").write_text("\n".join(lines) + "\n")
    frame, _ = sample_frame(lines)
    char = frame["characteristics"]
    keys = sorted({key for item in char for key in item})
    (OUT / "GSE295758_characteristic_keys.txt").write_text("\n".join(keys) + "\n")
    subject = frame["Sample_title"].str.extract(r"(subject\d+)", expand=False)
    # Prefer a subject characteristic when present.
    if any("subject" in item for item in char):
        subject = char.map(lambda d: d.get("subject", d.get("subjectid", "")))
    timepoint = char.map(lambda d: d.get("timepoint", ""))
    age = char.map(lambda d: to_num(d.get("age", "")))
    sex = char.map(lambda d: d.get("Sex", d.get("sex", "")))
    out = pd.DataFrame(
        {
            "dataset": "GSE295758",
            "gsm": frame["Sample_geo_accession"],
            "sample_title": frame["Sample_title"],
            "subject_id": subject,
            "timepoint": timepoint,
            "age": age,
            "sex": sex,
            "diagnosis_or_group": char.map(lambda d: d.get("opioids", d.get("diagnosis", d.get("group", "")))),
            "batch": char.map(lambda d: d.get("batch", "")),
            "library_strategy": frame.get("Sample_library_strategy", pd.Series([""] * len(frame))),
            "instrument": frame.get("Sample_instrument_model", pd.Series([""] * len(frame))),
            "raw_characteristics": frame["raw_characteristics"],
        }
    )
    # Keep every characteristic key as its own column.
    for key in keys:
        col = "char_" + key.replace(" ", "_")
        if col not in out.columns:
            out[col] = char.map(lambda d, key=key: d.get(key, ""))
    out.to_csv(OUT / "GSE295758_samples.tsv", sep="\t", index=False, float_format="%.10g")
    both = (
        out.groupby("subject_id")["timepoint"]
        .agg(lambda s: set(s.dropna()))
        .map(lambda values: "D" in values and "N" in values)
    )
    ages = out.drop_duplicates("subject_id")["age"]
    summary = pd.DataFrame(
        [
            {
                "n_samples": int(len(out)),
                "n_subjects": int(out["subject_id"].nunique()),
                "samples_per_subject_min": int(out.groupby("subject_id").size().min()),
                "samples_per_subject_median": float(out.groupby("subject_id").size().median()),
                "samples_per_subject_max": int(out.groupby("subject_id").size().max()),
                "n_subjects_with_D_and_N": int(both.sum()),
                "age_min": float(ages.min()),
                "age_q1": float(ages.quantile(0.25)),
                "age_median": float(ages.median()),
                "age_q3": float(ages.quantile(0.75)),
                "age_max": float(ages.max()),
                "n_sex": int(out.drop_duplicates("subject_id")["sex"].nunique()),
                "sex_counts_subjects": str(out.drop_duplicates("subject_id")["sex"].value_counts().to_dict()),
                "timepoint_counts": str(out["timepoint"].value_counts().to_dict()),
                "group_counts_subjects": str(out.drop_duplicates("subject_id")["diagnosis_or_group"].value_counts().to_dict()),
            }
        ]
    )
    summary.to_csv(OUT / "GSE295758_metadata_summary.tsv", sep="\t", index=False)
    processing = []
    for line in lines:
        if line.startswith(("!Sample_data_processing", "!Series_overall_design", "!Series_summary", "!Sample_source_name_ch1")):
            _, vals = split_geo(line)
            text = vals[0] if vals else ""
            if text not in [p.split("\t", 1)[-1] for p in processing]:
                processing.append(line.split("\t", 1)[0] + "\t" + text)
    (OUT / "GSE295758_processing_lines.txt").write_text("\n".join(processing) + "\n")
    suppl = pd.DataFrame(
        [
            {
                "file": "GSE295758_combatVstTranscriptCounts.csv.gz",
                "listing_size": "80M",
                "listing_modified": "2025-04-28 08:05",
                "downloaded": False,
                "source": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295758/suppl/",
            },
            {
                "file": "GSE295758_rawTranscriptCounts.csv.gz",
                "listing_size": "29M",
                "listing_modified": "2025-04-28 08:06",
                "downloaded": False,
                "source": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE295nnn/GSE295758/suppl/",
            },
        ]
    )
    suppl.to_csv(OUT / "GSE295758_supplement_listing.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


def module_clock_search() -> None:
    rows = [
        ("/home/jkang", "home directory, depth 4 name search and depth 5 git search", "no v0.3.0 repository"),
        ("/HDD8T1/WORK/GeroSci/revision", "folder name revision", "not a git repository; contains JAAA-D-26-02318.pdf and major_revision_1.txt only"),
        ("/HDD8T1/WORK/GeroSci/sc-module-aging", "git repository; scripts/00_download.sh mentions Zenodo 22166800 and tAge", "pyproject version 0.2.0; no git tag; no v0.3.0 commit message; download script does not record a checksum"),
        ("/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tAge", "nested git checkout", "tag 1.5.0, not the revision repository v0.3.0"),
        ("/HDD8T1/WORK/GeroSci/clock-specificity", "git repository", "no v0.3.0 tag"),
        ("/SSD1T /HDD6T /HDD8T2 /ssd240", "git repositories to depth 5", "no repository matching v0.3.0"),
    ]
    pd.DataFrame(rows, columns=["location", "what_was_checked", "result"]).to_csv(
        OUT / "module_clock_search.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        [
            {
                "module_clock_dir": "",
                "module_clock_code_dir": "",
                "status": "not_found",
                "part_g": "skipped",
                "reason": "No repository with tag or version v0.3.0 and a checksum-verifying setup script for Zenodo 22166800 was found. No new download.",
            }
        ]
    ).to_csv(OUT / "module_clock_status.tsv", sep="\t", index=False)


def main() -> None:
    append_checksums(
        [
            ("data/raw/GSE113883/GSE113883_series_matrix.txt.gz", RAW / "GSE113883" / "GSE113883_series_matrix.txt.gz"),
            ("data/raw/GSE113883/GSE113883_TPM_processed_data.txt.gz", RAW / "GSE113883" / "GSE113883_TPM_processed_data.txt.gz"),
            ("data/raw/GSE196728/GSE196728_series_matrix.txt.gz", RAW / "GSE196728" / "GSE196728_series_matrix.txt.gz"),
            ("data/raw/GSE196728/GSE196728_Peru_raw.xlsx", RAW / "GSE196728" / "GSE196728_Peru_raw.xlsx"),
            ("data/locked/GSE295758/GSE295758_series_matrix.txt.gz", LOCKED / "GSE295758" / "GSE295758_series_matrix.txt.gz"),
        ]
    )
    circadianphase()
    build_gse113883()
    build_gse196728()
    build_gse295758()
    module_clock_search()


if __name__ == "__main__":
    main()
