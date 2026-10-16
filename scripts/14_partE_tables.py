#!/usr/bin/env python3
"""Build Part E candidate and inspected tables from header files already saved."""

import csv
import json
from collections import Counter
from pathlib import Path

import importlib.util

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
OUT = ROOT / "results" / "stage2"
spec = importlib.util.spec_from_file_location("meta", ROOT / "scripts" / "13_partE_metadata.py")
meta = importlib.util.module_from_spec(spec)
spec.loader.exec_module(meta)


def processing_format(lines):
    text = " ".join(lines).lower()
    hits = []
    for label, token in (
        ("TPM", "tpm"),
        ("FPKM", "fpkm"),
        ("raw counts", "raw count"),
        ("counts", "count"),
    ):
        if token in text and label not in hits:
            hits.append(label)
    return "|".join(hits)


def genome_build(lines):
    for line in lines:
        if line.lower().startswith("genome_build") or line.lower().startswith("assembly"):
            return line
    return ""


def subject_key(chars):
    for name in ("subjectid", "subject", "subject id", "donor", "study subject id #", "individual"):
        if name in chars:
            return name
    return ""


def summarize(acc):
    buckets = meta.series_lines(acc)
    chars = meta.char_map(buckets)
    geos = (buckets.get("!Sample_geo_accession") or [[]])[0]
    proc = meta.unique_processing(buckets)
    titles = (buckets.get("!Sample_title") or [[]])[0]
    row = {
        "accession": acc,
        "title": meta.first(buckets, "!Series_title"),
        "organism": "|".join(sorted({v for row in buckets.get("!Sample_organism_ch1", []) for v in row if v})),
        "library_strategy": "|".join(sorted({v for row in buckets.get("!Sample_library_strategy", []) for v in row if v})),
        "source_name_unique": "|".join(sorted({v for row in buckets.get("!Sample_source_name_ch1", []) for v in row if v})),
        "platform_id": meta.first(buckets, "!Series_platform_id"),
        "pubmed_id": meta.first(buckets, "!Series_pubmed_id"),
        "n_samples": len([g for g in geos if g]),
        "characteristic_names": "|".join(chars),
        "processed_format_tokens": processing_format(proc),
        "genome_or_assembly_line": genome_build(proc),
        "overall_design": meta.first(buckets, "!Series_overall_design"),
        "title_examples": " | ".join(titles[:3]),
    }
    key = subject_key(chars)
    row["subject_field"] = key
    if key:
        counts = Counter(chars[key])
        per = Counter(counts.values())
        row["n_subjects"] = len(counts)
        row["samples_per_subject"] = json.dumps(dict(sorted(per.items())), ensure_ascii=False)
    else:
        row["n_subjects"] = ""
        row["samples_per_subject"] = ""
    # keep compact value examples for every characteristic
    examples = {}
    for name, vals in chars.items():
        examples[name] = {
            "n": len(vals),
            "n_unique": len(set(vals)),
            "top": Counter(vals).most_common(6),
        }
    row["characteristic_examples"] = json.dumps(examples, ensure_ascii=False)
    if "age" in chars:
        ages = []
        for value in chars["age"]:
            try:
                ages.append(float(value))
            except ValueError:
                pass
        if ages:
            row["age_min"] = min(ages)
            row["age_max"] = max(ages)
            row["age_n"] = len(ages)
        else:
            row["age_min"] = row["age_max"] = row["age_n"] = ""
    else:
        row["age_min"] = row["age_max"] = row["age_n"] = ""
    titles = (buckets.get("!Sample_title") or [[]])[0]
    src, subj_vals = title_subjects(acc, titles, chars)
    row["subject_id_source"] = src
    if subj_vals:
        counts = Counter(subj_vals)
        row["n_subjects_resolved"] = len(counts)
        row["samples_per_subject_resolved"] = json.dumps(dict(sorted(Counter(counts.values()).items())))
    else:
        row["n_subjects_resolved"] = ""
        row["samples_per_subject_resolved"] = ""
    row["candidate"], row["candidate_rule"] = classify(row, chars)
    return row


def title_subjects(acc, titles, chars):
    if acc == "GSE295758":
        subjects = [t.split("_")[0] for t in titles if t]
        return "sample_title", subjects
    if acc == "GSE150337":
        subjects = [t.split("_")[0] for t in titles if "_" in t]
        return "sample_title", subjects
    key = subject_key(chars)
    if key:
        return "characteristic:" + key, chars[key]
    return "", []


def classify(row, chars):
    acc = row["accession"]
    lib = row["library_strategy"]
    src = row["source_name_unique"].lower()
    names = set(row["characteristic_names"].split("|")) if row["characteristic_names"] else set()
    if lib in ("miRNA-Seq", "ncRNA-Seq"):
        return "no", "library_strategy_not_mrna_rnaseq"
    if "ctc" in src or "u2os" in src or "hns" in src or "cell line" in src:
        return "no", "source_not_blood_cell_transcriptome"
    bloodish = any(token in src for token in ("blood", "pbmc", "leukocyte"))
    if not bloodish:
        return "no", "source_not_blood"
    if acc == "GSE127530":
        return "no", "time_field_is_fasting_or_postprandial"
    if acc == "GSE346256":
        return "no", "time_field_is_hours_after_serum_shock"
    clock_fields = {"hour", "local_time", "sampling time", "time", "timepoint"}
    if names & clock_fields and acc in {"GSE113883", "GSE150337", "GSE196728", "GSE295758"}:
        return "yes", "repeated_samples_with_collection_time"
    return "no", "no_sample_level_clock_time"


def main():
    rows = [summarize(path.name.split("_")[0]) for path in sorted((OUT / "rnaseq_headers").glob("*_header.txt"))]
    fields = [
        "accession", "title", "organism", "library_strategy", "source_name_unique",
        "platform_id", "pubmed_id", "n_samples", "n_subjects", "subject_field",
        "samples_per_subject", "age_min", "age_max", "age_n", "characteristic_names",
        "processed_format_tokens", "genome_or_assembly_line", "title_examples",
        "overall_design", "characteristic_examples", "subject_id_source",
        "n_subjects_resolved", "samples_per_subject_resolved", "candidate", "candidate_rule",
    ]
    meta.write_tsv(OUT / "partE_header_table.tsv", rows, fields)
    print("rows", len(rows))


if __name__ == "__main__":
    main()
