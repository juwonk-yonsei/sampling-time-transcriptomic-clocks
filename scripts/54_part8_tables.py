#!/usr/bin/env python3
"""Manuscript tables for stage 8, as TSV and docx."""

from __future__ import annotations

import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
OUT = ROOT / "results" / "stage8" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def write_tsv(name: str, frame: pd.DataFrame) -> None:
    frame.to_csv(OUT / f"{name}.tsv", sep="\t", index=False)


def docx_table(name: str, title: str, frame: pd.DataFrame) -> None:
    rows = [list(frame.columns)] + frame.fillna("").astype(str).values.tolist()
    body = [
        "<w:p><w:r><w:t>",
        escape(title),
        "</w:t></w:r></w:p><w:tbl>",
        "<w:tblPr><w:tblW w:w=\"0\" w:type=\"auto\"/></w:tblPr>",
    ]
    for row in rows:
        body.append("<w:tr>")
        for cell in row:
            text = escape(cell).replace("\n", " ")
            body.append(f"<w:tc><w:p><w:r><w:t>{text}</w:t></w:r></w:p></w:tc>")
        body.append("</w:tr>")
    body.append("</w:tbl><w:sectPr/>")
    document = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<w:document xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\"><w:body>"
        + "".join(body)
        + "</w:body></w:document>"
    )
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
</Types>"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>"""
    path = OUT / f"{name}.docx"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("word/document.xml", document)


def emit(name: str, title: str, frame: pd.DataFrame) -> None:
    write_tsv(name, frame)
    docx_table(name, title, frame)


def main() -> None:
    table1 = pd.DataFrame([
        ["GSE39445", "development", "Agilent microarray, leukocytes, extended wakefulness", "leukocytes", "24", "221", "repeated measures", "clock_time_h"],
        ["GSE56931", "development", "Affymetrix microarray, whole blood, baseline control", "whole blood", "14", "130", "repeated measures", "clock_time_h"],
        ["GSE113883", "development", "RNA-seq TPM, whole blood, sampling every 2 h", "whole blood", "11", "153", "repeated measures", "clock_time_h"],
        ["GSE48113", "confirmation L1", "Microarray GPL15331, forced desynchrony, Day 1", "leukocytes", "22", "147", "repeated measures", "clock_time_h"],
        ["GSE122725", "confirmation L2, gate fail", "Microarray GPL24324, period 1 day-oriented schedule", "PBMC", "11", "70", "repeated measures", "wake_ref_h %% 24"],
        ["GSE295758", "confirmation L3", "RNA-seq raw counts, morning and evening", "whole blood", "58", "116", "paired", "morning 8.0, evening 20.0"],
        ["GSE196728", "excluded", "RNA-seq. Sea-level control subset. Decision 4.3 excludes it from later analyses.", "whole blood", "8", "48", "repeated measures", "clock_time_h: 2, 6, 10, 14, 18, 22"],
        ["GSE150337", "held", "RNA-seq. Decision 3.2.GSE150337: not opened.", "확인 못 함", "확인 못 함", "확인 못 함", "확인 못 함", "확인 못 함"],
    ], columns=["dataset", "role", "description", "tissue", "n_participants", "n_samples", "design", "time_axis"])
    table1.insert(0, "title", "Table 1. Datasets (post hoc, exploratory)")
    emit("Table1", "Table 1. Datasets (post hoc, exploratory)", table1)

    gates = []
    for dataset, path in [("L1", "L1"), ("L2", "L2"), ("L3", "L3")]:
        g = pd.read_csv(ROOT / f"results/confirm/{path}/gate.tsv", sep="\t", dtype=str)
        g.insert(0, "dataset", dataset)
        g.insert(0, "block", "gate")
        gates.append(g)
    fam_p = pd.read_csv(ROOT / "results/confirm/family_P.tsv", sep="\t", dtype=str)
    fam_s = pd.read_csv(ROOT / "results/confirm/family_S1.tsv", sep="\t", dtype=str)
    p_cols = ["dataset", "clock", "n_subjects", "n_samples", "two_a_over_sigma", "ratio_p025", "ratio_p975", "r2_within", "r2_p025", "r2_p975", "peak_h", "peak_center", "peak_halfwidth_95", "p", "q_holm", "decision", "mean_delta", "ci_low", "ci_high", "d"]
    fam_p = fam_p[p_cols].copy()
    fam_p.insert(0, "block", "family_P")
    s_cols = ["dataset", "clock", "n_subjects", "n_samples", "two_a_over_sigma", "ratio_p025", "ratio_p975", "r2_within", "r2_p025", "r2_p975", "peak_h", "peak_center", "peak_halfwidth_95", "p", "q_holm", "decision", "mean_delta", "ci_low", "ci_high", "d"]
    fam_s = fam_s[s_cols].copy()
    fam_s.insert(0, "block", "family_S1")
    rep = pd.read_csv(ROOT / "results/confirm/replication.tsv", sep="\t", dtype=str)
    rep.insert(0, "block", "replication")
    h1 = pd.read_csv(ROOT / "results/confirm/h1_secondary.tsv", sep="\t", dtype=str)
    h1.insert(0, "block", "h1_secondary")
    table2 = pd.concat(gates + [fam_p, rep, h1, fam_s], ignore_index=True)
    table2.insert(0, "title", "Table 2. Confirmation results (post hoc, exploratory)")
    emit("Table2", "Table 2. Confirmation results (post hoc, exploratory)", table2)

    blocks = []
    for label, path, cols in [
        ("H2", "h2.tsv", None),
        ("H3", "h3.tsv", None),
        ("H4", "h4.tsv", None),
        ("H5", "h5.tsv", None),
    ]:
        tab = pd.read_csv(ROOT / "results/confirm" / path, sep="\t", dtype=str)
        if cols:
            tab = tab[cols]
        tab.insert(0, "block", label)
        blocks.append(tab)
    s4 = pd.read_csv(ROOT / "results/confirm/s4.tsv", sep="\t", dtype=str)
    s4 = s4.loc[s4["in_f3"] == "TRUE", ["dataset", "clock", "n_subjects", "n_samples", "two_a", "peak_h", "p", "confirmed", "circ_dist", "phi_bar", "mean_delta", "observed_sign", "predicted_sign", "d_exp_median", "sign_agrees"]].copy()
    s4.insert(0, "block", "S4")
    blocks.append(s4)
    count = pd.read_csv(ROOT / "results/confirm/s4_l1_confirmed.tsv", sep="\t", dtype=str)
    count.insert(0, "block", "S4_L1_count")
    blocks.append(count)
    table3 = pd.concat(blocks, ignore_index=True)
    table3.insert(0, "title", "Table 3. Secondary hypotheses (post hoc, exploratory)")
    emit("Table3", "Table 3. Secondary hypotheses (post hoc, exploratory)", table3)

    decisions = pd.read_csv(ROOT / "config/decisions.tsv", sep="\t", dtype=str)
    keep = decisions["id"].astype(str).str.contains(r"\.D(?:[1-9]|1[0-2])$")
    table_s1 = decisions.loc[keep, ["date", "id", "content", "reason", "timing"]].copy()
    table_s1.insert(0, "title", "Table S1. Deviations D1-D12 (post hoc, exploratory)")
    emit("TableS1", "Table S1. Deviations D1-D12 (post hoc, exploratory)", table_s1)

    sens = []
    for label, path in [
        ("D1_complete", "results/confirm/L1/sensitivity_d1_complete_case.tsv"),
        ("D2", "results/confirm/L1/sensitivity_d2_imputed.tsv"),
        ("D6", "results/confirm/L3/sensitivity_d6_comprehensive.tsv"),
        ("D9", "results/confirm/L1/sensitivity_d9_imputed_h2.tsv"),
        ("D11", "results/confirm/L1/sensitivity_d11_participant_mean.tsv"),
        ("drop_flags_L1", "results/confirm/L1/sensitivity_drop_flags.tsv"),
        ("random_slopes", "results/confirm/L1/sensitivity_random_slopes.tsv"),
        ("combat", "results/confirm/L3/sensitivity_combat.tsv"),
        ("drop_flags_L3", "results/confirm/L3/sensitivity_drop_flags.tsv"),
        ("batch2", "results/confirm/L3/sensitivity_batch2.tsv"),
        ("intercept", "results/confirm/sensitivity_intercept.tsv"),
    ]:
        tab = pd.read_csv(ROOT / path, sep="\t", dtype=str)
        tab.insert(0, "block", label)
        sens.append(tab)
    table_s2 = pd.concat(sens, ignore_index=True)
    table_s2.insert(0, "title", "Table S2. Sensitivities (post hoc, exploratory)")
    emit("TableS2", "Table S2. Sensitivities (post hoc, exploratory)", table_s2)

    post = []
    for label, path in [
        ("l2_sleep_start", "results/stage8/l2_sleep_start_search.tsv"),
        ("arntl_alias", "results/stage8/gene_peak_circular_summary_arntl_alias.tsv"),
        ("effects", "results/stage8/effect_summary.tsv"),
        ("l3_prediction", "results/stage8/l3_d_vs_prediction.tsv"),
        ("d11_modules", "results/stage8/l1_f3_d11_model_m.tsv"),
        ("f3_features", "results/stage8/f3_feature_counts.tsv"),
        ("l2_h1", "results/posthoc_stage7c/h1_exploratory.tsv"),
        ("l2_h4", "results/posthoc_stage7c/h4_exploratory.tsv"),
    ]:
        tab = pd.read_csv(ROOT / path, sep="\t", dtype=str)
        if "title" in tab.columns:
            tab = tab.rename(columns={"title": "row_title"})
        tab.insert(0, "block", label)
        post.append(tab)
    table_s3 = pd.concat(post, ignore_index=True)
    table_s3.insert(0, "title", "Table S3. Post-hoc analyses (post hoc, exploratory)")
    emit("TableS3", "Table S3. Post-hoc analyses (post hoc, exploratory)", table_s3)


if __name__ == "__main__":
    main()
