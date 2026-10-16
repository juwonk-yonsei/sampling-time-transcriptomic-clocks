#!/usr/bin/env python3
"""QC for stage-3 RNA-seq gene matrices. Same rules as instruction 1 §6.4."""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage3"
EXPR = OUT / "expr"


def r2_ols(y, x):
    if len(y) < 3 or np.nanstd(y) == 0:
        return float("nan")
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    pred = x @ beta
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    ss_res = float(np.sum((y - pred) ** 2))
    if ss_tot == 0:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def categorical_r2(y, labels):
    labels = pd.Series(labels).astype(str)
    mask = np.isfinite(y) & labels.notna().to_numpy() & ~labels.isin(["NA", "nan", ""])
    if mask.sum() < 3:
        return float("nan")
    y_use = y[mask]
    labs = labels.to_numpy()[mask]
    if len(set(labs)) < 2:
        return 0.0
    dummies = pd.get_dummies(pd.Series(labs), drop_first=True)
    x = np.column_stack([np.ones(len(y_use)), dummies.to_numpy(dtype=float)])
    return r2_ols(y_use, x)


def cosinor_r2(y, clock):
    clock_num = pd.to_numeric(clock, errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(y) & np.isfinite(clock_num)
    if mask.sum() < 3:
        return float("nan")
    ang = 2.0 * np.pi * clock_num[mask] / 24.0
    x = np.column_stack([np.ones(mask.sum()), np.cos(ang), np.sin(ang)])
    return r2_ols(y[mask], x)


def qc(dataset, matrix_name, sample_name):
    mat = pd.read_csv(EXPR / matrix_name, sep="\t", index_col=0)
    samples = pd.read_csv(OUT / sample_name, sep="\t")
    samples = samples.set_index("gsm").loc[mat.columns]
    values = mat.to_numpy(dtype=float)
    sample_median = np.nanmedian(values, axis=0)
    sample_iqr = np.nanquantile(values, 0.75, axis=0) - np.nanquantile(values, 0.25, axis=0)
    center = float(np.median(sample_median))
    mad = float(np.median(np.abs(sample_median - center)))
    flagged = np.zeros(len(sample_median), dtype=bool) if mad == 0 else np.abs(sample_median - center) > (3.0 * mad)
    pd.DataFrame(
        {
            "gsm": mat.columns,
            "median_log2": sample_median,
            "iqr_log2": sample_iqr,
            "flag_3mad": flagged,
        }
    ).to_csv(OUT / f"{dataset}_sample_qc.tsv", sep="\t", index=False, float_format="%.10g")
    gene_var = np.nanvar(values, axis=1)
    complete = np.isfinite(values).all(axis=1)
    order = np.argsort(-np.where(complete, gene_var, -np.inf), kind="mergesort")
    top = values[order[:1000]]
    top = top - np.nanmean(top, axis=1, keepdims=True)
    u, s, vt = np.linalg.svd(np.nan_to_num(top), full_matrices=False)
    scores = vt[:5].T * s[:5]
    rows = []
    predictors = {
        "subject_id": ("cat", samples["subject_id"]),
        "condition": ("cat", samples["condition"]),
        "period": ("cat", samples["period"]),
        "batch_fields": ("cat", samples["batch_fields"]),
        "cos_sin_clock_time_h": ("clock", samples["clock_time_h"]),
    }
    for pc in range(5):
        for name, (kind, series) in predictors.items():
            y = scores[:, pc]
            stat = categorical_r2(y, series) if kind == "cat" else cosinor_r2(y, series)
            rows.append({"dataset": dataset, "pc": f"PC{pc+1}", "predictor": name, "r2": stat, "n_samples": int(mat.shape[1]), "n_genes_pca": 1000})
    pd.DataFrame(rows).to_csv(OUT / f"{dataset}_pca_r2.tsv", sep="\t", index=False, float_format="%.10g")
    pd.DataFrame(
        [
            {
                "dataset": dataset,
                "n_genes": int(mat.shape[0]),
                "n_samples": int(mat.shape[1]),
                "median_of_sample_medians": center,
                "mad_of_sample_medians": mad,
                "n_flagged_3mad": int(flagged.sum()),
                "flagged_gsm": ",".join(mat.columns[flagged]),
            }
        ]
    ).to_csv(OUT / f"{dataset}_qc_summary.tsv", sep="\t", index=False, float_format="%.10g")
    print(dataset, "flagged", int(flagged.sum()), "genes", mat.shape[0])


if __name__ == "__main__":
    qc("GSE113883", "GSE113883_log2tpm1.tsv.gz", "GSE113883_samples.tsv")
    if (EXPR / "GSE196728_log2fpkm1.tsv.gz").exists():
        qc("GSE196728", "GSE196728_log2fpkm1.tsv.gz", "GSE196728_samples.tsv")
