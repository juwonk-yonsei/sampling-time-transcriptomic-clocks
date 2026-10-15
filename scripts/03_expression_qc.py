#!/usr/bin/env python3
"""Expression processing and QC for the development series (Part C).

Does not open locked series matrices.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RES = ROOT / "results" / "stage1"

DATASETS = {
    "GSE39445": {
        "matrix": RAW / "GSE39445_series_matrix.txt.gz",
        "platform": RAW / "gpl" / "GPL15331_platform.txt",
        "platform_kind": "soft",
        "symbol_column": "GENE_SYMBOL",
    },
    "GSE56931": {
        "matrix": RAW / "GSE56931_series_matrix.txt.gz",
        "platform": RAW / "gpl" / "GPL10379.annot.gz",
        "platform_kind": "annot",
        "symbol_column": "Gene symbol",
    },
}

BLANK = {"", "NA", "N/A", "NULL", "NONE", "---"}


def load_symbol_map(path: Path, kind: str, column: str) -> pd.Series:
    if kind == "annot":
        handle = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        handle = path.open("r", encoding="utf-8", errors="replace")
    with handle as src:
        header = None
        for line in src:
            if line.startswith("!platform_table_begin"):
                header = next(src).rstrip("\n").split("\t")
                break
        if header is None:
            raise SystemExit(f"no platform table in {path}")
        idx = header.index(column)
        ids = []
        symbols = []
        for line in src:
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            ids.append(parts[0].strip().strip('"'))
            symbols.append(parts[idx].strip().strip('"') if len(parts) > idx else "")
    return pd.Series(symbols, index=pd.Index(ids, name="probe_id"), name="symbol_raw")


def classify_symbol(raw: str) -> tuple[str, str]:
    text = raw.strip().strip('"')
    if text.upper() in BLANK:
        return "NA", "dropped_no_symbol"
    parts = []
    for piece in text.split("///"):
        piece = piece.strip()
        if piece and piece.upper() not in BLANK:
            parts.append(piece)
    uniq = list(dict.fromkeys(parts))
    if len(uniq) == 1:
        return uniq[0], "candidate"
    if len(uniq) == 0:
        return "NA", "dropped_no_symbol"
    return "NA", "dropped_multi"


def load_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        comment="!",
        quotechar='"',
        compression="gzip",
        index_col=0,
        low_memory=False,
    )
    frame.index = frame.index.astype(str).str.strip().str.strip('"')
    frame.columns = frame.columns.astype(str).str.strip().str.strip('"')
    values = frame.apply(pd.to_numeric, errors="coerce")
    values.index.name = "probe_id"
    return values


def quantiles(values: np.ndarray) -> dict[str, float]:
    flat = values.ravel()
    flat = flat[np.isfinite(flat)]
    q0, q25, q99, q100 = np.quantile(flat, [0, 0.25, 0.99, 1])
    return {"q0": float(q0), "q25": float(q25), "q99": float(q99), "q100": float(q100), "n_numeric": int(flat.size)}


def call_scale(q: dict[str, float]) -> str:
    if q["q99"] > 100 or ((q["q100"] - q["q0"]) > 50 and q["q25"] > 0):
        return "linear"
    return "log"


def r2_ols(y: np.ndarray, x: np.ndarray) -> float:
    if len(y) < 3 or np.nanstd(y) == 0:
        return float("nan")
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    pred = x @ beta
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    ss_res = float(np.sum((y - pred) ** 2))
    if ss_tot == 0:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def categorical_r2(y: np.ndarray, labels: pd.Series) -> float:
    mask = np.isfinite(y) & labels.notna().to_numpy() & (labels.astype(str) != "NA")
    if mask.sum() < 3:
        return float("nan")
    y_use = y[mask]
    labs = labels.to_numpy()[mask]
    if len(set(labs)) < 2:
        return 0.0
    dummies = pd.get_dummies(pd.Series(labs), drop_first=True)
    x = np.column_stack([np.ones(len(y_use)), dummies.to_numpy(dtype=float)])
    return r2_ols(y_use, x)


def cosinor_r2(y: np.ndarray, clock: pd.Series) -> float:
    clock_num = pd.to_numeric(clock, errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(y) & np.isfinite(clock_num)
    if mask.sum() < 3:
        return float("nan")
    t = clock_num[mask]
    ang = 2.0 * np.pi * t / 24.0
    x = np.column_stack([np.ones(mask.sum()), np.cos(ang), np.sin(ang)])
    return r2_ols(y[mask], x)


def process(gse: str, spec: dict) -> None:
    print(f"loading {gse}", flush=True)
    expr = load_matrix(spec["matrix"])
    symbols = load_symbol_map(spec["platform"], spec["platform_kind"], spec["symbol_column"])
    values = expr.to_numpy(dtype=float)
    n_input_na = int(np.isnan(values).sum())
    q = quantiles(values)
    scale = call_scale(q)
    n_nonpositive = 0
    log2 = values.copy()
    if scale == "linear":
        nonpos = np.isfinite(log2) & (log2 <= 0)
        n_nonpositive = int(nonpos.sum())
        log2[nonpos] = np.nan
        positive = np.isfinite(log2)
        log2[positive] = np.log2(log2[positive])
    n_na_after = int(np.isnan(log2).sum())

    classified = symbols.reindex(expr.index).fillna("")
    gene = []
    status = []
    for raw in classified.tolist():
        symbol, state = classify_symbol(str(raw))
        gene.append(symbol)
        status.append(state)
    probe_map = pd.DataFrame(
        {
            "probe_id": expr.index.to_numpy(),
            "symbol_raw": classified.to_numpy(),
            "gene_symbol": gene,
            "status": status,
        }
    )
    unmapped = ~symbols.index.isin(expr.index)
    n_platform_only = int((~pd.Index(expr.index).isin(symbols.index)).sum())
    mean_log2 = np.nanmean(log2, axis=1)
    probe_map["mean_log2"] = mean_log2
    candidates = probe_map["status"].eq("candidate")
    # Highest mean log2; ties break by probe_id lexicographic order.
    ranked = probe_map.loc[candidates].sort_values(
        ["gene_symbol", "mean_log2", "probe_id"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    selected_idx = ranked.groupby("gene_symbol", sort=False).head(1).index
    probe_map.loc[candidates, "status"] = "not_selected"
    probe_map.loc[selected_idx, "status"] = "selected"

    selected = probe_map.loc[probe_map["status"].eq("selected")]
    gene_matrix = pd.DataFrame(log2[selected.index.to_numpy(), :], index=selected["gene_symbol"].to_numpy(), columns=expr.columns)
    gene_matrix.index.name = "gene"
    gene_matrix = gene_matrix.sort_index()
    out_expr = PROC / f"{gse}_expr_gene_log2.tsv.gz"
    gene_matrix.to_csv(out_expr, sep="\t", compression="gzip", float_format="%.10g")
    probe_map.to_csv(PROC / f"{gse}_probe_map.tsv", sep="\t", index=False, float_format="%.10g")

    # QC on gene-level log2 values.
    gene_values = gene_matrix.to_numpy(dtype=float)
    sample_median = np.nanmedian(gene_values, axis=0)
    q75 = np.nanquantile(gene_values, 0.75, axis=0)
    q25 = np.nanquantile(gene_values, 0.25, axis=0)
    sample_iqr = q75 - q25
    center = float(np.median(sample_median))
    mad = float(np.median(np.abs(sample_median - center)))
    if mad == 0:
        flagged = np.zeros(len(sample_median), dtype=bool)
    else:
        flagged = np.abs(sample_median - center) > (3.0 * mad)
    meta = pd.read_csv(PROC / f"{gse}_samples.tsv", sep="\t", dtype=str).set_index("gsm")
    qc = pd.DataFrame(
        {
            "gsm": gene_matrix.columns,
            "median_log2": sample_median,
            "iqr_log2": sample_iqr,
            "flag_3mad": flagged,
        }
    )
    qc.to_csv(RES / f"{gse}_sample_qc.tsv", sep="\t", index=False, float_format="%.10g")

    gene_var = np.nanvar(gene_values, axis=1)
    complete = np.isfinite(gene_values).all(axis=1)
    order = np.argsort(-np.where(complete, gene_var, -np.inf), kind="mergesort")
    complete_order = [i for i in order if complete[i]]
    n_top = min(1000, len(complete_order))
    top = complete_order[:n_top]
    top_mat = gene_values[top, :].T  # samples x genes
    top_mat = top_mat - top_mat.mean(axis=0, keepdims=True)
    _u, _s, vt = np.linalg.svd(top_mat, full_matrices=False)
    scores = top_mat @ vt[:5].T
    meta = meta.reindex(gene_matrix.columns)
    rows = []
    predictors = {
        "subject_id": meta["subject_id"],
        "condition": meta["condition"],
        "period": meta["period"],
        "batch_fields": meta["batch_fields"],
    }
    for pc in range(5):
        y = scores[:, pc]
        for name, labels in predictors.items():
            rows.append({"dataset": gse, "pc": f"PC{pc+1}", "predictor": name, "r2": categorical_r2(y, labels), "n": int(np.isfinite(y).sum())})
        rows.append({"dataset": gse, "pc": f"PC{pc+1}", "predictor": "cos_sin_clock_time_h", "r2": cosinor_r2(y, meta["clock_time_h"]), "n": int(pd.to_numeric(meta["clock_time_h"], errors="coerce").notna().sum())})
    pd.DataFrame(rows).to_csv(RES / f"{gse}_pca_r2.tsv", sep="\t", index=False, float_format="%.10g")

    summary = {
        "dataset": gse,
        "symbol_column": spec["symbol_column"],
        "n_probes_matrix": int(expr.shape[0]),
        "n_samples": int(expr.shape[1]),
        "n_input_na": n_input_na,
        "n_numeric": q["n_numeric"],
        "q0": q["q0"],
        "q25": q["q25"],
        "q99": q["q99"],
        "q100": q["q100"],
        "scale_call": scale,
        "n_nonpositive_set_na": n_nonpositive,
        "n_na_after_transform": n_na_after,
        "n_matrix_probes_absent_from_platform": n_platform_only,
        "n_dropped_no_symbol": int(probe_map["status"].eq("dropped_no_symbol").sum()),
        "n_dropped_multi": int(probe_map["status"].eq("dropped_multi").sum()),
        "n_not_selected": int(probe_map["status"].eq("not_selected").sum()),
        "n_selected_genes": int(probe_map["status"].eq("selected").sum()),
        "n_platform_rows_not_in_matrix": int(unmapped.sum()) if False else int((~symbols.index.isin(expr.index)).sum()),
        "median_of_sample_medians": center,
        "mad_of_sample_medians": mad,
        "n_flagged_3mad": int(flagged.sum()),
        "flagged_gsm": ",".join(qc.loc[flagged, "gsm"].tolist()),
        "pca_n_genes": n_top,
        "pca_complete_genes_available": len(complete_order),
    }
    print(summary, flush=True)
    return summary


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    rows = [process(gse, spec) for gse, spec in DATASETS.items()]
    pd.DataFrame(rows).to_csv(RES / "expression_qc_summary.tsv", sep="\t", index=False, float_format="%.10g")


if __name__ == "__main__":
    main()
