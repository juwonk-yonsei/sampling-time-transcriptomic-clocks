#!/usr/bin/env python3
"""Microarray preprocessing for the confirmation pipeline.

Locked series matrices are not opened here. The caller passes a path.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

BLANK = {"", "NA", "N/A", "NULL", "NONE", "---"}


def read_series_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", comment="!", quotechar='"', compression="gzip", index_col=0, low_memory=False)
    frame.index = frame.index.astype(str).str.strip().str.strip('"')
    frame.columns = frame.columns.astype(str).str.strip().str.strip('"')
    values = frame.apply(pd.to_numeric, errors="coerce")
    values.index.name = "probe_id"
    return values


def geo2r_scale(values: pd.DataFrame) -> str:
    flat = values.to_numpy(dtype=float).ravel()
    flat = flat[np.isfinite(flat)]
    q0, q25, q99, q100 = np.quantile(flat, [0, 0.25, 0.99, 1])
    if q99 > 100 or ((q100 - q0) > 50 and q25 > 0):
        return "linear"
    return "log"


def to_log2(values: pd.DataFrame, scale: str) -> pd.DataFrame:
    out = values.astype(float).copy()
    if scale == "linear":
        out = out.mask(out <= 0)
        out = np.log2(out)
    return out


def classify_symbol(raw: str) -> tuple[str, str]:
    text = str(raw).strip().strip('"')
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


def load_soft_symbols(path: Path, column: str) -> pd.Series:
    with path.open("r", encoding="utf-8", errors="replace") as src:
        header = None
        for line in src:
            if line.startswith("!platform_table_begin"):
                header = next(src).rstrip("\n").split("\t")
                break
        if header is None:
            raise SystemExit(f"no platform table in {path}")
        idx = header.index(column)
        ids, symbols = [], []
        for line in src:
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            ids.append(parts[0].strip().strip('"'))
            symbols.append(parts[idx].strip().strip('"') if len(parts) > idx else "")
    return pd.Series(symbols, index=pd.Index(ids, name="probe_id"))


def collapse_highest_mean(log2: pd.DataFrame, symbol_raw: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    classified = symbol_raw.reindex(log2.index).fillna("")
    genes, status = [], []
    for raw in classified.tolist():
        gene, state = classify_symbol(str(raw))
        genes.append(gene)
        status.append(state)
    probe_map = pd.DataFrame({
        "probe_id": log2.index.to_numpy(),
        "symbol_raw": classified.to_numpy(),
        "gene_symbol": genes,
        "status": status,
        "mean_log2": np.nanmean(log2.to_numpy(dtype=float), axis=1),
    })
    candidates = probe_map["status"].eq("candidate")
    ranked = probe_map.loc[candidates].sort_values(
        ["gene_symbol", "mean_log2", "probe_id"], ascending=[True, False, True], kind="mergesort"
    )
    selected_idx = ranked.groupby("gene_symbol", sort=False).head(1).index
    probe_map.loc[candidates, "status"] = "not_selected"
    probe_map.loc[selected_idx, "status"] = "selected"
    selected = probe_map.loc[probe_map["status"].eq("selected")]
    matrix = pd.DataFrame(
        log2.to_numpy(dtype=float)[selected.index.to_numpy(), :],
        index=selected["gene_symbol"].to_numpy(),
        columns=log2.columns,
    )
    matrix.index.name = "gene"
    matrix = matrix.sort_index()
    matrix = matrix.loc[~matrix.index.duplicated(keep="first")]
    return matrix, probe_map


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--symbol-column", default="GENE_SYMBOL")
    parser.add_argument("--samples", required=True, help="TSV with gsm column; columns are subset to these GSM ids")
    parser.add_argument("--out-expr", required=True)
    parser.add_argument("--out-audit", required=True)
    args = parser.parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    expr = read_series_matrix(Path(args.matrix))
    keep = [gsm for gsm in samples["gsm"].astype(str) if gsm in expr.columns]
    if len(keep) != samples["gsm"].nunique():
        missing = sorted(set(samples["gsm"].astype(str)) - set(keep))
        raise SystemExit(f"matrix is missing {len(missing)} samples")
    expr = expr.loc[:, keep]
    scale = geo2r_scale(expr)
    log2 = to_log2(expr, scale)
    symbols = load_soft_symbols(Path(args.platform), args.symbol_column)
    matrix, probe_map = collapse_highest_mean(log2, symbols)
    matrix.to_csv(args.out_expr, sep="\t", compression="gzip")
    audit = pd.DataFrame([{
        "scale_call": scale,
        "n_probes": int(expr.shape[0]),
        "n_samples": int(expr.shape[1]),
        "n_genes": int(matrix.shape[0]),
        "n_dropped_multi": int(probe_map["status"].eq("dropped_multi").sum()),
        "n_dropped_no_symbol": int(probe_map["status"].eq("dropped_no_symbol").sum()),
    }])
    audit.to_csv(args.out_audit, sep="\t", index=False)
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
