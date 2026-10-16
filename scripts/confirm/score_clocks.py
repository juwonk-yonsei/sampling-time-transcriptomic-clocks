#!/usr/bin/env python3
"""Score TAms clocks and module outcomes on a log or count matrix.

The caller passes an already extracted matrix. This script does not open locked files.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from clockspec import clocks
from clockspec.preprocess import array_to_diff, counts_to_diff
from scmodage import preprocess as pp
from scmodage.clocks import load_model


def symbols_to_human_entrez(log_expr: pd.DataFrame) -> pd.DataFrame:
    table = pp._gene_table("human")
    symbol = table["Gene.Symbol"].dropna().astype(str)
    lookup = pd.Series(table.loc[symbol.index, "Entrez"].astype(str).to_numpy(), index=symbol.to_numpy())
    lookup = lookup[~lookup.index.duplicated(keep="first")]
    entrez = pd.Index(log_expr.index.astype(str)).map(lookup)
    keep = pd.notna(entrez)
    mapped = log_expr.loc[keep].copy()
    mapped.index = pd.Index(pd.Index(entrez[keep]).astype(str))
    return mapped.loc[~mapped.index.duplicated(keep="first")]


def score_weights(diff: pd.DataFrame, weight: pd.Series) -> tuple[np.ndarray, list[str]]:
    universe = [gene for gene in diff.columns if bool(diff[gene].notna().all())]
    present = set(universe)
    genes = [gene for gene in weight.index.astype(str) if gene in present]
    if not genes:
        raise SystemExit("clock has no detected genes")
    values = diff.loc[:, genes].to_numpy(dtype=float) @ weight.loc[genes].to_numpy(dtype=float)
    return values, genes


def module_gene_sets(outcome: str) -> dict[str, set[str]]:
    sets = {}
    for clock in clocks.modules("rodent", outcome):
        model = load_model(str(clock.path))
        genes = {str(gene) for gene in model.feature_names_in_}
        if clock.module in sets:
            raise SystemExit(f"duplicate module {outcome} {clock.module}")
        sets[clock.module] = genes
    overlap = 0
    seen: set[str] = set()
    for genes in sets.values():
        if seen & genes:
            overlap += 1
        seen |= genes
    if overlap:
        raise SystemExit(f"rodent module genes overlap for {outcome}")
    return sets


def partial_sums(diff: pd.DataFrame, weight: pd.Series, modules: dict[str, set[str]]) -> dict[str, np.ndarray]:
    values, genes = score_weights(diff, weight)
    gene_index = {gene: i for i, gene in enumerate(genes)}
    used = weight.loc[genes].to_numpy(dtype=float)
    block = diff.loc[:, genes].to_numpy(dtype=float)
    assigned = np.zeros(len(genes), dtype=bool)
    out = {}
    for name, members in modules.items():
        idx = [gene_index[gene] for gene in genes if gene in members]
        if idx:
            assigned[idx] = True
            out[name] = block[:, idx] @ used[idx]
        else:
            out[name] = np.zeros(block.shape[0], dtype=float)
    rest = np.where(~assigned)[0]
    out["outside"] = block[:, rest] @ used[rest] if len(rest) else np.zeros(block.shape[0], dtype=float)
    total = out["outside"] + sum(out[name] for name in modules)
    if float(np.max(np.abs(total - values))) > 1e-8:
        raise SystemExit("component partial sums do not rebuild the weight score")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expr", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--mode", choices=["array", "counts"], required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--components", required=True)
    args = parser.parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    expr = pd.read_csv(args.expr, sep="\t", index_col=0)
    expr.columns = expr.columns.astype(str)
    order = samples["gsm"].astype(str).tolist()
    missing = [gsm for gsm in order if gsm not in expr.columns]
    if missing:
        raise SystemExit(f"expression is missing {len(missing)} samples")
    expr = expr.loc[:, order]
    reference = pd.Series(True, index=pd.Index(order, name="gsm"))
    if args.mode == "array":
        mapped = symbols_to_human_entrez(expr)
        diff, _level = array_to_diff(mapped, "human", reference)
    else:
        diff, _level = counts_to_diff(expr, "human", reference)
    diff.columns = diff.columns.astype(str)
    catalog = [
        ("TAms_chrono", clocks.composite("multispecies", "chronological"), "chronological"),
        ("TAms_mort", clocks.composite("multispecies", "mortality"), "mortality"),
    ]
    for clock in clocks.modules("multispecies", "mortality"):
        catalog.append((f"TM_{clock.module}_mortality", clock, "mortality"))
    rows = []
    components = []
    module_sets = {
        "chronological": module_gene_sets("chronological"),
        "mortality": module_gene_sets("mortality"),
    }
    meta = samples.set_index(samples["gsm"].astype(str))
    for clock_id, clock, outcome in catalog:
        weight = clocks.weights(clock).astype(float)
        weight.index = weight.index.astype(str)
        values, _genes = score_weights(diff, weight)
        for gsm, value in zip(diff.index.astype(str), values):
            rows.append({
                "gsm": gsm,
                "subject_id": meta.loc[gsm, "subject_id"],
                "clock": clock_id,
                "value": value,
            })
        if clock_id in {"TAms_chrono", "TAms_mort"}:
            parts = partial_sums(diff, weight, module_sets[outcome])
            for name, component in parts.items():
                comp_id = f"CMP_{clock_id}_{name}"
                for gsm, value in zip(diff.index.astype(str), component):
                    components.append({
                        "gsm": gsm,
                        "subject_id": meta.loc[gsm, "subject_id"],
                        "clock": comp_id,
                        "value": value,
                    })
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    pd.DataFrame(components).to_csv(args.components, sep="\t", index=False)
    print(f"clocks {len(catalog)} component_rows {len(components)}")


if __name__ == "__main__":
    main()
