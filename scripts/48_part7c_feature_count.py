#!/usr/bin/env python3
"""Count TAms genes with a non-zero weight that enter the score."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from clockspec import clocks

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "confirm"))
from score_clocks import score_weights, symbols_to_human_entrez  # noqa: E402
from clockspec.preprocess import array_to_diff  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expr", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    expr = pd.read_csv(args.expr, sep="\t", index_col=0)
    expr.columns = expr.columns.astype(str)
    order = samples["gsm"].astype(str).tolist()
    expr = expr.loc[:, order]
    mapped = symbols_to_human_entrez(expr)
    reference = pd.Series(True, index=pd.Index(order, name="gsm"))
    diff, _level = array_to_diff(mapped, "human", reference)
    diff.columns = diff.columns.astype(str)
    rows = []
    catalog = [
        ("TAms_chrono", clocks.composite("multispecies", "chronological"), 1839),
        ("TAms_mort", clocks.composite("multispecies", "mortality"), 10487),
    ]
    for clock_id, clock, stated in catalog:
        weight = clocks.weights(clock).astype(float)
        weight.index = weight.index.astype(str)
        _values, genes = score_weights(diff, weight)
        entered = [gene for gene in genes if float(weight.loc[gene]) != 0.0]
        rows.append({
            "version": args.version,
            "clock": clock_id,
            "n_weight": int((weight != 0).sum()),
            "stated_total": stated,
            "n_entered": len(entered),
        })
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
