#!/usr/bin/env python3
"""Count genes that enter each F3 outcome on three L1 matrices."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from clockspec import clocks
from clockspec.preprocess import array_to_diff

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "confirm"))
from score_clocks import module_gene_sets, score_weights, symbols_to_human_entrez  # noqa: E402

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
F3 = [
    "TM_turquoise_mortality",
    "CMP_TAms_chrono_turquoise",
    "CMP_TAms_chrono_yellow",
    "CMP_TAms_chrono_orange",
    "CMP_TAms_chrono_brown",
    "CMP_TAms_chrono_red",
    "CMP_TAms_mort_turquoise",
    "CMP_TAms_mort_darkred",
]


def participant_mean(expr: pd.DataFrame, subjects: pd.Series) -> tuple[pd.DataFrame, dict]:
    values = expr.to_numpy(dtype=float)
    subjects = subjects.astype(str).to_numpy()
    missing = ~np.isfinite(values)
    n_missing = int(missing.sum())
    filled_participant = 0
    out = values.copy()
    original = np.where(np.isfinite(values), values, np.nan)
    gene_mean = np.nanmean(original, axis=1)
    for subject in pd.unique(subjects):
        cols = np.where(subjects == subject)[0]
        block = out[:, cols]
        observed = np.isfinite(block)
        with np.errstate(all="ignore"):
            person = np.divide(np.where(observed, block, 0).sum(axis=1), observed.sum(axis=1))
        need = ~observed
        can = need & np.isfinite(person)[:, None]
        out[:, cols] = np.where(can, np.broadcast_to(person[:, None], block.shape), block)
        filled_participant += int(can.sum())
    need = ~np.isfinite(out)
    can = need & np.isfinite(gene_mean)[:, None]
    out = np.where(can, np.broadcast_to(gene_mean[:, None], out.shape), out)
    audit = {
        "n_missing": n_missing,
        "n_filled_participant": filled_participant,
        "n_filled_gene_mean": int(need.sum()),
        "n_nonfinite_after": int((~np.isfinite(out)).sum()),
    }
    filled = pd.DataFrame(out, index=expr.index, columns=expr.columns)
    return filled, audit


def entered_counts(expr: pd.DataFrame, samples: pd.DataFrame) -> dict[str, int]:
    order = samples["gsm"].astype(str).tolist()
    mapped = symbols_to_human_entrez(expr.loc[:, order])
    reference = pd.Series(True, index=pd.Index(order, name="gsm"))
    diff, _level = array_to_diff(mapped, "human", reference)
    diff.columns = diff.columns.astype(str)
    counts = {}
    module_sets = {
        "chronological": module_gene_sets("chronological"),
        "mortality": module_gene_sets("mortality"),
    }
    composites = {
        "TAms_chrono": ("chronological", clocks.composite("multispecies", "chronological")),
        "TAms_mort": ("mortality", clocks.composite("multispecies", "mortality")),
    }
    for parent, (outcome, clock) in composites.items():
        weight = clocks.weights(clock).astype(float)
        weight.index = weight.index.astype(str)
        _values, genes = score_weights(diff, weight)
        entered = {gene for gene in genes if float(weight.loc[gene]) != 0.0}
        for module, members in module_sets[outcome].items():
            counts[f"CMP_{parent}_{module}"] = len(entered & members)
    for clock in clocks.modules("multispecies", "mortality"):
        weight = clocks.weights(clock).astype(float)
        weight.index = weight.index.astype(str)
        _values, genes = score_weights(diff, weight)
        counts[f"TM_{clock.module}_mortality"] = len([gene for gene in genes if float(weight.loc[gene]) != 0.0])
    return {clock: counts[clock] for clock in F3}


def main() -> None:
    samples = pd.read_csv(ROOT / "results/confirm/L1/samples_scoring.tsv", sep="\t", dtype=str)
    frozen = pd.read_csv(ROOT / "results/confirm/L1/log2_genes.tsv.gz", sep="\t", index_col=0)
    frozen.columns = frozen.columns.astype(str)
    d2 = pd.read_csv(ROOT / "results/confirm/L1/log2_genes_d1_imputed.tsv", sep="\t", index_col=0)
    d2.columns = d2.columns.astype(str)
    d11, audit = participant_mean(frozen, samples.set_index("gsm").loc[frozen.columns, "subject_id"])
    pd.DataFrame([audit]).to_csv(ROOT / "results/stage8/d11_matrix_audit_recount.tsv", sep="\t", index=False)
    rows = []
    for version, expr in [("frozen", frozen), ("D2_knn", d2), ("D11_participant_mean", d11)]:
        counts = entered_counts(expr, samples)
        for clock, n_entered in counts.items():
            rows.append({"version": version, "clock": clock, "n_entered": n_entered, "title": "post hoc, exploratory"})
    pd.DataFrame(rows).to_csv(ROOT / "results/stage8/f3_feature_counts.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
