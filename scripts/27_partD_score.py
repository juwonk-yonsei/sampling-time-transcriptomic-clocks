#!/usr/bin/env python3
"""Score rodent module clocks on development control samples.

Input follows instruction 3 §7.1. Gene identifiers and coefficients are not written.
"""

from pathlib import Path

import joblib
import pandas as pd
from scmodage.clocks import Clock, predict
from scmodage.preprocess import map_to_mouse_entrez

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage4" / "module_clocks"
CLOCK_DIR = Path("/HDD8T1/WORK/GeroSci/zenodo_22166800_v2")
MODULE_DIR = CLOCK_DIR / "Rodent module clocks"

DATASETS = {
    "GSE39445": {
        "expr": ROOT / "data" / "processed" / "GSE39445_expr_gene_log2.tsv.gz",
        "samples": ROOT / "data" / "processed" / "GSE39445_samples.tsv",
        "value": "log2 microarray",
    },
    "GSE56931": {
        "expr": ROOT / "data" / "processed" / "GSE56931_expr_gene_log2.tsv.gz",
        "samples": ROOT / "data" / "processed" / "GSE56931_samples.tsv",
        "value": "log2 microarray",
    },
    "GSE113883": {
        "expr": ROOT / "results" / "stage3" / "expr" / "GSE113883_log2tpm1.tsv.gz",
        "samples": ROOT / "results" / "stage3" / "GSE113883_samples.tsv",
        "value": "log2(TPM+1)",
    },
}


def clocks():
    rows = [
        ("TA_chrono", "chronological", "composite", CLOCK_DIR / "EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl"),
        ("TA_mort", "mortality", "composite", CLOCK_DIR / "EN_Mortality_Rodents_Multitissue_scaleddiff.pkl"),
    ]
    for path in sorted(MODULE_DIR.glob("EN_Moduleclock_*.pkl")):
        parts = path.name.split("_")
        outcome_file = parts[2]
        module = parts[5]
        if module == "allmodulegenes":
            continue
        outcome = "chronological" if outcome_file == "Chronoage" else "mortality"
        rows.append((f"TM_{module}_{outcome}", outcome, module, path))
    return [Clock(name=name, outcome=outcome, module=module, function="", path=path) for name, outcome, module, path in rows]


def control_matrix(spec):
    expr = pd.read_csv(spec["expr"], sep="\t", index_col=0)
    samples = pd.read_csv(spec["samples"], sep="\t")
    use = samples["condition"].eq("control") & samples["clock_time_h"].notna()
    gsms = samples.loc[use, "gsm"].tolist()
    missing = [gsm for gsm in gsms if gsm not in expr.columns]
    if missing:
        raise SystemExit(f"missing expression columns: {len(missing)}")
    mat = expr.loc[:, gsms]
    values = mat.to_numpy(dtype=float)
    finite = pd.Series(pd.notna(values).all(axis=1), index=mat.index)
    mat = mat.loc[finite]
    mapped = map_to_mouse_entrez(mat, "human")
    center = mapped.mean(axis=1)
    relative = mapped.sub(center, axis=0).T
    relative.index = mapped.columns
    return relative, int(expr.shape[0]), int(finite.sum())


def main():
    catalog = clocks()
    if len([c for c in catalog if c.module != "composite"]) != 46:
        raise SystemExit(f"expected 46 module clocks, got {sum(c.module != 'composite' for c in catalog)}")
    prep = []
    score_rows = []
    cover_rows = []
    for dataset, spec in DATASETS.items():
        print(dataset, flush=True)
        relative, n_symbols, n_finite = control_matrix(spec)
        prep.append({
            "dataset": dataset,
            "value": spec["value"],
            "log_base_conversion": "none",
            "ortholog": "scmodage.preprocess.map_to_mouse_entrez species=human",
            "relative_expression": "mapped value minus mean of control samples in this dataset",
            "n_control_samples": relative.shape[0],
            "n_symbols_input": n_symbols,
            "n_symbols_finite": n_finite,
            "n_mouse_entrez": relative.shape[1],
        })
        for clock in catalog:
            model = joblib.load(clock.path)
            feats = [str(gene) for gene in model.feature_names_in_]
            aligned = relative.reindex(columns=feats)
            n_present = int(aligned.notna().all(axis=0).sum())
            scored = predict(clock, aligned)
            fraction_present = n_present / len(feats)
            for gsm, value in scored.items():
                score_rows.append({
                    "dataset": dataset,
                    "gsm": gsm,
                    "clock": clock.name,
                    "outcome": clock.outcome,
                    "module": clock.module,
                    "value": value,
                })
            cover_rows.append({
                "dataset": dataset,
                "clock": clock.name,
                "outcome": clock.outcome,
                "module": clock.module,
                "n_features": len(feats),
                "n_features_present": n_present,
                "fraction_present": fraction_present,
                "fraction_imputed": 1 - fraction_present,
                "imputer": "pipeline SimpleImputer",
            })
            print(" ", clock.name, n_present, "/", len(feats), flush=True)
    pd.DataFrame(prep).to_csv(OUT / "input_prep.tsv", sep="\t", index=False)
    pd.DataFrame(score_rows).to_csv(OUT / "clock_scores.tsv", sep="\t", index=False)
    pd.DataFrame(cover_rows).to_csv(OUT / "feature_coverage.tsv", sep="\t", index=False)
    print("scores", len(score_rows), flush=True)


if __name__ == "__main__":
    main()
