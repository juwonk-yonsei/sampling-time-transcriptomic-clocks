#!/usr/bin/env python3
"""Part B and Part C. Rescore tAge with the revision preprocessing and decompose the multispecies composites.

The wrapper only maps symbols to human Entrez. array_to_diff is called unchanged.
Gene identifiers and coefficients are not written.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from clockspec import clocks
from clockspec.preprocess import array_to_diff
from scmodage import preprocess as pp
from scmodage.clocks import Clock, load_model

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage5"
OUT.mkdir(parents=True, exist_ok=True)
ZENODO = Path("/HDD8T1/WORK/GeroSci/zenodo_22166800_v2")
PAPER = {"n_features": 10487, "n_in_module": 1292, "n_outside": 9195}

DATASETS = [
    {
        "dataset": "GSE39445",
        "expr": ROOT / "data/processed/GSE39445_expr_gene_log2.tsv.gz",
        "samples": ROOT / "data/processed/GSE39445_samples.tsv",
        "input_scale": "stage1_log2",
        "reference": "all_control_samples",
    },
    {
        "dataset": "GSE56931",
        "expr": ROOT / "data/processed/GSE56931_expr_gene_log2.tsv.gz",
        "samples": ROOT / "data/processed/GSE56931_samples.tsv",
        "input_scale": "stage1_log2",
        "reference": "all_control_samples",
    },
    {
        "dataset": "GSE113883",
        "expr": ROOT / "results/stage3/expr/GSE113883_log2tpm1.tsv.gz",
        "samples": ROOT / "results/stage3/GSE113883_samples.tsv",
        "input_scale": "log2_TPM_plus_1",
        "reference": "all_153_samples_all_control",
    },
]


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def symbols_to_human_entrez(log_expr: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Match the gene-table symbol lookup. One row per human Entrez, first row kept."""
    table = pp._gene_table("human")
    symbol = table["Gene.Symbol"].dropna().astype(str)
    lookup = pd.Series(table.loc[symbol.index, "Entrez"].astype(str).to_numpy(), index=symbol.to_numpy())
    lookup = lookup[~lookup.index.duplicated(keep="first")]
    entrez = pd.Index(log_expr.index.astype(str)).map(lookup)
    keep = pd.notna(entrez)
    mapped = log_expr.loc[keep].copy()
    mapped.index = pd.Index(pd.Index(entrez[keep]).astype(str))
    n_mapped = int(mapped.shape[0])
    mapped = mapped.loc[~mapped.index.duplicated(keep="first")]
    info = {
        "n_input_rows": int(log_expr.shape[0]),
        "n_symbol_mapped": n_mapped,
        "n_unique_human_entrez": int(mapped.shape[0]),
    }
    return mapped, info


def step_counts(mapped: pd.DataFrame) -> dict:
    """Counts of the same filters array_to_diff applies. The scored matrix comes from array_to_diff."""
    mean = mapped.mean(axis=1)
    kept = mean >= mean.quantile(0.25)
    mouse = pd.Series(mapped.index.astype(str), index=mapped.index).map(pp._ortholog_map("human"))
    mouse_kept = mouse.loc[kept]
    mapped_mouse = mouse_kept.dropna()
    unique_mouse = mapped_mouse.loc[~mapped_mouse.index.duplicated(keep="first")]
    return {
        "n_after_quantile": int(kept.sum()),
        "n_with_mouse_ortholog": int(mapped_mouse.shape[0]),
        "n_unique_mouse_entrez": int((~unique_mouse.duplicated(keep="first")).sum()),
        "quantile": 0.25,
    }


def rodent_module_sets() -> dict[str, dict[str, set[str]]]:
    sets: dict[str, dict[str, set[str]]] = {"Chronoage": {}, "Mortality": {}}
    for clock in clocks.modules("rodent", "chronological") + clocks.modules("rodent", "mortality"):
        token = "Chronoage" if clock.outcome == "chronological" else "Mortality"
        model = load_model(str(clock.path))
        genes = {str(gene) for gene in model.feature_names_in_}
        if clock.module in sets[token]:
            raise SystemExit(f"duplicate rodent module {token} {clock.module}")
        sets[token][clock.module] = genes
    return sets


def feature_ids(path: Path) -> set[str]:
    model = load_model(str(path))
    return {str(gene) for gene in model.feature_names_in_}


def paper_check(module_sets: dict[str, dict[str, set[str]]]) -> pd.DataFrame:
    rows = []
    for spec in ("multispecies", "rodent"):
        for outcome, token in (("chronological", "Chronoage"), ("mortality", "Mortality")):
            clock = clocks.composite(spec, outcome)
            feats = feature_ids(clock.path)
            union: set[str] = set()
            overlap_modules = 0
            for genes in module_sets[token].values():
                if union & genes:
                    overlap_modules += 1
                union |= genes
            n_in = len(feats & union)
            n_out = len(feats - union)
            match = (
                len(feats) == PAPER["n_features"]
                and n_in == PAPER["n_in_module"]
                and n_out == PAPER["n_outside"]
            )
            rows.append(
                {
                    "clock_set": spec,
                    "outcome": outcome,
                    "filename": clock.path.name,
                    "code_path_exists": clock.path.exists(),
                    "n_features": len(feats),
                    "n_rodent_modules": len(module_sets[token]),
                    "n_rodent_module_genes": len(union),
                    "n_modules_with_shared_gene": overlap_modules,
                    "n_in_rodent_modules": n_in,
                    "n_outside_rodent_modules": n_out,
                    "paper_n_features": PAPER["n_features"],
                    "paper_n_in_module": PAPER["n_in_module"],
                    "paper_n_outside": PAPER["n_outside"],
                    "all_three_match": match,
                    "passed_to_next_stage": match,
                }
            )
    return pd.DataFrame(rows)


def file_identity() -> pd.DataFrame:
    rows = []
    clocks_to_check = [
        clocks.composite("multispecies", "chronological"),
        clocks.composite("multispecies", "mortality"),
        clocks.composite("rodent", "chronological"),
        clocks.composite("rodent", "mortality"),
        *clocks.modules("multispecies", "mortality"),
        *clocks.modules("rodent", "chronological"),
        *clocks.modules("rodent", "mortality"),
    ]
    for clock in clocks_to_check:
        zenodo = ZENODO / clock.path.name
        if not zenodo.exists():
            parent = ZENODO / clock.path.parent.name
            zenodo = parent / clock.path.name
        code_md5 = md5(clock.path)
        zenodo_md5 = md5(zenodo) if zenodo.exists() else ""
        rows.append(
            {
                "filename": clock.path.name,
                "code_bytes": clock.path.stat().st_size,
                "zenodo_bytes": zenodo.stat().st_size if zenodo.exists() else 0,
                "md5_equal": code_md5 == zenodo_md5 and zenodo_md5 != "",
            }
        )
    return pd.DataFrame(rows)


def catalog() -> list[dict]:
    rows = [
        spec("TAms_chrono", clocks.composite("multispecies", "chronological"), "primary", "human_catalog"),
        spec("TAms_mort", clocks.composite("multispecies", "mortality"), "primary", "human_catalog"),
        spec("TArd_chrono", clocks.composite("rodent", "chronological"), "sensitivity", "rodent_composite"),
        spec("TArd_mort", clocks.composite("rodent", "mortality"), "sensitivity", "rodent_composite"),
    ]
    for clock in clocks.modules("multispecies", "mortality"):
        rows.append(spec(f"TM_{clock.module}_mortality", clock, "primary", "human_catalog"))
    for clock in clocks.modules("rodent", "chronological"):
        rows.append(spec(f"TMrd_{clock.module}_chronological", clock, "sensitivity", "rodent_module_rescore"))
    for clock in clocks.modules("rodent", "mortality"):
        rows.append(spec(f"TMrd_{clock.module}_mortality", clock, "sensitivity", "rodent_module_rescore"))
    return rows


def spec(clock_id: str, clock: Clock, role: str, catalog_name: str) -> dict:
    return {
        "clock": clock_id,
        "role": role,
        "catalog": catalog_name,
        "outcome": clock.outcome,
        "module": clock.module,
        "filename": clock.path.name,
        "object": clock,
    }


def weight_vector(clock: Clock) -> pd.Series:
    """Selected nonzero genes, the vector one_clock dots with the scaled-diff matrix."""
    return clocks.weights(clock).astype(float)


def score_weights(diff: pd.DataFrame, weight: pd.Series) -> tuple[np.ndarray, int, int]:
    universe = [gene for gene in diff.columns if diff[gene].notna().all()]
    universe_set = set(universe)
    genes = [gene for gene in weight.index if gene in universe_set]
    if not genes:
        raise SystemExit("clock has no detected genes")
    values = diff.loc[:, genes].to_numpy(dtype=float) @ weight.loc[genes].to_numpy(dtype=float)
    return values, len(genes), int(len(weight) - len(genes))


def feature_frame(model, diff: pd.DataFrame) -> pd.DataFrame:
    names = model.feature_names_in_
    frame = diff.reindex(columns=[str(gene) for gene in names])
    frame.columns = names
    return frame


def decompose(model, diff: pd.DataFrame, module_of: dict[str, str]) -> tuple[np.ndarray, dict[str, np.ndarray], float]:
    frame = feature_frame(model, diff)
    features = [str(gene) for gene in frame.columns]
    predicted = np.asarray(model.predict(frame), dtype=float).reshape(-1)
    transformed = model[:-1].transform(frame)
    estimator = model.steps[-1][1]
    coefficient = np.asarray(estimator.coef_, dtype=float).reshape(-1)
    intercept = float(np.asarray(estimator.intercept_).reshape(-1)[0])
    selector = next(step for _, step in model.steps if hasattr(step, "get_support"))
    support = np.asarray(selector.get_support(), dtype=bool)
    selected = [gene for gene, keep in zip(features, support) if keep]
    if transformed.shape[1] != coefficient.shape[0] or len(selected) != coefficient.shape[0]:
        raise SystemExit("pipeline width does not match the coefficient vector")
    terms = transformed * coefficient
    reconstructed = intercept + terms.sum(axis=1)
    max_abs = float(np.max(np.abs(reconstructed - predicted)))
    grouped: dict[str, np.ndarray] = {}
    for index, gene in enumerate(selected):
        name = module_of.get(gene, "outside")
        if name not in grouped:
            grouped[name] = np.zeros(frame.shape[0], dtype=float)
        grouped[name] += terms[:, index]
    return predicted, grouped, max_abs


def within_pearson(frame: pd.DataFrame, left: str, right: str) -> float:
    a = frame[left] - frame.groupby("subject_id")[left].transform("mean")
    b = frame[right] - frame.groupby("subject_id")[right].transform("mean")
    if a.std(ddof=1) == 0 or b.std(ddof=1) == 0:
        return float("nan")
    return float(a.corr(b))


def main() -> None:
    module_sets = rodent_module_sets()
    identity = file_identity()
    identity.to_csv(OUT / "code_vs_zenodo_md5.tsv", sep="\t", index=False)
    if not bool(identity["md5_equal"].all()):
        raise SystemExit("a revision-code clock file differs from the md5-checked Zenodo copy")
    check = paper_check(module_sets)
    check.to_csv(OUT / "paper_feature_check.tsv", sep="\t", index=False)
    if not bool(check.loc[check.clock_set == "multispecies", "all_three_match"].all()):
        raise SystemExit("multispecies composite failed the paper check")
    module_of = {}
    for token, outcome in (("Chronoage", "chronological"), ("Mortality", "mortality")):
        mapping: dict[str, str] = {}
        for color, genes in module_sets[token].items():
            for gene in genes:
                if gene in mapping:
                    raise SystemExit(f"rodent module genes overlap for {token}")
                mapping[gene] = color
        module_of[outcome] = mapping

    rows = catalog()
    pd.DataFrame([{k: v for k, v in row.items() if k != "object"} for row in rows]).to_csv(
        OUT / "clock_catalog.tsv", sep="\t", index=False
    )
    spec_rows = [
        {"item": "human_script", "value": "clock-specificity/scripts/10_score_human.py"},
        {"item": "human_composite_files", "value": "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl;EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl"},
        {"item": "human_modules", "value": "multispecies mortality modules returned by clocks.modules"},
        {"item": "human_module_count", "value": str(sum(row["catalog"] == "human_catalog" and row["module"] != "composite" for row in rows))},
        {"item": "scoring_method", "value": "dot product of effective weights of selected nonzero genes with scaled-diff columns present in every sample"},
        {"item": "intercept_in_weight_score", "value": "false"},
        {"item": "module_score_definition", "value": "weights of that module clock, not a split of the composite coefficient vector"},
        {"item": "predict_also_used_by_human_script", "value": "false"},
        {"item": "gene_filter", "value": "drop genes whose mean is below the 0.25 quantile of gene means"},
        {"item": "standardization_gene_set", "value": "genes remaining after the quantile filter and the ortholog map, before alignment to the clock gene list"},
        {"item": "standardization", "value": "z-score each sample across those genes, ddof=1"},
        {"item": "reference", "value": "subtract the per-gene median of the reference samples after alignment"},
        {"item": "log_transform", "value": "none inside array_to_diff; the input is already on a log scale"},
        {"item": "ortholog", "value": "human Entrez to mouse Entrez through the ortholog map; unmapped rows dropped; duplicate mouse ids keep the first row; values are not summed"},
        {"item": "wrapper", "value": "scripts/31_partB_score.py symbols_to_human_entrez"},
        {"item": "gse113883_input", "value": "log2(TPM+1) because counts are not stored; counts_to_diff was not used"},
        {"item": "code_files_md5", "value": "equal to the stage-4 Zenodo copies"},
    ]
    pd.DataFrame(spec_rows).to_csv(OUT / "preprocess_spec.tsv", sep="\t", index=False)

    score_parts = []
    count_rows = []
    predict_rows = []
    component_parts = []
    check_rows = []
    for item in DATASETS:
        print(f"scoring {item['dataset']}", flush=True)
        samples = pd.read_csv(item["samples"], sep="\t", dtype=str)
        expr = pd.read_csv(item["expr"], sep="\t", index_col=0)
        expr.columns = expr.columns.astype(str)
        missing = [gsm for gsm in samples.gsm.astype(str) if gsm not in expr.columns]
        if missing:
            raise SystemExit(f"{item['dataset']} is missing {len(missing)} sample columns")
        expr = expr.loc[:, samples.gsm.astype(str).tolist()]
        reference = pd.Series(samples.condition.eq("control").to_numpy(), index=samples.gsm.astype(str))
        mapped, map_info = symbols_to_human_entrez(expr)
        counts = step_counts(mapped)
        diff, _level = array_to_diff(mapped, "human", reference)
        diff.columns = diff.columns.astype(str)
        finite = int(diff.notna().all(axis=0).sum())
        count_rows.append(
            {
                "dataset": item["dataset"],
                "input_scale": item["input_scale"],
                "reference": item["reference"],
                "n_samples": int(expr.shape[1]),
                "n_reference": int(reference.sum()),
                **map_info,
                **counts,
                "n_clock_gene_list": int(diff.shape[1]),
                "n_clock_genes_finite_in_all_samples": finite,
            }
        )
        order = pd.Index(diff.index.astype(str))
        lookup = samples.loc[:, ["gsm", "subject_id", "condition"]].copy()
        lookup["gsm"] = lookup.gsm.astype(str)
        meta = lookup.set_index("gsm").loc[order, ["subject_id", "condition"]].copy()
        meta.insert(0, "gsm", order.to_numpy())
        for row in rows:
            weight = weight_vector(row["object"])
            values, n_used, n_dropped = score_weights(diff, weight)
            model = load_model(str(row["object"].path))
            predicted = np.asarray(model.predict(feature_frame(model, diff)), dtype=float).reshape(-1)
            block = meta.copy()
            block["dataset"] = item["dataset"]
            block["clock"] = row["clock"]
            block["role"] = row["role"]
            block["catalog"] = row["catalog"]
            block["outcome"] = row["outcome"]
            block["module"] = row["module"]
            block["value"] = values
            block["n_weight_genes_used"] = n_used
            block["n_weight_genes_absent"] = n_dropped
            score_parts.append(block.loc[:, ["dataset", "gsm", "subject_id", "condition", "clock", "role", "catalog", "outcome", "module", "value", "n_weight_genes_used", "n_weight_genes_absent"]])
            paired = pd.DataFrame({"subject_id": meta["subject_id"].to_numpy(), "weight": values, "predict": predicted})
            paired = paired.loc[meta["condition"].eq("control").to_numpy()].copy()
            predict_rows.append(
                {
                    "dataset": item["dataset"],
                    "clock": row["clock"],
                    "n_control": int(len(paired)),
                    "within_pearson_weight_vs_predict": within_pearson(paired, "weight", "predict"),
                }
            )
            if row["clock"] in {"TAms_chrono", "TAms_mort"}:
                predicted_d, grouped, max_abs = decompose(model, diff, module_of[row["outcome"]])
                check_rows.append(
                    {
                        "dataset": item["dataset"],
                        "clock": row["clock"],
                        "n_samples": int(len(predicted_d)),
                        "max_abs_component_sum_plus_intercept_minus_predict": max_abs,
                        "n_components": int(len(grouped)),
                    }
                )
                for name, component in grouped.items():
                    comp = meta.copy()
                    comp["dataset"] = item["dataset"]
                    comp["clock"] = row["clock"]
                    comp["component"] = name
                    comp["value"] = component
                    component_parts.append(comp.loc[:, ["dataset", "gsm", "subject_id", "condition", "clock", "component", "value"]])
        del expr, mapped, diff

    scores = pd.concat(score_parts, ignore_index=True)
    scores.to_csv(OUT / "clock_scores.tsv", sep="\t", index=False)
    pd.DataFrame(count_rows).to_csv(OUT / "preprocess_counts.tsv", sep="\t", index=False)
    pd.DataFrame(predict_rows).to_csv(OUT / "weight_vs_predict.tsv", sep="\t", index=False)
    components = pd.concat(component_parts, ignore_index=True)
    components.to_csv(OUT / "component_scores.tsv", sep="\t", index=False)
    pd.DataFrame(check_rows).to_csv(OUT / "decomposition_check.tsv", sep="\t", index=False)

    stage4 = pd.read_csv(ROOT / "results/stage4/module_clocks/clock_scores.tsv", sep="\t", dtype={"gsm": str})
    pair_map = {"TArd_chrono": "TA_chrono", "TArd_mort": "TA_mort"}
    for row in rows:
        if row["catalog"] != "rodent_module_rescore":
            continue
        pair_map[row["clock"]] = f"TM_{row['module']}_{row['outcome']}"
    corr_rows = []
    control = scores.loc[scores.condition.eq("control")]
    for dataset, left in control.groupby(["dataset", "clock"], sort=True):
        dataset_name, clock_id = dataset
        old_name = pair_map.get(clock_id)
        if old_name is None:
            corr_rows.append(
                {
                    "dataset": dataset_name,
                    "clock": clock_id,
                    "stage4_clock": "",
                    "n_control": int(len(left)),
                    "within_pearson": np.nan,
                    "paired": False,
                }
            )
            continue
        old = stage4.loc[(stage4.dataset == dataset_name) & (stage4.clock == old_name), ["gsm", "value"]]
        merged = left.loc[:, ["gsm", "subject_id", "value"]].merge(old, on="gsm", suffixes=("_new", "_old"))
        corr_rows.append(
            {
                "dataset": dataset_name,
                "clock": clock_id,
                "stage4_clock": old_name,
                "n_control": int(len(merged)),
                "within_pearson": within_pearson(merged.rename(columns={"value_new": "new", "value_old": "old"}), "new", "old"),
                "paired": True,
            }
        )
    correlations = pd.DataFrame(corr_rows)
    correlations.to_csv(OUT / "within_person_vs_stage4.tsv", sep="\t", index=False)
    paired = correlations.loc[correlations.paired]
    print("min within pearson", float(paired.within_pearson.min()))
    print("max decomposition diff", float(pd.DataFrame(check_rows).iloc[:, 3].max()))
    print("scores", scores.shape, "components", components.shape)


if __name__ == "__main__":
    main()
