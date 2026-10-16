#!/usr/bin/env python3
"""Part B structure, coverage, and scoring-code record.

Writes counts only. Module genes and coefficients are not written.
"""

from __future__ import annotations

import gzip
import hashlib
import subprocess
import warnings
import zipfile
from pathlib import Path

import joblib
import pandas as pd
from scmodage.clocks import UNITS
from scmodage.preprocess import map_to_mouse_entrez

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage4" / "module_clocks"
DEST = Path("/HDD8T1/WORK/GeroSci/zenodo_22166800_v2")
ZENODO = pd.read_csv(OUT / "zenodo_v2_files.tsv", sep="\t")
ENV = Path("/HDD8T1/WORK/GeroSci/tod4_clock_env")

NEEDED = [
    "EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
    "EN_Mortality_Rodents_Multitissue_scaleddiff.pkl",
    "EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl",
    "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
    "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl",
    "Rodent module clocks.zip",
    "Multispecies module clocks.zip",
]


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo: str, *args: str) -> str:
    return subprocess.check_output(["git", "-C", repo, *args], text=True).strip()


def tokens(raw: str) -> list[str]:
    blank = {"", "NA", "N/A", "NULL", "NONE", "---"}
    parts = []
    for piece in raw.split("///"):
        piece = piece.strip().strip('"')
        if piece and piece.upper() not in blank:
            parts.append(piece)
    return list(dict.fromkeys(parts))


def symbols_from_table(path: Path, column: str, gzipped: bool) -> list[str]:
    opener = gzip.open if gzipped else open
    found: list[str] = []
    seen: set[str] = set()
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("!platform_table_begin"):
                header = next(handle).rstrip("\n").split("\t")
                break
        if header is None:
            raise SystemExit(f"no platform table in {path}")
        idx = header.index(column)
        for line in handle:
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= idx:
                continue
            for token in tokens(parts[idx]):
                if token not in seen:
                    seen.add(token)
                    found.append(token)
    return found


def gse113883_symbols() -> list[str]:
    path = ROOT / "results" / "stage3" / "expr" / "GSE113883_log2tpm1.tsv.gz"
    symbols = []
    with gzip.open(path, "rt") as handle:
        next(handle)
        for line in handle:
            symbols.append(line.split("\t", 1)[0])
    return symbols


def mouse_ids(symbols: list[str]) -> set[str]:
    dummy = pd.DataFrame({"s": 1.0}, index=pd.Index(symbols, dtype=str))
    mapped = map_to_mouse_entrez(dummy, "human")
    return set(mapped.index.astype(str))


def load_features(path: Path, warning_rows: list[dict]) -> list[str]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = joblib.load(path)
    for item in caught:
        warning_rows.append(
            {
                "file": path.name,
                "category": item.category.__name__,
                "message": str(item.message),
            }
        )
    names = [str(gene) for gene in model.feature_names_in_]
    return names, model


def main() -> None:
    by_name = {row["name"]: row for _, row in ZENODO.iterrows()}
    placed = []
    for name in NEEDED:
        path = DEST / name
        digest = md5(path)
        placed.append(
            {
                "name": name,
                "path": str(path),
                "size": path.stat().st_size,
                "md5": digest,
                "zenodo_md5": by_name[name]["md5"],
                "md5_match": digest == by_name[name]["md5"],
                "source": "downloaded" if name.endswith(".zip") else "local_hardlink",
            }
        )
    placed_df = pd.DataFrame(placed)
    placed_df.to_csv(OUT / "placed_files.tsv", sep="\t", index=False)
    if not placed_df["md5_match"].all():
        raise SystemExit("placed file md5 does not match Zenodo")

    member_rows = []
    for archive_name in ("Rodent module clocks.zip", "Multispecies module clocks.zip"):
        archive_path = DEST / archive_name
        with zipfile.ZipFile(archive_path) as archive:
            archive.extractall(DEST)
        folder = DEST / archive_name.replace(".zip", "")
        for path in sorted(folder.glob("*.pkl")):
            member_rows.append(
                {
                    "archive": archive_name,
                    "name": path.name,
                    "size": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    pd.DataFrame(member_rows).to_csv(OUT / "module_archive_members.tsv", sep="\t", index=False)

    warning_rows: list[dict] = []
    rodent = DEST / "Rodent module clocks"
    module_rows = []
    feature_sets: dict[tuple[str, str], set[str]] = {}
    step_rows = []
    for path in sorted(rodent.glob("EN_Moduleclock_*.pkl")):
        parts = path.name.split("_")
        outcome = parts[2]
        module = parts[5]
        names, model = load_features(path, warning_rows)
        if not step_rows:
            step_rows.append(
                {
                    "example_file": path.name,
                    "pipeline_class": type(model).__name__,
                    "steps": ";".join(f"{label}:{type(step).__name__}" for label, step in model.steps),
                    "feature_id_pattern": "decimal_string" if all(name.isdigit() for name in names) else "other",
                    "feature_id_species": "mouse",
                    "feature_id_kind": "Entrez",
                }
            )
        feature_sets[(outcome, module)] = set(names)
        if module == "allmodulegenes":
            continue
        module_rows.append(
            {
                "outcome": outcome,
                "module": module,
                "n_features": len(names),
                "file": path.name,
            }
        )
    pd.DataFrame(module_rows).to_csv(OUT / "module_sizes.tsv", sep="\t", index=False)
    pd.DataFrame(step_rows).to_csv(OUT / "pipeline_format.tsv", sep="\t", index=False)

    composite_rows = []
    for outcome, filename in (
        ("Chronoage", "EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl"),
        ("Mortality", "EN_Mortality_Rodents_Multitissue_scaleddiff.pkl"),
    ):
        names, model = load_features(DEST / filename, warning_rows)
        feature_sets[(outcome, "composite")] = set(names)
        composite_rows.append(
            {
                "outcome": outcome,
                "n_features": len(names),
                "paper_n_features": 10487,
                "file": filename,
            }
        )
    pd.DataFrame(composite_rows).to_csv(OUT / "composite_sizes.tsv", sep="\t", index=False)

    structure = []
    for outcome in ("Chronoage", "Mortality"):
        sets = [genes for (got, module), genes in feature_sets.items() if got == outcome and module != "composite" and module != "allmodulegenes"]
        counts: dict[str, int] = {}
        for genes in sets:
            for gene in genes:
                counts[gene] = counts.get(gene, 0) + 1
        structure.append(
            {
                "outcome": outcome,
                "n_modules": len(sets),
                "paper_n_modules": 23,
                "n_unique_genes": len(counts),
                "paper_n_unique_genes": 2007,
                "n_genes_in_more_than_one_module": sum(value > 1 for value in counts.values()),
                "paper_overlap": 0,
                "output_unit": UNITS["chronological" if outcome == "Chronoage" else "mortality"],
            }
        )
    pd.DataFrame(structure).to_csv(OUT / "structure_vs_paper.tsv", sep="\t", index=False)

    if warning_rows:
        pd.DataFrame(warning_rows).drop_duplicates().to_csv(OUT / "version_warnings.tsv", sep="\t", index=False)
    else:
        pd.DataFrame([{"file": "", "category": "", "message": "none"}]).to_csv(
            OUT / "version_warnings.tsv", sep="\t", index=False
        )

    lists = {
        "GPL15331": symbols_from_table(ROOT / "data" / "raw" / "gpl" / "GPL15331_platform.txt", "GENE_SYMBOL", False),
        "GPL10379": symbols_from_table(ROOT / "data" / "raw" / "gpl" / "GPL10379.annot.gz", "Gene symbol", True),
        "GPL24324": (ROOT / "results" / "stage2" / "gpl24324_symbols.txt").read_text().splitlines(),
        "GSE113883": gse113883_symbols(),
    }
    mapped = {name: mouse_ids(symbols) for name, symbols in lists.items()}
    coverage = []
    for (outcome, module), genes in sorted(feature_sets.items()):
        if module == "allmodulegenes":
            continue
        for gene_list, present_ids in mapped.items():
            n_present = sum(gene in present_ids for gene in genes)
            coverage.append(
                {
                    "outcome": outcome,
                    "module": module,
                    "gene_list": gene_list,
                    "n_list_symbols": len(lists[gene_list]),
                    "n_list_mouse_entrez": len(present_ids),
                    "n_features": len(genes),
                    "n_features_present": n_present,
                    "fraction_present": n_present / len(genes) if genes else "",
                }
            )
    coverage_df = pd.DataFrame(coverage)
    coverage_df.to_csv(OUT / "coverage.tsv", sep="\t", index=False)
    composite_cov = coverage_df[coverage_df["module"] == "composite"]
    composite_cov.to_csv(OUT / "coverage_composite.tsv", sep="\t", index=False)
    module_cov = coverage_df[coverage_df["module"] != "composite"]
    summary_rows = []
    for (outcome, gene_list), sub in module_cov.groupby(["outcome", "gene_list"]):
        summary_rows.append(
            {
                "outcome": outcome,
                "gene_list": gene_list,
                "n_modules": len(sub),
                "fraction_min": sub["fraction_present"].min(),
                "fraction_median": sub["fraction_present"].median(),
                "fraction_max": sub["fraction_present"].max(),
                "n_features_present_min": int(sub["n_features_present"].min()),
                "n_features_present_max": int(sub["n_features_present"].max()),
            }
        )
    pd.DataFrame(summary_rows).to_csv(OUT / "coverage_module_summary.tsv", sep="\t", index=False)
    chrono = feature_sets[("Chronoage", "composite")]
    mort = feature_sets[("Mortality", "composite")]
    pd.DataFrame(
        [
            {
                "n_chronoage_features": len(chrono),
                "n_mortality_features": len(mort),
                "n_intersection": len(chrono & mort),
                "n_union": len(chrono | mort),
                "paper_n_features": 10487,
            }
        ]
    ).to_csv(OUT / "composite_feature_overlap.tsv", sep="\t", index=False)

    code_rows = []
    for repo, path, function in (
        ("/HDD8T1/WORK/GeroSci/clock-specificity", "scripts/10_score_human.py", "main"),
        ("/HDD8T1/WORK/GeroSci/clock-specificity", "src/clockspec/preprocess.py", "array_to_diff;counts_to_diff"),
        ("/HDD8T1/WORK/GeroSci/clock-specificity", "src/clockspec/clocks.py", "composite;modules"),
        ("/HDD8T1/WORK/GeroSci/sc-module-aging", "scripts/16_score_bulk.py", "module top level"),
        ("/HDD8T1/WORK/GeroSci/sc-module-aging", "src/scmodage/preprocess.py", "map_to_mouse_entrez;preprocess_stratum"),
        ("/HDD8T1/WORK/GeroSci/sc-module-aging", "src/scmodage/clocks.py", "predict;predict_all"),
    ):
        commit, date = git(repo, "log", "-1", "--format=%H %cI", "--", path).split(" ", 1)
        code_rows.append({"repo": repo, "file": path, "function": function, "commit": commit, "committer_date": date})
    pd.DataFrame(code_rows).to_csv(OUT / "scoring_code.tsv", sep="\t", index=False)

    method = pd.DataFrame(
        [
            {
                "chosen_repo": "/HDD8T1/WORK/GeroSci/clock-specificity",
                "reason": "both repositories contain human-to-rodent ortholog mapping and predict; scripts/10_score_human.py committer date 2026-10-14 is later than scripts/16_score_bulk.py 2026-10-12",
                "ortholog_source": "scmodage.preprocess._ortholog_map",
                "ortholog_table": "tAge tag 1.5.0 metadata/Table_of_orthologs.csv",
                "ortholog_rule": "Entrez human to Entrez mouse; drop_duplicates keep first",
                "code_relative_expression": "per-sample z-score then subtract reference-sample median; microarray log base treated as irrelevant because of that scaling; RNA-seq uses log10(RLE+1)",
                "code_filters": "counts: at least 10 counts in at least 20 percent of samples; microarray: drop genes below the 0.25 quantile of mean expression",
                "part_d_input": "instruction 3 section 7.1: log2 microarray or log2(TPM+1); same ortholog method; relative expression = log minus control-sample gene mean; then predict",
            }
        ]
    )
    method.to_csv(OUT / "scoring_method.tsv", sep="\t", index=False)

    repo_rows = []
    for repo in ("/HDD8T1/WORK/GeroSci/clock-specificity", "/HDD8T1/WORK/GeroSci/sc-module-aging"):
        repo_rows.append(
            {
                "repo": repo,
                "commit": git(repo, "rev-parse", "HEAD"),
                "describe": git(repo, "describe", "--tags", "--always"),
            }
        )
    tage = "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tAge"
    repo_rows.append({"repo": tage, "commit": git(tage, "rev-parse", "HEAD"), "describe": git(tage, "describe", "--tags", "--always")})
    pd.DataFrame(repo_rows).to_csv(OUT / "repo_versions.tsv", sep="\t", index=False)

    versions = subprocess.check_output(
        [str(ENV / "bin" / "pip"), "freeze"], text=True
    )
    (OUT / "venv_freeze.txt").write_text(
        "environment: /HDD8T1/WORK/GeroSci/tod4_clock_env\n"
        "base_python: /home/jkang/anaconda3/bin/python3\n"
        "pin_file: /HDD8T1/WORK/GeroSci/sc-module-aging/requirements.txt\n"
        "installed_pins: joblib==1.6.0 numpy==1.26.4 pandas==2.2.3 scikit-learn==1.3.2 scipy==1.13.1\n"
        "not_installed_from_pin_file: anndata h5py matplotlib openpyxl pyarrow pytest PyYAML scanpy seaborn statsmodels tqdm\n"
        "reason: those packages are not imported by the clock loader or the human scoring functions\n"
        "clock-specificity pin file: none\n\n"
        + versions
    )

    unused_composite = ZENODO[(ZENODO["class"] == "composite_clock") & (ZENODO["used_for_scoring"] == False)]
    training = ZENODO[ZENODO["class"] == "other"]
    judgment = pd.DataFrame(
        [
            {
                "item": "clock_file_set",
                "choice": "module archives plus EN multitissue scaleddiff composites that were present locally or are the module-clock zips",
                "n_composite_pkl_on_record": int((ZENODO["class"] == "composite_clock").sum()),
                "n_composite_placed": 5,
                "n_composite_not_downloaded": int(len(unused_composite)),
                "bytes_composite_not_downloaded": int(unused_composite["size"].sum()),
                "n_training_files": int(len(training)),
                "bytes_training_not_downloaded": int(training["size"].sum()),
                "training_files_downloaded": 0,
                "basis": "tAge README names 60 composite clocks and module-clock archives. Training files are Expression_data and Data_annotation. Part D uses the two rodent multitissue scaled-diff chronological and mortality composites and the 23 rodent module clocks. The other composite pkl files and the training files were not downloaded.",
            }
        ]
    )
    judgment.to_csv(OUT / "clock_file_judgment.tsv", sep="\t", index=False)

    (ROOT / "config" / "paths.yaml").write_text(
        "module_clock_dir: /HDD8T1/WORK/GeroSci/zenodo_22166800_v2\n"
        "module_clock_code_dir: /HDD8T1/WORK/GeroSci/clock-specificity\n"
    )
    print(pd.DataFrame(structure).to_string(index=False))
    print(pd.DataFrame(composite_rows).to_string(index=False))
    print("warnings", len(warning_rows), "coverage rows", len(coverage))


if __name__ == "__main__":
    main()
