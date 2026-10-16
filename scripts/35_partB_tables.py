#!/usr/bin/env python3
"""Part B. Place the preregistration text and fill tables F1-F6.

Existing sentences are not edited. Tables are appended.
"""

from __future__ import annotations

import csv
import hashlib
import math
import shutil
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECEIVED = Path("/HDD8T1/WORK/GeroSci/TOD_preregistration_v1.md")
DEST = ROOT / "prereg" / "TOD_preregistration_v1.md"
FROZEN = ROOT / "prereg" / "frozen"
RES = ROOT / "results" / "stage5"
STAGE4 = ROOT / "results" / "stage4" / "rnaseq"
ZENODO = Path("/HDD8T1/WORK/GeroSci/zenodo_22166800_v2")

PRIMARY = ["TAms_chrono", "TAms_mort", "RA_default"]
SECONDARY = ["RA_Peters", "RA_GTEx"]
CLOCKS = PRIMARY + SECONDARY
DATASETS = ["GSE39445", "GSE56931", "GSE113883"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(row[col] for col in columns) + "\n")


def markdown(rows: list[dict[str, str]], columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row[col] for col in columns) + " |")
    return "\n".join(lines)


def sign_of(text: str) -> str:
    value = Decimal(text)
    if value > 0:
        return "+"
    if value < 0:
        return "-"
    return "0"


def constant_c(beta_cos: str, beta_sin: str) -> str:
    """c = f(8) - f(20) = -beta_cos + beta_sin * sqrt(3)."""
    with localcontext() as ctx:
        ctx.prec = 50
        value = -Decimal(beta_cos) + Decimal(beta_sin) * Decimal(3).sqrt()
    return format(value, "f")


def main() -> None:
    FROZEN.mkdir(parents=True, exist_ok=True)
    received_hash = sha256(RECEIVED)
    DEST.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(RECEIVED, DEST)
    if sha256(DEST) != received_hash:
        raise SystemExit("placed preregistration bytes differ from the received file")

    sources = {
        "clock_summary_long": RES / "clock_summary_long.tsv",
        "clock_summary_phase": RES / "clock_summary_phase.tsv",
        "expected_d": RES / "gse295758_expected_d.tsv",
        "reproducible": RES / "reproducible_modules.tsv",
        "h1_stage5": RES / "h1_point_estimates.tsv",
        "h1_stage4_rnaseq": STAGE4 / "h1_point_estimates.tsv",
        "g2_clocks": STAGE4 / "g2_clocks.tsv",
        "g2_tams": RES / "g2_tams.tsv",
        "power": RES / "locked_design_power.tsv",
    }
    source_hash = {name: sha256(path) for name, path in sources.items()}
    long = {(r["clock"], r["dataset"]): r for r in load(sources["clock_summary_long"])}
    phase = {r["clock"]: r for r in load(sources["clock_summary_phase"])}
    expected = {r["clock"]: r for r in load(sources["expected_d"])}
    h1 = {(r["dataset"], r["clock"]): r for r in load(sources["h1_stage5"])}
    h1_ra = {(r["dataset"], r["clock"], r["method"]): r for r in load(sources["h1_stage4_rnaseq"])}

    f1_cols = ["clock", "role", "phi_bar", "max_circ_dist", "dataset", "two_a", "sigma", "sigma_u", "peak_h", "lrt_p"]
    f1 = []
    for clock in CLOCKS:
        role = "primary" if clock in PRIMARY else "secondary"
        for dataset in DATASETS:
            src = long[(clock, dataset)]
            ph = phase[clock]
            f1.append({
                "clock": clock,
                "role": role,
                "phi_bar": ph["phi_bar"],
                "max_circ_dist": ph["max_circ_dist"],
                "dataset": dataset,
                "two_a": src["two_a"],
                "sigma": src["sigma_reml"],
                "sigma_u": src["sigma_u"],
                "peak_h": src["peak_h"],
                "lrt_p": src["lrt_p"],
            })
    write_tsv(FROZEN / "F1.tsv", f1, f1_cols)

    f2_cols = ["clock", "role", "d_GSE39445", "d_GSE56931", "d_GSE113883", "d_exp_median", "predicted_sign"]
    f2 = []
    for clock in CLOCKS:
        src = expected[clock]
        f2.append({
            "clock": clock,
            "role": "primary" if clock in PRIMARY else "secondary",
            "d_GSE39445": src["d_GSE39445"],
            "d_GSE56931": src["d_GSE56931"],
            "d_GSE113883": src["d_GSE113883"],
            "d_exp_median": src["d_exp_median"],
            "predicted_sign": sign_of(src["d_exp_median"]),
        })
    write_tsv(FROZEN / "F2.tsv", f2, f2_cols)

    f3_cols = [
        "clock", "family", "module", "parent_clock", "phi_bar",
        "q_GSE39445", "q_GSE56931", "q_GSE113883",
        "peak_GSE39445", "peak_GSE56931", "peak_GSE113883",
        "d_exp_median", "predicted_sign",
    ]
    f3 = []
    for row in load(sources["reproducible"]):
        if row["reproducible"] != "TRUE":
            continue
        peaks = []
        for dataset in DATASETS:
            hit = h1.get((dataset, row["clock"]))
            if hit is None:
                raise SystemExit(f"missing peak for {row['clock']} {dataset}")
            peaks.append(hit["peak_h"])
        exp = expected[row["clock"]]
        f3.append({
            "clock": row["clock"],
            "family": row["family"],
            "module": row["module"],
            "parent_clock": row["parent_clock"],
            "phi_bar": row["phi_bar"],
            "q_GSE39445": row["q_GSE39445"],
            "q_GSE56931": row["q_GSE56931"],
            "q_GSE113883": row["q_GSE113883"],
            "peak_GSE39445": peaks[0],
            "peak_GSE56931": peaks[1],
            "peak_GSE113883": peaks[2],
            "d_exp_median": exp["d_exp_median"],
            "predicted_sign": sign_of(exp["d_exp_median"]),
        })
    if len(f3) != 8:
        raise SystemExit(f"expected 8 reproducible outcomes, found {len(f3)}")
    write_tsv(FROZEN / "F3.tsv", f3, f3_cols)

    f4_cols = ["clock", "dataset", "beta_cos", "beta_sin", "formula", "c"]
    formula = "c = f(8) - f(20); f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24); c = -beta_cos + beta_sin*sqrt(3)"
    f4 = []
    for clock in ("TAms_chrono", "TAms_mort"):
        src = h1[("GSE113883", clock)]
        f4.append({
            "clock": clock,
            "dataset": "GSE113883",
            "beta_cos": src["beta_cos"],
            "beta_sin": src["beta_sin"],
            "formula": formula,
            "c": constant_c(src["beta_cos"], src["beta_sin"]),
        })
    ra = h1_ra[("GSE113883", "RA_default", "lin")]
    f4.append({
        "clock": "RA_default",
        "dataset": "GSE113883",
        "beta_cos": ra["beta_cos"],
        "beta_sin": ra["beta_sin"],
        "formula": formula,
        "c": constant_c(ra["beta_cos"], ra["beta_sin"]),
    })
    write_tsv(FROZEN / "F4.tsv", f4, f4_cols)

    clock_files = [
        ("TAms_chrono", ZENODO / "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl", "7dc8a470f5db20a75d7052829efbd37a"),
        ("TAms_mort", ZENODO / "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl", "87307d7f5572d9932382912d3fd9571c"),
    ]
    f5_cols = ["item", "value"]
    f5 = [
        {"item": "R", "value": "4.5.3"},
        {"item": "lme4", "value": "2.0.6"},
        {"item": "RNAAgeCalc", "value": "1.22.0"},
        {"item": "MCPcounter", "value": "1.2.0"},
        {"item": "impute", "value": "1.84.0"},
        {"item": "limma", "value": "3.66.0"},
        {"item": "clariomshumanhttranscriptcluster.db", "value": "8.8.0"},
        {"item": "org.Hs.eg.db", "value": "3.22.0"},
        {"item": "python", "value": "3.12.7"},
        {"item": "numpy", "value": "1.26.4"},
        {"item": "pandas", "value": "2.2.3"},
        {"item": "scikit-learn", "value": "1.3.2"},
        {"item": "scipy", "value": "1.13.1"},
        {"item": "joblib", "value": "1.6.0"},
        {"item": "clock_python", "value": "/HDD8T1/WORK/GeroSci/tod4_clock_env"},
        {"item": "clock_specificity_commit", "value": "9d109499f0bf47a7869cd7b4bcffe311b8be9b75"},
        {"item": "tAge_tag", "value": "1.5.0"},
    ]
    for name, clock_path, expected_md5 in clock_files:
        found = md5(clock_path)
        if found != expected_md5:
            raise SystemExit(f"md5 mismatch for {clock_path.name}: {found}")
        f5.append({"item": f"md5_{name}", "value": found})
    for name, path in sources.items():
        f5.append({"item": f"sha256_{name}", "value": source_hash[name]})
        f5.append({"item": f"path_{name}", "value": str(path)})
    write_tsv(FROZEN / "F5.tsv", f5, f5_cols)

    f6_cols = ["block", "clock", "design", "param_dataset", "d_exp_median", "power_alpha_0.05", "power_alpha_0.01", "frac_p_lt_0.05", "frac_p_lt_0.01", "n_rep", "n_converged"]
    f6 = []
    for clock in CLOCKS:
        src = expected[clock]
        f6.append({
            "block": "L3",
            "clock": clock,
            "design": "GSE295758",
            "param_dataset": "median_of_three",
            "d_exp_median": src["d_exp_median"],
            "power_alpha_0.05": src["power_alpha_0.05"],
            "power_alpha_0.01": src["power_alpha_0.01"],
            "frac_p_lt_0.05": "",
            "frac_p_lt_0.01": "",
            "n_rep": "",
            "n_converged": "",
        })
    for row in load(sources["power"]):
        if row["clock"] not in CLOCKS:
            continue
        f6.append({
            "block": "L1L2",
            "clock": row["clock"],
            "design": row["design"],
            "param_dataset": row["param_dataset"],
            "d_exp_median": "",
            "power_alpha_0.05": "",
            "power_alpha_0.01": "",
            "frac_p_lt_0.05": row["frac_p_lt_0.05"],
            "frac_p_lt_0.01": row["frac_p_lt_0.01"],
            "n_rep": row["n_rep"],
            "n_converged": row["n_converged"],
        })
    write_tsv(FROZEN / "F6.tsv", f6, f6_cols)

    source_rows = []
    table_sources = {
        "F1": ["clock_summary_long", "clock_summary_phase"],
        "F2": ["expected_d"],
        "F3": ["reproducible", "expected_d", "h1_stage5"],
        "F4": ["h1_stage5", "h1_stage4_rnaseq"],
        "F5": list(sources),
        "F6": ["expected_d", "power"],
    }
    for table, names in table_sources.items():
        for name in names:
            source_rows.append({"table": table, "path": str(sources[name]), "sha256": source_hash[name]})
    write_tsv(FROZEN / "source_sha256.tsv", source_rows, ["table", "path", "sha256"])

    received_record = [{"item": "received_preregistration_sha256", "value": received_hash}]
    write_tsv(FROZEN / "received_preregistration_sha256.tsv", received_record, ["item", "value"])

    blocks = ["", "", "## Appendix. Tables F1-F6", ""]
    labels = {
        "F1": (f1, f1_cols, ["clock_summary_long", "clock_summary_phase"]),
        "F2": (f2, f2_cols, ["expected_d"]),
        "F3": (f3, f3_cols, ["reproducible", "expected_d", "h1_stage5"]),
        "F4": (f4, f4_cols, ["h1_stage5", "h1_stage4_rnaseq"]),
        "F5": (f5, f5_cols, list(sources)),
        "F6": (f6, f6_cols, ["expected_d", "power"]),
    }
    for table, (rows, cols, names) in labels.items():
        blocks.append(f"### Table {table}")
        blocks.append("")
        blocks.append(markdown(rows, cols))
        blocks.append("")
        blocks.append("Source files:")
        blocks.append("")
        for name in names:
            blocks.append(f"- `{sources[name]}` SHA-256 `{source_hash[name]}`")
        blocks.append("")
    with DEST.open("a") as handle:
        if not DEST.read_text().endswith("\n"):
            handle.write("\n")
        handle.write("\n".join(blocks))
        if not blocks[-1].endswith("\n"):
            handle.write("\n")
    print("received", received_hash)
    print("filled", sha256(DEST))
    print("F3", len(f3), "F4 c", [row["c"] for row in f4])


if __name__ == "__main__":
    main()
