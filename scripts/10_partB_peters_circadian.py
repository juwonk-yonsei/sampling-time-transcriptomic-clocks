#!/usr/bin/env python3
"""Part B. Peters coefficient subset check and GSE39445 circadianphase summary.

Does not refit models and does not read locked expression.
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage2"
OUT.mkdir(parents=True, exist_ok=True)

raw = pd.read_excel(
    ROOT / "data" / "raw" / "peters2015" / "ncomms9570-s2.xlsx",
    sheet_name="Data.1",
    header=2,
)
repl = raw.loc[raw["sign-REPL"].astype(str).str.lower() == "y", "NEW-Gene-ID"].astype(str)
repl_set = set(repl)
coefs = pd.read_csv(OUT / "rnaage_coef_genes.tsv", sep="\t")
peters_genes = set(coefs.loc[coefs.signature == "Peters", "gene"])
missing = sorted(peters_genes - repl_set)
audit = pd.DataFrame(
    [
        {
            "signature": "Peters",
            "n_coef_genes": len(peters_genes),
            "candidate_source": "Peters 2015 Supplementary Data 1",
            "n_candidate_symbols": len(repl_set),
            "candidate_filter": "sign-REPL == y",
            "candidate_column": "NEW-Gene-ID",
            "n_in_candidate": len(peters_genes & repl_set),
            "is_subset": len(missing) == 0,
            "n_not_in_candidate": len(missing),
            "not_in_candidate": ",".join(missing),
            "pass_to_next_stage": len(missing) == 0,
        }
    ]
)
audit.to_csv(OUT / "peters2015_subset_check.tsv", sep="\t", index=False)

absent = []
for signature in ("GTExAge", "DESeq2", "Peters"):
    absent.append(
        {
            "signature": signature,
            "separate_candidate_object_in_package": False,
            "namespace_objects": "count2FPKM,genelist_all,genelist_cau,lengthDB,makeplot,predict_age,predict_age_fromse",
            "stored_object": "genelist_all$blood[[signature]] numeric coefficients, including (Intercept)",
            "n_coef_genes_excluding_intercept": int((coefs.signature == signature).sum()),
        }
    )
pd.DataFrame(absent).to_csv(OUT / "rnaage_candidate_objects.tsv", sep="\t", index=False)

symbols = set((OUT / "gpl24324_symbols.txt").read_text().splitlines())
symbols.discard("")
cov_rows = []
for signature, clock in (("GTExAge", "RA_GTEx"), ("Peters", "RA_Peters"), ("DESeq2", "RA_default")):
    genes = set(coefs.loc[coefs.signature == signature, "gene"])
    n_on = len(genes & symbols)
    cov_rows.append(
        {
            "reference": clock,
            "platform": "GPL24324",
            "n_reference_genes": len(genes),
            "n_on_platform": n_on,
            "fraction": n_on / len(genes),
            "note": "clariomshumanhttranscriptcluster.db 8.8.0 PROBEID to SYMBOL",
        }
    )
mcp = pd.read_csv(
    ROOT / "data" / "processed" / "ref" / "MCPcounter_genes.txt", sep="\t"
)
symbol_col = [c for c in mcp.columns if "HUGO" in c][0]
pop_col = [c for c in mcp.columns if "population" in c.lower() or "Cell" in c][0]
for pop, block in mcp.groupby(pop_col, sort=False):
    genes = set(block[symbol_col].astype(str))
    n_on = len(genes & symbols)
    cov_rows.append(
        {
            "reference": f"MCPcounter:{pop}",
            "platform": "GPL24324",
            "n_reference_genes": len(genes),
            "n_on_platform": n_on,
            "fraction": n_on / len(genes) if genes else None,
            "note": "clariomshumanhttranscriptcluster.db 8.8.0 PROBEID to SYMBOL",
        }
    )
all_genes = set(mcp[symbol_col].astype(str))
cov_rows.append(
    {
        "reference": "MCPcounter:ALL",
        "platform": "GPL24324",
        "n_reference_genes": len(all_genes),
        "n_on_platform": len(all_genes & symbols),
        "fraction": len(all_genes & symbols) / len(all_genes),
        "note": "clariomshumanhttranscriptcluster.db 8.8.0 PROBEID to SYMBOL",
    }
)
pd.DataFrame(cov_rows).to_csv(OUT / "gpl24324_coverage.tsv", sep="\t", index=False)

file_row = pd.DataFrame(
    [
        {
            "path": "data/raw/peters2015/ncomms9570-s2.xlsx",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/instance/4639797/bin/ncomms9570-s2.xlsx",
            "paper": "Peters et al. 2015 Nature Communications Supplementary Data 1",
            "sheet": "Data.1",
            "header_row_excel": 3,
            "symbol_column": "NEW-Gene-ID",
            "filter_column": "sign-REPL",
            "filter_value": "y",
            "n_rows_in_sheet": int(raw.shape[0]),
            "n_sign_repl_y": int((raw["sign-REPL"].astype(str).str.lower() == "y").sum()),
            "n_unique_symbols": int(repl.nunique()),
            "used_as_coefficients": False,
        }
    ]
)
file_row.to_csv(OUT / "peters2015_file_record.tsv", sep="\t", index=False)

samples = pd.read_csv(ROOT / "data" / "processed" / "GSE39445_samples.tsv", sep="\t")
phase = pd.to_numeric(samples["circadian_ref_h"], errors="coerce")
dist_rows = []
for label, mask in (
    ("all", pd.Series(True, index=samples.index)),
    ("control", samples["condition"] == "control"),
):
    values = phase.loc[mask]
    observed = values.dropna()
    dist_rows.append(
        {
            "subset": label,
            "n": int(mask.sum()),
            "n_missing": int(values.isna().sum()),
            "min": observed.min(),
            "median": observed.median(),
            "max": observed.max(),
            "n_unique": int(observed.nunique()),
        }
    )
pd.DataFrame(dist_rows).to_csv(OUT / "gse39445_circadianphase_distribution.tsv", sep="\t", index=False)

control = samples.loc[samples["condition"] == "control"].copy()
subjects = sorted(control["subject_id"].dropna().unique())[:2]
example = control.loc[control["subject_id"].isin(subjects)].copy()
example = example.sort_values(["subject_id", "hours_awake", "clock_time_h", "gsm"])
example[
    ["gsm", "subject_id", "clock_time_h", "hours_awake", "circadian_ref_h"]
].to_csv(OUT / "gse39445_circadianphase_examples.tsv", sep="\t", index=False)
print(audit.to_string(index=False))
print(pd.DataFrame(dist_rows).to_string(index=False))
print("example subjects", subjects, "rows", len(example))
