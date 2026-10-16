#!/usr/bin/env python3
"""Build synthetic GSE295758-shaped transcript counts from GSE113883 TPM."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
OUT = ROOT / "results" / "stage6" / "rehearsal" / "R3" / "inputs"
OUT.mkdir(parents=True, exist_ok=True)

tpm = pd.read_csv(ROOT / "results" / "stage3" / "expr" / "GSE113883_tpm.tsv.gz", sep="\t", index_col=0)
samples = pd.read_csv(ROOT / "results" / "stage3" / "GSE113883_samples.tsv", sep="\t")
samples["clock_time_h"] = pd.to_numeric(samples["clock_time_h"])
picked = []
for hour, phase in ((8, "morning"), (20, "evening")):
    block = samples.loc[samples["clock_time_h"] == hour].copy()
    block = block.drop_duplicates("subject_id", keep="first")
    block["phase"] = phase
    picked.append(block)
meta = pd.concat(picked, ignore_index=True)
subjects = sorted(set(meta.loc[meta["phase"] == "morning", "subject_id"]) & set(meta.loc[meta["phase"] == "evening", "subject_id"]))
meta = meta.loc[meta["subject_id"].isin(subjects)].copy()
meta["subject_id"] = meta["subject_id"].astype(str)
meta["gsm"] = meta["gsm"].astype(str)

gencode = pd.read_csv(
    ROOT / "prereg" / "frozen" / "gencode_v38_basic_transcripts.tsv",
    sep="\t",
    usecols=["transcript_id", "gene_name"],
)
gencode = gencode.drop_duplicates("gene_name", keep="first")
common = [gene for gene in tpm.index.astype(str) if gene in set(gencode["gene_name"])]
gencode = gencode.set_index("gene_name").loc[common]
counts = np.rint(tpm.loc[common, meta["gsm"].tolist()].to_numpy(dtype=float) * 20.0)
counts = np.where(np.isfinite(counts), counts, 0.0)
counts = np.clip(counts, 0, None)

gene_table = pd.DataFrame(counts, index=common, columns=meta["gsm"].tolist())
gene_table.insert(0, "gene_name", gene_table.index)
gene_table.to_csv(OUT / "counts_gene_column.tsv.gz", sep="\t", index=False)

transcript_table = pd.DataFrame(counts, index=gencode["transcript_id"].to_numpy(), columns=meta["gsm"].tolist())
transcript_table.insert(0, "transcript_id", transcript_table.index)
transcript_table.to_csv(OUT / "counts_transcript_column.tsv.gz", sep="\t", index=False)

rng = np.random.default_rng(20261015)
cov = pd.DataFrame({
    "subject_id": subjects,
    "age": rng.integers(20, 76, size=len(subjects)),
    "sex": rng.choice(["F", "M"], size=len(subjects)),
    "opioid_group": rng.integers(0, 3, size=len(subjects)),
})
meta = meta.merge(cov, on="subject_id", how="left")
meta.to_csv(OUT / "samples.tsv", sep="\t", index=False)
cov.to_csv(OUT / "covariates.tsv", sep="\t", index=False)
print(f"subjects {len(subjects)} genes {len(common)} samples {meta.shape[0]}")
