# Time of blood sampling biases transcriptomic aging clocks

Code, frozen analysis plan, decision and deviation logs, result tables and figure source data for the article "Time of blood sampling biases transcriptomic aging clocks: development, prespecified confirmation, and practical guidance". Within-day variation of two tAge models and three RNAAgeCalc signatures was estimated in three public blood time-series datasets (GSE39445, GSE56931, GSE113883). The resulting predictions were frozen together with the analysis plan and code, and then tested in three datasets (GSE48113, GSE122725, GSE295758) whose expression data stayed locked until the freeze. Post hoc analyses cover a time-axis correction for GSE122725, missing-data handling in GSE48113, a bias simulation for common study designs and mitigation by within-person adjustment for blood-cell composition.

- Article: [Authors]. Time of blood sampling biases transcriptomic aging clocks: development, prespecified confirmation, and practical guidance. GeroScience. [Year]. doi:[DOI]
- This archive: Zenodo doi:[DOI]

The instructions and reports in `audit/` and `reports/` are written in Korean. Everything else is in English.

## Folder structure

| Folder | Content |
|---|---|
| `audit/` | Instructions given to the AI coding agent, the agent's reports, `INDEX.md`, and the lock-order table |
| `config/` | Paths, decision log, checksums of downloaded files, unlock approval (manifest hash), environment record, OSF status |
| `environment/` | Installed R and Python package versions, listed on 2026-10-22 |
| `figures/manuscript/` | Final figures (PDF, EPS, PNG, TIFF) and `source_data/` with the plotted values |
| `figures/submission/` | Journal files: Figs. 1–8 (EPS, TIFF), Figs. S1–S3 (PDF), `esm/ESM_1.pdf`, `esm/ESM_2.xlsx` |
| `figures/` (top level) | Earlier figure drafts from instructions 8–9 |
| `prereg/` | Prespecified plan, frozen parameter tables (`frozen/`), GENCODE 38 reference (`ref/`), `manifest.tsv`, `deviations.tsv`, `implementation_map.tsv`, OSF package (prepared, not uploaded) |
| `reports/` | One report per instruction |
| `results/stage1/`–`results/stage10/` | Outputs of each development, rehearsal and post hoc stage |
| `results/confirm/` | Confirmation run 5, the prespecified confirmatory output |
| `results/confirm_run2_partial/`, `_run3_partial/`, `_run4_partial/` | Outputs of interrupted confirmation runs 2–4 (see "Processed data included") |
| `results/posthoc_stage7c/` | Post hoc GSE122725 time-axis analysis |
| `scripts/` | Numbered analysis scripts, `confirm/` (locked confirmatory pipeline), `00_download_data.sh`, `verify_manifest.py`, `check_lock_order.py` |

## Record of events

Commit hashes, dates and times (Korea Standard Time, UTC+9) and messages from `git log`.

| Event | Commit | Date and time (KST) | Message |
|---|---|---|---|
| Rehearsal of the pipeline on development and synthetic data (R1–R3) | `967723e47b76d4094634bd52d9b056e56d9d210a` | 2026-10-16 16:20:54 | TOD6 Part D |
| Analysis plan freeze (manifest SHA-256 f19f79c1…) | `8a203ffdc5a5f1e933cd4f64154c0d5972480e32` | 2026-10-16 16:21:17 | TOD6 freeze |
| Pre-unlock record | `082ed49ed6f4733c28dcdcce9b502f7fb02c0862` | 2026-10-16 16:46:36 | TOD7 pre-unlock record |
| Unlock | `dc2bed1a64adb374884128eb91720ffbed3ad3e2` | 2026-10-16 16:49:28 | TOD7 unlock |
| Runs 1–3 stopped; bug fixes 7.3.1–7.3.5 recorded (later confirmed as D7) | `a868eee221933a57196936fcacbada4798485c82` | 2026-10-16 17:21:12 | TOD7 confirm stopped before family P |
| Deviations D1–D8 recorded before family-P results | `7ffde8edc6f8ffbb395482fcf430ad91457d16e1` | 2026-10-16 17:29:12 | TOD7b deviation plan (before results) |
| Code for D1–D8 | `8e38cffef0cc475f2267b176e8aaac9705a98445` | 2026-10-16 18:34:56 | TOD7b code fixes |
| Run 4 partial output | `eac683709ee4a63ce302c01b05fd3d629b02801f` | 2026-10-16 18:54:17 | TOD7b report |
| Deviations D9–D12 recorded before H2 and later results | `6904715f53f7d694ab5f710ec53b9e04ce07da34` | 2026-10-17 09:49:01 | TOD7c deviation plan (before H2+ results) |
| Code for D9–D11 | `bd728bacd95ee487679d0c0380cc99f00d780b88` | 2026-10-17 10:52:31 | TOD7c code fixes |
| Confirmation run 5 (prespecified output) | `e56a1c489a74453ca33a832835fdaa312e3e8371` | 2026-10-17 11:12:06 | TOD7c confirm output run5 |
| Post hoc, instruction 7c: GSE122725 time axis | `f6d8879c19bc710a7353dd9835dedee5efce5010` | 2026-10-17 11:20:50 | TOD7c posthoc L2 |
| Instruction 7c report and independent recomputation | `248b37ec68635bc6527a8da092f29daf70f680b2` | 2026-10-17 11:23:42 | TOD7c report |
| Post hoc, instruction 8: Part A | `67d4c2d7b3570f2b72aba06242706057cc4e635d` | 2026-10-19 11:09:06 | TOD8 Part A |
| Post hoc, instruction 8: Part B | `9317da205567e399907d98a2ae1ef21e50f4aab8` | 2026-10-19 12:03:27 | TOD8 Part B |
| Post hoc, instruction 8: report | `b97d159b3f1e86ce55c6020c8b91b29629f6c7da` | 2026-10-19 12:06:57 | TOD8 report |
| Post hoc, instruction 9: simulation and mitigation | `50e5f00580c580d77fa67215b9dd097528ca8268` | 2026-10-19 16:33:32 | TOD9 Part A–C |
| Post hoc, instruction 9: figures | `c31ddf74158a59a282e8f3546780b50081d74d55` | 2026-10-19 16:33:39 | TOD9 figures |
| Post hoc, instruction 9: report | `226949bc5cbe49c74dfd444f8d33bb680bcfd09b` | 2026-10-19 16:35:21 | TOD9 report |
| Instruction 10: figures | `122d9876599b2651709a8ebc718e75240e81b26b` | 2026-10-19 18:33:36 | TOD10 figures |
| Instruction 11: figures | `653eb14b0b751b45868b7132466076995ee08ac5` | 2026-10-21 16:28:52 | TOD11 figures |
| Instruction 12: figures | `64b482b816479fb03633fffdce9575f1527c6404` | 2026-10-21 18:00:58 | TOD12 figures |
| Instruction 13: figures | `694678a0913bcccb66a7979c97c68259ba7b0dbf` | 2026-10-21 18:33:01 | TOD13 figures |
| Instruction 14: figures | `927a5f35f6dd1ce4db9387fc2a93ce451fc338e6` | 2026-10-21 20:02:14 | TOD14 figures |
| Instruction 14b: figures and supplementary files | `d2dbad34eb5c21e03a2b6581da668fe99cfbf99c` | 2026-10-21 21:12:11 | TOD14b figures and ESM |
| Instruction 15: release check | `85486cce0e97cf0eadd7df88254b6e086d11d8bf` | 2026-10-22 08:49:51 | TOD15 release prep |

## Lock order

No expression values from the locked datasets were committed before the unlock commit.

`scripts/check_lock_order.py` reads the sample identifiers of GSE48113 (285), GSE122725 (141) and GSE295758 (116) from the stored confirmation matrices and searches every blob of every commit for them. A path whose content holds at least 10 identifiers of one locked dataset is listed with the first commit that added the path (`git log --all --diff-filter=A`). The unlock commit is `dc2bed1a64adb374884128eb91720ffbed3ad3e2` (2026-10-16 16:49:28 KST). The full table is `audit/locked_expression_first_commits.tsv`.

| First commit | Time (KST) | Message | Expression matrices | Per-sample value files |
|---|---|---|---|---|
| `a868eee` | 2026-10-16 17:21:12 | TOD7 confirm stopped before family P | 8 | 17 |
| `7ffde8e` | 2026-10-16 17:29:12 | TOD7b deviation plan (before results) | 4 | 9 |
| `eac6837` | 2026-10-16 18:54:17 | TOD7b report | 5 | 10 |
| `e56a1c4` | 2026-10-17 11:12:06 | TOD7c confirm output run5 | 2 | 10 |
| `f6d8879` | 2026-10-17 11:20:50 | TOD7c posthoc L2 | 0 | 2 |
| `122d987` | 2026-10-19 18:33:36 | TOD10 figures | 0 | 3 |

All 70 files (19 expression matrices and 51 files of per-sample values such as clock scores, module components, quality flags and figure source data) were first committed after the unlock commit. None of their blobs appears in any of the 42 commits up to and including the unlock commit.

The 42 pre-unlock trees contain 4 blobs with locked-dataset sample identifiers, none with expression values:

- `results/stage2/GSE122725_wake_ref.tsv` and `results/stage2/rnaseq_headers/GSE295758_header.txt` (same blob as `results/stage3/GSE295758_bang_lines.txt`): sample annotations read from the GEO series headers, which the plan allowed before the freeze.
- `results/stage6/rehearsal/R2/gene_log2.tsv.gz` and `results/stage6/rehearsal/R2/samples_input.tsv`: rehearsal R2, a matrix of N(0, 1) values (seed 20261015) labeled with the GSE122725 sample identifiers. Regenerating it with the same code and seed gives identical content (SHA-256 of the uncompressed table `cbcf0c2db76c092224ee0177d729d8e5a58258c87c12d9df1021721101976887`).

`git log --follow` reports 15:37:46 for 18 of the matrices. That time belongs to `prereg/ref/gencode.v38.basic.annotation.gtf.gz`: rename detection paired the newly added matrices with the GTF file as copies (for example `C068` in commit `a868eee`), and `--follow` then continued along the GTF history. None of the matrix paths exists in that commit. Both columns are in the TSV.

<details>
<summary>All 70 files with locked-dataset expression values or per-sample values</summary>

| Path | Content | Dataset (sample IDs found) | First commit | Time (KST) |
|---|---|---|---|---|
| `results/confirm/L1/flags_3mad.tsv` | per-sample values | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L1/log2_genes.tsv.gz` | expression matrix | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L1/ra_scores.tsv` | per-sample values | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L1/tage_components.tsv` | per-sample values | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L1/tage_scores.tsv` | per-sample values | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L2/flags_3mad.tsv` | per-sample values | GSE122725 (141) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L2/log2_genes.tsv` | expression matrix | GSE122725 (141) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L3/flags_3mad.tsv` | per-sample values | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L3/gene_counts.tsv` | expression matrix | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L3/log2_fpkm.tsv` | expression matrix | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L3/ra_scores.tsv` | per-sample values | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L3/tage_components.tsv` | per-sample values | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm/L3/tage_scores.tsv` | per-sample values | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L1/flags_3mad.tsv` | per-sample values | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L1/log2_genes.tsv.gz` | expression matrix | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L1/ra_scores.tsv` | per-sample values | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L1/tage_components.tsv` | per-sample values | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L1/tage_scores.tsv` | per-sample values | GSE48113 (285) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L2/flags_3mad.tsv` | per-sample values | GSE122725 (141) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L2/log2_genes.tsv` | expression matrix | GSE122725 (141) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L3/flags_3mad.tsv` | per-sample values | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L3/gene_counts.tsv` | expression matrix | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L3/log2_fpkm.tsv` | expression matrix | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L3/tage_components.tsv` | per-sample values | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run2_partial/L3/tage_scores.tsv` | per-sample values | GSE295758 (116) | `a868eee` | 2026-10-16 17:21:12 |
| `results/confirm_run3_partial/L1/flags_3mad.tsv` | per-sample values | GSE48113 (285) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L1/log2_genes.tsv.gz` | expression matrix | GSE48113 (285) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L1/ra_scores.tsv` | per-sample values | GSE48113 (285) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L1/tage_components.tsv` | per-sample values | GSE48113 (285) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L1/tage_scores.tsv` | per-sample values | GSE48113 (285) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L2/flags_3mad.tsv` | per-sample values | GSE122725 (141) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L2/log2_genes.tsv` | expression matrix | GSE122725 (141) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L3/flags_3mad.tsv` | per-sample values | GSE295758 (116) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L3/gene_counts.tsv` | expression matrix | GSE295758 (116) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L3/log2_fpkm.tsv` | expression matrix | GSE295758 (116) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L3/ra_scores.tsv` | per-sample values | GSE295758 (116) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L3/tage_components.tsv` | per-sample values | GSE295758 (116) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run3_partial/L3/tage_scores.tsv` | per-sample values | GSE295758 (116) | `7ffde8e` | 2026-10-16 17:29:12 |
| `results/confirm_run4_partial/L1/flags_3mad.tsv` | per-sample values | GSE48113 (285) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L1/log2_genes.tsv.gz` | expression matrix | GSE48113 (285) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L1/log2_genes_d1_imputed.tsv` | expression matrix | GSE48113 (285) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L1/ra_scores.tsv` | per-sample values | GSE48113 (285) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L1/ra_scores_complete_case.tsv` | per-sample values | GSE48113 (285) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L1/tage_components.tsv` | per-sample values | GSE48113 (285) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L1/tage_scores.tsv` | per-sample values | GSE48113 (285) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L2/flags_3mad.tsv` | per-sample values | GSE122725 (141) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L2/log2_genes.tsv` | expression matrix | GSE122725 (141) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L3/flags_3mad.tsv` | per-sample values | GSE295758 (116) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L3/gene_counts.tsv` | expression matrix | GSE295758 (116) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L3/log2_fpkm.tsv` | expression matrix | GSE295758 (116) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L3/ra_scores.tsv` | per-sample values | GSE295758 (116) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L3/tage_components.tsv` | per-sample values | GSE295758 (116) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm_run4_partial/L3/tage_scores.tsv` | per-sample values | GSE295758 (116) | `eac6837` | 2026-10-16 18:54:17 |
| `results/confirm/L1/log2_genes_d1_imputed.tsv` | expression matrix | GSE48113 (285) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L1/ra_scores_complete_case.tsv` | per-sample values | GSE48113 (285) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L1/ra_scores_d11.tsv` | per-sample values | GSE48113 (285) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L1/tage_components_d11.tsv` | per-sample values | GSE48113 (285) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L1/tage_components_d1_imputed.tsv` | per-sample values | GSE48113 (285) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L1/tage_scores_d11.tsv` | per-sample values | GSE48113 (285) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L1/tage_scores_d1_imputed.tsv` | per-sample values | GSE48113 (285) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L3/combat_gene_counts.tsv` | expression matrix | GSE295758 (116) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L3/combat_tage_components.tsv` | per-sample values | GSE295758 (116) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L3/combat_tage_scores.tsv` | per-sample values | GSE295758 (116) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L3/tage_components_comprehensive.tsv` | per-sample values | GSE295758 (116) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/confirm/L3/tage_scores_comprehensive.tsv` | per-sample values | GSE295758 (116) | `e56a1c4` | 2026-10-17 11:12:06 |
| `results/posthoc_stage7c/l2_tage_components.tsv` | per-sample values | GSE122725 (141) | `f6d8879` | 2026-10-17 11:20:50 |
| `results/posthoc_stage7c/l2_tage_scores.tsv` | per-sample values | GSE122725 (141) | `f6d8879` | 2026-10-17 11:20:50 |
| `figures/manuscript/source_data/Fig2_D.tsv` | per-sample values | GSE48113 (147) | `122d987` | 2026-10-19 18:33:36 |
| `figures/manuscript/source_data/Fig2_E.tsv` | per-sample values | GSE48113 (147) | `122d987` | 2026-10-19 18:33:36 |
| `figures/manuscript/source_data/Fig2_F.tsv` | per-sample values | GSE48113 (147) | `122d987` | 2026-10-19 18:33:36 |
</details>

To rerun the check (about 3 min):

```bash
python scripts/check_lock_order.py
```

## Verifying the frozen manifest

`prereg/manifest.tsv` lists the 26 frozen files with their sizes and SHA-256 hashes. Its own SHA-256 (`f19f79c1988d75e7085a3f570b80a96a2e3350b368335ac80b0ae4fb162120ba`) is the value in `config/unlock_approved.txt` that the locked pipeline required before it would open any locked file.

```bash
python scripts/verify_manifest.py
```

Each file is reported as `match`, `mismatch` or `missing (not redistributed)`. In this repository 25 files match and 1 differs: `scripts/confirm/run_pipeline.R` was changed after the unlock by the bug fixes 7.3.1–7.3.5 (commit `a868eee`, confirmed as deviation D7). For a mismatch, the script also hashes the copy stored in the freeze commit `8a203ff`, which matches the manifest, and lists the commits that changed the file.

Files that are not redistributed can be checked with the same script after downloading them from their source. For the tAge model files:

```bash
python scripts/verify_manifest.py --tage-dir /path/to/zenodo_22166800_v2
```

This compares the 7 Zenodo files used for scoring with the md5 values in `results/stage4/module_clocks/zenodo_v2_files.tsv`, and the 78 files unpacked from the two module archives with the SHA-256 values in `results/stage4/module_clocks/module_archive_members.tsv`.

## Reproducing the analyses

### Environment

The analyses ran on Ubuntu 22.04.5 LTS (AMD Ryzen 5 3600, 12 threads, 64 GB RAM) without GPU and without thread limits.

- R 4.5.3 (conda environment `r_env`). Installed packages: `environment/R_packages.tsv`. Main packages: lme4 2.0-6, RNAAgeCalc 1.22.0, MCPcounter 1.2.0, impute 1.84.0, limma 3.66.0, data.table 1.18.6.1, circular 0.5-2, clariomshumanhttranscriptcluster.db 8.8.0, org.Hs.eg.db 3.22.0.
- Python 3.12.7 for analysis, independent checks and figures: `environment/python_analysis_requirements.txt` (numpy 1.26.4, pandas 2.3.3, scipy 1.13.1, statsmodels 0.14.6, matplotlib 3.9.2, openpyxl 3.1.5, pypdf 6.21.1).
- Python 3.12.7 for clock scoring (`scripts/confirm/score_clocks.py`): `environment/python_clock_env_requirements.txt` (numpy 1.26.4, pandas 2.2.3, scikit-learn 1.3.2, scipy 1.13.1, joblib 1.6.0). It also needs two packages installed from local repositories: `clockspec` 0.1.0 from clock-specificity commit `9d109499f0bf47a7869cd7b4bcffe311b8be9b75` ([clock-specificity repository URL]) and `scmodage` 0.2.0 from sc-module-aging commit `a0e9d090cbabcdde441a599bb6443ed3f5ba0350` ([sc-module-aging repository URL]).
- The versions frozen before confirmation are in `prereg/frozen/environment.tsv`; the run record is `config/env_versions.txt`. No `renv.lock` or pinned requirements file was used during the analyses; the lists in `environment/` were generated afterwards from the same environments.

Paths: 30 scripts contain absolute paths of the analysis machine (`/HDD8T1/WORK/GeroSci/TOD`, `/HDD8T1/WORK/GeroSci/tod4_clock_env`, `/HDD8T1/WORK/GeroSci/zenodo_22166800_v2`, `/home/jkang/anaconda3/...`), and `config/paths.yaml` points to the tAge files. Recreate these paths (for example with symbolic links) or edit them. The frozen confirmatory scripts were not edited for this release because their hashes are part of the manifest.

### Data

```bash
bash scripts/00_download_data.sh
```

The script downloads the series matrices and supplementary files of the six analyzed GEO datasets (GSE39445, GSE56931, GSE113883, GSE48113, GSE122725, GSE295758), the excluded dataset GSE196728, and the platform tables GPL10379, GPL15331 and GPL24324 into `data/`. It then compares every file with `config/raw_checksums.tsv`. Supplementary Data 1 of Peters et al. (2015) is listed for manual download because PMC serves a browser check to scripted requests. GEO may regenerate series matrices and platform tables, so a new download can differ in bytes from the copies hashed on 2026-10-16. GSE150337 was held in reserve and never opened.

### tAge model files

Download the tAge model files from Zenodo, doi [10.5281/zenodo.22166800](https://doi.org/10.5281/zenodo.22166800), version 2 (2026-08-29; MGB Open Access License 1.0; source code https://github.com/Gladyshev-Lab/tAge). The files used are `EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl`, `EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl`, `EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl`, `EN_Mortality_Rodents_Multitissue_scaleddiff.pkl`, `EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl`, `Multispecies module clocks.zip` and `Rodent module clocks.zip`. Put them in the folder named in `config/paths.yaml` (`module_clock_dir`), unpack each archive into a folder of the same name, and check them with `scripts/verify_manifest.py --tage-dir`.

### Run order

Scripts run in numeric order within each stage. Elapsed times are those recorded with `/usr/bin/time -v` in the reports; blank cells were not recorded.

| Stage (instruction) | Scripts | Main outputs | Elapsed (recorded) |
|---|---|---|---|
| 1 | `01`–`07` | `results/stage1/`, `data/processed/` | `05_rnaage.R` 0:44.56 |
| 2 | `08`–`14`, `lib_cosinor.R` | `results/stage2/` | Part B 0:13.36; Part D 20:24.99 |
| 3 | `15`–`22` | `results/stage3/` | `20` 7:47.93; `21` 3:54.91; others under 0:25 |
| 4 | `23`–`29` | `results/stage4/` | `26` 11:39.96; `28` 18:10.03 |
| 5 | `30`–`33` | `results/stage5/` | Part B 2:56.93; Part D 26:17.35; Part E 29:00.35 |
| 6 (freeze, rehearsal) | `34`–`36`, `confirm/`, `41`–`45` | `prereg/`, `results/stage6/` | score rehearsal 8:37.07; R1 12:17.70; R2 0:06.96; R3 4:22.69 |
| 7, 7b, 7c (unlock, confirmation) | `confirm/run_confirm.sh`, `47`–`50` | `results/confirm/`, `results/stage7*/`, `results/posthoc_stage7c/` | run 5 19:03.48; D1 check 58:26.91; D11 check 56:38.20; `50` 0:32.01 |
| 8 (post hoc) | `51`–`55` | `results/stage8/` | `51` 46:00.17 |
| 9 (simulation, mitigation) | `56`–`58` | `results/stage9/` | `56` 1:48.72 |
| 10 | `59`–`61` | `results/stage10/`, figures | `59` 12:37.60 |
| 11–14b (figures, supplementary files) | `62`–`64`, `66` | `figures/manuscript/`, `figures/submission/` | |

The confirmatory pipeline opens files in `data/locked/` only through `load_locked_expression()`, and only when `config/unlock_approved.txt` holds the manifest hash.

## Processed data included

The repository contains 26 gene-by-sample tables (796 MB). All are processed from public GEO data or are synthetic; none is a file as downloaded from GEO. Two of them are 79.8 MB (one blob); none exceeds 100 MB.

| Path | Source GEO | Content | Run | Written by | Size (MB) |
|---|---|---|---|---|---|
| `results/confirm/L1/log2_genes.tsv.gz` | GSE48113 | log2 gene matrix, single-symbol Agilent probes | run 5 (prespecified output) | `scripts/confirm/run_locked_analysis.R` → `scripts/confirm/preprocess_array.py` | 34.3 |
| `results/confirm/L1/log2_genes_d1_imputed.tsv` | GSE48113 | log2 gene matrix with kNN imputation (deviation D1) | run 5 (prespecified output) | `scripts/confirm/run_locked_analysis.R` | 79.8 |
| `results/confirm/L2/log2_genes.tsv` | GSE122725 | log2 gene matrix, Clariom S probe sets collapsed to symbols | run 5 (prespecified output) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l2.R` | 25.5 |
| `results/confirm/L3/gene_counts.tsv` | GSE295758 | transcript counts summed to genes (frozen GENCODE 38 basic table) | run 5 (prespecified output) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l3.R` | 18.9 |
| `results/confirm/L3/log2_fpkm.tsv` | GSE295758 | log2(FPKM + 1) from the gene counts | run 5 (prespecified output) | `scripts/confirm/run_locked_analysis.R` (`counts_to_log2_fpkm`) | 43.0 |
| `results/confirm/L3/combat_gene_counts.tsv` | GSE295758 | ComBat-VST transcript counts summed to genes (sensitivity) | run 5 (prespecified output) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l3.R` | 35.4 |
| `results/confirm_run2_partial/L1/log2_genes.tsv.gz` | GSE48113 | log2 gene matrix, single-symbol Agilent probes | run 2 (interrupted) | `scripts/confirm/run_locked_analysis.R` → `scripts/confirm/preprocess_array.py` | 34.3 |
| `results/confirm_run2_partial/L2/log2_genes.tsv` | GSE122725 | log2 gene matrix, Clariom S probe sets collapsed to symbols | run 2 (interrupted) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l2.R` | 25.5 |
| `results/confirm_run2_partial/L3/gene_counts.tsv` | GSE295758 | transcript counts summed to genes (frozen GENCODE 38 basic table) | run 2 (interrupted) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l3.R` | 18.9 |
| `results/confirm_run2_partial/L3/log2_fpkm.tsv` | GSE295758 | log2(FPKM + 1) from the gene counts | run 2 (interrupted) | `scripts/confirm/run_locked_analysis.R` (`counts_to_log2_fpkm`) | 43.0 |
| `results/confirm_run3_partial/L1/log2_genes.tsv.gz` | GSE48113 | log2 gene matrix, single-symbol Agilent probes | run 3 (interrupted) | `scripts/confirm/run_locked_analysis.R` → `scripts/confirm/preprocess_array.py` | 34.3 |
| `results/confirm_run3_partial/L2/log2_genes.tsv` | GSE122725 | log2 gene matrix, Clariom S probe sets collapsed to symbols | run 3 (interrupted) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l2.R` | 25.5 |
| `results/confirm_run3_partial/L3/gene_counts.tsv` | GSE295758 | transcript counts summed to genes (frozen GENCODE 38 basic table) | run 3 (interrupted) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l3.R` | 18.9 |
| `results/confirm_run3_partial/L3/log2_fpkm.tsv` | GSE295758 | log2(FPKM + 1) from the gene counts | run 3 (interrupted) | `scripts/confirm/run_locked_analysis.R` (`counts_to_log2_fpkm`) | 43.0 |
| `results/confirm_run4_partial/L1/log2_genes.tsv.gz` | GSE48113 | log2 gene matrix, single-symbol Agilent probes | run 4 (interrupted) | `scripts/confirm/run_locked_analysis.R` → `scripts/confirm/preprocess_array.py` | 34.3 |
| `results/confirm_run4_partial/L1/log2_genes_d1_imputed.tsv` | GSE48113 | log2 gene matrix with kNN imputation (deviation D1) | run 4 (interrupted) | `scripts/confirm/run_locked_analysis.R` | 79.8 |
| `results/confirm_run4_partial/L2/log2_genes.tsv` | GSE122725 | log2 gene matrix, Clariom S probe sets collapsed to symbols | run 4 (interrupted) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l2.R` | 25.5 |
| `results/confirm_run4_partial/L3/gene_counts.tsv` | GSE295758 | transcript counts summed to genes (frozen GENCODE 38 basic table) | run 4 (interrupted) | `scripts/confirm/run_locked_analysis.R` + `scripts/confirm/preprocess_l3.R` | 18.9 |
| `results/confirm_run4_partial/L3/log2_fpkm.tsv` | GSE295758 | log2(FPKM + 1) from the gene counts | run 4 (interrupted) | `scripts/confirm/run_locked_analysis.R` (`counts_to_log2_fpkm`) | 43.0 |
| `results/stage6/rehearsal/R1/gene_log2.tsv.gz` | GSE39445 (development) | rehearsal R1: GSE39445 processed in the GSE48113 format | rehearsal (before freeze) | `scripts/43_partD_pipeline.R` (`run_r1`) → `scripts/confirm/preprocess_array.py` | 43.9 |
| `results/stage6/rehearsal/R1/scores_input.tsv.gz` | GSE39445 (development) | rehearsal R1: scoring input written from the same matrix | rehearsal (before freeze) | `scripts/43_partD_pipeline.R` (`finish_cosinor`) | 44.3 |
| `results/stage6/rehearsal/R2/gene_log2.tsv.gz` | synthetic | rehearsal R2: N(0, 1) values, seed 20261015, Clariom S probe IDs, GSE122725 sample IDs from the public annotation | rehearsal (before freeze) | `scripts/43_partD_pipeline.R` (`run_r2`) | 23.5 |
| `results/stage6/rehearsal/R3/inputs/counts_gene_column.tsv.gz` | GSE113883 (development), synthetic | rehearsal R3 input: round(20 × TPM) as GSE295758-shaped counts, gene-name column | rehearsal (before freeze) | `scripts/42_partD_make_r3.py` | 0.7 |
| `results/stage6/rehearsal/R3/inputs/counts_transcript_column.tsv.gz` | GSE113883 (development), synthetic | rehearsal R3 input: same counts, transcript-ID column | rehearsal (before freeze) | `scripts/42_partD_make_r3.py` | 0.8 |
| `results/stage6/rehearsal/R3/gene_column/gene_counts.tsv.gz` | GSE113883 (development), synthetic | rehearsal R3: gene counts parsed from the gene-name input | rehearsal (before freeze) | `scripts/43_partD_pipeline.R` (`run_r3_one`) | 0.7 |
| `results/stage6/rehearsal/R3/transcript_column/gene_counts.tsv.gz` | GSE113883 (development), synthetic | rehearsal R3: gene counts parsed from the transcript-ID input | rehearsal (before freeze) | `scripts/43_partD_pipeline.R` (`run_r3_one`) | 0.7 |

`results/confirm_run2_partial/`, `results/confirm_run3_partial/` and `results/confirm_run4_partial/` keep the outputs of confirmation runs that stopped before completion. They are the record behind the deviations. Runs 1 and 2 stopped on code errors (run 2 on renamed GSE295758 columns) that were fixed as 7.3.1–7.3.5 and confirmed as bug fixes by deviation D7; run 1 left no output folder. Run 3 stopped on empty RNAAgeCalc scores, which deviations D1 and D3 addressed. Run 4 stopped at H2 on a missing composition value, which led to deviations D9–D11, and one of its log lines disclosed the GSE48113 tAge-chronological family-P decision (deviation D12). Only run 5 (`results/confirm/`) is the prespecified confirmatory output.

Rehearsal inputs (before the freeze): R1 is GSE39445 processed into the GSE48113 format. R2 is a seeded matrix of normal random numbers (seed 20261015). R3 is synthetic transcript counts made from GSE113883 TPM values, with synthetic metadata (age, sex and opioid group drawn with seed 20261015).

## Not included

- tAge model files: the MGB Open Access License 1.0 does not allow redistribution; download them from Zenodo (above).
- Module gene memberships and coefficients: they are part of the tAge model files and fall under the same license.
- Raw GEO files (series matrices, supplementary files, platform tables): they are public; `scripts/00_download_data.sh` downloads and checks them.
- Run logs, caches and the manuscript draft: they are machine-specific or belong to the journal submission.

## Use of AI tools

From Section 2.11 of the manuscript:

> [Confirm before submission.] An AI coding agent (Grok, xAI, running in the Cursor editor) executed the analyses by following written, versioned instructions prepared by the authors. An AI assistant (Claude, Anthropic) helped with study design, instruction writing, review of analysis reports and drafting of this manuscript. All instructions, reports and decision logs were archived, and the confirmatory statistics were recomputed independently. The authors reviewed all outputs and take full responsibility for the content.

All instructions and reports are in `audit/` (see `audit/INDEX.md`).

## License

Code is released under the MIT License (`LICENSE`). Result tables, figures and documents are released under CC BY 4.0. Third-party data and files (GEO datasets, GENCODE annotation, tAge model files, Peters et al. 2015 supplementary data) remain under their original licenses.
