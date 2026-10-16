#!/usr/bin/env python3
"""Part A. Record instruction-3 decisions. Do not refit or open locked expression."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = ROOT / "config" / "decisions.tsv"
OUT = ROOT / "results" / "stage3"
OUT.mkdir(parents=True, exist_ok=True)

DATE = "2026-10-16"

ROWS = [
    (
        "3.1.1",
        "S1 comparison method is approved. LRT and R2_within drop only cos and sin. The null keeps hours_awake.",
        "Instruction 3 §1.1 item 1. Report 2 §9.",
    ),
    (
        "3.1.2",
        "S2 null and ICC are approved. The null keeps the random slopes. ICC stays blank.",
        "Instruction 3 §1.1 item 2. Report 2 §9.",
    ),
    (
        "3.1.3",
        "S3 definition is approved. S3 is descriptive. No mixed-model LRT. The time range is max(clock_time_h) minus min(clock_time_h). The R2 column is the median of per-subject OLS R2.",
        "Instruction 3 §1.1 item 3. Report 2 §9.",
    ),
    (
        "3.1.4",
        "The seed reset method is approved. This stage uses the same method. Bootstrap seeds are reset to 20261015 per dataset. Random-signature seeds are reset to 20261015 per dataset by clock.",
        "Instruction 3 §1.1 item 4. Report 2 §9.",
    ),
    (
        "3.1.5",
        "The composition-analysis sample is approved: control samples with non-missing clock_time_h.",
        "Instruction 3 §1.1 item 5. Report 2 §9.",
    ),
    (
        "3.1.6",
        "The reduction-ratio interval rule is approved. The metric changes to the rhythm-vector projection ratio in instruction 3 §4.",
        "Instruction 3 §1.1 item 6. Report 2 §9.",
    ),
    (
        "3.1.7",
        "GSE196728 is designated a development dataset. Gene IDs in the processed file are checked in instruction 3 §2.3. It is used only if that check passes.",
        "Instruction 3 §1.1 item 7 and §1.2.",
    ),
    (
        "3.1.8",
        "GSE127530 and GSE346256 are approved as not candidates.",
        "Instruction 3 §1.1 item 8. Report 2 §9.",
    ),
    (
        "3.1.9",
        "miRNA-Seq and ncRNA-Seq are approved as not candidates.",
        "Instruction 3 §1.1 item 9. Report 2 §9.",
    ),
    (
        "3.1.10",
        "The unit of GSE39445 circadianphase is checked from the paper and the data in instruction 3 §2.1. It is not decided in this part.",
        "Instruction 3 §1.1 item 10. Report 2 §9.",
    ),
    (
        "3.1.11",
        "The GTEx judgment is approved. H3 uses GSE295758.",
        "Instruction 3 §1.1 item 11. Report 2 §9.",
    ),
    (
        "3.1.12",
        "Skipping the official bootstrap when all within-person correlations are at least 0.99 is approved. RNA-seq datasets follow instruction 3 §3.2.",
        "Instruction 3 §1.1 item 12. Report 2 §9.",
    ),
    (
        "3.2.GSE39445",
        "Role development. Technology microarray. Uses: H1, H2, H4 development.",
        "Instruction 3 §1.2.",
    ),
    (
        "3.2.GSE56931",
        "Role development. Technology microarray. Uses: H1, H2, H4 development.",
        "Instruction 3 §1.2.",
    ),
    (
        "3.2.GSE113883",
        "Role development, designated in this stage. Technology RNA-seq. Uses: year-scale H1, H2, H4 (DLMO) development, and the gate test.",
        "Instruction 3 §1.2.",
    ),
    (
        "3.2.GSE196728",
        "Role development, designated in this stage, conditional. Technology RNA-seq. Uses: sea-level condition H1 and H2 development, only if instruction 3 §2.3 passes.",
        "Instruction 3 §1.2.",
    ),
    (
        "3.2.GSE48113",
        "Role locked. Technology microarray. Uses: confirmation H1 and H4 (sleep phase misalignment). Expression values stay unread.",
        "Instruction 3 §1.2. Instruction 1 §1.3.",
    ),
    (
        "3.2.GSE122725",
        "Role locked. Technology microarray. Uses: confirmation H1 and H4 (schedule shift, light). Expression values stay unread.",
        "Instruction 3 §1.2. Instruction 1 §1.3.",
    ),
    (
        "3.2.GSE295758",
        "Role locked, designated in this stage. Technology RNA-seq. Uses: confirmation H1 (morning-evening) and H3 (age interaction). From this decision, instruction 1 §1.3 applies. Expression values and supplementary files are not opened.",
        "Instruction 3 §1.2.",
    ),
    (
        "3.2.GSE150337",
        "Role held. Technology RNA-seq. Not opened.",
        "Instruction 3 §1.2.",
    ),
    (
        "3.3",
        "RA_Peters is the provisional primary clock. The final choice is made at preregistration. This stage analyzes RA_GTEx, RA_Peters, and RA_default.",
        "Instruction 3 §1.3. Reason 1 stated in the instruction: coverage on all four platforms is at least 0.95 and is the highest of the three clocks. Coverages already computed in this project: GPL15331 280/289 = 0.9688581315 (results/stage1/rnaage_platform_coverage.tsv); GPL10379 278/289 = 0.9619377163 (same file); GPL24324 276/289 = 0.9550173010380623 (results/stage2/gpl24324_coverage.tsv). No fourth platform coverage row is in those files. Reason 2 stated in the instruction: the rhythm was strongest on development dataset GSE56931.",
    ),
]

header = "date\tid\tcontent\treason\ttiming\n"
existing = DECISIONS.read_text()
if not existing.endswith("\n"):
    existing += "\n"
    DECISIONS.write_text(existing)

present = set()
for line in existing.splitlines()[1:]:
    parts = line.split("\t")
    if len(parts) >= 2:
        present.add(parts[1])

new_lines = []
for item_id, content, reason in ROWS:
    if item_id in present:
        continue
    new_lines.append(f"{DATE}\t{item_id}\t{content}\t{reason}\tprior")

if new_lines:
    with DECISIONS.open("a") as handle:
        handle.write("\n".join(new_lines) + "\n")

role_path = OUT / "partA_decisions.tsv"
with role_path.open("w") as handle:
    handle.write(header)
    for item_id, content, reason in ROWS:
        handle.write(f"{DATE}\t{item_id}\t{content}\t{reason}\tprior\n")

print(f"appended {len(new_lines)} rows")
print(f"wrote {role_path}")
