#!/usr/bin/env python3
"""Part A. Record instruction-6 approvals of report-5 judgments."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = ROOT / "config" / "decisions.tsv"
OUT = ROOT / "results" / "stage6"
OUT.mkdir(parents=True, exist_ok=True)
DATE = "2026-10-16"

ROWS = [
    (
        "6.1.1",
        "The human-data module list is the 14 multispecies mortality modules. The rodent module rescores stay under the TMrd_ prefix.",
        "Instruction 6 section 1 item 1. Report 5 section 9. Approved.",
        "prior",
    ),
    (
        "6.1.2",
        "Composite components use the 23 rodent module gene sets. That is the same set as the 1292-gene paper check.",
        "Instruction 6 section 1 item 2. Report 5 section 9. Approved.",
        "prior",
    ),
    (
        "6.1.3",
        "sigma is the REML residual SD of model M. sigma_u is the REML random-intercept SD of the null model, recorded when that SD has length 1.",
        "Instruction 6 section 1 item 3. Report 5 section 9. Approved.",
        "prior",
    ),
    (
        "6.1.4",
        "The 3-hour reproducible-rhythm rule is circular distance at most 3 hours.",
        "Instruction 6 section 1 item 4. Report 5 section 9. Approved.",
        "prior",
    ),
    (
        "6.1.5",
        "The TMrd_ clock names are approved.",
        "Instruction 6 section 1, deviation row. Report 5 section 10. Approved.",
        "prior",
    ),
]


def append_decisions() -> None:
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
    for item_id, content, reason, timing in ROWS:
        if item_id in present:
            continue
        new_lines.append(f"{DATE}\t{item_id}\t{content}\t{reason}\t{timing}")
    if new_lines:
        with DECISIONS.open("a") as handle:
            handle.write("\n".join(new_lines) + "\n")
    with (OUT / "partA_decisions.tsv").open("w") as handle:
        handle.write("date\tid\tcontent\treason\ttiming\n")
        for item_id, content, reason, timing in ROWS:
            handle.write(f"{DATE}\t{item_id}\t{content}\t{reason}\t{timing}\n")
    print(f"appended {len(new_lines)} decision rows")


if __name__ == "__main__":
    append_decisions()
