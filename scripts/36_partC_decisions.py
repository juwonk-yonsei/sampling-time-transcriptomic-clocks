#!/usr/bin/env python3
"""Part C. Record implementation judgments that do not change frozen tables."""

from pathlib import Path

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
DECISIONS = ROOT / "config" / "decisions.tsv"
DATE = "2026-10-16"
ROWS = [
    (
        "6.2.1",
        "The L2 gate window of 18-26 h after habitual wake, modulo 24, is peak_mod in [0, 24) with peak_mod >= 18 or peak_mod <= 2, inclusive.",
        "Preregistration section 6. 26 modulo 24 is 2. The sentence is not edited.",
        "prior",
    ),
    (
        "6.2.2",
        "GENCODE 38 transcript and gene identifiers are joined after removing the version suffix following the first dot. The versioned identifiers stay in the table.",
        "Preregistration section 5 names GENCODE release 38 and does not state the version suffix. The join key is recorded in gencode_v38_basic_transcripts.tsv.",
        "prior",
    ),
    (
        "6.2.3",
        "L2 probesets with more than one SYMBOL are dropped. The dropped count is written to the preprocess audit. L1 already drops multi-symbol probes. The L2 sentence does not state that rule.",
        "Preregistration section 5. Conservative reading of a one-symbol map. The sentence is not edited.",
        "prior",
    ),
    (
        "6.2.4",
        "Holm correction is applied to the hypothesis tests that the run actually performs. A rehearsal is one dataset. A failed gate removes that dataset's tests.",
        "Preregistration section 7 states families of 9, 6, 3 and 6 tests across the locked datasets. A run that has not unlocked those datasets cannot include their tests.",
        "prior",
    ),
    (
        "6.2.5",
        "The L2 confirmation time axis is wake_ref_h modulo 24.",
        "Preregistration section 3. Stage 5 power used the unwrapped wake_ref_h. Frozen table F6 is not recomputed.",
        "prior",
    ),
    (
        "6.2.6",
        "Confirmation component scores are partial sums of the weight score over the 23 rodent module gene sets plus outside. Table F3 keeps the stage-5 values.",
        "Preregistration section 4. Within-person, the stage-5 pipeline terms and the weight score differ by a per-sample constant, so the frozen cosinor totals are unchanged.",
        "prior",
    ),
    (
        "6.2.7",
        "The scoring reference is every analysed sample of that dataset.",
        "Preregistration section 4. Stage 5 used control samples as the reference. Confirmation follows the preregistration.",
        "prior",
    ),
    (
        "6.2.8",
        "MCPcounter on the RNA-seq path uses log2(FPKM+1), the same matrix as the RNAAgeCalc linear score. Gate genes use log2(CPM+1).",
        "Preregistration section 5 specifies log2(FPKM+1) for the RA score and log2(CPM+1) for gate genes. It does not name the MCPcounter input for counts.",
        "prior",
    ),
    (
        "6.2.9",
        "confirm_root walks parent directories until it finds scripts/lib_cosinor.R.",
        "The lock and pipeline functions are sourced from more than one entry script.",
        "prior",
    ),
    (
        "6.2.10",
        "Decision 6.2.3 is withdrawn. L2 keeps every SYMBOL annotation and selects the highest-mean probeset per symbol. A probeset with more than one SYMBOL can be selected for each of those symbols. The count of multi-symbol probesets is still written to the audit.",
        "clariomshumanhttranscriptcluster.db has 24351 PROBEID keys. PER1 has one probeset and that probeset has more than one SYMBOL. Dropping it would remove PER1, so the L2 gate in preregistration section 6 could not be applied. The preregistration sentence says highest-mean probeset per symbol and does not say to drop multi-symbol probesets.",
        "prior",
    ),
]


def main() -> None:
    text = DECISIONS.read_text()
    extra = []
    for row in ROWS:
        if f"\t{row[0]}\t" in text:
            continue
        extra.append("\t".join((DATE, *row)))
    if extra:
        with DECISIONS.open("a") as handle:
            if not text.endswith("\n"):
                handle.write("\n")
            handle.write("\n".join(extra) + "\n")
    print(f"added {len(extra)}")


if __name__ == "__main__":
    main()
