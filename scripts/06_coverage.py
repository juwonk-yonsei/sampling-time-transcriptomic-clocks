#!/usr/bin/env python3
"""Platform inclusion of RNAAgeCalc blood signatures and MCPcounter markers.

EPIC is not used: its LICENSE is academic non-commercial and forbids
redistribution. GPL24324 has no gene-symbol column.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage1"
BLANK = {"", "NA", "N/A", "NULL", "NONE", "---"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tokens(raw: str) -> list[str]:
    parts = []
    for piece in raw.split("///"):
        piece = piece.strip().strip('"')
        if piece and piece.upper() not in BLANK:
            parts.append(piece)
    return list(dict.fromkeys(parts))


def symbols_from_table(path: Path, column: str, gzipped: bool) -> set[str]:
    opener = gzip.open if gzipped else open
    found: set[str] = set()
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
            found.update(tokens(parts[idx]))
    return found


def main() -> None:
    gene_path = ROOT / "data" / "processed" / "rnaage_genes" / "blood_signature_genes.tsv"
    by_clock: dict[str, set[str]] = {}
    with gene_path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_clock.setdefault(row["clock"], set()).add(row["gene"])

    platforms = {
        "GPL15331": symbols_from_table(ROOT / "data" / "raw" / "gpl" / "GPL15331_platform.txt", "GENE_SYMBOL", False),
        "GPL10379": symbols_from_table(ROOT / "data" / "raw" / "gpl" / "GPL10379.annot.gz", "Gene symbol", True),
    }
    platform_note = {
        "GPL15331": "column GENE_SYMBOL; symbol present if it is a /// token",
        "GPL10379": "column Gene symbol; symbol present if it is a /// token",
        "GPL24324": "no gene-symbol column; second SPOT_ID holds free text; coverage not computed",
    }

    rows = []
    for clock, genes in by_clock.items():
        n = len(genes)
        for gpl in ("GPL15331", "GPL10379", "GPL24324"):
            if gpl == "GPL24324":
                rows.append({
                    "reference": clock,
                    "platform": gpl,
                    "n_reference_genes": n,
                    "n_on_platform": "NA",
                    "fraction": "NA",
                    "note": platform_note[gpl],
                })
                continue
            hit = len(genes & platforms[gpl])
            rows.append({
                "reference": clock,
                "platform": gpl,
                "n_reference_genes": n,
                "n_on_platform": hit,
                "fraction": f"{hit / n:.10g}" if n else "NA",
                "note": platform_note[gpl],
            })

    mcp_path = ROOT / "data" / "processed" / "ref" / "MCPcounter_genes.txt"
    mcp_desc = Path("/ssd240/anaconda3/envs/r_env/lib/R/library/MCPcounter/DESCRIPTION")
    populations: dict[str, set[str]] = {}
    with mcp_path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            symbol = row["HUGO symbols"].strip().strip('"')
            pop = row["Cell population"].strip().strip('"')
            if symbol:
                populations.setdefault(pop, set()).add(symbol)
    all_markers = set().union(*populations.values())
    mcp_rows = []
    for pop, genes in sorted(populations.items()):
        for gpl in ("GPL15331", "GPL10379", "GPL24324"):
            if gpl == "GPL24324":
                mcp_rows.append({
                    "population": pop,
                    "platform": gpl,
                    "n_markers": len(genes),
                    "n_on_platform": "NA",
                    "fraction": "NA",
                })
                continue
            hit = len(genes & platforms[gpl])
            mcp_rows.append({
                "population": pop,
                "platform": gpl,
                "n_markers": len(genes),
                "n_on_platform": hit,
                "fraction": f"{hit / len(genes):.10g}",
            })
    for gpl in ("GPL15331", "GPL10379", "GPL24324"):
        if gpl == "GPL24324":
            mcp_rows.append({
                "population": "ALL",
                "platform": gpl,
                "n_markers": len(all_markers),
                "n_on_platform": "NA",
                "fraction": "NA",
            })
            continue
        hit = len(all_markers & platforms[gpl])
        mcp_rows.append({
            "population": "ALL",
            "platform": gpl,
            "n_markers": len(all_markers),
            "n_on_platform": hit,
            "fraction": f"{hit / len(all_markers):.10g}",
        })

    def write(path: Path, records: list[dict]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)

    write(OUT / "rnaage_platform_coverage.tsv", rows)
    write(OUT / "mcpcounter_platform_coverage.tsv", mcp_rows)
    meta = [
        f"mcpcounter_package\tMCPcounter",
        f"mcpcounter_version\t1.2.0",
        f"mcpcounter_license\tGPL-3",
        f"mcpcounter_description_sha256\t{sha256(mcp_desc)}",
        f"mcpcounter_description_bytes\t{mcp_desc.stat().st_size}",
        f"mcpcounter_genes_path\tdata/processed/ref/MCPcounter_genes.txt",
        f"mcpcounter_genes_sha256\t{sha256(mcp_path)}",
        f"mcpcounter_genes_bytes\t{mcp_path.stat().st_size}",
        f"mcpcounter_signature_license\tCC0 1.0",
        f"mcpcounter_n_markers\t{len(all_markers)}",
        f"mcpcounter_n_populations\t{len(populations)}",
        "epic\tnot_downloaded",
        "epic_reason\tLICENSE is academic non-commercial and forbids redistribution",
        "gpl24324\tno gene-symbol column; coverage not computed",
    ]
    (OUT / "cell_reference_files.tsv").write_text("\n".join(meta) + "\n", encoding="utf-8")
    print("coverage rows", len(rows), "mcp rows", len(mcp_rows), "markers", len(all_markers))


if __name__ == "__main__":
    main()
