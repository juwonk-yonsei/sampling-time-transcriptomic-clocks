#!/usr/bin/env python3
"""Extract GEO series-matrix metadata lines only.

For every series matrix, read line by line. Keep lines that start with '!'.
Stop at '!series_matrix_table_begin' and do not read further. Non-'!' lines
before that marker are not stored.
"""

from __future__ import annotations

import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data" / "processed" / "geo_headers"

SOURCES = [
    ROOT / "data" / "raw" / "GSE39445_series_matrix.txt.gz",
    ROOT / "data" / "raw" / "GSE56931_series_matrix.txt.gz",
    ROOT / "data" / "locked" / "GSE48113_series_matrix.txt.gz",
    ROOT / "data" / "locked" / "GSE122725_series_matrix.txt.gz",
]

STOP = "!series_matrix_table_begin"


def extract(path: Path, out_path: Path) -> dict:
    n_meta = 0
    n_blank = 0
    stopped = False
    non_bang = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        with out_path.open("w", encoding="utf-8") as out:
            for line in handle:
                if line.startswith("!"):
                    n_meta += 1
                    # The stop line is itself a '!' line. Record that it was
                    # reached, then stop before any expression row.
                    if line.startswith(STOP):
                        out.write(line if line.endswith("\n") else line + "\n")
                        stopped = True
                        break
                    out.write(line if line.endswith("\n") else line + "\n")
                elif line.strip() == "":
                    n_blank += 1
                else:
                    non_bang += 1
                    stopped = True
                    break
    return {
        "path": str(path.relative_to(ROOT)),
        "out": str(out_path.relative_to(ROOT)),
        "meta_lines": n_meta,
        "blank_lines": n_blank,
        "stopped_at_table": stopped and non_bang == 0,
        "non_bang_before_stop": non_bang,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in SOURCES:
        if not path.exists():
            raise SystemExit(f"missing {path}")
        gse = path.name.split("_")[0]
        out_path = OUT_DIR / f"{gse}_header.txt"
        rows.append(extract(path, out_path))
    report = ROOT / "results" / "stage1" / "header_extract_log.tsv"
    report.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    with report.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(keys) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[k]) for k in keys) + "\n")
    for row in rows:
        print(
            f"{row['path']}\tmeta_lines={row['meta_lines']}\t"
            f"stopped_at_table={row['stopped_at_table']}\t"
            f"non_bang_before_stop={row['non_bang_before_stop']}"
        )


if __name__ == "__main__":
    main()
