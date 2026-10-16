#!/usr/bin/env python3
"""Build a GENCODE release 38 basic transcript-to-gene table.

Identifier versions after the first dot are kept in transcript_id and gene_id.
transcript_id_base and gene_id_base drop that suffix for joins.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

ATTR = re.compile(r'(transcript_id|gene_id|gene_name) "([^"]+)"')


def base_id(value: str) -> str:
    return value.split(".", 1)[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    n = 0
    with gzip.open(args.gtf, "rt") as src, Path(args.out).open("w", encoding="utf-8") as dst:
        dst.write("transcript_id\tgene_id\tgene_name\ttranscript_id_base\tgene_id_base\n")
        for line in src:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            attrs = dict(ATTR.findall(parts[8]))
            transcript_id = attrs["transcript_id"]
            gene_id = attrs["gene_id"]
            gene_name = attrs.get("gene_name", "")
            dst.write(
                f"{transcript_id}\t{gene_id}\t{gene_name}\t{base_id(transcript_id)}\t{base_id(gene_id)}\n"
            )
            n += 1
    print(f"transcripts {n}")


if __name__ == "__main__":
    main()
