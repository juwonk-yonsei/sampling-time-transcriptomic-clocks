#!/usr/bin/env python3
"""Record byte size and SHA-256 for downloaded raw files.

Locked supplementary files are not downloaded. Their names and sizes come
from HTTP Content-Length.
"""

from __future__ import annotations

import hashlib
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DOWNLOADED = [
    ROOT / "data" / "raw" / "GSE39445_series_matrix.txt.gz",
    ROOT / "data" / "raw" / "GSE56931_series_matrix.txt.gz",
    ROOT / "data" / "locked" / "GSE48113_series_matrix.txt.gz",
    ROOT / "data" / "locked" / "GSE122725_series_matrix.txt.gz",
    ROOT / "data" / "raw" / "gpl" / "GPL10379.annot.gz",
    ROOT / "data" / "raw" / "gpl" / "GPL15331_platform.txt",
    ROOT / "data" / "raw" / "gpl" / "GPL24324_platform.txt",
]

LOCKED_SUPPL = [
    ("GSE48113", "GSE48113_RAW.tar", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48113/suppl/GSE48113_RAW.tar"),
    ("GSE48113", "filelist.txt", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48113/suppl/filelist.txt"),
    ("GSE122725", "GSE122725_RAW.tar", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122725/suppl/GSE122725_RAW.tar"),
    ("GSE122725", "filelist.txt", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122725/suppl/filelist.txt"),
]

DEV_SUPPL_NOT_DOWNLOADED = [
    ("GSE39445", "GSE39445_RAW.tar", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39445/suppl/GSE39445_RAW.tar"),
    ("GSE39445", "filelist.txt", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39445/suppl/filelist.txt"),
    ("GSE56931", "GSE56931_RAW.tar", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56931/suppl/GSE56931_RAW.tar"),
    ("GSE56931", "filelist.txt", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56931/suppl/filelist.txt"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def content_length(url: str) -> str:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.headers.get("Content-Length", "NA")


def main() -> None:
    lines = ["path\tbytes\tsha256"]
    for path in DOWNLOADED:
        if not path.exists():
            raise SystemExit(f"missing {path}")
        rel = path.relative_to(ROOT).as_posix()
        lines.append(f"{rel}\t{path.stat().st_size}\t{sha256(path)}")
    out = ROOT / "config" / "raw_checksums.tsv"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(out)

    locked_lines = ["dataset\tname\tbytes\tdownloaded"]
    for dataset, name, url in LOCKED_SUPPL:
        locked_lines.append(f"{dataset}\t{name}\t{content_length(url)}\tno")
    locked_path = ROOT / "results" / "stage1" / "locked_suppl_listing.tsv"
    locked_path.write_text("\n".join(locked_lines) + "\n", encoding="utf-8")
    print(locked_path)

    dev_lines = ["dataset\tname\tbytes\tdownloaded\treason"]
    reason = "series matrix contains the processed expression table"
    for dataset, name, url in DEV_SUPPL_NOT_DOWNLOADED:
        dev_lines.append(f"{dataset}\t{name}\t{content_length(url)}\tno\t{reason}")
    dev_path = ROOT / "results" / "stage1" / "development_suppl_not_downloaded.tsv"
    dev_path.write_text("\n".join(dev_lines) + "\n", encoding="utf-8")
    print(dev_path)


if __name__ == "__main__":
    main()
