#!/usr/bin/env python3
"""Recompute SHA-256 for every file in the frozen manifest and compare with the stored value.

Each file is reported as one of
  match                      the recomputed hash equals the manifest hash
  mismatch                   the file exists but the hash differs
  missing (not redistributed) the file is not in this copy of the repository
For a mismatch inside a git checkout, the copy stored in the commit that added the
manifest is also hashed, and the later commits that changed the file are listed.

With --tage-dir, the tAge model files used for scoring (Zenodo 10.5281/zenodo.22166800,
version 2) are checked the same way against the md5 values in
results/stage4/module_clocks/zenodo_v2_files.tsv and, for files unpacked from the module
archives, the SHA-256 values in results/stage4/module_clocks/module_archive_members.tsv.

Usage:
  python scripts/verify_manifest.py
  python scripts/verify_manifest.py --tage-dir /path/to/zenodo_22166800_v2
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "prereg" / "manifest.tsv"
FREEZE_HASH = ROOT / "config" / "unlock_approved.txt"
ZENODO_FILES = ROOT / "results" / "stage4" / "module_clocks" / "zenodo_v2_files.tsv"
ARCHIVE_MEMBERS = ROOT / "results" / "stage4" / "module_clocks" / "module_archive_members.tsv"
MISSING = "missing (not redistributed)"


def digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def status(path: Path, expected: str, algorithm: str) -> str:
    if not path.is_file():
        return MISSING
    return "match" if digest(path, algorithm) == expected else "mismatch"


def git(*args: str) -> bytes | None:
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None


def check_manifest() -> dict[str, int]:
    counts = {"match": 0, "mismatch": 0, MISSING: 0}
    print("# frozen manifest: prereg/manifest.tsv")
    manifest_hash = digest(MANIFEST, "sha256")
    approved = FREEZE_HASH.read_text().strip() if FREEZE_HASH.is_file() else ""
    print(f"manifest_sha256\t{manifest_hash}\t{'equals' if manifest_hash == approved else 'differs from'} config/unlock_approved.txt")
    freeze = git("log", "--diff-filter=A", "--format=%H", "--", "prereg/manifest.tsv")
    freeze_commit = freeze.decode().split()[-1] if freeze else ""
    for row in read_tsv(MANIFEST):
        result = status(ROOT / row["path"], row["sha256"], "sha256")
        counts[result] += 1
        note = ""
        if result == "mismatch" and freeze_commit:
            frozen = git("show", f"{freeze_commit}:{row['path']}")
            same = frozen is not None and hashlib.sha256(frozen).hexdigest() == row["sha256"]
            later = git("log", "--format=%h %s", f"{freeze_commit}..HEAD", "--", row["path"])
            changes = "; ".join(later.decode().strip().splitlines()) if later else ""
            note = f"\tcopy at freeze commit {freeze_commit[:7]} {'matches' if same else 'does not match'}; changed after freeze in: {changes}"
        print(f"{result}\t{row['path']}{note}")
    return counts


def check_tage(tage_dir: Path) -> dict[str, int]:
    counts = {"match": 0, "mismatch": 0, MISSING: 0}
    print(f"# tAge model files: {tage_dir}")
    used = [r for r in read_tsv(ZENODO_FILES) if r["used_for_scoring"] == "True"]
    for row in used:
        result = status(tage_dir / row["name"], row["md5"], "md5")
        counts[result] += 1
        print(f"{result}\tmd5\t{row['name']}")
    used_archives = {r["name"] for r in used if r["name"].endswith(".zip")}
    for row in read_tsv(ARCHIVE_MEMBERS):
        if row["archive"] not in used_archives:
            continue
        folder = tage_dir / row["archive"].removesuffix(".zip")
        result = status(folder / row["name"], row["sha256"], "sha256")
        counts[result] += 1
        print(f"{result}\tsha256\t{folder.name}/{row['name']}")
    return counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tage-dir", type=Path, help="folder with the Zenodo 22166800 v2 files and the unpacked module archives")
    args = parser.parse_args()
    summaries = [("manifest", check_manifest())]
    if args.tage_dir:
        summaries.append(("tage", check_tage(args.tage_dir)))
    print("# summary")
    for name, counts in summaries:
        print(f"{name}\t" + "\t".join(f"{k}={v}" for k, v in counts.items()))
    return 1 if any(c["mismatch"] for _, c in summaries) else 0


if __name__ == "__main__":
    sys.exit(main())
