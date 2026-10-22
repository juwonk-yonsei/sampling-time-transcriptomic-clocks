#!/usr/bin/env python3
"""Check that no locked-dataset expression values entered git before the unlock commit.

Locked datasets: GSE48113, GSE122725, GSE295758. Their sample identifiers are
read from the headers of the stored confirmation matrices. Every blob in every
commit is searched for those identifiers. A path whose content holds at least
10 identifiers of one dataset is listed with the first commit that added it.
Every tree up to and including the unlock commit is also searched for any
identifier at all.

Usage: python scripts/check_lock_order.py [--out audit/locked_expression_first_commits.tsv]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import sys
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADERS = {
    "GSE48113": "results/confirm/L1/log2_genes.tsv.gz",
    "GSE122725": "results/confirm/L2/log2_genes.tsv",
    "GSE295758": "results/confirm/L3/gene_counts.tsv",
}
MIN_IDS = 10
GSM_RE = re.compile(rb"GSM\d{6,8}")
L3_RE = re.compile(rb"(?<![A-Za-z0-9_])X?(\d{1,2}_[DN](?:_batch\d)?)(?![A-Za-z0-9])")
METADATA_NAMES = ("samples_", "wake_ref", "_header.txt", "bang_lines", "column_names", "layout")
# Rehearsal R2 is N(0, 1) noise (seed 20261015, scripts/43_partD_pipeline.R run_r2)
# that only borrows the GSE122725 sample identifiers from the public annotation.
SYNTHETIC_PREFIX = "results/stage6/rehearsal/"
DATE = "--date=format:%Y-%m-%d %H:%M:%S %z"


def git(*args: str, data: bytes | None = None) -> str:
    out = subprocess.run(["git", "-c", "core.quotePath=false", *args], cwd=ROOT, input=data,
                         capture_output=True, check=True)
    return out.stdout.decode("utf-8")


def header_ids(path: str) -> set[str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(ROOT / path, "rt") as handle:
        return set(handle.readline().rstrip("\n").split("\t")[1:])


def tree_blobs(commits: list[str]) -> dict[str, set[str]]:
    blobs: dict[str, set[str]] = {}
    for commit in commits:
        for entry in git("ls-tree", "-r", "-z", commit).split("\0"):
            if not entry:
                continue
            meta, path = entry.split("\t", 1)
            _, kind, sha = meta.split()
            if kind == "blob":
                blobs.setdefault(sha, set()).add(path)
    return blobs


TRUNCATED: dict[str, set[str]] = {}


def gunzip(data: bytes, sha: str = "", paths: set[str] | None = None) -> bytes:
    """Decompress gzip data. A truncated stream is decompressed as far as it goes."""
    try:
        return gzip.decompress(data)
    except (OSError, EOFError):
        try:
            partial = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(data)
        except zlib.error:
            return data
        if sha:
            TRUNCATED[sha] = set(paths or ())
        return partial


def count_ids(blobs: dict[str, set[str]], gsm: dict[str, str], l3: set[str]) -> dict[str, dict[str, int]]:
    proc = subprocess.Popen(["git", "cat-file", "--batch"], cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    counts: dict[str, dict[str, int]] = {}
    for sha, paths in blobs.items():
        proc.stdin.write((sha + "\n").encode())
        proc.stdin.flush()
        size = int(proc.stdout.readline().split()[2])
        data = proc.stdout.read(size)
        proc.stdout.read(1)
        if any(p.endswith(".gz") for p in paths):
            data = gunzip(data, sha, paths)
        found: dict[str, set[bytes]] = {ds: set() for ds in HEADERS}
        for match in GSM_RE.findall(data):
            ds = gsm.get(match.decode())
            if ds:
                found[ds].add(match)
        for match in L3_RE.findall(data):
            if match.decode() in l3:
                found["GSE295758"].add(match)
        hit = {ds: len(v) for ds, v in found.items() if v}
        if hit:
            counts[sha] = hit
    proc.stdin.close()
    proc.wait()
    return counts


def first_add(path: str, follow: bool) -> tuple[str, str]:
    args = ["log", "--all", "--diff-filter=A", DATE, "--format=%H\t%ad"]
    if follow:
        args.append("--follow")
    lines = git(*args, "--", path).strip().splitlines()
    if not lines:
        return "", ""
    sha, when = lines[-1].split("\t")
    return sha, when


def kind_of(path: str, first_line: str) -> str:
    if path.startswith("reports/"):
        return "report"
    if path.startswith(SYNTHETIC_PREFIX):
        return "synthetic_rehearsal"
    if first_line.split("\t", 1)[0] in ("gene", "gene_name", "transcript_id"):
        return "expression_matrix"
    if any(name in path for name in METADATA_NAMES):
        return "sample_metadata"
    return "derived_per_sample_values"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="audit/locked_expression_first_commits.tsv")
    args = parser.parse_args()

    ids = {ds: header_ids(path) for ds, path in HEADERS.items()}
    gsm = {s: ds for ds in ("GSE48113", "GSE122725") for s in ids[ds]}
    l3 = ids["GSE295758"]

    unlock_lines = git("log", "--diff-filter=A", DATE, "--format=%H\t%ad\t%s", "--", "config/unlock_approved.txt").strip().splitlines()
    unlock_sha, unlock_time, unlock_msg = unlock_lines[-1].split("\t")
    unlock_epoch = int(git("show", "-s", "--format=%ct", unlock_sha).strip())

    all_commits = git("rev-list", "--all").split()
    blobs = tree_blobs(all_commits)
    counts = count_ids(blobs, gsm, l3)

    rows = []
    for sha, hit in counts.items():
        strong = {ds: n for ds, n in hit.items() if n >= MIN_IDS}
        if not strong:
            continue
        for path in sorted(blobs[sha]):
            rows.append((path, sha, strong))
    by_path: dict[str, dict] = {}
    for path, sha, strong in rows:
        entry = by_path.setdefault(path, {"blobs": set(), "datasets": {}})
        entry["blobs"].add(sha)
        for ds, n in strong.items():
            entry["datasets"][ds] = max(n, entry["datasets"].get(ds, 0))

    out_rows = []
    for path in sorted(by_path):
        sha_any = sorted(by_path[path]["blobs"])[0]
        raw = subprocess.run(["git", "cat-file", "-p", sha_any], cwd=ROOT, capture_output=True, check=True).stdout
        if path.endswith(".gz"):
            raw = gunzip(raw)
        first_line = raw.split(b"\n", 1)[0].decode("utf-8", "replace")
        add_sha, add_time = first_add(path, follow=False)
        follow_sha, follow_time = first_add(path, follow=True)
        add_epoch = int(git("show", "-s", "--format=%ct", add_sha).strip())
        out_rows.append({
            "path": path,
            "kind": kind_of(path, first_line),
            "datasets": ";".join(f"{ds}:{n}" for ds, n in sorted(by_path[path]["datasets"].items())),
            "first_commit": add_sha,
            "first_commit_time_kst": add_time,
            "after_unlock": str(add_epoch > unlock_epoch),
            "follow_commit": follow_sha,
            "follow_time_kst": follow_time,
        })

    pre_commits = git("rev-list", unlock_sha).split()
    pre_blobs = tree_blobs(pre_commits)
    pre_counts = count_ids(pre_blobs, gsm, l3)

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    values = [r for r in out_rows if r["kind"] in ("expression_matrix", "derived_per_sample_values")]
    print(f"unlock_commit\t{unlock_sha}\t{unlock_time}\t{unlock_msg}")
    print(f"commits_scanned\t{len(all_commits)}\tblobs\t{len(blobs)}\tpaths_with_locked_ids\t{len(out_rows)}")
    for kind in ("expression_matrix", "derived_per_sample_values", "synthetic_rehearsal", "sample_metadata", "report"):
        sel = [r for r in out_rows if r["kind"] == kind]
        if sel:
            earliest = min(sel, key=lambda r: r["first_commit_time_kst"])
            print(f"kind\t{kind}\tfiles\t{len(sel)}\tearliest_first_commit\t{earliest['first_commit_time_kst']}\t{earliest['first_commit'][:7]}\tafter_unlock_all\t{all(r['after_unlock'] == 'True' for r in sel)}")
    print(f"preunlock_commits\t{len(pre_commits)}\tblobs\t{len(pre_blobs)}\tblobs_with_any_locked_id\t{len(pre_counts)}")
    for sha, hit in sorted(pre_counts.items(), key=lambda kv: sorted(pre_blobs[kv[0]])):
        print(f"preunlock_blob\t{sha[:12]}\t{';'.join(sorted(pre_blobs[sha]))}\t{hit}")
    for sha, paths in sorted(TRUNCATED.items()):
        print(f"truncated_gzip_blob\t{sha[:12]}\t{';'.join(sorted(paths))}")
    value_blobs = {sha for r in values for sha in by_path[r["path"]]["blobs"]}
    leaked = value_blobs & set(pre_blobs)
    print(f"value_blobs_in_preunlock_trees\t{len(leaked)}")
    ok = all(r["after_unlock"] == "True" for r in values) and not leaked
    print(f"verdict\t{'all after unlock' if ok else 'files before unlock'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
