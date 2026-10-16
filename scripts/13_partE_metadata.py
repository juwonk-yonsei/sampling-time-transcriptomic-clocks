#!/usr/bin/env python3
"""Part E. Metadata-only search log and header facts. No expression values."""

import csv
import json
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
OUT = ROOT / "results" / "stage2"
HEAD = OUT / "rnaseq_headers"
DATE = "2026-10-16"


def get_bytes(url, timeout=60, headers=None):
    req = urllib.request.Request(url, headers=headers or {"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.status, resp.read()


def esearch(db, term, retmax=0):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "term": term, "retmode": "json", "retmax": retmax}
    )
    for attempt in range(5):
        try:
            _, raw = get_bytes(url)
            data = json.loads(raw)
            return data["esearchresult"]
        except Exception as exc:
            time.sleep(2 + attempt)
            last = exc
    raise last


def split_geo(line):
    return [part.strip().strip('"') for part in line.rstrip("\n").split("\t")]


def series_lines(acc):
    path = HEAD / f"{acc}_header.txt"
    buckets = {}
    for line in path.read_text(errors="replace").splitlines():
        parts = split_geo(line)
        buckets.setdefault(parts[0], []).append(parts[1:])
    return buckets


def first(buckets, key):
    rows = buckets.get(key) or []
    if not rows or not rows[0]:
        return ""
    return rows[0][0]


def char_map(buckets):
    found = {}
    for row in buckets.get("!Sample_characteristics_ch1", []):
        if not row or ":" not in row[0]:
            continue
        name = row[0].split(":", 1)[0].strip()
        if name not in found:
            found[name] = [part.split(":", 1)[1].strip() if ":" in part else part for part in row]
    return found


def unique_processing(buckets):
    seen = []
    for row in buckets.get("!Sample_data_processing", []):
        if row and row[0] not in seen:
            seen.append(row[0])
    return seen


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def header_facts():
    rows = []
    for path in sorted(HEAD.glob("*_header.txt")):
        acc = path.name.split("_")[0]
        buckets = series_lines(acc)
        chars = char_map(buckets)
        titles = (buckets.get("!Sample_title") or [[]])[0]
        geos = (buckets.get("!Sample_geo_accession") or [[]])[0]
        proc = unique_processing(buckets)
        rows.append({
            "accession": acc,
            "n_sample_geo_accession": len([g for g in geos if g]),
            "n_sample_title": len([t for t in titles if t]),
            "title": first(buckets, "!Series_title"),
            "organism": "|".join(sorted({v for row in buckets.get("!Sample_organism_ch1", []) for v in row if v})),
            "library_strategy": "|".join(sorted({v for row in buckets.get("!Sample_library_strategy", []) for v in row if v})),
            "library_source": "|".join(sorted({v for row in buckets.get("!Sample_library_source", []) for v in row if v})),
            "source_name": "|".join(sorted({v for row in buckets.get("!Sample_source_name_ch1", []) for v in row if v})),
            "platform_id": first(buckets, "!Series_platform_id"),
            "pubmed_id": first(buckets, "!Series_pubmed_id"),
            "characteristic_names": "|".join(chars.keys()),
            "overall_design": first(buckets, "!Series_overall_design"),
            "processing_lines": " || ".join(proc),
            "title_examples": " | ".join(titles[:3]),
        })
        for name, vals in chars.items():
            counts = Counter(vals)
            rows[-1][f"char_n::{name}"] = len(vals)
            rows[-1][f"char_unique::{name}"] = len(counts)
            rows[-1][f"char_top::{name}"] = "; ".join(f"{k}:{v}" for k, v in counts.most_common(8))
    # Flatten only the stable columns plus a long JSON of char stats.
    flat = []
    for row in rows:
        extra = {k: v for k, v in row.items() if k.startswith("char_")}
        base = {k: v for k, v in row.items() if not k.startswith("char_")}
        base["characteristic_stats"] = json.dumps(extra, ensure_ascii=False)
        flat.append(base)
    fields = [
        "accession", "n_sample_geo_accession", "n_sample_title", "title", "organism",
        "library_strategy", "library_source", "source_name", "platform_id", "pubmed_id",
        "characteristic_names", "title_examples", "overall_design", "processing_lines",
        "characteristic_stats",
    ]
    write_tsv(OUT / "partE_header_facts.tsv", flat, fields)
    print("wrote header facts", len(flat))


def main():
    search_rows = []
    queries = [
        ("GEO", "gds", 'human[Organism] AND (blood[All Fields] OR PBMC[All Fields] OR leukocyte[All Fields]) AND ("RNA-seq"[All Fields] OR "RNA sequencing"[All Fields]) AND (circadian[All Fields] OR "time of day"[All Fields] OR "clock time"[All Fields])'),
        ("GEO", "gds", '"whole blood"[All Fields] AND (RNA-seq[All Fields] OR "high throughput sequencing"[All Fields]) AND (circadian[All Fields] OR "time of day"[All Fields] OR "time-of-day"[All Fields]) AND human[Organism]'),
        ("GEO", "gds", 'PBMC[All Fields] AND "RNA-seq"[All Fields] AND (circadian[All Fields] OR "time of day"[All Fields]) AND human[Organism]'),
        ("GEO", "gds", "TimeSignature[All Fields] AND human[Organism]"),
        ("GEO", "gds", "BodyTime[All Fields] AND human[Organism]"),
        ("PubMed", "pubmed", "TimeSignature blood RNA-seq"),
        ("PubMed", "pubmed", "TimeSignature circadian blood"),
        ("PubMed", "pubmed", "Wittenbrink TimeSignature"),
        ("PubMed", "pubmed", "BodyTime blood transcriptome circadian"),
        ("PubMed", "pubmed", "BodyTime circadian transcriptome"),
        ("PubMed", "pubmed", '"whole blood" RNA-seq "time of day" circadian'),
        ("PubMed", "pubmed", "blood RNA-seq circadian older adults"),
    ]
    for database, db, term in queries:
        time.sleep(0.4)
        result = esearch(db, term, retmax=5 if db == "pubmed" else 0)
        search_rows.append({
            "date": DATE,
            "database": database,
            "query": term,
            "hit_count": result["count"],
            "id_list": ";".join(result.get("idlist", [])),
            "note": "",
        })
        print(database, result["count"], term[:70])

    biostudies = [
        "blood AND RNA-seq AND circadian",
        '"whole blood" AND RNA-seq AND "time of day"',
        "PBMC AND RNA-seq AND circadian",
        "leukocyte AND RNA-seq AND circadian",
        "BodyTime AND blood",
        "TimeSignature AND blood",
    ]
    for query in biostudies:
        time.sleep(0.3)
        url = "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urllib.parse.urlencode(
            {"query": query, "pageSize": 1}
        )
        _, raw = get_bytes(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
        total = json.loads(raw).get("totalHits")
        search_rows.append({
            "date": DATE,
            "database": "BioStudies",
            "query": query,
            "hit_count": total,
            "id_list": "",
            "note": "",
        })
        print("BioStudies", total, query)

    for query in ["blood RNA-seq circadian", "whole blood time of day"]:
        url = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?" + urllib.parse.urlencode(
            {"keywords": query}
        )
        note = ""
        count = ""
        try:
            status, raw = get_bytes(url)
            note = f"HTTP {status}; body_start={raw[:80]!r}"
        except Exception as exc:
            note = f"{type(exc).__name__}: {exc}"
        search_rows.append({
            "date": DATE,
            "database": "ArrayExpress",
            "query": query,
            "hit_count": count,
            "id_list": "",
            "note": note,
        })
        print("ArrayExpress", note[:120])

    write_tsv(
        OUT / "partE_search_log.tsv",
        search_rows,
        ["date", "database", "query", "hit_count", "id_list", "note"],
    )
    print("wrote search log", len(search_rows))
    header_facts()


if __name__ == "__main__":
    main()
