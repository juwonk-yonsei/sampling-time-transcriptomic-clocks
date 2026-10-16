#!/usr/bin/env python3
"""Part B. Confirm Zenodo 22166800 version 2 and the local clock files.

Clock coefficients and gene lists are not written.
"""

import hashlib
import json
import os
import subprocess
import warnings
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "stage4" / "module_clocks"
OUT.mkdir(parents=True, exist_ok=True)
DEST = Path("/HDD8T1/WORK/GeroSci/zenodo_22166800_v2")
DEST.mkdir(parents=True, exist_ok=True)
ZENODO_JSON = Path("/tmp/zenodo_rec.json")
ZENODO_VERSIONS = Path("/tmp/zenodo_versions.json")

# Elastic-net multitissue scaled-diff composites and the module archives.
# The other  composite pkl files are listed but not downloaded (see the report).
NEEDED = [
    "EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
    "EN_Mortality_Rodents_Multitissue_scaleddiff.pkl",
    "EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl",
    "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
    "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl",
    "Rodent module clocks.zip",
    "Multispecies module clocks.zip",
]
LOCAL_CANDIDATES = {
    "EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl": [
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_zenodo/EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl",
    ],
    "EN_Mortality_Rodents_Multitissue_scaleddiff.pkl": [
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_Mortality_Rodents_Multitissue_scaleddiff.pkl",
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_zenodo/EN_Mortality_Rodents_Multitissue_scaleddiff.pkl",
    ],
    "EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl": [
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl",
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_zenodo/EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl",
    ],
    "EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl": [
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl",
    ],
    "EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl": [
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl",
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_zenodo/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl",
    ],
    "Rodent module clocks.zip": [
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_modules/rodent_modules.zip",
    ],
    "Multispecies module clocks.zip": [
        "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_modules/multispecies_modules.zip",
    ],
}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify(name: str) -> str:
    if "module clocks" in name and name.endswith(".zip"):
        return "module_clock"
    if name.endswith(".pkl"):
        return "composite_clock"
    return "other"


def main() -> None:
    rec = json.loads(ZENODO_JSON.read_text())
    versions = json.loads(ZENODO_VERSIONS.read_text())
    hits = versions["hits"]["hits"]
    version_rows = []
    for index, hit in enumerate(sorted(hits, key=lambda h: h["metadata"]["publication_date"]), start=1):
        version_rows.append(
            {
                "order_by_publication_date": index,
                "record_id": hit["id"],
                "doi": hit.get("doi"),
                "publication_date": hit["metadata"].get("publication_date"),
                "created": hit.get("created"),
                "version_label": hit["metadata"].get("version") or "",
                "is_2026_08_29": hit["metadata"].get("publication_date") == "2026-08-29",
            }
        )
    pd.DataFrame(version_rows).to_csv(OUT / "zenodo_versions.tsv", sep="\t", index=False)

    files = []
    for item in rec["files"]:
        files.append(
            {
                "record_id": rec["id"],
                "publication_date": rec["metadata"]["publication_date"],
                "version_label": rec["metadata"].get("version") or "",
                "name": item["key"],
                "size": item["size"],
                "md5": item["checksum"].split(":", 1)[1],
                "class": classify(item["key"]),
                "used_for_scoring": item["key"] in NEEDED,
            }
        )
    file_table = pd.DataFrame(files)
    file_table.to_csv(OUT / "zenodo_v2_files.tsv", sep="\t", index=False)

    license_row = pd.DataFrame(
        [
            {
                "record_id": rec["id"],
                "conceptrecid": rec["conceptrecid"],
                "publication_date": rec["metadata"]["publication_date"],
                "version_label": rec["metadata"].get("version") or "",
                "access_right": rec["metadata"].get("access_right"),
                "license_field": rec["metadata"].get("license") or "",
                "description_repository": "https://github.com/Gladyshev-Lab/tAge",
                "local_license_title": "MGB OPEN ACCESS LICENSE 1.0",
                "local_license_path": "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tAge/LICENSE",
                "n_versions": len(version_rows),
                "newer_than_2026_08_29": any(
                    row["publication_date"] > "2026-08-29" for row in version_rows
                ),
            }
        ]
    )
    license_row.to_csv(OUT / "zenodo_record.tsv", sep="\t", index=False)

    by_name = {row["name"]: row for row in files}
    local_rows = []
    accepted = {}
    for name, paths in LOCAL_CANDIDATES.items():
        zenodo = by_name[name]
        for raw in paths:
            path = Path(raw)
            if not path.exists():
                local_rows.append(
                    {
                        "zenodo_name": name,
                        "local_path": raw,
                        "exists": False,
                        "size": "",
                        "md5": "",
                        "name_match": path.name == name,
                        "size_match": False,
                        "md5_match": False,
                        "accepted": False,
                    }
                )
                continue
            digest = md5(path)
            name_match = path.name == name
            size_match = path.stat().st_size == zenodo["size"]
            md5_match = digest == zenodo["md5"]
            ok = name_match and size_match and md5_match
            local_rows.append(
                {
                    "zenodo_name": name,
                    "local_path": raw,
                    "exists": True,
                    "size": path.stat().st_size,
                    "md5": digest,
                    "name_match": name_match,
                    "size_match": size_match,
                    "md5_match": md5_match,
                    "accepted": ok,
                }
            )
            if ok and name not in accepted:
                accepted[name] = path
    pd.DataFrame(local_rows).to_csv(OUT / "local_md5.tsv", sep="\t", index=False)

    search_rows = [
        {
            "location": "/HDD8T1/WORK/GeroSci/clock-specificity",
            "what": "repository tree",
            "result": "no Zenodo clock pkl or module zip",
        },
        {
            "location": "/HDD8T1/WORK/GeroSci/sc-module-aging",
            "what": "data/external",
            "result": "EN multitissue scaleddiff pkl files and renamed module zips; see local_md5.tsv",
        },
        {
            "location": "/home/jkang/.cache",
            "what": "name search for 22166800, scaleddiff, module clock",
            "result": "no match",
        },
        {
            "location": "tAge download_clocks",
            "what": "default destination in the package",
            "result": "dest_dir is an argument; vignette uses tempdir(); no fixed cache path",
        },
    ]
    pd.DataFrame(search_rows).to_csv(OUT / "local_search.tsv", sep="\t", index=False)

    # Name-matched files are hardlinked. Zips whose content md5 matches but whose
    # filename differs are downloaded under the Zenodo name.
    for name, path in accepted.items():
        target = DEST / name
        if target.exists():
            continue
        os.link(path, target)
    for name in NEEDED:
        target = DEST / name
        if target.exists() and md5(target) == by_name[name]["md5"]:
            continue
        url = f"https://zenodo.org/api/records/22166800/files/{name.replace(' ', '%20')}/content"
        subprocess.run(
            ["curl", "-fL", "--retry", "5", "--max-time", "120", "-o", str(target), url],
            check=True,
        )
        digest = md5(target)
        if digest != by_name[name]["md5"]:
            raise SystemExit(f"md5 mismatch after download: {name} {digest}")

    placed = []
    for name in NEEDED:
        path = DEST / name
        placed.append(
            {
                "name": name,
                "path": str(path),
                "size": path.stat().st_size,
                "md5": md5(path),
                "zenodo_md5": by_name[name]["md5"],
                "md5_match": md5(path) == by_name[name]["md5"],
                "source": "local_hardlink" if name in accepted and not name.endswith(".zip") else "downloaded_or_replaced",
            }
        )
    # Recompute source more carefully.
    for row in placed:
        if row["name"] in accepted and row["name"].endswith(".pkl"):
            row["source"] = "local_hardlink"
        elif row["name"].endswith(".zip"):
            row["source"] = "downloaded"
    pd.DataFrame(placed).to_csv(OUT / "placed_files.tsv", sep="\t", index=False)

    module_dir = DEST / "Rodent module clocks"
    if not module_dir.exists():
        with zipfile.ZipFile(DEST / "Rodent module clocks.zip") as archive:
            archive.extractall(DEST)
    members = []
    for path in sorted(module_dir.glob("*.pkl")):
        members.append(
            {
                "name": path.name,
                "size": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(members).to_csv(OUT / "rodent_module_members.tsv", sep="\t", index=False)
    print("placed", len(placed), "module members", len(members))


if __name__ == "__main__":
    main()
