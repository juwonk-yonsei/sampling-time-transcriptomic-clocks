"""Project paths. tAge resources (MGB non-commercial licence) are read from an external
directory and are never copied into this repository."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
CURATION = DATA / "curation"
RESULTS = ROOT / "results"

TAGE = Path(os.environ.get("CLOCKSPEC_TAGE", "/HDD8T1/WORK/GeroSci/sc-module-aging/data/external"))
TAGE_META = TAGE / "tAge" / "inst" / "extdata"
TAGE_TRAINING = TAGE / "tage_zenodo" / "Data_annotation_relative_rodents.xlsx"
GTEX = Path(os.environ.get("CLOCKSPEC_GTEX", "/HDD8T1/WORK/GeroSci/sc-module-aging/data/raw/GTEx"))


def ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path
