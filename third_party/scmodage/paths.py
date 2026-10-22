from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
EXTERNAL = DATA / "external"
PROCESSED = DATA / "processed"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

TAGE_REPO = EXTERNAL / "tAge"
TAGE_META = TAGE_REPO / "inst" / "extdata"
TAGE_ZENODO = EXTERNAL / "tage_zenodo"
TAGE_MODULES = EXTERNAL / "tage_modules" / "rodent" / "Rodent module clocks"
RETRAINED_CLOCKS = PROCESSED / "retrained_clocks"


def ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path
