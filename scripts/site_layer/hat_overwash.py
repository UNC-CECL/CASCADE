"""
Where does the observed overwash record live, and its figures and tables?

    from site_layer.hat_overwash import OVERWASH_ROOT

Resolves data/hatteras_init/8-overwash-analysis/ once: the record, what it shows,
and the comparisons against the footprint and the model. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""

from __future__ import annotations

from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

OVERWASH_ROOT = INIT_ROOT / "8-overwash-analysis"
CAPTIONS = OVERWASH_ROOT / "CAPTIONS.md"

OBSERVATIONS = OVERWASH_ROOT / "1-observations"
# The hand-digitised record. Edit it here; everything else is drawn from it.
WORKBOOK = OBSERVATIONS / "Hatteras_Overwash_Data.xlsx"

RECORD = OVERWASH_ROOT / "2-record"
HEATMAPS = RECORD / "heatmaps"
MAP = RECORD / "map"

VS_FOOTPRINT = OVERWASH_ROOT / "3-vs-footprint"
VS_FOOTPRINT_TABLES = VS_FOOTPRINT / "tables"

VS_MODEL = OVERWASH_ROOT / "4-vs-model"
VS_MODEL_TABLES = VS_MODEL / "tables"

ARCHIVE = OVERWASH_ROOT / "archive"

# The label each figure's CAPTIONS.md entry carries, in the order entries are sorted
CAPTION_FOLDERS = {
    "heatmaps": "2-record/heatmaps",
    "map": "2-record/map",
    "vs-footprint": "3-vs-footprint",
    "vs-model": "4-vs-model",
}


if __name__ == "__main__":
    for name in ("OVERWASH_ROOT", "CAPTIONS", "OBSERVATIONS", "WORKBOOK",
                 "RECORD", "HEATMAPS", "MAP", "VS_FOOTPRINT",
                 "VS_FOOTPRINT_TABLES", "VS_MODEL", "VS_MODEL_TABLES", "ARCHIVE"):
        path = globals()[name]
        print(f"{'ok' if path.exists() else 'MISSING':8} {name:20} "
              f"{path.relative_to(PROJECT_ROOT).as_posix()}")
