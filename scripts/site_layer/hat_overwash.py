# ==============================================================================
# hat_overwash.py
#
# WHERE DOES THE OBSERVED OVERWASH RECORD LIVE, AND ITS FIGURES AND TABLES?
#
# WHY THIS EXISTS
#   overwash_data.py already held the root for its three sibling scripts, but
#   the two model-comparison scripts in scripts/analyze_output/overwash/ typed
#   the workbook's path themselves, and both named a location it had left
#   (scripts/input_prep/8-overwash-analysis/), so neither could run. Same
#   answer as hat_observed_rates.py for 5-scr: resolved ONCE (2026-09-18).
#
# THE LAYOUT, grouped by job (2026-09-18). Before that it was grouped by file
# type -- figures/ and tables/ -- which split the footprint comparison across
# figures/vs-footprint/ and a top-level vs-footprint/.
#
#     data/hatteras_init/8-overwash-analysis/
#         README.md
#         CAPTIONS.md               one entry per figure, headed by its folder
#         1-observations/           THE RECORD, and its long-form tables
#             Hatteras_Overwash_Data.xlsx
#             overwash_observations.csv  storms_by_image.csv
#         2-record/                 what the record shows
#             heatmaps/  map/
#         3-vs-footprint/           the record against the 1984 footprint:
#                                   figures here, their tables in tables/
#         archive/                  superseded_20260910/
# ==============================================================================

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

ARCHIVE = OVERWASH_ROOT / "archive"

# The label each figure's CAPTIONS.md entry carries, keyed by the short
# folder name the scripts pass to overwash_data.upsert_caption. Its order is
# the order the entries are sorted into.
CAPTION_FOLDERS = {
    "heatmaps": "2-record/heatmaps",
    "map": "2-record/map",
    "vs-footprint": "3-vs-footprint",
}


if __name__ == "__main__":
    for name in ("OVERWASH_ROOT", "CAPTIONS", "OBSERVATIONS", "WORKBOOK",
                 "RECORD", "HEATMAPS", "MAP", "VS_FOOTPRINT",
                 "VS_FOOTPRINT_TABLES", "ARCHIVE"):
        path = globals()[name]
        print(f"{'ok' if path.exists() else 'MISSING':8} {name:20} "
              f"{path.relative_to(PROJECT_ROOT).as_posix()}")
