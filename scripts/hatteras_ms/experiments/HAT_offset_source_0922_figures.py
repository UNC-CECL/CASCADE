"""
The 2026-09-22 shoreline-offset study, drawn in the island-offset house form.

    python scripts/hatteras_ms/experiments/HAT_offset_source_0922_figures.py

Draws the study's run and its archived matrix control with the shared
figures; runs nothing. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-28
"""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_offset_source_comparison as base  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = base.grid.RAW_RUNS
STUDY_DIR = RAW / "experiments" / "island-offset" / "2026-09-22-div10-offset-shoreline-trial-original"
RUN = "HAT_1996_2010_edgeBE_road_bdm_nogroin"
RUNS = {"shoreline": STUDY_DIR / "1996_2010" / "edgeBE" / RUN,
        "duneline": RAW / "archive" / "2026-09-24-pre-metres" / "matrix" / "1996_2010"
                    / "edgeBE" / RUN}
# -----------------------------------------------------------------------------


# Run: the figures
def main():
    import pandas as pd
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    # 2010-2024 was not run in this study: panel (b) shows the observation alone
    panels = [dict(label="Full management 1996–2010", start=1996, end=2010,
                   rates={src: pd.read_csv(d / "tables" / "shoreline_change_rate.csv")
                          .set_index("gis_domain") for src, d in RUNS.items()}),
              dict(label="Full management 2010–2024: not run in this study",
                   start=2010, end=2024, rates={})]
    note = ("(a) full management 1996-2010, as run on 2026-09-22; (b) 2010-2024 was not "
            "run in this study and shows the observation only. (a): island offset ÷10 (offset_mode "
            "asrun) on superseded_20260924_pre-metres/v1 of each source; waves Hs 2.5 m, Tp 8 s, "
            "asymmetry 0.7, high-angle 0.1 (the ÷10-era defaults); edge source/sink "
            "correction at the ÷10-era ends, GIS 1 +32.2 and GIS 90 +10.0 m/yr (edgeBE); no "
            "groin; Barrier3D before the route_overwash fix. The dune-line arm is the matrix "
            "control, run 2026-09-18 (7e492c3) and since archived under "
            "archive/2026-09-24-pre-metres/; the shoreline arm ran 2026-09-22 (c8ed4a1). Drawn "
            "2026-09-28 in the island-offset house form.")
    legend = ("÷10 island offset, as run 2026-09-22; edge source/sink correction "
              "(÷10-era ends +32.2 / +10.0 m/yr)")
    for p in base.house_figures(panels, STUDY_DIR / "figures", legend, note, "full_management"):
        print(p.relative_to(STUDY_DIR))
    return 0


if __name__ == "__main__":
    sys.exit(main())
