"""The 2026-09-22 shoreline-offset study, drawn in the island-offset house form (2026-09-28).

Hannah, 2026-09-28: make the figures across the island-offset experiments
consistent with the most recent version, so they are easier to compare. The
09-22 study had only the runner's per-run figures. This draws its runs with
the shared house_figures (HAT_offset_source_comparison.py) and runs nothing.

    shoreline  the study's own run,
               island-offset/2026-09-22-div10-offset-shoreline-trial-original/
               1996_2010/edgeBE/HAT_1996_2010_edgeBE_road_bdm_nogroin
    dune line  its matrix control, since archived:
               archive/2026-09-24-pre-metres/matrix/1996_2010/edgeBE/
               HAT_1996_2010_edgeBE_road_bdm_nogroin
    both       ÷10 offset (asrun) on the pre-metres builds, Hs 2.5 / Tp 8 /
               asymmetry 0.7 / high-angle 0.1, edgeBE +32.2 / +10.0 m/yr,
               full management, no groin, the unfixed Barrier3D. The control
               ran 2026-09-18 (7e492c3), the shoreline arm 2026-09-22 (c8ed4a1)

The road_reloc_bdm arm is left out: no relocation falls inside 1996-2010, so
it is the same run (NOTE.md), and it has no archived control.

    python scripts/hatteras_ms/experiments/HAT_offset_source_0922_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_offset_source_comparison as base  # noqa: E402

RAW = base.grid.RAW_RUNS
STUDY_DIR = RAW / "experiments" / "island-offset" / "2026-09-22-div10-offset-shoreline-trial-original"
RUN = "HAT_1996_2010_edgeBE_road_bdm_nogroin"
RUNS = {"shoreline": STUDY_DIR / "1996_2010" / "edgeBE" / RUN,
        "duneline": RAW / "archive" / "2026-09-24-pre-metres" / "matrix" / "1996_2010"
                    / "edgeBE" / RUN}


def main():
    import pandas as pd
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    # (b) 2010-2024 was never run in this study: the panel carries the
    # observation alone so the layout matches the others (full management x
    # both periods, Hannah 2026-09-28); no run is added to a record
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
