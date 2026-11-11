"""
Do the survey crest windows sit on the dune Hannah traced from the aerial imagery? The check of the window rule.

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/6-toes-vs-dunelines/HAT_crest_windows_vs_dunelines.py

The survey crests are read inside a window per 10 m cell (HAT_survey_dune_crests.py). Since 2026-11-03 the window
is the start's moved landward by the observed shoreline retreat since the start, never seaward (shoreline_retreat).
Hannah's traced dune lines are an independent record of where the dune was: for each survey with a line from the
same years (PAIRS in HAT_dune_toes_vs_dunelines.py: 1997, 2004 and 2009 on the 1996 start, the 2023 line with the
2020 lidar on the 2009 start), each cell's window is mapped to eastings through the start's straightened frame and
the line's crossing of the cell's profile is tested against it: inside the window (within HALF_CELL_M of its
edges), seaward of it or landward of it. The same for the start's fixed window and, where it exists, the
follow-the-dune trial (survey-crests/trial-follow/). Geometry only: no lidar is read.

Writes, in the workspace's dune-collision/1-observations/6-toes-vs-dunelines/supporting/:
crest_windows_vs_dunelines.csv (per survey, cell and window) and crest_windows_vs_dunelines_summary.csv (per survey,
window and zone: the share of cells whose line is inside, seaward, landward).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_toes_vs_dunelines as V  # noqa: E402

S, T = V.S, V.T
HALF_CELL_M = 5.0
# the hotspot domains, reported beside the zones
HOTSPOTS = (("hot GIS 9-13", 9, 13), ("hot GIS 22-25", 22, 25), ("hot GIS 51-55", 51, 55),
            ("hot GIS 75-77", 75, 77), ("hot GIS 84-88", 84, 88))


def shifts(start, folder=""):
    p = S.tv.product_dir(start) / "survey-crests" / folder / S.SHIFTS_FILE
    if not p.exists():
        return None
    t = pd.read_csv(p)
    return {(int(s), int(g)): d.sort_values("cell").shift_cells.to_numpy(int) for (s, g), d in t.groupby(["survey", "gis"])}


def main():
    rows = []
    for start, year, vintage, _ in V.PAIRS:
        windows = json.loads(S.tv.picks_dir(start).joinpath(
            f"HAT_dune_search_windows_{S.tv.topo_dirs(start)[2]}.json").read_text(encoding="utf-8"))
        prod, _ = S.STARTS[start]
        clips = S.product(prod).gapfill_1m
        grids = {g: S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in T.GIS}
        line = V.duneline(vintage, grids[1][2].to_wkt())
        modes = {"fixed": None, "retreat": shifts(start), "follow": shifts(start, "trial-follow")}
        for g in T.GIS:
            shear, c0, i0, i1 = S.frame(start, g, windows)
            xx, yy = V.frame_coords(grids[g], shear, c0)
            for cell in range(T.CELLS):
                y = np.nanmean(yy[cell])
                xs = xx[cell]
                if not np.isfinite(y) or not np.isfinite(xs).any():
                    continue
                xl, n = V.crossing(line, y, np.nanmax(xs), np.nanmin(xs), np.nanmean(xs[i0:i1]))
                if not np.isfinite(xl):
                    continue
                for mode, sh in modes.items():
                    if mode != "fixed" and sh is None:
                        continue
                    s = 0 if mode == "fixed" else int(sh.get((year, g), np.zeros(T.CELLS, int))[cell])
                    a, b = max(i0 + s, 0), min(i1 + s, xs.size)
                    if b <= a or not np.isfinite(xs[a:b]).any():
                        continue
                    # eastings fall landward along the profile: the seaward edge is the larger
                    sea, land = np.nanmax(xs[a:b]), np.nanmin(xs[a:b])
                    where = ("inside" if land - HALF_CELL_M <= xl <= sea + HALF_CELL_M
                             else "seaward" if xl > sea else "landward")
                    rows.append(dict(start=start, survey=year, line=vintage, gis=g, cell=cell, window=mode,
                                     shift_m=s * 10, line_x=xl, window_sea_x=sea, window_land_x=land, position=where))
        print(f"{start} {year} done", flush=True)
    t = pd.DataFrame(rows)
    out = V.OUT_DIR / "supporting"
    out.mkdir(parents=True, exist_ok=True)
    t.round(1).to_csv(out / "crest_windows_vs_dunelines.csv", index=False)
    summ = []
    for (start, year, mode), d in t.groupby(["start", "survey", "window"]):
        for name, lo, hi in (("island", 1, 90), *V.ZONES, *HOTSPOTS):
            x = d[d.gis.between(lo, hi)]
            if len(x):
                summ.append(dict(start=start, survey=year, window=mode, zone=name, cells=len(x),
                                 inside=(x.position == "inside").mean(), seaward=(x.position == "seaward").mean(),
                                 landward=(x.position == "landward").mean()))
    s = pd.DataFrame(summ)
    s.round(3).to_csv(out / "crest_windows_vs_dunelines_summary.csv", index=False)
    pd.set_option("display.width", 200)
    print(s.pivot_table(index=["survey", "zone"], columns="window", values="inside").round(2).to_string())


if __name__ == "__main__":
    main()
