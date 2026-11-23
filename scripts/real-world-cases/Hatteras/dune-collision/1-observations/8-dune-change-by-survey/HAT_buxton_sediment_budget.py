"""
Was the Buxton erosion overwash-driven? A subaerial sediment budget from the lidar (2026-11-23).

    python scripts/real-world-cases/Hatteras/dune-collision/1-observations/8-dune-change-by-survey/HAT_buxton_sediment_budget.py

If the barrier rolled over by overwash, the sand lost from the beach and dune reappears behind the dune (on the road,
the interior, or pushed back by NCDOT's bulldozers) and the profile's volume above mean high water stays roughly
constant. If the sand went offshore or alongshore, the total falls. On each 1 m alongshore line of the coast unrolled
along NC-12's 2008 alignment, every survey on the 2009 datum (hat_lidar_strips.surveys), the volume above MHW
(0.36 m NAVD88) per m of shore:
- landward start: the most seaward of the surveys' first valid points, so every survey covers the same ground
  (no further landward than LAND_M behind the road)
- split at 10 m landward of the line's fall-1997 dune crest (the highest point at least 10 m seaward of the road):
  "seaward" (dune front, beach) and "landward" (old dune back, road, interior)
- seaward of each survey's last valid point the ground is taken as below MHW (water); gaps inside are interpolated
Per domain, the median over its lines; GIS 9-15 (the NCDOT hot spot and its edges) against GIS 16-22 (control).

Writes, in the workspace's dune-collision/1-observations/8-dune-change-by-survey/supporting/ (D:):
- buxton_sediment_budget.csv: per domain and survey, the landward, seaward and total volume change from 1997 (m3/m),
  and the share of lines with data

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-23
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_change_by_survey as D  # noqa: E402

L = D.L
OUT = D.OUT
# --- CONFIG ------------------------------------------------------------------
GIS = range(9, 23)
LAND_M = 250.0                 # at most this far landward of NC-12 (2008)
SEA_M = 320.0
REF = 1997                     # the pre-storm survey
SPLIT_BEHIND_CREST_M = 10.0
STEP = 1.0
# -----------------------------------------------------------------------------


def line_volumes(prof_by_year, d):
    """One line: {year: (landward, seaward)} m3/m above MHW, on the ground every covering survey has; a survey with no
    data seaward of the road on this line is left out (NaN)."""
    have = {y: p for y, p in prof_by_year.items() if np.isfinite(p[d >= 0]).sum() > 20}
    if REF not in have:
        return None
    i0 = max(np.where(np.isfinite(p))[0][0] for p in have.values())
    ref = have[REF]
    sea = (d >= 10) & (np.arange(d.size) >= i0)
    if not np.isfinite(ref[sea]).any():
        return None
    k = np.where(sea)[0][np.nanargmax(ref[sea])]
    split = d[k] - SPLIT_BEHIND_CREST_M
    out = {y: (np.nan, np.nan) for y in prof_by_year}
    for y, p in have.items():
        q = p[i0:].copy()
        dd = d[i0:]
        ok = np.where(np.isfinite(q))[0]
        if ok.size < 20:
            continue
        q[ok[-1] + 1:] = L.MHW_M - 1.0                     # water seaward of the last return
        bad = ~np.isfinite(q)
        if bad.any():
            q[bad] = np.interp(np.where(bad)[0], np.where(~bad)[0], q[~bad])
        h = np.clip(q - L.MHW_M, 0, None) * STEP
        out[y] = (float(h[dd < split].sum()), float(h[dd >= split].sum()))
    return out


def measure():
    r = L.Region("budget", "GIS 9-22", GIS, (), sea_m=SEA_M, land_m=LAND_M)
    surv, g, dom = L.surveys(r)
    st = L.Strip(r, g, L.road_lines(g[2]), dom)
    d = np.arange(-LAND_M, SEA_M, STEP)
    P = {}
    for y, (z, _, meta) in surv.items():
        zs = st.sample(z)
        P[y] = np.array([np.interp(d, st.n - st.road08[j], zs[:, j], left=np.nan, right=np.nan)
                         if np.isfinite(st.road08[j]) else np.full(d.size, np.nan) for j in range(len(st.s))])
    years = sorted(P)
    rows = []
    for gis, (a, b) in sorted(st.spans.items()):
        lines = np.where((st.s >= a) & (st.s <= b))[0]
        res = [line_volumes({y: P[y][j] for y in years}, d) for j in lines]
        res = [x for x in res if x is not None]
        for y in years:
            land = np.array([x[y][0] - x[REF][0] for x in res])
            sea = np.array([x[y][1] - x[REF][1] for x in res])
            ok = np.isfinite(land)
            if ok.sum() < 50:
                continue
            rows.append(dict(gis=gis, survey=y, lines=int(ok.sum()), share_lines=ok.sum() / max(len(lines), 1),
                             landward_m3_per_m=float(np.median(land[ok])), seaward_m3_per_m=float(np.median(sea[ok])),
                             total_m3_per_m=float(np.median(land[ok] + sea[ok]))))
    t = pd.DataFrame(rows)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(2).to_csv(OUT / "supporting" / "buxton_sediment_budget.csv", index=False)
    print(f"wrote {OUT / 'supporting' / 'buxton_sediment_budget.csv'}")
    return t


def report(t):
    pd.set_option("display.width", 250)
    for col in ("seaward_m3_per_m", "landward_m3_per_m", "total_m3_per_m"):
        print(f"\n{col} (change from {REF}, median over lines):")
        print(t.pivot(index="gis", columns="survey", values=col).round(0).to_string())
    print("\nlines per domain:", t.groupby("gis").lines.first().to_dict())


if __name__ == "__main__":
    report(pd.read_csv(OUT / "supporting" / "buxton_sediment_budget.csv") if "report" in sys.argv[1:] else measure())
