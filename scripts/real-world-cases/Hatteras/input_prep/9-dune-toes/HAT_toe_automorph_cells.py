"""
Automorph on every 10 m model cell of a few domains (test, 2026-11-09): is a toe per cell any different from one per
domain? Not hand-checked.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_automorph_cells.py [--all]

--all: every domain, written to automorph_cell_toes_all_1996.csv.

Steps, each cell of GIS_RANGE:
1. The profile: the 1996 lidar every 1 m along the cell's centre row, east to west (the model's own cross-shore), from
   OCEAN_PAD_M offshore of the domain's clip.
2. The shoreline as the transects (HAT_toe_method_line.pick).
3. The crest: Automorph's automated crest with the 3 m floor (HAT_toe_method_automorph.auto_crest).
4. The toe: the adapted version on that crest, with its toe screen.
5. A cell with no toe takes the median of its domain's kept cells; a domain with fewer than MIN_KEPT kept cells keeps
   its hand-checked domain toe.

Writes, in data/hatteras_init/9-dune-toes/1-transects-1996/1-toes/ (D:): automorph_cell_toes_1996.csv, with toe_m_navd88
(this cell's toe, filled) and toe_domain_median_m_navd88 (the median of the domain's kept cells). Read by the
cal_automorph_cells_* variants of HAT_dune_collision.py.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-09
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import LineString

sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_toe_method_automorph as A  # noqa: E402
from HAT_toe_common import MIN_KEPT, PICKS, S, STEP_M, TOES, Windows, line_points, sample, survey  # noqa: E402
from HAT_toe_method_line import pick  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
ALL = "--all" in sys.argv     # every domain, in its own file, instead of the GIS 9-15 test
GIS_RANGE = range(1, 91) if ALL else range(9, 16)
OUT_NAME = f"automorph_cell_toes{'_all' if ALL else ''}_1996.csv"
CELLS = 50
OCEAN_PAD_M = 200.0        # the profile starts this far offshore of the clip, so that it crosses the shoreline
# -----------------------------------------------------------------------------


def cell_row(s, z):
    p = pick(s, z, mode="front")
    shore = p["shore"]
    out = dict(crest=None, toe=None, reason="no shoreline")
    if shore is None:
        return out
    c = A.auto_crest(z, shore)
    if c is None:
        return out | dict(reason="no automatic crest")
    a = A.automorph(s, z, shore, pick=dict(crest_x=float(s[c] - s[shore])), **A.VERSIONS["adapted_"][0])
    zs, g = A.smooth(z), lambda i: None if i is None else float(zs[i])
    return dict(crest=g(a["crest"]), toe=g(a["toe"]), reason=a["reason"] or "kept",
                crest_x=None if a["crest"] is None else float(s[a["crest"]] - s[shore]))


def main():
    src, scale, offset = survey()
    win = Windows()
    rows = []
    for g in GIS_RANGE:
        t, h, w, *_ = win.domain(g)
        _, _, crs = S.clip_grid(win.clips / f"clip_domain_{g}_filled.tif")
        east, west, bottom = t.c + w * t.a + OCEAN_PAD_M, t.c, t.f - h
        for c in range(CELLS):
            y = bottom + (c + 0.5) * S.BLOCK
            s, xy, *_ = line_points(LineString([(east, y), (west, y)]))
            rows.append(dict(gis=g, cell=c, northing=y, **cell_row(s, sample(src, scale, offset, xy, crs))))
    t = pd.DataFrame(rows).rename(columns=dict(crest="crest_m_navd88", toe="raw_toe_m_navd88"))
    t["kept"] = t.raw_toe_m_navd88.notna()
    dom = pd.read_csv(TOES / "automorph_domain_toes_1996.csv").set_index("gis").toe_m_navd88
    med = t[t.kept].groupby("gis").raw_toe_m_navd88.agg(["median", "count"])
    t["toe_domain_median_m_navd88"] = [med["median"][g] if g in med.index and med["count"][g] >= MIN_KEPT else dom[g]
                                       for g in t.gis]
    t["toe_m_navd88"] = t.raw_toe_m_navd88.fillna(t.toe_domain_median_m_navd88)
    t["hand_domain_toe_m_navd88"] = t.gis.map(dom)
    t.round(3).to_csv(TOES / OUT_NAME, index=False)

    print(t.groupby("gis").agg(kept=("kept", "sum"), low=("raw_toe_m_navd88", "min"),
                               cell_median=("toe_domain_median_m_navd88", "first"),
                               high=("raw_toe_m_navd88", "max"), sd=("raw_toe_m_navd88", "std"),
                               hand=("hand_domain_toe_m_navd88", "first")).round(2).to_string())
    print(t.reason.value_counts().to_string())

    # Against your checked transects: the cell nearest each one
    tr = pd.read_csv(TOES / "automorph_transect_toes_1996.csv")
    tr = tr[tr.kept & tr.gis.isin(GIS_RANGE)]
    from HAT_toe_common import transects
    geo = transects()[0].set_index("ObjectID").geometry
    rows = []
    for r in tr.itertuples():
        y = float(np.mean(np.asarray(geo[r.transect].coords)[:, 1]))
        c = t[t.gis == r.gis].iloc[(t[t.gis == r.gis].northing - y).abs().argmin()]
        rows.append(dict(transect=r.transect, gis=r.gis, cell=int(c.cell), hand_toe=r.toe_m_navd88,
                         cell_toe=c.raw_toe_m_navd88, hand_crest_x=r.crest_x, cell_crest_x=c.crest_x))
    v = pd.DataFrame(rows)
    d = (v.cell_toe - v.hand_toe).dropna()
    dx = (v.cell_crest_x - v.hand_crest_x).dropna()
    print(f"\nnearest cell against your checked transect toe, {len(d)} of {len(v)}: median {d.median():+.2f} m, RMSE "
          f"{np.sqrt((d ** 2).mean()):.2f} m, within 0.5 m {100 * (d.abs() <= 0.5).mean():.0f}%; crest position within "
          f"10 m {100 * (dx.abs() <= 10).mean():.0f}%")
    print(v.round(2).to_string(index=False))
    print(f"wrote {TOES / OUT_NAME}")


if __name__ == "__main__":
    main()
