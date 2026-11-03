"""
A year's storms on the Buxton dunes, storm by storm, through the dune collision module.

    python collision_year_gif.py [--year 1999] [--gis 8 15]
        ->  figures/collision_4_year_of_storms.gif   (and a still of the last frame, collision_4_year_of_storms.png)

Every 10 m dune cell of the chosen domains starts the year at its lidar crest (the survey flown before the
year, 1998-09 for 1999), with the toes and beach widths the hindcast reads (dune_toe_cells.csv). Each storm of
the year in the hindcast storm series is then put through the module's own functions (collision_cells,
larson_volume, lower_dunes) at the hindcast's settings, as DuneCollisionManager.update does: the cells it
leaves in collision are lowered, and the next storm sees the lowered dunes. Two frames a storm: the water
arriving, then the dunes after. Cells the storm overtops are Barrier3D's own overwash, marked but not
lowered here; dune growth (once a year, before the storms) is left out. Nothing is fitted or run; seconds.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import collision_explained_figures as E  # noqa: E402  (paths, settings, units, the storm series)
from cascade.dune_collision import beach_factor, collision_cells, larson_volume, lower_dunes  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _title, apply_style, caption, figsize, open_frame, save)

# --- CONFIG ------------------------------------------------------------------
INIT = E.WORKSPACE / "data" / "hatteras_init" / "1-barrier3d-domains" / "1996-start"
CRESTS = INIT / "survey-crests" / "survey_dune_crests.csv"
TOES = INIT / "dune-toes" / "dune_toe_cells.csv"
YEAR = 1999
GIS = (8, 15)                   # first and last domain drawn (Buxton and the road north of it)
CELLS = 50                      # 10 m dune cells per domain
FPS = 1.2
HOLD_END = 4                    # extra frames on the last one
DPI = 110
C_START = "0.55"                # the crest at the start of the year
# -----------------------------------------------------------------------------


# The cells: start crest (the survey before the year), toe and beach width, alongshore position
def cells(year, gis):
    c = pd.read_csv(CRESTS)
    c = c[(c.survey == year - 1) & c.gis.between(*gis)]
    t = pd.read_csv(TOES)
    d = t[t.gis.between(*gis)].merge(c[["gis", "cell", "crest_m_navd88", "status", "flown"]],
                                    on=["gis", "cell"], how="left")
    d.loc[d.status != "ok", "crest_m_navd88"] = np.nan
    d = d.sort_values(["gis", "cell"]).reset_index(drop=True)
    return d.assign(x=d.gis - 0.5 + (d.cell + 0.5) / CELLS)


# Every storm of the year on every cell, as DuneCollisionManager.update: the crest before and after each
def run_year(d, st):
    crest0 = d.crest_m_navd88.to_numpy()
    ok = np.isfinite(crest0)
    h = np.where(ok, crest0 - E.BERM, 0.0) / 10.0
    rows = np.column_stack([h, h])                     # both rows at the crest: lowered alike
    pre = np.where(ok, E.dam(crest0), np.inf)          # the pre-storm crest, kept after dune growth
    toe = E.dam(d.toe_m_navd88.to_numpy())
    beach = beach_factor(d.beach_width_m.to_numpy() / 10.0, E.W0 / 10.0)
    steps = []
    for _, s in st.iterrows():
        r = E.dam(s.r_m)
        crest = rows.max(axis=1) + E.dam(E.BERM)
        hit = collision_cells(pre, crest, toe, r, beach) & ok
        over = ok & (pre < r)
        before = np.where(ok, E.m(crest), np.nan)
        if hit.any():
            v = larson_volume(r, toe[hit], s.period, s.duration, E.CS, beach[hit])
            rows[hit] = lower_dunes(rows[hit], v, toe[hit], E.dam(E.BERM), E.DUNE_W / 10.0, E.RESTART / 10.0)
        after = np.where(ok, E.m(rows.max(axis=1) + E.dam(E.BERM)), np.nan)
        steps.append(dict(storm=s, hit=hit, over=over, before=before, after=after))
    return steps


# One frame: the dunes alongshore, the beach under them, the year's storms below
def draw(fig, axes, d, st, steps, k, phase, year):
    a, w, b = axes
    for ax in axes:
        ax.clear()
    x = d.x.to_numpy()
    crest0 = d.crest_m_navd88.to_numpy()
    s = steps[k]
    now = s["before"] if phase == "arrive" else s["after"]
    storm = s["storm"]

    # (a) crests: the start, now, what is lost; the toe; the water while the storm is on
    a.fill_between(x, now, crest0, where=crest0 > now + 1e-6, color=E.C_MOD_FILL, lw=0, step="mid",
                   label="lowered since the start of the year")
    a.plot(x, crest0, color=C_START, lw=0.8, drawstyle="steps-mid", label=f"crest at the start ({d.flown.dropna().iloc[0]} lidar)")
    a.plot(x, now, color=E.C_MOD, lw=1.4, drawstyle="steps-mid", label="crest now")
    a.plot(x, d.toe_m_navd88, color=INK_MUTED, lw=0.7, ls=":", drawstyle="steps-mid", label="toe")
    if phase == "arrive":
        a.axhspan(0, storm.r_m, color=E.C_WATER, alpha=0.45, lw=0)
        a.axhline(storm.r_m, color=C["LATE"], lw=1.2)
    # which cells this storm acts on: a tick rug along the top, so the crest line stays readable
    a.vlines(x[s["hit"]], 10.55, 10.95, color=C["LATE"], lw=0.9, label="in collision this storm (ticks)")
    a.vlines(x[s["over"]], 10.55, 10.95, color=C_START, lw=0.9, label="overtopped: Barrier3D's overwash (ticks)")
    for g in range(GIS[0], GIS[1]):
        a.axvline(g + 0.5, color="0.9", lw=0.6, zorder=0)
    a.set(xlim=(GIS[0] - 0.5, GIS[1] + 0.5), ylim=(0, 11), ylabel="Elevation (m NAVD88)")
    a.set_yticks(range(0, 10, 2))
    a.set_xticks(range(GIS[0], GIS[1] + 1))
    a.tick_params(labelbottom=False)
    a.legend(loc="upper left", bbox_to_anchor=(0, 0.95), ncol=3, frameon=False, fontsize=7, columnspacing=1.0, handlelength=1.4)
    open_frame(a)
    n_hit = int(s["hit"].sum())
    lost = np.nanmean(crest0 - now)
    when = storm.StartTime.strftime("%d %b %Y")
    verb = "arrives" if phase == "arrive" else "after"
    a.set_title(f"Storm {k + 1} of {len(st)} {verb}: {when}, R_high {storm.r_m:.2f} m for {storm.duration:.0f} h"
                f"  ·  {n_hit} cells in collision  ·  mean lowering since January {lost:.2f} m",
                fontsize=8.5, loc="center", pad=4)

    # (b) the beach in front of each dune: the module's scale on the loss
    bw = d.beach_width_m.to_numpy()
    w.bar(x, bw, width=1 / CELLS, color=np.where(bw >= E.W0, "0.8", E.C_MOD_FILL), lw=0)
    w.axhline(E.W0, color=INK_MUTED, lw=0.7, ls=":")
    w.text(GIS[1] + 0.48, E.W0 + 3, f"{E.W0:.0f} m: no collision loss", ha="right", va="bottom", fontsize=7,
           color=INK_MUTED)
    w.set(xlim=a.get_xlim(), ylim=(0, 100), ylabel="Beach (m)", xlabel="GIS domain (south → north)")
    w.set_xticks(range(GIS[0], GIS[1] + 1))
    open_frame(w)

    # (c) the year's storms, this one lit
    t = st.StartTime
    col = np.where(np.arange(len(st)) == k, C["LATE"], E.C_WATER)
    b.vlines(t, 0, st.r_m, colors=col, lw=np.where(np.arange(len(st)) == k, 2.4, 1.4))
    b.axhline(np.nanmedian(d.toe_m_navd88), color=INK_MUTED, lw=0.7, ls=":")
    b.text(pd.Timestamp(f"{year}-12-31"), np.nanmedian(d.toe_m_navd88) + 0.1, "median toe", ha="right",
           va="bottom", fontsize=7, color=INK_MUTED)
    b.set(xlim=(pd.Timestamp(f"{year}-01-01"), pd.Timestamp(f"{year}-12-31")), ylim=(0, 4.2),
          ylabel="R_high (m)")
    b.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
    open_frame(b)


# Run: cells, the year, the frames, the still
def main() -> None:
    global GIS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--year", type=int, default=YEAR)
    ap.add_argument("--gis", type=int, nargs=2, default=GIS)
    args = ap.parse_args()
    GIS = tuple(args.gis)

    apply_style()
    d = cells(args.year, GIS)
    st = E.storms()
    st = st[st.calendar_year == args.year].sort_values("StartTime").reset_index(drop=True)
    steps = run_year(d, st)

    fig, axes = plt.subplots(3, 1, figsize=figsize("double", 0.66), layout="constrained",
                             gridspec_kw=dict(height_ratios=(3.2, 0.9, 1.1)))
    out = E.OUT / "collision_4_year_of_storms.gif"
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=DPI):
        for k in range(len(steps)):
            for phase in ("arrive", "after"):
                draw(fig, axes, d, st, steps, k, phase, args.year)
                writer.grab_frame()
        for _ in range(HOLD_END):
            writer.grab_frame()

    crest0 = d.crest_m_navd88.to_numpy()
    end = steps[-1]["after"]
    lowered = np.nansum((crest0 - end) > 0.01)
    acted = sum(int(s["hit"].any()) for s in steps)
    big = max(range(len(steps)), key=lambda i: int(steps[i]["hit"].sum()))
    caption(fig, f"**Tests:** nothing new; this figure (the last frame of `collision_4_year_of_storms.gif`) "
            f"shows the dune collision module acting on a real stretch of coast through a year of storms. "
            f"**How to read:** (a) every 10 m dune cell of GIS {GIS[0]}-{GIS[1]} ({len(d)} cells, "
            f"{int(np.isfinite(crest0).sum())} with a lidar crest), starting the year at its {args.year - 1} "
            f"lidar crest (grey), with the toes and beach widths the hindcast reads (dotted, the toe; "
            f"(b) the dry beach, grey where {E.W0:.0f} m or wider). Each of the {len(st)} storms of {args.year} in the "
            f"hindcast storm series is put through the module's own functions at the hindcast's settings "
            f"(Cs {E.CS:g}, W0 {E.W0:.0f} m): blue ticks along the top, the cells a storm leaves in collision; purple line "
            f"and fill, the crest now and what has been lost since January; grey ticks, cells a storm overtops, "
            f"which Barrier3D's own overwash lowers (not drawn here). Dune growth, once a year before the "
            f"storms, is left out. (c) the year's storms, the current one dark. In the GIF each storm has two "
            f"frames: the water arriving (blue), then the dunes after. **Shows:** {acted} of the {len(st)} "
            f"storms lower some dune; storm {big + 1} ({st.StartTime[big]:%d %b}, R_high "
            f"{st.r_m[big]:.2f} m, {st.duration[big]:.0f} h) puts the most cells in collision "
            f"({int(steps[big]['hit'].sum())}). By December {lowered} cells are lower, by "
            f"{np.nanmean(crest0 - end):.2f} m on average and up to {np.nanmax(crest0 - end):.2f} m, "
            f"most where the beach is narrow; behind dry beaches of {E.W0:.0f} m or more nothing is lost.")
    save(fig, E.OUT / "collision_4_year_of_storms", close=True)
    print(f"{out}\n{acted} storms act; {lowered} cells lowered, mean {np.nanmean(crest0 - end):.2f} m, "
          f"max {np.nanmax(crest0 - end):.2f} m")


if __name__ == "__main__":
    main()
