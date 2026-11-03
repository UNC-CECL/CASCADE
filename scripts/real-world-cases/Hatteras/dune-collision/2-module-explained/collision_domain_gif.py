"""
A year's storms on one domain's dunes, storm by storm: the simple version of collision_year_gif.py.

    python collision_domain_gif.py [--year 1999] [--gis 12]
        ->  figures/collision_5_one_domain.gif   (and its last frame, collision_5_one_domain.png)

The domain's 50 dune cells as bars, each at its lidar crest from the survey before the year, with the toes
the hindcast reads; each storm of the year in the hindcast storm series lowers them through the module's own
functions, as collision_year_gif.py does (whose cells() and run_year() this uses). GIS 12 by default: a full
lidar record, a 0-20 m dry beach everywhere, no 1999 storm overtopping it, so all it shows is collision.

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
import collision_explained_figures as E  # noqa: E402  (paths, settings, the storm series)
import collision_year_gif as G  # noqa: E402  (the cells and the year's storms through the module)
from site_layer.hat_figure_style import C, INK_MUTED, apply_style, caption, figsize, open_frame, save  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
YEAR = 1999
GIS = 12
FPS = 1.0
DPI = 110
C_DUNE = E.C_SAND
C_LOST = E.C_MOD
# -----------------------------------------------------------------------------


# One frame: the dunes as bars, the water while the storm is on, the year's storms below
def draw(axes, d, st, steps, k, phase, year, gis):
    a, b = axes
    a.clear(); b.clear()
    x = d.cell.to_numpy() * 10 + 5                      # m alongshore within the domain
    crest0 = d.crest_m_navd88.to_numpy()
    s = steps[k]
    now = s["before"] if phase == "arrive" else s["after"]
    storm = s["storm"]

    a.bar(x, now - E.BERM, bottom=E.BERM, width=9, color=C_DUNE, lw=0)
    a.bar(x, crest0 - now, bottom=now, width=9, color=C_LOST, lw=0, label="lost since January")
    a.step(x, d.toe_m_navd88, where="mid", color="0.25", lw=1, ls=":", label="dune toe")
    if phase == "arrive":
        a.axhspan(0, storm.r_m, color=E.C_WATER, alpha=0.55, lw=0, zorder=3)
        a.text(497, storm.r_m + 0.12, f"storm water level {storm.r_m:.2f} m", ha="right", va="bottom",
               fontsize=8, color=C["LATE"], zorder=4)
    a.set(xlim=(0, 500), ylim=(0, 9), xlabel=f"Alongshore in GIS {gis} (m, south → north)",
          ylabel="Elevation (m NAVD88)")
    a.legend(loc="upper left", ncol=2, frameon=False, fontsize=8)
    open_frame(a)
    lost = np.nanmean(crest0 - now)
    head = (f"{storm.StartTime:%d %b %Y}: water at {storm.r_m:.2f} m for {storm.duration:.0f} h"
            if phase == "arrive" else f"After the storm: {int(s['hit'].sum())} of {len(d)} dunes lowered")
    a.set_title(f"Storm {k + 1} of {len(st)}  ·  {head}  ·  average lost so far {lost:.2f} m",
                fontsize=9, loc="center")

    on = np.arange(len(st)) == k
    b.vlines(st.StartTime, 0, st.r_m, colors=np.where(on, C["LATE"], E.C_WATER), lw=np.where(on, 2.6, 1.4))
    b.axhline(np.nanmedian(d.toe_m_navd88), color=INK_MUTED, lw=0.7, ls=":")
    b.text(pd.Timestamp(f"{year}-12-31"), np.nanmedian(d.toe_m_navd88) + 0.1, "typical toe", ha="right",
           va="bottom", fontsize=7, color=INK_MUTED)
    b.set(xlim=(pd.Timestamp(f"{year}-01-01"), pd.Timestamp(f"{year}-12-31")), ylim=(0, 4.2),
          ylabel="Storm water\nlevel (m)")
    b.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
    open_frame(b)


# Run: the domain, the year, the frames, the still
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--year", type=int, default=YEAR)
    ap.add_argument("--gis", type=int, default=GIS)
    args = ap.parse_args()

    apply_style()
    d = G.cells(args.year, (args.gis, args.gis))
    st = E.storms()
    st = st[st.calendar_year == args.year].sort_values("StartTime").reset_index(drop=True)
    steps = G.run_year(d, st)

    fig, axes = plt.subplots(2, 1, figsize=figsize("double", 0.55), layout="constrained",
                             gridspec_kw=dict(height_ratios=(3, 1)))
    out = E.OUT / "collision_5_one_domain.gif"
    out.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=DPI):
        for k in range(len(steps)):
            for phase in ("arrive", "after"):
                draw(axes, d, st, steps, k, phase, args.year, args.gis)
                writer.grab_frame()
        for _ in range(3):                 # identical frames merge into one held ~4 s
            writer.grab_frame()

    crest0 = d.crest_m_navd88.to_numpy()
    end = steps[-1]["after"]
    big = max(range(len(steps)), key=lambda i: np.nansum(steps[i]["before"] - steps[i]["after"]))
    big_lost = np.nanmean(steps[big]["before"] - steps[big]["after"])
    caption(fig, f"**Tests:** nothing new; the last frame of `collision_5_one_domain.gif`, the dune collision "
            f"module on one domain through a year of storms (the simple version of "
            f"`collision_4_year_of_storms.gif`). **How to read:** (a) the {len(d)} dune cells of GIS {args.gis} "
            f"(10 m each) as bars from the berm to the crest, starting at the {d.flown.dropna().iloc[0]} lidar "
            f"crest, with the hindcast's toes (dotted) and beach widths ({d.beach_width_m.min():.0f}-"
            f"{d.beach_width_m.max():.0f} m). Each of the {len(st)} storms of {args.year} in the hindcast "
            f"storm series lowers them through the module's own functions at the hindcast's settings "
            f"(Cs {E.CS:g}, W0 {E.W0:.0f} m); purple, what has been lost since January. In the GIF each storm "
            f"has two frames, the water arriving (blue) and the dunes after. (b) the year's storms, the current "
            f"one dark. Dune growth (once a year, before the storms) is left out. **Shows:** storm {big + 1} "
            f"({st.StartTime[big]:%d %b}, {st.r_m[big]:.2f} m for {st.duration[big]:.0f} h) does most of the "
            f"year's work, {big_lost:.2f} m on average; by December the dunes are {np.nanmean(crest0 - end):.2f} m "
            f"lower on average, the low dunes and those with low toes most.")
    save(fig, E.OUT / "collision_5_one_domain", close=True)
    print(f"{out}\nGIS {args.gis}: storm {big + 1} -{big_lost:.2f} m; year -{np.nanmean(crest0 - end):.2f} m mean")


if __name__ == "__main__":
    main()
