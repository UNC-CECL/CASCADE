"""
One storm along the whole island: which dunes Dennis leaves in collision, overwashes or leaves alone.

    python collision_dennis_island.py [--year 1999]
        ->  figures/collision_8_dennis_island.png

Every 10 m dune cell of GIS 1-90 at its lidar crest from the survey before the year (1998-09 for 1999), with
the toes and beach widths the hindcast reads, against the year's largest storm in the hindcast storm series
(Dennis, one water level along the whole island, as Barrier3D takes it). Each cell's regime: overtopped
(Barrier3D's overwash), in collision and lowered by the module, in collision but behind a dry beach of 25 m or
more (the module takes nothing), or below the toe. Then the module's crest lowering per domain, through its
own functions (collision_year_gif.py's cells() and run_year()). Barrier3D's overwash lowering is not computed.

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
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import collision_explained_figures as E  # noqa: E402  (settings, the storm series)
import collision_year_gif as G  # noqa: E402  (the cells and a storm through the module)
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, caption, figsize, open_frame, save, town_bands)

# --- CONFIG ------------------------------------------------------------------
YEAR = 1999
GIS = (1, 90)
REGIMES = (("overwash", "overtopped: Barrier3D's overwash", E.C_B3D),
           ("collision", "collision: lowered by the module", E.C_MOD),
           ("wide", f"collision, dry beach ≥ {E.W0:.0f} m: nothing taken", E.C_MOD_FILL),
           ("below", "water below the toe", "0.85"),
           ("nodata", "no lidar crest", "white"))
# -----------------------------------------------------------------------------


# Each cell's regime in the storm
def regimes(d, r):
    crest, toe, w = d.crest_m_navd88.to_numpy(), d.toe_m_navd88.to_numpy(), d.beach_width_m.to_numpy()
    out = np.full(len(d), "below", dtype=object)
    out[np.isnan(crest)] = "nodata"
    ok = ~np.isnan(crest)
    out[ok & (crest <= r)] = "overwash"
    reach = ok & (crest > r) & (toe < r)
    out[reach & (w < E.W0)] = "collision"
    out[reach & (w >= E.W0)] = "wide"
    return out


# Run: cells, the storm, regimes and lowering by domain, the figure
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--year", type=int, default=YEAR)
    args = ap.parse_args()

    d = G.cells(args.year, GIS)
    st = E.storms()
    yr = st[st.calendar_year == args.year]
    storm = yr.loc[[yr.r_m.idxmax()]].reset_index(drop=True)
    s = storm.iloc[0]
    d["regime"] = regimes(d, s.r_m)
    step = G.run_year(d, storm)[0]
    d["lowered"] = d.crest_m_navd88 - step["after"]

    share = d.groupby("gis").regime.value_counts(normalize=True).unstack(fill_value=0)
    share = share.reindex(columns=[k for k, _, _ in REGIMES], fill_value=0)
    low = d.groupby("gis").lowered.mean()

    apply_style()
    fig, (a, b) = plt.subplots(2, 1, figsize=figsize("double", 0.62), layout="constrained", sharex=True,
                               gridspec_kw=dict(height_ratios=(1.3, 1)))
    bottom = np.zeros(len(share))
    for key, label, col in REGIMES:
        a.bar(share.index, share[key] * 100, bottom=bottom, width=1.0, color=col, lw=0,
              edgecolor="0.7" if col == "white" else None, label=label)
        bottom += share[key].to_numpy() * 100
    a.set(xlim=(GIS[0] - 0.5, GIS[1] + 0.5), ylim=(0, 100), ylabel="Dune cells (%)")
    open_frame(a)
    _title(a, 0, f"Each dune's regime in the {s.StartTime:%d %b %Y} storm ({s.r_m:.2f} m for {s.duration:.0f} h)")

    b.bar(low.index, low, width=0.85, color=E.C_MOD, lw=0)
    b.set(ylabel="Crest lowered by\nthe module (m)", xlabel=DOMAIN_AXIS_LABEL, ylim=(0, max(1.0, low.max() * 1.1)))
    b.set_xticks(range(0, 91, 10))
    town_bands(b, where="top")
    open_frame(b)
    _title(b, 1, "Mean crest lowering per domain, this storm alone")
    handles = [Patch(facecolor=c, edgecolor="0.7" if c == "white" else c, label=l) for _, l, c in REGIMES]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=7.5)

    n = d.regime.value_counts()
    tot = int((d.regime != "nodata").sum())
    pct = {k: 100 * n.get(k, 0) / tot for k, _, _ in REGIMES if k != "nodata"}
    top = low.nlargest(3)
    caption(fig, f"**Tests:** nothing new; this figure shows where along Hatteras Island one storm puts the dunes "
            f"in each regime, and where the dune collision module acts. **How to read:** every 10 m dune cell of "
            f"GIS 1-90 (1 at Cape Point, 90 at Pea Island; {len(d)} cells, {tot} with a lidar crest) at its "
            f"{args.year - 1} lidar crest, with the toes and beach widths the hindcast reads, against the "
            f"largest {args.year} storm of the hindcast storm series ({s.StartTime:%d %b}, {s.r_m:.2f} m held "
            f"{s.duration:.0f} h; one water level along the island, as Barrier3D takes it). (a) the share of each "
            f"domain's cells in each regime: overtopped (crest at or below the water; Barrier3D's own overwash), "
            f"in collision (toe below the water, crest above it) and lowered by the module, in collision behind "
            f"a beach of {E.W0:.0f} m or more (beach factor 0: the module takes nothing), or with the water below "
            f"the toe. (b) the module's mean crest lowering per domain from this "
            f"storm alone, through its own functions at the hindcast's settings (Cs {E.CS:g}); Barrier3D's "
            f"overwash lowering is not computed; villages shaded. **Shows:** {pct['collision']:.0f}% of cells are in collision "
            f"and lowered, {pct['wide']:.0f}% in collision behind a wide beach, {pct['overwash']:.0f}% "
            f"overtopped and {pct['below']:.0f}% below the toe. The module acts on a different stretch from "
            f"overwash; its largest lowering is at GIS "
            + ", ".join(f"{g} ({v:.2f} m)" for g, v in top.items()) + ".")
    save(fig, E.OUT / "collision_8_dennis_island", close=True)
    print(f"regimes % {({k: round(v) for k, v in pct.items()})}; top lowering {top.round(2).to_dict()}")


if __name__ == "__main__":
    main()
