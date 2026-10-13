#!/usr/bin/env python3
"""
The three calendar-year mean shorelines of the advisor's 2026-10-12 proposal, as shoreline positions.

    python end_line_positions.py  ->  figures/calendar_year_shorelines_1996_2009_2025.png, supporting/

Advisor's proposal: each period starts and ends on a calendar-year CoastSat mean, and storms run
in the years between. Calibration 1996 -> 2009 (storms 1997-2009), test 2009 -> 2025 (storms
2010-2025). Each line is the distance from the offshore datum per GIS domain, by the
island-offset build's own intersection (coastsat_mean_shoreline_compared.datum_stations), in the
six-section layout of mean_shoreline/compared/. No runs, no smoothing.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
MS = REPO / "scripts" / "input_prep" / "5-scr" / "1-observations" / "mean_shoreline"
sys.path[:0] = [str(REPO / "scripts"), str(MS)]
from coastsat_mean_shoreline_compared import _town_bands_y, datum_stations  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import mean_shoreline_csv  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# year -> (window, role in the proposal, colour); 1996 red, 2009 amber, 2025 blue, as on every
# figure that draws all three (Hannah, 2026-10-12). All built at min 5 positions per transect. Thin solid lines, the
# last drawn on top (Hannah, 2026-10-12).
LINES = {
    1996: (("1996-01-01", "1996-12-31"), "start of calibration", fs.C_1984),
    2009: (("2009-01-01", "2009-12-31"), "end of calibration and start of test", fs.C["ADDED"]),
    2025: (("2025-01-01", "2025-12-31"), "end of test", fs.C_1997),
}
SECTIONS = ((1, 15), (16, 30), (31, 45), (46, 60), (61, 75), (76, 90))
LW = 1.0
FIG_SIZE = (13.0, 6.0)
TITLE = "CoastSat mean shoreline positions, calendar years 1996, 2009 and 2025: the start and end lines of the proposed calibration (1996→2009) and test (2009→2025) periods"
STEM = "calendar_year_shorelines_1996_2009_2025"
# -----------------------------------------------------------------------------


# Shoal zones against a vertical alongshore axis, hatched as the net-change figures draw them;
# the name sits at the top of the band's visible part, on the landward side
def shoal_bands_y(ax, lo, hi):
    a = HATTERAS_ANNOTATIONS
    for name, (s_lo, s_hi) in a.shoal_zones.items():
        if s_hi + 0.5 <= lo - 0.5 or s_lo - 0.5 >= hi + 0.5:
            continue
        ax.axhspan(s_lo - 0.5, s_hi + 0.5, color=a.color_shoal, alpha=0.10, zorder=0,
                   hatch="///", edgecolor=a.color_shoal, linewidth=0)
        top = min(s_hi + 0.5, hi + 0.5) - 0.15
        ax.text(0.015, top, name, transform=ax.get_yaxis_transform(), ha="left", va="top",
                fontsize=6.5, color="#7A5800", style="italic", zorder=1, clip_on=True,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))


# Piers and groins as horizontal lines at their GIS position, styled as the net-change figures
# draw them; call after the shorelines are drawn, so each name goes on the side away from them
def feature_lines_y(ax, lo, hi, st):
    a = HATTERAS_ANNOTATIONS
    lines = [(n, d, a.color_pier, "-.", 1.0) for n, (d, _) in a.piers.items()]
    lines += [(n, d, a.color_groin, ":", 1.1) for n, d in a.groins.items()]
    for name, dom, c, ls, lw in lines:
        if not lo - 0.5 <= dom <= hi + 0.5:
            continue
        ax.axhline(dom, color=c, ls=ls, lw=lw, alpha=0.85, zorder=2)
        x = np.interp(dom, st.index, st.mean(axis=1))
        frac = ax.transAxes.inverted().transform(ax.transData.transform((x, dom)))[0]
        right = frac < 0.5
        # Village names sit at the left, mid-span (_town_bands_y); go below the line if one is just above
        mids = [(max(t_lo - 0.5, lo - 0.5) + min(t_hi + 0.5, hi + 0.5)) / 2
                for t_lo, t_hi in a.town_spans.values()
                if not (t_hi + 0.5 < lo - 0.5 or t_lo - 0.5 > hi + 0.5)]
        below = not right and any(0 < m - dom < 1.0 for m in mids)
        ax.text(0.985 if right else 0.015, dom - 0.1 if below else dom + 0.1, name,
                transform=ax.get_yaxis_transform(),
                ha="right" if right else "left", va="top" if below else "bottom", fontsize=6.5, color=c, zorder=2,
                clip_on=True, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))


def figure(st, n_obs, png):
    fs.apply_style()
    fig, axes = plt.subplots(1, len(SECTIONS), figsize=FIG_SIZE, layout="constrained")
    for i, ((lo, hi), ax) in enumerate(zip(SECTIONS, axes)):
        sec = st.loc[lo:hi]
        _town_bands_y(ax, lo, hi)
        shoal_bands_y(ax, lo, hi)
        for z, (y, (_, role, c)) in enumerate(LINES.items()):
            ax.plot(sec[y], sec.index, color=c, lw=LW, zorder=3 + z,
                    label=f"{y} mean shoreline, 1 Jan–31 Dec {y}\n({role})")
        ax.set_ylim(lo - 0.5, hi + 0.5)
        ax.invert_xaxis()          # landward left, ocean right
        ax.xaxis.set_major_locator(mticker.MaxNLocator(3))
        ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.tick_params(axis="x", labelsize=7)
        feature_lines_y(ax, lo, hi, st)
        fs._title(ax, i, f"GIS {lo}–{hi}")
    fig.supylabel(fs.DOMAIN_AXIS_LABEL, fontsize=9)
    fig.supxlabel("Distance from the offshore datum (m)   —   landward ←   |   → ocean", fontsize=9)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside upper center", ncol=3, fontsize=7.5, frameon=False,
               columnspacing=3.0, title=TITLE, title_fontsize=10.5)

    d1, d2 = st[1996] - st[2009], st[2009] - st[2025]
    fs.caption(fig, (
        "CoastSat mean shorelines for calendar years 1996 (red), 2009 (amber) and 2025 "
        "(blue, drawn on top): each transect's mean over 1 Jan to 31 Dec of that year, min 5 positions. "
        "Under the advisor's 2026-10-12 proposal, 1996 starts calibration (storms 1997-2009), 2009 "
        "ends calibration and starts the test (storms 2010-2025), and 2025 ends the test. Six "
        "consecutive sections of Hatteras, south to north from left to right, alongshore up the page. "
        "Each line is the distance from the shared offshore datum along the 100 m model transects, "
        "averaged per 500 m domain, ocean on the right, computed with the island-offset build's own "
        f"intersection. Median satellite positions per transect: 1996 {n_obs[1996]}, 2009 "
        f"{n_obs[2009]}, 2025 {n_obs[2025]}. Change per domain, seaward +, island medians 1996→2009 {d1.median():+.1f} m, 2009→2025 "
        f"{d2.median():+.1f} m. No smoothing. Village spans shaded grey; shoal zones (Avon Shoals GIS 26-38, Wimble Shoals GIS 61-74) hatched; piers dash-dot blue, Buxton groin dotted red. Per-domain numbers in {STEM}.csv."))
    return fs.save(fig, png, close=True)


def main():
    st = datum_stations([(y, w, None) for y, (w, _, _) in LINES.items()])
    n_obs = {y: int(pd.read_csv(mean_shoreline_csv(*w)).query("included").n_obs.median())
             for y, (w, _, _) in LINES.items()}
    out = figure(st, n_obs, HERE / "figures" / "calendar_year_shorelines" / f"{STEM}.png")
    tab = st.rename(columns=lambda y: f"station_{y}_m")
    tab["change_1996_2009_m"] = st[1996] - st[2009]
    tab["change_2009_2025_m"] = st[2009] - st[2025]
    tab.to_csv(fs.support_dir(HERE / "figures" / "calendar_year_shorelines") / f"{STEM}.csv", float_format="%.2f")
    print(f"{len(st)} domains -> {out[0]}  (median positions {n_obs})")
    for c in ("change_1996_2009_m", "change_2009_2025_m"):
        d = tab[c]
        print(f"  {c}: median {d.median():+.1f} m, range {d.min():+.1f} to {d.max():+.1f}")


if __name__ == "__main__":
    main()
