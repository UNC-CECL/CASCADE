#!/usr/bin/env python3
"""
The 2025 end lines of the 10-07 test figure (three of its original windows plus calendar 2025;
the DEM-date window left out at Hannah's request), as shoreline positions rather than net change (the advisor's request, 2026-10-12).

    python oct07_end_lines_positions.py  ->  figures/end_lines_2025_positions.png, supporting/

Each line is the distance from the offshore datum per GIS domain, by the island-offset build's
own intersection (coastsat_mean_shoreline_compared.datum_stations), in the six-section layout of
mean_shoreline/compared/. Labels and colours are the 10-07 test figure's. No runs, no smoothing.

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
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
MS = REPO / "scripts" / "input_prep" / "5-scr" / "1-observations" / "mean_shoreline"
sys.path[:0] = [str(REPO / "scripts"), str(MS), str(HERE)]
from coastsat_mean_shoreline_compared import _town_bands_y, datum_stations  # noqa: E402
from end_line_positions import LW, feature_lines_y, shoal_bands_y  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import mean_shoreline_csv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# key -> (window, legend label, colour), thin solid lines, calendar 2025 drawn last, on top: the 10-07 test figure's end lines with its
# labels and colours; the first three are windows of the original 10-07 figure (its DEM-date
# window, 2025-02-17 to 2026-02-17, left out), calendar 2025 the advisor's addition
LINES = {
    "cal_2024": (("2023-07-01", "2025-06-30"), "1 Jul 2024 ±1 yr", "#b2182b"),
    "jan_2025": (("2024-01-01", "2026-01-01"), "1 Jan 2025 (model end) ±1 yr", "#e08214"),
    "jan_2025_6mo": (("2024-07-01", "2025-07-01"), "1 Jan 2025 (model end) ±6 months", "#762a83"),
    "year_2025": (("2025-01-01", "2025-12-31"), "Calendar year 2025", "#1b7837"),
}
SECTIONS = ((1, 15), (16, 30), (31, 45), (46, 60), (61, 75), (76, 90))
FIG_SIZE = (13.0, 6.0)
TITLE = "CoastSat mean shoreline positions for four averaging windows around the end of the test period (2025): how much the choice of end window moves the shoreline"
STEM = "end_lines_2025_positions"
# -----------------------------------------------------------------------------


def dates(w):
    a, b = map(pd.Timestamp, w)
    return f"{a.day} {a:%b %Y}–{b.day} {b:%b %Y}"


def figure(st, n_obs, png):
    fs.apply_style()
    fig, axes = plt.subplots(1, len(SECTIONS), figsize=FIG_SIZE, layout="constrained")
    for i, ((lo, hi), ax) in enumerate(zip(SECTIONS, axes)):
        sec = st.loc[lo:hi]
        _town_bands_y(ax, lo, hi)
        shoal_bands_y(ax, lo, hi)
        for z, (k, (w, lab, c)) in enumerate(LINES.items()):
            ax.plot(sec[k], sec.index, color=c, lw=LW, zorder=3 + z,
                    label=f"Mean shoreline, {dates(w)}\n({lab})")
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
    fig.legend(h, l, loc="outside upper center", ncol=len(LINES), fontsize=7.5, frameon=False,
               columnspacing=3.0, title=TITLE, title_fontsize=10.5)

    spread = st.max(axis=1) - st.min(axis=1)
    fs.caption(fig, (
        "CoastSat mean shoreline positions for four test end lines of the 2026-10-07 figure: three "
        "windows around the end of the test plus calendar year 2025, in six consecutive sections of "
        "Hatteras, south to north from left to right, alongshore up the page. Each line is the distance "
        "from the shared offshore datum along the 100 m model transects, averaged per 500 m domain, ocean "
        "on the right, computed with the island-offset build's own intersection. CoastSat ends "
        "2026-01-13, so windows reaching past it hold data up to then. Median satellite positions per "
        "transect: " + ", ".join(f"{LINES[k][1]} {n}" for k, n in n_obs.items()) +
        f". The four lines differ by at most {spread.max():.0f} m in any domain (island median "
        f"{spread.median():.1f} m). No smoothing. Village spans shaded grey; shoal zones (Avon Shoals GIS 26-38, Wimble Shoals GIS 61-74) hatched; piers dash-dot blue, Buxton groin dotted red. Per-domain numbers in "
        f"{STEM}.csv."))
    return fs.save(fig, png, close=True)


def main():
    st = datum_stations([(k, w, None) for k, (w, _, _) in LINES.items()])
    n_obs = {k: int(pd.read_csv(mean_shoreline_csv(*w)).query("included").n_obs.median())
             for k, (w, _, _) in LINES.items()}
    out = figure(st, n_obs, HERE / "figures" / "end_lines_2025" / f"{STEM}.png")
    tab = st.rename(columns=lambda k: f"station_{k}_m")
    tab["spread_m"] = st.max(axis=1) - st.min(axis=1)
    tab.to_csv(fs.support_dir(HERE / "figures" / "end_lines_2025") / f"{STEM}.csv", float_format="%.2f")
    print(f"{len(st)} domains -> {out[0]}")
    print(f"  spread of the four lines: median {tab.spread_m.median():.1f} m, "
          f"max {tab.spread_m.max():.1f} m (GIS {tab.spread_m.idxmax()})")


if __name__ == "__main__":
    main()
