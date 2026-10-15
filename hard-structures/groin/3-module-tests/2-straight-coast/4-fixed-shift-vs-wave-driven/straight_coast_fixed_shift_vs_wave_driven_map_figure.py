#!/usr/bin/env python3
"""
The fixed-shift groin and the wave-driven groin drawn as maps: the coast at its own angle, true scale.

    python straight_coast_fixed_shift_vs_wave_driven_map_figure.py   ->  figures/straight_coast_fixed_shift_vs_wave_driven_map.png

The comparison figure plots shoreline change against cell number, which hides the coast itself. This draws
the same runs (the year-0, 10 and 50 profiles in each version's tables/) in plan view, 1:1, north up: the
straight start coast at its orientation, the shoreline after 10 and 50 years, the sand built out or cut back
between the start and year 50, the groin, and the net drift. Reads the tables straight_coast_test.py writes;
nothing is rerun.

Geometry, from straight_coast_test.run(): cells are 500 m alongshore with the higher index to the north; the
start coast's landward position grows by tan(theta0) x 500 m per cell, so a positive angle trends west going
north and a negative one east (Buxton's coast at GIS 5-6 runs about 21 degrees east of north, near -20 here);
the net drift runs south, so the north side is updrift. Each cell's shoreline moves as one, as in BRIE.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984_FILL, INK, INK_MUTED, MAP_TEXT_DARK, apply_style, figsize, north_dart, record_caption, save,
    scale_bar_km, spines_for_image)

# --- CONFIG ------------------------------------------------------------------
RUN = HERE.parent / "1-run"                    # the test's tables
OUT = HERE / "figures" / "straight_coast_fixed_shift_vs_wave_driven_map.png"
VERSIONS = (("fixed_shift", "Fixed-shift groin", C["ADDED"]),
            ("wave_driven", "Wave-driven groin", INK))
ANGLES = ((-20, "coast turned against the drift (−20°)"), (10, "coast turned with the drift (+10°)"))
DY_M = 500.0                         # cell length alongshore, as straight_coast_test.DY_M
CELLS = 3                            # cells drawn either side of the groin
HALF_WIDTH_M = 800.0                 # map half-width east-west about the groin, the same on every map
YEARS = (10, 50)                     # shorelines drawn besides the start
GROIN_DRAWN_M = 250.0                # length drawn; in the model the groin is a face between two cells
BALANCED_M = 30.0                    # the sand-balance band used in the comparison figure
# -----------------------------------------------------------------------------


def profiles(folder):
    t = pd.read_csv(RUN / "tables" / f"straight_coast_{folder}_profiles.csv")
    return t


def by_year(folder):
    t = pd.read_csv(RUN / "tables" / f"straight_coast_{folder}_by_year.csv")
    return t


# One shoreline as a north-up polyline: cell k (offset from the groin, + north) spans northing
# (k - 0.5) x DY to (k + 0.5) x DY minus half a cell, so the groin face (between offsets 0 and 1) sits at 0
def shoreline(theta, seaward):
    xs, ys = [], []
    tan = np.tan(np.deg2rad(theta))
    for k, s in seaward.items():
        y0, y1 = (k - 1) * DY_M, k * DY_M
        xs += [-tan * y0 + s, -tan * y1 + s]
        ys += [y0, y1]
    return np.array(xs), np.array(ys)


def main():
    apply_style()
    fig = plt.figure(figsize=figsize("double", height=5.0), constrained_layout=True)
    subs = fig.subfigures(1, 2, wspace=0.03)
    axes = [ax for s in subs for ax in s.subplots(1, 2)]
    shades = {10: "0.55", 50: INK}
    results = {}
    for i, (theta, label) in enumerate(ANGLES):
        for j, (folder, title, col) in enumerate(VERSIONS):
            ax = axes[2 * i + j]
            p = profiles(folder)
            p = p[(p["theta0_deg"] == theta) & p["domain_from_groin"].between(-CELLS + 1, CELLS)]
            start = pd.Series(0.0, index=sorted(p["domain_from_groin"].unique()))
            x0, y0 = shoreline(theta, start)
            lines = {yr: shoreline(theta, p[p["year"] == yr].set_index("domain_from_groin")["seaward_m"].sort_index())
                     for yr in YEARS}
            x50, y50 = lines[YEARS[-1]]
            xmin, xmax = -HALF_WIDTH_M, HALF_WIDTH_M
            # land west of the year-50 shoreline, sea east; built out (gold) or cut back (red) since the start
            ax.fill_betweenx(y50, xmin, x50, color=C["ADDED_FILL"], alpha=0.55, lw=0, zorder=0)
            ax.fill_betweenx(y50, x50, xmax, color=C["WATER"], alpha=0.45, lw=0, zorder=0)
            ax.fill_betweenx(y0, x0, x50, where=x50 >= x0, color=C["ADDED"], alpha=0.75, lw=0, zorder=1)
            ax.fill_betweenx(y0, x50, x0, where=x50 < x0, color=C_1984_FILL, alpha=0.95, lw=0, zorder=1)
            ax.plot(x0, y0, color=INK_MUTED, lw=1.0, ls=(0, (4, 2)), zorder=3)
            for yr, (x, y) in lines.items():
                ax.plot(x, y, color=shades[yr], lw=1.0 if yr != YEARS[-1] else 1.5, zorder=4)
            # the groin, drawn shore-normal from the start shoreline at the face
            nx, ny = np.cos(np.deg2rad(theta)), np.sin(np.deg2rad(theta))
            ax.plot([-60 * nx, GROIN_DRAWN_M * nx], [-60 * ny, GROIN_DRAWN_M * ny], color=C["GROIN"], lw=2.6,
                    solid_capstyle="butt", zorder=5)
            # the net drift, offshore, running south along the coast
            tan = np.tan(np.deg2rad(theta))
            ya, yb = 2.3 * DY_M, 0.7 * DY_M
            off = 420.0
            ax.annotate("", xy=(-tan * yb + off, yb), xytext=(-tan * ya + off, ya),
                        arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.2), zorder=6)
            ax.text(-tan * 1.5 * DY_M + off + 60, 1.5 * DY_M, "net drift", rotation=(90 - theta) - (180 if theta < 0 else 0), fontsize=6.5,
                    color=INK_MUTED, va="center", ha="left")
            b = by_year(folder)
            b = b[(b["theta0_deg"] == theta) & (b["year"] == YEARS[-1])].iloc[0]
            results[folder, theta] = b
            tot = float(b["reach_change_m"])
            ok = abs(tot) <= BALANCED_M
            ax.text(0.04, 0.985, f"year 50:\nnorth side {b['updrift_seaward_m']:+.0f} m\n"
                    f"south side {b['downdrift_seaward_m']:+.0f} m", transform=ax.transAxes, fontsize=6.6, va="top",
                    color=INK, bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85), zorder=8)
            ax.text(0.04, 0.83, f"sand made or lost:\n{tot:+,.0f} m" + (" (balanced)" if ok else ""),
                    transform=ax.transAxes, fontsize=6.6, va="top", fontweight="bold",
                    color=C["REF"] if ok else C["GROIN"],
                    bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85), zorder=8)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(-(CELLS - 0.0) * DY_M, CELLS * DY_M)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            spines_for_image(ax)
            if j == 0:
                scale_bar_km(ax, length_m=500, segments=1, unit="m", x=0.08, y=0.07, text_kw=MAP_TEXT_DARK)
            if j == 1:
                north_dart(ax, (xmax - 170, -(CELLS - 0.6) * DY_M), arrow_m=220, text_kw=MAP_TEXT_DARK)
            ax.set_title(title, fontsize=8.5, color=col if j == 0 else INK, fontweight="bold")
        subs[i].suptitle(label[0].upper() + label[1:], fontsize=9)
    handles = [Line2D([], [], color=INK_MUTED, lw=1.0, ls=(0, (4, 2)), label="start: a straight coast"),
               Line2D([], [], color=shades[10], lw=1.0, label="year 10"),
               Line2D([], [], color=shades[50], lw=1.5, label="year 50"),
               Patch(facecolor=C["ADDED"], alpha=0.75, label="built out since the start"),
               Patch(facecolor=C_1984_FILL, alpha=0.95, label="cut back since the start"),
               Line2D([], [], color=C["GROIN"], lw=2.6, label="groin (drawn 250 m)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=7)
    save(fig, OUT, close=True)
    d20, c20 = results["fixed_shift", -20], results["wave_driven", -20]
    d10, c10 = results["fixed_shift", 10], results["wave_driven", 10]
    record_caption(OUT, (
        "**Tests:** nothing new; this draws the comparison figure's runs as maps, so the coast's shape can be "
        "seen. **How to read:** plan view at true scale (1:1), north up, ocean to the east, 1.5 km of coast "
        "either side of the groin. Each coast starts straight (dashed) at its orientation: −20 degrees trends "
        "east going north, close to Buxton's own coast at GIS 5-6; +10 degrees trends west going north. Grey: the "
        "shoreline after 10 years; black: after 50. Gold: sand built out between the start and year 50; red: "
        "sand cut back. Red bar: the groin, drawn shore-normal and 250 m long for visibility; in the model it is "
        "a face between two 500 m cells, and each cell's shoreline moves as one, which is why the shorelines "
        "step. The net longshore drift runs south, so the north side is updrift. Top left of each map: the cells "
        "either side of the groin at year 50; bottom: the sand gained or lost over the 61-cell reach, green when "
        f"within ±{BALANCED_M:.0f} m. BRIE's alongshore step alone under option A waves; no storms, no "
        "Barrier3D, no sea-level rise. **Shows:** the wave-driven groin builds a fillet against its "
        f"north face and cuts the south side back by the same amount on both coasts ({c20['updrift_seaward_m']:+.0f} "
        f"/ {c20['downdrift_seaward_m']:+.0f} m at −20 degrees, {c10['updrift_seaward_m']:+.0f} / "
        f"{c10['downdrift_seaward_m']:+.0f} m at +10), the shape of a real groin, with the sand balanced. The "
        "fixed-shift groin does not: on the coast turned against the drift both sides build out "
        f"({d20['updrift_seaward_m']:+.0f} / {d20['downdrift_seaward_m']:+.0f} m) and the reach gains "
        f"{d20['reach_change_m']:+,.0f} m of sand from nowhere; on the coast turned with the drift both sides "
        f"retreat ({d10['updrift_seaward_m']:+.0f} / {d10['downdrift_seaward_m']:+.0f} m) and "
        f"{abs(d10['reach_change_m']):,.0f} m is lost. That is BRIE's solve leaking at the bend the fixed-shift groin "
        "makes (`brie_sand_balance_at_a_step.png`), not the groin."))
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
