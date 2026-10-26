#!/usr/bin/env python3
"""
The fixed-shift groin and the wave-driven groin drawn as maps: the coast at its own angle, true scale.

    python straight_coast_fixed_shift_vs_wave_driven_map_figure.py   ->  figures/straight_coast_fixed_shift_vs_wave_driven_map.png
    python straight_coast_fixed_shift_vs_wave_driven_map_figure.py --module fixed-shift
                                                 ->  ../5-fixed-shift/figures/straight_coast_fixed_shift_map.png
    python straight_coast_fixed_shift_vs_wave_driven_map_figure.py --module wave-driven
                                                 ->  ../6-wave-driven/figures/straight_coast_wave_driven_map.png

The comparison figure plots shoreline change against cell number, which hides the coast itself. This draws
the same runs (the year-0, 10 and 50 profiles in each version's tables/) in plan view, 1:1, north up: the
straight start coast at its orientation, the shoreline after 10 and 50 years, the sand built out or cut back
between the start and year 50, the groin, and the net drift, with the wave climate every coast shares.
--module maps one version alone. Reads the tables straight_coast_test.py writes; nothing is rerun.

Geometry, from straight_coast_test.run(): cells are 500 m alongshore with the higher index to the north; the
start coast's landward position grows by tan(theta0) x 500 m per cell, so a positive angle trends west going
north and a negative one east (Buxton's coast at GIS 5-6 runs about 21 degrees east of north, near -20 here);
the net drift runs south, so the north side is updrift. Each cell's shoreline moves as one, as in BRIE, so each
is drawn as a rectangle along the coast moved shore-normal (map_geometry.py).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-16
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import argparse
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
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(HERE.parent / "1-run"))   # the test, for its wave climate
import straight_coast_test as sc  # noqa: E402
from map_geometry import draw_coast, shoreline  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984_FILL, INK, INK_MUTED, MAP_TEXT_DARK, apply_style, figsize, north_dart, record_caption, save,
    scale_bar_km, spines_for_image)

# --- CONFIG ------------------------------------------------------------------
RUN = workspace_twin(HERE.parent) / "1-run"                    # the test's tables
OUT = workspace_twin(HERE) / "figures" / "straight_coast_fixed_shift_vs_wave_driven_map.png"
VERSIONS = (("fixed_shift", "Fixed-shift groin", C["ADDED"]),
            ("wave_driven", "Wave-driven groin", INK))
ALONE = {                            # --module: the version mapped alone, and its own folder
    "fixed-shift": ("fixed_shift", workspace_twin(HERE.parent) / "5-fixed-shift" / "figures" / "straight_coast_fixed_shift_map.png"),
    "wave-driven": ("wave_driven", workspace_twin(HERE.parent) / "6-wave-driven" / "figures" / "straight_coast_wave_driven_map.png"),
}
ANGLES = ((-20, "coast turned against the drift (−20°)"), (10, "coast turned with the drift (+10°)"))
DY_M = 500.0                         # cell length alongshore, as straight_coast_test.DY_M
CELLS = 3                            # cells in view either side of the groin
DRAWN = CELLS + 3                    # cells drawn, so a turned coast still fills the view
FILLS = dict(land=dict(facecolor=C["ADDED_FILL"], alpha=0.55), sea=dict(facecolor=C["WATER"], alpha=0.45),
             built=dict(facecolor=C["ADDED"], alpha=0.75), cut=dict(facecolor=C_1984_FILL, alpha=0.95))
HALF_WIDTH_M = 800.0                 # map half-width east-west about the groin, the same on every map
YEARS = (10, 50)                     # shorelines drawn besides the start
GROIN_DRAWN_M = 250.0                # length drawn; in the model the groin is a face between two cells
WAVES = ("Waves on every coast: Hs {Hs:g} m, Tp {Tp:g} s, asymmetry {asym:g}, "
         "high-angle fraction {ahf:g}").format(**sc.CLIMATE)
# -----------------------------------------------------------------------------


def profiles(folder):
    t = pd.read_csv(workspace_twin(RUN) / "tables" / f"straight_coast_{folder}_profiles.csv")
    return t


def by_year(folder):
    t = pd.read_csv(workspace_twin(RUN) / "tables" / f"straight_coast_{folder}_by_year.csv")
    return t


def main(module="both"):
    apply_style()
    if module == "both":
        versions, out = VERSIONS, OUT
    else:
        versions, out = tuple(v for v in VERSIONS if v[0] == ALONE[module][0]), ALONE[module][1]
    n = len(versions)
    fig = plt.figure(figsize=figsize("double" if n == 2 else "single", height=5.0), constrained_layout=True)
    subs = fig.subfigures(1, 2, wspace=0.03)
    axes = [ax for s in subs for ax in np.atleast_1d(s.subplots(1, n))]
    shades = {10: "0.55", 50: INK}
    results = {}
    for i, (theta, label) in enumerate(ANGLES):
        for j, (folder, title, col) in enumerate(versions):
            ax = axes[n * i + j]
            p = profiles(folder)
            p = p[(p["theta0_deg"] == theta) & p["domain_from_groin"].between(-DRAWN + 1, DRAWN)]
            seaward = {yr: p[p["year"] == yr].set_index("domain_from_groin")["seaward_m"].sort_index()
                       for yr in YEARS}
            offs = seaward[YEARS[-1]].index.to_numpy()
            x0, y0 = shoreline(theta, offs, np.zeros(len(offs)), DY_M)
            lines = {yr: shoreline(theta, offs, v.to_numpy(), DY_M) for yr, v in seaward.items()}
            xmin, xmax = -HALF_WIDTH_M, HALF_WIDTH_M
            # land behind the year-50 shoreline, sea in front; built out (gold) or cut back (red) since the start
            draw_coast(ax, theta, offs, seaward[YEARS[-1]].to_numpy(), DY_M, **FILLS)
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
            ax.text(0.04, 0.985, f"year 50:\nnorth side {b['updrift_seaward_m']:+.0f} m\n"
                    f"south side {b['downdrift_seaward_m']:+.0f} m", transform=ax.transAxes, fontsize=6.6, va="top",
                    color=INK, bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85), zorder=8)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(-(CELLS - 0.0) * DY_M, CELLS * DY_M)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            spines_for_image(ax)
            if j == 0:
                scale_bar_km(ax, length_m=500, segments=1, unit="m", x=0.08, y=0.07, text_kw=MAP_TEXT_DARK)
            if j == n - 1:
                north_dart(ax, (xmax - 170, -(CELLS - 0.6) * DY_M), arrow_m=220, text_kw=MAP_TEXT_DARK)
            if n == 2:
                ax.set_title(title, fontsize=8.5, color=col if j == 0 else INK, fontweight="bold")
        head = label[0].upper() + label[1:]
        subs[i].suptitle(head if n == 2 else head.replace(" (", "\n("), fontsize=9)
    handles = [Line2D([], [], color=INK_MUTED, lw=1.0, ls=(0, (4, 2)), label="start: a straight coast"),
               Line2D([], [], color=shades[10], lw=1.0, label="year 10"),
               Line2D([], [], color=shades[50], lw=1.5, label="year 50"),
               Patch(facecolor=C["ADDED"], alpha=0.75, label="built out since the start"),
               Patch(facecolor=C_1984_FILL, alpha=0.95, label="cut back since the start"),
               Line2D([], [], color=C["GROIN"], lw=2.6, label="groin (drawn 250 m)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3 if n == 2 else 2, frameon=False, fontsize=7,
               title=WAVES if n == 2 else WAVES.replace(": ", ":\n"), title_fontsize=7)
    if n == 1:
        fig.suptitle(versions[0][1], fontsize=10, fontweight="bold")
    save(fig, out, close=True)
    if n == 1:
        record_caption(out, alone_caption(versions[0][0], results))
        print(f"wrote {out.relative_to(REPO)}")
        return
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
        "either side of the groin at year 50. Below the key: the wave climate, the same on every coast. BRIE's "
        "alongshore step alone; no storms, no Barrier3D, no sea-level rise. **Shows:** the wave-driven groin builds a fillet against its "
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


# The caption of one version mapped alone (--module)
def alone_caption(folder, results):
    a, b = results[folder, -20], results[folder, 10]
    how = ("**How to read:** plan view at true scale (1:1), north up, ocean to the east, 1.5 km of coast either "
           "side of the groin; left the coast at −20 degrees (trends east going north, close to Buxton's own coast "
           "at GIS 5-6), right +10 degrees (trends west going north). Dashed: the straight start coast. Grey: the "
           "shoreline after 10 years; black: after 50. Gold: sand built out between the start and year 50; red: "
           "sand cut back. Red bar: the groin, drawn shore-normal and 250 m long for visibility; in the model it "
           "is a face between two 500 m cells, and each cell's shoreline moves as one, which is why the shorelines "
           "step. The net longshore drift runs south, so the north side is updrift. Top left of each map: the "
           "cells either side of the groin at year 50. Below the key: the wave climate, the same on both coasts. "
           "BRIE's alongshore step alone; no storms, no Barrier3D, no sea-level rise. The two versions side by side are in "
           "`straight_coast_fixed_shift_vs_wave_driven_map.png`. ")
    if folder == "wave_driven":
        return ("**Tests:** nothing new; the wave-driven groin's straight-coast runs drawn as maps, so the coast's "
                "shape can be seen. " + how + "**Shows:** a fillet against the groin's north face and the south "
                f"side cut back by the same amount on both coasts ({a['updrift_seaward_m']:+.0f} / "
                f"{a['downdrift_seaward_m']:+.0f} m at −20 degrees, {b['updrift_seaward_m']:+.0f} / "
                f"{b['downdrift_seaward_m']:+.0f} m at +10), the shape of a real groin, with the sand balanced.")
    return ("**Tests:** nothing new; the fixed-shift groin's straight-coast runs drawn as maps, so the coast's "
            "shape can be seen. " + how + "**Shows:** not the shape of a groin: on the coast turned against the "
            f"drift both sides build out ({a['updrift_seaward_m']:+.0f} / {a['downdrift_seaward_m']:+.0f} m) and "
            f"the reach gains {a['reach_change_m']:+,.0f} m of sand from nowhere; on the coast turned with the "
            f"drift both sides retreat ({b['updrift_seaward_m']:+.0f} / {b['downdrift_seaward_m']:+.0f} m) and "
            f"{abs(b['reach_change_m']):,.0f} m is lost. That is BRIE's solve leaking at the bend the groin makes "
            "(`brie_sand_balance_at_a_step.png`), not the groin.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--module", choices=("both", "fixed-shift", "wave-driven"), default="both")
    main(ap.parse_args().module)
