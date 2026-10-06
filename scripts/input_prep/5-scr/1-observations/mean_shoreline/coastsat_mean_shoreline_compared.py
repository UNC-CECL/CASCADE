"""
The three mean shorelines the model starts from and is graded against, drawn as the island.

    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_compared.py

The 1996, 2009 and 2025 window means of NET_CHANGE_WINDOWS, as distance from
the offshore datum per GIS domain, six sections in one row (the panels of
2-brie-offset's duneline_vs_shoreline figure). Writes mean_shoreline/compared/
and a copy to output/figures/2-observations/mean_shoreline/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-06
"""

from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))

from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    NET_CHANGE_CENTRES, NET_CHANGE_WINDOWS, mean_shoreline_csv, mean_shoreline_geojson,
)

# --- CONFIG ------------------------------------------------------------------
CALIBRATION, TEST = (1996, 2009), (2009, 2025)
# 1996 red, 2009 blue, 2025 amber (Hannah, 2026-10-06); purple stays the dune line's
STYLES = ((fs.C_1984, "-"), (fs.C_1997, (0, (4, 2))), (fs.C["ADDED"], (0, (1.2, 1.2))))
OFFSET_PRODUCER = (_REPO / "scripts" / "input_prep" / "2-brie-offset" / "1-produce"
                   / "duneline_to_raw_offsets.py")
SECTIONS = ((1, 15), (16, 30), (31, 45), (46, 60), (61, 75), (76, 90))
# One row, the island laid out left to right (Hannah, 2026-10-06)
GRID_COLS = 6
FIG_SIZE = (13.0, 5.8)
OUT_NAME = "compared"
STEM = "mean_shoreline_compared_1996_2009_2025"
# -----------------------------------------------------------------------------


# (year, window, centre) for the calibration start, the shared calibration end / test start, the test end
def shorelines():
    (c0, c1), (t0, t1) = NET_CHANGE_WINDOWS[CALIBRATION], NET_CHANGE_WINDOWS[TEST]
    if c1 != t0:
        raise RuntimeError("the calibration end window is not the test start window")
    centres = (NET_CHANGE_CENTRES[CALIBRATION][0], NET_CHANGE_CENTRES[CALIBRATION][1],
               NET_CHANGE_CENTRES[TEST][1])
    return list(zip((CALIBRATION[0], CALIBRATION[1], TEST[1]), (c0, c1, t1), centres))


# duneline_to_raw_offsets.py, loaded as a module
def _offset_producer():
    spec = importlib.util.spec_from_file_location("duneline_to_raw_offsets", OFFSET_PRODUCER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Per domain, each line's mean station from the offshore datum, by the offset build's own intersection
def datum_stations(lines):
    import geopandas as gpd
    prod = _offset_producer()
    tr = prod.load_transects()
    out = {}
    for year, (lo, hi), _ in lines:
        line = gpd.read_file(mean_shoreline_geojson(lo, hi)).to_crs(tr.crs).iloc[0].geometry
        out[year] = prod.intersect(tr, line).groupby("domain_id")["ORIG_LEN"].mean()
    st = pd.DataFrame(out)
    st.index.name = "gis_domain"
    return st


# Village spans against a vertical alongshore axis, as the duneline_vs_shoreline figure draws them
def _town_bands_y(ax, lo, hi):
    try:
        from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS
        spans = HATTERAS_ANNOTATIONS.town_spans
    except ImportError:
        return
    for name, (s_lo, s_hi) in spans.items():
        if s_hi + 0.5 < lo - 0.5 or s_lo - 0.5 > hi + 0.5:
            continue
        ax.axhspan(s_lo - 0.5, s_hi + 0.5, color="0.94", lw=0, zorder=0)
        mid = (max(s_lo - 0.5, lo - 0.5) + min(s_hi + 0.5, hi + 0.5)) / 2
        ax.text(0.015, mid, name, transform=ax.get_yaxis_transform(), ha="left",
                va="center", fontsize=6.5, color=fs.INK_MUTED, zorder=1, clip_on=True)


def figure(st, lines, n_obs, folder):
    fs.apply_style()
    years = [y for y, _, _ in lines]
    n_rows = int(np.ceil(len(SECTIONS) / GRID_COLS))
    fig, axes = plt.subplots(n_rows, GRID_COLS, figsize=FIG_SIZE, layout="constrained", squeeze=False)
    for i, ((lo, hi), ax) in enumerate(zip(SECTIONS, axes.flat)):
        sec = st.loc[lo:hi]
        _town_bands_y(ax, lo, hi)
        for (y, (wlo, whi), c), (col, ls) in zip(lines, STYLES):
            ax.plot(sec[y], sec.index, color=col, ls=ls, lw=1.4, zorder=3 + years.index(y),
                    label=f"{y} (±{'6 months' if y == years[2] else '1 yr'} of {c})")
        ax.set_ylim(lo - 0.5, hi + 0.5)
        ax.invert_xaxis()          # landward left, ocean right
        ax.xaxis.set_major_locator(mticker.MaxNLocator(3))
        ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.tick_params(axis="x", labelsize=7)
        # Each target's range of change in this section, seaward +
        d1, d2 = sec[years[0]] - sec[years[1]], sec[years[1]] - sec[years[2]]
        ax.text(0.5, 0.015, f"{years[0]}→{years[1]}: {d1.min():+.0f} to {d1.max():+.0f} m\n"
                f"{years[1]}→{years[2]}: {d2.min():+.0f} to {d2.max():+.0f} m",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=6.5, color=fs.INK_MUTED,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0), zorder=8)
        fs._title(ax, i, f"GIS {lo}–{hi}")
    fig.supylabel(fs.DOMAIN_AXIS_LABEL, fontsize=9)
    fig.supxlabel("Distance from the offshore datum (m)   —   landward ←   |   → ocean", fontsize=9)
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside upper center", ncol=3, fontsize=7, frameon=False)

    d1, d2 = st[years[0]] - st[years[1]], st[years[1]] - st[years[2]]
    fs.caption(fig, (
        f"The three CoastSat mean shorelines the hindcast starts from and is graded against, as the "
        f"island-offset build sees them, in six consecutive sections of Hatteras in one row, south to north from left to right, "
        f"alongshore up the page, south at the bottom. Each line is the distance from the shared "
        f"offshore datum along the 100 m model transects, averaged per 500 m domain, with the ocean on "
        f"the right, computed with the offset build's own intersection (duneline_to_raw_offsets.intersect). "
        f"{years[0]} (red, solid) and {years[1]} (blue, dashed) are means over ±1 yr of the middle of the "
        f"lidar flights of the start DEMs (1996 NOAA/NASA ALACE, {lines[0][1][0]} to {lines[0][1][1]}; "
        f"2009 USACE NCMP, {lines[1][1][0]} to {lines[1][1][1]}). {years[2]} (amber, dotted) has no DEM: "
        f"it is the mean over ±6 months of {lines[2][2]}, 16 yr after the 2009 centre, and CoastSat ends "
        f"2026-01-13, so it spans about 11 months; median satellite positions per transect "
        f"{n_obs[years[0]]}, {n_obs[years[1]]} and {n_obs[years[2]]}. The text in each panel gives the "
        f"range of change per domain over the calibration ({years[0]}→{years[1]}) and test "
        f"({years[1]}→{years[2]}) periods, seaward +; island medians {d1.median():+.1f} m and "
        f"{d2.median():+.1f} m. No smoothing. Village spans shaded. Per-domain numbers in "
        f"{STEM}.csv."))
    return fs.save(fig, folder / f"{STEM}.png", close=True)


def main():
    lines = shorelines()
    st = datum_stations(lines)
    years = [y for y, _, _ in lines]
    n_obs = {y: int(pd.read_csv(mean_shoreline_csv(lo, hi)).query("included").n_obs.median())
             for y, (lo, hi), _ in lines}
    folder = mean_shoreline_csv(*lines[0][1]).parent.parent / OUT_NAME
    out = figure(st, lines, n_obs, folder)
    # Stations grow landward, so earlier minus later is + where the shoreline moved seaward
    tab = st.rename(columns=lambda y: f"station_{y}_m")
    tab[f"change_{years[0]}_{years[1]}_m"] = st[years[0]] - st[years[1]]
    tab[f"change_{years[1]}_{years[2]}_m"] = st[years[1]] - st[years[2]]
    tab.to_csv(fs.support_dir(folder) / f"{STEM}.csv", float_format="%.2f")
    pub = fs.figure_dir("observations", "mean_shoreline")
    shutil.copy2(out[0], pub / out[0].name)
    print(f"{len(st)} domains -> {out[0]}\n  copy -> {pub / out[0].name}")
    for a, b in ((years[0], years[1]), (years[1], years[2])):
        d = tab[f"change_{a}_{b}_m"]
        print(f"  {a} -> {b}: median {d.median():+.1f} m, range {d.min():+.1f} to {d.max():+.1f}")


if __name__ == "__main__":
    main()
