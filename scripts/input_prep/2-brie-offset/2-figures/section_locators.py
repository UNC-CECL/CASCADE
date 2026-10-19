"""
A locator map for each section figure: the modelled reach laid out left to right, the section boxed.

    python scripts/input_prep/2-brie-offset/2-figures/section_locators.py --year 1996

One small map per 15-domain section (GIS 1-15 ... 76-90), the same sixths as the dune line vs
shoreline figures, for putting beside them on a slide. Drawn from the island outline and the
model's domain boxes, a quarter turn with north to the right, so no network is needed.
Written to domains/locators/ inside the dune line vs shoreline comparison.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-07
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "figure_making" / "island"))

from site_layer import hat_topo_version as _tv  # noqa: E402
from site_layer.hat_figure_style import C, FIG_W_DOUBLE, INK, INK_MUTED, apply_style, save  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402
from study_area_figures import (  # noqa: E402
    HALO, LAND, LAND_EDGE, WATER_MAP, Frame, load_layers, map_axes)

# --- CONFIG ------------------------------------------------------------------
# The same sixths as compare_offset_sources.py and offset_sources_on_domains.py
SECTIONS = tuple((lo, lo + 14) for lo in range(1, 91, 15))
BOX_PAD_M = 350              # the section box stands this far off its domains
PAD_ALONG_KM, PAD_SEA_KM, PAD_SOUND_KM = 2.5, 2.6, 3.0
LABEL_SEA_M = 2300          # ocean-side labels sit this far off the reach axis
LABEL_CLEAR_M = 650         # and at least this far beyond the section box
# -----------------------------------------------------------------------------


# One locator: the reach in outline, every domain faint, the section boxed and named
def _locator(frame, dom, outline, lo, hi, window):
    fig, ax = plt.subplots(figsize=(FIG_W_DOUBLE, 1.9), constrained_layout=True)
    map_axes(ax, window)
    ax.set_facecolor(WATER_MAP)
    frame.geoms(outline.geometry).plot(ax=ax, facecolor=LAND, edgecolor=LAND_EDGE, lw=0.4, zorder=1)
    frame.geoms(dom.geometry).plot(ax=ax, facecolor="none", edgecolor="0.7", lw=0.25, zorder=2)

    sec = dom[(dom.ID >= lo) & (dom.ID <= hi)]
    frame.geoms(sec.geometry).plot(ax=ax, facecolor=C["LOCATOR"], alpha=0.25, edgecolor="none", zorder=3)
    x0, y0, x1, y1 = frame.geoms(sec.geometry).total_bounds
    ax.add_patch(Rectangle((x0 - BOX_PAD_M, y0 - BOX_PAD_M), (x1 - x0) + 2 * BOX_PAD_M,
                           (y1 - y0) + 2 * BOX_PAD_M, facecolor="none", edgecolor=C["LOCATOR"],
                           lw=2.0, zorder=6))

    # Section name on the sound side of the box, villages on the ocean side
    sea = frame.seaward
    sound_y = (y1 + BOX_PAD_M if sea[1] < 0 else y0 - BOX_PAD_M)
    ax.text((x0 + x1) / 2, sound_y + (-sea[1]) * 250, f"GIS {lo}–{hi}", ha="center",
            va="bottom" if sea[1] < 0 else "top", fontsize=11, fontweight="bold",
            color=C["LOCATOR"], zorder=8, path_effects=HALO)
    box_sea_edge = (y0 - BOX_PAD_M) if sea[1] < 0 else (y1 + BOX_PAD_M)

    # An ocean-side label, pushed past the section box where it would land under it
    def ocean_label(gis, text, colour):
        p = frame.along(gis)[0] + sea * LABEL_SEA_M
        if x0 - BOX_PAD_M - 1500 < p[0] < x1 + BOX_PAD_M + 1500:
            clear = box_sea_edge + sea[1] * LABEL_CLEAR_M
            p[1] = min(p[1], clear) if sea[1] < 0 else max(p[1], clear)
        ax.text(p[0], p[1], text, ha="center", va="center", fontsize=8, color=colour,
                zorder=7, path_effects=HALO)

    for name, (s_lo, s_hi) in HATTERAS_ANNOTATIONS.town_spans.items():
        ocean_label((s_lo + s_hi) / 2, name, INK_MUTED)
    # The two ends named on the ocean side too
    for gis in (1, 90):
        ocean_label(gis, f"GIS {gis}", INK)
    ax.text(0.995, 0.96 if sea[1] < 0 else 0.04, "north →", transform=ax.transAxes, ha="right",
            va="top" if sea[1] < 0 else "bottom", fontsize=8, color=INK)
    ax.text(0.5, 0.03 if sea[1] < 0 else 0.97, "ATLANTIC OCEAN", transform=ax.transAxes,
            ha="center", va="bottom" if sea[1] < 0 else "top", fontsize=7, color=INK_MUTED)
    return fig


# Run: one locator per section
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--year", type=int, default=1996)
    ap.add_argument("--shoreline-version", default=None, help="default: the shoreline CURRENT")
    args = ap.parse_args(argv)
    year = args.year
    shore_v = args.shoreline_version or _tv.offset_version(year, "shoreline")
    out_dir = _tv.offset_source_comparison_dir(year, "duneline_vs_shoreline", shore_v) / "domains" / "locators"

    apply_style()
    dom, outline, _roads = load_layers()
    frame = Frame(dom, outline)
    window = frame.window(dom, PAD_ALONG_KM, PAD_SEA_KM, PAD_SOUND_KM)
    for lo, hi in SECTIONS:
        fig = _locator(frame, dom, outline, lo, hi, window)
        stem = f"locator_GIS{lo:02d}-{hi:02d}"
        for path in save(fig, out_dir / stem, close=True):
            print(f"  wrote {path}")


if __name__ == "__main__":
    main()
