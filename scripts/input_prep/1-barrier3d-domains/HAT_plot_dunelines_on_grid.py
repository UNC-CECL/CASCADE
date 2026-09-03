#!/usr/bin/env python3
r"""
HAT_plot_dunelines_on_grid.py
==============================================================================
The digitized dune lines drawn ON the Barrier3D grid, and the same measurement
across all 90 domains.

WHY PUT THE LINES ON THE GRID
    N is a difference between two digitized lines, expressed in model cells.
    Every earlier figure showed that as numbers or as a profile. Drawing the
    lines on the cells they are measured in makes three things checkable by eye:

      * that the lines fall where the topography says a dune line should --
        the 1997 line should track the seaward face of the surveyed dune;
      * that interior row 0 sits LANDWARD of both, which is why the raw
        measurement carries a definitional offset at all;
      * that the band between the two lines -- the date term, which is N -- is
        a coherent alongshore feature and not per-profile noise.

    Both lines are drawn per profile, at the fractional cell where the geometry
    actually crosses that profile's raster row, so the sawtooth is real: it is
    the per-profile shear of the north-up clip, and it is present in both lines
    identically, which is why it cancels in the difference.

USAGE
    python HAT_plot_dunelines_on_grid.py [--domain 85]
==============================================================================
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import array_name, dune_topo_root            # noqa: E402
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_topo_version import duneline_shift_dir  # noqa: E402
from hat_figure_style import (apply_style, C, caption,             # noqa: E402
                              elevation_cmap, panel_title,
                              spines_for_image)

# Resolved through hat_topo_version.duneline_shift_dir - ONE definition
# of a path that eight scripts used to build by hand. Moved under
# row-insert-scope/ on 2026-09-03.
S = duneline_shift_dir("1984-start")
BASE_V = "v3"
BERM_EL_M = 1.7
DUNE_ROWS = 2
L84, L97, LROW0 = C["ACCENT"], "#1b6ca8", "#111111"


def per_profile(fname, D):
    out = {}
    for r in csv.DictReader(open(S / fname)):
        if int(r["domain"]) == D:
            out[int(r["profile"])] = (float(r["duneline_cell"]),
                                      int(r["interior_row0_cell"]))
    return out


def domain_medians(fname):
    out = {}
    for r in csv.DictReader(open(S / fname)):
        out[int(r["domain"])] = (float(r["duneline_cell_median"]),
                                 float(r["row0_cell_median"]),
                                 float(r["shift_m_median"]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--rows", type=int, default=26)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    D = args.domain
    apply_style()

    R = dune_topo_root("1984-start")
    topo = np.load(R / BASE_V / "topography" / array_name("topography", D)) * 10.0
    dune = np.load(R / BASE_V / "dunes" / array_name("dune", D)) * 10.0

    p84, p97 = per_profile("duneline_shift_1984_profiles.csv", D), \
        per_profile("duneline_shift_1997_profiles.csv", D)
    common = sorted(set(p84) & set(p97))
    prof = np.array(common, dtype=float)
    c84 = np.array([p84[k][0] for k in common])
    c97 = np.array([p97[k][0] for k in common])
    r0 = np.array([p84[k][1] for k in common], dtype=float)

    m84, m97 = domain_medians("duneline_shift_1984.csv"), \
        domain_medians("duneline_shift_1997.csv")
    mdate = domain_medians("duneline_retreat_1984_1997.csv")
    gis = np.array(sorted(set(m84) & set(m97) & set(mdate)))

    fig = plt.figure(figsize=(11.0, 10.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.35, 1.0, 1.0])

    # ---- (a) the lines on the grid --------------------------------------
    ax = fig.add_subplot(gs[0])
    cmap, norm, bounds = elevation_cmap()
    strip = np.tile(BERM_EL_M + dune[None, :topo.shape[1]], (DUNE_ROWS, 1))
    ax.imshow(np.vstack([strip, topo[:args.rows, :]]), cmap=cmap, norm=norm,
              aspect="auto", interpolation="nearest", zorder=1,
              extent=[-0.5, topo.shape[1] - 0.5,
                      args.rows - 0.5, -DUNE_ROWS - 0.5])
    ax.axhline(-0.5, color="#333333", lw=0.9, zorder=4)
    ax.fill_between(prof, c84, c97, color=C["ADDED"], alpha=0.30, zorder=5,
                    label="date term = N")
    ax.plot(prof, c84, "-", color=L84, lw=2.0, zorder=6, label="1984 dune line")
    ax.plot(prof, c97, "-", color=L97, lw=2.0, zorder=6, label="1997 dune line")
    ax.plot(prof, r0, "-", color=LROW0, lw=2.2, zorder=6,
            label="interior row 0")
    ax.set_xlim(-0.5, topo.shape[1] - 0.5)
    ax.set_ylim(args.rows - 0.5, -DUNE_ROWS - 0.5)
    ax.set_xlabel("alongshore cell")
    ax.set_ylabel("cross-shore cell\n(negative = dune rows)")
    ax.set_title(panel_title("a", "GIS {} — the two dune lines drawn on the "
                                  "Barrier3D grid".format(D)))
    ax.legend(fontsize=7.4, loc="upper center", ncol=4,
              bbox_to_anchor=(0.5, -0.20), frameon=False)
    spines_for_image(ax)

    # ---- (b) the same three references, all 90 domains ------------------
    ax2 = fig.add_subplot(gs[1])
    a84 = np.array([m84[g][0] for g in gis])
    a97 = np.array([m97[g][0] for g in gis])
    ar0 = np.array([m84[g][1] for g in gis])
    ax2.fill_between(gis, a84, a97, color=C["ADDED"], alpha=0.30, zorder=2,
                     label="date term = N")
    ax2.plot(gis, a84, "-", color=L84, lw=1.5, label="1984 dune line")
    ax2.plot(gis, a97, "-", color=L97, lw=1.5, label="1997 dune line")
    ax2.plot(gis, ar0, "-", color=LROW0, lw=1.8, label="interior row 0")
    for lo, hi in ((9, 14), (84, 87)):
        ax2.axvspan(lo - .5, hi + .5, color="#ffe9b0", alpha=.5, zorder=0)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 91)
    ax2.set_xlabel("GIS domain")
    ax2.set_ylabel("cross-shore cell\n(median per domain)")
    ax2.set_title(panel_title("b", "All 90 domains, same three references "
                                   "(shaded = the two relocation blocks)"))
    ax2.legend(fontsize=7.6, ncol=4, loc="upper left")
    ax2.grid(alpha=.3)

    # ---- (c) the decomposition, all 90 domains --------------------------
    ax3 = fig.add_subplot(gs[2])
    tot = np.array([m84[g][2] for g in gis])
    fea = np.array([m97[g][2] for g in gis])
    dat = np.array([mdate[g][2] for g in gis])
    ax3.plot(gis, tot, "-", color=L84, lw=1.4,
             label="total  row 0 − 1984  (median {:+.1f} m)".format(
                 float(np.median(tot))))
    ax3.plot(gis, fea, "-", color=L97, lw=1.4,
             label="feature  row 0 − 1997  (median {:+.1f} m)".format(
                 float(np.median(fea))))
    ax3.plot(gis, dat, "-", color=C["REF"], lw=1.9,
             label="DATE = N  (median {:+.1f} m)".format(float(np.median(dat))))
    ax3.axhline(0, color="k", lw=0.9)
    for lo, hi in ((9, 14), (84, 87)):
        ax3.axvspan(lo - .5, hi + .5, color="#ffe9b0", alpha=.5, zorder=0)
    ax3.set_xlim(0, 91)
    ax3.set_xlabel("GIS domain   (1 = Cape Point  →  90 = north Pea Island)")
    ax3.set_ylabel("metres")
    # Do not overstate this. The feature term is TIGHT over most of the island
    # (IQR +14.5 to +26.2 m) but it is not constant: it spikes to 130-145 m
    # around GIS 35 and 63-68, the reaches where the date term is strongly
    # negative -- i.e. where the shoreline prograded and the two lines are on
    # opposite sides of row 0. Those are the domains where the differencing
    # argument is weakest, and the figure should say so rather than average
    # them away.
    ax3.set_title(panel_title("c", "The decomposition island-wide — the feature "
                                   "term is tight over most of the island, "
                                   "with excursions at GIS 35 and 63–68"))
    ax3.legend(fontsize=7.6, ncol=3, loc="upper left")
    ax3.grid(alpha=.3)

    fig_h = fig.get_figheight()
    # Attached to panel (a), vertically. It describes ONLY that panel, and a
    # figure-level bar at the bottom kept landing on the caption while implying
    # it applied to (b) and (c) as well.
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax,
                      orientation="vertical", boundaries=bounds[1:],
                      ticks=bounds[1:-1], pad=0.012, fraction=0.030,
                      aspect=18)
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(labelsize=6.8, length=2)
    cb.set_label("elevation (m MHW)\npanel (a) only", fontsize=7.0,
                 labelpad=3)

    fig.suptitle("The digitized dune lines in model cells, and the "
                 "decomposition across the island",
                 fontsize=11, fontweight="bold", x=0.055, ha="left",
                 y=1 - 0.24 / fig_h)
    caption(fig,
            "Lines are drawn at the fractional cell where each geometry crosses "
            "that profile's raster row. The sawtooth in (a) is the per-profile "
            "shear of the north-up clip;\n"
            "it appears identically in both lines, which is why it cancels in "
            "the difference. (b) and (c) are medians per domain over the 50 "
            "profiles.",
            y=0.06 / fig_h, size=7.4)
    fig.subplots_adjust(top=1 - 0.58 / fig_h, bottom=0.62 / fig_h,
                        left=0.085, right=0.985, hspace=0.62)

    out = Path(args.out) if args.out else (
        insert_figures_dir("1984-start") / "HAT_dunelines_on_grid_GIS{}.png".format(D))
    fig.savefig(out)
    print("wrote {}".format(out))
    print("  island-wide medians:  total {:+.1f}  feature {:+.1f}  date {:+.1f} m"
          .format(float(np.median(tot)), float(np.median(fea)),
                  float(np.median(dat))))
    print("  feature IQR {:+.1f} to {:+.1f} m   date IQR {:+.1f} to {:+.1f} m"
          .format(float(np.percentile(fea, 25)), float(np.percentile(fea, 75)),
                  float(np.percentile(dat, 25)), float(np.percentile(dat, 75))))


if __name__ == "__main__":
    main()
