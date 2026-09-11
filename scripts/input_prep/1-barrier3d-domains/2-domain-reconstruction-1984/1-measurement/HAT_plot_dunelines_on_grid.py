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
from hat_topo_version import insert_figures_dir_for_domain  # noqa: E402
from hat_topo_version import duneline_shift_dir  # noqa: E402
from hat_figure_style import (apply_style, C, C_1984, C_1984_FILL,   # noqa: E402
                              C_1997, DOMAIN_AXIS_LABEL, INK, caption,
                              elevation_cmap, figsize, open_frame, save,
                              spines_for_image, town_bands, _title)

# Resolved through hat_topo_version.duneline_shift_dir - ONE definition
# of a path that eight scripts used to build by hand. Moved under
# 2-domain-reconstruction-1984/ on 2026-09-03.
S = duneline_shift_dir("1984-start")
BASE_V = "v2"   # the re-pick base; was "v3" until the 2026-09-04 renumber
BERM_EL_M = 1.7
DUNE_ROWS = 2
L84, L97, LROW0 = C_1984, C_1997, C["ROAD"]
L_DATE, L_BAND = C["REF"], C_1984_FILL
# the two relocation blocks, GIS 9-14 and 84-87: the modification under test
BLOCKS = ((9, 14), (84, 87))


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


def _blocks(ax):
    for lo, hi in BLOCKS:
        ax.axvspan(lo - .5, hi + .5, color=C["ACCENT_FILL"], alpha=.45,
                   lw=0, zorder=0)


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

    fig = plt.figure(figsize=figsize("double", height=8.2),
                     constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.3, 1.0, 1.0])

    # ---- (a) the lines on the grid --------------------------------------
    ax = fig.add_subplot(gs[0])
    cmap, norm, bounds = elevation_cmap()
    strip = np.tile(BERM_EL_M + dune[None, :topo.shape[1]], (DUNE_ROWS, 1))
    im = ax.imshow(np.vstack([strip, topo[:args.rows, :]]), cmap=cmap, norm=norm,
                   aspect="auto", interpolation="nearest", zorder=1,
                   extent=[-0.5, topo.shape[1] - 0.5,
                           args.rows - 0.5, -DUNE_ROWS - 0.5])
    ax.axhline(-0.5, color=INK, lw=0.6, zorder=4)
    ax.fill_between(prof, c84, c97, color=L_BAND, alpha=0.55, lw=0, zorder=5,
                    label="band between the lines (the date term, N)")
    ax.plot(prof, c84, "-", color=L84, lw=1.5, zorder=6, label="1984 dune line")
    ax.plot(prof, c97, "-", color=L97, lw=1.5, zorder=6, label="1997 dune line")
    ax.plot(prof, r0, "-", color=LROW0, lw=1.6, zorder=6,
            label="interior row 0")
    ax.set_xlim(-0.5, topo.shape[1] - 0.5)
    ax.set_ylim(args.rows - 0.5, -DUNE_ROWS - 0.5)
    ax.set_xlabel("alongshore cell")
    ax.set_ylabel("cross-shore cell\n(negative = dune rows)")
    _title(ax, 0, "the two dune lines on the model grid, GIS {}".format(D))
    ax.legend(loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.20),
              frameon=False, columnspacing=1.2)
    spines_for_image(ax)
    # The bar is attached to (a) and describes only that panel; (b) and (c)
    # are charts.
    cax = ax.inset_axes([1.012, 0.0, 0.016, 1.0])
    cb = fig.colorbar(im, cax=cax, boundaries=bounds[1:], ticks=bounds[1:-1])
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(length=2)
    cb.set_label("elevation (m MHW)")

    # ---- (b) the same three references, all 90 domains ------------------
    ax2 = fig.add_subplot(gs[1])
    a84 = np.array([m84[g][0] for g in gis])
    a97 = np.array([m97[g][0] for g in gis])
    ar0 = np.array([m84[g][1] for g in gis])
    town_bands(ax2)
    _blocks(ax2)
    ax2.fill_between(gis, a84, a97, color=L_BAND, alpha=0.7, lw=0, zorder=2,
                     label="band between the lines (the date term, N)")
    ax2.plot(gis, a84, "-", color=L84, lw=1.1, label="1984 dune line")
    ax2.plot(gis, a97, "-", color=L97, lw=1.1, label="1997 dune line")
    ax2.plot(gis, ar0, "-", color=LROW0, lw=1.3, label="interior row 0")
    ax2.set_ylim(50, 0)
    ax2.set_xlim(0, 91)
    ax2.set_xlabel(DOMAIN_AXIS_LABEL)
    ax2.set_ylabel("cross-shore cell\n(median per domain)")
    _title(ax2, 1, "the same three references, all domains")
    # no legend: the four handles are those of (a), whose legend sits
    # directly above this panel
    ax2.grid(axis="y")
    open_frame(ax2)

    # ---- (c) the decomposition, all 90 domains --------------------------
    ax3 = fig.add_subplot(gs[2])
    tot = np.array([m84[g][2] for g in gis])
    fea = np.array([m97[g][2] for g in gis])
    dat = np.array([mdate[g][2] for g in gis])
    town_bands(ax3, where="bottom")
    _blocks(ax3)
    ax3.plot(gis, tot, "-", color=L84, lw=1.1,
             label="row 0 − 1984 line (total)")
    ax3.plot(gis, fea, "-", color=L97, lw=1.1,
             label="row 0 − 1997 line (feature)")
    ax3.plot(gis, dat, "-", color=L_DATE, lw=1.5,
             label="1997 line − 1984 line (date, N)")
    ax3.axhline(0, color=INK, lw=0.6)
    ax3.set_xlim(0, 91)
    ax3.set_xlabel(DOMAIN_AXIS_LABEL)
    ax3.set_ylabel("offset (m)")
    # Do not overstate this. The feature term is TIGHT over most of the island
    # (IQR +14.5 to +26.2 m) but it is not constant: it spikes to 130-145 m
    # around GIS 35 and 63-68, the reaches where the date term is strongly
    # negative -- i.e. where the shoreline prograded and the two lines are on
    # opposite sides of row 0. Those are the domains where the differencing
    # argument is weakest, and the caption says so rather than averaging
    # them away.
    _title(ax3, 2, "the decomposition, all domains")
    fig.legend(*ax3.get_legend_handles_labels(), loc="outside lower center",
               ncol=3, frameon=False)
    ax3.grid(axis="y")
    open_frame(ax3)

    med_tot, med_fea, med_dat = (float(np.median(tot)), float(np.median(fea)),
                                 float(np.median(dat)))
    q_fea = (float(np.percentile(fea, 25)), float(np.percentile(fea, 75)))
    q_dat = (float(np.percentile(dat, 25)), float(np.percentile(dat, 75)))
    caption(fig,
            "The digitized dune lines in model cells, and the decomposition "
            "across the island. (a) GIS {}: the 1984 (red) and 1997 (blue) dune "
            "lines and interior row 0 (black) drawn on the Barrier3D grid as "
            "extracted from the 1996 surface (elevation in classes, m MHW; the "
            "two rows above the rule are the dune rows). Lines are drawn at the "
            "fractional cell where each geometry crosses that profile's raster "
            "row; the sawtooth is the per-profile shear of the north-up clip and "
            "appears identically in both lines, which is why it cancels in the "
            "difference. The band between the lines is the date term, N. "
            "(b) The same three references (colours as in (a)) as medians per "
            "domain over the 50 profiles, along the island (1 at Cape Point, 90 at north Pea "
            "Island). (c) The decomposition per domain: total (row 0 − 1984 "
            "line, island median {:+.1f} m), feature (row 0 − 1997 line, "
            "median {:+.1f} m, IQR {:+.1f} to {:+.1f} m) and date (1997 − 1984 "
            "line, median {:+.1f} m, IQR {:+.1f} to {:+.1f} m). The feature "
            "term is tight over most of the island but excursions at GIS 35 "
            "and 63–68, where the shoreline prograded and the two lines lie on "
            "opposite sides of row 0, are where the differencing argument is "
            "weakest. Grey bands are the villages, purple bands the two "
            "relocation blocks (GIS 9–14 and 84–87)."
            .format(D, med_tot, med_fea, q_fea[0], q_fea[1], med_dat,
                    q_dat[0], q_dat[1]))

    out = Path(args.out) if args.out else (
        insert_figures_dir_for_domain("1984-start", "1-measurement", D)
        / "HAT_dunelines_on_grid_GIS{}.png".format(D))
    save(fig, out)
    print("wrote {}".format(out))
    print("  island-wide medians:  total {:+.1f}  feature {:+.1f}  date {:+.1f} m"
          .format(med_tot, med_fea, med_dat))
    print("  feature IQR {:+.1f} to {:+.1f} m   date IQR {:+.1f} to {:+.1f} m"
          .format(q_fea[0], q_fea[1], q_dat[0], q_dat[1]))


if __name__ == "__main__":
    main()
