#!/usr/bin/env python3
r"""
HAT_plot_fill_options_grid.py
==============================================================================
The candidate interior fills as BARRIER3D DOMAIN VIEWS -- the grid the model
would be handed under each, not a profile line through it.

A median profile hides the thing that decides overwash: whether the added rows
are uniform across the domain or carry alongshore structure. The flat fill has
none by construction; the matched backdune carries today's backdune texture;
the measured fill inherits the 1996 dune face and is as variable as that dune
was. The grid shows that directly.

    reference          v3, no rows added. NC-12 at its MEASURED offset.
    flat backdune      median of interior rows 1-3, per column
    matched backdune   today's near-dune profile copied to the 1984 position
    measured + floor   the real DEM cell where dry, floored at the platform.
                       THE SHIPPED RULE -- v4 at the ten block domains, v5
                       island-wide.
    measured + median  every dry cell kept as measured; only the cells at or
                       below MHW filled, with the median of the block's own
                       dry cells. `--fill median`.
    raw DEM, no floor  A CONTROL, NOT A CANDIDATE.

`taper` was drawn here until 2026-09-03 and is gone: it was fully invented and
it anchored on row 0, which at GIS 85 IS the mis-picked 1996 crest, so it
inherited a known-bad endpoint. `--fill taper` still exists in
HAT_insert_seaward_rows.py - removing a build capability is a different
decision from removing a figure panel.

`matched backdune` is NOT a `--fill` choice in HAT_insert_seaward_rows.py. It
is drawn here as a candidate; building it would need a new fill rule.

THE ROAD IS DRAWN AT ITS MEASURED OFFSET, NOT THE FLOORED ONE
    int() truncates toward zero, so a measured -15 m gives road_start = -1: the
    road lands SEAWARD of interior row 0, drawn overlapping the dune strip and
    hatched. That is not a plotting artefact, it is the failure --
    roadway_manager would evaluate xyz_interior_grid[-1:1, :], valid Python
    indexing from the LANDWARD end, and bulldoze the sound-side marsh. Flooring
    it in the figure is what made the problem invisible in earlier versions.

USAGE
    python HAT_plot_fill_options_grid.py [--domain 85] [--rows 26]
==============================================================================
"""

from __future__ import annotations

import argparse
import csv
import importlib.util as _iu
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import array_name, dune_topo_root            # noqa: E402
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_figure_style import (apply_style, C, caption,             # noqa: E402
                              elevation_cmap, panel_title,
                              spines_for_image)

# INS_V supplies N and the post-insert setback ONLY; the figure draws
# blocks it builds itself. Repointed v4 -> v5 on 2026-09-03: N is
# identical at all ten block domains (verified), and v5 is the version
# taken forward, so v4 no longer has to exist for this figure to build.
BASE_V, INS_V = "v3", "v5"
OFFSET_SCRIPT = (REPO / "scripts/input_prep/4-mgmt-forcings/road_offset"
                 / "1-produce/HAT_road_offset_from_dune_start.py")
BERM_EL_M = 1.7
DUNE_ROWS = 2
ROAD_CELLS = 2
BACKDUNE_ROWS = 3

# THE ROAD ELEVATION CASCADE ACTUALLY USES, m MHW.
#
# hatteras_site_config.HATTERAS_ROAD_ELEVATION_FILE resolves to THIS file.
# There is a second one, road_offset/dunestart_offset/1984/
# RoadElevation_1984_dunestart.csv, and it must NOT be used here: it samples
# along the 1984 alignment, which at the relocated domains (GIS 9-15, 84-87)
# now lies UNDER the foredune, so it returns dune rather than roadbed. At
# GIS 85 the two read 0.807 m and 1.833 m - a metre of difference that is
# entirely the abandoned corridor being buried. See the long note beside
# HATTERAS_ROAD_ELEVATION_FILE in hatteras_site_config.py.
ROAD_ELEV_FILE = (REPO / "data" / "hatteras_init" / "4-mgmt-forcing"
                  / "road_elevation" / "RoadElevation.csv")


def road_elevation(domain):
    """Per-domain road elevation, m MHW. None if this domain has no NC-12."""
    rows = list(csv.reader(open(ROAD_ELEV_FILE)))
    ids = [int(float(x)) for x in rows[0]]
    vals = [float(x) for x in rows[1]]
    return dict(zip(ids, vals)).get(domain)


def load_ext():
    spec = _iu.spec_from_file_location("hat_off", OFFSET_SCRIPT)
    m = _iu.module_from_spec(spec)
    sys.modules["hat_off"] = m
    spec.loader.exec_module(m)
    return m.load_extractor("1984-start")


def real_cells(dom, row0, n, n_along):
    """The DEM cell at each added position, per profile. NaN where off-array."""
    z = dom["z"]
    out = np.full((n, n_along), np.nan)
    for i in range(n_along):
        for k in range(n):
            src = int(row0[i]) - n + k
            if 0 <= src < z.shape[1]:
                out[k, i] = z[i, src]
    return out


def draw(ax, topo, dune, setback, n_added, nrows, letter, title, note,
         reference=False, road_ele=None):
    """`road_ele` paints the road block at a SINGLE elevation, which is what
    the model holds: bulldoze() does `np.zeros(...) + road_ele`, so after the
    first pass every road cell carries the same value regardless of the ground
    that was there. Left as None for the reference panel, whose setback is
    negative and whose road therefore falls outside the interior array."""
    cmap, norm, _ = elevation_cmap()
    n_along = topo.shape[1]
    strip = np.tile(BERM_EL_M + dune[None, :n_along], (DUNE_ROWS, 1))
    shown = np.vstack([strip, topo[:nrows, :]])
    if road_ele is not None and setback >= 0:
        r0 = int(setback / 10.0) + DUNE_ROWS
        shown[r0:r0 + ROAD_CELLS, :] = road_ele
    ax.imshow(shown, cmap=cmap, norm=norm,
              aspect="auto", interpolation="nearest", zorder=1)

    ax.axhline(DUNE_ROWS - 0.5, color="#333333", lw=0.9, zorder=4)
    if n_added:
        ax.add_patch(Rectangle((-0.5, DUNE_ROWS - 0.5), n_along, n_added,
                               fill=False, ec=C["ADDED"], lw=1.4,
                               ls=(0, (4, 2)), zorder=5))

    rs = int(setback / 10.0) + DUNE_ROWS
    ax.add_patch(Rectangle((-0.5, rs - 0.5), n_along, ROAD_CELLS, fill=False,
                           ec=C["ROAD"], lw=1.5, zorder=8))
    if setback < 0:
        ax.add_patch(Rectangle((-0.5, rs - 0.5), n_along, ROAD_CELLS,
                               facecolor="none", hatch="////", ec=C["ROAD"],
                               lw=0.0, zorder=7))

    # One quantitative note per panel, same corner every time, so the panels
    # are comparable without the reader going to a table. BOTTOM right, not
    # top: the added rows, the dune strip and the road band all sit in the top
    # third of the frame, which is the part being compared, and a box up there
    # covers it. The back-barrier below is uniform and carries nothing.
    ax.text(0.985, 0.028, note, transform=ax.transAxes, fontsize=7.2,
            ha="right", va="bottom", zorder=9, linespacing=1.35,
            bbox=dict(fc="white", ec="#bbbbbb", lw=0.5, alpha=0.92, pad=2.2))

    ax.set_xlim(-0.5, n_along - 0.5)
    ax.set_ylim(nrows + DUNE_ROWS - 0.5, -0.5)
    ax.set_xticks([0, 25, 49])
    ax.set_yticks([0, 5, 10, 15, 20, 25])
    ax.set_title(panel_title(letter, title), fontsize=9)
    spines_for_image(ax)
    if reference:
        for sp in ax.spines.values():
            sp.set_linewidth(1.6)
            sp.set_edgecolor(C["ROAD"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--rows", type=int, default=26)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    D = args.domain
    apply_style()

    R = dune_topo_root("1984-start")
    v3 = np.load(R / BASE_V / "topography" / array_name("topography", D)) * 10.0
    dune = np.load(R / BASE_V / "dunes" / array_name("dune", D)) * 10.0
    aud = {int(r["domain"]): r for r in csv.DictReader(
        open(R / INS_V / "HAT_seaward_row_insert_audit.csv"))}
    n = int(aud[D]["n_rows_inserted"])
    setb = float(aud[D]["setback_model_after_m"])
    setb_meas = float(aud[D]["setback_raw_before_m"])

    ext = load_ext()
    dom = ext.load_profiles(ext.LOAD_PATH / "domain_{}.npy".format(D))
    prof = ext.masked_profiles(dom["z"])
    w = json.load(open(ext.PICKS_DIR / "HAT_dune_search_windows_{}.json"
                       .format(BASE_V)))["domain_{}".format(D)]
    _e, dl = ext.find_dunes(prof, dom["start_beach"], int(w["i0"]), int(w["i1"]))
    row0, _ = ext.interior_row0_line(prof, dl)
    n_along = v3.shape[1]

    plat = np.median(v3[1:1 + BACKDUNE_ROWS, :], axis=0)
    real = real_cells(dom, row0, n, n_along)
    dry = np.isfinite(real) & (real > 0.0)
    road_ele = road_elevation(D)
    if road_ele is None:
        raise SystemExit(
            "domain {} has no entry in {} - it carries no NC-12, so there is "
            "no road elevation to draw.".format(D, ROAD_ELEV_FILE.name))

    blocks = [
        # MATCHED BACKDUNE. The existing first n interior rows copied in front
        # of themselves, so the near-dune PROFILE is reproduced at the 1984
        # position: "the 1984 backdune looked like the present one, just
        # further seaward". Every value is a measured cell, so unlike the flat
        # fill it carries real alongshore texture - but it is counted as 0%
        # measured because those cells are measurements of the WRONG PLACE,
        # copied, not of the ground being filled.
        ("matched backdune  ·  profile moved seaward",
         v3[:n, :].copy(), 0.0),
        # KEEP DRY, FILL WET WITH THE BLOCK'S OWN MEDIAN. `--fill median`.
        # Same dry-land test as the shipped rule, but no second step: a
        # measurement is never raised, and the one invented number comes from
        # the ground being filled rather than from interior rows 1-3.
        ("measured + median fill  ·  --fill median",
         np.where(dry, real, np.median(real[dry]) if dry.any() else plat.mean()),
         100.0 * dry.mean()),
        ("raw DEM, no floor  ·  control, not a candidate",
         np.where(np.isfinite(real), real, plat[None, :]),
         100.0 * np.isfinite(real).mean()),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(7.6, 7.6),
                             sharex=True, sharey=True)
    axf = axes.ravel()
    rs_ref = int(setb_meas / 10.0)
    draw(axf[0], v3, dune, setb_meas, 0, args.rows, "a",
         "reference — {} as extracted".format(BASE_V),
         "NC-12 measured {:+.0f} m\nrow {} — SEAWARD of row 0\nwould not "
         "initialise".format(setb_meas, rs_ref), reference=True)

    rs = int(setb / 10.0)
    for k, (lab, blk, pct) in enumerate(blocks, start=1):
        full = np.vstack([blk, v3])
        under = np.median(full[rs:rs + ROAD_CELLS, :], axis=1)
        note = ("NC-12 bulldozed to {:.2f} m\nground beneath was "
                "{:.2f} / {:.2f} m\nadded block mean {:.2f} m\n"
                "{:.0f}% taken from the DEM".format(
                    road_ele, under[0], under[1], float(np.mean(blk)), pct))
        draw(axf[k], full, dune, setb, n, args.rows, "abcdef"[k], lab, note,
             road_ele=road_ele)

    for ax in axes[-1]:
        ax.set_xlabel("alongshore cell")
    for ax in axes[:, 0]:
        ax.set_ylabel("cross-shore row\n(0 = first interior cell)")

    fig_h = fig.get_figheight()
    # CAP_IN grew from 0.44 when the caption went from two lines to four
    # (the road-elevation sentence); at 0.44 it ran up into the legend row.
    CAP_IN, LEG_IN, CBAR_IN, XLAB_IN = 0.92, 0.26, 0.50, 0.46
    bottom_in = CAP_IN + LEG_IN + CBAR_IN + XLAB_IN

    cmap, norm, bounds = elevation_cmap()
    cax = fig.add_axes([0.345, (CAP_IN + LEG_IN + 0.34) / fig_h, 0.31,
                        0.095 / fig_h])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax,
                      orientation="horizontal", boundaries=bounds[1:],
                      ticks=bounds[1:-1])
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("elevation (m MHW); leftmost class is water (≤ 0)",
                 fontsize=7.4, labelpad=2)

    fig.legend(handles=[
        Line2D([0], [0], color=C["ROAD"], lw=1.5,
               label="NC-12, measured offset, at road elevation"),
        Line2D([0], [0], color=C["ROAD"], lw=6, alpha=0.35,
               label="hatched = seaward of row 0"),
        Line2D([0], [0], color=C["ADDED"], lw=1.4, ls=(0, (4, 2)),
               label="the {} added rows".format(n)),
    ], loc="lower center", bbox_to_anchor=(0.5, CAP_IN / fig_h), ncol=3,
        fontsize=7.4)

    fig.suptitle("GIS {}  ·  candidate interior fills, as the Barrier3D domain"
                 .format(D), fontsize=11, fontweight="bold", x=0.055,
                 ha="left", y=1 - 0.24 / fig_h)
    caption(fig,
            "Two dune rows over the interior domain, row 0 first. Panels "
            "(b)-(d) differ ONLY in the {} added rows; everything landward of "
            "the dashed band is identical in all three.\n"
            "NC-12 is drawn at its measured offset AND at the elevation "
            "bulldoze() gives it — one value for every road cell, from "
            "RoadElevation.csv, independent of the fill.\n"
            "The road does not move; the added rows move interior row 0 out "
            "from under it. (d) is a control, not a candidate.".format(n),
            y=0.075 / fig_h, size=7.4)
    fig.subplots_adjust(top=1 - 0.58 / fig_h, bottom=bottom_in / fig_h,
                        left=0.135, right=0.985, hspace=0.30, wspace=0.10)

    out = Path(args.out) if args.out else (
        insert_figures_dir("1984-start") / "HAT_fill_options_grid_GIS{}.png".format(D))
    fig.savefig(out)
    print("wrote {}".format(out))


if __name__ == "__main__":
    main()
