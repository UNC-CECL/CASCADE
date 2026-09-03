#!/usr/bin/env python3
r"""
HAT_plot_b3d_grid.py
==============================================================================
The Barrier3D grid itself -- the arrays CASCADE is handed -- for v1 and v2.

WHAT IS DRAWN
    Every cell the model receives for a domain, in the model's own indexing:

        dune rows      DuneDomain, 2 rows, height above berm -> m MHW
        interior rows  InteriorDomain, row 0 seaward, increasing landward
        NC-12          the 2-cell block bulldoze() writes, at
                       road_start = int(setback / dy)

    Cross-shore runs DOWN the page (row 0 at the top, behind the dune), and
    alongshore runs across -- 50 cells, 500 m.

    The beach and shoreface are NOT in these arrays. Barrier3D carries them as
    parameters, not cells, so they are named in the caption but not drawn; the
    top of the dune strip is where the model's grid begins.

WHY BOTH VERSIONS SIDE BY SIDE
    v2 prepends N interior rows behind the dune so the 1984 roadway starts where
    it historically did. In the grid that shows as the road block moving DOWN the
    page: the dune does not move relative to the array, the ground between the
    dune and the road grows.

    At GIS 85 the difference is stark. In v1 the road occupies rows 0-1 -- on
    interior row 0, which is the dune crest itself -- because the setback
    measured -10 m and was floored to 0. In v2 it sits at rows 5-6 on backdune.

USAGE
    python HAT_plot_b3d_grid.py
    python HAT_plot_b3d_grid.py --domains 85 --rows 40
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
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import array_name, dune_topo_root          # noqa: E402
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_figure_style import (apply_style, C, caption,           # noqa: E402
                              elevation_cmap, panel_title,
                              spines_for_image)

ROAD_DIR = (REPO / "data/hatteras_init/4-mgmt-forcing/road_offset"
            / "dunestart_offset/1984")
BERM_EL_M = 1.7      # BermEl, Hatteras-CASCADE-parameters.yaml
DUNE_ROWS = 2        # DuneWidth 20 m / dy 10 m; row 1 is a copy of row 0
ROAD_CELLS = 2       # road_width 20 m / dy 10 m
# Defaults, overridable per run. HARDCODING THESE IS WHAT WENT WRONG: the
# script was named _v1_v2 and defaulted its output to that name, then was
# repointed at v3/v4 -- so a default run would have written v3/v4 content
# into a file called v1_v2 and quietly replaced a correct figure.
# Default v3 -> v5 since 2026-09-03. It was v3 -> v4, which meant a bare
# run wrote HAT_b3d_grid_v3_v4.png - a figure deliberately deleted as
# superseded, so the default recreated the thing the cleanup removed.
# Pass --insert v4 for the scope comparison; v4 is still on disk.
BASE_V, INS_V = "v3", "v5"


def load_version(version, domain):
    root = dune_topo_root("1984-start") / version
    topo = np.load(root / "topography" / array_name("topography", domain)) * 10.0
    dune = np.load(root / "dunes" / array_name("dune", domain)) * 10.0
    aud = {}
    p = root / "HAT_seaward_row_insert_audit.csv"
    if p.is_file():
        aud = {int(r["domain"]): r for r in csv.DictReader(open(p))}
    return topo, dune, aud.get(domain)


def baseline_setback(domain):
    for r in csv.DictReader(open(ROAD_DIR / "RoadOffset_1984_domains.csv")):
        if int(r["domain"]) == domain and r["setback_dunestart_m"] not in ("", "nan"):
            return float(r["setback_dunestart_m"])
    return np.nan


def draw(ax, topo, dune, setback_model, n_inserted, nrows, title, ylabel=True):
    """One domain's grid: dune strip on top, interior below, NC-12 outlined."""
    cmap, norm, _ = elevation_cmap()
    n_along = topo.shape[1]
    dune_strip = np.tile(BERM_EL_M + dune[None, :n_along], (DUNE_ROWS, 1))
    grid = np.vstack([dune_strip, topo[:nrows, :]])

    ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto",
              interpolation="nearest", zorder=1)

    # dune strip, ruled off from the interior
    ax.axhline(DUNE_ROWS - 0.5, color="#333333", lw=1.1, zorder=4)
    ax.text(n_along - 0.8, DUNE_ROWS / 2 - 0.5, "dune rows", fontsize=7.5,
            va="center", ha="right", zorder=6,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))

    if n_inserted:
        ax.add_patch(Rectangle((-0.5, DUNE_ROWS - 0.5), n_along, n_inserted,
                               fill=False, ec=C["ADDED"], lw=1.6, ls=(0, (4, 2)),
                               zorder=5))
        ax.text(n_along * 0.5, DUNE_ROWS + n_inserted / 2 - 0.5,
                "+{} rows".format(n_inserted), fontsize=7.5, ha="center",
                va="center", color="#7a5200", zorder=6,
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))

    rs = int(setback_model / 10.0) + DUNE_ROWS
    ax.add_patch(Rectangle((-0.5, rs - 0.5), n_along, ROAD_CELLS, fill=False,
                           ec=C["ROAD"], lw=1.6, zorder=6))
    ax.text(0.8, rs + ROAD_CELLS / 2 - 0.5, "NC-12", fontsize=7.5, ha="left",
            va="center", color=C["ROAD"], zorder=7,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))

    ax.set_xlim(-0.5, n_along - 0.5)
    ax.set_ylim(nrows + DUNE_ROWS - 0.5, -0.5)
    ax.set_xticks([0, 10, 20, 30, 40, 49])
    ax.set_xlabel("alongshore cell")
    if ylabel:
        ax.set_ylabel("cross-shore row\n(0 = first interior cell)")
    ax.set_title(title)
    spines_for_image(ax)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domains", default="84,85,86")
    ap.add_argument("--rows", type=int, default=32)
    global BASE_V, INS_V
    ap.add_argument("--base", default=BASE_V, help="extraction version")
    ap.add_argument("--insert", default=INS_V, help="version with rows added")
    ap.add_argument("--out", default=None,
                    help="default is named for the two versions actually drawn")
    args = ap.parse_args()
    BASE_V, INS_V = args.base, args.insert
    domains = [int(x) for x in args.domains.split(",")]
    apply_style()

    nd = len(domains)
    fig, axes = plt.subplots(nd, 2, squeeze=False,
                             figsize=(9.2, 2.45 * nd + 2.5))
    letters = "abcdefgh"
    for r, D in enumerate(domains):
        t1, d1, _ = load_version(BASE_V, D)
        t2, d2, aud = load_version(INS_V, D)
        n = int(aud["n_rows_inserted"]) if aud else 0
        raw1 = baseline_setback(D)
        sb1 = max(raw1, 0.0)                       # the model floors it
        sb2 = float(aud["setback_model_after_m"]) if aud else sb1

        left = panel_title(letters[2 * r], "GIS {} · " + BASE_V + " · setback {:+.0f} m{}"
                           .format(D, raw1,
                                   " (floored to 0)" if raw1 <= 0 else ""))
        right = panel_title(letters[2 * r + 1],
                            ("GIS {} · " + INS_V + " · setback +{:.0f} m").format(D, sb2))
        draw(axes[r][0], t1, d1, sb1, 0, args.rows, left)
        draw(axes[r][1], t2, d2, sb2, n, args.rows, right, ylabel=False)

    # LAYOUT IS COMPUTED IN INCHES, then converted to figure fractions.
    #
    # Fractional margins are only correct for the figure height they were tuned
    # on. This figure's height scales with the number of domains, so a bottom
    # margin of "0.195" reserved 1.9 in at three domains and 0.87 in at one --
    # and the one-domain version put the caption through the panels. Inches are
    # the quantity the furniture actually needs.
    # Stacked from the bottom up: caption, legend, colourbar label, colourbar
    # bar, then the axes' own x-label, then the panels. Each gets its height in
    # inches and the fractions fall out of that.
    fig_h = fig.get_figheight()
    CAP_IN = 0.62      # three wrapped lines at 7.8 pt
    LEG_IN = 0.26
    CBAR_IN = 0.50     # bar plus its label, which sits BELOW a horizontal bar
    XLABEL_IN = 0.42   # the panels' own "alongshore cell"
    TITLE_IN = 0.80    # suptitle plus the panel titles
    bottom_in = CAP_IN + LEG_IN + CBAR_IN + XLABEL_IN

    cmap, norm, bounds = elevation_cmap()
    cbar_y = (CAP_IN + LEG_IN + 0.34) / fig_h
    cax = fig.add_axes([0.36, cbar_y, 0.30, 0.10 / fig_h])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax,
                      orientation="horizontal", boundaries=bounds[1:],
                      ticks=bounds[1:-1])
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("elevation (m MHW); leftmost class is water (≤ 0)",
                 fontsize=7.6, labelpad=2)

    fig.legend(handles=[
        Line2D([0], [0], color=C["ROAD"], lw=1.6,
               label="NC-12, the 2-cell block bulldoze() writes"),
        Line2D([0], [0], color=C["ADDED"], lw=1.6, ls=(0, (4, 2)),
               label="interior rows added in " + INS_V),
    ], loc="lower center", bbox_to_anchor=(0.5, CAP_IN / fig_h), ncol=2,
        fontsize=7.6)

    fig.suptitle("The Barrier3D grid CASCADE receives, before and after the "
                 "row insert", fontsize=11, fontweight="bold", x=0.055,
                 ha="left", y=1 - 0.26 / fig_h)

    caption(fig,
            "Cells are the arrays the model is handed: two dune rows (crest "
            "height above berm, as m MHW) over the interior domain, row 0 "
            "first. Beach and\n"
            "shoreface are Barrier3D parameters, not cells, and are not drawn. "
            "v4 adds N rows behind the dune, N = (1984 − 1997 dune-line "
            "difference) / 10 m,\n"
            "so interior row 0 sits at the 1984 dune position. The added rows are "
            "the ONLY change — no existing cell is modified — so island\n"
            "width grows by exactly N at the 8 domains that got rows "
            "(GIS 85, 353 → 413 m).",
            y=0.09 / fig_h)

    fig.subplots_adjust(top=1 - TITLE_IN / fig_h, bottom=bottom_in / fig_h,
                        hspace=0.50, wspace=0.09, left=0.105, right=0.985)
    # Named for what it DRAWS, so it cannot silently replace another pair's
    # figure.
    out = Path(args.out) if args.out else (
        insert_figures_dir("1984-start")
        / "HAT_b3d_grid_{}_{}{}.png".format(
            BASE_V, INS_V,
            "_GIS{}".format(domains[0]) if len(domains) == 1 else ""))
    fig.savefig(out)
    print("wrote {}".format(out))


if __name__ == "__main__":
    main()
