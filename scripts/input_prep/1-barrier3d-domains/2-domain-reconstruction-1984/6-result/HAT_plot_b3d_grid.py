#!/usr/bin/env python3
r"""
HAT_plot_b3d_grid.py
==============================================================================
The Barrier3D grid itself -- the arrays CASCADE is handed -- for the
extraction and the layer built from it, side by side (--base/--insert).

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
from site_layer.hat_topo_version import array_name, dune_topo_root          # noqa: E402
from site_layer.hat_topo_version import insert_figures_dir  # noqa: E402
from site_layer.hat_topo_version import require_version  # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, INK, caption,      # noqa: E402
                              elevation_cmap, figsize, save,
                              spines_for_image, _title)

from site_layer import hat_topo_version as _tv  # noqa: E402
ROAD_DIR = _tv.road_setback_dir(1984)
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
# DELETED 2026-09-07 with every layer (only unmodified topography is kept);
# the literal is kept as the name of what this drew. require_version() in
# main() says so before any array is opened.
BASE_V, INS_V = "v2", "v4"   # base was "v3" and the layer "v5" until the 2026-09-04 renumber


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


def draw(ax, topo, dune, setback_model, n_inserted, nrows, idx, title,
         ylabel=True, xlabel=True):
    """One domain's grid: dune strip on top, interior below, NC-12 outlined.

    The per-panel "+N rows" count came off the canvas 2026-09-10 with the
    title's setback figures: the house rule puts statistics in the caption,
    and main() writes one caption clause per domain from the same numbers."""
    cmap, norm, _ = elevation_cmap()
    n_along = topo.shape[1]
    dune_strip = np.tile(BERM_EL_M + dune[None, :n_along], (DUNE_ROWS, 1))
    grid = np.vstack([dune_strip, topo[:nrows, :]])

    ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto",
              interpolation="nearest", zorder=1)

    # dune strip, ruled off from the interior
    ax.axhline(DUNE_ROWS - 0.5, color=INK, lw=1.0, zorder=4)
    ax.text(n_along - 0.8, DUNE_ROWS / 2 - 0.5, "dune rows", fontsize=7,
            va="center", ha="right", zorder=6, color=INK,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))

    if n_inserted:
        ax.add_patch(Rectangle((-0.5, DUNE_ROWS - 0.5), n_along, n_inserted,
                               fill=False, ec=C["ADDED"], lw=1.4, ls=(0, (4, 2)),
                               zorder=5))

    rs = int(setback_model / 10.0) + DUNE_ROWS
    ax.add_patch(Rectangle((-0.5, rs - 0.5), n_along, ROAD_CELLS, fill=False,
                           ec=C["ROAD"], lw=1.4, zorder=6))
    ax.text(0.8, rs + ROAD_CELLS / 2 - 0.5, "NC-12", fontsize=7, ha="left",
            va="center", color=C["ROAD"], zorder=7,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))

    ax.set_xlim(-0.5, n_along - 0.5)
    ax.set_ylim(nrows + DUNE_ROWS - 0.5, -0.5)
    ax.set_xticks([0, 10, 20, 30, 40, 49])
    if xlabel:
        ax.set_xlabel("alongshore cell")
    if ylabel:
        ax.set_ylabel("cross-shore row\n(0 = first interior cell)")
    spines_for_image(ax)
    _title(ax, idx, title)


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
    require_version("1984-start", INS_V, "the version with rows added; pass --insert")
    domains = [int(x) for x in args.domains.split(",")]
    apply_style()

    nd = len(domains)
    fig, axes = plt.subplots(nd, 2, squeeze=False,
                             figsize=figsize("double",
                                             height=min(2.1 * nd + 0.9, 9.4)),
                             constrained_layout=True)
    notes = []
    for r, D in enumerate(domains):
        t1, d1, _ = load_version(BASE_V, D)
        t2, d2, aud = load_version(INS_V, D)
        n = int(aud["n_rows_inserted"]) if aud else 0
        raw1 = baseline_setback(D)
        sb1 = max(raw1, 0.0)                       # the model floors it
        sb2 = float(aud["setback_model_after_m"]) if aud else sb1

        draw(axes[r][0], t1, d1, sb1, 0, args.rows, 2 * r,
             "GIS {} as extracted".format(D), xlabel=(r == nd - 1))
        draw(axes[r][1], t2, d2, sb2, n, args.rows, 2 * r + 1,
             "GIS {} with rows added".format(D), ylabel=False,
             xlabel=(r == nd - 1))
        notes.append(
            "At GIS {D} the setback measures {raw:+.0f} m{fl}; {n} rows are "
            "added, and the road ends up {aft:+.0f} m behind interior row 0."
            .format(D=D, raw=raw1,
                    fl=", which the model floors to 0, putting NC-12 on interior "
                       "row 0 itself — the dune crest" if raw1 <= 0 else "",
                    n=n, aft=sb2))

    cmap, norm, bounds = elevation_cmap()
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm),
                      ax=axes.ravel().tolist(), orientation="horizontal",
                      location="bottom", boundaries=bounds[1:],
                      ticks=bounds[1:-1], shrink=0.42, aspect=34, pad=0.015)
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("elevation (m MHW); leftmost class is water (≤ 0)",
                 fontsize=7.6, labelpad=2)

    fig.legend(handles=[
        Line2D([0], [0], color=C["ROAD"], lw=1.4,
               label="NC-12, the 2-cell block bulldoze() writes"),
        Line2D([0], [0], color=C["ADDED"], lw=1.4, ls=(0, (4, 2)),
               label="interior rows added behind the dune"),
    ], loc="outside lower center", ncol=2, frameon=False, fontsize=7.6)

    caption(fig,
            "The arrays CASCADE is handed, in the model’s own indexing: two "
            "dune rows (crest height above the berm, as m above mean high "
            "water) over the interior domain, row 0 first, cross-shore down "
            "the page and alongshore across (50 cells, 500 m). Beach and "
            "shoreface are Barrier3D parameters rather than cells and are not "
            "drawn, so the top of the dune strip is where the grid begins. "
            "The left column of each row is the domain as extracted; the right "
            "column is the same domain with N interior rows added behind the "
            "dune, N = (1984 − 1997 dune-line difference) / 10 m, so that "
            "interior row 0 sits at the 1984 dune position. The added rows are "
            "the ONLY change — no existing cell is modified — so the island "
            "grows by exactly N and the road block moves DOWN the page: the "
            "dune does not move relative to the array, the ground between the "
            "dune and the road grows. {notes}"
            .format(notes=" ".join(notes)))

    # Named for what it DRAWS, so it cannot silently replace another pair's
    # figure.
    out = Path(args.out) if args.out else (
        # NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
        # guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
        insert_figures_dir("1984-start", "6-result")
        / "HAT_b3d_grid_{}_{}{}.png".format(
            BASE_V, INS_V,
            "_GIS{}".format(domains[0]) if len(domains) == 1 else ""))
    written = save(fig, out, vector=False)
    print("wrote {}".format(written[0]))


if __name__ == "__main__":
    main()
