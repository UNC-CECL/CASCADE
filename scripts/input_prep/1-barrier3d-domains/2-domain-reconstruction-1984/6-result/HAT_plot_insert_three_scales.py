#!/usr/bin/env python3
r"""
HAT_plot_insert_three_scales.py
==============================================================================
The seaward row insert at three scales: GIS 85, the two relocation blocks, and
the whole island.

WHAT THE INSERT IS
    The 1984-start DEM is a 1996 ALACE beach on a 2009 backdune, so its dune has
    already migrated landward past the 1984 NC-12 alignment. At GIS 85 that puts
    the 1984 roadbed SEAWARD of interior row 0 -- setback -10 m, floored to 0,
    and a road that relocates in model year 1 by construction.

    The fix measures how far the dune line moved between 1984 and 1997 (both
    digitized from imagery, same feature) and prepends that many interior rows
    behind the dune, so row 0 sits at the 1984 dune position and NC-12 lands its
    true distance behind it.

WHY THREE SCALES
    Row 1 shows the mechanism on the domain the work was for.
    Row 2 shows every domain with a documented historical relocation, which is
        where the correction is actually applied.
    Row 3 shows the measured retreat for all 90 domains -- context for whether
        GIS 85 is exceptional or typical. NOTE it is MEASURED island-wide but
        APPLIED only to the ten block domains; the row 3 bars outside the shaded
        blocks are what an island-wide version WOULD insert, not what it did.

USAGE
    python HAT_plot_insert_three_scales.py
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
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import array_name, dune_topo_root  # noqa: E402
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_topo_version import require_version  # noqa: E402
from hat_topo_version import duneline_shift_dir  # noqa: E402
from hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              DOMAIN_AXIS_LABEL, caption, figsize,
                              open_frame, save, town_bands, _title)

# Resolved through hat_topo_version.duneline_shift_dir - ONE definition
# of a path that eight scripts used to build by hand. Moved under
# 2-domain-reconstruction-1984/ on 2026-09-03.
SHIFT = duneline_shift_dir("1984-start")
ROAD = (REPO / "data/hatteras_init/4-mgmt-forcing/road_offset/dunestart_offset/1984")
BASE_VERSION = "v2"   # the re-picked extraction (was "v3")
# DELETED 2026-09-07 with every layer (only unmodified topography is kept);
# the literal is kept as the name of what this drew. require_version() in
# main() says so before any array is opened.
VERSION = "v4"   # island scope, measured + floor (was "v5"). v2 + rows. Was "v4" (block
                      # scope) until 2026-09-03. N is identical at the ten
                      # BLOCK domains, so panels (a)-(c) are unchanged - but
                      # panel (d) is NOT: it goes from 8 red bars to 38,
                      # because v5 applies rows wherever the measurement
                      # selects them. Every label below is derived from
                      # VERSION and from the data, never written literally,
                      # so switching the version cannot leave a stale caption.
BLOCK_A, BLOCK_B = list(range(9, 15)), list(range(84, 88))
# Shared vocabulary: BASE = the unmodified input, ACCENT = the change
# under test, ADDED = fabricated ground. Same meanings in every figure.
GREY, RED = C["BASE"], C["ACCENT"]
LAND = 0.0


def read(path, key="shift_m_median", idx="domain"):
    return {int(r[idx]): float(r[key]) for r in csv.DictReader(open(path))
            if r[key] not in ("", "nan")}


def main() -> None:
    ap = argparse.ArgumentParser()
    # Named for the versions drawn. The old default was a fixed
    # "HAT_insert_three_scales.png", so repointing BASE_VERSION/VERSION from
    # v1/v2 to v3/v4 would have overwritten the v1/v2 figure in place.
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    require_version("1984-start", VERSION, "VERSION, the layer drawn")
    apply_style()

    date = read(SHIFT / "duneline_retreat_1984_1997.csv")
    before = read(ROAD / "RoadOffset_1984_domains.csv", "setback_dunestart_m")
    audit = {int(r["domain"]): r for r in csv.DictReader(
        open(dune_topo_root("1984-start") / VERSION
             / "HAT_seaward_row_insert_audit.csv"))}
    after = {d: float(r["setback_raw_after_m"]) for d, r in audit.items()}
    nrows = {d: int(r["n_rows_inserted"]) for d, r in audit.items()}

    fig = plt.figure(figsize=figsize("double", aspect=1.12),
                     constrained_layout=True)
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.05, 0.92, 1.0])

    # ---------------------------------------------------- ROW 1: GIS 85 alone
    #
    # The monospace block of GIS 85 numbers that used to sit at gs[0, 1] came
    # off the canvas 2026-09-10: a statistics table is caption text under the
    # house rules, and it is now the last third of caption() below. Panel (a)
    # takes the full width in its place.
    ax = fig.add_subplot(gs[0, :])
    v1 = np.load(dune_topo_root("1984-start") / BASE_VERSION / "topography"
                 / array_name("topography", 85))
    ins = np.load(dune_topo_root("1984-start") / VERSION / "topography"
                  / array_name("topography", 85))
    n = nrows[85]
    x1 = np.arange(v1.shape[0]) * 10.0
    x2 = (np.arange(ins.shape[0]) - n) * 10.0        # common frame: v1's row 0
    m1, m2 = np.median(v1, axis=1) * 10, np.median(ins, axis=1) * 10
    ax.axvspan(-n * 10 - 5, 0, color=C["ADDED_FILL"], zorder=0)
    k1, k2 = x1 <= 420, x2 <= 420
    ax.plot(x1[k1], m1[k1], "-", color=GREY, lw=3.2, label="as extracted")
    ax.plot(x2[k2], m2[k2], "-", color=RED, lw=1.4,
            label="rows added, nothing else changed")
    ax.add_patch(plt.Rectangle((before[85], -0.85), 20, 0.5, color=C["ROAD"],
                               zorder=6))
    ax.annotate("NC-12, which does not move",
                xy=(before[85] + 10, -0.35), xytext=(70, -0.78), fontsize=7.5,
                color=INK, va="center",
                arrowprops=dict(arrowstyle="->", lw=0.8, color=INK))
    ax.axvline(0, color=INK, lw=1.0)
    ax.text(6, 5.4, "row 0 as extracted", fontsize=7.5, rotation=90, va="top",
            color=INK_MUTED)
    ax.axvline(-n * 10, color=RED, lw=1.0, ls=(0, (4, 2)))
    ax.text(-n * 10 + 6, 5.4, "row 0 after the insert", fontsize=7.5,
            rotation=90, va="top", color=RED)
    ax.axhline(0, color=C["WATER"], lw=0.9, ls=":")
    ax.set_xlim(-n * 10 - 25, 420)
    ax.set_ylim(-1.0, 5.6)
    ax.set_xlabel("metres landward of interior row 0 as extracted")
    ax.set_ylabel("median elevation (m MHW)")
    ax.grid(axis="y")
    open_frame(ax)
    _title(ax, 0, "one domain: GIS 85")

    # ------------------------------------- ROW 2: the two relocation blocks
    for j, (blk, name, yr) in enumerate(
            ((BLOCK_A, "GIS 9–14, inter-village", 1999),
             (BLOCK_B, "GIS 84–87, Pea Island", 1989))):
        axb = fig.add_subplot(gs[1, j])
        d = [g for g in blk if g in after]
        xi = np.arange(len(d))
        axb.bar(xi - 0.2, [before[g] for g in d], 0.4, color=GREY,
                label="as extracted")
        axb.bar(xi + 0.2, [after[g] for g in d], 0.4, color=RED,
                label="after the rows are added")
        for k, g in enumerate(d):
            if nrows[g]:
                axb.text(k + 0.2, after[g] + 2, "+{}".format(nrows[g]),
                         ha="center", fontsize=7, color=RED)
            if before[g] <= 0:
                axb.text(k - 0.2, 2, "floored", ha="center", fontsize=7,
                         rotation=90, color=INK_MUTED)
        axb.axhline(0, color=INK, lw=0.8)
        axb.set_xticks(xi)
        axb.set_xticklabels(d)
        axb.set_xlabel(DOMAIN_AXIS_LABEL)
        axb.set_ylabel("road setback (m)")
        axb.grid(axis="y")
        open_frame(axb)
        _title(axb, 1 + j, name)

    # ------------------------------------------------ ROW 3: the whole island
    axi = fig.add_subplot(gs[2, :])
    gis = np.array(sorted(date))
    val = np.array([date[g] for g in gis])
    applied = np.array([g in after and nrows[g] > 0 for g in gis])
    axi.bar(gis[~applied], val[~applied], 0.85, color=C["BASE_FILL"],
            label="measured, no rows added")
    axi.bar(gis[applied], val[applied], 0.85, color=RED,
            label="rows added")
    for lo, hi in ((9, 14), (84, 87)):
        axi.axvspan(lo - .5, hi + .5, color=C["ADDED_FILL"], alpha=.55,
                    zorder=0)
    axi.axhline(0, color=INK, lw=0.8)
    axi.axhline(np.median(val), color=C["REF"], ls=(0, (4, 2)), lw=1.1,
                label="island median")
    axi.annotate("GIS 85", xy=(85, date[85]), xytext=(76, date[85] + 20),
                 fontsize=8, color=INK,
                 arrowprops=dict(arrowstyle="->", lw=0.8, color=INK))
    axi.set_xlim(0, 91)
    axi.set_xlabel(DOMAIN_AXIS_LABEL)
    axi.set_ylabel("measured 1984→1997\ndune-line retreat (m)")
    axi.grid(axis="y")
    open_frame(axi)
    town_bands(axi, strip=0.07)
    _title(axi, 3, "the whole island")

    fig.legend(handles=[
        Patch(facecolor=GREY, label="as extracted"),
        Patch(facecolor=RED, label="rows added behind the dune"),
        Patch(facecolor=C["BASE_FILL"], label="measured, no rows added"),
        Line2D([0], [0], color=C["REF"], ls=(0, (4, 2)), lw=1.1,
               label="island median"),
        Patch(facecolor=C["ADDED_FILL"], label="relocation block"),
    ], loc="outside lower center", ncol=5, frameon=False, fontsize=7.5)

    caption(fig,
            "Interior rows inserted behind the dune so that the 1984 roadway "
            "starts where it historically did. N = (1984 − 1997 dune-line "
            "difference) / 10 m, both lines digitized from imagery to the same "
            "feature, so the definitional offset between a digitized line and "
            "the model’s interior row 0 cancels and what remains is date. "
            "(a) is the domain the work was for, in profile: the road does not "
            "move, row 0 moves seaward by {n} cells. (b) and (c) are the two "
            "blocks with a documented historical relocation, which is where "
            "the correction is applied; “floored” marks a domain whose "
            "measured setback is negative and which the model therefore starts "
            "with the road on interior row 0, and “+N” is the number of rows "
            "added. (d) is the measurement for all 90 domains, context for "
            "whether GIS 85 is exceptional: domain 1 is at Cape Point in the "
            "south and domain 90 at Pea Island in the north, the light strips "
            "along the top name the villages, and the two shaded bands are the "
            "relocation blocks of (b) and (c). Rows were added at {na} of the "
            "90; the pale bars are the domains where the measurement rounds to "
            "N = 0, so no 1984 land is missing and nothing is inserted. The "
            "island median retreat is {med:+.1f} m (IQR {q1:+.1f} to "
            "{q3:+.1f} m) against GIS 85’s {g85:+.1f} m, its "
            "{pc:.0f}th percentile. At GIS 85 that is {n} cells, so the road "
            "setback goes from {b85:+.0f} m — floored to 0, the road on the "
            "dune crest itself — to {a85:+.0f} m, and the Barrier3D island "
            "grows by {w} m, because nothing is retired to pay for the added "
            "rows. Recorded with the figure when it was built: as extracted "
            "the road sat on 2.94 m and relocated in model year 1985; with the "
            "rows added it sits on real DEM cells, 47% of the block, and "
            "relocates in 1995. The observed relocation was 1989."
            .format(n=n, na=int(applied.sum()), med=np.median(val),
                    q1=np.percentile(val, 25), q3=np.percentile(val, 75),
                    g85=date[85], pc=100 * (val < date[85]).mean(),
                    b85=before[85], a85=after[85], w=int(10 * n)))

    out = Path(args.out) if args.out else (
        # NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
        # guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
        insert_figures_dir("1984-start", "6-result")
        / "HAT_insert_three_scales_{}_{}.png".format(BASE_VERSION, VERSION))
    written = save(fig, out)
    print("wrote {}".format(written[0]))
    print("\nisland-wide 1984->1997 dune retreat: median {:+.1f} m, "
          "IQR {:+.1f} to {:+.1f}, GIS 85 {:+.1f} m ({:.0f}th percentile)".format(
              np.median(val), np.percentile(val, 25), np.percentile(val, 75),
              date[85], 100 * (val < date[85]).mean()))


if __name__ == "__main__":
    main()
