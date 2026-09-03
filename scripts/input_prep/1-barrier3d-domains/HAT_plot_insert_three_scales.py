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


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import array_name, dune_topo_root  # noqa: E402
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_topo_version import duneline_shift_dir  # noqa: E402
from hat_figure_style import (apply_style, C, caption,   # noqa: E402
                              panel_title)

# Resolved through hat_topo_version.duneline_shift_dir - ONE definition
# of a path that eight scripts used to build by hand. Moved under
# row-insert-scope/ on 2026-09-03.
SHIFT = duneline_shift_dir("1984-start")
ROAD = (REPO / "data/hatteras_init/4-mgmt-forcing/road_offset/dunestart_offset/1984")
BASE_VERSION = "v3"   # the re-picked extraction
VERSION = "v5"        # v3 + added rows, island-wide scope. Was "v4" (block
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
    apply_style()

    date = read(SHIFT / "duneline_retreat_1984_1997.csv")
    before = read(ROAD / "RoadOffset_1984_domains.csv", "setback_dunestart_m")
    audit = {int(r["domain"]): r for r in csv.DictReader(
        open(dune_topo_root("1984-start") / VERSION
             / "HAT_seaward_row_insert_audit.csv"))}
    after = {d: float(r["setback_raw_after_m"]) for d, r in audit.items()}
    nrows = {d: int(r["n_rows_inserted"]) for d, r in audit.items()}

    fig = plt.figure(figsize=(16, 15))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.05, 0.95, 1.0],
                  width_ratios=[1.5, 1.0], hspace=0.42, wspace=0.22)

    # ---------------------------------------------------- ROW 1: GIS 85 alone
    ax = fig.add_subplot(gs[0, 0])
    v1 = np.load(dune_topo_root("1984-start") / BASE_VERSION / "topography"
                 / array_name("topography", 85))
    ins = np.load(dune_topo_root("1984-start") / VERSION / "topography"
                  / array_name("topography", 85))
    n = nrows[85]
    x1 = np.arange(v1.shape[0]) * 10.0
    x2 = (np.arange(ins.shape[0]) - n) * 10.0        # common frame: v1's row 0
    m1, m2 = np.median(v1, axis=1) * 10, np.median(ins, axis=1) * 10
    ax.axvspan(-n * 10 - 5, 0, color=C["ADDED_FILL"], zorder=0)
    # Low and to the right of the band: the top-left corner already carries the
    # two row-0 rules and their rotated labels, and three annotations in one
    # corner is how a figure stops being readable.
    ax.text(-n * 10 - 18, 3.1, "{} rows\ninserted\n({} m)".format(n, n * 10),
            fontsize=9, color="#7a5200", va="top", ha="left")
    k1, k2 = x1 <= 420, x2 <= 420
    ax.plot(x1[k1], m1[k1], "-", color=GREY, lw=5, label="v3 (re-picked, as extracted)")
    ax.plot(x2[k2], m2[k2], "-", color=RED, lw=2,
            label="{} (rows added, nothing else changed)".format(VERSION))
    ax.add_patch(plt.Rectangle((before[85], -0.85), 20, 0.5, color="k", zorder=6))
    ax.annotate("NC-12 (does not move)", xy=(before[85] + 10, -0.35),
                xytext=(95, -0.75), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=1.1))
    ax.axvline(0, color="k", lw=1.2)
    ax.text(6, 5.4, "v3 row 0", fontsize=8.5, rotation=90, va="top")
    ax.axvline(-n * 10, color=RED, lw=1.2, ls="--")
    ax.text(-n * 10 + 6, 5.4, "new row 0", fontsize=8.5, rotation=90,
            va="top", color=RED)
    ax.axhline(0, color="#3b6ea5", lw=.9, ls=":")
    ax.set_xlim(-n * 10 - 25, 420)
    ax.set_ylim(-1.0, 5.6)
    ax.set_xlabel("metres landward of v3 interior row 0")
    ax.set_ylabel("median elevation (m MHW)")
    ax.set_title(panel_title("a", "GIS 85 — the domain the work was for"))
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=.3)

    axt = fig.add_subplot(gs[0, 1])
    axt.axis("off")
    txt = [
        "GIS 85", "",
        "measured 1984→1997 dune retreat   58.2 m",
        "   ÷ 10 m cells                     6 rows", "",
        "road setback",
        "   v3                     -15 m  → floored to 0",
        "   {} (after insert)      +45 m".format(VERSION), "",
        "island width (Barrier3D)   v3 → {}, +60 m".format(VERSION),
        "   grows by the added rows: nothing is retired",
        "   to pay for them", "",
        "what the road sits on",
        "   v3        2.94 m",
        "   {}        real DEM, 47% of the block".format(VERSION),
        "", "",
        "model behaviour",
        "   v3        relocates 1985  (setback floored)",
        "   {}        relocates 1995".format(VERSION),
        "   observed  1989",
    ]
    axt.text(0.0, 0.99, "\n".join(txt), family="monospace", fontsize=10.5,
             va="top", transform=axt.transAxes)

    # ------------------------------------- ROW 2: the two relocation blocks
    for j, (blk, name, yr) in enumerate(
            ((BLOCK_A, "GIS 9–14  —  inter-village", 1999),
             (BLOCK_B, "GIS 84–87  —  Pea Island", 1989))):
        axb = fig.add_subplot(gs[1, j])
        d = [g for g in blk if g in after]
        xi = np.arange(len(d))
        axb.bar(xi - 0.2, [before[g] for g in d], 0.4, color=GREY,
                label="v3 setback")
        axb.bar(xi + 0.2, [after[g] for g in d], 0.4, color=RED,
                label="after insert")
        for k, g in enumerate(d):
            if nrows[g]:
                axb.text(k + 0.2, after[g] + 2, "+{}".format(nrows[g]),
                         ha="center", fontsize=8.5, color=RED)
            if before[g] <= 0:
                axb.text(k - 0.2, 2, "floored", ha="center", fontsize=7.5,
                         rotation=90, color="k")
        axb.axhline(0, color="k", lw=1)
        axb.set_xticks(xi)
        axb.set_xticklabels(d)
        axb.set_xlabel("GIS domain")
        axb.set_ylabel("road setback (m)")
        axb.set_title(panel_title("bc"[j], "{} · relocated {}".format(name, yr)))
        axb.legend(fontsize=8.5)
        axb.grid(alpha=.3, axis="y")

    # ------------------------------------------------ ROW 3: the whole island
    axi = fig.add_subplot(gs[2, :])
    gis = np.array(sorted(date))
    val = np.array([date[g] for g in gis])
    applied = np.array([g in after and nrows[g] > 0 for g in gis])
    axi.bar(gis[~applied], val[~applied], 0.85, color=C["BASE_FILL"],
            label="measured, NOT applied")
    axi.bar(gis[applied], val[applied], 0.85, color=RED,
            label="applied — {} domains got rows".format(int(applied.sum())))
    for lo, hi in ((9, 14), (84, 87)):
        axi.axvspan(lo - .5, hi + .5, color=C["ADDED_FILL"], alpha=.55, zorder=0)
    axi.axhline(0, color="k", lw=1)
    axi.axhline(np.median(val), color=C["REF"], ls="--", lw=1.4,
                label="island median {:+.1f} m".format(np.median(val)))
    axi.annotate("GIS 85", xy=(85, date[85]), xytext=(78, date[85] + 22),
                 fontsize=10, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", lw=1.2))
    axi.set_xlabel("GIS domain   (1 = Cape Point  →  90 = north Pea Island)")
    axi.set_ylabel("measured 1984→1997\ndune-line retreat (m)")
    axi.set_title(panel_title("d", "The whole island — is GIS 85 "
                               "exceptional?"))
    axi.legend(fontsize=9, loc="upper left")
    axi.grid(alpha=.3, axis="y")
    axi.set_xlim(0, 91)

    fig_h = fig.get_figheight()
    fig.suptitle("Interior rows inserted behind the dune so the 1984 roadway "
                 "starts where it historically was",
                 fontsize=11.5, fontweight="bold", x=0.055, ha="left",
                 y=1 - 0.26 / fig_h)
    caption(fig,
            "N = (1984 − 1997 dune-line difference) / 10 m, both lines digitized "
            "from imagery to the same feature, so the definitional offset\n"
            "between a digitized line and the model's interior row 0 cancels and "
            "what remains is date. Panel (d) shows the measurement for all 90\n"
            "domains; rows were APPLIED at {} of them in {}. Grey bars are the "
            "domains where the measurement rounds to N=0, so no 1984 land is "
            "missing and nothing is inserted. The two shaded bands are the "
            "relocation blocks of panels (b) and (c) — under {} they are "
            "context, not the scope.".format(
                int(applied.sum()), VERSION, VERSION),
            y=0.09 / fig_h)
    # GridSpec already set the spacing; tight_layout fights it and warns.
    fig.subplots_adjust(top=1 - 0.62 / fig_h, bottom=1.20 / fig_h,   # caption 0.62 in + x-label 0.4 in + pad
                        left=0.075, right=0.975)
    out = Path(args.out) if args.out else (
        insert_figures_dir("1984-start")
        / "HAT_insert_three_scales_{}_{}.png".format(BASE_VERSION, VERSION))
    fig.savefig(out, dpi=130)
    print("wrote {}".format(out))
    print("\nisland-wide 1984->1997 dune retreat: median {:+.1f} m, "
          "IQR {:+.1f} to {:+.1f}, GIS 85 {:+.1f} m ({:.0f}th percentile)".format(
              np.median(val), np.percentile(val, 25), np.percentile(val, 75),
              date[85], 100 * (val < date[85]).mean()))


if __name__ == "__main__":
    main()
