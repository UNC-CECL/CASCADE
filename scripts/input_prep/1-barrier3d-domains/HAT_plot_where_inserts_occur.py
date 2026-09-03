#!/usr/bin/env python3
r"""
HAT_plot_where_inserts_occur.py
==============================================================================
Where the insert changes the domain and where it does not -- zoomed on the two
relocation blocks, with island-wide context.

THREE REASONS A DOMAIN IS UNCHANGED, and they are not the same thing
    1. OUT OF SCOPE. 80 domains were never targeted. Many of them have a
       measured date term well over a cell, so they are unchanged by a decision
       about scope, not by a measurement.
    2. IN SCOPE, N ROUNDS TO 0. GIS 9 (-4.9 m) and GIS 14 (+2.7 m) were
       targeted and measured; the measurement says no cell is needed. That is
       the method declining to act, which is different from not being asked.
    3. MEASURED PROGRADATION. GIS 9's date term is NEGATIVE -- the 1997 line is
       seaward of the 1984 one. Nothing is inserted, and nothing is removed
       either: N is floored at 0. Whether progradation should ever REMOVE rows
       is an open question, not a settled one.

THE INTERVAL CAVEAT, DRAWN
    The dune lines are 1984 and 1997 -- 13 years. The DEM surface at row 0 and
    at every inserted cell is 1996 ALACE (verified: 100% of cells, all ten
    domains, no 2009 or 2014). So the interval wanted is 12 years and the one
    measured is 13. The 1996 and 1997 aerials are 1996-10-14 and 1997-10-12,
    363 days apart, so the overstatement is very close to one year of retreat.
    The x12/13 bar shows which domains' N that moves.

USAGE
    python HAT_plot_where_inserts_occur.py
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
from hat_topo_version import dune_topo_root                       # noqa: E402
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_topo_version import duneline_shift_dir  # noqa: E402
from hat_figure_style import (apply_style, C, caption,            # noqa: E402
                              panel_title)

# Resolved through hat_topo_version.duneline_shift_dir - ONE definition
# of a path that eight scripts used to build by hand. Moved under
# row-insert-scope/ on 2026-09-03.
S = duneline_shift_dir("1984-start")
BLOCK_A, BLOCK_B = list(range(9, 15)), list(range(84, 88))
SCOPE = BLOCK_A + BLOCK_B
YEARS_MEASURED, YEARS_WANTED = 13.0, 12.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="v5", help="the version to describe")
    ap.add_argument("--compare", default="v4",
                    help="the earlier, block-scoped version to contrast with")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    apply_style()

    date = {int(r["domain"]): float(r["shift_m_median"])
            for r in csv.DictReader(open(S / "duneline_retreat_1984_1997.csv"))}
    def read_audit(v):
        return {int(r["domain"]): r for r in csv.DictReader(
            open(dune_topo_root("1984-start") / v
                 / "HAT_seaward_row_insert_audit.csv"))}

    aud = read_audit(args.version)
    prev = read_audit(args.compare)
    nrow = {d: int(r["n_rows_inserted"]) for d, r in aud.items()}
    nprev = {d: int(r["n_rows_inserted"]) for d, r in prev.items()}
    road = {d: r.get("has_road", "1") == "1" for d, r in aud.items()}

    def fnum(x):
        try:
            return float(x)
        except ValueError:
            return float("nan")

    before = {d: fnum(r["setback_raw_before_m"]) for d, r in aud.items()}
    after = {d: fnum(r["setback_raw_after_m"]) for d, r in aud.items()}

    def n_of(m):
        return max(0, int(round(m / 10.0)))

    fig = plt.figure(figsize=(11.4, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])

    # ---- (a) and (b): the two blocks, zoomed -----------------------------
    for j, (blk, name, yr, letter) in enumerate(
            ((BLOCK_A, "GIS 9–14 · inter-village", 1999, "a"),
             (BLOCK_B, "GIS 84–87 · Pea Island", 1989, "b"))):
        ax = fig.add_subplot(gs[0, j])
        d = [g for g in blk if g in date]
        xi = np.arange(len(d))
        dm = np.array([date[g] for g in d])
        n_now = np.array([nrow.get(g, 0) for g in d])
        n_12 = np.array([n_of(date[g] * YEARS_WANTED / YEARS_MEASURED) for g in d])

        cols = [C["ACCENT"] if nrow.get(g, 0) > 0 else "#9a9a9a" for g in d]
        ax.bar(xi, dm, 0.62, color=cols, zorder=3)
        for k, g in enumerate(d):
            lab = "+{}".format(n_now[k]) if n_now[k] else "no change"
            ax.text(k, dm[k] + (1.5 if dm[k] >= 0 else -4.5), lab,
                    ha="center", va="bottom" if dm[k] >= 0 else "top",
                    fontsize=8, color=C["ACCENT"] if n_now[k] else "#5a5a5a",
                    fontweight="bold" if n_now[k] else "normal")
            if n_12[k] != n_now[k]:
                ax.plot([k], [date[g] * YEARS_WANTED / YEARS_MEASURED], "v",
                        color=C["REF"], ms=7, zorder=5)
                ax.text(k, date[g] * YEARS_WANTED / YEARS_MEASURED - 5.5,
                        "→+{}".format(n_12[k]), ha="center", va="top",
                        fontsize=7.5, color=C["REF"])
        for t in (5, 15, 25, 35, 45, 55):
            ax.axhline(t, color="#dddddd", lw=0.6, zorder=0)
        ax.axhline(0, color="k", lw=1.0, zorder=2)
        ax.axhline(5, color="#888888", ls=":", lw=1.1, zorder=2)
        # its own left margin, so the label never lands on a bar
        ax.set_xlim(-1.35, len(d) - 0.4)
        ax.text(-1.28, 5, " round()\n threshold", fontsize=6.6, ha="left",
                va="center", color="#666666")
        ax.set_xticks(xi)
        ax.set_xticklabels(d)
        ax.set_xlabel("GIS domain")
        ax.set_ylabel("measured 1984→1997\ndune-line retreat (m)")
        ax.set_ylim(min(-12, dm.min() - 10), max(dm.max() + 14, 20))
        ax.set_title(panel_title(letter, "{} · relocated {}".format(name, yr)))
        ax.grid(False)

    # ---- (c) setback before/after, the two blocks ------------------------
    ax3 = fig.add_subplot(gs[1, 0])
    d = [g for g in SCOPE if g in aud]
    xi = np.arange(len(d))
    ax3.bar(xi - 0.2, [before[g] for g in d], 0.4, color="#9a9a9a",
            label="before insert (measured)")
    ax3.bar(xi + 0.2, [after[g] for g in d], 0.4, color=C["ACCENT"],
            label="{}, after insert".format(args.version))
    for k, g in enumerate(d):
        if before[g] < 0:
            ax3.text(k - 0.2, before[g] - 1.5, "would fail", ha="center",
                     va="top", fontsize=6.4, rotation=90, color="k")
        if nrow.get(g, 0) == 0:
            ax3.text(k, max(before[g], after[g]) + 1.5, "unchanged",
                     ha="center", va="bottom", fontsize=6.2, rotation=90,
                     color="#5a5a5a")
    ax3.axhline(0, color="k", lw=1.0)
    lo3 = min(0.0, min(before[g] for g in d))
    hi3 = max(after[g] for g in d)
    ax3.set_ylim(lo3 - 6, hi3 + 22)          # headroom for the rotated labels
    ax3.set_xticks(xi)
    ax3.set_xticklabels(d, fontsize=7.5)
    ax3.set_xlabel("GIS domain")
    ax3.set_ylabel("road setback (m)")
    ax3.set_title(panel_title("c", "What the insert does to the setback"))
    ax3.legend(fontsize=7.6)
    ax3.grid(alpha=.3, axis="y")

    # ---- (d) island-wide: the scope decision, removed ---------------------
    # Under the block scope, "unchanged" meant two different things: a domain
    # whose measurement said no cell was needed, and a domain that was never
    # asked. 30 domains sat in the second category with N >= 1, several of them
    # a larger measured retreat than GIS 86's. Island-wide scope collapses that
    # into one category, and this panel is what shows the two apart.
    ax4 = fig.add_subplot(gs[1, 1])
    gis = np.array(sorted(date))
    dm = np.array([date[g] for g in gis])
    now = np.array([nrow.get(g, 0) > 0 for g in gis])
    was = np.array([nprev.get(g, 0) > 0 for g in gis])

    m_kept = now & was
    m_new = now & (~was)
    m_zero = ~now
    ax4.bar(gis[m_kept], dm[m_kept], 0.9, color=C["ACCENT"], zorder=3,
            label="already in {} ({})".format(args.compare, int(m_kept.sum())))
    ax4.bar(gis[m_new], dm[m_new], 0.9, color="#c8a06a", zorder=3,
            label="added by island-wide scope ({})".format(int(m_new.sum())))
    ax4.bar(gis[m_zero], dm[m_zero], 0.9, color="#dcdcdc", zorder=1,
            label="N rounds to 0, no land missing ({})".format(int(m_zero.sum())))
    ax4.axhline(0, color="k", lw=1.0)
    ax4.axhline(5, color="#888888", ls=":", lw=1.0)
    for lo, hi in ((9, 14), (84, 87)):
        ax4.axvspan(lo - .5, hi + .5, color="#ffe9b0", alpha=.45, zorder=0)
    noroad = [g for g in gis[now] if not road.get(g, True)]
    lo4, hi4 = float(dm.min()), float(dm.max())
    ax4.set_ylim(lo4 - 6, hi4 + 26)
    if noroad:
        ax4.plot(noroad, [hi4 + 6] * len(noroad), "v", color="#444444", ms=4,
                 zorder=6, clip_on=False,
                 label="no NC-12 in the domain ({})".format(len(noroad)))
    ax4.set_xlim(0, 91)
    ax4.set_xlabel("GIS domain   (1 = Cape Point → 90 = north Pea Island)")
    ax4.set_ylabel("measured retreat (m)")
    ax4.set_title(panel_title("d", "Island-wide: every domain with a measured "
                                   "shift now gets its rows"))
    ax4.legend(fontsize=6.6, loc="lower left", ncol=2, framealpha=.92)
    ax4.grid(alpha=.3, axis="y")

    fig_h = fig.get_figheight()
    fig.suptitle("Where the insert changes the domain, and where it does not ({})".format(args.version),
                 fontsize=11, fontweight="bold", x=0.055, ha="left",
                 y=1 - 0.24 / fig_h)
    caption(fig,
            "Green triangles in (a) and (b): N if the measurement were scaled "
            "13 → 12 years. The dune lines are 1984/1997 while the DEM surface "
            "at row 0 and at every inserted cell is 1996 ALACE\n"
            "(verified 100% of cells); the 1996 and 1997 aerials are 363 days "
            "apart. That bias is recorded, not corrected — scaling would assume "
            "steady retreat across 13 storm-driven years. It moves N\n"
            "at GIS 12, 84 and 85. N is floored at 0: where the shoreline "
            "prograded no 1984 land is missing from the 1996 survey, so nothing "
            "is inserted and nothing is removed.",
            y=0.075 / fig_h, size=7.4)
    fig.subplots_adjust(top=1 - 0.60 / fig_h, bottom=0.78 / fig_h,
                        left=0.075, right=0.98, hspace=0.46, wspace=0.24)

    out = Path(args.out) if args.out else (
        insert_figures_dir("1984-start") / "HAT_where_inserts_occur.png")
    fig.savefig(out)
    print("wrote {}".format(out))
    print("  {}: changed {} domains ({} carried over from {}, {} added by "
          "island-wide scope); {} unchanged because N rounds to 0"
          .format(args.version, int(now.sum()), int(m_kept.sum()),
                  args.compare, int(m_new.sum()), int(m_zero.sum())))
    print("  of the changed domains, {} have no NC-12 in them"
          .format(len(noroad)))


if __name__ == "__main__":
    main()
