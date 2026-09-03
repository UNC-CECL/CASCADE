#!/usr/bin/env python3
r"""
HAT_plot_how_N_is_determined.py
==============================================================================
How the number of inserted rows, N, is measured. One domain, the whole chain.

THE QUANTITY WANTED
    How far the dune line moved between 1984 and the surveyed surface, in 10 m
    cells. That distance is how far interior row 0 has to move seaward for the
    1984 roadway to sit its true distance behind the dune.

WHY IT IS A DIFFERENCE OF TWO LINES AND NOT ONE MEASUREMENT
    The obvious measurement -- 1984 dune line against the extractor's interior
    row 0 -- confounds two things:

        row 0 - line_1984  =  (how far the island moved)          DATE
                           +  (digitized line vs the model's row 0)  FEATURE

    The feature term is not small. The 1997 line, measured the same way against
    the same row 0, sits +16.2 m seaward of it island-wide (IQR +12.8 to +21.0)
    -- a near-constant offset, which is what a definitional difference looks
    like. Island-wide it accounts for ~85% of the naive number.

    Differencing two digitized lines cancels it exactly:

        (row0 - line_1984) - (row0 - line_1997) = line_1997 - line_1984

    Row 0 drops out algebraically, so no assumption about where row 0 sits
    survives into N. And because the same person digitized the same feature from
    the same kind of imagery at both dates, the definitional term cancels too.

    N = round( median over profiles / 10 m ),  floored at 0.

USAGE
    python HAT_plot_how_N_is_determined.py [--domain 85]
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
from hat_topo_version import duneline_shift_dir  # noqa: E402
from hat_topo_version import insert_figures_dir  # noqa: E402
from hat_figure_style import (apply_style, C, caption,   # noqa: E402
                              panel_title)

# Resolved through hat_topo_version.duneline_shift_dir - ONE definition
# of a path that eight scripts used to build by hand. Moved under
# row-insert-scope/ on 2026-09-03.
S = duneline_shift_dir("1984-start")


def per_profile(fname, D):
    out = {}
    for r in csv.DictReader(open(S / fname)):
        if int(r["domain"]) == D:
            out[int(r["profile"])] = (float(r["duneline_cell"]),
                                      int(r["interior_row0_cell"]),
                                      float(r["shift_m"]))
    return out


def domain_row(fname, D):
    for r in csv.DictReader(open(S / fname)):
        if int(r["domain"]) == D:
            return r
    return {}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=int, default=85)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    D = args.domain
    apply_style()

    p84 = per_profile("duneline_shift_1984_profiles.csv", D)
    p97 = per_profile("duneline_shift_1997_profiles.csv", D)
    common = sorted(set(p84) & set(p97))
    prof = np.array(common)
    line84 = np.array([p84[k][0] for k in common])
    line97 = np.array([p97[k][0] for k in common])
    row0 = np.array([p84[k][1] for k in common], dtype=float)
    s84 = np.array([p84[k][2] for k in common])
    s97 = np.array([p97[k][2] for k in common])

    med84, med97 = float(np.median(s84)), float(np.median(s97))
    date = med84 - med97
    n_rows = max(0, int(round(date / 10.0)))
    per_prof_date = float(np.median(s84 - s97))

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.0))
    (ax1, ax2), (ax3, ax4) = axes

    # ---- (a) the geometry, per profile ---------------------------------
    ax1.plot(prof, row0, "-", color="#333333", lw=2.0,
             label="interior row 0 (the model's reference)")
    ax1.plot(prof, line84, "-", color=C["ACCENT"], lw=1.6,
             label="1984 dune line")
    ax1.plot(prof, line97, "-", color="#1b6ca8", lw=1.6,
             label="1997 dune line")
    ax1.fill_between(prof, line84, line97, color=C["ADDED_FILL"], alpha=0.85,
                     zorder=0, label="the date term (what N measures)")
    ax1.invert_yaxis()
    ax1.set_xlabel("alongshore profile")
    ax1.set_ylabel("cross-shore cell\n(0 = ocean)")
    ax1.set_title(panel_title("a", "Where the three references fall, "
                                   "profile by profile"))
    ax1.legend(fontsize=7.2, loc="upper center", ncol=2,
               bbox_to_anchor=(0.5, -0.16), frameon=False)
    ax1.grid(alpha=.3)

    # ---- (b) the two shifts, and their difference ----------------------
    ax2.plot(prof, s84, "-", color=C["ACCENT"], lw=1.5,
             label="row 0 − 1984 line  (total)")
    ax2.plot(prof, s97, "-", color="#1b6ca8", lw=1.5,
             label="row 0 − 1997 line  (feature)")
    ax2.plot(prof, s84 - s97, "-", color=C["REF"], lw=1.8,
             label="difference  (date)")
    for v, col in ((med84, C["ACCENT"]), (med97, "#1b6ca8"), (date, C["REF"])):
        ax2.axhline(v, color=col, ls="--", lw=1.0, alpha=0.8)
    ax2.axhline(0, color="k", lw=0.9)
    ax2.set_xlabel("alongshore profile")
    ax2.set_ylabel("metres")
    ax2.set_title(panel_title("b", "The two shifts and what is left when "
                                   "one is subtracted"))
    ax2.legend(fontsize=7.6, loc="upper right")
    ax2.grid(alpha=.3)

    # ---- (c) the decomposition as a bar ---------------------------------
    ax3.barh([2], [med84], color=C["ACCENT"], height=0.55)
    ax3.barh([1], [med97], color="#1b6ca8", height=0.55)
    ax3.barh([0], [date], color=C["REF"], height=0.55)
    for y, v in ((2, med84), (1, med97), (0, date)):
        ax3.text(v + 1.5, y, "{:.1f} m".format(v), va="center",
                 fontsize=8.5, fontweight="bold")
    # The tick labels carry the meaning, so a legend here would repeat itself
    # and was colliding with the caption.
    ax3.set_yticks([0, 1, 2])
    ax3.set_yticklabels(["DATE\n= N",
                         "feature\nrow 0 − 1997",
                         "total\nrow 0 − 1984"])
    ax3.set_xlim(0, max(med84, date) * 1.28)
    ax3.set_xlabel("metres (median over {} profiles)".format(len(common)))
    ax3.set_title(panel_title("c", "Row 0 cancels; the definitional offset "
                                   "cancels with it"))
    ax3.grid(alpha=.3, axis="x")

    # ---- (d) from metres to rows ----------------------------------------
    ax4.axis("off")
    lines = [
        "GIS {}".format(D), "",
        "median over {} profiles".format(len(common)),
        "  total   row 0 − 1984 line      {:6.1f} m".format(med84),
        "  feature row 0 − 1997 line      {:6.1f} m".format(med97),
        "  " + "-" * 40,
        "  DATE    1997 line − 1984 line  {:6.1f} m".format(date), "",
        "  ÷ 10 m cell size               {:6.2f} cells".format(date / 10.0),
        "  round()                        {:6d} rows".format(n_rows), "",
        "robustness",
        "  median of per-profile diffs    {:6.1f} m".format(per_prof_date),
        "    (vs {:.1f} m for the difference of".format(date),
        "     medians — {:+.1f} m, under a".format(per_prof_date - date),
        "     tenth of a cell)", "",
        "  p10–p90 of the total           {:.0f}–{:.0f} m".format(
            float(np.percentile(s84, 10)), float(np.percentile(s84, 90))),
        "  p10–p90 of the feature         {:.0f}–{:.0f} m".format(
            float(np.percentile(s97, 10)), float(np.percentile(s97, 90))),
    ]
    ax4.text(0.0, 0.99, "\n".join(lines), family="monospace", fontsize=9.2,
             va="top", transform=ax4.transAxes)
    ax4.set_title(panel_title("d", "From metres to rows"))

    fig_h = fig.get_figheight()
    fig.suptitle("GIS {}  ·  how the number of inserted rows is determined"
                 .format(D), fontsize=11, fontweight="bold", x=0.055,
                 ha="left", y=1 - 0.24 / fig_h)
    caption(fig,
            "The 1984 line measured against interior row 0 mixes island "
            "movement with the offset between a digitized line and the model's "
            "row 0. Differencing two\n"
            "digitized lines removes both row 0 and that definitional offset "
            "algebraically, leaving date. N is the median of the per-profile "
            "difference, in 10 m cells.",
            y=0.075 / fig_h, size=7.4)
    fig.subplots_adjust(top=1 - 0.62 / fig_h, bottom=0.62 / fig_h,
                        left=0.085, right=0.98, hspace=0.58, wspace=0.22)

    out = Path(args.out) if args.out else (
        insert_figures_dir("1984-start") / "HAT_how_N_determined_GIS{}.png".format(D))
    fig.savefig(out)
    print("wrote {}".format(out))
    print("  total {:.1f} m = feature {:.1f} + date {:.1f}  ->  N = {} rows"
          .format(med84, med97, date, n_rows))


if __name__ == "__main__":
    main()
