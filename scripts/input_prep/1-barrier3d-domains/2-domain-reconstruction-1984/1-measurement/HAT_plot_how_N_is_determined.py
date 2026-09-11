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
from hat_topo_version import insert_figures_dir_for_domain  # noqa: E402
from hat_figure_style import (apply_style, C, C_1984, C_1984_FILL,   # noqa: E402
                              C_1997, INK, caption, figsize, open_frame,
                              save, _title)

# Resolved through hat_topo_version.duneline_shift_dir - ONE definition
# of a path that eight scripts used to build by hand. Moved under
# 2-domain-reconstruction-1984/ on 2026-09-03.
S = duneline_shift_dir("1984-start")
L84, L97, LROW0, L_DATE = C_1984, C_1997, C["ROAD"], C["REF"]


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
    p_tot = (float(np.percentile(s84, 10)), float(np.percentile(s84, 90)))
    p_fea = (float(np.percentile(s97, 10)), float(np.percentile(s97, 90)))

    # (a) and (b) side by side, (c) the bar below them. The arithmetic that
    # used to be a monospace panel (d) is in the caption.
    fig = plt.figure(figsize=figsize("double", aspect=0.74),
                     constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.62])
    ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    # ---- (a) the geometry, per profile ---------------------------------
    ax1.plot(prof, row0, "-", color=LROW0, lw=1.6, label="interior row 0")
    ax1.plot(prof, line84, "-", color=L84, lw=1.3, label="1984 dune line")
    ax1.plot(prof, line97, "-", color=L97, lw=1.3, label="1997 dune line")
    ax1.fill_between(prof, line84, line97, color=C_1984_FILL, alpha=0.7,
                     lw=0, zorder=0, label="band between the lines (the date term, N)")
    ax1.invert_yaxis()
    ax1.set_xlabel("alongshore profile")
    ax1.set_ylabel("cross-shore cell (0 = ocean)")
    _title(ax1, 0, "the three references")
    ax1.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.22),
               frameon=False, columnspacing=1.0)
    ax1.grid(axis="y")
    open_frame(ax1)

    # ---- (b) the two shifts, and their difference ----------------------
    ax2.plot(prof, s84, "-", color=L84, lw=1.2,
             label="row 0 − 1984 line (total)")
    ax2.plot(prof, s97, "-", color=L97, lw=1.2,
             label="row 0 − 1997 line (feature)")
    ax2.plot(prof, s84 - s97, "-", color=L_DATE, lw=1.5,
             label="difference (date)")
    for v, col in ((med84, L84), (med97, L97), (date, L_DATE)):
        ax2.axhline(v, color=col, ls="--", lw=0.8, alpha=0.8)
    ax2.plot([], [], color=INK, ls="--", lw=0.8, label="median over the profiles")
    ax2.axhline(0, color=INK, lw=0.6)
    ax2.set_xlabel("alongshore profile")
    ax2.set_ylabel("offset (m)")
    _title(ax2, 1, "the two shifts, differenced")
    ax2.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.22),
               frameon=False, columnspacing=1.0)
    ax2.grid(axis="y")
    open_frame(ax2)

    # ---- (c) the decomposition as a bar ---------------------------------
    ax3.barh([2], [med84], color=L84, height=0.55)
    ax3.barh([1], [med97], color=L97, height=0.55)
    ax3.barh([0], [date], color=L_DATE, height=0.55)
    for y, v in ((2, med84), (1, med97), (0, date)):
        ax3.text(v + 0.8, y, "{:.1f} m".format(v), va="center", ha="left",
                 fontsize=8, fontweight="bold", color=INK)
    # The tick labels carry the meaning, so a legend here would repeat itself.
    ax3.set_yticks([0, 1, 2])
    ax3.set_yticklabels(["date (N)\n1997 − 1984 line",
                         "feature\nrow 0 − 1997 line",
                         "total\nrow 0 − 1984 line"])
    ax3.set_xlim(0, max(med84, date) * 1.15)
    ax3.set_xlabel("offset (m, median over {} profiles)".format(len(common)))
    _title(ax3, 2, "row 0 cancels, and the offset with it")
    ax3.grid(axis="x")
    open_frame(ax3)

    caption(fig,
            "How the number of inserted rows, N, is determined at GIS {}. The "
            "1984 line measured against interior row 0 mixes island movement "
            "with the offset between a digitized line and the model's row 0; "
            "differencing two digitized lines removes both row 0 and that "
            "definitional offset algebraically, leaving the date term. (a) The "
            "1984 (red) and 1997 (blue) dune lines and interior row 0 (black) "
            "per profile, in grid cells; the band between the lines is what N "
            "measures. (b) Each line's offset from row 0 per profile and the "
            "difference, with the medians dashed. (c) The medians over the {} "
            "profiles: total (row 0 − 1984 line) {:.1f} m, feature (row 0 − "
            "1997 line) {:.1f} m, date (1997 line − 1984 line) {:.1f} m. "
            "Divided by the 10 m cell that is {:.2f} cells, rounded to N = {} "
            "rows. Robustness: the median of the per-profile differences is "
            "{:.1f} m, {:+.1f} m from the difference of medians, under a tenth "
            "of a cell; the p10–p90 range of the total is {:.0f}–{:.0f} m and "
            "of the feature {:.0f}–{:.0f} m."
            .format(D, len(common), med84, med97, date, date / 10.0, n_rows,
                    per_prof_date, per_prof_date - date, p_tot[0], p_tot[1],
                    p_fea[0], p_fea[1]))

    out = Path(args.out) if args.out else (
        insert_figures_dir_for_domain("1984-start", "1-measurement", D)
        / "HAT_how_N_determined_GIS{}.png".format(D))
    save(fig, out)
    print("wrote {}".format(out))
    print("  total {:.1f} m = feature {:.1f} + date {:.1f}  ->  N = {} rows"
          .format(med84, med97, date, n_rows))


if __name__ == "__main__":
    main()
