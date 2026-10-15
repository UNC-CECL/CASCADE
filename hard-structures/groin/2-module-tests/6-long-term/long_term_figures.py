#!/usr/bin/env python3
"""
The long-run figures, drawn from the tables long_term_test.py writes.

    python long_term_figures.py   ->  figures/long_term_<n>_<what>.png

Captions go to figures/supporting/CAPTIONS.md. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / d) for d in ("1-straight-coast/1-run", "4-groin-field", "5-fine-grid-reference")]
import straight_coast_test as sc  # noqa: E402
import groin_field_test as gf  # noqa: E402
import long_term_test as lt  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
TABLES = HERE / "tables"
FIGS = HERE / "figures"
REF_KW = dict(color=C["BASE"], lw=4.0, alpha=0.55, zorder=1)
# -----------------------------------------------------------------------------


def module_kw(m):
    return dict(color=sc.MODULES[m]["col"], lw=1.2, marker=gf.MARKERS[m], ms=2.8)


# Figure 1: the gap across the groin over the whole run, one panel per orientation
def gap_over_time(frame):
    ths = list(lt.ORIENTATIONS)
    fig, axes = plt.subplots(2, 4, figsize=figsize("double", height=5.0), constrained_layout=True)
    last = int(frame.year.max())
    for i, th in enumerate(ths):
        ax = axes.flat[i]
        f = frame[frame.theta0_deg == th]
        r = f[f.source == "reference"]
        ax.plot(r.year, r.gap_m, **REF_KW)
        for m in lt.MODULES:
            s_ = f[f.source == m]
            ax.plot(s_.year, s_.gap_m, markevery=25, **module_kw(m))
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_xlim(0, last)
        ax.grid(axis="y")
        open_frame(ax)
        _letter_inside(ax, i)
        ax.set_title(f"Coast at {th:+d}$^\\circ$", loc="center", fontsize=8)
        if i % 4 == 0:
            ax.set_ylabel("Gap across the groin (m)")
        if i >= 4:
            ax.set_xlabel("Year")
    ax = axes.flat[-1]
    ax.axis("off")
    handles = [Line2D([], [], label="reference, 110 m groins\n(25 m cells)", **REF_KW)]
    handles += [Line2D([], [], label=sc.MODULES[m]["label"], **module_kw(m)) for m in lt.MODULES]
    ax.legend(handles=handles, loc="center", fontsize=7, frameon=False)
    fig.suptitle(f"Does each version's groin hold steady over {last} years?")
    png = FIGS / "long_term_1_gap_over_time.png"
    save(fig, png, dpi=250, close=True)
    return png


# Figure 2: the end state against orientation: net sand, growth in the last quarter, steepest angle
def end_state(summ, last):
    q = last // 4
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=2.9), constrained_layout=True)
    panels = (("net_change_end_m", f"Net sand change at year {last}\n(m of summed shoreline)", "Sand conserved?", True),
              ("gap_growth_last_quarter_m", f"Gap change, years {last - q}-{last} (m)", "Still growing?", True),
              ("max_angle_deg", "Steepest shoreline angle,\nany face, any year ($^\\circ$)", "Stable?", False))
    for i, (col, ylab, title, sym) in enumerate(panels):
        ax = axes[i]
        for m in lt.MODULES + ("reference",):
            s_ = summ[summ.source == m].sort_values("theta0_deg")
            if s_[col].isna().all():
                continue
            kw = dict(color=C["BASE"], lw=2.4, alpha=0.6) if m == "reference" else module_kw(m)
            ax.plot(s_.theta0_deg, s_[col], **kw)
        if sym:
            ax.set_yscale("symlog", linthresh=10)
            ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlabel("Shoreline orientation ($^\\circ$)")
        ax.set_ylabel(ylab, fontsize=7)
        ax.grid(axis="y")
        open_frame(ax)
        _letter_inside(ax, i)
        ax.set_title(title, loc="center")
    handles = [Line2D([], [], color=C["BASE"], lw=2.4, alpha=0.6, label="reference")]
    handles += [Line2D([], [], label=sc.MODULES[m]["label"], **module_kw(m)) for m in lt.MODULES]
    fig.legend(handles=handles, loc="outside lower center", ncol=6, fontsize=6.8, frameon=False)
    fig.suptitle(f"After {last} years: is sand conserved, has the groin levelled off, and is it stable?")
    png = FIGS / "long_term_2_end_state.png"
    save(fig, png, dpi=250, close=True)
    return png


# Run: read the tables, draw both figures, write their captions
def main():
    apply_style()
    frame = pd.read_csv(TABLES / "long_term_by_year.csv").dropna(subset=["year"])
    summ = pd.read_csv(TABLES / "long_term_summary.csv")
    FIGS.mkdir(parents=True, exist_ok=True)
    last = int(frame.year.max())
    s = summ.set_index(["source", "theta0_deg"])
    png = gap_over_time(frame)
    record_caption(png, (
        f"**Tests:** whether each version of the groin module holds a steady gap across the groin over {last} "
        "years, or keeps growing, collapses or runs away, and whether it keeps tracking a real-sized groin field. "
        "**How to read:** one groin on a straight coast at seven orientations, BRIE's alongshore solve alone "
        f"(option A waves; no storms, no Barrier3D, no sea-level rise), {lt.NY} domains. Lines: the gap between "
        "the updrift and downdrift domains, the quantity measured at Buxton. Broad grey: the 25 m reference with "
        "three 110 m groins, averaged to domains. **Shows:** see the README beside this figure for the numbers; "
        f"at -20 degrees the corrected version's gap is {s.loc[('corrected', -20), f'gap_{last}_m']:.0f} m at year "
        f"{last} against the reference's {s.loc[('reference', -20), f'gap_{last}_m']:.0f} m (tables/)."))
    png = end_state(summ, last)
    record_caption(png, (
        f"**Tests:** at the end of a {last}-year run, whether each version of the groin module still conserves "
        "sand, has stopped growing, and has kept the shoreline within the angles BRIE can handle. **How to read:** "
        "against the shoreline orientation: (a) the net change in sand over the coast at the last year, as summed "
        f"shoreline (symmetric log scale; zero is conserved); (b) the change in the gap over the last {last // 4} "
        "years (near zero is levelled off); (c) the steepest shoreline angle across any face in any year (BRIE's "
        "diffusivity table ends at 90 degrees, where a runaway stops the model). Grey: the 25 m reference. "
        "**Shows:** the numbers are in tables/long_term_summary.csv and the README beside this figure."))


if __name__ == "__main__":
    main()
