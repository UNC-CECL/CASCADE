#!/usr/bin/env python3
"""
The fine-grid reference figures, drawn from the tables fine_grid_reference.py writes.

    python fine_grid_reference_figures.py   ->  figures/fine_reference_<n>_<what>.png

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
sys.path[:0] = [str(HERE.parent / "1-straight-coast"), str(HERE.parent / "4-groin-field")]
import straight_coast_test as sc  # noqa: E402
import groin_field_test as gf  # noqa: E402
import fine_grid_reference as fr  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
TABLES = HERE / "tables"
FIGS = HERE / "figures"
SHOWN = (-20, 0, 20)
REF_KW = dict(color=C["BASE"], lw=4.0, alpha=0.55, zorder=1)   # a broad grey band behind the modules
REF_FINE_KW = dict(color=C["BASE"], lw=0.7, ls=(0, (1, 1.5)), zorder=1)
# -----------------------------------------------------------------------------


def module_kw(m):
    return dict(color=sc.MODULES[m]["col"], lw=1.3, marker=gf.MARKERS[m], ms=3)


# Figure 1: the gap across the groin over time (top) and the shoreline at the last year (bottom)
def gap_and_profiles(frame, prof, score):
    fig, axes = plt.subplots(2, len(SHOWN), figsize=figsize("double", height=6.2), constrained_layout=True)
    last = int(frame.year.max())
    for c_, th in enumerate(SHOWN):
        ax = axes[0, c_]
        f = frame[frame.theta0_deg == th]
        ref = f[f.source == "reference"]
        for length, kw in ((fr.LENGTH_SENSITIVITY_M[0], dict(ls=(0, (2, 2)))), (fr.LENGTH_SENSITIVITY_M[1], dict(ls=(0, (6, 2))))):
            r = ref[ref.groin_length_m == length]
            ax.plot(r.year, r.gap_m, color=C["BASE"], lw=1.0, zorder=1, **kw)
        r = ref[ref.groin_length_m == fr.GROIN_LENGTH_M]
        ax.plot(r.year, r.gap_m, **REF_KW)
        for m in fr.MODULES:
            s_ = f[f.source == m]
            ax.plot(s_.year, s_.gap_m, markevery=7, **module_kw(m))
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_xlabel("Year")
        ax.grid(axis="y")
        open_frame(ax)
        _letter_inside(ax, c_)
        ax.set_title(f"Coast at {th:+d}$^\\circ$", loc="center")
        ax = axes[1, c_]
        p = prof[prof.theta0_deg == th]
        fine = p[p.source == "reference_fine"]
        ax.plot(fine.domain - fr.LO - 1, fine.seaward_m, **REF_FINE_KW)
        r = p[p.source == "reference"]
        ax.step(r.domain - fr.LO - 1, r.seaward_m, where="mid", **REF_KW)
        for m in fr.MODULES:
            s_ = p[p.source == m]
            ax.plot(s_.domain - fr.LO - 1, s_.seaward_m, **module_kw(m))
        ax.axvline(0, color=C["GROIN"], lw=0.8)
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_xlim(-8, 8)
        ax.grid(axis="y")
        open_frame(ax)
        _letter_inside(ax, len(SHOWN) + c_)
    axes[0, 0].set_ylabel("Gap across the groin\n(updrift minus downdrift\ndomain, m)")
    axes[1, 0].set_ylabel(f"Shoreline change after\n{last} yr (m, seaward +)")
    axes[1, len(SHOWN) // 2].set_xlabel("Domains from the groin face (updrift right)")
    handles = [Line2D([], [], label=f"reference, {fr.GROIN_LENGTH_M:g} m groins", **REF_KW),
               Line2D([], [], color=C["BASE"], lw=1.0, ls=(0, (2, 2)), label=f"reference, {fr.LENGTH_SENSITIVITY_M[0]:g} m"),
               Line2D([], [], color=C["BASE"], lw=1.0, ls=(0, (6, 2)), label=f"reference, {fr.LENGTH_SENSITIVITY_M[1]:g} m"),
               Line2D([], [], label="reference, 25 m cells", **REF_FINE_KW)]
    handles += [Line2D([], [], label=sc.MODULES[m]["label"], **module_kw(m)) for m in fr.MODULES]
    handles.append(Line2D([], [], color=C["GROIN"], lw=0.8, label="groin face"))
    fig.legend(handles=handles, loc="outside lower center", ncol=5, fontsize=6.8, frameon=False)
    fig.suptitle("Does a 500 m groin module reproduce a real-sized groin field?")
    png = FIGS / "fine_reference_1_gap_and_profiles.png"
    save(fig, png, dpi=250, close=True)
    sm = score.set_index(["module", "theta0_deg"])
    ref_end = frame[(frame.source == "reference") & (frame.year == last)].pivot_table(
        index="theta0_deg", columns="groin_length_m", values="gap_m")
    record_caption(png, (
        "**Tests:** whether each version of the groin module, with one groin on BRIE's 500 m grid, reproduces "
        "what a real-sized groin field does to the two domains either side of it: the gap between them over time, "
        "and the shoreline after "
        f"{last} years. **How to read:** the reference is a one-line model on 25 m cells with the same wave-climate "
        f"drift as BRIE, and three groins {fr.GROIN_OFFSETS_M[1]:g} m apart inside the updrift domain, as at Buxton; "
        f"each groin first holds {sc.BLOCKING_FRACTION} of the drift and lets more sand round its tip as the "
        "fillet reaches it (GENESIS-style bypassing). Broad grey: the reference with "
        f"{fr.GROIN_LENGTH_M:g} m groins, averaged to domains; dashed grey: the same with {fr.LENGTH_SENSITIVITY_M[0]:g} "
        f"and {fr.LENGTH_SENSITIVITY_M[1]:g} m groins; dotted grey in (d-f): the reference on its own 25 m cells. "
        "Coloured: the four modules. (a-c) the updrift domain minus the downdrift domain, the quantity measured "
        "at Buxton. (d-f) shoreline change, seaward positive, updrift side on the right. **Shows:** the reference "
        f"field fills to its tips within 5 to 20 years at -30 to +20 degrees and then holds a gap of about twice the groin length "
        f"({ref_end.loc[-20, fr.LENGTH_SENSITIVITY_M[0]]:.0f}, {ref_end.loc[-20, fr.GROIN_LENGTH_M]:.0f} and "
        f"{ref_end.loc[-20, fr.LENGTH_SENSITIVITY_M[1]]:.0f} m at -20 degrees for the three lengths); at +30 degrees, where "
        "the drift is a tenth of its grid-aligned value, it is still filling at year 50. No module has a groin length, so "
        "none levels off at a length. At -20 degrees the corrected version "
        f"ends at {sm.loc[('corrected', -20), 'module_gap_end_m']:.0f} m against the reference's "
        f"{sm.loc[('corrected', -20), 'ref_gap_end_m']:.0f} m; at +20 degrees it holds only "
        f"{sm.loc[('corrected', 20), 'module_gap_end_m']:.0f} m, the drift version "
        f"{sm.loc[('drift', 20), 'module_gap_end_m']:.0f} m, and the pinned version builds its gap on the wrong side "
        f"({sm.loc[('pinned', 20), 'module_gap_end_m']:.0f} m) (tables/)."))


# Figure 2: the error against the reference, by orientation
def errors(score):
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=2.9), constrained_layout=True)
    for i, (col, ylab, title) in enumerate((
            ("gap_rmse_m", "RMS error in the gap,\nyears 1-50 (m)", "The gap over time"),
            ("profile_rmse_m", f"RMS error in the shoreline,\n{2 * fr.WINDOW + 2} domains at year 50 (m)",
             "The shoreline at the end"))):
        ax = axes[i]
        for m in fr.MODULES:
            s_ = score[score.module == m].sort_values("theta0_deg")
            ax.plot(s_.theta0_deg, s_[col], **module_kw(m))
        ax.set_xlabel("Shoreline orientation ($^\\circ$)")
        ax.set_ylabel(ylab)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y")
        open_frame(ax)
        _letter_inside(ax, i)
        ax.set_title(title, loc="center")
    fig.legend(handles=[Line2D([], [], label=sc.MODULES[m]["label"], **module_kw(m)) for m in fr.MODULES],
               loc="outside lower center", ncol=4, fontsize=7, frameon=False)
    fig.suptitle("How far is each version from the real-sized field?")
    png = FIGS / "fine_reference_2_errors.png"
    save(fig, png, dpi=250, close=True)
    mean = score.groupby("module")[["gap_rmse_m", "profile_rmse_m"]].mean().round(0)
    record_caption(png, (
        "**Tests:** across shoreline orientations from -30 to +30 degrees, how far each version of the groin module, "
        f"on BRIE's 500 m grid, sits from the 25 m reference field with {fr.GROIN_LENGTH_M:g} m groins. **How to "
        "read:** (a) the root-mean-square difference in the gap between the updrift and downdrift domains over "
        f"years 1 to 50; (b) the same for the shoreline after 50 years, over the {2 * fr.WINDOW + 2} domains "
        "around the groin. Lower is closer. **Shows:** averaged over the seven orientations, the gap error is "
        f"{mean.loc['corrected', 'gap_rmse_m']:.0f} m for the corrected version, {mean.loc['drift', 'gap_rmse_m']:.0f} m "
        f"drift, {mean.loc['dipole', 'gap_rmse_m']:.0f} m dipole and {mean.loc['pinned', 'gap_rmse_m']:.0f} m pinned; "
        f"the shoreline error {mean.loc['corrected', 'profile_rmse_m']:.0f}, {mean.loc['drift', 'profile_rmse_m']:.0f}, "
        f"{mean.loc['dipole', 'profile_rmse_m']:.0f} and {mean.loc['pinned', 'profile_rmse_m']:.0f} m. The corrected "
        "version is closest on coasts turned against the drift and drifts away from the reference on coasts turned "
        "with it, where the reference still fills to its tips and the module traps little "
        "(tables/fine_reference_scores.csv)."))
    return mean


# Run: read the tables, draw both figures
def main():
    apply_style()
    frame = pd.read_csv(TABLES / "fine_reference_by_year.csv")
    prof = pd.read_csv(TABLES / "fine_reference_profiles.csv")
    score = pd.read_csv(TABLES / "fine_reference_scores.csv")
    FIGS.mkdir(parents=True, exist_ok=True)
    gap_and_profiles(frame, prof, score)
    print(errors(score).to_string())


if __name__ == "__main__":
    main()
