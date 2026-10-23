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

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / "2-straight-coast" / "1-run"), str(HERE.parent / "3-groin-field")]
import straight_coast_test as sc  # noqa: E402
import groin_field_test as gf  # noqa: E402
import fine_grid_reference as fr  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
TABLES = workspace_twin(HERE) / "tables"
FIGS = workspace_twin(HERE) / "figures"
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
            s_ = f[(f.source == m) & (f.groin_length_m.isna() | (f.groin_length_m == fr.GROIN_LENGTH_M))]
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
    png = workspace_twin(FIGS) / "fine_reference_1_gap_and_profiles.png"
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
        "none levels off at a length. At -20 degrees the wave-driven groin "
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
    png = workspace_twin(FIGS) / "fine_reference_2_errors.png"
    save(fig, png, dpi=250, close=True)
    mean = score.groupby("module")[["gap_rmse_m", "profile_rmse_m"]].mean().round(0)
    record_caption(png, (
        "**Tests:** across shoreline orientations from -30 to +30 degrees, how far each version of the groin module, "
        f"on BRIE's 500 m grid, sits from the 25 m reference field with {fr.GROIN_LENGTH_M:g} m groins. **How to "
        "read:** (a) the root-mean-square difference in the gap between the updrift and downdrift domains over "
        f"years 1 to 50; (b) the same for the shoreline after 50 years, over the {2 * fr.WINDOW + 2} domains "
        "around the groin. Lower is closer. **Shows:** averaged over the seven orientations, the gap error is "
        f"{mean.loc['corrected', 'gap_rmse_m']:.0f} m for the wave-driven groin, {mean.loc['drift', 'gap_rmse_m']:.0f} m "
        f"drift, {mean.loc['dipole', 'gap_rmse_m']:.0f} m dipole and {mean.loc['pinned', 'gap_rmse_m']:.0f} m pinned; "
        f"the shoreline error {mean.loc['corrected', 'profile_rmse_m']:.0f}, {mean.loc['drift', 'profile_rmse_m']:.0f}, "
        f"{mean.loc['dipole', 'profile_rmse_m']:.0f} and {mean.loc['pinned', 'profile_rmse_m']:.0f} m. The corrected "
        "version is closest on coasts turned against the drift and drifts away from the reference on coasts turned "
        "with it, where the reference still fills to its tips and the module traps little "
        "(tables/fine_reference_scores.csv)."))
    return mean


# Figure 3: the gap at the last year against groin length, the reference and module 7
def length_check(lengths):
    fig, ax = plt.subplots(figsize=figsize("single", height=3.0), constrained_layout=True)
    ths = sorted(lengths.theta0_deg.unique())
    ramp = plt.get_cmap("Greys")(np.linspace(0.35, 0.95, len(ths)))
    for c, th in zip(ramp, ths):
        l_ = lengths[lengths.theta0_deg == th].sort_values("groin_length_m")
        ax.plot(l_.groin_length_m, l_.ref_gap_end_m, color=c, lw=1.4, marker="o", ms=3)
        ax.plot(l_.groin_length_m, l_.module_gap_end_m, color=c, lw=1.0, ls=(0, (3, 2)),
                marker=gf.MARKERS["length"], ms=3, mfc=sc.MODULES["length"]["col"], mec=sc.MODULES["length"]["col"])
    x = np.array([0, lengths.groin_length_m.max()])
    ax.plot(x, 2 * x, color=C["GROIN"], lw=0.8, ls=(0, (1, 1.5)))
    ax.text(20, 75, "gap = 2L", color=C["GROIN"], fontsize=6.5, ha="left", va="bottom", rotation=0)
    ax.set_xlabel("Groin length L (m)")
    ax.set_ylabel("Gap after 50 yr (m)")
    ax.set_xlim(0, x[1] * 1.05)
    ax.set_ylim(0, None)
    ax.grid(axis="y")
    open_frame(ax)
    handles = [Line2D([], [], color=INK, lw=1.4, marker="o", ms=3, label="reference, 25 m cells"),
               Line2D([], [], color=INK, lw=1.0, ls=(0, (3, 2)), marker=gf.MARKERS["length"], ms=3,
                      mfc=sc.MODULES["length"]["col"], mec=sc.MODULES["length"]["col"], label=sc.MODULES["length"]["label"]),
               Line2D([], [], color=ramp[0], lw=1.4, label=f"{ths[0]:+d}$^\\circ$ (light) to {ths[-1]:+d}$^\\circ$ (dark)")]
    ax.legend(handles=handles, loc="upper left", fontsize=6.3, frameon=False)
    fig.suptitle("Does the gap scale with groin length?")
    png = workspace_twin(FIGS) / "fine_reference_3_length.png"
    save(fig, png, dpi=250, close=True)
    sub = lengths[lengths.theta0_deg <= 10]
    record_caption(png, (
        "**Tests:** whether the wave-driven groin with a groin length (module 7) makes the gap across the groin "
        "scale with the groin's length the way a real-sized field does. **How to read:** the gap between the "
        "updrift and downdrift domains after 50 years against the groin length L, for the 25 m reference (solid, "
        "circles) and module 7 on BRIE's 500 m grid with the same L (dashed, pink triangles), one shade per "
        "orientation from -30 (light) to +30 degrees (dark); the dotted red line is a gap of twice the groin "
        "length. **Shows:** the reference follows gap = 2L at -30 to +20 degrees "
        f"({sub[sub.groin_length_m == fr.LENGTH_SENSITIVITY_M[0]].ref_gap_end_m.min():.0f} to "
        f"{sub[sub.groin_length_m == fr.LENGTH_SENSITIVITY_M[1]].ref_gap_end_m.max():.0f} m). Module 7 rises far "
        "more slowly: from "
        f"{sub[sub.groin_length_m == fr.LENGTH_SENSITIVITY_M[0]].module_gap_end_m.min():.0f}-"
        f"{sub[sub.groin_length_m == fr.LENGTH_SENSITIVITY_M[0]].module_gap_end_m.max():.0f} m at "
        f"{fr.LENGTH_SENSITIVITY_M[0]:g} m to "
        f"{sub[sub.groin_length_m == fr.LENGTH_SENSITIVITY_M[1]].module_gap_end_m.min():.0f}-"
        f"{sub[sub.groin_length_m == fr.LENGTH_SENSITIVITY_M[1]].module_gap_end_m.max():.0f} m at "
        f"{fr.LENGTH_SENSITIVITY_M[1]:g} m at -30 to +10 degrees, and near zero at +20 and +30. Its blocking reads "
        "the build-out of the whole 500 m updrift cell, which BRIE's diffusion keeps well short of the tip, so "
        "the groin starts passing sand long before a real one would (tables/fine_reference_length_check.csv)."))


# Run: read the tables, draw the figures
def main():
    apply_style()
    frame = pd.read_csv(workspace_twin(TABLES) / "fine_reference_by_year.csv")
    prof = pd.read_csv(workspace_twin(TABLES) / "fine_reference_profiles.csv")
    score = pd.read_csv(workspace_twin(TABLES) / "fine_reference_scores.csv")
    FIGS.mkdir(parents=True, exist_ok=True)
    gap_and_profiles(frame, prof, score)
    print(errors(score).to_string())
    length_check(pd.read_csv(workspace_twin(TABLES) / "fine_reference_length_check.csv"))


if __name__ == "__main__":
    main()
