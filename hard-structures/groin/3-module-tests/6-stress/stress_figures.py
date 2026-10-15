#!/usr/bin/env python3
"""
The stress-test figures, drawn from the tables stress_test.py writes.

    python stress_figures.py   ->  figures/stress_<n>_<what>.png

Captions go to figures/supporting/CAPTIONS.md. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import SymLogNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stress_test as st  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
TABLES = HERE / "tables"
FIGS = HERE / "figures"
CMAP = "RdBu"                     # diverging: sand lost red, gained blue, a grey-white midpoint
LIMIT_M = 10_000                  # colour scale ends; larger values sit at the end colour
WAVE_ROWS = [f"{p}={v}" for p, vals in st.WAVE_SWEEP.items() for v in vals]
FAIL_MARK = {"crash": "X", "conserve": "c", "angle": "a", "steady": "s"}
MODE_CASES = (("waves", ("ahf=0.1", "ahf=0.3", "ahf=0.5")), ("waves", ("asym=0.6", "asym=0.7", "asym=0.8")))
MODE_THETA, MODE_S = -20, 0.25
MODE_COLS = (C["ADDED"], C["GROIN"])           # option A is grey; the stronger stresses gold then red
OPTION_A_CASES = ("ahf=0.5", "asym=0.6")
# -----------------------------------------------------------------------------


# The first rule a run breaks, as one letter; blank if it passes
def fail_letter(r):
    if r.passes:
        return ""
    if pd.notna(r.crashed_year):
        return FAIL_MARK["crash"]
    if not r.conserves:
        return FAIL_MARK["conserve"]
    if not r.angle_ok:
        return FAIL_MARK["angle"]
    return FAIL_MARK["steady"]


# One heat panel: rows x orientations of the groin's net sand against the no-groin run, failures lettered
def heat(ax, d, rows, norm):
    cols = list(st.ORIENTATIONS)
    grid = np.full((len(rows), len(cols)), np.nan)
    marks = {}
    for _, r in d.iterrows():
        if r.row in rows and r.theta0 in cols:
            i, j = rows.index(r.row), cols.index(r.theta0)
            grid[i, j] = r.net_vs_nogroin_m
            marks[(i, j)] = fail_letter(r)
    im = ax.imshow(np.clip(grid, -LIMIT_M, LIMIT_M), cmap=CMAP, norm=norm, aspect="auto")
    for (i, j), m in marks.items():
        if m:
            dark = abs(norm(np.clip(grid[i, j], -LIMIT_M, LIMIT_M)) - 0.5) > 0.3   # deep red or blue cell
            ax.text(j, i, m, ha="center", va="center", fontsize=6.5, fontweight="bold",
                    color="white" if dark else INK)
    ax.set_xticks(range(len(cols)), [f"{c:+d}" for c in cols])
    ax.set_yticks(range(len(rows)), rows, fontsize=6.5)
    ax.set_xlabel("Shoreline orientation ($^\\circ$)")
    ax.tick_params(length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    return im


# Figure 1: where it breaks, across wave climates
def waves_figure(t):
    d = t[t.family == "waves"].assign(row=lambda x: x.case)
    norm = SymLogNorm(linthresh=st.CONSERVE_ABS_M, vmin=-LIMIT_M, vmax=LIMIT_M, base=10)
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=5.2), constrained_layout=True, sharey=True)
    for k, (ax, s) in enumerate(zip(axes, st.STRENGTHS)):
        im = heat(ax, d[d.strength == s], WAVE_ROWS, norm)
        for row in ("ahf=0.5", "asym=0.6", "Hs=2.0"):           # option A, one row in each sweep
            ax.axhline(WAVE_ROWS.index(row) - 0.5, color=INK_MUTED, lw=0.4)
            ax.axhline(WAVE_ROWS.index(row) + 0.5, color=INK_MUTED, lw=0.4)
        for edge in (WAVE_ROWS.index("ahf=0.1") - 0.5, WAVE_ROWS.index("asym=0.5") - 0.5):
            ax.axhline(edge, color=INK, lw=1.0)
        _letter_inside(ax, k)
        ax.set_title(f"Strength {s:g}", loc="center")
    cb = fig.colorbar(im, ax=axes, shrink=0.7, pad=0.01)
    cb.set_label("Net sand the groin adds, against the no-groin run\n(m of summed shoreline; red lost, blue made)",
                 fontsize=7)
    fig.suptitle("Where does the runner's groin break, across wave climates?")
    png = FIGS / "stress_1_where_it_breaks.png"
    save(fig, png, dpi=250, close=True)
    opt = d[(d.case == "asym=0.6")]
    record_caption(png, (
        "**Tests:** whether the runner's groin (`cascade.groin.CorrectedGroinCallback`, unchanged, 10-cell "
        "correction zone) keeps the sand balance and stays stable when one wave parameter at a time is moved far "
        "from option A (Hs 2.0, Tp 7.5, asymmetry 0.6, high-angle fraction 0.5). **How to read:** BRIE's "
        "alongshore step alone on a straight coast, one groin, 100 years; each cell is one run, scored against the "
        "same run with the groin off. Colour: the sand the groin made (blue) or lost (red) by year 100, summed "
        f"over the coast, on a log scale that is linear inside ±{st.CONSERVE_ABS_M:g} m. Letters: the first rule "
        f"broken. X the run crashed (an angle past BRIE's diffusivity table); c sand not conserved (more than "
        f"{st.CONSERVE_ABS_M:g} m or {st.CONSERVE_FRAC:.0%} of the groin's gross change); a a face steeper than "
        f"{st.MAX_ANGLE_DEG:g}$^\\circ$; s the gap still moving over the last fifth of the run. Grey-ruled rows are "
        "option A. Buxton's groin cells start at about −12$^\\circ$ and −25$^\\circ$. **Shows:** at option A "
        f"the groin holds the balance at every orientation but +10$^\\circ$ "
        f"({opt[(opt.strength == 0.25) & (opt.theta0 == 10)].net_vs_nogroin_m.iloc[0]:+.0f} m at strength 0.25). "
        "It breaks two ways: where BRIE's diffusion is strong (high-angle fraction 0.3 or less, Hs 2.5 and more "
        "on coasts turned with the drift) the leak correction overshoots and runs away, and where BRIE's "
        "diffusivity is clipped to zero at the groin (asymmetry 0.7 and more on coasts turned against the "
        "drift) nothing carries sand back across the face, so the gap grows without limit."))


# A row label: the family's word, then the case's value without its own name
def row_label(labels, fam, case):
    value = case.split("=", 1)[1] if "=" in case else case
    if fam == "notch":
        value = f"{int(value):+d} m"
    return f"{labels[fam]} {value}".strip() if labels[fam] else case


# Figure 2: the other families, as the same scorecard
def families_figure(t):
    labels = {"strength": "strength", "notch": "start gap", "switch-off": "", "cycling": "",
              "background": "background", "placement": "groin at", "long": "", "variable": "waves vary:"}
    d = t[t.family.isin(labels)].copy()
    d["row"] = [row_label(labels, f, c) for f, c in zip(d.family, d.case)]
    v = d[d.family == "variable"].copy()
    v["row"] = "waves vary: " + v.case.str.split(", ").str[1]
    v = (v.groupby(["row", "theta0", "strength"], as_index=False)
          .agg(net_vs_nogroin_m=("net_vs_nogroin_m", "mean"), passes=("passes", "all"),
               crashed_year=("crashed_year", "max"), conserves=("conserves", "all"), angle_ok=("angle_ok", "all")))
    d = pd.concat([d[d.family != "variable"], v], ignore_index=True)
    order = []
    for fam in labels:
        sub = t[t.family == fam]
        if fam == "variable":
            order += ["waves vary: as written", "waves vary: refreshed"]
        else:
            order += list(dict.fromkeys(row_label(labels, fam, c) for c in sub.case))
    norm = SymLogNorm(linthresh=st.CONSERVE_ABS_M, vmin=-LIMIT_M, vmax=LIMIT_M, base=10)
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=6.4), constrained_layout=True, sharey=True)
    strength_rows = [r for r in order if r.startswith("strength")]
    for k, (ax, s) in enumerate(zip(axes, st.STRENGTHS)):
        sub = d[(d.strength == s) | d.row.isin(strength_rows)]
        im = heat(ax, sub, order, norm)
        _letter_inside(ax, k)
        ax.set_title(f"Strength {s:g} (top rows: strength 0 to 1)", loc="center")
    cb = fig.colorbar(im, ax=axes, shrink=0.6, pad=0.01)
    cb.set_label("Net sand the groin adds, against the no-groin run\n(m of summed shoreline; red lost, blue made)",
                 fontsize=7)
    fig.suptitle("What else breaks it, at option A waves?")
    png = FIGS / "stress_2_other_stresses.png"
    save(fig, png, dpi=250, close=True)
    p = t[t.family == "placement"]
    record_caption(png, (
        "**Tests:** at option A waves, whether the runner's groin keeps the sand balance and stays stable under "
        "the other stresses: its strength from 0 to 1; a sharp step already at the groin when it starts (the gap "
        "in metres, updrift seaward positive; Buxton's is 256 m); on for 50 years then off for 100; on and off "
        "every 5 years; cross-shore change passed in as Barrier3D passes it (2 m/yr retreat, ±10 m/yr noise per "
        "domain, a 30 m cut every 10 years); the groin face near the end of the array (lo is the lower domain of "
        "the face, of 61); 500 years; and a new wave climate every year (asymmetry 0.4-0.8, Hs 1.5-2.5, the "
        "callback as written and with its drift table refreshed each year; mean of three seeds, three "
        "orientations). **How to read:** as figure 1: colour is the net sand against the matched no-groin run, "
        "letters the first rule broken. **Shows:** strength 0 is the no-groin run exactly, and up to 0.5 the "
        "groin conserves at every orientation; the switch-off, the background change and 500 years leave it "
        "conserving at strength 0.25 except at +10$^\\circ$; it leaks when switched on and off every 5 years, "
        f"when its face is within {st.HALF_WIDTH} cells of the end of the array (up to "
        f"{p.net_vs_nogroin_m.abs().max():.0f} m: the correction zone is cut off there and does not wrap), and "
        "in almost every run where the waves change from year to year."))


# The year rows of the groin run at strength s: the year table holds both strengths back to back, each from year 1
def runs_of(rows, s):
    starts = list(np.flatnonzero(rows.year.to_numpy() == 1)) + [len(rows)]
    for a, b in zip(starts[:-1], starts[1:]):
        if abs(rows.strength.iloc[a] - s) < 1e-9:
            return rows.iloc[a:b]
    raise KeyError(f"no run at strength {s}")


# Figure 3: the two ways it breaks, year by year
def modes_figure(t, by):
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", height=5.0), constrained_layout=True, sharex=True)
    g = by[(by.family == "waves") & (by.theta0 == MODE_THETA)]
    for k, (fam, cases) in enumerate(MODE_CASES):
        stressed = iter(MODE_COLS)
        for case in sorted(cases, key=lambda c: (c not in OPTION_A_CASES, abs(float(c.split("=")[1]) - 0.55))):
            col = C["BASE"] if case in OPTION_A_CASES else next(stressed)
            label = f"{case} (option A)" if case in OPTION_A_CASES else case
            seg = runs_of(g[g.groin & (g.case == case)], MODE_S)
            n = g[~g.groin & (g.case == case)].drop_duplicates("year").set_index("year").reindex(seg.year)
            axes[0, k].plot(seg.year, seg.gap_m.to_numpy() - n.gap_m.to_numpy(), color=col, lw=1.4, label=label)
            axes[1, k].plot(seg.year, seg.net_m.to_numpy() - n.net_m.to_numpy(), color=col, lw=1.4)
    titles = ("Diffusion too strong: the correction overshoots", "Diffusion clipped to zero: the coast is frozen")
    for k in range(2):
        axes[0, k].set_title(titles[k], loc="center", fontsize=8)
        axes[0, k].legend(fontsize=6.5, frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.12))
        axes[1, k].set_xlabel("Model year")
        for i in range(2):
            ax = axes[i, k]
            ax.axhline(0, color=INK_MUTED, lw=0.5)
            ax.grid(axis="y")
            open_frame(ax)
            _letter_inside(ax, 2 * i + k)
    axes[0, 0].set_ylabel("Gap across the groin,\nagainst no groin (m)")
    axes[1, 0].set_ylabel("Net sand the groin adds\n(m of summed shoreline)")
    axes[1, 0].set_yscale("symlog", linthresh=10)
    axes[1, 1].set_yscale("symlog", linthresh=10)
    axes[0, 0].set_yscale("symlog", linthresh=100)
    fig.suptitle(f"The two ways it breaks, coast at {MODE_THETA:+d}$^\\circ$, strength {MODE_S:g}")
    png = FIGS / "stress_3_two_failure_modes.png"
    save(fig, png, dpi=250, close=True)
    record_caption(png, (
        "**Tests:** how the runner's groin fails when it does, year by year. **How to read:** BRIE's alongshore "
        f"step alone, a straight coast at {MODE_THETA:+d}$^\\circ$ (Buxton's groin cells are at about −12$^\\circ$ "
        f"and −25$^\\circ$), strength {MODE_S:g}, against the same run with no groin. Top: the gap across the groin "
        "(updrift seaward of downdrift); bottom: the sand the groin made or lost, summed over the coast; both on "
        "log scales that are linear near zero (the top right is linear). Grey is option A in each column. (a, c) Lower high-angle fraction "
        "makes BRIE's diffusion stronger; at 0.3 and 0.1 the leak correction, which is an explicit term in front "
        "of BRIE's implicit step, overshoots and grows each year; at 0.1 the run ends where an angle runs off "
        "BRIE's diffusivity table. (b, d) Higher asymmetry clips BRIE's diffusivity to zero at the groin's angles; "
        "nothing carries sand back across the face, the groin keeps moving its share of the drift every year, and "
        "the gap grows without limit. **Shows:** neither happens at option A; both happen inside the range of the "
        "groin-on wave sweep (high-angle 0.3, asymmetry 0.7 and 0.8)."))


def main():
    apply_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(TABLES / "stress_runs.csv")
    by = pd.read_csv(TABLES / "stress_by_year.csv.gz")
    waves_figure(t)
    families_figure(t)
    modes_figure(t, by)
    print("wrote", FIGS.relative_to(st.REPO))


if __name__ == "__main__":
    main()
