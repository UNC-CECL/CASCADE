#!/usr/bin/env python3
"""Every M value as its own panel, with all six f curves drawn on it.

The top-N overlay (fig_top_profiles.png) shows only the best cells, and they sit
so close together that nothing about the parameter response is visible. This
draws the whole grid instead: one panel per M, six f curves inside it, the same
observed and no-groin reference on every panel, and a shared y axis so panels
can be read against each other.

WHAT TO LOOK FOR
    * within a panel: how much f moves the profile at fixed M. Period 1 mostly
      PRECEDES the 1996-2003 deterioration ramp, so f should move it little --
      period-1 cumulative trapping is M(15.5 + 4.5f), which f changes by only
      29% across its whole range.
    * across panels: M lifts the whole curve rather than building a local
      fillet at D5/D6. That is the finding the per-domain decomposition in
      fig_d4d7_window.png makes numerically -- the groin's gain comes from D4,
      outside the dipole, while D5 (downdrift) gets worse.

Writes output/groin_sweep/figures/profiles_by_M/
    fig_all_M_profiles.png     the grid, for comparison across M
    fig_M<value>.png           one file per M, for detail
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

_H = Path(__file__).resolve(); BASE = _H.parents[3]
for p in (BASE/"scripts", _H.parent):
    if str(p) not in sys.path: sys.path.insert(0, str(p))
from HAT_fullperiod_target import observed_change_profile
from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED, caption,
                              figsize, open_frame, save, _title)

SH = list(range(1, 13)); FIT = list(range(4, 9))
fi = [SH.index(k) for k in FIT]
PINNED_BE1, YRS = -42.6, 20.0
OUT = BASE/"output"/"groin_sweep"/"figures"/"profiles_by_M"
apply_style()

d = pd.DataFrame([json.loads(l) for l in
                  open(BASE/"output"/"groin_sweep"/"1984_2004_edgeBE"/"sweep_results.jsonl")
                  if l.strip()])
d = d[(d.be1 == PINNED_BE1) & d.differential_err.notna()].copy()
obs = np.array([observed_change_profile(1984, 2004, SH)[k] for k in SH])
c = lambda v: np.asarray(v, float) - np.asarray(v, float)[fi].mean()
prof = lambda r: np.array([r[f"rate_D{k}"] for k in SH]) * YRS
obs_c = c(obs)
rmse = lambda v: float(np.sqrt(((c(v) - obs_c)[fi] ** 2).mean()))
d["rmse"] = d.apply(lambda r: rmse(prof(r)), axis=1)

nog_row = d[d.M == 0].iloc[0]
nog_c, nog_rmse = c(prof(nog_row)), rmse(prof(nog_row))
MS = sorted(m for m in d.M.unique() if m > 0)
FS = sorted(d.fraction.unique())
# The f family as one colour family, so a panel reads as "one M, six f" and
# not as six unrelated series; darkest is the largest f. viridis was used here
# until 2026-09-11 and shared its green with the reference marks elsewhere.
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
CMAP = LinearSegmentedColormap.from_list(
    "hat_accent_ramp", [C["ACCENT_FILL"], C["ACCENT"]])(
        np.linspace(0.15, 1.0, len(FS)))
# Top widened on the 2026-08-30 flip to landward-positive: the high-M cells
# (M >= 110) reach past +70 at D12 and were clipping.
YLIM = (-60, 88)


def draw(ax, M, compact, letter=None):
    ax.axvspan(3.5, 8.5, color="0.94", zorder=0)
    ax.axvline(5.5, color=INK_MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.plot(SH, nog_c, ":", lw=1.4, color=C["BASE"], zorder=3,
            label=f"no groin, {nog_rmse:.1f} m")
    sub = d[d.M == M].sort_values("fraction")
    for colour, (_, r) in zip(CMAP, sub.iterrows()):
        # The per-cell error is panel-specific: in the grid these handles feed
        # ONE shared legend, so quoting a number there would attach panel (a)'s
        # errors to every panel.
        ax.plot(SH, c(prof(r)), "-", lw=1.2, color=colour, zorder=4,
                label=(f"f {r.fraction:g}" if compact
                       else f"f {r.fraction:g}, {r.rmse:.1f} m"))
    ax.plot(SH, obs_c, "s--", ms=3.0, lw=1.6, color=INK, zorder=6,
            label="observed")
    ax.set_ylim(*YLIM); ax.set_xlim(0.5, 12.5)
    ax.set_xticks(SH if not compact else SH[::2])
    ax.grid(axis="y"); ax.set_axisbelow(True)
    open_frame(ax)
    best = sub.loc[sub.rmse.idxmin()]
    if letter is None:
        ax.set_title(f"M = {M:g}, best f {best.fraction:g} at "
                     f"{best.rmse:.2f} m", loc="left")
    else:
        # Four to a row: the letter and a full sentence collide, so the panel
        # carries the M and the caption carries the best f per panel.
        _title(ax, letter, f"M = {M:g}")
    return best


OUT.mkdir(parents=True, exist_ok=True)

# ---- the grid --------------------------------------------------------------
ncol = 4; nrow = int(np.ceil(len(MS) / ncol))
fig, axes = plt.subplots(nrow, ncol,
                         figsize=figsize("double", height=1.55 * nrow),
                         sharex=True, sharey=True, constrained_layout=True)
bests = {}
for i, (ax, M) in enumerate(zip(axes.flat, MS)):
    bests[M] = draw(ax, M, compact=True, letter=i)
for ax in axes.flat[len(MS):]:
    ax.axis("off")
for ax in axes[-1]: ax.set_xlabel("GIS domain")
# One shared y label: repeated on three rows it took a third of the canvas.
fig.supylabel("change 1984 to 2004 (m), demeaned over D4−D8;"
              " positive is landward", fontsize=9)
handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="outside lower center", ncol=len(labels),
           frameon=False, fontsize=7)
caption(fig,
        "Every M in the period-1 grid as its own panel, with all {nf} f curves "
        "on each, at be1 = {be:g}. The observed profile and the no-groin "
        "baseline are repeated on every panel and the y axis is shared, so "
        "panels can be read against one another; the shaded band is the D4−D8 "
        "fit window and the dashed vertical is the structure between D5 and "
        "D6. Two things to look for. Within a panel, f moves the profile far "
        "less than M does, because period 1 mostly PRECEDES the 1996 to 2003 "
        "deterioration ramp: period-1 cumulative trapping is M(15.5 + 4.5f), "
        "which f changes by only 29% across its whole range. Across panels, M "
        "shifts the ENTIRE curve rather than building a local fillet at D5 and "
        "D6 — which is the finding the per-domain decomposition in the D4−D7 "
        "window figure makes numerically, where the groin's gain comes from "
        "D4, outside the dipole, while the downdrift domain gets worse. "
        "Everything is demeaned over the fit window and drawn "
        "landward-positive, so erosion is up. The best f in each panel, with "
        "its error: {best}."
        .format(nf=len(FS), be=PINNED_BE1,
                best="; ".join(
                    "M {:g}, f {:g} at {:.2f} m".format(
                        m, b.fraction, b.rmse)
                    for m, b in bests.items())))
p = save(fig, OUT/"fig_all_M_profiles.png", close=True)[0]
print(f"  {p.name}")

# ---- one per M -------------------------------------------------------------
for M in MS:
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.44),
                           constrained_layout=True)
    draw(ax, M, compact=False)
    ax.set_xlabel("GIS domain")
    ax.set_ylabel("shoreline change 1984 to 2004 (m)\ndemeaned over D4−D8;"
                  " positive is landward")
    fig.legend(loc="outside lower center", ncol=4, frameon=False, fontsize=7)
    caption(fig,
            "One M from the period-1 grid, {M:g} m/yr, with all {nf} "
            "deterioration floors drawn on it against the observed profile and "
            "the no-groin baseline, at be1 = {be:g}. Darkest is the largest f. "
            "The shaded band is the D4−D8 fit window, the dashed vertical the "
            "structure between D5 and D6, and everything is demeaned over the "
            "fit window and drawn landward-positive so erosion is up. The "
            "companion grid figure puts every M side by side; f moves these "
            "curves little because period 1 mostly precedes the 1996 to 2003 "
            "deterioration ramp."
            .format(M=M, nf=len(FS), be=PINNED_BE1))
    save(fig, OUT/f"fig_M{M:g}.png", close=True)
print(f"  + {len(MS)} per-M files")
