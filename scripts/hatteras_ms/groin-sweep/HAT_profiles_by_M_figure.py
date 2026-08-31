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

SH = list(range(1, 13)); FIT = list(range(4, 9))
fi = [SH.index(k) for k in FIT]
PINNED_BE1, YRS = -42.6, 20.0
OUT = BASE/"output"/"groin_sweep"/"figures"/"profiles_by_M"
INK, MUTED, GRID = "#1a1a2e", "#5c6068", "#d5d8dd"

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.7,
                     "legend.frameon": False, "pdf.fonttype": 42})

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
CMAP = plt.cm.viridis(np.linspace(0.08, 0.92, len(FS)))
# Top widened on the 2026-08-30 flip to landward-positive: the high-M cells
# (M >= 110) reach past +70 at D12 and were clipping.
YLIM = (-60, 88)


def draw(ax, M, compact):
    ax.axvspan(3.5, 8.5, color="#fdf6d8", zorder=0)
    ax.axvline(5.5, color="#b0b4ba", lw=1.1, ls=(0, (3, 2)), zorder=1)
    ax.plot(SH, nog_c, ":", lw=1.6, color=MUTED, zorder=3,
            label=f"no groin ({nog_rmse:.1f} m)")
    sub = d[d.M == M].sort_values("fraction")
    for colour, (_, r) in zip(CMAP, sub.iterrows()):
        ax.plot(SH, c(prof(r)), "-", lw=1.5, color=colour, zorder=4,
                label=f"f={r.fraction:g}  ({r.rmse:.1f})")
    ax.plot(SH, obs_c, "s--", ms=4.5, lw=1.8, color=INK, zorder=6,
            label="observed")
    ax.set_ylim(*YLIM); ax.set_xlim(0.5, 12.5)
    ax.set_xticks(SH if not compact else SH[::2])
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    best = sub.loc[sub.rmse.idxmin()]
    ax.set_title(f"M = {M:g}      best f = {best.fraction:g}, "
                 f"RMSE {best.rmse:.2f} m", loc="left",
                 fontsize=10.5, color=INK, pad=5)


OUT.mkdir(parents=True, exist_ok=True)

# ---- the grid --------------------------------------------------------------
ncol = 4; nrow = int(np.ceil(len(MS) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.3*ncol, 3.5*nrow),
                         sharex=True, sharey=True)
for ax, M in zip(axes.flat, MS):
    draw(ax, M, compact=True)
for ax in axes.flat[len(MS):]:
    ax.axis("off")
axes.flat[0].legend(fontsize=8, loc="upper left", ncol=2)
for ax in axes[-1]: ax.set_xlabel("GIS domain")
for ax in axes[:, 0]:
    ax.set_ylabel("change 1984→2004 (m)\ndemeaned over D4−D8\n[+ = landward, erosion]")
fig.suptitle("Every M, every f — period 1, be1 = −42.6.  Shaded = fit window "
             "D4−D8;  dashed vertical = the groin at D5/D6",
             fontsize=13.5, fontweight="bold", color=INK, y=0.998)
fig.text(0.5, 0.005,
         "f moves the profile far less than M does, because period 1 mostly "
         "precedes the 1996−2003 deterioration ramp: period-1 cumulative "
         "trapping is M(15.5 + 4.5f), which f changes by only 29% across its "
         "whole range.\nNote what M actually does — it shifts the entire curve "
         "rather than building a local fillet at D5/D6.",
         ha="center", fontsize=9, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 0.975))
p = OUT/"fig_all_M_profiles.png"
fig.savefig(p, dpi=180, bbox_inches="tight", facecolor="white")
plt.close(fig); print(f"  {p.name}")

# ---- one per M -------------------------------------------------------------
for M in MS:
    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    draw(ax, M, compact=False)
    ax.set_xlabel("GIS domain")
    ax.set_ylabel("Shoreline change 1984→2004 (m)\ndemeaned over D4−D8   [+ = landward, erosion]")
    ax.legend(fontsize=9, loc="upper left", ncol=2)
    fig.savefig(OUT/f"fig_M{M:g}.png", dpi=180, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
print(f"  + {len(MS)} per-M files")
