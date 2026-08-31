#!/usr/bin/env python3
"""The D4-D7 window -- largest RMSE gain of any window, and why that is misleading.

D4-D7 returns M = 70 with a 3.09 m gain, the largest of the eight windows
tested on 2026-08-30. This figure exists to show that a larger gain is not a
better fit: it decomposes the gain per domain, and the decomposition says the
groin improves the domain OUTSIDE the dipole and degrades the downdrift domain
the structure actually acts on.

Writes output/groin_sweep/figures/fig_d4d7_window.png
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

SH = list(range(1, 13)); FIT = list(range(4, 8))     # D4-D7
fi = [SH.index(k) for k in FIT]
OUT = BASE/"output"/"groin_sweep"/"figures"
INK, MUTED, GRID = "#1a1a2e", "#5c6068", "#d5d8dd"
GOOD, BAD, ACC = "#2a9d5c", "#c2185b", "#c2185b"

d = pd.DataFrame([json.loads(l) for l in
                  open(BASE/"output"/"groin_sweep"/"1984_2004_edgeBE"/"sweep_results.jsonl")
                  if l.strip()])
d = d[(d.be1 == -42.6) & d.differential_err.notna()]
# LANDWARD-POSITIVE from here on, so erosion is UP and panel (a) reads as a
# plan view, matching the gifs. rate_D* and observed_change_profile are
# both SEAWARD-positive at source, so both are negated -- at the PLOTTING
# layer only. The scoring pipeline (FLIP_SIGN_MODEL, the sweep worker) is
# untouched, and panel (b) is unchanged either way because it plots
# |residual|.
prof = lambda r: -np.array([r[f"rate_D{k}"] for k in SH]) * 20.0
obs = -np.array([observed_change_profile(1984, 2004, SH)[k] for k in SH])
c = lambda v: np.asarray(v, float) - np.asarray(v, float)[fi].mean()
nog = prof(d[d.M == 0].iloc[0])
best = d[(d.M == 70) & (d.fraction == 1.0)].iloc[0]

plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.7,
                     "legend.frameon": False, "pdf.fonttype": 42})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14.5, 5.4),
                             gridspec_kw={"width_ratios": [1.35, 1]})

# ---- (a) the profile -------------------------------------------------------
a1.axvspan(3.5, 7.5, color="#fdf6d8", zorder=0)
a1.annotate("FIT WINDOW  D4−D7", xy=(5.5, 0.965), xycoords=("data", "axes fraction"),
            ha="center", va="top", fontsize=9.5, color="#8a7b1f",
            fontweight="bold")
a1.axvline(5.5, color=ACC, lw=1.4, ls=(0, (4, 2)), zorder=2)
a1.annotate("Buxton groin\nD6 updrift / D5 downdrift", xy=(5.56, 0.03),
            xycoords=("data", "axes fraction"), rotation=90, ha="left",
            va="bottom", fontsize=8.5, color=ACC)
a1.plot(SH, c(obs), "s--", ms=7, lw=2.2, color=INK, label="observed 1984→2004", zorder=6)
a1.plot(SH, c(nog), ":", lw=1.8, color=MUTED, label="no groin", zorder=4)
a1.plot(SH, c(prof(best)), "-", lw=2.0, color=ACC,
        label="M=70, f=1.0  (best on D4−D7)", zorder=5)
a1.set_xticks(SH); a1.set_xlabel("GIS domain")
a1.set_ylabel("Shoreline change 1984→2004 (m)\ndemeaned over D4−D7   [+ = landward, erosion]")
a1.grid(axis="y", color=GRID, lw=0.6); a1.set_axisbelow(True)
for s in ("top", "right"): a1.spines[s].set_visible(False)
a1.legend(fontsize=9, loc="lower left")
a1.set_title("(a)  D4−D7 gives the largest gain of any window (3.09 m)",
             loc="left", fontsize=11.5, color=INK, pad=8)

# ---- (b) where the gain comes from ----------------------------------------
rn = np.abs(c(nog) - c(obs))[fi]
rg = np.abs(c(prof(best)) - c(obs))[fi]
x = np.arange(len(FIT)); w = 0.38
a2.bar(x - w/2, rn, w, color=MUTED, label="no groin")
a2.bar(x + w/2, rg, w, color=ACC, label="M=70, f=1.0")
for i, (n, g) in enumerate(zip(rn, rg)):
    better = g < n
    a2.annotate(f"{n-g:+.1f}", xy=(i, max(n, g) + 0.6), ha="center", fontsize=10,
                fontweight="bold", color=GOOD if better else BAD)
a2.set_xticks(x)
a2.set_xticklabels([f"D{k}" + ("\ndowndrift" if k == 5 else
                               "\nupdrift" if k == 6 else "") for k in FIT])
a2.set_ylabel("|residual| (m)")
a2.grid(axis="y", color=GRID, lw=0.6); a2.set_axisbelow(True)
for s in ("top", "right"): a2.spines[s].set_visible(False)
a2.legend(fontsize=9)
a2.set_title("(b)  …but it degrades the domain the groin acts on",
             loc="left", fontsize=11.5, color=INK, pad=8)

fig.text(0.5, -0.06,
         "The gain is carried by D4, OUTSIDE the dipole. D5 — the downdrift domain — "
         "gets nearly 3× worse (3.4 → 9.2 m); D7 also worsens. RMSE squares residuals, "
         "so D4's large\nimprovement outweighs the degradations and the window scores well. "
         "This is the volume-neutral dipole failing as GROIN_PLAN.md predicts: observed "
         "downdrift extent is 0 m, the model's is 2,500 m.",
         ha="center", fontsize=8.8, color=MUTED)
OUT.mkdir(parents=True, exist_ok=True)
p = OUT/"fig_d4d7_window.png"
fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
print(f"  {p}")
