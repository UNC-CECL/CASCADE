#!/usr/bin/env python3
"""The D4-D7 window -- largest RMSE gain of any window, and why that is misleading.

D4-D7 returns M = 70 with a 3.09 m gain, the largest of the eight windows
tested on 2026-08-30. This figure exists to show that a larger gain is not a
better fit: it decomposes the gain per domain, and the decomposition says the
groin improves the domain OUTSIDE the dipole and degrades the downdrift domain
the structure actually acts on.

STYLE, 2026-09-11
    Under the house style (`scripts/hat_figure_style.py`). The footnote
    paragraph is now the caption in output/groin_sweep/figures/CAPTIONS.md, the
    canvas is a 190 mm printed column instead of 14.5 in, and the bars are the
    house BASE grey for the baseline against ACCENT purple for the run under
    test. The per-domain deltas keep a green/red split, which is the one place
    in this figure set where those two mean better and worse rather than the
    1984/1997 vintages -- nothing on this figure is a vintage, and the whole
    point of panel (b) is which domains moved the wrong way.

Writes output/groin_sweep/figures/fig_d4d7_window.png (and .pdf)
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
from hat_figure_style import (apply_style, C, C_1984, INK, INK_MUTED, caption,
                              figsize, open_frame, save, _title)

SH = list(range(1, 13)); FIT = list(range(4, 8))     # D4-D7
fi = [SH.index(k) for k in FIT]
OUT = BASE/"output"/"groin_sweep"/"figures"
# BETTER / WORSE on panel (b) only; ACC is the run under test, FOIL the
# baseline it is scored against.
BETTER, WORSE = C["REF"], C_1984
ACC, FOIL, BAND = C["ACCENT"], C["BASE"], "0.94"

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

apply_style()
fig, (a1, a2) = plt.subplots(1, 2, figsize=figsize("double", aspect=0.42),
                             gridspec_kw={"width_ratios": [1.35, 1]},
                             constrained_layout=True)

# ---- (a) the profile -------------------------------------------------------
a1.axvspan(3.5, 7.5, color=BAND, zorder=0)
a1.annotate("fit window, D4−D7", xy=(5.5, 0.965), xycoords=("data", "axes fraction"),
            ha="center", va="top", fontsize=7.5, color=INK_MUTED)
a1.axvline(5.5, color=INK_MUTED, lw=0.8, ls=(0, (4, 2)), zorder=2)
a1.annotate("Buxton groin: D6 updrift, D5 downdrift", xy=(5.56, 0.03),
            xycoords=("data", "axes fraction"), rotation=90, ha="left",
            va="bottom", fontsize=7, color=INK_MUTED)
a1.plot(SH, c(obs), "s--", ms=4.0, lw=1.8, color=INK,
        label="observed, 1984 to 2004", zorder=6)
a1.plot(SH, c(nog), ":", lw=1.4, color=FOIL, label="no groin", zorder=4)
a1.plot(SH, c(prof(best)), "-", lw=1.6, color=ACC,
        label="M 70, f 1.0, the best cell on D4−D7", zorder=5)
a1.set_xticks(SH); a1.set_xlabel("GIS domain")
a1.set_ylabel("shoreline change 1984 to 2004 (m)\ndemeaned over D4−D7;"
              " positive is landward")
a1.grid(axis="y"); a1.set_axisbelow(True)
open_frame(a1)
_title(a1, 0, "the window with the largest gain")

# ---- (b) where the gain comes from ----------------------------------------
rn = np.abs(c(nog) - c(obs))[fi]
rg = np.abs(c(prof(best)) - c(obs))[fi]
x = np.arange(len(FIT)); w = 0.38
a2.bar(x - w/2, rn, w, color=FOIL, label="no groin")
a2.bar(x + w/2, rg, w, color=ACC, label="M 70, f 1.0")
for i, (n, g) in enumerate(zip(rn, rg)):
    a2.annotate(f"{n-g:+.1f}", xy=(i, max(n, g) + 0.6), ha="center",
                fontsize=8, color=BETTER if g < n else WORSE)
a2.set_xticks(x)
a2.set_xticklabels([f"D{k}" + ("\ndowndrift" if k == 5 else
                               "\nupdrift" if k == 6 else "") for k in FIT],
                   fontsize=7.5)
a2.set_ylabel("|residual| (m)")
a2.grid(axis="y"); a2.set_axisbelow(True)
open_frame(a2)
_title(a2, 1, "the gain, domain by domain")

# Panel (a)'s handles only: panel (b)'s bars are the same two things in the
# same two colours, and a figure-wide legend listed each of them twice.
fig.legend(*a1.get_legend_handles_labels(), loc="outside lower center",
           ncol=3, frameon=False)

caption(fig,
        "D4−D7 returns M = 70 with a 3.09 m error gain, the largest of the "
        "eight windows tested on 2026-08-30, and this figure is why a larger "
        "gain is not a better fit. (a) The observed profile against the best "
        "cell on this window and the no-groin baseline, demeaned over D4−D7 "
        "and drawn landward-positive so erosion is up and the panel reads as a "
        "plan view; both sources are seaward-positive and are negated at the "
        "plotting layer only, which leaves the scoring pipeline untouched. "
        "(b) The same fit decomposed per domain, as the absolute residual with "
        "the groin off and on, with the change annotated: green where the "
        "groin improves the domain, red where it degrades it. The gain is "
        "carried by D4, OUTSIDE the dipole. D5, the downdrift domain the "
        "structure actually acts on, gets nearly three times worse, from 3.4 "
        "to 9.2 m, and D7 also worsens. Because RMSE squares residuals, D4's "
        "large improvement outweighs the two degradations and the window "
        "scores well. This is the volume-neutral dipole failing exactly as "
        "GROIN_PLAN.md predicts: the observed downdrift extent is 0 m and the "
        "model's is 2,500 m. Green and red mean better and worse on panel (b) "
        "only; elsewhere in this figure set they are the 1984 and 1997 "
        "vintages, and no vintage is drawn here.")

OUT.mkdir(parents=True, exist_ok=True)
for path in save(fig, OUT/"fig_d4d7_window.png", close=True):
    print(f"  {path}")
