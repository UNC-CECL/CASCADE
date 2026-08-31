#!/usr/bin/env python3
"""Animated D2-D12 shoreline for a selection of (M, f) cells, for eyeballing the fit.

The existing run GIFs cover D1-D15 and exist only for the pair that was run as a
hindcast. Every SWEEP cell carries its own shoreline_matrix.npy, so any (M, f)
can be animated -- which is what is needed to judge whether M = 60 / f = 0.6 is
actually the best pair rather than just the best-scoring one.

Each frame shows shoreline CHANGE since 1984, demeaned over D4-D8 so the
alongshore shape is what the eye compares (a uniform level offset belongs to the
source/sink term, not the groin). The no-groin cell at the same be1 is drawn on
every frame as the reference, and the observed 1984->2004 change is drawn as a
fixed target so the endpoint can be judged.

    THE OBSERVED TARGET IS FIXED, AND THAT IS A LIMIT OF THIS WINDOW, NOT A
    CHOICE. Inside 1984-2004 the observation IS a single endpoint (an OLS fit
    to CoastSat chainage evaluated at both ends). The full-life companion,
    HAT_groin_full_life_gif.py, animates the observations too, because the
    1967 wet/dry record carries 19 dated surveys inside its window.

ORIENTATION: EROSION IS UP, SO THE PANEL READS AS A PLAN VIEW
    Changed 2026-08-30. These gifs were SEAWARD-positive and the full-life gif
    was LANDWARD-positive, so two animations in the same folder had opposite y
    axes. They are all landward-positive now: a retreating shoreline moves UP
    and the reader is looking down on the island with the ocean below the axis.

    Sign handling, which differs per source and is the easy thing to get wrong:
      `shoreline_matrix.npy`        Barrier3D's x_s, ALREADY landward-positive
          and already in METRES. Used as-is -- the negation that used to be
          here is gone.
      `observed_change_profile()`   chainage, SEAWARD-positive (verified in
          HAT_fullperiod_target.py against the published CoastSat LRR). It is
          NEGATED here.
    So the flip did not simply move a minus sign; it moved it from one series
    to the other.

DOMAINS D2-D12, matching the extent of the 1967 wet/dry survey and centred on
the groin at D5/D6. D1 is excluded: the cape's change over period 1 is 81-104 m,
about five times the groin's signal, and it swamps the axis.

Writes output/groin_sweep/figures/zoom_gifs_D2_D12/
"""
from __future__ import annotations
import sys
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
import numpy as np

_H = Path(__file__).resolve(); BASE = _H.parents[3]
for p in (BASE/"scripts", _H.parent):
    if str(p) not in sys.path: sys.path.insert(0, str(p))
from HAT_fullperiod_target import observed_change_profile

SWEEP = BASE/"output"/"groin_sweep"/"1984_2004_edgeBE"
OUT = BASE/"output"/"groin_sweep"/"figures"/"zoom_gifs_D2_D12"
BE, BUF = "-42.6", 15                 # real GIS n -> padded index BUF + (n-1)
DOM = list(range(2, 13)); FIT = list(range(4, 9))
PAD = [BUF + (n - 1) for n in DOM]
fi = [DOM.index(k) for k in FIT]
START = 1984
HOLD_FRAMES = 5                       # frames held on the last year
# NO dam->m CONVERSION. shoreline_matrix.npy is written in METRES already.
# Multiplying by 10 on 2026-08-30 made every curve ten times too large and put
# it off a +/-70 m axis. Checked against the cell's own shoreline_change_rate.csv:
# (m[-1] - m[0]) gives 82.9 m at D4 where the CSV reports 84.1, agreeing to
# the endpoint-vs-LRR estimator difference. Anything that rescales this must be
# re-checked against that CSV.

# Okabe-Ito, colour-vision-safe and muted enough to print. Shared with
# HAT_groin_full_life_gif.py so the two animations read as one set.
MODEL_C = "#D55E00"              # vermillion -- the swept cell
BASE_C = "#8C8C8C"               # grey       -- the no-groin reference
OBS_C = "#111111"                # near-black -- the observed endpoint
MARK_C = "#0072B2"               # blue       -- the structure
INK, SUBTLE, GRID = "#1A1A1A", "#6A6F76", "#E3E6EA"

# The pairs worth comparing: the chosen value, its f-neighbours, and M values
# either side -- plus a high-M cell where the module actually draws a dipole.
CELLS = [(60, 0.60), (60, 1.00), (60, 0.00),
         (40, 0.60), (95, 0.60), (160, 0.60)]

# What each cell is in the folder for. Shown as the subtitle so a gif opened on
# its own still says why it exists.
WHY = {
    (60, 0.60): "the selected pair",
    (60, 1.00): "f = 1.0 — the structure never deteriorates",
    (60, 0.00): "f = 0 — no trapping at all after 2003",
    (40, 0.60): "M too low — the fillet barely forms",
    (95, 0.60): "M too high — ~2x the littoral drift",
    (160, 0.60): "what the railed joint fit would have run",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10.5,
    "axes.linewidth": 0.8, "axes.edgecolor": "#4A4F55",
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 4.0, "ytick.major.size": 4.0,
    "xtick.color": "#4A4F55", "ytick.color": "#4A4F55",
    "text.color": INK, "axes.labelcolor": INK,
    "legend.frameon": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})


def load(M, f):
    """Change since 1984 per domain, LANDWARD-positive, metres."""
    p = SWEEP/f"M{M:g}_be{BE}_f{f:.2f}"/"shoreline_matrix.npy"
    if not p.is_file(): return None
    m = np.load(p)[:, PAD]                  # metres, landward-+, both already
    return m - m[0]


# observed_change_profile is SEAWARD-positive; negate it onto the landward-
# positive axis these panels now use.
obs = -np.array([observed_change_profile(1984, 2004, DOM)[k] for k in DOM])
c = lambda v: np.asarray(v, float) - np.asarray(v, float)[fi].mean()
nog = load(0, 0.0)
if nog is None:
    nog = load(0, 0.00)
OUT.mkdir(parents=True, exist_ok=True)

made = []
for M, f in CELLS:
    ch = load(M, f)
    if ch is None:
        print(f"  [skip] M={M:g} f={f:g} -- cell absent"); continue
    nyr = ch.shape[0]
    fig = plt.figure(figsize=(10.2, 6.4))
    ax = fig.add_axes([0.135, 0.175, 0.845, 0.615])

    fig.text(0.135, 0.940, f"M = {M:g} m/yr,  f = {f:g}",
             fontsize=16.5, fontweight="bold", color=INK, ha="left",
             va="center")
    fig.text(0.135, 0.898,
             f"{WHY.get((M, f), '')}   ·   period 1, 1984–2004   ·   "
             f"be1 = {BE}",
             fontsize=10.5, color=SUBTLE, ha="left", va="center")

    fig.legend(handles=[
        Line2D([], [], color=OBS_C, marker="s", markersize=6, linestyle="--",
               linewidth=1.8, label="Observed 1984→2004 (endpoint)"),
        Line2D([], [], color=MODEL_C, marker="o", markersize=5.5,
               linewidth=2.2, label=f"Model, M = {M:g}, f = {f:g}"),
        Line2D([], [], color=BASE_C, linestyle=":", linewidth=2.0,
               label="Model, no groin (same year)"),
    ], loc="upper left", bbox_to_anchor=(0.132, 0.862), ncol=3,
        fontsize=9.5, handlelength=2.4, columnspacing=2.6)

    fig.text(
        0.135, 0.030,
        "Landward-positive: erosion moves UP, so the panel reads as a plan view "
        "with the ocean below and the island above. Demeaned over D4–D8 — a "
        "uniform alongshore offset belongs to the source/sink calibration, not "
        "the groin, so only the SHAPE is being compared here. The observed "
        "target is a fixed endpoint because that is all this window has; the "
        "full-life gif animates the surveys.",
        fontsize=7.9, color=SUBTLE, ha="left", va="bottom", wrap=True)

    def frame(index, M=M, f=f, ch=ch, nyr=nyr, ax=ax):
        t = min(index, nyr - 1)
        ax.clear()
        ax.axvspan(3.5, 8.5, color="#FBF6E4", zorder=0)
        ax.annotate("fit window D4–D8", xy=(6.0, 0.015),
                    xycoords=("data", "axes fraction"), ha="center",
                    va="bottom", fontsize=8.2, color="#93842A")
        ax.axhline(0.0, color="#B9BEC4", linewidth=0.9, zorder=1)
        ax.axvline(5.5, color=MARK_C, lw=1.3, ls=(0, (5, 3)), zorder=2,
                   alpha=0.85)
        ax.annotate("Buxton groin field   D6 updrift | D5 downdrift",
                    xy=(5.5, 1.012), xycoords=("data", "axes fraction"),
                    ha="center", va="bottom", fontsize=8.2, color=MARK_C)

        ax.plot(DOM, c(obs), "s--", ms=6, lw=1.8, color=OBS_C, zorder=6)
        if nog is not None:
            ax.plot(DOM, c(nog[t]), ":", lw=2.0, color=BASE_C, zorder=4)
        ax.plot(DOM, c(ch[t]), "-o", ms=5.5, lw=2.2, color=MODEL_C, zorder=5)

        # -98 not -70: M = 160 dives past -75 at D12 and ran into the
        # fillet readout. All six cells share one limit so they stay
        # directly comparable.
        ax.set_xlim(1.6, 12.4); ax.set_ylim(-98, 70)
        ax.set_xticks(DOM)
        ax.set_xlabel("Alongshore GIS domain      "
                      "(D2 = toward Cape Point   →   D12 = north)", labelpad=7)
        ax.set_ylabel("Shoreline change since 1984 (m)\ndemeaned over D4−D8",
                      labelpad=26)
        ax.grid(axis="y", color=GRID, lw=0.7); ax.set_axisbelow(True)
        for s in ("top", "right"): ax.spines[s].set_visible(False)

        # Orientation cues, so "up = erosion" needs no caption to decode.
        ax.annotate("erosion\nlandward", xy=(-0.105, 0.94),
                    xycoords="axes fraction", ha="center", va="top",
                    fontsize=8.4, color=SUBTLE)
        ax.annotate("▲", xy=(-0.105, 0.955), xycoords="axes fraction",
                    ha="center", va="bottom", fontsize=7.5, color=SUBTLE)
        ax.annotate("accretion\nseaward", xy=(-0.105, 0.06),
                    xycoords="axes fraction", ha="center", va="bottom",
                    fontsize=8.4, color=SUBTLE)
        ax.annotate("▼", xy=(-0.105, 0.045), xycoords="axes fraction",
                    ha="center", va="top", fontsize=7.5, color=SUBTLE)

        ax.annotate(str(START + t), xy=(0.988, 0.965),
                    xycoords="axes fraction", ha="right", va="top",
                    fontsize=31, color="#DADEE3", fontweight="bold", zorder=0)

        # The fillet, model against the fixed observed endpoint.
        up_i, down_i = DOM.index(6), DOM.index(5)
        model_fillet = ch[t][down_i] - ch[t][up_i]
        observed_fillet = obs[down_i] - obs[up_i]
        ax.annotate(f"fillet  D5 − D6        model  {model_fillet:+.0f} m"
                    f"        observed  {observed_fillet:+.0f} m",
                    xy=(0.988, 0.035), xycoords="axes fraction", ha="right",
                    va="bottom", fontsize=9.4, color=INK)

    anim = FuncAnimation(fig, frame, frames=nyr + HOLD_FRAMES, interval=420)
    p = OUT/f"shoreline_D2-D12_M{M:g}_f{f:.2f}.gif"
    anim.save(p, writer=PillowWriter(fps=2.4))
    plt.close(fig); made.append(p.name); print(f"  {p.name}")

print(f"\n{len(made)} gifs -> {OUT}")
print("  y axis is LANDWARD-positive: erosion up, matching the full-life gif")
