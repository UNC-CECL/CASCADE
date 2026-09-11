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
from hat_figure_style import (apply_style, C, INK, INK_MUTED, figsize,
                              open_frame, record_caption)

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
# House palette (2026-09-11), replacing an Okabe-Ito set chosen in this file.
MODEL_C = C["ACCENT"]            # the swept cell
BASE_C = C["BASE"]               # the no-groin reference
OBS_C = INK                      # the observed endpoint
MARK_C = INK_MUTED               # the structure

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

# The house style is the whole of it now; the local rcParams block that used
# to sit here set its own type stack, ink and tick sizes.
apply_style()


def load(M, f):
    """Change since 1984 per domain, LANDWARD-positive, metres.

    Two directory spellings, because the no-groin cell has no deterioration
    floor to name: the swept cells are `M<M>_be<be>_f<f>` while the baseline on
    disk is `M0_be<be>`. Only the suffixed name was tried until 2026-09-11, so
    the baseline never loaded and the legend advertised a dotted line these
    gifs never drew.
    """
    for name in (f"M{M:g}_be{BE}_f{f:.2f}", f"M{M:g}_be{BE}"):
        p = SWEEP/name/"shoreline_matrix.npy"
        if p.is_file():
            m = np.load(p)[:, PAD]          # metres, landward-+, both already
            return m - m[0]
    return None


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
    fig = plt.figure(figsize=figsize("double", aspect=0.60))
    ax = fig.add_axes([0.135, 0.165, 0.845, 0.615])

    # Title, subtitle and legend stay ON the canvas: a GIF is watched
    # standalone, with no caption file beside it in a viewer.
    fig.text(0.135, 0.952, f"M = {M:g} m/yr, f = {f:g}",
             fontsize=11, fontweight="bold", color=INK, ha="left",
             va="center")
    fig.text(0.135, 0.905,
             f"{WHY.get((M, f), '')} · period 1, 1984 to 2004",
             fontsize=8, color=INK_MUTED, ha="left", va="center")

    fig.legend(handles=[
        Line2D([], [], color=OBS_C, marker="s", markersize=3.6,
               linestyle="--", linewidth=1.6,
               label="observed 1984 to 2004, the endpoint"),
        Line2D([], [], color=MODEL_C, marker="o", markersize=3.4,
               linewidth=1.8, label=f"modelled, M {M:g}, f {f:g}"),
        Line2D([], [], color=BASE_C, linestyle=":", linewidth=1.4,
               label="modelled, no groin, same year"),
    ], loc="upper left", bbox_to_anchor=(0.132, 0.872), ncol=3,
        fontsize=8, handlelength=2.2, columnspacing=2.2, frameon=False)

    # The footnote paragraph goes to CAPTIONS.md beside the GIFs.

    def frame(index, M=M, f=f, ch=ch, nyr=nyr, ax=ax):
        t = min(index, nyr - 1)
        ax.clear()
        ax.axvspan(3.5, 8.5, color="0.94", zorder=0)
        ax.annotate("fit window, D4–D8", xy=(6.0, 0.015),
                    xycoords=("data", "axes fraction"), ha="center",
                    va="bottom", fontsize=7, color=INK_MUTED)
        ax.axhline(0.0, color=INK_MUTED, linewidth=0.8, ls=(0, (4, 3)),
                   zorder=1)
        ax.axvline(5.5, color=MARK_C, lw=0.8, ls=(0, (5, 3)), zorder=2)
        ax.annotate("the groin field: D6 updrift, D5 downdrift",
                    xy=(5.5, 1.012), xycoords=("data", "axes fraction"),
                    ha="center", va="bottom", fontsize=7, color=MARK_C)

        ax.plot(DOM, c(obs), "s--", ms=3.6, lw=1.6, color=OBS_C, zorder=6)
        if nog is not None:
            ax.plot(DOM, c(nog[t]), ":", lw=1.4, color=BASE_C, zorder=4)
        ax.plot(DOM, c(ch[t]), "-o", ms=3.4, lw=1.8, color=MODEL_C, zorder=5)

        # -98 not -70: M = 160 dives past -75 at D12 and ran into the
        # fillet readout. All six cells share one limit so they stay
        # directly comparable.
        ax.set_xlim(1.6, 12.4); ax.set_ylim(-98, 70)
        ax.set_xticks(DOM)
        ax.set_xlabel("GIS domain (D2 toward Cape Point → D12 north)",
                      labelpad=6)
        # 26, not 40: the axes start at 0.135 of the figure width, so a
        # bigger pad pushes the label off the canvas. The orientation
        # cues moved further out instead.
        ax.set_ylabel("shoreline change since 1984 (m)\n"
                      "demeaned over D4−D8", labelpad=8)
        ax.grid(axis="y"); ax.set_axisbelow(True)
        open_frame(ax)

        # Orientation cues, so "up = erosion" needs no caption to decode.
        ax.annotate("▲ erosion, landward", xy=(0.012, 0.975),
                      xycoords="axes fraction", ha="left", va="top",
                      fontsize=7.5, color=INK_MUTED,
                      bbox=dict(facecolor="white", alpha=0.85,
                                edgecolor="none",
                                boxstyle="square,pad=0.12"))
        ax.annotate("▼ accretion, seaward", xy=(0.012, 0.025),
                      xycoords="axes fraction", ha="left",
                      va="bottom", fontsize=7.5, color=INK_MUTED,
                      bbox=dict(facecolor="white", alpha=0.85,
                                edgecolor="none",
                                boxstyle="square,pad=0.12"))

        ax.annotate(str(START + t), xy=(0.988, 0.965),
                    xycoords="axes fraction", ha="right", va="top",
                    fontsize=24, color="0.88", fontweight="bold", zorder=0)

        # The fillet, model against the fixed observed endpoint.
        up_i, down_i = DOM.index(6), DOM.index(5)
        model_fillet = ch[t][down_i] - ch[t][up_i]
        observed_fillet = obs[down_i] - obs[up_i]
        ax.annotate(f"fillet, D5 − D6:  modelled {model_fillet:+.0f} m,"
                    f"  observed {observed_fillet:+.0f} m",
                    xy=(0.988, 0.035), xycoords="axes fraction", ha="right",
                    va="bottom", fontsize=8, color=INK)

    anim = FuncAnimation(fig, frame, frames=nyr + HOLD_FRAMES, interval=420)
    p = OUT/f"shoreline_D2-D12_M{M:g}_f{f:.2f}.gif"
    anim.save(p, writer=PillowWriter(fps=2.4))
    plt.close(fig); made.append(p.name); print(f"  {p.name}")

    record_caption(
        p,
        "One period-1 sweep cell animated year by year, M = {M:g} m/yr with "
        "f = {f:g}, against the observed 1984 to 2004 change and the no-groin "
        "run of the same year, at be1 = {be}. {why} The frames keep their own "
        "title, legend and year clock because a GIF is watched standalone; "
        "everything else about how to read them is here. The y axis is "
        "LANDWARD-POSITIVE, so erosion moves up and the panel reads as a plan "
        "view with the ocean below and the island above. Every curve is "
        "DEMEANED over the D4–D8 fit window, because a uniform alongshore "
        "offset belongs to the source/sink calibration rather than to the "
        "groin, so only the SHAPE is being compared. The observed target is a "
        "fixed endpoint because that is all this window has; the full-life GIF "
        "is the one that animates the surveys. All cells share one y limit so "
        "they stay directly comparable."
        .format(M=M, f=f, be=BE, why=WHY.get((M, f), "").strip()))

print(f"\n{len(made)} gifs -> {OUT}")
print("  y axis is LANDWARD-positive: erosion up, matching the full-life gif")
