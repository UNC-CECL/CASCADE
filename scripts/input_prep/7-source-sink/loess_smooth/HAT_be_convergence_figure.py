#!/usr/bin/env python3
"""Did the source/sink calibration converge, and was its zone set fixed in advance?

Those are the two questions a reader has to be able to answer, because the
calibration is a FIXED-POINT SOLVE rather than a closed-form one. Where it
stopped is a scientific claim -- "this is the model's limit" -- and a stopping
point is only meaningful if the sequence was contracting and the target was not
moving while it ran.

WHY ITERATE AT ALL
    The ordinary calibration measures the residual of a base run and imposes it
    as the background-erosion field, which assumes that giving a domain X m/yr
    moves that domain's shoreline rate by X m/yr. It does not: BRIE diffuses an
    imposed rate alongshore and the domain keeps only a fraction of it. Measured
    here, one pass closes 42% of the misfit in period 1 and 57% in period 2 --
    so a one-shot residual mixes "the model cannot reproduce this" with "the
    correction was only half applied", and no reader can separate them.
    Iterating removes the second, leaving a residual that means one thing.

    Iterating also needs no estimate of the surviving fraction g, which matters
    because g is not a constant: a contiguous same-signed block of corrections
    passes at ~0.8-1.2 while a pattern alternating at the grid scale is damped
    to ~0.1. Dividing by g instead would amplify narrow features roughly tenfold
    into rates that are indefensible read as sediment fluxes.

WHY THE ZONES ARE FROZEN, AND WHY THE LEFT PANEL SHOWS THE RUN THAT SCORED
BETTER
    Zone membership is the scientific step: it says this stretch of coast has a
    real sediment-budget deficit and here is the process. Magnitude is
    arithmetic. Iterating both lets the arithmetic rewrite the science, because
    each pass re-derives zones from a NEW residual -- so as coherent features
    are satisfied, less coherent ones cross the threshold. And since adding BE
    at a domain pushes sediment into its neighbours, later passes partly correct
    the spillover of earlier ones, which never terminates.

    The unmasked run is drawn because it scored BETTER (dashed, right of the
    converged points). Hiding it would be the wrong kind of tidy: the gap is the
    fit available only by correcting outside justifiable zones, and the argument
    for this calibration is that the gap was declined deliberately, which the
    reader can only weigh by seeing its size.

WHAT THE RIGHT PANEL IS FOR
    To show the zone set was fixed BEFORE the iteration ran, not grown to fit.
    D5-D7 are marked separately: they are the groin's own footprint, reserved so
    the source/sink field cannot absorb the groin's shortfall and double-count
    against the M/f fit. D6 carries the largest residual in both periods and is
    deliberately never corrected.

Usage:
    python HAT_be_convergence_figure.py

Reads  output/convergence_history.json, and the live FROZEN_ZONE_DOMAINS /
       GROIN_RESERVED_DOMAINS / HATTERAS_BE_RATES_CALIBRATED, so the figure
       cannot drift from the calibration it documents.
Writes output/fig_be_convergence.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import importlib.util
import json
import pathlib
import sys

import numpy as np

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
OUTPUT_DIR = _HERE.parent / "output"
HISTORY = OUTPUT_DIR / "convergence_history.json"

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

PERIOD_LABEL = {"1984_2004": "1984–2004", "2004_2024": "2004–2024"}
PERIOD_KEY = {"1984_2004": 1984, "2004_2024": 2004}
COLOUR = {"1984_2004": "#1565C0", "2004_2024": "#B71C1C"}
RESERVED_COLOUR = "#FF8C00"


def load_calibration():
    """The live constants, imported rather than copied.

    A figure that hardcodes the zone set would keep rendering happily after
    someone edited the calibration, which is the failure mode this exists to
    guard against.
    """
    spec = importlib.util.spec_from_file_location(
        "_loess_analysis", _HERE.parent / "HAT_be_zone_LOESS_analysis.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from hatteras_site_config import HATTERAS_BE_RATES_CALIBRATED as rates
    return module.FROZEN_ZONE_DOMAINS, module.GROIN_RESERVED_DOMAINS, rates


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    history = json.loads(HISTORY.read_text(encoding="utf-8"))
    frozen, reserved, rates = load_calibration()
    periods = list(history["passes"])

    figure, (left, right) = plt.subplots(
        1, 2, figsize=(15.5, 6.8), gridspec_kw=dict(width_ratios=[1, 1.35]))

    # ---- LEFT: the convergence sequence ---------------------------------
    for period in periods:
        colour = COLOUR[period]
        passes = history["passes"][period]
        x = [p["pass"] for p in passes]
        y = [p["rmse"] for p in passes]

        left.plot(x, y, marker="o", markersize=8, color=colour, linewidth=2.4,
                  zorder=5, label=f"{PERIOD_LABEL[period]}  (frozen zones)")
        # Labels sit at SEGMENT MIDPOINTS, not on the markers. A gain belongs to
        # the step, not the endpoint, and at the markers the two periods' labels
        # collided with each other and with the lines.
        for step in range(1, len(x)):
            xm = (x[step - 1] + x[step]) / 2.0
            ym = (y[step - 1] + y[step]) / 2.0
            dy = 15 if period == periods[0] else -21
            left.annotate(f"{passes[step]['gain_pct']:.1f}%",
                          xy=(xm, ym), xytext=(0, dy),
                          textcoords="offset points", ha="center",
                          fontsize=9, color=colour, weight="bold")
        left.plot(x[-1], y[-1], marker="o", markersize=15, markerfacecolor="none",
                  markeredgecolor=colour, markeredgewidth=2.0, zorder=6)

        # The abandoned run, drawn because it scored better -- see the docstring.
        unmasked = history["_abandoned_unmasked"][period]
        left.plot([x[-1] + 0.55], [unmasked], marker="x", markersize=11,
                  color=colour, markeredgewidth=2.2, linestyle="none", zorder=6)
        left.plot([x[-1], x[-1] + 0.55], [y[-1], unmasked], color=colour,
                  linestyle="--", linewidth=1.2, alpha=0.55, zorder=4)

    left.set_xlabel("iteration pass")
    left.set_ylabel("shoreline-rate RMSE vs CoastSat target, D2–D89 (m/yr)")
    left.set_xticks(sorted({p["pass"] for pp in history["passes"].values()
                            for p in pp}))
    left.set_title("CONVERGENCE\neach pass adds the residual the last one left; "
                   "% is what that pass bought", fontsize=11.5)
    left.grid(alpha=0.25)
    # Scaled to the SEQUENCE. edgeBE and zeroBE are 2-4x these values and drawing
    # them as lines squashed the whole iteration into the bottom fifth of the
    # panel, which defeats the point of the figure; they are stated instead.
    left.set_ylim(0.44, 0.82)
    baselines = "   |   ".join(
        f"{PERIOD_LABEL[p]}: edgeBE {history['baselines']['edgeBE'][p]:.2f}, "
        f"zeroBE {history['baselines']['zeroBE'][p]:.2f}" for p in periods)
    left.text(0.015, 0.02, f"off-scale above — {baselines}", fontsize=7.6,
              color="#444444", ha="left", va="bottom", style="italic",
              transform=left.transAxes)
    handles = [Line2D([], [], color=COLOUR[p], marker="o", linewidth=2.4,
                      label=PERIOD_LABEL[p]) for p in periods]
    handles += [
        Line2D([], [], color="#555555", marker="x", linestyle="none",
               markeredgewidth=2.2, label="unmasked run — abandoned"),
        Line2D([], [], color="#555555", linestyle="--", linewidth=1.0,
               alpha=0.5, label="edgeBE / zeroBE baselines")]
    left.legend(handles=handles, loc="upper right", fontsize=8.5)

    # ---- RIGHT: the frozen zone set --------------------------------------
    for row, period in enumerate(periods):
        gis = PERIOD_KEY[period]
        members = set(frozen[gis])
        y0 = row * 1.0
        for domain in range(1, 91):
            if domain in reserved:
                face, alpha = RESERVED_COLOUR, 0.95
            elif domain in members:
                face, alpha = COLOUR[period], 0.75
            else:
                face, alpha = "#DDDDDD", 0.7
            right.add_patch(plt.Rectangle((domain - 0.5, y0), 1.0, 0.62,
                                          facecolor=face, alpha=alpha,
                                          edgecolor="none"))
        right.text(-1.5, y0 + 0.31, PERIOD_LABEL[period], ha="right",
                   va="center", fontsize=10, color=COLOUR[period], weight="bold")
        right.text(91.5, y0 + 0.31, f"{len(members)} domains", ha="left",
                   va="center", fontsize=8.5, color=COLOUR[period])

    # the final BE field, on a shared axis below the zone bars
    scale = 0.075
    base = -0.62
    right.axhline(base, color="#999999", linewidth=0.9, zorder=2)
    for period in periods:
        gis = PERIOD_KEY[period]
        values = [rates[gis].get(d, 0.0) for d in range(2, 90)]
        clipped = np.clip(values, -4.0, 4.0)   # interior only; GIS 1/90 dwarf it
        right.plot(range(2, 90), base + np.array(clipped) * scale,
                   color=COLOUR[period], linewidth=1.6, alpha=0.9, zorder=3,
                   label=f"final BE, {PERIOD_LABEL[period]}")
    right.text(91.5, base, "final BE\n(m/yr)", ha="left", va="center",
               fontsize=8, color="#444444")

    right.set_xlim(-8, 100)
    right.set_ylim(base - 0.40, 2.05)
    right.set_yticks([])
    right.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    right.set_xlabel("GIS domain (south → north)")
    right.set_title("THE ZONE SET, FIXED BEFORE ITERATION\n"
                    "coloured = correctable; grey = withheld however large its "
                    "residual; orange = reserved for the groin", fontsize=11.5)
    right.legend(loc="upper right", fontsize=8)

    figure.tight_layout(rect=(0, 0.115, 1, 1))
    figure.text(
        0.01, 0.012,
        "WHY ITERATE. Imposing X m/yr of background erosion does not move a domain's rate by X -- BRIE diffuses most of it alongshore -- so the "
        "one-shot solve closes only 42% (P1) and 57% (P2) of the misfit and its residual conflates 'the model cannot do this' with 'the correction "
        "was half applied'. Each pass re-measures and adds what is left, which needs no estimate of the surviving fraction; that fraction is not "
        "constant anyway, running ~0.8-1.2 for a contiguous block of corrections and ~0.1 for one alternating at the grid scale.\n"
        "WHY THE ZONES ARE FROZEN. Zone membership is the science, magnitude is arithmetic, and iterating both lets the arithmetic rewrite the "
        "science: re-deriving zones each pass let 19 domains (P1) and 12 (P2) outside the original identification pick up corrections, and part of "
        "what later passes 'find' is the alongshore spillover of earlier passes. The abandoned unmasked run is plotted BECAUSE it scored better "
        "(0.4931 / 0.4843) -- that gap is the fit available only outside justifiable zones, and it was declined deliberately.\n"
        "WHAT IS LEFT. D6 carries the largest residual in both periods (2.00 and 2.59 m/yr) and is never corrected: it is the groin's own shortfall "
        "-- too little fillet built in period 1, no release in period 2 -- and absorbing it here would double-count against the M/f fit.",
        fontsize=7.3, color="#333333", wrap=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / "fig_be_convergence.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)

    print(f"wrote {path}")
    for period in periods:
        passes = history["passes"][period]
        edge = history["baselines"]["edgeBE"][period]
        print(f"  {PERIOD_LABEL[period]}: {passes[0]['rmse']:.4f} -> "
              f"{passes[-1]['rmse']:.4f} in {len(passes)-1} pass(es); "
              f"closed {(edge - passes[-1]['rmse'])/edge*100:.1f}% of edgeBE; "
              f"{len(frozen[PERIOD_KEY[period]])} domains correctable")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
