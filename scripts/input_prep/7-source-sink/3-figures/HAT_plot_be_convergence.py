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
    python 3-figures/HAT_plot_be_convergence.py

Reads  2-calibrate/1984_2004__2004_2024/convergence_history.json, and the live FROZEN_ZONE_DOMAINS /
       GROIN_RESERVED_DOMAINS / HATTERAS_BE_RATES_CALIBRATED, so the figure
       cannot drift from the calibration it documents.
Writes data/hatteras_init/7-source-sink/3-figures/1984_2004__2004_2024/2-method/fig_be_convergence.png (and the
       PDF beside it); the caption is written to CAPTIONS.md in that folder.

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
sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))
from site_layer import hat_source_sink as _be  # noqa: E402

# Reads and writes the data tree, where the fit now puts its products: the
# default pair's folder (resolved by hat_source_sink.py since 2026-09-18).
OUTPUT_DIR = _be.calibrate_dir()
HISTORY = OUTPUT_DIR / "convergence_history.json"
# The figure belongs with the rest of the section 7 figures, in the data tree;
# the iteration's own record stays beside the calibration that wrote it.
FIG_DIR = _be.figures_dir()

from site_layer.hat_figure_style import (                                   # noqa: E402
    apply_style, figsize, save, caption, town_bands, open_frame,
    DOMAIN_AXIS_LABEL, C, C_1984, C_1997, INK, INK_MUTED, _title)

PERIOD_LABEL = {"1984_2004": "1984–2004", "2004_2024": "2004–2024"}
PERIOD_KEY = {"1984_2004": 1984, "2004_2024": 2004}
# The earlier period is the red of the house vintage pair, the later the blue.
COLOUR = {"1984_2004": C_1984, "2004_2024": C_1997}


def load_calibration():
    """The live constants, imported rather than copied.

    A figure that hardcodes the zone set would keep rendering happily after
    someone edited the calibration, which is the failure mode this exists to
    guard against.
    """
    spec = importlib.util.spec_from_file_location(
        "_loess_analysis",
        _HERE.parent.parent / "2-calibrate" / "HAT_be_zone_residual_fit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from site_layer.hatteras_site_config import HATTERAS_BE_RATES_CALIBRATED as rates
    return module.FROZEN_ZONE_DOMAINS, module.GROIN_RESERVED_DOMAINS, rates


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    apply_style()

    history = json.loads(HISTORY.read_text(encoding="utf-8"))
    frozen, reserved, rates = load_calibration()
    periods = list(history["passes"])

    figure, (left, right) = plt.subplots(
        1, 2, figsize=figsize("double", aspect=0.46),
        constrained_layout=True, gridspec_kw=dict(width_ratios=[1, 1.45]))

    # ---- LEFT: the convergence sequence ---------------------------------
    for period in periods:
        colour = COLOUR[period]
        passes = history["passes"][period]
        x = [p["pass"] for p in passes]
        y = [p["rmse"] for p in passes]

        left.plot(x, y, marker="o", markersize=4.5, color=colour, linewidth=1.6,
                  zorder=5, label=PERIOD_LABEL[period])
        # Labels sit at SEGMENT MIDPOINTS, not on the markers. A gain belongs to
        # the step, not the endpoint, and at the markers the two periods' labels
        # collided with each other and with the lines.
        for step in range(1, len(x)):
            xm = (x[step - 1] + x[step]) / 2.0
            ym = (y[step - 1] + y[step]) / 2.0
            # The earlier period runs BELOW its line and the later one above,
            # so the two sets of gains cannot meet in the middle.
            dy = -15 if period == periods[0] else 11
            left.annotate(f"{passes[step]['gain_pct']:.1f}%",
                          xy=(xm, ym), xytext=(0, dy),
                          textcoords="offset points", ha="center",
                          fontsize=7, color=colour)
        left.plot(x[-1], y[-1], marker="o", markersize=10,
                  markerfacecolor="none", markeredgecolor=colour,
                  markeredgewidth=1.2, zorder=6)

        # The abandoned run, drawn because it scored better -- see the docstring.
        unmasked = history["_abandoned_unmasked"][period]
        left.plot([x[-1] + 0.55], [unmasked], marker="x", markersize=7,
                  color=colour, markeredgewidth=1.6, linestyle="none", zorder=6)
        left.plot([x[-1], x[-1] + 0.55], [y[-1], unmasked], color=colour,
                  linestyle=(0, (3, 2)), linewidth=0.9, zorder=4)

    left.set_xlabel("iteration pass")
    left.set_ylabel("shoreline-rate RMSE against the\nCoastSat target, D2\u2013D89 (m/yr)")
    left.set_xticks(sorted({p["pass"] for pp in history["passes"].values()
                            for p in pp}))
    # Scaled to the SEQUENCE. edgeBE and zeroBE are 2-4x these values and drawing
    # them as lines squashed the whole iteration into the bottom fifth of the
    # panel, which defeats the point of the figure; they are in the caption.
    left.set_ylim(0.44, 0.82)
    left.grid(axis="y")
    open_frame(left)
    _title(left, 0, "convergence")
    handles = [Line2D([], [], color=COLOUR[p], marker="o", markersize=4.5,
                      linewidth=1.6, label=PERIOD_LABEL[p]) for p in periods]
    handles += [
        Line2D([], [], color=INK_MUTED, marker="x", linestyle="none",
               markeredgewidth=1.6, label="zone set not imposed")]
    left.legend(handles=handles, loc="upper right", frameon=False, fontsize=7)

    # ---- RIGHT: the frozen zone set --------------------------------------
    for row, period in enumerate(periods):
        gis = PERIOD_KEY[period]
        members = set(frozen[gis])
        y0 = row * 1.0
        for domain in range(1, 91):
            kw = dict(facecolor=C["BASE_FILL"], edgecolor="none")
            if domain in reserved:
                kw = dict(facecolor="none", edgecolor=C["BASE"], hatch="///",
                          linewidth=0.0)
            elif domain in members:
                kw = dict(facecolor=COLOUR[period], edgecolor="none")
            right.add_patch(plt.Rectangle((domain - 0.5, y0), 1.0, 0.62,
                                          zorder=3, **kw))
        right.text(-1.5, y0 + 0.31, PERIOD_LABEL[period], ha="right",
                   va="center", fontsize=8, color=COLOUR[period])
        right.text(91.5, y0 + 0.31, f"{len(members)} domains", ha="left",
                   va="center", fontsize=7, color=INK_MUTED)

    # the final background-erosion field, on a shared axis below the zone bars
    scale = 0.075
    base = -0.62
    right.axhline(base, color=INK_MUTED, linewidth=0.6, zorder=2)
    for period in periods:
        gis = PERIOD_KEY[period]
        values = [rates[gis].get(d, 0.0) for d in range(2, 90)]
        clipped = np.clip(values, -4.0, 4.0)   # interior only; GIS 1/90 dwarf it
        right.plot(range(2, 90), base + np.array(clipped) * scale,
                   color=COLOUR[period], linewidth=1.0, zorder=3)
    right.text(91.5, base, "calibrated\nfield (m/yr)", ha="left", va="center",
               fontsize=7, color=INK_MUTED)

    right.set_xlim(-10, 104)
    right.set_ylim(base - 0.40, 2.85)
    right.set_yticks([])
    right.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    right.set_xlabel(DOMAIN_AXIS_LABEL)
    for side in ("top", "right", "left"):
        right.spines[side].set_visible(False)
    town_bands(right, where="top", strip=0.07, fontsize=7)
    _title(right, 1, "the zone set, fixed before the first pass")
    right.legend(handles=[
        Patch(facecolor=COLOUR[periods[0]],
              label=f"correctable, {PERIOD_LABEL[periods[0]]}"),
        Patch(facecolor=C["BASE_FILL"], label="withheld, left at zero"),
        Patch(facecolor=COLOUR[periods[1]],
              label=f"correctable, {PERIOD_LABEL[periods[1]]}"),
        Patch(facecolor="none", edgecolor=C["BASE"], hatch="///",
              label="reserved for the groin")],
        loc="upper center", bbox_to_anchor=(0.5, 0.93), ncol=2, frameon=False,
        fontsize=7)

    baselines = "; ".join(
        f"{PERIOD_LABEL[p]} {history['baselines']['edgeBE'][p]:.2f} and "
        f"{history['baselines']['zeroBE'][p]:.2f}" for p in periods)
    unmasked = ", ".join(f"{history['_abandoned_unmasked'][p]:.4f}"
                         for p in periods)
    caption(figure, (
        "The source/sink calibration is a fixed-point solve, so where it stopped "
        "is a claim about the model's limit, and that claim only means something "
        "if the sequence was contracting and the target was not moving while it "
        "ran. (a) the shoreline-rate RMSE against the CoastSat target after each "
        "pass, for each period; the percentage on a segment is what that pass "
        "bought, and the ringed marker is the pass the calibration stopped at. "
        "Imposing X m/yr of background erosion does not move a domain's rate by "
        "X -- BRIE diffuses most of it alongshore -- so a single pass closes only "
        "42 per cent of the misfit in the first period and 57 per cent in the "
        "second, and its residual conflates 'the model cannot do this' with 'the "
        "correction was half applied'. Each further pass re-measures and adds "
        "what is left, which needs no estimate of the surviving fraction; that "
        "fraction is not constant anyway, running about 0.8-1.2 for a contiguous "
        "block of corrections and about 0.1 for one alternating at the grid "
        "scale. The crosses are the same iteration run without the zone set "
        f"imposed; it scored better ({unmasked}) and was declined anyway, and "
        "the gap is the size of the fit available only by correcting outside "
        "justifiable zones. The axis is scaled to the sequence: the edgeBE and "
        f"zeroBE baselines are off the top at {baselines} m/yr. (b) the zone "
        "set, identified once from the first residual and then held for every "
        "pass, because zone membership is the scientific step and magnitude is "
        "arithmetic -- re-deriving the zones each pass would let less coherent "
        "features cross the threshold as real ones were satisfied, and let later "
        "passes correct the alongshore spillover of earlier ones. Domains 5-7 "
        "are the Buxton groin's own footprint, reserved so the source/sink field "
        "cannot absorb the groin's shortfall and double-count against the "
        "trapping fit; domain 6 carries the largest residual in both periods "
        "(2.00 and 2.59 m/yr) and is deliberately never corrected. The line "
        "below each pair of rows is the calibrated field over the interior "
        "domains, clipped to plus or minus 4 m/yr. Domain 1 is at Cape Point and "
        "domain 90 at Pea Island."))

    path = FIG_DIR / "2-method" / "fig_be_convergence.png"
    save(figure, path)
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
