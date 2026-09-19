#!/usr/bin/env python3
"""How the module's logic produces (and fails to produce) the observed gap.

The two-shoreline figure shows WHAT happened. This one shows WHY the module can
follow part of it and not the rest, by putting the mechanics next to the
consequence.

THE MECHANICS (left panel)
    `GroinCallback` is four lines of arithmetic, applied once a year just
    before BRIE's alongshore solve:

        M_eff        = M, tapering to M*f between 1996 and 2003
        dx_updrift   = -M_eff      seaward advance at D6
        dx_downdrift = +M_eff      landward retreat at D5

    So each year the module pushes the two sides APART by 2*M_eff, and BRIE's
    diffusion immediately starts spreading that dipole back out. The modelled
    gap is the balance of those two.

WHY THAT MATTERS (right panel)
    M_eff is bounded at >= 0. There is no value of M or f that makes the module
    pull the two sides together. It can widen the gap, or -- at f = 0 -- stop
    widening it and let diffusion slowly close it. It can never actively close
    it.

    The observations do close it after 2004. Diffusion alone manages a small
    fraction of that; the fraction is computed here from the groin-off run and
    quoted in the caption. So the post-2004 narrowing is outside the
    parameterisation, and no choice of (M, f) reaches it. That is the single
    fact that determines how the module should be applied.

DATA
    Observed gap from the wet/dry surveys; modelled gap from the continuous
    1984-2024 sweep cells, which carry a full annual trajectory each.

STYLE, 2026-09-11
    Drawn under the project house style (`scripts/site_layer/hat_figure_style.py`): a
    190 mm column rather than a 15 in canvas, the sheltered side in ACCENT
    purple and the unprotected side in BASE grey to match the two-shoreline
    figure, and the schematic's two shouted sentences and the three footnote
    paragraphs moved into CAPTIONS.md beside the image.

    The "about a tenth" that footnote asserted is now measured off the
    groin-off run each time the figure is drawn, and it does NOT survive the
    measurement as it was stated. Diffusion alone closes the gap by about a
    sixth of the observed NET change across the window, which is where that
    number came from; its post-2004 RATE is about a sixty-sixth of the observed
    closure rate, which is what the sentence claimed. The caption carries both,
    because they disagree by a factor of ten and the rate is the one the
    argument rests on.

Usage:
    python HAT_groin_module_logic_figure.py

Writes groin_module_logic.png (and .pdf) beside this file.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import pathlib
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
GROIN_DIR = HERE.parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)

WETDRY_TABLE = (GROIN_DIR / "HAT-groin-buxton-output" / "shoreline_position_output"
                / "Change_from_wetdry_1967_D2_D12.csv")
SWEEP_ROOT = REPO / "output" / "calibration" / "groin" / "fullperiod_1984_2024"

UP_COLOR, DOWN_COLOR = C["ACCENT"], C["BASE"]
STORM_YEAR = 2003


def observed_gap():
    """{year: gap in metres} = D5 retreat minus D6 retreat, from 1967 datum."""
    frame = pd.read_csv(WETDRY_TABLE).set_index("Domain_ID")
    out = {}
    for column in frame.columns:
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", column)
        if not match:
            continue
        up, down = frame.loc[6, column], frame.loc[5, column]
        if pd.isna(up) or pd.isna(down):
            continue
        out[int(match.group(1))] = float(down - up)
    return dict(sorted(out.items()))


def modelled_gap(combo):
    """Modelled gap per year for one sweep cell, referenced to its own year 0."""
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as geometry
    path = SWEEP_ROOT / combo / "shoreline_matrix.npy"
    if not path.exists():
        return None, None
    matrix = np.load(path)
    up, down = geometry.gis_to_pad(6), geometry.gis_to_pad(5)
    gap = matrix[:, down] - matrix[:, up]
    return 1984 + np.arange(matrix.shape[0]), gap - gap[0]


def rate(years, series, start, end):
    """Endpoint rate of `series` over the samples inside [start, end]."""
    years = np.asarray(years)
    inside = np.where((years >= start) & (years <= end))[0]
    if inside.size < 2:
        return float("nan")
    first, last = inside[0], inside[-1]
    return (series[last] - series[first]) / (years[last] - years[first])


def draw_mechanics(axis):
    """Schematic: what the module adds to the two domains each year.

    Laid out for a HALF of a 190 mm column, which is about a quarter of the
    canvas this panel had before 2026-09-11. Every label was rewritten short
    and moved OUTWARD from the structure: at the printed width the old
    full-width text ran through the groin and through its own arrows, and a
    label touching an arrow reads as a label FOR that arrow."""
    axis.set_xlim(0, 10)
    axis.set_ylim(0, 10)
    axis.axis("off")

    # ocean / land reference
    axis.axhspan(0.4, 2.4, color=C["WATER"], alpha=0.6, zorder=0)
    axis.text(0.25, 0.75, "ocean", fontsize=7.5, color=INK_MUTED)
    axis.axhspan(8.2, 9.8, color="0.93", zorder=0)
    axis.text(0.25, 9.4, "land", fontsize=7.5, color=INK_MUTED)

    # the structure
    axis.plot([5, 5], [2.2, 8.0], color=INK, linewidth=3.4,
              solid_capstyle="butt", zorder=5)
    axis.text(5, 8.15, "groin", ha="center", fontsize=8, color=INK)

    # the two shorelines, each labelled at its OUTER end and on the side away
    # from the arrows
    axis.plot([0.8, 5], [6.2, 6.2], color=DOWN_COLOR, linewidth=2.4, zorder=4)
    axis.text(0.8, 6.45, "downdrift (GIS 5)", fontsize=7.5, color=DOWN_COLOR,
              va="bottom")
    axis.plot([5, 9.2], [4.6, 4.6], color=UP_COLOR, linewidth=2.4, zorder=4)
    axis.text(9.2, 4.3, "updrift (GIS 6)", fontsize=7.5, color=UP_COLOR,
              ha="right", va="top")

    # What the module applies each year.
    # Each label sits beyond its own arrowhead and on the SAME SIDE of the
    # groin as its arrow, on one line: beside the shaft it landed at the same
    # height as the shoreline label next to it, and on the far side it ran
    # straight through the structure.
    axis.annotate("", xy=(4.3, 7.6), xytext=(4.3, 6.35),
                  arrowprops=dict(arrowstyle="-|>", color=DOWN_COLOR, lw=1.8))
    axis.text(4.05, 7.55, "+M_eff, retreat", fontsize=7.5, color=DOWN_COLOR,
              ha="right", va="center")
    axis.annotate("", xy=(5.7, 3.1), xytext=(5.7, 4.45),
                  arrowprops=dict(arrowstyle="-|>", color=UP_COLOR, lw=1.8))
    axis.text(5.95, 3.05, "-M_eff, advance", fontsize=7.5, color=UP_COLOR,
              ha="left", va="center")

    # diffusion pushing back, across the structure
    axis.annotate("", xy=(4.1, 5.4), xytext=(5.9, 5.4),
                  arrowprops=dict(arrowstyle="<|-|>", color=INK_MUTED, lw=1.4,
                                  linestyle="--"))
    axis.text(5.0, 1.5, "alongshore diffusion spreads the dipole\n"
                        "and works to close the gap",
              ha="center", va="center", fontsize=7.5, color=INK_MUTED)


def main():
    apply_style()

    obs = observed_gap()
    window = {y: v for y, v in obs.items() if 1984 <= y <= 2024}
    base_year = min(window)
    years = np.array(sorted(window))
    values = np.array([window[y] - window[base_year] for y in years])

    fig, (left, right) = plt.subplots(
        1, 2, figsize=figsize("double", aspect=0.44),
        gridspec_kw=dict(width_ratios=[1, 1.25]), constrained_layout=True)

    draw_mechanics(left)
    _title(left, 0, "what the module applies each year")

    right.axhline(0.0, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 3)),
                  zorder=1)
    right.plot(years, values, marker="o", markersize=3.4, color=INK,
               linewidth=1.6, zorder=6, label="observed gap")

    runs = {}
    for combo, style, label in (
            ("M0", dict(color=DOWN_COLOR, linestyle=(0, (1, 2)), linewidth=1.6),
             "modelled, groin off: diffusion only"),
            # The chosen pair is (60, 0.6), but the continuous-window sweep's
            # f grid runs 0.1/0.3/0.5/0.7/0.9, so f=0.6 was never run here.
            # M60_f0.50 is its nearest neighbour and the label says so.
            ("M60_f0.50", dict(color=UP_COLOR, linewidth=1.8),
             "modelled, groin on at M 60, f 0.5, the nearest cell run to the"
             " chosen f 0.6")):
        span, gap = modelled_gap(combo)
        if span is None:
            continue
        runs[combo] = (span, gap)
        right.plot(span, gap, zorder=5, label=label, **style)

    right.axvline(STORM_YEAR, color=INK_MUTED, linestyle=(0, (1, 2)),
                  linewidth=0.8, zorder=2)
    right.text(STORM_YEAR + 0.6, 0.03, "storm damage 2003", rotation=90,
               fontsize=7, color=INK_MUTED, va="bottom", zorder=6,
               transform=right.get_xaxis_transform())

    right.set_xlabel("year")
    right.set_ylabel("change in the gap since 1984 (m)")
    right.grid(axis="y")
    right.set_axisbelow(True)
    open_frame(right)
    _title(right, 1, "what that produces, against the surveys")

    fig.legend(loc="outside lower center", ncol=3, frameon=False)

    last = int(years[-1])
    obs_release = rate(years, values, 2004, last)
    caption_parts = [
        "(a) The module's arithmetic, applied once a year just before BRIE's "
        "alongshore solve: the sheltered updrift domain is advanced seaward by "
        "the effective trapping rate and the downdrift domain retreated "
        "landward by the same amount, so the two shorelines are pushed apart by "
        "twice that rate every year. Alongshore diffusion then spreads the "
        "dipole and works to close the gap again, and the modelled gap is the "
        "balance of the two, saturating where they cancel. The rate is bounded "
        "at zero, so there is no (M, f) that makes the module close the gap: it "
        "can widen it, or at f = 0 stop widening it and let diffusion close it "
        "slowly. (b) The consequence over the net 1984-2024 window. The "
        "observed gap changes {obs:+.0f} m across it, closing at "
        "{obs_rate:+.1f} m/yr after 2004".format(
            obs=values[-1], obs_rate=obs_release),
    ]
    if "M0" in runs:
        span0, gap0 = runs["M0"]
        m0_release = rate(span0, gap0, 2004, 2024)
        rate_share = abs(m0_release / obs_release) if obs_release else float("nan")
        net_share = abs(gap0[-1] / values[-1]) if values[-1] else float("nan")
        # BOTH comparisons, because they disagree by a factor of ten and the
        # footnote this caption replaces quoted only the flattering one. It
        # said diffusion manages "about a tenth" of the observed closure; that
        # holds for the net change across the window ({net_share}), not for the
        # post-2004 rate ({rate_share}), which is what the sentence claimed.
        caption_parts.append(
            "; with the groin off, diffusion alone closes it by "
            "{m0_net:+.0f} m across the window, {net:.0%} of the observed "
            "change, but at only {m0:+.2f} m/yr after 2004, {rt:.1%} of the "
            "observed closure rate. That rate is the most the parameterisation "
            "can deliver at f = 0, short of the observed closure by a factor "
            "of about {factor:.0f}"
            .format(m0_net=gap0[-1], net=net_share, m0=m0_release,
                    rt=rate_share, factor=1.0 / rate_share))
    if "M60_f0.50" in runs:
        caption_parts.append(
            ". Turning the groin on moves the modelled line the other way, to "
            "{on:+.0f} m over the window, away from the surveys after 2004"
            .format(on=runs["M60_f0.50"][1][-1]))
    caption_parts.append(
        ". Read this as the module's BOUND, not as a test of the groin: the net "
        "1984-2024 window cancels period 1's build against period 2's collapse, "
        "so a module that can only widen the gap cannot fit it however it is "
        "scored. Within period 1 alone the groin does help, closing 23% of the "
        "shape misfit on GIS domains 4-8 (recorded from "
        "output/calibration/groin/figures/why_M60_f06.png when this figure was built, "
        "not computed here). Fit on period 1; period 2's closure is carried by "
        "the source/sink calibration. The modelled cell drawn here is M 60, "
        "f 0.5: the continuous-window sweep's f grid is 0.1/0.3/0.5/0.7/0.9, so "
        "the chosen f 0.6 was never run in it.")
    caption(fig, "".join(caption_parts))

    written = save(fig, HERE / "groin_module_logic.png", close=True)
    for path in written:
        print("wrote {}".format(path))
    print("  observed gap change 1984-{}: {:+.1f} m".format(last, values[-1]))
    for combo, (span, gap) in runs.items():
        print("  {:<12} modelled: {:+.1f} m   post-2004 {:+.2f} m/yr"
              .format(combo, gap[-1], rate(span, gap, 2004, 2024)))
    print("  observed post-2004: {:+.2f} m/yr".format(obs_release))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
