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

    The observations do close it, at -4.0 m/yr after 2004. Diffusion alone
    manages about a tenth of that. So the post-2004 narrowing is outside the
    parameterisation, and no choice of (M, f) reaches it. That is the single
    fact that determines how the module should be applied.

DATA
    Observed gap from the wet/dry surveys; modelled gap from the continuous
    1984-2024 sweep cells, which carry a full annual trajectory each.

Usage:
    python HAT_groin_module_logic_figure.py

Writes groin_module_logic.png beside this file.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import pathlib
import re

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
PROJECT_BASE_DIR = next(p for p in HERE.parents if (p / "pyproject.toml").exists())

WETDRY_TABLE = (HERE / "HAT-groin-test-output" / "shoreline_position_output"
                / "Change_from_wetdry_1967_D2_D12.csv")
SWEEP_ROOT = (PROJECT_BASE_DIR / "output" / "groin_sweep"
              / "fullperiod_1984_2024")

UP_COLOR, DOWN_COLOR, MODEL_COLOR = "#B71C1C", "#1565C0", "#FF8C00"


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
    import sys
    if str(PROJECT_BASE_DIR / "scripts") not in sys.path:
        sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))
    from hatteras_site_config import HATTERAS_DOMAINS as geometry
    path = SWEEP_ROOT / combo / "shoreline_matrix.npy"
    if not path.exists():
        return None, None
    matrix = np.load(path)
    up, down = geometry.gis_to_pad(6), geometry.gis_to_pad(5)
    gap = matrix[:, down] - matrix[:, up]
    return 1984 + np.arange(matrix.shape[0]), gap - gap[0]


def draw_mechanics(axis):
    """Schematic: what the module adds to the two domains each year."""
    axis.set_xlim(0, 10)
    axis.set_ylim(0, 10)
    axis.axis("off")

    # ocean / land reference
    axis.axhspan(0, 3.2, color="#BBDEFB", alpha=0.5)
    axis.text(0.3, 0.4, "OCEAN", fontsize=9, color="#0D47A1", weight="bold")
    axis.axhspan(7.4, 10, color="#FFF3E0", alpha=0.7)
    axis.text(0.3, 9.3, "LAND", fontsize=9, color="#E65100", weight="bold")

    # the groin
    axis.plot([5, 5], [3.0, 7.2], color="#333333", linewidth=7,
              solid_capstyle="butt", zorder=5)
    axis.text(5, 7.5, "GROIN", ha="center", fontsize=10, weight="bold")

    # the two shorelines
    axis.plot([1.2, 5], [5.6, 5.6], color=DOWN_COLOR, linewidth=3.4, zorder=4)
    axis.text(1.2, 5.85, "D5  DOWNDRIFT", fontsize=9.5, color=DOWN_COLOR,
              weight="bold", va="bottom")
    axis.plot([5, 8.8], [4.3, 4.3], color=UP_COLOR, linewidth=3.4, zorder=4)
    axis.text(8.8, 4.5, "D6  UPDRIFT", fontsize=9.5, color=UP_COLOR,
              weight="bold", ha="right", va="bottom")

    # What the module applies each year. Arrows and their labels sit clear of
    # the domain labels above -- an overlap here reads as one being a caption
    # for the other, which is exactly the wrong reading.
    axis.annotate("", xy=(3.9, 7.0), xytext=(3.9, 5.7),
                  arrowprops=dict(arrowstyle="-|>", color=DOWN_COLOR, lw=3))
    axis.text(4.1, 6.5, "+M_eff\nRETREAT", fontsize=9, color=DOWN_COLOR,
              weight="bold", va="center")
    axis.annotate("", xy=(6.1, 3.0), xytext=(6.1, 4.2),
                  arrowprops=dict(arrowstyle="-|>", color=UP_COLOR, lw=3))
    axis.text(6.3, 3.4, "-M_eff\nADVANCE", fontsize=9, color=UP_COLOR,
              weight="bold", va="center")

    # diffusion pushing back
    axis.annotate("", xy=(4.3, 5.0), xytext=(5.8, 5.0),
                  arrowprops=dict(arrowstyle="<|-|>", color="#2E7D32", lw=2.4,
                                  linestyle="--"))
    axis.text(5.05, 2.05, "BRIE diffusion spreads the dipole\n"
                          "and works to CLOSE the gap",
              ha="center", fontsize=9, color="#2E7D32", weight="bold")

    axis.text(5.0, 9.0,
              "each year the module pushes the two sides APART by 2 x M_eff",
              ha="center", fontsize=10, weight="bold")
    axis.text(5.0, 0.95,
              "M_eff >= 0 ALWAYS  ->  the module can only WIDEN the gap,\n"
              "or stop widening it. It can NEVER pull the sides together.",
              ha="center", fontsize=9.5, color="#B71C1C", weight="bold",
              bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFEBEE",
                        edgecolor="#B71C1C"))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    obs = observed_gap()
    window = {y: v for y, v in obs.items() if 1984 <= y <= 2024}
    base_year = min(window)
    years = np.array(sorted(window))
    values = np.array([window[y] - window[base_year] for y in years])

    figure, (left, right) = plt.subplots(1, 2, figsize=(15, 7),
                                         gridspec_kw=dict(width_ratios=[1, 1.25]))
    draw_mechanics(left)
    left.set_title("THE LOGIC\nwhat GroinCallback does every year", fontsize=12)

    right.axhline(0.0, color="#999999", linewidth=1.0, zorder=1)
    right.plot(years, values, marker="o", markersize=6, color="#1A1A1A",
               linewidth=2.4, zorder=6, label="OBSERVED gap")

    for combo, style, label in (
            ("M0", dict(color="#777777", linestyle=":", linewidth=2.2),
             "model, groin OFF (diffusion only)"),
            # The chosen pair is (60, 0.6), but the continuous-window sweep's
            # f grid runs 0.1/0.3/0.5/0.7/0.9, so f=0.6 was never run here.
            # M60_f0.50 is its nearest neighbour and the label says so.
            ("M60_f0.50", dict(color=MODEL_COLOR, linewidth=2.6),
             "model, groin ON (M=60, f=0.5 -- nearest cell to the chosen f=0.6)")):
        span, gap = modelled_gap(combo)
        if span is None:
            continue
        right.plot(span, gap, zorder=5, label=label, **style)

    right.axvline(2003, color="#B71C1C", linestyle=":", linewidth=1.6, zorder=2)
    right.text(2003.6, 0.03, "storm damage 2003", rotation=90, fontsize=8.5,
               color="#B71C1C", va="bottom",
               transform=right.get_xaxis_transform())

    right.annotate("the module widens the gap;\nobserved CLOSES it",
                   xy=(2018, -18), xytext=(2006, 4),
                   fontsize=10, color="#B71C1C", weight="bold",
                   arrowprops=dict(arrowstyle="->", color="#B71C1C", lw=1.8))

    right.set_xlabel("year")
    right.set_ylabel("change in the gap since 1984 (m)")
    right.set_title("THE CONSEQUENCE\nover the NET 1984-2024 window the module "
                    "cannot follow the closure", fontsize=12)
    right.grid(alpha=0.25)
    right.legend(loc="lower left", fontsize=9)

    figure.tight_layout(rect=(0, 0.10, 1, 1))
    figure.text(
        0.01, 0.015,
        "LEFT: the module adds -M_eff at D6 and +M_eff at D5 once a year, so it pushes the two shorelines apart by 2 x M_eff; "
        "BRIE's diffusion then spreads that dipole and works to close the gap again. The modelled gap is the balance of the two, "
        "and it saturates where they cancel.\n"
        "RIGHT: because M_eff is bounded at >= 0 there is no (M, f) that makes the module CLOSE the gap. Turning the groin on "
        "moves the orange line UP, away from the observations after 2004. The best the module can do is f = 0 -- stop widening "
        "and let diffusion close the gap slowly -- which is about a tenth of the observed -4.0 m/yr. That is why f -> 0 is the "
        "correct answer for period 2, and why the remaining closure is carried by the source/sink calibration.\n"
        "DO NOT READ THIS AS 'the groin never helps'. This panel is the NET 1984-2024 window, which cancels period 1's build "
        "(+52 m) against period 2's collapse (-76 m) -- so a module that can only widen the gap cannot win here, and this window "
        "cannot fit a groin however it is scored. Within PERIOD 1 alone the groin DOES help: on D4-D8 it closes 23% of the shape "
        "misfit (see why_M60_f06.png). Fit on period 1; use this panel to understand the module's BOUND, not to test the groin.",
        fontsize=8, color="#333333", wrap=True)

    path = HERE / "groin_module_logic.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    print(f"wrote {path}")
    print(f"  observed gap change 1984-2024: {values[-1]:+.1f} m")
    for combo in ("M0", "M60_f0.50"):
        span, gap = modelled_gap(combo)
        if span is not None:
            print(f"  {combo:<12} modelled: {gap[-1]:+.1f} m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
