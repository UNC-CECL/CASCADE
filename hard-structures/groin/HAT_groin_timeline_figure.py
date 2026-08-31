#!/usr/bin/env python3
"""The Buxton groin's observed history, and how it maps onto the hindcast.

One figure answering three questions that kept getting tangled:

    WHEN did the groin affect the shoreline?
        Continuously from 1969, but in two opposite regimes. It TRAPPED through
        2004 and has been RELEASING since. The turning point coincides with the
        2003 storm damage, not with anything in the model.

    IS THERE A DIFFERENT SIGNAL IN EACH HINDCAST PERIOD?
        Yes, and they have opposite sign. Period 1 (1984-2004) is still
        accumulating at +2.7 m/yr. Period 2 (2004-2024) is draining at
        -4.0 m/yr.

    HOW DOES THE MODULE REPRESENT THAT?
        `GroinCallback` carries an absolute calendar schedule -- install 1969,
        deterioration onset 1996 (last repair), linear ramp to 2003 (storm),
        then hold at M*f. The lower panel draws that schedule against the
        observations, which is the clearest way to see that the module's
        built-in timeline already matches the measured history, and where it
        cannot follow.

WHAT THE MODULE CANNOT DO, DRAWN RATHER THAN FOOTNOTED
    Trapping is bounded at >= 0, so the groin can stop adding sand but cannot
    actively drain the fillet. Period 2's observed -76 m release is therefore
    outside what the parameterisation can produce at any (M, f), and the lower
    panel shades that region so the limitation is visible next to the data.

DATA
    Fillet is x_s[D5] - x_s[D6] against a fixed 1967 datum, from
    `Change_from_wetdry_1967_D2_D12.csv` -- 24 dated wet/dry surveys produced by
    the GIS analysis in HAT-groin-gis-analysis. Landward-positive, so a rising
    curve means the updrift side is holding while the downdrift side retreats,
    which is what a groin builds.

Usage:
    python HAT_groin_timeline_figure.py

Writes groin_timeline_and_hindcast.png beside this file.

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

# The structure's documented history. These are the values GroinCallback is
# configured with, not fitted quantities.
INSTALL_YEAR = 1969
LAST_REPAIR_YEAR = 1996
STORM_YEAR = 2003

# The chosen parameters. M is an EFFECTIVE, grid-specific rate -- see the
# accompanying GROIN_PLAN.md for why it is not a sediment flux.
CHOSEN_M, CHOSEN_F = 60.0, 0.6

PERIODS = ((1984, 2004, "#1565C0"), (2004, 2024, "#2E7D32"))
UPDRIFT_GIS, DOWNDRIFT_GIS = 6, 5


def observed_fillet():
    """{year: fillet in metres} against the fixed 1967 datum."""
    frame = pd.read_csv(WETDRY_TABLE).set_index("Domain_ID")
    out = {}
    for column in frame.columns:
        # SECOND year is the survey; the first is the 1967 datum.
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", column)
        if not match:
            continue
        up, down = frame.loc[UPDRIFT_GIS, column], frame.loc[DOWNDRIFT_GIS, column]
        if pd.isna(up) or pd.isna(down):
            continue
        out[int(match.group(1))] = float(down - up)
    out.setdefault(1967, 0.0)
    return dict(sorted(out.items()))


def effective_trapping(year, M=CHOSEN_M, f=CHOSEN_F):
    """M_eff for a calendar year -- mirrors GroinCallback._effective_trapping_rate."""
    if year < INSTALL_YEAR:
        return 0.0
    if year < LAST_REPAIR_YEAR:
        return M
    ramp_years = STORM_YEAR - LAST_REPAIR_YEAR
    floor = M * f
    if year >= STORM_YEAR:
        return floor
    taper = (year - LAST_REPAIR_YEAR) / ramp_years
    return M - taper * (M - floor)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    obs = observed_fillet()
    years = np.array(sorted(obs))
    values = np.array([obs[y] for y in years])

    figure, (top, bottom) = plt.subplots(
        2, 1, figsize=(14, 9.5), sharex=True,
        gridspec_kw=dict(height_ratios=[3, 2]))

    # ---- hindcast period bands, on both panels --------------------------
    for axis in (top, bottom):
        for start, end, colour in PERIODS:
            axis.axvspan(start, end, color=colour, alpha=0.07, zorder=0)
        for year, label, colour in ((INSTALL_YEAR, "installed 1969", "#333333"),
                                    (LAST_REPAIR_YEAR, "last repair 1996", "#E65100"),
                                    (STORM_YEAR, "storm damage 2003", "#B71C1C")):
            axis.axvline(year, color=colour, linestyle=":", linewidth=1.6, zorder=2)

    # ---- TOP: the observed fillet ---------------------------------------
    top.plot(years, values, marker="o", markersize=6, color="#1A1A1A",
             linewidth=2.0, zorder=5, label="observed fillet (wet/dry surveys)")

    phases = ((1967, 1978, "#C8E6C9", "BUILD\n+10.7 m/yr"),
              (1978, 2004, "#FFF9C4", "slow growth\n+1.3 m/yr"),
              (2004, 2023, "#FFCDD2", "RELEASE\n-4.0 m/yr"))
    for start, end, colour, label in phases:
        top.axvspan(start, end, color=colour, alpha=0.45, zorder=1)
        top.text((start + end) / 2, 0.06, label, ha="center", va="bottom",
                 fontsize=9, transform=top.get_xaxis_transform(), zorder=6)

    for year, colour, label in ((INSTALL_YEAR, "#333333", "installed"),
                                (LAST_REPAIR_YEAR, "#E65100", "last repair"),
                                (STORM_YEAR, "#B71C1C", "storm damage")):
        top.text(year + 0.5, 0.965, label, rotation=90, fontsize=8, color=colour,
                 va="top", transform=top.get_xaxis_transform(), zorder=6)

    top.set_ylabel("FILLET (m)\n"
                   "how much LESS the updrift side (D6) has retreated\n"
                   "than the downdrift side (D5), since 1967")
    top.set_title(
        "Buxton groin -- observed fillet history and how it maps onto the "
        "1984-2024 hindcast\n"
        "FILLET = differential protection, NOT accumulation: both sides eroded, "
        "the updrift side just eroded far less", fontsize=12)
    top.grid(alpha=0.25)
    top.legend(loc="upper left", fontsize=9)

    # period annotations with the numbers that matter
    for (start, end, colour), text in zip(
            PERIODS,
            ("PERIOD 1  1984-2004\n+52 m  (+2.74 m/yr)\ngroin still TRAPPING",
             "PERIOD 2  2004-2024\n-76 m  (-4.02 m/yr)\ngroin RELEASING")):
        top.text((start + end) / 2, 0.97, text, ha="center", va="top",
                 fontsize=9.5, color=colour, weight="bold",
                 transform=top.get_xaxis_transform(), zorder=7,
                 bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                           edgecolor=colour, alpha=0.85))

    # ---- BOTTOM: what the module actually applies ------------------------
    span = np.arange(1967, 2025)
    m_eff = np.array([effective_trapping(y) for y in span])
    bottom.plot(span, m_eff, color="#FF8C00", linewidth=2.6, zorder=5,
                label=f"M_eff applied by GroinCallback (M={CHOSEN_M:g}, f={CHOSEN_F:g})")
    bottom.fill_between(span, 0, m_eff, color="#FF8C00", alpha=0.18, zorder=3)

    # The floor the module cannot go below.
    bottom.axhline(0.0, color="#B71C1C", linewidth=1.4, zorder=4)
    bottom.fill_between([2004, 2024], -12, 0, color="#B71C1C", alpha=0.13, zorder=2)
    bottom.text(2014, -6, "the module CANNOT go here:\ntrapping is bounded at >= 0, so it can stop\n"
                          "adding sand but cannot DRAIN the fillet",
                ha="center", va="center", fontsize=8.5, color="#B71C1C", zorder=6)

    bottom.set_ylim(-12, CHOSEN_M * 1.15)
    bottom.set_ylabel("effective trapping rate M_eff (m/yr)")
    bottom.set_xlabel("year")
    bottom.set_title(
        "The module's calendar schedule already encodes this history -- "
        "no period-specific configuration is needed", fontsize=10.5)
    bottom.grid(alpha=0.25)
    bottom.legend(loc="upper right", fontsize=9)

    bottom.annotate("full trapping", xy=(1982, CHOSEN_M), xytext=(1972, CHOSEN_M * 0.72),
                    fontsize=8.5, color="#E65100",
                    arrowprops=dict(arrowstyle="->", color="#E65100", lw=1.2))
    bottom.annotate(f"ramp to floor\nM x f = {CHOSEN_M * CHOSEN_F:g} m/yr",
                    xy=(2003, CHOSEN_M * CHOSEN_F),
                    xytext=(2008, CHOSEN_M * 0.62), fontsize=8.5, color="#E65100",
                    arrowprops=dict(arrowstyle="->", color="#E65100", lw=1.2))

    figure.tight_layout(rect=(0, 0.055, 1, 1))
    figure.text(
        0.01, 0.008,
        "Fillet is the D5-D6 shoreline offset against a fixed 1967 datum, from 24 dated wet/dry surveys (GIS analysis). It is a "
        "RELATIVE measure: by 2004 the downdrift side had retreated 175 m since 1967 while the sheltered updrift side had "
        "retreated only 25 m, and that 150 m gap is the fillet -- not 150 m of new beach. The post-2004 collapse is ~85% D6 "
        "eroding once the structure failed, not impounded sand draining downdrift. "
        "The observed turning point is 2004 -- immediately after the 2003 storm damage -- and the module's built-in ramp "
        "(onset at the 1996 repair, complete by the 2003 storm) reproduces that timing without period-specific settings. "
        "PERIOD 1 is the only window in the hindcast where the groin does something the module can reproduce; at M=60, f=0.6 "
        "it supplies about 33% of the observed +52 m. PERIOD 2's release is outside the parameterisation entirely and is "
        "carried by the source/sink calibration.",
        fontsize=7.5, color="#444444", wrap=True)

    path = HERE / "groin_timeline_and_hindcast.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    print(f"wrote {path}")

    print("\nphase summary")
    for a, b, label in ((1967, 1978, "BUILD"), (1978, 2004, "slow growth"),
                        (1984, 2004, "PERIOD 1"), (2004, 2023, "PERIOD 2")):
        ya = min(y for y in years if y >= a)
        yb = max(y for y in years if y <= b)
        delta = obs[yb] - obs[ya]
        print(f"  {label:<14} {ya}-{yb}  {delta:+7.1f} m  ({delta / (yb - ya):+.2f} m/yr)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
