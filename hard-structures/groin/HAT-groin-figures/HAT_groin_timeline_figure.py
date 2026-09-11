#!/usr/bin/env python3
"""The Buxton groin's observed history, and how it maps onto the hindcast.

One figure answering three questions that kept getting tangled:

    WHEN did the groin affect the shoreline?
        Continuously from 1969, but in two opposite regimes. It TRAPPED through
        2004 and has been RELEASING since. The turning point coincides with the
        2003 storm damage, not with anything in the model.

    IS THERE A DIFFERENT SIGNAL IN EACH HINDCAST PERIOD?
        Yes, and they have opposite sign. Period 1 (1984-2004) is still
        accumulating; period 2 (2004-2024) is draining. Both rates are computed
        here and quoted in the caption.

    HOW DOES THE MODULE REPRESENT THAT?
        `GroinCallback` carries an absolute calendar schedule -- install 1969,
        deterioration onset 1996 (last repair), linear ramp to 2003 (storm),
        then hold at M*f. The lower panel draws that schedule against the
        observations, which is the clearest way to see that the module's
        built-in timeline already matches the measured history, and where it
        cannot follow.

WHAT THE MODULE CANNOT DO, DRAWN RATHER THAN FOOTNOTED
    Trapping is bounded at >= 0, so the groin can stop adding sand but cannot
    actively drain the fillet. Period 2's observed release is therefore outside
    what the parameterisation can produce at any (M, f), and the lower panel
    hatches that region so the limitation is visible next to the data.

DATA
    Fillet is x_s[D5] - x_s[D6] against a fixed 1967 datum, from
    `Change_from_wetdry_1967_D2_D12.csv` -- 24 dated wet/dry surveys produced by
    the GIS analysis in HAT-groin-gis-analysis. Landward-positive, so a rising
    curve means the updrift side is holding while the downdrift side retreats,
    which is what a groin builds.

STYLE, 2026-09-11
    Drawn under the project house style (`scripts/hat_figure_style.py`). Three
    things changed beyond type and colour. The canvas is a 190 mm column
    instead of 14 in, so its type survives a page. The phase washes are gone:
    four overlapping shades on one panel (three phases plus two hindcast
    windows) could not be told apart, so the phases are named with their rates
    in the caption and the hindcast windows are a strip against the top edge.
    And the title sentences and the footnote paragraph are in CAPTIONS.md
    beside the image, with every number in them computed here rather than
    written in.

Usage:
    python HAT_groin_timeline_figure.py

Writes groin_timeline_and_hindcast.png (and .pdf) beside this file.

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
from hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)

WETDRY_TABLE = (GROIN_DIR / "HAT-groin-buxton-output" / "shoreline_position_output"
                / "Change_from_wetdry_1967_D2_D12.csv")

# The structure's documented history. These are the values GroinCallback is
# configured with, not fitted quantities.
INSTALL_YEAR = 1969
LAST_REPAIR_YEAR = 1996
STORM_YEAR = 2003

# The chosen parameters. M is an EFFECTIVE, grid-specific rate -- see the
# accompanying GROIN_PLAN.md for why it is not a sediment flux.
CHOSEN_M, CHOSEN_F = 60.0, 0.6

PERIODS = ((1984, 2004, "hindcast period 1"), (2004, 2024, "hindcast period 2"))
UPDRIFT_GIS, DOWNDRIFT_GIS = 6, 5

# The schedule is the thing under test; the observations are the record it is
# read against. The RdBu vintage pair is not used here -- nothing on this
# figure is a 1984/1997 vintage.
SCHEDULE_COLOR, SCHEDULE_FILL = C["ACCENT"], C["ACCENT_FILL"]


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


def retreat_at(year):
    """(updrift, downdrift) retreat since 1967, positive landward, one survey.

    The caption's "175 m against 25 m" pair: read from the table rather than
    written into the prose. The column is matched on its survey year, not
    built by format -- the real names carry a trailing unit ("_m") and two of
    them a month ("1972_july"), so a formatted name misses."""
    frame = pd.read_csv(WETDRY_TABLE).set_index("Domain_ID")
    pattern = re.compile(r"change_from_wetdry_1967_wetdry_{}(?:_|$)".format(year))
    column = next(c for c in frame.columns if pattern.match(c))
    return (float(frame.loc[UPDRIFT_GIS, column]),
            float(frame.loc[DOWNDRIFT_GIS, column]))


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


def phase(obs, start, end):
    """(first survey, last survey, change, rate) inside [start, end]."""
    years = sorted(y for y in obs if start <= y <= end)
    first, last = years[0], years[-1]
    delta = obs[last] - obs[first]
    return first, last, delta, delta / (last - first)


def period_strip(axis, spans, frac=0.07, shade="0.94"):
    """The hindcast windows as a band against the top edge, named once.

    A strip rather than a full-height wash, for the same reason the phase
    washes were dropped: several overlapping greys on one panel read as one.
    The twin of this helper is in HAT_groin_two_shorelines_figure.py."""
    for start, end, label in spans:
        axis.add_patch(plt.Rectangle(
            (start, 1.0 - frac), end - start, frac,
            transform=axis.get_xaxis_transform(), facecolor=shade,
            edgecolor="white", linewidth=0.8, zorder=0, clip_on=True))
        axis.text((start + end) / 2, 1.0 - frac / 2, label,
                  transform=axis.get_xaxis_transform(), ha="center",
                  va="center", fontsize=7, color=INK_MUTED, zorder=1)


def main():
    apply_style()

    obs = observed_fillet()
    years = np.array(sorted(obs))
    values = np.array([obs[y] for y in years])

    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=figsize("double", aspect=0.78), sharex=True,
        gridspec_kw=dict(height_ratios=[3, 2]), constrained_layout=True)

    # ---- the structure's dated events, on both panels -------------------
    for axis in (top, bottom):
        for year in (INSTALL_YEAR, LAST_REPAIR_YEAR, STORM_YEAR):
            axis.axvline(year, color=INK_MUTED, linestyle=(0, (1, 2)),
                         linewidth=0.8, zorder=2)

    # ---- (a) the observed fillet ----------------------------------------
    top.plot(years, values, marker="o", markersize=3.4, color=INK,
             linewidth=1.6, zorder=5, label="observed fillet, wet/dry surveys")

    for year, label in ((INSTALL_YEAR, "built"),
                        (LAST_REPAIR_YEAR, "last repair"),
                        (STORM_YEAR, "storm damage")):
        # Clear of the curve: at 0.025 the install label ran through the 1970
        # survey, which reads as a label for the marker.
        top.text(year + 0.7, 0.10, "{} {}".format(label, year), rotation=90,
                 fontsize=7, color=INK_MUTED, va="bottom", zorder=6,
                 transform=top.get_xaxis_transform())

    top.set_ylabel("fillet (m)\nhow much less the sheltered side\nhas retreated"
                   " since 1967")
    top.grid(axis="y")
    top.set_axisbelow(True)
    open_frame(top)
    _title(top, 0, "the measured history")

    # ---- (b) what the module actually applies ----------------------------
    span = np.arange(1967, 2025)
    m_eff = np.array([effective_trapping(y) for y in span])
    bottom.plot(span, m_eff, color=SCHEDULE_COLOR, linewidth=1.8, zorder=5,
                label="trapping rate applied by the module,"
                      " M {:g} m/yr falling to {:g}".format(
                          CHOSEN_M, CHOSEN_M * CHOSEN_F))
    bottom.fill_between(span, 0, m_eff, color=SCHEDULE_FILL, alpha=0.45,
                        linewidth=0, zorder=3)

    # The floor the module cannot go below, hatched rather than washed so it
    # cannot be mistaken for another shaded period.
    bottom.axhline(0.0, color=INK, linewidth=0.8, zorder=4)
    bottom.fill_between([2004, 2024], -14, 0, facecolor="none",
                        edgecolor=C["BASE"], hatch="////", linewidth=0.0,
                        zorder=2)
    bottom.text(2014, -7.2, "not reachable: trapping is bounded at zero",
                ha="center", va="center", fontsize=7, color=INK_MUTED,
                zorder=6,
                bbox=dict(facecolor="white", alpha=0.9, edgecolor="none",
                          boxstyle="square,pad=0.25"))

    bottom.annotate("full trapping", xy=(1980, CHOSEN_M),
                    xytext=(1971, CHOSEN_M * 0.70), fontsize=7.5,
                    color=INK_MUTED,
                    arrowprops=dict(arrowstyle="->", color=INK_MUTED, lw=0.8))
    bottom.annotate("ramp to the floor, M x f",
                    xy=(STORM_YEAR, CHOSEN_M * CHOSEN_F),
                    xytext=(2006, CHOSEN_M * 0.72), fontsize=7.5,
                    color=INK_MUTED,
                    arrowprops=dict(arrowstyle="->", color=INK_MUTED, lw=0.8))

    bottom.set_ylim(-14, CHOSEN_M * 1.15)
    bottom.set_ylabel("effective trapping rate (m/yr)")
    bottom.set_xlabel("year")
    bottom.grid(axis="y")
    bottom.set_axisbelow(True)
    open_frame(bottom)
    _title(bottom, 1, "the module's calendar schedule")

    bottom.set_xlim(years[0] - 1, 2025)
    period_strip(top, PERIODS)
    fig.legend(loc="outside lower center", ncol=2, frameon=False)

    build = phase(obs, 1967, 1978)
    slow = phase(obs, 1978, 2004)
    one = phase(obs, 1984, 2004)
    two = phase(obs, 2004, int(years[-1]))
    up04, down04 = retreat_at(2004)
    caption(fig,
            "(a) The Buxton groin's fillet against a fixed 1967 datum, from "
            "{n} dated wet/dry surveys: the downdrift shoreline position minus "
            "the updrift one, landward-positive, so a rising curve means the "
            "sheltered side is holding while the exposed side retreats. It is a "
            "RELATIVE measure, not accumulation. By 2004 the downdrift side had "
            "retreated {d:.0f} m since 1967 while the sheltered updrift side "
            "had retreated {u:.0f} m, and that {g:.1f} m difference is the "
            "fillet, not {g:.1f} m of new beach. The record runs in three "
            "phases: a build at {rb:+.1f} m/yr to {b1}, slow growth at "
            "{rs:+.1f} m/yr to 2004, then release at {rr:+.1f} m/yr. The "
            "post-2004 collapse is about 85% the updrift side eroding once the "
            "structure failed, not impounded sand draining downdrift (recorded "
            "from the per-domain attribution when this figure was built, not "
            "computed here). Across the two hindcast windows the signal has "
            "opposite sign: {p1a}-{p1b} gains {d1:+.0f} m ({r1:+.2f} m/yr), "
            "{p2a}-{p2b} loses {d2:+.0f} m ({r2:+.2f} m/yr); all rates are "
            "endpoint differences between the surveys bracketing the window, "
            "not regressions. (b) The trapping rate the module applies, from "
            "the same absolute calendar the structure has: full rate from the "
            "1969 install, deterioration from the 1996 last repair, linear to "
            "the 2003 storm, then held at M x f. The observed turning point is "
            "2004, immediately after that storm, so the built-in ramp "
            "reproduces the timing with no period-specific configuration. The "
            "hatched region is what the parameterisation cannot reach: the rate "
            "is bounded at zero, so the module can stop adding sand but never "
            "drain the fillet, which puts period 2's release outside it at any "
            "(M, f). Period 1 is the only window in the hindcast where the "
            "groin does something the module can reproduce; period 2's closure "
            "is carried by the source/sink calibration."
            .format(n=len(years), d=down04, u=up04, g=down04 - up04,
                    rb=build[3], b1=build[1], rs=slow[3], rr=two[3],
                    p1a=one[0], p1b=one[1], d1=one[2], r1=one[3],
                    p2a=two[0], p2b=two[1], d2=two[2], r2=two[3]))

    written = save(fig, HERE / "groin_timeline_and_hindcast.png", close=True)
    for path in written:
        print("wrote {}".format(path))

    print("\nphase summary")
    for label, a, b in (("build", 1967, 1978), ("slow growth", 1978, 2004),
                        ("period 1", 1984, 2004),
                        ("period 2", 2004, int(years[-1]))):
        first, last, delta, rate = phase(obs, a, b)
        print("  {:<14} {}-{}  {:+7.1f} m  ({:+.2f} m/yr)"
              .format(label, first, last, delta, rate))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
