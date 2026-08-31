#!/usr/bin/env python3
"""The groin in one picture: two shorelines, and the gap between them.

WHY THIS EXISTS ALONGSIDE THE TIMELINE FIGURE
    The fillet is a DERIVED quantity -- a difference between two domains -- so
    it needs a sentence of explanation before anyone can read it, and it invites
    the wrong reading ("150 m of new beach") when it actually means "150 m of
    relative position".

    Plotting the two shorelines themselves needs no explanation. Both retreat.
    One retreats far less. The gap between the lines IS the groin's effect, and
    its width over time is the whole story:

        widening gap  ->  the groin is trapping
        closing gap   ->  the groin has stopped

    Everything the timeline figure says is visible here without a definition.

WHAT IT IS FOR
    Deciding how to apply the module. The groin's job in the model is to hold
    the updrift line above the downdrift one. Reading the gap tells you directly
    what the module has to do in each hindcast period, and where it cannot help.

DATA
    Shoreline change since a fixed 1967 datum at the two domains flanking the
    structure -- D5 downdrift, D6 updrift -- from 24 dated wet/dry surveys
    (`Change_from_wetdry_1967_D2_D12.csv`, GIS analysis).

    Plotted SEAWARD-POSITIVE so the lines fall as the shoreline retreats, which
    is how a reader expects an eroding coast to look. The source table is
    landward-positive, so it is negated here; the timeline figure keeps the
    source convention, which is why its curve rises where these fall.

Usage:
    python HAT_groin_two_shorelines_figure.py

Writes groin_two_shorelines.png beside this file.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import pathlib
import re

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
WETDRY_TABLE = (HERE / "HAT-groin-test-output" / "shoreline_position_output"
                / "Change_from_wetdry_1967_D2_D12.csv")

INSTALL_YEAR, LAST_REPAIR_YEAR, STORM_YEAR = 1969, 1996, 2003
UPDRIFT_GIS, DOWNDRIFT_GIS = 6, 5

UP_COLOR, DOWN_COLOR = "#B71C1C", "#1565C0"


def two_shorelines():
    """{year: (updrift, downdrift)} change since 1967, SEAWARD-positive."""
    frame = pd.read_csv(WETDRY_TABLE).set_index("Domain_ID")
    out = {}
    for column in frame.columns:
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", column)
        if not match:
            continue
        up, down = frame.loc[UPDRIFT_GIS, column], frame.loc[DOWNDRIFT_GIS, column]
        if pd.isna(up) or pd.isna(down):
            continue
        # Source table is landward-positive; negate so retreat plots downward.
        out[int(match.group(1))] = (-float(up), -float(down))
    out.setdefault(1967, (0.0, 0.0))
    return dict(sorted(out.items()))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = two_shorelines()
    years = np.array(sorted(data))
    updrift = np.array([data[y][0] for y in years])
    downdrift = np.array([data[y][1] for y in years])

    figure, axis = plt.subplots(figsize=(13, 7.5))

    # The gap IS the groin. Drawn first so the lines sit on top of it.
    axis.fill_between(years, downdrift, updrift, color="#FFC107", alpha=0.35,
                      zorder=1, label="the GAP = the groin's effect")

    axis.plot(years, updrift, marker="o", markersize=5, color=UP_COLOR,
              linewidth=2.6, zorder=4,
              label=f"D{UPDRIFT_GIS} UPDRIFT  (behind the groin -- protected)")
    axis.plot(years, downdrift, marker="s", markersize=5, color=DOWN_COLOR,
              linewidth=2.6, zorder=4,
              label=f"D{DOWNDRIFT_GIS} DOWNDRIFT  (starved of sand)")

    for year, colour, label in ((INSTALL_YEAR, "#333333", "groin built 1969"),
                                (LAST_REPAIR_YEAR, "#E65100", "last repair 1996"),
                                (STORM_YEAR, "#B71C1C", "storm damage 2003")):
        axis.axvline(year, color=colour, linestyle=":", linewidth=1.6, zorder=2)
        axis.text(year + 0.6, 0.02, label, rotation=90, fontsize=8.5,
                  color=colour, va="bottom", transform=axis.get_xaxis_transform())

    # The two hindcast periods, and what the gap does in each.
    for (start, end), text, colour in (
            ((1984, 2004), "PERIOD 1\ngap still WIDENING\ngroin working",
             "#1B5E20"),
            ((2004, 2024), "PERIOD 2\ngap CLOSING\ngroin failed", "#B71C1C")):
        axis.axvspan(start, end, color=colour, alpha=0.05, zorder=0)
        axis.text((start + end) / 2, 0.96, text, ha="center", va="top",
                  fontsize=10.5, weight="bold", color=colour,
                  transform=axis.get_xaxis_transform(), zorder=6,
                  bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                            edgecolor=colour, alpha=0.9))

    # Annotate the gap at its widest and at the end.
    for year, note in ((2004, "widest gap: 150 m"), (2023, "gap now 74 m")):
        if year not in data:
            continue
        up, down = data[year]
        axis.annotate("", xy=(year, up), xytext=(year, down),
                      arrowprops=dict(arrowstyle="<->", color="#F57F17", lw=2.2),
                      zorder=5)
        axis.text(year + 1.0, (up + down) / 2, note, fontsize=9.5,
                  color="#F57F17", weight="bold", va="center", zorder=6)

    axis.axhline(0.0, color="#999999", linewidth=1.0, zorder=1)
    axis.set_xlabel("year")
    axis.set_ylabel("shoreline change since 1967 (m)   [down = eroding]")
    axis.set_title(
        "The Buxton groin in one picture: both shorelines eroded, "
        "the protected one eroded far less\n"
        "The gap between the lines is what the groin module has to reproduce",
        fontsize=12.5)
    axis.grid(alpha=0.25)
    axis.legend(loc="lower left", fontsize=9.5)

    figure.tight_layout(rect=(0, 0.085, 1, 1))
    figure.text(
        0.01, 0.012,
        "HOW TO READ IT.  Both lines fall: the whole reach is eroding. The red line falls more slowly because the groin shelters "
        "it, so the gap opens -- that gap is the groin's entire effect, and it is RELATIVE protection, not new beach.\n"
        "FOR THE HINDCAST.  In PERIOD 1 the gap is still opening (+2.7 m/yr), so the module has something to reproduce and does "
        "about 33% of it at M=60, f=0.6. In PERIOD 2 the gap closes (-4.0 m/yr) as the unprotected updrift side catches up -- the "
        "module CANNOT do this, because trapping is bounded at zero and it can stop adding sand but never remove it. That closure "
        "is carried by the source/sink calibration.",
        fontsize=8, color="#333333", wrap=True)

    path = HERE / "groin_two_shorelines.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    print(f"wrote {path}")
    for year in (1967, 1978, 2004, 2023):
        if year in data:
            up, down = data[year]
            print(f"  {year}:  D6 {up:+7.1f}   D5 {down:+7.1f}   gap {up - down:6.1f} m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
