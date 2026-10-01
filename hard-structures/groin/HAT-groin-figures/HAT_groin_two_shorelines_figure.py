#!/usr/bin/env python3
"""
The groin in one picture: the two shorelines either side of it, and the gap between them.

    python HAT_groin_two_shorelines_figure.py

Reads the wet/dry change table (GIS 5 and 6 since 1967); writes
groin_two_shorelines.png (and .pdf) beside this file, caption in CAPTIONS.md.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
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
from site_layer.hat_figure_style import (apply_style, C, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save)

# --- CONFIG ------------------------------------------------------------------
WETDRY_TABLE = (GROIN_DIR / "HAT-groin-buxton-output" / "shoreline_position_output"
                / "Change_from_wetdry_1967_D2_D12.csv")

INSTALL_YEAR, LAST_REPAIR_YEAR, STORM_YEAR = 1969, 1996, 2003
UPDRIFT_GIS, DOWNDRIFT_GIS = 6, 5

# Sheltered side under test, unprotected side as baseline, the band between is the effect
UP_COLOR, DOWN_COLOR, GAP_FILL = C["ACCENT"], C["BASE"], C["ACCENT_FILL"]

PERIODS = ((1984, 2004, "hindcast period 1"), (2004, 2024, "hindcast period 2"))
# -----------------------------------------------------------------------------


# {year: (updrift, downdrift)} change since 1967, seaward-positive
def two_shorelines():
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


# The hindcast windows as a strip against the top edge, named once (README)
def period_strip(axis, spans, frac=0.062, shade="0.94"):
    for start, end, label in spans:
        axis.add_patch(plt.Rectangle(
            (start, 1.0 - frac), end - start, frac,
            transform=axis.get_xaxis_transform(), facecolor=shade,
            edgecolor="white", linewidth=0.8, zorder=0, clip_on=True))
        axis.text((start + end) / 2, 1.0 - frac / 2, label,
                  transform=axis.get_xaxis_transform(), ha="center",
                  va="center", fontsize=7, color=INK_MUTED, zorder=1)


# Endpoint rate of `series` between the first and last survey inside [start, end]
def endpoint_rate(years, series, start, end):
    inside = [y for y in years if start <= y <= end]
    if len(inside) < 2:
        return float("nan")
    first, last = min(inside), max(inside)
    order = list(years)
    return ((series[order.index(last)] - series[order.index(first)])
            / (last - first))


# Run: read the two shorelines, draw them and the gap, write the figure and caption
def main():
    apply_style()

    # The two shorelines and the gap between them
    data = two_shorelines()
    years = np.array(sorted(data))
    updrift = np.array([data[y][0] for y in years])
    downdrift = np.array([data[y][1] for y in years])
    gap = updrift - downdrift
    order = list(years)

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.50),
                           constrained_layout=True)

    # The gap IS the groin. Drawn first so the lines sit on top of it.
    ax.fill_between(years, downdrift, updrift, color=GAP_FILL, alpha=0.75,
                    linewidth=0, zorder=1,
                    label="the gap between them: the groin's effect")

    ax.plot(years, downdrift, marker="s", markersize=3.4, color=DOWN_COLOR,
            linewidth=1.6, zorder=4,
            label="downdrift of the structure, GIS domain {}".format(DOWNDRIFT_GIS))
    ax.plot(years, updrift, marker="o", markersize=3.4, color=UP_COLOR,
            linewidth=1.6, zorder=4,
            label="updrift, sheltered by the structure, GIS domain {}".format(UPDRIFT_GIS))

    # The structure's dated events
    for year, label in ((INSTALL_YEAR, "built"),
                        (LAST_REPAIR_YEAR, "last repair"),
                        (STORM_YEAR, "storm damage")):
        ax.axvline(year, color=INK_MUTED, linestyle=(0, (1, 2)), linewidth=0.8,
                   zorder=2)
        ax.text(year + 0.7, 0.02, "{} {}".format(label, year), rotation=90,
                fontsize=7, color=INK_MUTED, va="bottom", zorder=6,
                transform=ax.get_xaxis_transform())

    ax.axhline(0.0, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 3)),
               zorder=1)

    # Gap at its widest and at the last survey, measured off the series (README)
    widest = int(years[int(np.argmax(gap))])
    final = int(years[-1])
    for year, side in ([(widest, -1.0)] if widest == final
                       else [(widest, -1.0), (final, 1.0)]):
        up, down = data[year]
        ax.annotate("", xy=(year, up), xytext=(year, down),
                    arrowprops=dict(arrowstyle="<->", color=UP_COLOR, lw=1.0,
                                    shrinkA=0, shrinkB=0), zorder=5)
        ax.text(year + side, (up + down) / 2, "{:.0f} m".format(up - down),
                fontsize=7.5, color=UP_COLOR, va="center", zorder=6,
                ha="right" if side < 0 else "left")

    # Axes, strip and legend
    ax.set_xlim(years[0] - 1, years[-1] + 5)
    period_strip(ax, PERIODS)
    ax.set_xlabel("year")
    ax.set_ylabel("shoreline change since 1967 (m)\nnegative is retreat")
    ax.set_title("Shoreline change either side of the Buxton groin", loc="left")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)
    fig.legend(loc="outside lower center", ncol=3, frameon=False)

    # Caption, every number computed here
    first = int(years[0])
    caption(fig,
            "Shoreline change at the two Barrier3D domains flanking the Buxton "
            "groin field, against a fixed {first} datum, from {n} dated wet/dry "
            "surveys extracted by the GIS analysis. Plotted seaward-positive, "
            "so both lines falling means the whole reach is eroding; the source "
            "table is landward-positive and is negated here, which is why the "
            "timeline figure's fillet curve rises where these fall. The "
            "sheltered updrift line falls more slowly, so the gap opens, and "
            "that gap is the groin's entire effect. It is RELATIVE protection, "
            "not new beach: at the end of period 1 the downdrift side had "
            "retreated {dr04:.0f} m since {first} while the sheltered side had "
            "retreated {ur04:.0f} m, a {gw04:.1f} m gap, and by {last} it had "
            "closed to {gl:.0f} m. The annotated widths are measured off the "
            "plotted series, which puts the widest gap in the record at {wide} "
            "({gw:.0f} m) rather than at 2004: the {wide} survey is a "
            "single-survey spike on the downdrift side, and the figure "
            "annotates where the data actually peaks. For the hindcast, the "
            "gap is still opening through period 1 at {r1:+.1f} m/yr, so the "
            "module has something to reproduce, and it supplies about a third "
            "of that at the chosen parameters (recorded from the period-1 sweep "
            "when this figure was built, not computed here). Through period 2 "
            "the gap closes at {r2:+.1f} m/yr as the unprotected side catches "
            "up, which the module cannot do: trapping is bounded at zero, so it "
            "can stop adding sand but never remove it, and that closure is "
            "carried by the source/sink calibration. Both period rates are "
            "endpoint differences between the surveys bracketing the window, "
            "not regressions."
            .format(first=first, last=final, n=len(years), wide=widest,
                    dr04=-downdrift[order.index(2004)],
                    ur04=-updrift[order.index(2004)],
                    gw04=gap[order.index(2004)],
                    gw=gap[order.index(widest)], gl=gap[-1],
                    r1=endpoint_rate(years, gap, 1984, 2004),
                    r2=endpoint_rate(years, gap, 2004, final)))

    # Write the figure and print the gap at four surveys
    written = save(fig, HERE / "groin_two_shorelines.png", close=True)
    for path in written:
        print("wrote {}".format(path))
    for year in (1967, 1978, 2004, 2023):
        if year in data:
            up, down = data[year]
            print("  {}:  D6 {:+7.1f}   D5 {:+7.1f}   gap {:6.1f} m"
                  .format(year, up, down, up - down))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
