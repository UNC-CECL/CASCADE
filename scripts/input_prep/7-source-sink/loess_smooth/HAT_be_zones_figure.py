#!/usr/bin/env python3
"""Which domains were eligible for correction, and how much each one got.

Two questions a reader of the calibrated field has to be able to answer, and
neither is visible in a table of 90 numbers:

    WHICH DOMAINS QUALIFIED
        A 0.0 in the field is ambiguous on its face -- it can mean "no residual
        here" or "this domain was never eligible". Those are opposite claims.
        The top panel separates them: coloured by physical zone where the
        domain was inside the frozen zone set, grey where it was withheld
        however large its residual, orange at D5-D7 where the groin owns the
        misfit, purple at the two locked ends.

    HOW MUCH CORRECTION IT GOT, AND FROM WHERE
        The lower panels split each final rate into the ONE-SHOT solve and what
        the ITERATION added on top. That split is the case for iterating at
        all: if the light bars were most of the height, the ordinary
        calibration was already converged and the extra passes bought nothing.
        They are not -- the iteration roughly doubles the field, because
        imposing X m/yr of background erosion moves a domain's rate by far less
        than X once BRIE has diffused it alongshore.

WHY ZONE MEMBERSHIP IS FIXED
    Zone identification is the scientific step -- this stretch has a real
    sediment-budget deficit, here is the process. Magnitude is arithmetic.
    Iterating both lets the arithmetic rewrite the science: each pass
    re-derives zones from a new residual, so as coherent features are satisfied
    less coherent ones cross the threshold, and since adding background erosion
    at one domain changes its neighbours' residuals, later passes start
    correcting the spillover of earlier ones. So zones are identified once and
    held (`FROZEN_ZONE_DOMAINS`).

Usage:
    python HAT_be_zones_figure.py

Reads  the live FROZEN_ZONE_DOMAINS / GROIN_RESERVED_DOMAINS / PHYSICAL_ZONES,
       the calibrated field from hatteras_site_config.py, and the pass-0 field
       from the masked iteration's first backup -- so the figure cannot drift
       from the calibration it documents.
Writes output/fig_be_zones_and_corrections.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import importlib.util
import pathlib
import re
import sys

import numpy as np

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
OUTPUT_DIR = _HERE.parent / "output"
CONFIG = PROJECT_BASE_DIR / "scripts" / "hatteras_site_config.py"
PASS0_BACKUP = (PROJECT_BASE_DIR / "scripts"
                / "hatteras_site_config_prebe_20260824_223143.py")

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

_ROW = re.compile(r"^\s*(\d+):\s*([+-]?\d+\.?\d*),", re.M)

PERIODS = ((1984, "1984–2004", "#1565C0"), (2004, "2004–2024", "#B71C1C"))
WITHHELD_COLOUR = "#C9C9C9"
RESERVED_COLOUR = "#FF8C00"
LOCKED_COLOUR = "#5E35B1"
ZONE_COLOURS = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974",
                "#64B5CD", "#937860"]


def rates_from(path, period):
    text = pathlib.Path(path).read_text(encoding="utf-8")
    block = text.split("HATTERAS_BE_RATES_CALIBRATED")[1]
    segment = block.split(f"{period}:")[1]
    segment = segment[:segment.find("},")]
    return {int(m.group(1)): float(m.group(2)) for m in _ROW.finditer(segment)}


def analysis_module():
    spec = importlib.util.spec_from_file_location(
        "_loess", _HERE.parent / "HAT_be_zone_LOESS_analysis.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    module = analysis_module()
    frozen = module.FROZEN_ZONE_DOMAINS
    reserved = set(module.GROIN_RESERVED_DOMAINS)
    locked = set(module.LOCKED_DOMAINS)
    zones = module.PHYSICAL_ZONES
    zone_colour = {name: ZONE_COLOURS[i % len(ZONE_COLOURS)]
                   for i, name in enumerate(zones)}

    final = {p: rates_from(CONFIG, p) for p, _, _ in PERIODS}
    pass0 = {p: rates_from(PASS0_BACKUP, p) for p, _, _ in PERIODS}

    gis = np.arange(1, 91)
    figure = plt.figure(figsize=(16, 9.6))
    grid = figure.add_gridspec(3, 1, height_ratios=[1.05, 1, 1], hspace=0.30)
    strip = figure.add_subplot(grid[0])
    bars = [figure.add_subplot(grid[1]), figure.add_subplot(grid[2])]

    # ---- TOP: eligibility, one row per period ----------------------------
    for row, (period, label, _) in enumerate(PERIODS):
        y0 = (1 - row) * 1.0
        for d in gis:
            if d in locked:
                face = LOCKED_COLOUR
            elif d in reserved:
                face = RESERVED_COLOUR
            elif d in frozen[period]:
                face = zone_colour[module.assign_physical_zone(d)]
            else:
                face = WITHHELD_COLOUR
            strip.add_patch(plt.Rectangle((d - 0.5, y0), 1.0, 0.66,
                                          facecolor=face, edgecolor="none"))
        n = len([d for d in frozen[period] if d not in reserved | locked])
        strip.text(-1.2, y0 + 0.33, label, ha="right", va="center",
                   fontsize=11.5, weight="bold")
        strip.text(91.2, y0 + 0.33, f"{n} correctable", ha="left", va="center",
                   fontsize=9.5, color="#444444")

    # zone names under the strip, so the colours are decodable without a
    # seven-entry legend competing with the status colours
    for name, (d0, d1, _) in zones.items():
        strip.plot([d0 - 0.5, d1 + 0.5], [-0.16, -0.16],
                   color=zone_colour[name], linewidth=3.5, solid_capstyle="butt")
        strip.text((d0 + d1) / 2.0, -0.30, name.replace(" / ", "/\n"),
                   ha="center", va="top", fontsize=7.5, color=zone_colour[name])

    strip.set_xlim(-9, 100)
    # Headroom above the 1984 row for the legend: at the previous limit it
    # was landing on top of the row it describes.
    strip.set_ylim(-0.95, 2.75)
    strip.set_yticks([])
    strip.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    strip.set_title("WHICH DOMAINS QUALIFIED\n"
                    "zone set identified once from the pass-0 residual, then "
                    "held fixed for every iteration pass", fontsize=12.5)
    strip.legend(handles=[
        Patch(facecolor="#777777", label="inside the frozen zone set "
                                         "(coloured by physical zone)"),
        Patch(facecolor=WITHHELD_COLOUR,
              label="withheld — outside the set, left at 0.0 however "
                    "large the residual"),
        Patch(facecolor=RESERVED_COLOUR, label="D5–D7 reserved for the groin"),
        Patch(facecolor=LOCKED_COLOUR, label="D1 / D90 locked "
                                             "(solved separately)")],
        loc="upper center", fontsize=8.5, ncol=2, framealpha=0.9)

    # ---- MIDDLE / BOTTOM: how much correction, and from which pass --------
    for axis, (period, label, colour) in zip(bars, PERIODS):
        p0 = np.array([pass0[period].get(d, 0.0) for d in gis])
        fin = np.array([final[period].get(d, 0.0) for d in gis])
        interior = (gis >= 2) & (gis <= 89)      # D1/D90 dwarf everything

        axis.bar(gis[interior], p0[interior], width=0.86, color=colour,
                 alpha=0.35, zorder=3, label="one-shot solve (pass 0)")
        # Stacked in the SAME direction as pass 0, so the bar reads as a total
        # rather than a difference; where the iteration reversed a sign the
        # segment simply crosses zero, which is itself worth seeing.
        axis.bar(gis[interior], (fin - p0)[interior], width=0.86,
                 bottom=p0[interior], color=colour, alpha=0.95, zorder=4,
                 label="added by iteration")

        for d in gis[interior]:
            if d in reserved:
                axis.axvspan(d - 0.5, d + 0.5, color=RESERVED_COLOUR,
                             alpha=0.20, zorder=0)
            elif d not in frozen[period]:
                axis.axvspan(d - 0.5, d + 0.5, color=WITHHELD_COLOUR,
                             alpha=0.45, zorder=0)

        axis.axhline(0.0, color="#333333", linewidth=0.9, zorder=5)
        moved = int(np.sum(np.abs(fin - p0)[interior] > 1e-9))
        added = np.abs(fin - p0)[interior]
        axis.set_title(
            f"{label}    final field, split by which pass produced it    "
            f"(iteration moved {moved} domains, mean "
            f"{added[added > 1e-9].mean():.2f}, max {added.max():.1f} m/yr)",
            fontsize=11.5, loc="left")
        axis.set_ylabel("background erosion\nrate (m/yr)", fontsize=10)
        axis.set_xlim(-9, 100)
        axis.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
        axis.grid(alpha=0.22, axis="y")
        axis.legend(fontsize=9, loc="lower right", ncol=2, framealpha=0.9)

    bars[1].set_xlabel("GIS domain  (south → north;  500 m per domain)",
                       fontsize=11.5)

    figure.suptitle("Source/sink calibration — eligible domains and the "
                    "correction each received", fontsize=14, y=0.985)
    figure.tight_layout(rect=(0, 0.075, 1, 0.962))
    figure.text(
        0.008, 0.010,
        "A 0.0 IN THE FIELD IS AMBIGUOUS on its face — it can mean 'no residual here' or 'never eligible'. The top panel separates them: grey domains "
        "were outside the frozen zone set and stay at 0.0 however large their residual, which is honest unexplained variance rather than a fitted "
        "constant. D1/D90 are boundary absorbers solved by buffer-cell reproduction and are excluded from the bars, where their ~10x rates would "
        "flatten everything else.\n"
        "WHY THE SPLIT MATTERS. Imposing X m/yr of background erosion does not move a domain's rate by X — BRIE diffuses most of it alongshore — so "
        "the one-shot solve closes only 42% (1984–2004) and 57% (2004–2024) of the misfit. The dark segments are what re-measuring and adding bought; "
        "that they are comparable to the light ones IS the case for iterating. Zone membership was NOT iterated: re-deriving it each pass let less "
        "coherent features cross the threshold as real ones were satisfied, and let later passes correct the alongshore spillover of earlier ones.",
        fontsize=7.5, color="#333333", wrap=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / "fig_be_zones_and_corrections.png"
    figure.savefig(path, dpi=170, facecolor="white")
    plt.close(figure)

    print(f"wrote {path}")
    for period, label, _ in PERIODS:
        eligible = [d for d in frozen[period] if d not in reserved | locked]
        print(f"  {label}: {len(eligible)} correctable, "
              f"{88 - len(eligible)} interior domains withheld or reserved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
