#!/usr/bin/env python3
"""
Which domains were eligible for correction, and how much each one got.

    python scripts/input_prep/7-source-sink/3-figures/plot_be_zones.py

Zones and the pass-0 and final fields per domain, both periods. Details: scripts/input_prep/7-source-sink/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import importlib.util
import pathlib
import re
import sys

import numpy as np

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
# The figure lives with the rest of the section 7 figures, in the data tree.
sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))
from site_layer import hat_source_sink as _be  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
FIG_DIR = _be.figures_dir()   # the default pair's (2026-09-18)
CONFIG = PROJECT_BASE_DIR / "scripts" / "site_layer" / "hatteras_site_config.py"

# THE PASS-0 FIELD, AND WHY IT IS THIS FILE AND NOT THE ONE BESIDE IT
PASS0_BACKUP = _be.PASS0_BACKUP

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

from site_layer.hat_figure_style import (                                   # noqa: E402
    apply_style, figsize, save, caption, town_bands, open_frame,
    DOMAIN_AXIS_LABEL, C, C_1984, C_1997, INK, INK_MUTED, _title)

_ROW = re.compile(r"^\s*(\d+):\s*([+-]?\d+\.?\d*),", re.M)

# The earlier period is the red of the house vintage pair, the later the blue.
WITHHELD_FILL = C["BASE_FILL"]
PERIODS = ((1984, "1984–2004", C_1984), (2004, "2004–2024", C_1997))
# -----------------------------------------------------------------------------


# One period's rates from a config file
def rates_from(path, period):
    path = pathlib.Path(path)
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} is missing. The lower panels split the final field into "
            f"the one-shot solve and what the iteration added, and the one-shot "
            f"half can only come from this lineage's pass-0 backup -- the file "
            f"the apply step wrote going INTO pass 1. A backup from another "
            f"lineage must not be substituted, and convergence_history.json "
            f"records only the per-pass RMSE, not the per-domain field. Restore "
            f"that file, or re-run the masked iteration from pass 0 and keep "
            f"every backup it writes, before drawing this figure.")
    text = path.read_text(encoding="utf-8")
    block = text.split("HATTERAS_BE_RATES_CALIBRATED")[1]
    segment = block.split(f"{period}:")[1]
    segment = segment[:segment.find("},")]
    return {int(m.group(1)): float(m.group(2)) for m in _ROW.finditer(segment)}


# be_zone_residual_fit.py, loaded as a module for its LOWESS helpers
def analysis_module():
    spec = importlib.util.spec_from_file_location(
        "_lowess",
        _HERE.parent.parent / "2-calibrate" / "be_zone_residual_fit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Run: the zones figure
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

    final = {p: rates_from(CONFIG, p) for p, _, _ in PERIODS}
    pass0 = {p: rates_from(PASS0_BACKUP, p) for p, _, _ in PERIODS}

    gis = np.arange(1, 91)
    apply_style()
    figure = plt.figure(figsize=figsize("double", height=6.2),
                        constrained_layout=True)
    grid = figure.add_gridspec(3, 1, height_ratios=[1.05, 1, 1])
    strip = figure.add_subplot(grid[0])
    bars = [figure.add_subplot(grid[1]), figure.add_subplot(grid[2])]

    # Top: eligibility, one row per period

    # The fill says what happened to the domain, not which zone it is in
    for row, (period, label, colour) in enumerate(PERIODS):
        y0 = (1 - row) * 1.0
        for d in gis:
            kw = dict(facecolor=WITHHELD_FILL, edgecolor="none")
            if d in locked:
                kw = dict(facecolor=C["BASE"], edgecolor="none")
            elif d in reserved:
                kw = dict(facecolor="none", edgecolor=C["BASE"], hatch="///",
                          linewidth=0.0)
            elif d in frozen[period]:
                kw = dict(facecolor=colour, edgecolor="none")
            strip.add_patch(plt.Rectangle((d - 0.5, y0), 1.0, 0.66, zorder=3,
                                          **kw))
        n = len([d for d in frozen[period] if d not in reserved | locked])
        strip.text(-1.2, y0 + 0.33, label, ha="right", va="center",
                   fontsize=8, color=colour)
        strip.text(91.2, y0 + 0.33, f"{n} correctable", ha="left", va="center",
                   fontsize=7, color=INK_MUTED)

    # Zone extents named once on a ruler under the strip
    short = getattr(module, "ZONE_DISPLAY_NAMES", {})
    for k, (name, (d0, d1, _mech)) in enumerate(zones.items()):
        strip.plot([d0 - 0.4, d1 + 0.4], [-0.14, -0.14], color=INK_MUTED,
                   linewidth=1.4, solid_capstyle="butt", zorder=3)
        strip.text((d0 + d1) / 2.0, -0.24 - 0.20 * (k % 2),
                   short.get(name, name), ha="center", va="top", fontsize=6.5,
                   color=INK_MUTED)

    strip.set_xlim(-10, 104)
    strip.set_ylim(-0.80, 2.60)      # headroom above the upper row for the key
    strip.set_yticks([])
    strip.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    for side in ("top", "right", "left"):
        strip.spines[side].set_visible(False)
    town_bands(strip, where="top", strip=0.07, fontsize=7)
    _title(strip, 0, "which domains were eligible")
    strip.legend(handles=[
        Patch(facecolor=C_1984, label="inside the zone set, 1984\u20132004"),
        Patch(facecolor=WITHHELD_FILL, label="withheld, left at zero"),
        Patch(facecolor=C_1997, label="inside the zone set, 2004\u20132024"),
        Patch(facecolor="none", edgecolor=C["BASE"], hatch="///",
              label="reserved for the groin"),
        Patch(facecolor=C["BASE"], label="boundary domain, solved separately")],
        loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=3, frameon=False,
        fontsize=7)

    # Middle / bottom: how much correction, and from which pass
    for i, (axis, (period, label, colour)) in enumerate(zip(bars, PERIODS)):
        p0 = np.array([pass0[period].get(d, 0.0) for d in gis])
        fin = np.array([final[period].get(d, 0.0) for d in gis])
        interior = (gis >= 2) & (gis <= 89)      # D1/D90 dwarf everything

        axis.bar(gis[interior], p0[interior], width=0.86, color=colour,
                 alpha=0.40, zorder=3, label="the one-shot solve")
        # Stacked in the SAME direction as pass 0, so the bar reads as a total rather than a difference
        axis.bar(gis[interior], (fin - p0)[interior], width=0.86,
                 bottom=p0[interior], color=colour, zorder=4,
                 label="what the further passes added")

        for d in gis[interior]:
            if d in reserved:
                axis.axvspan(d - 0.5, d + 0.5, facecolor="none",
                             edgecolor=C["BASE"], hatch="///", linewidth=0.0,
                             alpha=0.6, zorder=0)
            elif d not in frozen[period]:
                axis.axvspan(d - 0.5, d + 0.5, color=WITHHELD_FILL, alpha=0.7,
                             lw=0, zorder=0)

        axis.axhline(0.0, color=INK, linewidth=0.7, zorder=5)
        _title(axis, i + 1, label)
        axis.set_ylabel("background erosion\nrate (m/yr)")
        axis.set_xlim(-10, 104)
        axis.set_xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
        axis.grid(axis="y")
        open_frame(axis)
        # No village bands here: panel (a) above carries them
        axis.legend(loc="upper left", ncol=1, frameon=False, fontsize=7)

    bars[1].set_xlabel(DOMAIN_AXIS_LABEL)

    moved = {}
    for period, label, _c in PERIODS:
        p0 = np.array([pass0[period].get(d, 0.0) for d in gis])
        fin = np.array([final[period].get(d, 0.0) for d in gis])
        interior = (gis >= 2) & (gis <= 89)
        added = np.abs(fin - p0)[interior]
        moved[label] = (int(np.sum(added > 1e-9)),
                        added[added > 1e-9].mean() if np.any(added > 1e-9) else 0.0,
                        added.max())
    detail = "; ".join(
        f"{k}, {v[0]} domains moved, mean {v[1]:.2f} and at most {v[2]:.1f} m/yr"
        for k, v in moved.items())

    caption(figure, (
        "Which domains the source/sink calibration was allowed to correct, and "
        "how much each one received. A rate of zero in the field is ambiguous on "
        "its face -- it can mean 'no residual here' or 'this domain was never "
        "eligible', which are opposite claims. (a) separates them: a domain is "
        "coloured where it lay inside the zone set that was identified once from "
        "the first residual and then held for every pass, pale where it was "
        "withheld and stays at zero however large its residual, hatched at "
        "domains 5 to 7 where the Buxton groin owns the misfit, and dark grey at "
        "the two ends, which are boundary absorbers solved separately by "
        "buffer-cell reproduction. The physical zones are named on the ruler "
        "beneath. (b, c) split each final rate into the one-shot solve and what "
        "the further passes added on top, for the interior domains only -- the "
        "two ends carry rates about ten times larger and would flatten "
        "everything else. That split is the case for iterating at all: imposing "
        "X m/yr of background erosion does not move a domain's rate by X, "
        "because BRIE diffuses most of it alongshore, so the one-shot solve "
        "closes only 42 per cent of the misfit in the first period and 57 per "
        f"cent in the second ({detail}). Zone membership itself was not "
        "iterated: re-deriving it each pass would let less coherent features "
        "cross the threshold as the real ones were satisfied, and would let "
        "later passes correct the alongshore spillover of earlier ones, which "
        "never terminates. Domain 1 is at Cape Point and domain 90 at Pea "
        "Island, 500 m per domain."))

    path = FIG_DIR / "2-method" / "fig_be_zones_and_corrections.png"
    save(figure, path)
    plt.close(figure)

    print(f"wrote {path}")

    for period, label, _ in PERIODS:
        eligible = [d for d in frozen[period] if d not in reserved | locked]
        print(f"  {label}: {len(eligible)} correctable, "
              f"{88 - len(eligible)} interior domains withheld or reserved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
