#!/usr/bin/env python3
"""
The fixed-shift groin (the first version, the source/sink dipole) and the wave-driven groin side by side on straight coasts.

    python straight_coast_fixed_shift_vs_wave_driven_figure.py   ->  figures/straight_coast_fixed_shift_vs_wave_driven.png
    python straight_coast_fixed_shift_vs_wave_driven_figure.py --module fixed-shift
                                                 ->  ../5-fixed-shift/figures/straight_coast_fixed_shift.png
    python straight_coast_fixed_shift_vs_wave_driven_figure.py --module wave-driven
                                                 ->  ../6-wave-driven/figures/straight_coast_wave_driven.png

Left column module 1, right column module 6, the same rows and axes in each: the rule, the shoreline
around the groin over 50 years at -20 and +10 degrees, what each side has done by year 50 at every
orientation, and whether sand is conserved. --module draws one column alone, on its own axes. Reads the
tables straight_coast_test.py writes (run that first).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "1-run"))   # the test, in ../1-run/
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = HERE / "figures" / "straight_coast_fixed_shift_vs_wave_driven.png"
PAIR = ("dipole", "corrected")
ALONE = {                           # --module: the version drawn alone, and its own folder
    "fixed-shift": ("dipole", HERE.parent / "5-fixed-shift" / "figures" / "straight_coast_fixed_shift.png"),
    "wave-driven": ("corrected", HERE.parent / "6-wave-driven" / "figures" / "straight_coast_wave_driven.png"),
}
ANGLES = (-20, 10)                  # a coast turned against the drift, and one turned with it
CELLS = 12                          # cells either side of the groin drawn in the profiles
MAX_THETA = 30                      # the drift reverses beyond about +32 degrees, which swaps the sides
SKETCH = {                          # what each rule reads, and its labels on the two cells
    "dipole": dict(reads="reads nothing", up="+M", down="$-$M", drift=False),
    "corrected": dict(reads="reads the net drift, and\nreturns BRIE's leak\nnear the groin",
                      up="+f$_{block}$ Q$_{net}$", down="$-$f$_{block}$ Q$_{net}$", drift=True),
}
# -----------------------------------------------------------------------------


# Plan view of one rule: land above, ocean below, updrift on the right, the groin between two cells
def sketch(ax, name, spec, i, title=True):
    s = SKETCH[name]
    x = np.array([0, 4, 4, 5, 5, 6, 6, 10])
    y = np.array([0, 0, 1, 1, -1, -1, 0, 0])
    ax.fill_between(x, y, 3.5, color=C["ADDED_FILL"], lw=0)
    ax.fill_between(x, -4, y, color=C["WATER"], lw=0)
    ax.fill_between([5, 6], [-1] * 2, 0, color=spec["col"], alpha=0.45, lw=0)
    ax.fill_between([4, 5], 0, [1] * 2, color="white", alpha=0.75, lw=0)
    ax.plot([0, 10], [0, 0], color=INK_MUTED, lw=0.7, ls=(0, (3, 2)))
    ax.plot(x, y, color=INK, lw=1.2)
    for xc in (3, 4, 6, 7):
        ax.plot([xc, xc], [0.25, -0.25], color=INK_MUTED, lw=0.5)
    ax.plot([5, 5], [1.2, -2.2], color=C["GROIN"], lw=2.2, solid_capstyle="butt")
    ax.text(6.15, -0.5, s["up"], ha="left", va="center", fontsize=6.0, color=INK)
    ax.text(3.85, 0.5, s["down"], ha="right", va="center", fontsize=6.0, color=INK)
    if s["drift"]:
        ax.annotate("", xy=(5.15, -2.9), xytext=(9.6, -2.9), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3))
        ax.annotate("", xy=(1.0, -2.9), xytext=(4.85, -2.9),
                    arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=0.7, ls=(0, (2, 1.5))))
        ax.text(7.4, -3.15, "drift", ha="center", va="top", fontsize=6.5, color=INK)
        ax.text(2.9, -3.15, "1$-$f$_{block}$ passes", ha="center", va="top", fontsize=6.5, color=INK_MUTED)
    else:
        ax.annotate("", xy=(1.0, -2.9), xytext=(9.6, -2.9),
                    arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=0.7, ls=(0, (2, 1.5))))
        ax.text(7.4, -3.15, "drift (not read)", ha="center", va="top", fontsize=6.5, color=INK_MUTED)
    ax.text(5, 3.35, s["reads"], ha="center", va="top", fontsize=6.6, color=INK)
    ax.text(9.8, -1.8, "updrift", ha="right", va="top", fontsize=6.3, color=INK_MUTED, style="italic")
    ax.text(0.2, -1.8, "downdrift", ha="left", va="top", fontsize=6.3, color=INK_MUTED, style="italic")
    ax.set_xlim(0, 10)
    ax.set_ylim(-4, 3.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_color(INK_MUTED)
        side.set_linewidth(0.6)
    _letter_inside(ax, i)
    if title:
        ax.set_title(spec["label"], loc="center", color=spec["col"], fontweight="bold")


def tables(name):
    f = sc.MODULES[name]["folder"]
    base = sc.module_dir(sc.MODULES[name]) / "tables" / f"straight_coast_{f}"
    return pd.read_csv(f"{base}_by_year.csv"), pd.read_csv(f"{base}_profiles.csv")


def main(module="both"):
    apply_style()
    pair, out = (PAIR, OUT) if module == "both" else ((ALONE[module][0],), ALONE[module][1])
    data = {n: tables(n) for n in pair}
    fig, axes = plt.subplots(5, len(pair), figsize=figsize("double" if len(pair) == 2 else "single", height=9.4),
                             constrained_layout=True, gridspec_kw=dict(height_ratios=(0.9, 1, 1, 1, 1)),
                             squeeze=False)

    # profile and orientation ranges shared across the two columns
    prof = {a: pd.concat([data[n][1][(data[n][1].theta0_deg == a) & (data[n][1].domain_from_groin.abs() <= CELLS)]
                          for n in pair]) for a in ANGLES}
    y50 = {n: data[n][0][(data[n][0].year == 50) & (data[n][0].theta0_deg <= MAX_THETA)].sort_values("theta0_deg")
           for n in pair}
    k = 0
    for j, name in enumerate(pair):
        spec = sc.MODULES[name]
        by_year, profiles = data[name]
        sketch(axes[0, j], name, spec, k, title=len(pair) == 2)    # alone, the suptitle names it
        k += 1
        for r, a in enumerate(ANGLES, start=1):
            ax = axes[r, j]
            for c, yr in zip(sc.YEAR_RAMP, sc.PROFILE_YEARS):
                p = profiles[(profiles.theta0_deg == a) & (profiles.year == yr)
                             & (profiles.domain_from_groin.abs() <= CELLS)].sort_values("domain_from_groin")
                ax.plot(p.domain_from_groin, p.seaward_m, color=c, lw=1.3, drawstyle="steps-mid")
            ax.axvline(0.5, color=C["GROIN"], lw=1.6)      # the face between cells 0 and 1
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            lo, hi = prof[a].seaward_m.min(), prof[a].seaward_m.max()
            pad = 0.08 * (hi - lo)
            ax.set_ylim(lo - pad, hi + pad)
            ax.annotate("", xy=(1.5, 0.93), xytext=(CELLS - 0.5, 0.93), xycoords=("data", "axes fraction"),
                        arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.0))
            ax.text(CELLS / 2, 0.89, "drift", transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=6.5, color=INK_MUTED)
            end = by_year[(by_year.theta0_deg == a) & (by_year.year == 50)].iloc[0]
            ax.text(0.02, 0.05, f"year 50: total sand {end.reach_change_m:+,.0f} m", transform=ax.transAxes,
                    fontsize=6.8, color=INK)
            ax.set_xlabel("Cells from the groin (500 m); updrift →")
            if j == 0:
                ax.set_ylabel(f"Shoreline (m), + seaward\ncoast at {a:+d}°")
            else:
                ax.tick_params(labelleft=False)
            ax.grid(axis="y")
            open_frame(ax)
            _letter_inside(ax, k)
            ax.set_title(f"Coast turned {'against' if a < 0 else 'with'} the drift ({a:+d}°)", fontsize=8)
            k += 1

        # year 50 against orientation: each side, then the sand budget
        d = y50[name]
        ax = axes[3, j]
        ax.plot(d.theta0_deg, d.updrift_seaward_m, color=sc.MODULES[name]["col"], lw=1.6, marker="o", ms=2.5)
        ax.plot(d.theta0_deg, d.downdrift_seaward_m, color=sc.MODULES[name]["col"], lw=1.2, ls=(0, (3, 2)),
                marker="o", ms=2, mfc="white")
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_ylabel("Year 50 (m), + seaward" if j == 0 else "")
        ax.set_title("Each side of the groin at year 50", fontsize=8)
        ax = axes[4, j]
        ax.axhspan(-30, 30, color=C["BASE_FILL"], alpha=0.6, lw=0)
        ax.plot(d.theta0_deg, d.reach_change_m, color=sc.MODULES[name]["col"], lw=1.6, marker="o", ms=2.5)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_ylabel("Sand gained or lost\nby year 50 (m of shoreline)" if j == 0 else "")
        ax.set_title("Is sand conserved? (0 = yes)", fontsize=8)
        for ax in axes[3:, j]:
            ax.axvline(0, color=INK_MUTED, lw=0.5, ls=(0, (2, 2)))
            ax.set_xlabel("Coast orientation (°)")
            ax.grid(axis="y")
            open_frame(ax)
            _letter_inside(ax, k)
            k += 1
    for r in (3, 4):
        lo = min(axes[r, j].get_ylim()[0] for j in range(len(pair)))
        hi = max(axes[r, j].get_ylim()[1] for j in range(len(pair)))
        for j in range(len(pair)):
            axes[r, j].set_ylim(lo, hi)
        if len(pair) == 2:
            axes[r, 1].tick_params(labelleft=False)

    handles = ([Line2D([], [], color=c, lw=1.3, label=f"year {y}") for c, y in zip(sc.YEAR_RAMP, sc.PROFILE_YEARS)]
               + [Line2D([], [], color=INK, lw=1.6, label="updrift cell"),
                  Line2D([], [], color=INK, lw=1.2, ls=(0, (3, 2)), label="downdrift cell"),
                  Line2D([], [], color=C["GROIN"], lw=1.6, label="groin")])
    fig.legend(handles=handles, loc="outside lower center", ncol=8 if len(pair) == 2 else 4, frameon=False,
               fontsize=6.8)
    fig.suptitle("The fixed-shift groin against the wave-driven groin, on straight coasts" if len(pair) == 2
                 else f"The {sc.MODULES[pair[0]]['label'].lower()} groin on straight coasts")
    save(fig, out, dpi=250, close=True)

    e = {n: y50[n].set_index("theta0_deg") for n in pair}
    if module != "both":
        record_caption(out, alone_caption(pair[0], e[pair[0]]))
        print(f"wrote {out.relative_to(sc.REPO)}")
        return
    dp, cr = e["dipole"], e["corrected"]
    record_caption(OUT, (
        "**Tests:** how the fixed-shift groin (the first version, the source/sink dipole) and the current candidate "
        "(the wave-driven groin, the corrected drift groin) differ on coasts with nothing else going on. **How to read:** "
        "left column module 1, right column module 6; the same axes on each row. (a, f) The rule each year, "
        "in plan view (land above, ocean below, updrift on the right, red the groin): the fixed-shift groin moves a fixed "
        f"M = {sc.M:g} m from the downdrift cell to the updrift cell whatever the waves do; the wave-driven groin "
        f"holds back a fraction f$_{{block}}$ = {sc.BLOCKING_FRACTION} of the wave climate's net drift arriving "
        f"at the groin, and gives back what BRIE's row-scaled solve creates or destroys within "
        f"{sc.CORRECTION_HALF_WIDTH} cells of it. (b, c, g, h) The shoreline around the groin in years 1, 5, "
        "10, 25 and 50 (light to dark), on a coast turned against the drift (-20 degrees) and one turned with "
        "it (+10 degrees), with the total sand gained or lost across the 61-cell reach by year 50. (d, i) "
        "The updrift (solid) and downdrift (dashed) cells at year 50 against the coast's orientation. "
        "(e, j) The total change across the reach at year 50, in metres of shoreline summed over cells; grey "
        "band ±30 m, the pinned groin's range. Orientations stop at +30 degrees, before the drift reverses. "
        "BRIE's alongshore step alone under option A waves; no storms, no Barrier3D. **Shows:** the fixed-shift groin "
        "does not read the waves, so what it builds depends on BRIE, not the drift: at -20 degrees BRIE's solve "
        f"adds sand and both sides end seaward (updrift {dp.loc[-20, 'updrift_seaward_m']:+.0f} m, downdrift "
        f"{dp.loc[-20, 'downdrift_seaward_m']:+.0f} m); at +10 degrees it removes sand and both sides end "
        f"landward ({dp.loc[10, 'updrift_seaward_m']:+.0f} and {dp.loc[10, 'downdrift_seaward_m']:+.0f} m), "
        f"and the reach gains or loses up to {dp.reach_change_m.abs().max():,.0f} m. The wave-driven groin builds "
        f"the updrift side and cuts the downdrift side by the same amount at every orientation "
        f"({cr.loc[-20, 'updrift_seaward_m']:+.0f} / {cr.loc[-20, 'downdrift_seaward_m']:+.0f} m at -20 "
        f"degrees, {cr.loc[0, 'updrift_seaward_m']:+.0f} / {cr.loc[0, 'downdrift_seaward_m']:+.0f} m at 0), "
        "fading to nothing as the drift falls toward zero, and keeps the reach's sand within "
        f"{cr.reach_change_m.abs().max():.0f} m."))
    print(f"wrote {OUT.relative_to(sc.REPO)}")


# The caption of one version drawn alone (--module)
def alone_caption(name, d):
    how = ("**How to read:** (a) The rule each year, in plan view (land above, ocean below, updrift on the right, "
           "red the groin). (b, c) The shoreline around the groin in years 1, 5, 10, 25 and 50 (light to dark), on "
           "a coast turned against the drift (-20 degrees) and one turned with it (+10 degrees), with the total "
           "sand gained or lost across the 61-cell reach by year 50. (d) The updrift (solid) and downdrift "
           "(dashed) cells at year 50 against the coast's orientation. (e) The total change across the reach at "
           "year 50, in metres of shoreline summed over cells; grey band ±30 m, the pinned groin's range. "
           "Orientations stop at +30 degrees, before the drift reverses. BRIE's alongshore step alone under "
           "option A waves; no storms, no Barrier3D. The two versions side by side are in "
           "`straight_coast_fixed_shift_vs_wave_driven.png`. ")
    if name == "dipole":
        return ("**Tests:** what the fixed-shift groin (the first version, the source/sink dipole) does on coasts "
                f"with nothing else going on. It moves a fixed M = {sc.M:g} m from the downdrift cell to the "
                "updrift cell each year, whatever the waves do. " + how +
                "**Shows:** what it builds depends on BRIE, not the drift: at -20 degrees BRIE's solve adds sand "
                f"and both sides end seaward (updrift {d.loc[-20, 'updrift_seaward_m']:+.0f} m, downdrift "
                f"{d.loc[-20, 'downdrift_seaward_m']:+.0f} m); at +10 degrees it removes sand and both sides end "
                f"landward ({d.loc[10, 'updrift_seaward_m']:+.0f} and {d.loc[10, 'downdrift_seaward_m']:+.0f} m), "
                f"and the reach gains or loses up to {d.reach_change_m.abs().max():,.0f} m.")
    return ("**Tests:** what the wave-driven groin (the corrected drift groin, the current candidate) does on "
            f"coasts with nothing else going on. It holds back a fraction f$_{{block}}$ = {sc.BLOCKING_FRACTION} "
            "of the wave climate's net drift arriving at the groin, and gives back what BRIE's row-scaled solve "
            f"creates or destroys within {sc.CORRECTION_HALF_WIDTH} cells of it. " + how +
            "**Shows:** it builds the updrift side and cuts the downdrift side by the same amount at every "
            f"orientation ({d.loc[-20, 'updrift_seaward_m']:+.0f} / {d.loc[-20, 'downdrift_seaward_m']:+.0f} m "
            f"at -20 degrees, {d.loc[0, 'updrift_seaward_m']:+.0f} / {d.loc[0, 'downdrift_seaward_m']:+.0f} m at "
            "0), fading to nothing as the drift falls toward zero, and keeps the reach's sand within "
            f"{d.reach_change_m.abs().max():.0f} m.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--module", choices=("both", "fixed-shift", "wave-driven"), default="both")
    main(ap.parse_args().module)
