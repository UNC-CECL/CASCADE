#!/usr/bin/env python3
"""Intersects the two periods' sweep surfaces to fit M and f together.

WHY THIS IS A SEPARATE STEP
    Neither period identifies both knobs on its own. Cumulative trapping per
    unit M, measured from `cascade.groin.GroinCallback` with the documented
    1969 install / 1996 onset / 7-year ramp schedule:

        1984-2004    16 + 4f      f moves this only 16.0 -> 20.0, so holding
                                  cumulative trapping fixed while sliding f
                                  across its whole range needs just a 25%
                                  change in M. f is nearly free here.
        2004-2024    20f          the run sits entirely past the 2003 ramp,
                                  so M and f enter only as their product.

    One period gives one constraint on two unknowns. Two periods intersect.
    This script does that intersection, per source/sink preset.

HOW be1 IS TREATED
    be1 is a NUISANCE parameter, not a fitted result: it is swept in period 1
    under edgeBE and absent everywhere else. For each (M, f) the period-1
    contribution is MINIMISED over be1 -- a profile likelihood -- and the
    winning be1 is reported alongside. Summing over be1, or fixing it at one
    value, would both charge the groin for background erosion the model was
    free to place elsewhere.

WHAT THE ANSWER WILL LOOK LIKE, AND WHY
    Period 2's observed D6 - D5 is negative (-2.47 m/yr). The source/sink pair
    adds -M updrift and +M downdrift, so the modelled differential is
    non-negative at any M >= 0 and no cell can reach that target. Period 2
    therefore contributes a monotone penalty in M*f, pushing the joint
    solution toward f = 0, and period 1 then sets M through M*(16 + 4f).

    Algebraically, with A the period-1 constraint and B the period-2 one:

        f = 4*B / (5*A - B)

    B pinned near zero gives f near zero. That is a RESULT -- the structure
    stopped trapping after the 2003 storm -- but it is reached by railing to
    a grid edge, so this script flags every fitted value that lands on a
    bound rather than in the interior. A railed value is a bound, not an
    optimum, and must not be quoted as a fitted parameter.

Usage:
    python HAT_groin_joint_fit.py [--preset edgeBE] [--no-figures]

Reads   output/groin_sweep/<period>_<preset>/sweep_results.csv  (all four)
Writes  output/groin_sweep/joint_fit.json    fitted (M, f, be1) per preset
        output/groin_sweep/joint_fit.csv     the full joint surface
        output/groin_sweep/figures/joint_<preset>_surface.png
        output/groin_sweep/joint_<preset>_constraints.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
# The guard below is what makes a future move fail here, loudly, instead of
# resolving to scripts/scripts and surfacing as a missing data file several
# imports deeper.
PROJECT_BASE_DIR = _HERE.parents[3]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/hatteras_ms/groin-sweep/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from HAT_groin_sweep_config import (  # noqa: E402
    END_YEAR,
    F_VALUES,
    M_VALUES,
    OBSERVED_DIFFERENTIAL,
    PERIODS,
    PERIOD_DIFFERENTIAL_IS_REACHABLE,
    PRESETS,
    sweep_output_dir,
)

OUTPUT_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep"
# Figures live in a subdirectory; joint_fit.json does NOT. That file is read by
# HAT_run_all.py (stage 6 takes its fitted M and f from it) and by
# HAT_be_zone_LOESS_analysis.py (which uses it to find a groin-aware base run),
# both of which pin the top-level path. Moving it would break stage 6 silently.
FIGURE_DIR = OUTPUT_DIR / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
JOINT_JSON = OUTPUT_DIR / "joint_fit.json"
JOINT_CSV = OUTPUT_DIR / "joint_fit.csv"

MODEL_COLOR = "#FF8C00"
GROIN_COLOR = "#B71C1C"
OBSERVED_COLOR = "#1A1A1A"


def load_period(period, preset):
    """Loads one sweep's scored results.

    Returns:
        A DataFrame, or None if that sweep has not produced a CSV yet.
    """
    path = sweep_output_dir(period, preset) / "sweep_results.csv"
    if not path.exists():
        return None
    frame = pd.read_csv(path)
    return frame[frame["differential_err"].notna()].copy()


def period_surface(frame, period):
    """Reduces one period's rows to a value per (M, f) cell.

    `err` is the FILLET-SIZE error where the sweep recorded one, the
    differential error otherwise. The fillet saturates, so its level
    carries the information about M and its slope -- which is what the
    differential measures -- carries almost none.

    be1 is profiled out: for each (M, f) the best-scoring be1 is kept, and
    which one it was is carried along so the fitted background erosion can be
    reported with the fitted groin.

    The M = 0 baseline carries no f, so it is broadcast across every f -- with
    no groin attached, all f give the identical run, and leaving the row at a
    single f would punch a hole in the surface at M = 0.

    Returns:
        A DataFrame indexed by (M, fraction) with columns err, be1,
        differential.
    """
    rows = []
    zero = frame[frame["M"] == 0]
    for _, row in zero.iterrows():
        for fraction in F_VALUES:
            rows.append(dict(M=0.0, fraction=fraction,
                             err=row.get("fillet_err",
                                         row["differential_err"]),
                             be1=row.get("be1"),
                             differential=row["differential_m_yr"]))
    for _, row in frame[frame["M"] > 0].iterrows():
        rows.append(dict(M=row["M"], fraction=row["fraction"],
                         err=row.get("fillet_err",
                                     row["differential_err"]),
                         be1=row.get("be1"),
                         differential=row["differential_m_yr"]))

    expanded = pd.DataFrame(rows)
    # Profile out be1: keep the best-scoring one per (M, f).
    best = (expanded.sort_values("err")
            .drop_duplicates(subset=["M", "fraction"], keep="first")
            .set_index(["M", "fraction"]))
    return best.rename(columns={
        "err": f"err_{period}", "be1": f"be1_{period}",
        "differential": f"differential_{period}"})


def joint_fit(preset):
    """Builds the joint surface for one preset and picks its best cell.

    Returns:
        (surface, fit) where surface is a DataFrame over (M, f) and fit is a
        dict describing the winning cell, or (None, reason) if a period's
        sweep is missing.
    """
    surfaces = {}
    for period in PERIODS:
        frame = load_period(period, preset)
        if frame is None or frame.empty:
            return None, (f"no scored results for {period}-{END_YEAR[period]} "
                          f"{preset}; run its sweep first")
        surfaces[period] = period_surface(frame, period)

    surface = surfaces[PERIODS[0]].join(surfaces[PERIODS[1]], how="inner")
    if surface.empty:
        return None, (f"the two {preset} sweeps share no (M, f) cells -- one "
                      f"of them is incomplete")

    # RANKED ON PERIOD 1 ALONE. Cumulative trapping is M*(16 + 4f) in
    # 1984-2004 but 20*M*f in 2004-2024, so period 2 sees only the PRODUCT
    # M*f and cannot separate the two parameters. Every bit of information
    # distinguishing M from f lives in period 1, the one window straddling
    # the 1996-2003 ramp.
    #
    # Summing both legs let a window that cannot identify the parameters pull
    # them anyway, and in the wrong direction: holding period-1 trapping
    # fixed, driving f from 0.9 to 0 forces M from 50 to 61.3 m/yr -- away
    # from what the sediment budget can afford, not toward it.
    #
    # joint_err is still computed and reported, so the change is visible in
    # the output rather than only in this comment.
    err_columns = [f"err_{p}" for p in PERIODS]
    surface["joint_err"] = surface[err_columns].sum(axis=1)
    surface["fit_err"] = surface[f"err_{PERIODS[0]}"]
    surface = surface.sort_values("fit_err").reset_index()

    best = surface.iloc[0]
    fitted_M, fitted_f = float(best["M"]), float(best["fraction"])

    # A value sitting on the edge of its grid is a bound, not an optimum: the
    # search wanted to keep going and ran out of grid. Reported explicitly so
    # a railed result cannot be quoted as a fitted parameter.
    bounds = []
    if fitted_M in (min(M_VALUES), max(M_VALUES)):
        bounds.append("M")
    if fitted_f in (min(F_VALUES), max(F_VALUES)):
        bounds.append("f")

    fit = dict(
        preset=preset,
        M=fitted_M,
        fraction=fitted_f,
        fit_err=float(best["fit_err"]),
        fit_period=PERIODS[0],
        joint_err=float(best["joint_err"]),
        at_grid_bound=bounds,
        # be1 is period 1's; period 2 has none to fit.
        be1_1984=(None if pd.isna(best.get("be1_1984"))
                  else float(best["be1_1984"])),
        **{f"err_{p}": float(best[f"err_{p}"]) for p in PERIODS},
        **{f"differential_{p}": float(best[f"differential_{p}"])
           for p in PERIODS},
        **{f"observed_differential_{p}": float(OBSERVED_DIFFERENTIAL[p])
           for p in PERIODS},
        periods_reachable={str(p): bool(PERIOD_DIFFERENTIAL_IS_REACHABLE[p])
                           for p in PERIODS},
    )
    return surface, fit


def plot_surface(surface, fit, preset):
    """Draws the joint score over the M-f grid, with the ridge visible."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grid = surface.pivot(index="fraction", columns="M", values="fit_err")
    figure, axis = plt.subplots(figsize=(9, 5.5))
    mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                           shading="nearest", cmap="viridis_r")
    figure.colorbar(mesh, ax=axis,
                    label="joint |differential - observed|, both periods (m/yr)")
    axis.plot(fit["M"], fit["fraction"], marker="*", markersize=20,
              color=GROIN_COLOR, markeredgecolor="white", markeredgewidth=1.2,
              linestyle="none", label="best cell", zorder=5)
    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    title = f"Joint two-period fit -- {preset}"
    if fit["at_grid_bound"]:
        title += f"   (railed on {', '.join(fit['at_grid_bound'])})"
    axis.set_title(title)
    axis.legend(loc="upper right")
    figure.tight_layout()
    path = FIGURE_DIR / f"joint_{preset}_surface.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path


def plot_constraints(surface, fit, preset):
    """Draws each period's own best-fit valley in (M, f), and where they cross.

    This is the figure that shows WHY the joint fit lands where it does: each
    period contributes a valley, period 1's running along M*(16 + 4f) = const
    and period 2's along M*f = const, and the fitted point is their crossing.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(9, 5.5))
    colors = {PERIODS[0]: MODEL_COLOR, PERIODS[1]: GROIN_COLOR}
    for period in PERIODS:
        grid = surface.pivot(index="fraction", columns="M",
                             values=f"err_{period}")
        # The valley floor: for each f, the M that period likes best.
        valley_M = [grid.columns[int(np.nanargmin(grid.loc[f].values))]
                    for f in grid.index]
        label = f"{period}-{END_YEAR[period]} best M"
        if not PERIOD_DIFFERENTIAL_IS_REACHABLE[period]:
            label += " (bound: target unreachable)"
        axis.plot(valley_M, grid.index, marker="o", color=colors[period],
                  label=label)

    axis.plot(fit["M"], fit["fraction"], marker="*", markersize=20,
              color=OBSERVED_COLOR, linestyle="none",
              label="joint fit", zorder=5)
    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    axis.set_title(f"Where the two periods constrain (M, f) -- {preset}")
    axis.legend(loc="best", fontsize=9)
    axis.grid(alpha=0.3)
    figure.tight_layout()
    path = FIGURE_DIR / f"joint_{preset}_constraints.png"
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--preset", choices=PRESETS, action="append",
                        help="restrict to one preset (repeatable)")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()

    presets = args.preset or list(PRESETS)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fits, surfaces, problems = {}, [], {}
    for preset in presets:
        surface, outcome = joint_fit(preset)
        if surface is None:
            problems[preset] = outcome
            print(f"\n{preset}: SKIPPED -- {outcome}")
            continue
        fits[preset] = outcome
        surface.insert(0, "preset", preset)
        surfaces.append(surface)

        print("\n" + "=" * 72)
        print(f"JOINT FIT  {preset}")
        print("=" * 72)
        print(f"  M = {outcome['M']:g} m/yr,  f = {outcome['fraction']:g}")
        if outcome["be1_1984"] is not None:
            print(f"  be1 (1984, profiled) = {outcome['be1_1984']:+g} m/yr")
        print(f"  joint error {outcome['joint_err']:.3f} m/yr")
        for period in PERIODS:
            reach = ("" if PERIOD_DIFFERENTIAL_IS_REACHABLE[period]
                     else "   [target unreachable -- this leg is a bound]")
            print(f"    {period}-{END_YEAR[period]}: modelled "
                  f"{outcome[f'differential_{period}']:+.3f} vs observed "
                  f"{outcome[f'observed_differential_{period}']:+.3f}, "
                  f"err {outcome[f'err_{period}']:.3f}{reach}")
        if outcome["at_grid_bound"]:
            print(f"\n  WARNING: {', '.join(outcome['at_grid_bound'])} landed "
                  f"on a grid bound.\n"
                  f"  This is a BOUND, not a fitted optimum. Report it as "
                  f"such, and widen\n"
                  f"  the grid in HAT_groin_sweep_config.py if an interior "
                  f"solution is wanted.")

        if not args.no_figures:
            print(f"  figures: {plot_surface(surface, outcome, preset).name}, "
                  f"{plot_constraints(surface, outcome, preset).name}")

    if surfaces:
        pd.concat(surfaces, ignore_index=True).to_csv(JOINT_CSV, index=False)
        print(f"\n  surface written to {JOINT_CSV}")
    if fits:
        JOINT_JSON.write_text(json.dumps(fits, indent=2))
        print(f"  fitted values written to {JOINT_JSON}")

    return 0 if fits else 1


if __name__ == "__main__":
    sys.exit(main())
