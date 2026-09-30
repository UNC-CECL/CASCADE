#!/usr/bin/env python3
"""
Fit M and f together from the two periods' sweep surfaces.

    python scripts/hatteras_ms/groin-sweep/HAT_groin_joint_fit.py --preset edgeBE
    python scripts/hatteras_ms/groin-sweep/HAT_groin_joint_fit.py --no-figures

Ranked on period 1 (the only window that separates M and f); writes
joint_fit.json and the surface figures. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# parents[3]: this file is in hatteras_ms/groin-sweep/; the guard makes a move fail loudly here
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/hatteras_ms/groin-sweep/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997,  # noqa: E402
                              caption, error_cmap, figsize, open_frame, save)
from HAT_groin_sweep_config import (  # noqa: E402
    GROIN_SWEEP_ROOT,
    END_YEAR,
    F_VALUES,
    M_VALUES,
    OBSERVED_DIFFERENTIAL,
    PERIODS,
    PERIOD_DIFFERENTIAL_IS_REACHABLE,
    PRESETS,
    sweep_output_dir,
    joint_fit_paths,
)

# --- CONFIG ------------------------------------------------------------------
OUTPUT_DIR = GROIN_SWEEP_ROOT
# Figures live in a subdirectory
FIGURE_DIR = OUTPUT_DIR / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
JOINT_JSON, JOINT_CSV = joint_fit_paths()

# House colours (2026-09-11): observations in INK, the run under test the ACCENT, guides muted
EARLY_COLOR, LATE_COLOR = C_1984, C_1997
FIT_COLOR = C["ACCENT"]
# -----------------------------------------------------------------------------


# Loads one sweep's scored results
def load_period(period, preset):
    path = sweep_output_dir(period, preset) / "sweep_results.csv"
    if not path.exists():
        return None
    frame = pd.read_csv(path)
    return frame[frame["differential_err"].notna()].copy()


# Reduces one period's rows to a value per (M, f) cell
def period_surface(frame, period):
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


# Builds the joint surface for one preset and picks its best cell
def joint_fit(preset):
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

    # Ranked on period 1 alone: period 2 sees only M*f; joint_err still reported
    err_columns = [f"err_{p}" for p in PERIODS]
    surface["joint_err"] = surface[err_columns].sum(axis=1)
    surface["fit_err"] = surface[f"err_{PERIODS[0]}"]
    surface = surface.sort_values("fit_err").reset_index()

    best = surface.iloc[0]
    fitted_M, fitted_f = float(best["M"]), float(best["fraction"])

    # A value sitting on the edge of its grid is a bound, not an optimum
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


# Draws the joint score over the M-f grid, with the ridge visible
def plot_surface(surface, fit, preset):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_style()
    grid = surface.pivot(index="fraction", columns="M", values="fit_err")
    figure, axis = plt.subplots(figsize=figsize("double", aspect=0.46),
                                constrained_layout=True)
    mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                           shading="nearest", cmap=error_cmap())
    cb = figure.colorbar(mesh, ax=axis)
    cb.set_label("joint |modelled − observed| fillet trend,\n"
                 "both periods (m/yr)")
    cb.outline.set_linewidth(0.6)
    axis.plot(fit["M"], fit["fraction"], marker="*", markersize=14,
              color=FIT_COLOR, markeredgecolor="white", markeredgewidth=0.8,
              linestyle="none", label="best cell", zorder=5)
    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    axis.set_title("Joint two-period fit", loc="left")
    axis.legend(loc="upper right")

    railed = (" The fit is RAILED on {}, so it is a grid bound rather than an "
              "interior minimum.".format(", ".join(fit["at_grid_bound"]))
              if fit["at_grid_bound"] else
              " The fit is an interior minimum, not a grid bound.")
    caption(figure,
            "The joint two-period score over the (M, f) grid for the {p} "
            "source/sink preset: for each cell, how far the modelled fillet "
            "trend sits from the observed one in both hindcast periods at "
            "once. Dark is worse, and the marked cell is the best on this "
            "score.{railed} Read this beside the constraints figure for the "
            "same preset, which separates the two periods' own valleys and "
            "shows where they cross. The joint fit is recorded here because it "
            "was attempted, not because it is the answer: period 2's observed "
            "gap NARROWS, a groin with trapping at or above zero can only "
            "widen it, and so fitting the two periods together asks for "
            "something the parameterisation cannot produce. The production "
            "pair is fitted on period 1 alone."
            .format(p=preset, railed=railed))

    path = FIGURE_DIR / f"joint_{preset}_surface.png"
    return save(figure, path, close=True)[0]


# Draws each period's own best-fit valley in (M, f), and where they cross
def plot_constraints(surface, fit, preset):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_style()
    figure, axis = plt.subplots(figsize=figsize("double", aspect=0.46),
                                constrained_layout=True)
    colors = {PERIODS[0]: EARLY_COLOR, PERIODS[1]: LATE_COLOR}
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

    axis.plot(fit["M"], fit["fraction"], marker="*", markersize=14,
              color=FIT_COLOR, markeredgecolor="white", markeredgewidth=0.8,
              linestyle="none", label="joint fit", zorder=5)
    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    axis.set_title("Where each period constrains the pair", loc="left")
    axis.legend(loc="best", fontsize=7.5)
    axis.grid()
    axis.set_axisbelow(True)
    open_frame(axis)

    unreachable = [str(p) for p in PERIODS
                   if not PERIOD_DIFFERENTIAL_IS_REACHABLE[p]]
    caption(figure,
            "Why the joint fit for the {p} preset lands where it does. Each "
            "line is one period's own valley floor: for every deterioration "
            "floor f, the trapping rate M that period scores best. Period 1 "
            "runs along constant M(16 + 4f) because it mostly precedes the "
            "1996 to 2003 deterioration ramp, period 2 along constant M·f "
            "because it lies entirely after it, and the marked point is where "
            "the two cross. The earlier period is red and the later blue, as "
            "everywhere in this project. {un} A valley drawn against an "
            "unreachable target is a grid bound, not a constraint: the module "
            "can only widen the gap between the structure's flanks, and "
            "period 2's observed gap narrows, which is why the production pair "
            "is fitted on period 1 alone."
            .format(p=preset,
                    un=("The target is UNREACHABLE for {}, so that line is "
                        "railed.".format(" and ".join(unreachable))
                        if unreachable else
                        "Both periods' targets are reachable on this grid."))
            )

    path = FIGURE_DIR / f"joint_{preset}_constraints.png"
    return save(figure, path, close=True)[0]


# Presets in an existing joint_fit.json that were set by hand, not fitted
def _pinned_presets(path):
    if not path.exists():
        return []
    try:
        existing = json.loads(path.read_text())
    except (OSError, ValueError):
        return []
    return sorted(k for k, v in existing.items()
                  if isinstance(v, dict) and v.get("provenance"))


# Writes the ranking's answer, unless it would clobber a hand pin
def _write_fits(fits, force=False):
    pinned = _pinned_presets(JOINT_JSON)
    if pinned and not force:
        sidecar = JOINT_JSON.with_name("joint_fit_ranking.json")
        sidecar.write_text(json.dumps(fits, indent=2))
        print("")
        print(f"  {'=' * 66}")
        print(f"  REFUSED to overwrite {JOINT_JSON.name}: it holds hand-pinned "
              f"values")
        print(f"  pinned presets   {', '.join(pinned)}")
        for preset, fit in sorted(fits.items()):
            bound = fit.get("at_grid_bound") or []
            flag = f"   RAILED on {', '.join(bound)}" if bound else ""
            print(f"  this ranking     {preset:<8} M = {fit.get('M')}, "
                  f"f = {fit.get('fraction')}{flag}")
        print(f"  ranking written to {sidecar.name} instead")
        print(f"  re-run with --force to overwrite the pin deliberately")
        print(f"  {'=' * 66}")
        return
    JOINT_JSON.write_text(json.dumps(fits, indent=2))
    print(f"  fitted values written to {JOINT_JSON}")


# Run: load both periods' sweeps, rank, write the fit and figures
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--preset", choices=PRESETS, action="append",
                        help="restrict to one preset (repeatable)")
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument(
        "--force", action="store_true",
        help="overwrite joint_fit.json even if it holds hand-pinned "
             "values. Without this, a pinned file is preserved and the "
             "ranking is written beside it instead.")
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
        _write_fits(fits, force=args.force)

    return 0 if fits else 1


if __name__ == "__main__":
    sys.exit(main())
