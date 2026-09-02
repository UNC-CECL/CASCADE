#!/usr/bin/env python3
"""The calibrated hindcast: both periods, one page.

THE FIGURE
    The headline result, and only the headline result. One row per period,
    each showing the modelled shoreline change rate against the CoastSat
    target across all 90 real domains, with the misfit shaded between them.

    Everything on it is one configuration -- calibBE source/sink, full
    management, groin on at the fitted (M, f) -- because this figure answers
    "how well does the calibrated model reproduce the observed shoreline?"
    and nothing else. `HAT_scenario_grid.py` is the figure for comparing
    presets and scenarios; overlaying those here would turn a result into a
    contrast and bury it.

WHY THE RATE, NOT THE POSITION
    The rate is the calibrated quantity. Every skill number in the project --
    the source/sink residual, the groin fit, the convergence history -- is an
    RMSE on `lrr_m_yr`, an OLS slope through the run's annual states, scored
    against a target built the same way. Plotting position instead would put a
    quantity nobody calibrated on the headline figure and invite a reader to
    judge the model on it.

WHAT IS MARKED, AND WHY EACH ONE MATTERS
    D1 and D90    Locked. Solved separately by buffer-cell reproduction, not
                  fit from a residual, and they carry rates an order of
                  magnitude larger than the interior because only ~10% of an
                  imposed edge rate survives diffusion. They are boundary
                  absorbers, not sediment budgets, and reading them as skill
                  would be a mistake.
    D5-D7         Reserved for the groin. The residual there is the groin
                  module's shortfall and is deliberately never corrected, so
                  the largest misfit on this figure is expected and explained.
                  See fig_groin_reserved_residual.png.
    grey bands    Outside the frozen zone set: diagnosed but never corrected,
                  because zone membership was fixed before the calibration
                  iterated. The misfit there is honest unexplained variance.

    The point of marking all three is that a reader can see WHERE the model is
    allowed to miss and why, rather than reading every excursion as failure.

Usage:
    python HAT_hindcast_final_figure.py [--preset calibBE]

Writes output/comparisons/hindcast_calibrated_full_management.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import importlib.util
import pathlib
import sys

import numpy as np
import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
OUT_DIR = PROJECT_BASE_DIR / "output" / "comparisons"

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

from cascade_pipeline.run_registry import find_run_dir          # noqa: E402

LOESS_PATH = (PROJECT_BASE_DIR / "scripts" / "input_prep" / "7-source-sink"
              / "loess_smooth" / "HAT_be_zone_LOESS_analysis.py")

PERIODS = {
    "1984_2004": dict(label="1984–2004", start=1984, scenario="road_bdm",
                      colour="#1565C0"),
    "2004_2024": dict(label="2004–2024", start=2004,
                      scenario="road_bdm_nourish", colour="#B71C1C"),
}
LOCKED = (1, 90)
RESERVED_COLOUR = "#FF8C00"


def analysis_module():
    spec = importlib.util.spec_from_file_location("_loess", LOESS_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_rates(period, preset, scenario):
    name = f"HAT_{period}_{preset}_{scenario}_groin"
    # RESOLVED, NOT JOINED BY HAND. A run forced off the calibration wave
    # climate is filed under an arm component that a <period>/<preset>/<name>
    # join has no slot for, so the run would read as missing. find_run_dir
    # defaults to the calibration arm -- this figure's runs -- and names the
    # arms a run IS under when it is not there.
    try:
        run_dir = find_run_dir(RAW_RUNS, name, period, preset)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"{exc}"
            f"  Run HAT_run_all.py --stages 2,6 --presets {preset} first."
        ) from None
    path = run_dir / f"{name}_shoreline_change_rate.csv"
    return pd.read_csv(path).set_index("gis_domain")["lrr_m_yr"], name


def skill(model, target, domains):
    """RMSE, bias and correlation over the scored interior."""
    shared = [d for d in domains if d in model.index and d in target.index]
    m = np.array([model[d] for d in shared], dtype=float)
    t = np.array([target[d] for d in shared], dtype=float)
    good = ~(np.isnan(m) | np.isnan(t))
    m, t = m[good], t[good]
    return (float(np.sqrt(np.mean((m - t) ** 2))), float(np.mean(m - t)),
            float(np.corrcoef(m, t)[0, 1]), len(m))


def contiguous(domains):
    """[(first, last), ...] runs, for shading a domain set compactly."""
    out, run = [], []
    for d in sorted(domains):
        if run and d == run[-1] + 1:
            run.append(d)
        else:
            if run:
                out.append((run[0], run[-1]))
            run = [d]
    if run:
        out.append((run[0], run[-1]))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--preset", default="calibBE")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    module = analysis_module()
    reserved = list(module.GROIN_RESERVED_DOMAINS)
    frozen = module.FROZEN_ZONE_DOMAINS

    figure, axes = plt.subplots(2, 1, figsize=(16, 9.4), sharex=True)

    summary = []
    for axis, (key, meta) in zip(axes, PERIODS.items()):
        csv = (module.P1_COASTSAT_CSV if meta["start"] == 1984
               else module.P2_COASTSAT_CSV)
        target = module.load_observed(meta["start"], csv)[1]
        model, run_name = run_rates(key, args.preset, meta["scenario"])

        gis = np.arange(1, 91)
        t = np.array([target.get(g, np.nan) for g in gis], dtype=float)
        m = np.array([model.get(g, np.nan) for g in gis], dtype=float)

        # Withheld domains first, so every line and marker sits on top of them.
        for a, b in contiguous(set(range(2, 90)) - set(frozen[meta["start"]])):
            axis.axvspan(a - 0.5, b + 0.5, color="#C8C8C8", alpha=0.30, zorder=0)
        for a, b in contiguous(reserved):
            axis.axvspan(a - 0.5, b + 0.5, color=RESERVED_COLOUR, alpha=0.28,
                         zorder=1)
        for locked in LOCKED:
            axis.axvspan(locked - 0.5, locked + 0.5, color="#5E35B1", alpha=0.22,
                         zorder=1)

        axis.fill_between(gis, t, m, where=~(np.isnan(t) | np.isnan(m)),
                          color=meta["colour"], alpha=0.16, zorder=2,
                          label="misfit")
        axis.plot(gis, t, color="#1A1A1A", linewidth=2.6, linestyle="--",
                  marker="s", markersize=4.5, zorder=6, label="observed (CoastSat)")
        axis.plot(gis, m, color=meta["colour"], linewidth=2.4, marker="o",
                  markersize=4.5, zorder=5, label="CASCADE, calibrated")
        axis.axhline(0.0, color="#999999", linewidth=0.9, zorder=3)

        # Scored on the INTERIOR: D1 and D90 are boundary absorbers and their
        # rates are ~10x the interior, so including them would let two tuned
        # buffer cells dominate a 90-domain skill number.
        rmse, bias, corr, n = skill(model, target, range(2, 90))
        summary.append((meta["label"], run_name, rmse, bias, corr, n))
        axis.set_title(
            f"{meta['label']}    RMSE {rmse:.3f} m/yr   bias {bias:+.3f}   "
            f"r = {corr:.2f}   (D2–D89, n = {n})",
            fontsize=12, loc="left")
        axis.set_ylabel("shoreline change rate (m/yr)\n[+ = seaward]", fontsize=10)
        axis.grid(alpha=0.22)

    axes[1].set_xlabel("GIS domain  (south → north:  D1 = Cape Point, "
                       "D90 = Pea Island)", fontsize=11)
    axes[1].set_xlim(0, 91)
    axes[1].set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])

    handles, labels = axes[0].get_legend_handles_labels()
    handles += [
        Patch(facecolor=RESERVED_COLOUR, alpha=0.28,
              label="D5–D7 reserved for the groin"),
        Patch(facecolor="#C8C8C8", alpha=0.30,
              label="outside the frozen zone set — never corrected"),
        Patch(facecolor="#5E35B1", alpha=0.22,
              label="D1 / D90 locked (boundary absorbers)")]
    axes[0].legend(handles=handles, loc="upper left", fontsize=8.5, ncol=2)

    figure.suptitle(
        "Cape Hatteras hindcast, calibrated — full management, groin on, "
        "calibrated source/sink", fontsize=14, y=0.985)
    figure.tight_layout(rect=(0, 0.105, 1, 0.965))
    figure.text(
        0.01, 0.012,
        "ONE CONFIGURATION ONLY: calibBE source/sink, full management (roadway + beach/dune, with nourishment in 2004-2024), groin on at the "
        "fitted M = 60 m/yr, f = 0.6. This figure answers how well the calibrated model reproduces the observed shoreline; scenario_grid_by_preset.png "
        "is where presets and management scenarios are compared. The quantity is lrr_m_yr, an OLS slope through the run's annual states, scored "
        "against a target built the same way -- the rate is what was calibrated, so it is what the model is judged on.\n"
        "WHERE THE MODEL IS ALLOWED TO MISS. Orange: D5-D7, reserved for the groin -- that residual is the groin module's shortfall (too little "
        "fillet built in period 1, no release in period 2) and absorbing it into the source/sink field would double-count against the M/f fit. "
        "Grey: outside the frozen zone set, diagnosed but never corrected, because zone membership was fixed before the calibration iterated rather "
        "than grown to fit. Purple: D1 and D90, solved separately by buffer-cell reproduction and carrying rates ~10x the interior, so they are "
        "excluded from the skill numbers rather than allowed to dominate them.",
        fontsize=7.4, color="#333333", wrap=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "hindcast_calibrated_full_management.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)

    print(f"wrote {path}")
    for label, run_name, rmse, bias, corr, n in summary:
        print(f"  {label}  {run_name}")
        print(f"      RMSE {rmse:.4f} m/yr   bias {bias:+.4f}   r {corr:.3f}   n {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
