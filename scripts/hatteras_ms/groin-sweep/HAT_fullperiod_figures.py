#!/usr/bin/env python3
"""The four sweep outputs for the continuous 1984-2024 groin calibration.

Deliberately the same set the 1967 rig produced, because that set answered the
question the scalar-fillet figures could not: a heatmap showing whether the
optimum is INTERIOR, and profile plots showing whether the winning cell has the
right SHAPE and not merely the right magnitude at one point.

    heatmap.png            profile RMSE over the M-f grid, best cell marked,
                           cells the model refused drawn as gaps rather than
                           silently dropped
    best_fit_profile.png   the winning cell against the observed change profile
    top_n_profiles.png     the best N cells together, so the spread shows how
                           sharply the metric discriminates

SIGN CONVENTION, AND WHY IT DIFFERS FROM THE RIG FIGURE
    These plot SEAWARD-POSITIVE change, because that is what the CoastSat
    chainage target is measured in and converting for display would put two
    conventions in one workflow. The 1967 rig's figures are landward-positive
    ("+ = landward"), so a curve that rises here falls there. The axis label
    states it on every panel rather than relying on the reader to remember.

WHAT AN INTERIOR OPTIMUM WOULD MEAN
    Every earlier attempt at this calibration railed: the best cell sat on a
    grid edge, which means the search wanted to keep going and ran out of grid,
    not that it found a minimum. A best cell with neighbours on all four sides
    is the evidence that this window and this metric can actually identify the
    pair. The heatmap flags the outcome either way.

Usage:
    python HAT_fullperiod_figures.py [--top-n 5]

Reads  output/groin_sweep/fullperiod_1984_2024/results.csv
Writes output/groin_sweep/fullperiod_1984_2024/figures/

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from HAT_fullperiod_target import (  # noqa: E402
    END_YEAR,
    FIT_DOMAINS_GIS,
    START_YEAR,
    observed_change_profile,
)

OUT_ROOT = PROJECT_BASE_DIR / "output" / "groin_sweep" / "fullperiod_1984_2024"
RESULTS_CSV = OUT_ROOT / "results.csv"
FIGURE_DIR = OUT_ROOT / "figures"

GROIN_UPDRIFT_GIS, GROIN_DOWNDRIFT_GIS = 6, 5
OBSERVED_COLOR = "#1A1A1A"
GROIN_COLOR = "#B71C1C"


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def load_results():
    """Scored cells, ranked. Raises if the sweep has not collated yet."""
    if not RESULTS_CSV.exists():
        raise FileNotFoundError(
            f"{RESULTS_CSV} not found -- run HAT_fullperiod_sweep.py first "
            f"(or with --collate-only if the cells are already on disk).")
    frame = pd.read_csv(RESULTS_CSV)
    if frame.empty:
        raise ValueError(f"{RESULTS_CSV.name} has no scored cells.")
    return frame.sort_values("rmse_m").reset_index(drop=True)


def _profile_of(row):
    """The per-domain change profile carried on one results row."""
    return np.array([float(row[f"change_D{d}"]) for d in FIT_DOMAINS_GIS])


def _mark_groin(axis):
    """Groin line plus updrift / downdrift shading, on a domain axis."""
    axis.axvspan(GROIN_DOWNDRIFT_GIS - 0.5, GROIN_DOWNDRIFT_GIS + 0.5,
                 color="#1565C0", alpha=0.10, zorder=0)
    axis.axvspan(GROIN_UPDRIFT_GIS - 0.5, GROIN_UPDRIFT_GIS + 0.5,
                 color=GROIN_COLOR, alpha=0.10, zorder=0)
    axis.axvline((GROIN_DOWNDRIFT_GIS + GROIN_UPDRIFT_GIS) / 2.0,
                 color=GROIN_COLOR, linestyle="--", linewidth=1.8, zorder=2)
    axis.text((GROIN_DOWNDRIFT_GIS + GROIN_UPDRIFT_GIS) / 2.0 + 0.08,
              0.97, "Buxton groin", rotation=90, color=GROIN_COLOR,
              fontsize=8, va="top", transform=axis.get_xaxis_transform())


def _profile_axis(axis):
    axis.axhline(0.0, color="#BBBBBB", linewidth=0.8, zorder=1)
    axis.set_xlabel("GIS domain (south to north)")
    axis.set_ylabel(f"shoreline change {START_YEAR}-{END_YEAR} (m)  "
                    f"[+ = seaward]")
    axis.set_xticks(list(FIT_DOMAINS_GIS))
    axis.grid(alpha=0.25)


def fig_heatmap(frame):
    """Profile RMSE over the M-f grid, with the optimum's interiority stated."""
    plt = _plt()
    groin = frame[frame["M"] > 0]
    grid = groin.pivot_table(index="fraction", columns="M", values="rmse_m")

    figure, axis = plt.subplots(figsize=(10, 6))
    mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                           shading="nearest", cmap="viridis_r")
    figure.colorbar(mesh, ax=axis, label="profile RMSE, D1-D12 (m)")

    best = groin.iloc[0]
    axis.plot(best["M"], best["fraction"], marker="*", markersize=22,
              color=GROIN_COLOR, markeredgecolor="white", markeredgewidth=1.4,
              linestyle="none", zorder=5,
              label=f"best: M={best['M']:g}, f={best['fraction']:.2f} "
                    f"(RMSE {best['rmse_m']:.1f} m)")

    # Interiority is the point of this figure, so it is computed, not eyeballed.
    rails = []
    if best["M"] in (grid.columns.min(), grid.columns.max()):
        rails.append("M")
    if best["fraction"] in (grid.index.min(), grid.index.max()):
        rails.append("f")
    verdict = ("INTERIOR on both axes -- the grid contains the optimum"
               if not rails else
               f"RAILED on {', '.join(rails)} -- the search ran out of grid, "
               f"so this is a bound, not an optimum")

    baseline = frame[frame["M"] == 0]
    sub = (f"no-groin baseline RMSE {baseline.iloc[0]['rmse_m']:.1f} m"
           if not baseline.empty else "no M = 0 baseline")
    # The caveat belongs in the TITLE, not a footnote. This window nets period
    # 1's fillet build against period 2's collapse, so a module whose trapping
    # is bounded at >= 0 can never win here regardless of how it is scored. The
    # sweep bounds M from above; it does not test whether a groin operated.
    # Fitting is done on period 1 -- see why_M60_f06.png.
    axis.set_title(f"Continuous {START_YEAR}-{END_YEAR} groin sweep -- "
                   f"THIS WINDOW CANNOT FIT A GROIN BY CONSTRUCTION\n"
                   f"(it nets period 1's build against period 2's collapse; "
                   f"it bounds M, it does not test one)\n"
                   f"{verdict}   |   {sub}", fontsize=10)
    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    axis.legend(loc="upper right", fontsize=8)
    figure.tight_layout()
    path = FIGURE_DIR / "heatmap.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path, verdict


def fig_best_fit(frame, observed):
    """The winning cell against the observed change profile."""
    plt = _plt()
    gis = np.array(FIT_DOMAINS_GIS, dtype=float)
    obs = np.array([observed[d] for d in FIT_DOMAINS_GIS])
    best = frame[frame["M"] > 0].iloc[0]

    figure, axis = plt.subplots(figsize=(11, 6))
    _mark_groin(axis)
    axis.plot(gis, obs, marker="s", linestyle="--", color=OBSERVED_COLOR,
              linewidth=2.4, label=f"Observed {END_YEAR}", zorder=5)
    axis.plot(gis, _profile_of(best), marker="o", color="#4A148C",
              linewidth=2.2, zorder=4,
              label=f"Best fit (M={best['M']:g}, f={best['fraction']:.2f}, "
                    f"RMSE={best['rmse_m']:.1f} m)")
    baseline = frame[frame["M"] == 0]
    if not baseline.empty:
        axis.plot(gis, _profile_of(baseline.iloc[0]), color="#888888",
                  linestyle=":", linewidth=1.8, zorder=3,
                  label=f"No groin (RMSE={baseline.iloc[0]['rmse_m']:.1f} m)")
    _profile_axis(axis)
    axis.set_title(f"Best-fit sweep result vs observed shoreline change, "
                   f"{START_YEAR}-{END_YEAR}", fontsize=12)
    axis.legend(loc="best", fontsize=9)
    figure.tight_layout()
    path = FIGURE_DIR / "best_fit_profile.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


def fig_top_n(frame, observed, n):
    """The best N cells together; a tight bundle means weak discrimination."""
    plt = _plt()
    gis = np.array(FIT_DOMAINS_GIS, dtype=float)
    obs = np.array([observed[d] for d in FIT_DOMAINS_GIS])
    top = frame[frame["M"] > 0].head(n)

    figure, axis = plt.subplots(figsize=(11, 6))
    _mark_groin(axis)
    axis.plot(gis, obs, marker="s", linestyle="--", color=OBSERVED_COLOR,
              linewidth=2.6, label=f"Observed {END_YEAR}", zorder=10)
    colours = plt.get_cmap("viridis")(np.linspace(0.15, 0.85, len(top)))
    for colour, (_, row) in zip(colours, top.iterrows()):
        axis.plot(gis, _profile_of(row), color=colour, linewidth=1.8, zorder=4,
                  label=f"M={row['M']:g}, f={row['fraction']:.2f} "
                        f"(RMSE={row['rmse_m']:.1f})")
    _profile_axis(axis)
    spread = float(top["rmse_m"].max() - top["rmse_m"].min())
    axis.set_title(f"Top {len(top)} sweep results vs observed shoreline change"
                   f"\nRMSE spread across these cells: {spread:.1f} m",
                   fontsize=12)
    axis.legend(loc="best", fontsize=8)
    figure.tight_layout()
    path = FIGURE_DIR / "top_n_profiles.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--top-n", type=int, default=5)
    args = parser.parse_args()

    frame = load_results()
    observed = observed_change_profile()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    heatmap, verdict = fig_heatmap(frame)
    print(f"  {len(frame)} cells scored")
    print(f"  {verdict}")
    for path in (heatmap, fig_best_fit(frame, observed),
                 fig_top_n(frame, observed, args.top_n)):
        print(f"  wrote {path.name}")
    print(f"\n  -> {FIGURE_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
