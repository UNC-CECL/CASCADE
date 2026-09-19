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
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, error_cmap, figsize, open_frame,
                              save)
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
# House colours (2026-09-11): observations in INK; the structure is a
# muted guide line, not the vintage red it used to be drawn in.
OBSERVED_COLOR = INK
GROIN_COLOR = INK_MUTED


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
    # Two greys, not two hues: these bands say WHERE the structure is, and
    # the colour on this figure belongs to what is plotted.
    axis.axvspan(GROIN_DOWNDRIFT_GIS - 0.5, GROIN_DOWNDRIFT_GIS + 0.5,
                 color="0.94", zorder=0)
    axis.axvspan(GROIN_UPDRIFT_GIS - 0.5, GROIN_UPDRIFT_GIS + 0.5,
                 color="0.90", zorder=0)
    axis.axvline((GROIN_DOWNDRIFT_GIS + GROIN_UPDRIFT_GIS) / 2.0,
                 color=INK_MUTED, linestyle=(0, (4, 2)), linewidth=0.8,
                 zorder=2)
    axis.text((GROIN_DOWNDRIFT_GIS + GROIN_UPDRIFT_GIS) / 2.0 + 0.08,
              0.97, "Buxton groin", rotation=90, color=INK_MUTED,
              fontsize=7, va="top", transform=axis.get_xaxis_transform())


def _profile_axis(axis):
    axis.axhline(0.0, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 3)),
                 zorder=1)
    axis.set_xlabel("GIS domain (south → north)")
    axis.set_ylabel(f"shoreline change {START_YEAR} to {END_YEAR} (m)\n"
                    f"positive is seaward")
    axis.set_xticks(list(FIT_DOMAINS_GIS))
    axis.grid(axis="y")
    axis.set_axisbelow(True)
    open_frame(axis)


def fig_heatmap(frame):
    """Profile RMSE over the M-f grid, with the optimum's interiority stated."""
    plt = _plt()
    groin = frame[frame["M"] > 0]
    grid = groin.pivot_table(index="fraction", columns="M", values="rmse_m")

    apply_style()
    figure, axis = plt.subplots(figsize=figsize("double", aspect=0.48),
                                constrained_layout=True)
    mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                           shading="nearest", cmap=error_cmap())
    cb = figure.colorbar(mesh, ax=axis)
    cb.set_label("profile RMSE, D1−D12 (m)")
    cb.outline.set_linewidth(0.6)

    best = groin.iloc[0]
    axis.plot(best["M"], best["fraction"], marker="*", markersize=13,
              color=C["ACCENT"], markeredgecolor="white", markeredgewidth=0.8,
              linestyle="none", zorder=5,
              label=f"best: M {best['M']:g}, f {best['fraction']:.2f}, "
                    f"{best['rmse_m']:.1f} m")

    # Interiority is the point of this figure, so it is computed, not eyeballed.
    rails = []
    if best["M"] in (grid.columns.min(), grid.columns.max()):
        rails.append("M")
    if best["fraction"] in (grid.index.min(), grid.index.max()):
        rails.append("f")
    verdict = ("INTERIOR on both axes, so the grid contains the optimum"
               if not rails else
               f"RAILED on {', '.join(rails)} — the search ran out of grid, "
               f"so this is a bound and not an optimum")

    baseline = frame[frame["M"] == 0]
    sub = (f"no-groin baseline RMSE {baseline.iloc[0]['rmse_m']:.1f} m"
           if not baseline.empty else "no M = 0 baseline")
    # The caveat belongs in the TITLE, not a footnote. This window nets period
    # 1's fillet build against period 2's collapse, so a module whose trapping
    # is bounded at >= 0 can never win here regardless of how it is scored. The
    # sweep bounds M from above; it does not test whether a groin operated.
    # Fitting is done on period 1 -- see why_M60_f06.png.
    axis.set_title(f"Profile error over the (M, f) grid, "
                   f"{START_YEAR} to {END_YEAR}", loc="left")
    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    axis.legend(loc="upper right", fontsize=7)

    caption(figure,
            "The continuous {start} to {end} sweep's profile error over the "
            "(M, f) grid, dark worse, with the best cell marked. THIS WINDOW "
            "CANNOT FIT A GROIN BY CONSTRUCTION: it nets period 1's fillet "
            "build against period 2's collapse, so a module whose trapping is "
            "bounded at or above zero can never win here however it is "
            "scored. The sweep bounds M from above; it does not test whether a "
            "groin operated, and the fitting is done on period 1 alone. On "
            "this grid the optimum is {verdict}, against a {sub}."
            .format(start=START_YEAR, end=END_YEAR,
                    verdict=verdict[0].lower() + verdict[1:], sub=sub))

    return save(figure, FIGURE_DIR / "heatmap.png", close=True)[0], verdict


def fig_best_fit(frame, observed):
    """The winning cell against the observed change profile."""
    plt = _plt()
    gis = np.array(FIT_DOMAINS_GIS, dtype=float)
    obs = np.array([observed[d] for d in FIT_DOMAINS_GIS])
    best = frame[frame["M"] > 0].iloc[0]

    apply_style()
    figure, axis = plt.subplots(figsize=figsize("double", aspect=0.46),
                                constrained_layout=True)
    _mark_groin(axis)
    axis.plot(gis, obs, marker="s", markersize=3.4, linestyle="--",
              color=OBSERVED_COLOR, linewidth=1.8,
              label=f"observed at {END_YEAR}", zorder=5)
    axis.plot(gis, _profile_of(best), marker="o", markersize=3.4,
              color=C["ACCENT"], linewidth=1.6, zorder=4,
              label=f"best fit, M {best['M']:g}, f {best['fraction']:.2f}, "
                    f"{best['rmse_m']:.1f} m")
    baseline = frame[frame["M"] == 0]
    if not baseline.empty:
        axis.plot(gis, _profile_of(baseline.iloc[0]), color=C["BASE"],
                  linestyle=":", linewidth=1.4, zorder=3,
                  label=f"no groin, {baseline.iloc[0]['rmse_m']:.1f} m")
    _profile_axis(axis)
    axis.set_title(f"The best cell against the observed change, "
                   f"{START_YEAR} to {END_YEAR}", loc="left")
    axis.legend(loc="best", fontsize=7)

    caption(figure,
            "The best-scoring cell of the continuous {start} to {end} sweep "
            "against the observed change profile, with the no-groin baseline "
            "for reference. Remember what this window is: it nets period 1's "
            "fillet build against period 2's collapse, so it bounds M rather "
            "than testing whether a groin operated. The pair taken forward is "
            "fitted on period 1 alone."
            .format(start=START_YEAR, end=END_YEAR))

    return save(figure, FIGURE_DIR / "best_fit_profile.png", close=True)[0]


def fig_top_n(frame, observed, n):
    """The best N cells together; a tight bundle means weak discrimination."""
    plt = _plt()
    gis = np.array(FIT_DOMAINS_GIS, dtype=float)
    obs = np.array([observed[d] for d in FIT_DOMAINS_GIS])
    top = frame[frame["M"] > 0].head(n)

    apply_style()
    from matplotlib.colors import LinearSegmentedColormap
    figure, axis = plt.subplots(figsize=figsize("double", aspect=0.46),
                                constrained_layout=True)
    _mark_groin(axis)
    axis.plot(gis, obs, marker="s", markersize=3.4, linestyle="--",
              color=OBSERVED_COLOR, linewidth=1.8,
              label=f"observed at {END_YEAR}", zorder=10)
    # One colour family, darkest is best.
    ramp = LinearSegmentedColormap.from_list(
        "hat_accent_ramp", [C["ACCENT_FILL"], C["ACCENT"]])
    for colour, (_, row) in zip(ramp(np.linspace(1.0, 0.2, len(top))),
                                top.iterrows()):
        axis.plot(gis, _profile_of(row), color=colour, linewidth=1.4, zorder=4,
                  label=f"M {row['M']:g}, f {row['fraction']:.2f}, "
                        f"{row['rmse_m']:.1f} m")
    _profile_axis(axis)
    spread = float(top["rmse_m"].max() - top["rmse_m"].min())
    axis.set_title(f"The top {len(top)} cells against the observed change",
                   loc="left")
    axis.legend(loc="best", fontsize=7)

    caption(figure,
            "The {n} best-scoring cells of the continuous {start} to {end} "
            "sweep on one set of axes, best darkest, against the observed "
            "change profile. Their errors span {spread:.1f} m: a tight bundle "
            "means this window discriminates weakly between cells, which is a "
            "statement about identifiability rather than about fit. The window "
            "nets period 1's build against period 2's collapse, so it bounds M "
            "rather than testing a groin."
            .format(n=len(top), start=START_YEAR, end=END_YEAR,
                    spread=spread))

    return save(figure, FIGURE_DIR / "top_n_profiles.png", close=True)[0]


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
