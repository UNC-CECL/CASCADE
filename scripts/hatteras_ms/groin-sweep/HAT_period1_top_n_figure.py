#!/usr/bin/env python3
"""Top-N sweep results against observed change -- PERIOD 1, production geometry.

The direct counterpart to the 1967 rig's `HAT_groin_sweep_top_n_profiles.png`,
so the two calibrations can be read side by side. Same question: do the
best-scoring cells actually reproduce the observed alongshore SHAPE, or do they
only match a summary number?

THREE DIFFERENCES FROM THE RIG FIGURE, ALL DELIBERATE

    window      1984-2004, not 1967-2018. Period 1 is the only window in the
                hindcast where the observed gap between the groin's flanks
                WIDENS, which is the only behaviour a module with trapping
                >= 0 can produce.

    geometry    120-domain production grid, not the rig's 41. M is
                grid-specific -- a confined array preserves dipole amplitude
                that an open one diffuses away -- so a value fitted here
                transfers to the hindcast and one fitted on the rig does not.

    score       DEMEANED, and ranked on D4-D8 only. A uniform level offset in
                the groin's neighbourhood is absorbed by the source/sink
                calibration downstream, so it is not the groin's job; what the
                groin must get right is the shape. D1 is excluded from the
                ranking because the cape's change over period 1 is 81-104 m,
                roughly five times the groin's signal, and it swamps it.

WHAT TO LOOK FOR
    The no-groin baseline is drawn alongside the top cells. If the groin is
    doing real work the top cells should sit closer to the observations than
    that grey line does, in the shaded fit window. They do: 15.58 -> 11.69 m
    for the chosen cell (M = 60, f = 0.6), a 25% reduction.

    Watch also how tightly the top five bundle together. They span M = 40-95
    and f = 0.4-1.0 yet differ by under 0.5 m, which is the visual statement of
    the ridge in period-1 cumulative trapping, M(15.5 + 4.5f). An earlier
    version of this caption said "the metric identifies a product, not a
    pair"; fig_Mf_identifiability.png tested that and refuted it
    (corr(RMSE, M*f) = -0.07). See CALIBRATION_FIGURES.md and GROIN_PLAN.md.

Usage:
    python HAT_period1_top_n_figure.py [--top-n 5]

Writes output/groin_sweep/figures/period1_top_n_profiles.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402

from HAT_fullperiod_target import observed_change_profile  # noqa: E402

SWEEP = PROJECT_BASE_DIR / "output" / "groin_sweep" / "1984_2004_edgeBE"
FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

SHOW_DOMAINS = list(range(1, 13))     # plot the whole neighbourhood
FIT_DOMAINS = list(range(4, 9))       # but rank on D4-D8 only
# The CALIBRATED value, re-solved 2026-08-28. This was -40.0, which was the
# nearest grid point to the old -41.8 rather than the value itself; -42.6
# is now on the grid and is what production spends.
PINNED_BE1 = -42.6                    # production's edgeBE value for 1984
CHOSEN_M, CHOSEN_F = 60.0, 0.6
GROIN_COLOR = "#B71C1C"


def load():
    """Cells at the pinned be1, with their profiles and D4-D8 demeaned score."""
    # LANDWARD-POSITIVE, so erosion is UP and the panel reads as a plan view,
    # matching the gifs. observed_change_profile is SEAWARD-positive, so it is
    # negated; x_s is landward-positive already, so the negation that used to
    # sit on `change` below is gone. Both series flip together, so `score` is
    # unchanged, and the scoring pipeline itself is untouched.
    show_obs = -np.array([observed_change_profile(1984, 2004, SHOW_DOMAINS)[d]
                          for d in SHOW_DOMAINS])
    fit_index = [SHOW_DOMAINS.index(d) for d in FIT_DOMAINS]
    fit_obs = show_obs[fit_index]
    fit_shape = fit_obs - fit_obs.mean()

    rows = []
    for cell in sorted(SWEEP.iterdir()):
        path = cell / "shoreline_matrix.npy"
        if not path.exists():
            continue
        name = cell.name
        be1 = np.nan if "beNA" in name else float(name.split("_be")[1].split("_f")[0])
        if not (np.isnan(be1) or be1 == PINNED_BE1):
            continue
        change = np.load(path)[-1] - np.load(path)[0]
        profile = np.array([change[GEOMETRY.gis_to_pad(d)] for d in SHOW_DOMAINS])
        fit = profile[fit_index]
        rows.append(dict(
            M=0.0 if name.startswith("M0_") else float(name.split("_")[0][1:]),
            f=0.0 if name.startswith("M0_") else float(name.split("_f")[1]),
            score=float(np.sqrt(np.mean((fit - fit.mean() - fit_shape) ** 2))),
            profile=profile))
    return pd.DataFrame(rows), show_obs


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--top-n", type=int, default=5)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame, observed = load()
    groin = frame[frame.M > 0].sort_values("score").reset_index(drop=True)
    baseline = frame[frame.M == 0].sort_values("score").iloc[0]
    top = groin.head(args.top_n)

    # Demeaned over the FIT window, so every curve is centred the same way the
    # score centres it. Centring on the plotted window instead would show a
    # different quantity from the one that was ranked.
    fit_index = [SHOW_DOMAINS.index(d) for d in FIT_DOMAINS]
    centre = lambda v: np.asarray(v, float) - np.mean(np.asarray(v, float)[fit_index])

    x = np.array(SHOW_DOMAINS, dtype=float)
    figure, axis = plt.subplots(figsize=(13, 7))

    axis.axvspan(min(FIT_DOMAINS) - 0.5, max(FIT_DOMAINS) + 0.5,
                 color="#FFF9C4", alpha=0.55, zorder=0)
    axis.text((min(FIT_DOMAINS) + max(FIT_DOMAINS)) / 2, 0.02,
              f"FIT WINDOW D{min(FIT_DOMAINS)}-D{max(FIT_DOMAINS)}",
              ha="center", va="bottom", fontsize=9, color="#8D6E00",
              weight="bold", transform=axis.get_xaxis_transform(), zorder=6)
    axis.axvline(5.5, color=GROIN_COLOR, linestyle="--", linewidth=1.8, zorder=2)
    axis.text(5.58, 0.97, "Buxton groin", rotation=90, fontsize=8.5,
              color=GROIN_COLOR, va="top", transform=axis.get_xaxis_transform())

    axis.plot(x, centre(observed), marker="s", markersize=8, linestyle="--",
              color="#1A1A1A", linewidth=2.8, zorder=10, label="Observed 2004")
    axis.plot(x, centre(baseline.profile), color="#888888", linestyle=":",
              linewidth=2.2, zorder=4,
              label=f"no groin  (RMSE={baseline.score:.1f} m)")

    colours = plt.get_cmap("viridis")(np.linspace(0.12, 0.82, len(top)))
    for colour, (_, row) in zip(colours, top.iterrows()):
        marker = "o" if (row.M == CHOSEN_M and row.f == CHOSEN_F) else None
        axis.plot(x, centre(row.profile), color=colour, linewidth=2.0,
                  marker=marker, markersize=7, zorder=5,
                  label=f"M={row.M:g}, frac={row.f:.2f}  (RMSE={row.score:.1f})")

    chosen = frame[(frame.M == CHOSEN_M) & (frame.f == CHOSEN_F)]
    if not chosen.empty and not ((top.M == CHOSEN_M) & (top.f == CHOSEN_F)).any():
        row = chosen.iloc[0]
        axis.plot(x, centre(row.profile), color=GROIN_COLOR, linewidth=2.6,
                  marker="o", markersize=7, zorder=8,
                  label=f"CHOSEN M={CHOSEN_M:g}, f={CHOSEN_F:g}  "
                        f"(RMSE={row.score:.1f})")

    axis.set_xticks(SHOW_DOMAINS)
    axis.set_xlabel("GIS Domain ID (D1-D12)")
    axis.set_ylabel("Shoreline change 1984->2004 (m), demeaned over the fit window\n"
                    "[+ = landward, erosion]")
    spread = float(top.score.max() - top.score.min())
    axis.set_title(
        f"Top {len(top)} sweep results vs observed shoreline change -- "
        f"PERIOD 1, production geometry\n"
        f"ranked on D4-D8 shape; RMSE spread across these cells: {spread:.2f} m "
        f"(M {top.M.min():g}-{top.M.max():g}, f {top.f.min():g}-{top.f.max():g})",
        fontsize=12)
    axis.grid(alpha=0.25)
    # Lower left: after the 2026-08-30 flip to landward-positive the observed
    # curve occupies the upper left.
    axis.legend(loc="lower left", fontsize=9)

    figure.tight_layout(rect=(0, 0.075, 1, 1))
    figure.text(
        0.01, 0.012,
        "Counterpart to the 1967 rig's top-N figure, for comparison. Differences are deliberate: PERIOD 1 (the only hindcast "
        "window where the observed gap WIDENS, which is all a module with trapping >= 0 can produce), PRODUCTION geometry (M is "
        "grid-specific, so a value fitted here transfers and one fitted on the rig does not), and a DEMEANED score over D4-D8 "
        "(a uniform level offset is absorbed by the source/sink calibration downstream; D1 is excluded because the cape's 81-104 m "
        "change swamps the groin's ~17 m signal). The top cells bundle tightly despite spanning a wide M and f range -- that is "
        "the ridge in period-1 cumulative trapping, M(15.5 + 4.5f) -- NOT in "
        "M*f, which was tested and refuted (corr = -0.07).",
        fontsize=7.5, color="#333333", wrap=True)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURE_DIR / "period1_top_n_profiles.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    print(f"wrote {path}")
    print(f"  no groin {baseline.score:.2f} m")
    print(top[["M", "f", "score"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
