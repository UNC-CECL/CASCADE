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
from hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save)

SWEEP = PROJECT_BASE_DIR / "output" / "groin_sweep" / "1984_2004_edgeBE"
FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

SHOW_DOMAINS = list(range(1, 13))     # plot the whole neighbourhood
FIT_DOMAINS = list(range(4, 9))       # but rank on D4-D8 only
# The CALIBRATED value, re-solved 2026-08-28. This was -40.0, which was the
# nearest grid point to the old -41.8 rather than the value itself; -42.6
# is now on the grid and is what production spends.
PINNED_BE1 = -42.6                    # production's edgeBE value for 1984
CHOSEN_M, CHOSEN_F = 60.0, 0.6
# House semantics: the cell taken forward is the ACCENT, the no-groin baseline
# is BASE, the structure's position is a guide line in muted ink. GROIN_COLOR
# was a dark red here, which is the 1984 vintage colour elsewhere.
CHOSEN_COLOR, BASELINE_COLOR, BAND = C["ACCENT"], C["BASE"], "0.94"


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
    from matplotlib.colors import LinearSegmentedColormap

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
    apply_style()
    figure, axis = plt.subplots(figsize=figsize("double", aspect=0.48),
                                constrained_layout=True)

    axis.axvspan(min(FIT_DOMAINS) - 0.5, max(FIT_DOMAINS) + 0.5,
                 color=BAND, zorder=0)
    axis.text((min(FIT_DOMAINS) + max(FIT_DOMAINS)) / 2, 0.02,
              f"fit window, D{min(FIT_DOMAINS)}−D{max(FIT_DOMAINS)}",
              ha="center", va="bottom", fontsize=7.5, color=INK_MUTED,
              transform=axis.get_xaxis_transform(), zorder=6)
    axis.axvline(5.5, color=INK_MUTED, linestyle=(0, (4, 2)), linewidth=0.8,
                 zorder=2)
    axis.text(5.58, 0.97, "Buxton groin", rotation=90, fontsize=7,
              color=INK_MUTED, va="top",
              transform=axis.get_xaxis_transform())

    axis.plot(x, centre(observed), marker="s", markersize=4.0, linestyle="--",
              color=INK, linewidth=1.8, zorder=10,
              label="observed, 1984 to 2004")
    axis.plot(x, centre(baseline.profile), color=BASELINE_COLOR, linestyle=":",
              linewidth=1.4, zorder=4,
              label=f"no groin, {baseline.score:.1f} m")

    # One colour family for the top cells, darkest is best, so they read as
    # variations of one thing. viridis was used here until 2026-09-11.
    ramp = LinearSegmentedColormap.from_list(
        "hat_accent_ramp", [C["ACCENT_FILL"], C["ACCENT"]])
    for colour, (_, row) in zip(ramp(np.linspace(1.0, 0.2, len(top))),
                                top.iterrows()):
        marker = "o" if (row.M == CHOSEN_M and row.f == CHOSEN_F) else None
        axis.plot(x, centre(row.profile), color=colour, linewidth=1.4,
                  marker=marker, markersize=3.4, zorder=5,
                  label=f"M {row.M:g}, f {row.f:.2f}, {row.score:.1f} m")

    chosen = frame[(frame.M == CHOSEN_M) & (frame.f == CHOSEN_F)]
    if not chosen.empty and not ((top.M == CHOSEN_M) & (top.f == CHOSEN_F)).any():
        row = chosen.iloc[0]
        axis.plot(x, centre(row.profile), color=CHOSEN_COLOR, linewidth=1.8,
                  marker="o", markersize=3.4, zorder=8,
                  label=f"the pair taken forward, M {CHOSEN_M:g}, "
                        f"f {CHOSEN_F:g}, {row.score:.1f} m")

    axis.set_xticks(SHOW_DOMAINS)
    axis.set_xlabel("GIS domain")
    axis.set_ylabel("shoreline change 1984 to 2004 (m)\ndemeaned over the fit"
                    " window; positive is landward")
    spread = float(top.score.max() - top.score.min())
    axis.set_title(f"The top {len(top)} cells against the observed change,"
                   " period 1", loc="left")
    axis.grid(axis="y")
    axis.set_axisbelow(True)
    open_frame(axis)
    figure.legend(loc="outside lower center", ncol=4, frameon=False,
                  fontsize=7.5)

    caption(figure,
            "The direct counterpart to the 1967 rig's top-N figure, so the two "
            "calibrations can be read side by side, and the question is the "
            "same: do the best-scoring cells reproduce the observed alongshore "
            "SHAPE, or only a summary number? Three differences from the rig "
            "figure are deliberate. The window is 1984 to 2004, the only "
            "window in the hindcast where the observed gap between the groin's "
            "flanks widens, which is all a module with trapping at or above "
            "zero can produce. The geometry is the 120-domain production grid "
            "rather than the rig's 41, and M is grid-specific — a confined "
            "array preserves dipole amplitude that an open one diffuses away — "
            "so a value fitted here transfers to the hindcast and one fitted "
            "on the rig does not. And the score is demeaned and ranked on "
            "D4−D8 only, because a uniform level offset near the groin is "
            "absorbed by the source/sink calibration downstream and D1's "
            "81−104 m of cape change is about five times the groin's signal. "
            "The no-groin baseline is the grey dotted line: the groin is doing "
            "real work, taking the chosen cell from {base:.2f} to {chosen} m. "
            "The top five bundle within {spread:.2f} m while spanning M {mlo:g} "
            "to {mhi:g} and f {flo:g} to {fhi:g}, which is the visual "
            "statement of the ridge in period-1 cumulative trapping, "
            "M(15.5 + 4.5f) — NOT in M·f, which was tested and refuted at a "
            "correlation of −0.07. Everything is demeaned over the fit window, "
            "the same way the score centres it, and drawn landward-positive so "
            "erosion is up."
            .format(base=baseline.score,
                    chosen="{:.2f}".format(
                        float(chosen.iloc[0].score)) if not chosen.empty
                    else "{:.2f}".format(float(top.score.min())),
                    spread=spread, mlo=top.M.min(), mhi=top.M.max(),
                    flo=top.f.min(), fhi=top.f.max()))

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    written = save(figure, FIGURE_DIR / "period1_top_n_profiles.png",
                   close=True)
    print(f"wrote {written[0]}")
    print(f"  no groin {baseline.score:.2f} m")
    print(top[["M", "f", "score"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
