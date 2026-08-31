#!/usr/bin/env python3
"""Why M = 60, f = 0.6 -- the period-1 fit that supports it.

REPLACES AN EARLIER VERSION OF THIS FIGURE, and the reason matters.

    The previous version argued that the fillet metric was unsatisfiable, so M
    had to be bounded by physics (sediment budget, barrier stability) with the
    reach fit refining inside those bounds. Two things were wrong with it:

      1. The premise is superseded. M = 60 is now supported by a DIRECT FIT on
         period 1 over D4-D8 -- a 25% improvement on no groin -- rather than by
         an argument from constraints.

      2. It shaded "unstable M >= 70" and "barrier drowns M >= 100" on a
         PRODUCTION-geometry plot. Those thresholds were measured on the
         41-domain rig and do not transfer: all 36 production cells, including
         the full M = 70 and M = 80 rows, ran clean. Shading them was
         misleading.

WHY PERIOD 1, AND WHY D4-D8
    Period 1 is the only window in the hindcast where the observed gap between
    the groin's two flanks WIDENS (+52 m). That is the only behaviour the
    module can produce -- trapping is bounded at >= 0, so it can widen the gap
    or stop widening it, never close it. Period 2 closes, and the continuous
    1984-2024 window nets the two against each other, so neither can fit a
    groin however it is scored.

    D4-D8 excludes D1, where the cape's shoreline change over period 1 is
    81-104 m -- roughly five times the groin's ~17 m signal. On the full window
    with a raw score the cape swamps the groin and no-groin wins by 0.18 m; on
    D4-D8 the groin wins by 5.06 m. The signal was always there; the window was
    hiding it.

WHY THE SCORE IS DEMEANED
    A uniform level offset in the groin's neighbourhood is absorbed by the
    source/sink calibration that runs afterwards, so correcting it is not the
    groin's job. What the groin must get right is the SHAPE. Demeaning removes
    a constant and keeps every gradient -- unlike a linear detrend, which would
    partly absorb the dipole's own gradient and hide a working groin.

WHAT THE FIGURE DOES NOT CLAIM
    Not that (60, 0.6) is uniquely determined. Fourteen cells lie within 0.5 m
    of the best, spanning M = 40-95 and f = 0.4-1.0. The ridge is in
    PERIOD-1 CUMULATIVE TRAPPING, M(15.5 + 4.5f) -- not in M*f, which
    fig_Mf_identifiability.png tested and refuted (corr(RMSE, M*f) = -0.07).
    See CALIBRATION_FIGURES.md. The chosen pair is not the top-scoring one
    (M = 50, f = 1.0 scores 11.59 m, 0.10 m better) because f = 1.0 asserts the
    groin never deteriorated, which the GIS record contradicts outright. Within
    the tied band, (60, 0.6) is the cell that carries a real deterioration floor
    and still stays inside the affordable drift. The right panel draws the band rather
    than a single star, because a star would assert precision the data does not
    support.

Usage:
    python HAT_groin_choice_figure.py

Writes output/groin_sweep/figures/why_M60_f06.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline.hindcast import implied_interception_m3_yr  # noqa: E402
from hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402

from HAT_fullperiod_target import observed_change_profile  # noqa: E402

SWEEP = PROJECT_BASE_DIR / "output" / "groin_sweep" / "1984_2004_edgeBE"
FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

FIT_DOMAINS = list(range(4, 9))          # D4-D8: excludes the cape at D1
PINNED_BE1 = -40.0                       # nearest grid value to production's -41.8
CHOSEN_M, CHOSEN_F = 60.0, 0.6
BAND_M = 0.5                             # cells within this of the best are tied
PROFILE_HEIGHT_M = 1.7 + 22.25
DRIFT_LOW, DRIFT_HIGH = 5.0e5, 7.0e5


def load_cells():
    """Every period-1 cell at the pinned be1, scored on the demeaned profile."""
    # LANDWARD-POSITIVE, so erosion is UP in the left panel and it reads as a
    # plan view, matching the gifs and the other profile figures.
    # observed_change_profile is SEAWARD-positive, so it is negated; x_s is
    # landward-positive already, so the negation on `change` below is gone.
    # Both flip together, so every RMSE here is unchanged.
    observed = -np.array([observed_change_profile(1984, 2004, FIT_DOMAINS)[d]
                          for d in FIT_DOMAINS])
    observed_shape = observed - observed.mean()

    rows = []
    for cell in sorted(SWEEP.iterdir()):
        path = cell / "shoreline_matrix.npy"
        if not path.exists():
            continue
        name = cell.name
        be1 = np.nan if "beNA" in name else float(name.split("_be")[1].split("_f")[0])
        if not (np.isnan(be1) or be1 == PINNED_BE1):
            continue
        matrix = np.load(path)
        change = matrix[-1] - matrix[0]             # x_s is landward-positive
        model = np.array([change[GEOMETRY.gis_to_pad(d)] for d in FIT_DOMAINS])
        M = 0.0 if name.startswith("M0_") else float(name.split("_")[0][1:])
        fraction = 0.0 if name.startswith("M0_") else float(name.split("_f")[1])
        rows.append(dict(
            M=M, f=fraction,
            demeaned=float(np.sqrt(np.mean((model - model.mean() - observed_shape) ** 2))),
            profile=model))
    return pd.DataFrame(rows), observed


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame, observed = load_cells()
    groin = frame[frame.M > 0].sort_values("demeaned").reset_index(drop=True)
    baseline = frame[frame.M == 0].sort_values("demeaned").iloc[0]
    best = groin.iloc[0]
    chosen = frame[(frame.M == CHOSEN_M) & (frame.f == CHOSEN_F)].iloc[0]
    band = groin[groin.demeaned <= best.demeaned + BAND_M]

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure, (left, right) = plt.subplots(1, 2, figsize=(15, 6.5))

    # ---- LEFT: the profiles this fit is scored on ------------------------
    x = np.array(FIT_DOMAINS, dtype=float)
    centre = lambda v: np.asarray(v, float) - np.mean(v)
    left.axvspan(4.5, 5.5, color="#1565C0", alpha=0.10, zorder=0)
    left.axvspan(5.5, 6.5, color="#B71C1C", alpha=0.10, zorder=0)
    left.axvline(5.5, color="#B71C1C", linestyle="--", linewidth=1.8, zorder=2)
    left.text(5.58, 0.97, "groin", rotation=90, fontsize=8.5, color="#B71C1C",
              va="top", transform=left.get_xaxis_transform())

    left.plot(x, centre(observed), marker="s", markersize=8, linestyle="--",
              color="#1A1A1A", linewidth=2.6, zorder=6, label="OBSERVED 1984-2004")
    left.plot(x, centre(baseline.profile), marker="^", color="#777777",
              linewidth=2.0, linestyle=":", zorder=4,
              label=f"no groin  (RMSE {baseline.demeaned:.2f} m)")
    left.plot(x, centre(chosen.profile), marker="o", color="#FF8C00",
              linewidth=2.6, zorder=5,
              label=f"M={CHOSEN_M:g}, f={CHOSEN_F:g}  (RMSE {chosen.demeaned:.2f} m)")

    left.set_xticks(FIT_DOMAINS)
    left.set_xlabel("GIS domain")
    left.set_ylabel("shoreline change 1984-2004, demeaned (m)"
                    "   [+ = landward, erosion]")
    left.set_title(
        f"THE FIT\nperiod 1, D4-D8, shape only -- the groin closes "
        f"{(baseline.demeaned - chosen.demeaned) / baseline.demeaned * 100:.0f}% "
        f"of the misfit", fontsize=11.5)
    left.grid(alpha=0.25)
    left.legend(loc="best", fontsize=9)

    # ---- RIGHT: the M-f surface, with the indistinguishable band ---------
    grid = groin.pivot_table(index="f", columns="M", values="demeaned")
    mesh = right.pcolormesh(grid.columns, grid.index, grid.values,
                            shading="nearest", cmap="viridis_r")
    figure.colorbar(mesh, ax=right, label="demeaned profile RMSE, D4-D8 (m)")

    right.plot(band.M, band.f, marker="o", markersize=9, linestyle="none",
               markerfacecolor="none", markeredgecolor="white",
               markeredgewidth=2.0, zorder=5,
               label=f"{len(band)} cells within {BAND_M} m -- indistinguishable")
    right.plot(chosen.M, chosen.f, marker="*", markersize=26, color="#B71C1C",
               markeredgecolor="white", markeredgewidth=1.5, linestyle="none",
               zorder=7, label=f"CHOSEN M={CHOSEN_M:g}, f={CHOSEN_F:g}")

    # affordability, the one physical constraint that DOES apply on this grid
    m_axis = np.array(sorted(grid.columns), dtype=float)
    afford = [implied_interception_m3_yr(float(m), PROFILE_HEIGHT_M, GEOMETRY)
              for m in m_axis]
    m_lo = float(np.interp(DRIFT_LOW, afford, m_axis))
    m_hi = float(np.interp(DRIFT_HIGH, afford, m_axis))
    for edge in (m_lo, m_hi):
        right.axvline(edge, color="#1B5E20", linestyle="--", linewidth=1.6, zorder=4)
    right.text((m_lo + m_hi) / 2, 0.02,
               f"affordable\n{m_lo:.0f}-{m_hi:.0f} m/yr",
               ha="center", va="bottom", fontsize=8.5, color="#1B5E20",
               transform=right.get_xaxis_transform(), zorder=6)

    right.set_xlabel("groin trapping rate M (m/yr)")
    right.set_ylabel("deterioration floor f")
    right.set_title(
        "THE UNCERTAINTY\nM and f trade off along a ridge of constant "
        "period-1 trapping, M(15.5 + 4.5f)", fontsize=11.5)
    right.legend(loc="upper right", fontsize=8.5)

    figure.tight_layout(rect=(0, 0.10, 1, 1))
    figure.text(
        0.01, 0.012,
        f"PERIOD 1 is the only hindcast window where the observed gap between the groin's flanks WIDENS (+52 m), which is the only "
        f"behaviour a module with trapping >= 0 can produce. D4-D8 excludes D1, where the cape's change over period 1 is 81-104 m -- "
        f"about five times the groin's signal; on the full window with a raw score the cape swamps it and no-groin wins by 0.18 m. "
        f"The score is DEMEANED because a uniform level offset is absorbed by the source/sink calibration downstream, so only shape "
        f"is the groin's responsibility.\n"
        f"NO STABILITY SHADING HERE, unlike the previous version of this figure: the M>=70 instability and M>=100 drowning were "
        f"measured on the 41-domain rig and do NOT transfer -- all 36 production cells including M=70 and M=80 ran clean. "
        f"Affordability is the one physical constraint that does apply.",
        fontsize=7.5, color="#333333", wrap=True)

    path = FIGURE_DIR / "why_M60_f06.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    print(f"wrote {path}")
    print(f"  no groin        {baseline.demeaned:.2f} m")
    print(f"  chosen ({CHOSEN_M:g},{CHOSEN_F:g}) {chosen.demeaned:.2f} m  "
          f"({(baseline.demeaned - chosen.demeaned) / baseline.demeaned * 100:.0f}% better)")
    print(f"  best  ({best.M:g},{best.f:g})  {best.demeaned:.2f} m")
    print(f"  {len(band)} cells within {BAND_M} m: M {band.M.min():g}-{band.M.max():g}, "
          f"f {band.f.min():g}-{band.f.max():g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
