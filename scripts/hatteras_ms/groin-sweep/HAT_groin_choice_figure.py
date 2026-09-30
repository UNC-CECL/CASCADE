#!/usr/bin/env python3
"""
Why M = 60, f = 0.6: the period-1 fit that supports it.

    python scripts/hatteras_ms/groin-sweep/HAT_groin_choice_figure.py

The profile fit and the ridge the chosen pair sits on. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from HAT_groin_sweep_config import GROIN_SWEEP_ROOT  # noqa: E402
from cascade_pipeline.hindcast import implied_interception_m3_yr  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402

from HAT_fullperiod_target import observed_change_profile  # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, error_cmap, figsize, open_frame, save,
                              _title)

# --- CONFIG ------------------------------------------------------------------
SWEEP = GROIN_SWEEP_ROOT / "1984_2004_edgeBE"
FIGURE_DIR = GROIN_SWEEP_ROOT / "figures"

FIT_DOMAINS = list(range(4, 9))          # D4-D8: excludes the cape at D1
PINNED_BE1 = -40.0                       # nearest grid value to production's -41.8
CHOSEN_M, CHOSEN_F = 60.0, 0.6
BAND_M = 0.5                             # cells within this of the best are tied
PROFILE_HEIGHT_M = 1.7 + 22.25
DRIFT_LOW, DRIFT_HIGH = 5.0e5, 7.0e5
# -----------------------------------------------------------------------------


# Every period-1 cell at the pinned be1, scored on the demeaned profile
def load_cells():
    # LANDWARD-POSITIVE, so erosion is UP in the left panel and it reads as a plan view
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


# Run: the figure
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
    apply_style()
    figure, (left, right) = plt.subplots(
        1, 2, figsize=figsize("double", aspect=0.44), constrained_layout=True)

    # Left: the profiles this fit is scored on
    x = np.array(FIT_DOMAINS, dtype=float)
    centre = lambda v: np.asarray(v, float) - np.mean(v)
    # The two flanks of the structure, named rather than colour-coded
    left.axvspan(4.5, 5.5, color="0.94", zorder=0)
    left.axvspan(5.5, 6.5, color="0.90", zorder=0)
    left.axvline(5.5, color=INK_MUTED, linestyle=(0, (4, 2)), linewidth=0.8,
                 zorder=2)
    left.text(5.42, 0.97, "downdrift", rotation=90, fontsize=7,
              color=INK_MUTED, ha="right", va="top",
              transform=left.get_xaxis_transform())
    left.text(5.58, 0.97, "updrift", rotation=90, fontsize=7,
              color=INK_MUTED, ha="left", va="top",
              transform=left.get_xaxis_transform())

    left.plot(x, centre(observed), marker="s", markersize=4.0, linestyle="--",
              color=INK, linewidth=1.8, zorder=6,
              label="observed, 1984 to 2004")
    left.plot(x, centre(baseline.profile), marker="^", markersize=3.4,
              color=C["BASE"], linewidth=1.4, linestyle=":", zorder=4,
              label=f"no groin, {baseline.demeaned:.2f} m")
    left.plot(x, centre(chosen.profile), marker="o", markersize=3.4,
              color=C["ACCENT"], linewidth=1.8, zorder=5,
              label=f"M {CHOSEN_M:g}, f {CHOSEN_F:g}, {chosen.demeaned:.2f} m")

    left.set_xticks(FIT_DOMAINS)
    left.set_xlabel("GIS domain")
    left.set_ylabel("shoreline change 1984 to 2004, demeaned (m)\n"
                    "positive is landward")
    _title(left, 0, "the fit, on shape alone")
    left.grid(axis="y")
    left.set_axisbelow(True)
    open_frame(left)
    left.legend(loc="lower left", fontsize=7.5)

    # Right: the M-f surface, with the indistinguishable band
    grid = groin.pivot_table(index="f", columns="M", values="demeaned")
    mesh = right.pcolormesh(grid.columns, grid.index, grid.values,
                            shading="nearest", cmap=error_cmap())
    cb = figure.colorbar(mesh, ax=right)
    cb.set_label("demeaned profile RMSE, D4−D8 (m)")
    cb.outline.set_linewidth(0.6)

    right.plot(band.M, band.f, marker="o", markersize=5.5, linestyle="none",
               markerfacecolor="none", markeredgecolor=C["ACCENT"],
               markeredgewidth=1.2, zorder=5,
               label=f"{len(band)} cells within {BAND_M} m of the best")
    right.plot(chosen.M, chosen.f, marker="*", markersize=14,
               color=C["ACCENT"], markeredgecolor="white",
               markeredgewidth=0.8, linestyle="none", zorder=7,
               label=f"the pair taken forward, M {CHOSEN_M:g}, f {CHOSEN_F:g}")

    # affordability, the one physical constraint that DOES apply on this grid
    m_axis = np.array(sorted(grid.columns), dtype=float)
    afford = [implied_interception_m3_yr(float(m), PROFILE_HEIGHT_M, GEOMETRY)
              for m in m_axis]
    m_lo = float(np.interp(DRIFT_LOW, afford, m_axis))
    m_hi = float(np.interp(DRIFT_HIGH, afford, m_axis))
    for edge in (m_lo, m_hi):
        right.axvline(edge, color=C["REF"], linestyle=(0, (4, 2)),
                      linewidth=1.0, zorder=4)
    right.text((m_lo + m_hi) / 2, 0.02,
               f"affordable, {m_lo:.0f} to {m_hi:.0f} m/yr",
               ha="center", va="bottom", fontsize=7, color=C["REF"],
               transform=right.get_xaxis_transform(), zorder=6,
               bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                         boxstyle="square,pad=0.15"))

    right.set_xlabel("groin trapping rate M (m/yr)")
    right.set_ylabel("deterioration floor f")
    _title(right, 1, "the uncertainty in the pair")
    right.legend(loc="upper right", fontsize=7)

    closed = (baseline.demeaned - chosen.demeaned) / baseline.demeaned * 100
    caption(figure,
            "Why M = {M:g}, f = {f:g}. (a) The fit the choice is made on: "
            "period 1, D4−D8, shape only. The groin closes {closed:.0f}% of "
            "the misfit the no-groin baseline carries, from {base:.2f} to "
            "{ch:.2f} m. Period 1 is the only hindcast window where the "
            "observed gap between the structure's flanks WIDENS, by 52 m, "
            "which is the only behaviour a module with trapping at or above "
            "zero can produce. D4−D8 excludes D1, where the cape's change over "
            "period 1 is 81 to 104 m, about five times the groin's signal; on "
            "the full window with a raw score the cape swamps it and no groin "
            "wins by 0.18 m. The score is demeaned because a uniform level "
            "offset is absorbed by the source/sink calibration downstream, so "
            "only shape is the groin's responsibility. The two light bands are "
            "the downdrift and updrift domains either side of the structure. "
            "(b) What that fit does and does not pin down. {n} cells sit "
            "within {bandm} m of the best and are not distinguishable by this "
            "target; M and f trade off along a ridge of constant period-1 "
            "cumulative trapping, M(15.5 + 4.5f). The green lines are "
            "affordability, the one physical constraint that does apply on "
            "this grid: the drift a groin of that trapping rate would have to "
            "intercept. There is deliberately NO stability shading, unlike an "
            "earlier version of this figure — the instability above M = 70 and "
            "the drowning above M = 100 were measured on the 41-domain rig and "
            "do not transfer, and all 36 production cells including M = 70 and "
            "M = 80 ran clean."
            .format(M=CHOSEN_M, f=CHOSEN_F, closed=closed,
                    base=baseline.demeaned, ch=chosen.demeaned,
                    n=len(band), bandm=BAND_M))

    written = save(figure, FIGURE_DIR / "why_M60_f06.png", close=True)
    print(f"wrote {written[0]}")
    print(f"  no groin        {baseline.demeaned:.2f} m")
    print(f"  chosen ({CHOSEN_M:g},{CHOSEN_F:g}) {chosen.demeaned:.2f} m  "
          f"({(baseline.demeaned - chosen.demeaned) / baseline.demeaned * 100:.0f}% better)")
    print(f"  best  ({best.M:g},{best.f:g})  {best.demeaned:.2f} m")
    print(f"  {len(band)} cells within {BAND_M} m: M {band.M.min():g}-{band.M.max():g}, "
          f"f {band.f.min():g}-{band.f.max():g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
