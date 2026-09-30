#!/usr/bin/env python3
"""
What the 1967 rig can and cannot say: f is bracketed, M is railed.

    python scripts/hatteras_ms/groin-sweep/HAT_rig_f_bracket_figure.py

From the rig sweep's results. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")

# --- CONFIG ------------------------------------------------------------------
SWEEP_CSV = (PROJECT_BASE_DIR / "hard-structures" / "groin"
             / "HAT-buxton-hindcast-groin-test" / "sensitivity_sweep"
             / "HAT_groin_sweep_results.csv")
FIGURE_DIR = GROIN_SWEEP_ROOT / "figures"

if str(PROJECT_BASE_DIR / "scripts") not in sys.path:
    sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))
from HAT_groin_sweep_config import GROIN_SWEEP_ROOT  # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, C_1984, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)

CHOSEN_M, CHOSEN_F = 60.0, 0.6
# House colours (2026-09-11)
DEAD = C_1984
# -----------------------------------------------------------------------------


# Run: the figure
def main() -> None:
    if not SWEEP_CSV.exists():
        raise SystemExit(f"missing rig sweep results: {SWEEP_CSV}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame = pd.read_csv(SWEEP_CSV)
    frame["rmse"] = pd.to_numeric(frame["rmse"], errors="coerce")
    # A blank rmse is a cell that did not complete -- the barrier drowned or the solver diverged
    crashed = frame[frame["rmse"].isna()]
    ok = frame.dropna(subset=["rmse"])

    apply_style()
    from matplotlib.colors import LinearSegmentedColormap
    figure, (ax_f, ax_m) = plt.subplots(
        1, 2, figsize=figsize("double", aspect=0.44),
        gridspec_kw={"width_ratios": [1, 1.25]}, constrained_layout=True)

    # (a) f at m = 60
    series = ok[ok["M"] == CHOSEN_M].sort_values("fraction")
    ax_f.plot(series["fraction"], series["rmse"], "-o", color=C["ACCENT"],
              linewidth=1.6, markersize=3.4, zorder=4)
    best = series.loc[series["rmse"].idxmin()]
    ax_f.plot([best["fraction"]], [best["rmse"]], "o", markersize=9,
              markerfacecolor="none", markeredgecolor=INK,
              markeredgewidth=1.2, zorder=5)
    ax_f.annotate(f"f = {best['fraction']:g}, {best['rmse']:.2f} m",
                  xy=(best["fraction"], best["rmse"]),
                  xytext=(best["fraction"] + 0.06, best["rmse"] + 6.0),
                  fontsize=7.5, color=INK,
                  arrowprops=dict(arrowstyle="-", color=INK, linewidth=0.7))
    ax_f.set_xlabel("deterioration fraction f")
    ax_f.set_ylabel("rig RMSE, D2–D12 change profile (m)")
    _title(ax_f, 0, "f at M = 60, over 1967 to 2018")
    ax_f.grid()
    ax_f.set_axisbelow(True)
    open_frame(ax_f)

    # (b) m across the grid
    ceiling = ok["rmse"].max() * 2.4
    fractions = sorted(ok["fraction"].unique())
    ramp = LinearSegmentedColormap.from_list(
        "hat_accent_ramp", [C["ACCENT_FILL"], C["ACCENT"]])
    shades = dict(zip(fractions, ramp(np.linspace(0.2, 1.0, len(fractions)))))
    for fraction, group in ok.groupby("fraction"):
        group = group.sort_values("M")
        ax_m.plot(group["M"], group["rmse"], "-o", linewidth=1.2,
                  markersize=2.6, color=shades[fraction],
                  label=f"f = {fraction:g}", zorder=3)

    if not crashed.empty:
        ax_m.plot(crashed["M"], np.full(len(crashed), ceiling), "x",
                  color=DEAD, markersize=5, markeredgewidth=1.2,
                  label="did not complete", zorder=4)
        ax_m.axhline(ceiling, color=DEAD, linewidth=0.7,
                     linestyle=(0, (3, 3)), zorder=2)

    ax_m.axvspan(65, max(ok["M"].max(), crashed["M"].max() if not crashed.empty
                         else 0) + 8, color="0.94", zorder=0)
    ax_m.axvline(CHOSEN_M, color=INK_MUTED, linewidth=0.8,
                 linestyle=(0, (4, 2)), zorder=3)
    ax_m.annotate("M = 60, the last value\nthat runs on the rig",
                  xy=(CHOSEN_M, 0.80), xycoords=("data", "axes fraction"),
                  ha="right", va="center", fontsize=7, color=INK_MUTED)
    ax_m.annotate("cells here did not complete",
                  xy=(0.985, 0.16), xycoords="axes fraction", ha="right",
                  va="center", fontsize=7, color=INK_MUTED)

    ax_m.set_yscale("log")
    ax_m.set_xlabel("trapping rate M (m/yr)")
    ax_m.set_ylabel("rig RMSE (m, log scale)")
    _title(ax_m, 1, "M across the grid")
    ax_m.grid(which="both")
    ax_m.set_axisbelow(True)
    open_frame(ax_m)
    # Outside
    figure.legend(*ax_m.get_legend_handles_labels(),
                  loc="outside lower center", ncol=8, frameon=False,
                  fontsize=7)

    caption(figure,
            "What the 1967 rig does and does not settle. (a) f IS bracketed: "
            "at M = 60 over 1967 to 2018 the rig's error has an interior "
            "minimum at f = {bf:g}, and this is the only window that contains "
            "the 1996 to 2003 deterioration ramp, which is what makes f "
            "visible at all. (b) M is NOT. Every f family improves "
            "monotonically in M right up to the last cell that completes, and "
            "past M = 65 cells stop completing at all, the barrier drowning or "
            "the solver diverging, so the apparent optimum at M = 60 is a "
            "stability wall rather than a fit. Those failures are information, "
            "not missing data, so they are drawn on a ceiling line instead of "
            "being dropped. The wall is RIG-SPECIFIC and does not transfer: "
            "the production grid ran M = 70 and M = 80 clean. The rig "
            "therefore corroborates f = {cf:g} and is only consistent with "
            "M = {cm:g}, which is set instead by the production period-1 "
            "D4–D8 fit and reproduced independently on D3–D9. This corrects "
            "an earlier claim that the rig bracketed both."
            .format(bf=best["fraction"], cf=CHOSEN_F, cm=CHOSEN_M))

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    written = save(figure, FIGURE_DIR / "rig_f_bracket.png", close=True)
    print(f"wrote {written[0]}")
    print(f"  f at M=60: best f = {best['fraction']:g}, RMSE {best['rmse']:.2f}")
    print(f"  cells that did not complete: {len(crashed)} "
          f"(M >= {crashed['M'].min():g})" if not crashed.empty else "")


if __name__ == "__main__":
    main()
