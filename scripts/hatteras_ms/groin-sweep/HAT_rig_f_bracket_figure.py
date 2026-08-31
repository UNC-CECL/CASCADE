#!/usr/bin/env python3
"""What the 1967 rig can and cannot say: f is bracketed, M is railed.

WHY THIS FIGURE EXISTS
    The rig is quoted as corroborating BOTH parameters. It does not, and the
    difference matters enough to draw. `scripts/hatteras_site_config.py` said
    until 2026-08-30 that in the rig "f = 0.6 is a clean INTERIOR minimum ...
    So is M ... Neither railed." The rig's own sweep CSV refuses the second
    half: RMSE improves monotonically to M = 60 and then the model blows up.

    So the rig resolves f and is only CONSISTENT with M. M is set by the
    production period-1 fit instead. This figure is the evidence for both
    halves of that sentence, in one place, so the claim cannot drift back.

WHAT IS PLOTTED
    (a) f AT M = 60 -- a clean interior minimum at f = 0.6, bracketed on both
        sides with steep curvature. This is the parameter the rig owns: it is
        the only window containing the 1996-2003 deterioration ramp, because
        both hindcast windows begin 15 years after the structure went in.

    (b) M ACROSS THE WHOLE GRID, on a log axis because the failure spans three
        orders of magnitude. Every f-series improves monotonically to M = 60
        and then jumps ~13x at M = 70. Cells at M >= 100 do not complete at
        all and are drawn on the "crashed" rule at the top. M = 60 is the LAST
        VALUE THAT RUNS, not the value where the fit stops improving.

    The stability wall is RIG-SPECIFIC and does not transfer: all 36 production
    cells, including the full M = 70 and M = 80 rows, ran clean on the
    120-domain grid. Do not quote this ceiling for production runs.

Usage:
    python HAT_rig_f_bracket_figure.py

Writes output/groin_sweep/figures/rig_f_bracket.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")

SWEEP_CSV = (PROJECT_BASE_DIR / "hard-structures" / "groin"
             / "HAT-hindcast-groin-test" / "sensitivity_sweep"
             / "HAT_groin_sweep_results.csv")
FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

CHOSEN_M, CHOSEN_F = 60.0, 0.6
GOOD, BAD, INK, GRID = "#FF8C00", "#C2185B", "#1A1A2E", "#D5D8DD"


def main() -> None:
    if not SWEEP_CSV.exists():
        raise SystemExit(f"missing rig sweep results: {SWEEP_CSV}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frame = pd.read_csv(SWEEP_CSV)
    frame["rmse"] = pd.to_numeric(frame["rmse"], errors="coerce")
    # A blank rmse is a cell that did not complete -- the barrier drowned or
    # the solver diverged. Those are information, not missing data.
    crashed = frame[frame["rmse"].isna()]
    ok = frame.dropna(subset=["rmse"])

    plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.7,
                         "legend.frameon": False, "pdf.fonttype": 42})
    figure, (ax_f, ax_m) = plt.subplots(1, 2, figsize=(13.5, 5.4),
                                        gridspec_kw={"width_ratios": [1, 1.25]})

    # ---- (a) f at M = 60 -------------------------------------------------
    series = ok[ok["M"] == CHOSEN_M].sort_values("fraction")
    ax_f.plot(series["fraction"], series["rmse"], "-o", color=GOOD,
              linewidth=2.2, markersize=6, zorder=4)
    best = series.loc[series["rmse"].idxmin()]
    ax_f.plot([best["fraction"]], [best["rmse"]], "o", markersize=14,
              markerfacecolor="none", markeredgecolor=INK,
              markeredgewidth=2.0, zorder=5)
    ax_f.annotate(f"f = {best['fraction']:g}\nRMSE {best['rmse']:.2f} m",
                  xy=(best["fraction"], best["rmse"]),
                  xytext=(best["fraction"] + 0.06, best["rmse"] + 6.0),
                  fontsize=9.5, color=INK,
                  arrowprops=dict(arrowstyle="-", color=INK, linewidth=0.9))
    ax_f.set_xlabel("deterioration fraction f")
    ax_f.set_ylabel("Rig RMSE, D2–D12 change profile (m)")
    ax_f.set_title("(a)  f IS bracketed  —  at M = 60, 1967–2018",
                   fontsize=11, loc="left")
    ax_f.grid(color=GRID, linewidth=0.6)
    ax_f.set_axisbelow(True)
    for side in ("top", "right"):
        ax_f.spines[side].set_visible(False)
    ax_f.annotate("the only window containing the\n1996–2003 deterioration ramp",
                  xy=(0.03, 0.06), xycoords="axes fraction", ha="left",
                  va="bottom", fontsize=8.8, color="#555555")

    # ---- (b) M across the grid -------------------------------------------
    ceiling = ok["rmse"].max() * 2.4
    for fraction, group in ok.groupby("fraction"):
        group = group.sort_values("M")
        ax_m.plot(group["M"], group["rmse"], "-o", linewidth=1.5,
                  markersize=4.5, alpha=0.85, label=f"f = {fraction:g}",
                  zorder=3)

    if not crashed.empty:
        ax_m.plot(crashed["M"], np.full(len(crashed), ceiling), "x",
                  color=BAD, markersize=8, markeredgewidth=1.8,
                  label="did not complete", zorder=4)
        ax_m.axhline(ceiling, color=BAD, linewidth=0.9,
                     linestyle=(0, (3, 3)), alpha=0.6, zorder=2)
        ax_m.annotate("CRASHED — barrier drowned / solver diverged",
                      xy=(0.985, ceiling), xycoords=("axes fraction", "data"),
                      ha="right", va="bottom", fontsize=8.5, color=BAD)

    ax_m.axvspan(65, max(ok["M"].max(), crashed["M"].max() if not crashed.empty
                         else 0) + 8, color="#FDECEF", zorder=0)
    ax_m.axvline(CHOSEN_M, color=INK, linewidth=1.3, linestyle=(0, (4, 2)),
                 zorder=3)
    ax_m.annotate("M = 60\nlast value that runs",
                  xy=(CHOSEN_M, 0.42), xycoords=("data", "axes fraction"),
                  ha="right", va="center", fontsize=9.0, color=INK)
    ax_m.annotate("STABILITY WALL\n(rig-specific — production\nran M = 70 and 80 clean)",
                  xy=(0.985, 0.30), xycoords="axes fraction", ha="right",
                  va="center", fontsize=8.6, color=BAD)

    ax_m.set_yscale("log")
    ax_m.set_xlabel("trapping rate M (m/yr)")
    ax_m.set_ylabel("Rig RMSE (m, log scale)")
    ax_m.set_title("(b)  M is NOT — it rails against a stability wall",
                   fontsize=11, loc="left")
    ax_m.grid(color=GRID, linewidth=0.6, which="both")
    ax_m.set_axisbelow(True)
    for side in ("top", "right"):
        ax_m.spines[side].set_visible(False)
    ax_m.legend(fontsize=8.2, ncol=2, loc="lower left")

    figure.text(
        0.012, 0.015,
        "The rig therefore CORROBORATES f = 0.6 and is only CONSISTENT WITH "
        "M = 60: RMSE improves monotonically right up to the last cell that "
        "completes, so a larger M cannot be ruled out here. M is set instead "
        "by the production period-1 D4–D8 fit, and reproduced independently by "
        "D3–D9. Corrects an earlier claim that the rig bracketed both.",
        fontsize=7.8, color="#555555", wrap=True, va="bottom")

    figure.subplots_adjust(top=0.92, bottom=0.20, left=0.065, right=0.985,
                           wspace=0.22)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FIGURE_DIR / "rig_f_bracket.png"
    figure.savefig(out, dpi=200)
    print(f"wrote {out}")
    print(f"  f at M=60: best f = {best['fraction']:g}, RMSE {best['rmse']:.2f}")
    print(f"  cells that did not complete: {len(crashed)} "
          f"(M >= {crashed['M'].min():g})" if not crashed.empty else "")


if __name__ == "__main__":
    main()
