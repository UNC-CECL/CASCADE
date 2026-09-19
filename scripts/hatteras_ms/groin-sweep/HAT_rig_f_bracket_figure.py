#!/usr/bin/env python3
"""What the 1967 rig can and cannot say: f is bracketed, M is railed.

WHY THIS FIGURE EXISTS
    The rig is quoted as corroborating BOTH parameters. It does not, and the
    difference matters enough to draw. `scripts/site_layer/hatteras_site_config.py` said
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

SWEEP_CSV = (PROJECT_BASE_DIR / "hard-structures" / "groin"
             / "HAT-buxton-hindcast-groin-test" / "sensitivity_sweep"
             / "HAT_groin_sweep_results.csv")
FIGURE_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

if str(PROJECT_BASE_DIR / "scripts") not in sys.path:
    sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))
from site_layer.hat_figure_style import (apply_style, C, C_1984, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)

CHOSEN_M, CHOSEN_F = 60.0, 0.6
# House colours (2026-09-11). The series under test is the ACCENT; a cell
# that did not complete is drawn in the vintage red, the one warm colour in
# the palette, because "this run does not exist" has to be findable at a
# glance; the f families in panel (b) are a ramp of the accent so they read
# as one family. GOOD/BAD/INK/GRID were an orange, a pink and a navy chosen
# in this file.
DEAD = C_1984


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

    apply_style()
    from matplotlib.colors import LinearSegmentedColormap
    figure, (ax_f, ax_m) = plt.subplots(
        1, 2, figsize=figsize("double", aspect=0.44),
        gridspec_kw={"width_ratios": [1, 1.25]}, constrained_layout=True)

    # ---- (a) f at M = 60 -------------------------------------------------
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

    # ---- (b) M across the grid -------------------------------------------
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
    # Outside: eight entries inside panel (b) sat on top of the low-M cells,
    # which are the ones the panel is about.
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
