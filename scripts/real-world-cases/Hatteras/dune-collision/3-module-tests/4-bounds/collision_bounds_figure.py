"""
The dune collision module's tested range against every storm the Hatteras hindcast runs.

    python collision_bounds_figure.py
        ->  figures/collision_bounds.png, figures/supporting/storms_against_bounds.csv

Every storm of the 1996-2025 hindcast storm series (its summary file: peak water level R_high, peak period,
duration) against the limits DuneCollisionManager stops a run at (cascade/dune_collision.py: duration 200 h,
period 5-20 s, water 5 m above a toe it acts on; the Cs limit is a setting, not a storm property). The water
above the toe is taken at the lowest toe of the start each storm runs from (1996-start up to 2009, 2009-start
after, dune_toe_cells.csv as it stands), so it is the most any cell sees; the median toe is drawn too. Only
storms that can reach a toe are counted: with these toes, every storm. Takes seconds.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
_REPO = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
for _p in (_HAT_DIR, _REPO):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
from hatteras_paths import WORKSPACE, workspace_twin  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK_MUTED, _title, apply_style, caption, figsize, open_frame, save, support_dir)
from cascade.dune_collision import DuneCollisionManager as M  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
INIT = WORKSPACE / "data" / "hatteras_init"
STORMS = INIT / "3-env-forcings" / "3-storms" / "hindcast_storms" / "1996_2025" / "1996_2025_storms_v3_240_summary.csv"
TOES = {start: INIT / "1-barrier3d-domains" / start / "dune-toes" / "dune_toe_cells.csv"
        for start in ("1996-start", "2009-start")}
SPLIT_YEAR = 2009                 # storms up to this calendar year run from the 1996 start
MHW = 0.36                        # m NAVD88
MAX_DUR = M.MAX_STORM_DURATION_HR
T_LO, T_HI = M.STORM_PERIOD_RANGE_S
MAX_ABOVE_M = M.MAX_WATER_ABOVE_TOE_DAM * 10
OUT = workspace_twin(_HERE.parent) / "figures"
C_OUT = "#f2dede"                 # outside the tested range
# -----------------------------------------------------------------------------


# The storms, with the water above the lowest and the median toe of the start they run from
def storms():
    s = pd.read_csv(STORMS)
    s["r_m"] = s.Rhigh * 10 + MHW
    toe = {k: pd.read_csv(f).toe_m_navd88.dropna() for k, f in TOES.items()}
    start = np.where(s.calendar_year <= SPLIT_YEAR, "1996-start", "2009-start")
    s["start"] = start
    s["above_lowest_toe_m"] = s.r_m - np.array([toe[k].min() for k in start])
    s["above_median_toe_m"] = s.r_m - np.array([toe[k].median() for k in start])
    return s, toe


# Run: the storms, the table, the figure
def main() -> None:
    s, toe = storms()
    s.to_csv(support_dir(OUT) / "storms_against_bounds.csv", index=False)
    reach = s[s.above_lowest_toe_m > 0]
    near = dict(duration=reach.duration.max() / MAX_DUR, above=reach.above_lowest_toe_m.max() / MAX_ABOVE_M,
                t_hi=reach.period.max() / T_HI, t_lo=T_LO / reach.period.min())
    print(f"{len(s)} storms, {len(reach)} reach a toe; longest {reach.duration.max():.0f} h of {MAX_DUR:g} "
          f"({near['duration']:.0%}); most water above the lowest toe {reach.above_lowest_toe_m.max():.2f} m of "
          f"{MAX_ABOVE_M:g} ({near['above']:.0%}); periods {reach.period.min():.1f}-{reach.period.max():.1f} s "
          f"of {T_LO:g}-{T_HI:g}; outside any limit: "
          f"{int(((reach.duration > MAX_DUR) | (reach.above_lowest_toe_m > MAX_ABOVE_M) | (reach.period > T_HI) | (reach.period < T_LO)).sum())}")

    apply_style()
    fig, (a, b) = plt.subplots(1, 2, figsize=figsize("double", 0.42), layout="constrained")
    # (a) duration against the water above the toe
    a.axvspan(MAX_DUR, 260, color=C_OUT, lw=0)
    a.axhspan(MAX_ABOVE_M, 6.5, color=C_OUT, lw=0)
    for start, col, mk in (("1996-start", C["EARLY"], "o"), ("2009-start", C["LATE"], "s")):
        x = s[s.start == start]
        a.scatter(x.duration, x.above_lowest_toe_m, s=12, color=col, marker=mk, lw=0, alpha=0.8,
                  label=f"{'1996-2009' if start == '1996-start' else '2010-2025'} storms, lowest toe")
        a.scatter(x.duration, x.above_median_toe_m, s=10, facecolors="none", edgecolors=col, marker=mk, lw=0.6,
                  alpha=0.6)
    a.scatter([], [], s=10, facecolors="none", edgecolors="0.4", lw=0.6, label="the same, median toe")
    a.axvline(MAX_DUR, color=C["GROIN"], lw=1)
    a.axhline(MAX_ABOVE_M, color=C["GROIN"], lw=1)
    a.text(MAX_DUR + 4, 5.9, f"{MAX_DUR:g} h", color=C["GROIN"], fontsize=7.5)
    a.text(5, MAX_ABOVE_M + 0.12, f"{MAX_ABOVE_M:g} m above the toe", color=C["GROIN"], fontsize=7.5)
    big = reach.loc[reach.above_lowest_toe_m.idxmax()]
    a.annotate(f"Isabel {big.r_m:.2f} m", (big.duration, big.above_lowest_toe_m),
               xytext=(big.duration + 25, big.above_lowest_toe_m + 0.05), fontsize=7.5,
               arrowprops=dict(arrowstyle="-", lw=0.6))
    a.set(xlim=(0, 260), ylim=(-1.6, 6.5), xlabel="Storm duration (h)", ylabel="Peak water above the dune toe (m)")
    a.axhline(0, color=INK_MUTED, lw=0.6, ls=":")
    open_frame(a)
    _title(a, 0, "Duration and water level")

    # (b) duration against peak period
    b.axvspan(MAX_DUR, 260, color=C_OUT, lw=0)
    b.axhspan(T_HI, 24, color=C_OUT, lw=0)
    b.axhspan(0, T_LO, color=C_OUT, lw=0)
    for start, col, mk in (("1996-start", C["EARLY"], "o"), ("2009-start", C["LATE"], "s")):
        x = s[s.start == start]
        b.scatter(x.duration, x.period, s=12, color=col, marker=mk, lw=0, alpha=0.8)
    for y in (T_LO, T_HI):
        b.axhline(y, color=C["GROIN"], lw=1)
    b.axvline(MAX_DUR, color=C["GROIN"], lw=1)
    b.text(5, T_HI + 0.4, f"{T_HI:g} s", color=C["GROIN"], fontsize=7.5)
    b.text(5, T_LO - 1.4, f"{T_LO:g} s", color=C["GROIN"], fontsize=7.5)
    b.set(xlim=(0, 260), ylim=(0, 24), xlabel="Storm duration (h)", ylabel="Peak wave period (s)")
    open_frame(b)
    _title(b, 1, "Duration and wave period")
    h, lab = a.get_legend_handles_labels()
    fig.legend(h, lab, loc="outside lower center", ncol=3, frameon=False, fontsize=7.5)

    caption(fig, f"**Tests:** how close the Hatteras hindcast's storms come to the limits the dune collision "
            f"module stops a run at (`DuneCollisionOutsideTestedRange`). **How to read:** every storm of the "
            f"1996-2025 hindcast storm series ({len(s)}; red circles 1996-2009, run from the 1996 start; blue squares "
            f"2010-2025, from the 2009 start). (a) duration against the peak water level above the lowest dune toe of "
            f"that start (filled; {toe['1996-start'].min():.2f} and {toe['2009-start'].min():.2f} m NAVD88), the most "
            f"any cell sees, and above the median toe (open); (b) duration against the peak wave period. Red lines "
            f"and shading: the module's limits ({MAX_DUR:g} h, {T_LO:g}-{T_HI:g} s, {MAX_ABOVE_M:g} m above a toe). "
            f"The module checks a storm only where it lowers a dune, so a storm past a limit on cells it overtops "
            f"does not stop a run. **Shows:** every storm is inside the tested range, by little on two of the "
            f"limits: the longest is {reach.duration.max():.0f} h ({near['duration']:.0%} of the limit), the "
            f"highest puts the water {reach.above_lowest_toe_m.max():.2f} m above the lowest toe "
            f"({near['above']:.0%}, Isabel, which overtops most dunes); the periods span "
            f"{reach.period.min():.1f}-{reach.period.max():.1f} s, well inside {T_LO:g}-{T_HI:g} s. A longer "
            f"storm series (a forecast) is the likeliest to reach a limit, first on duration.")
    save(fig, OUT / "collision_bounds", close=True)


if __name__ == "__main__":
    main()
