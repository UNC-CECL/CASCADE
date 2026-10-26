"""
Can the dune ceiling be the measured dune height, if the growth rate is slowed to the slow-growth candidate?

    python scripts/hatteras_ms/experiments/HAT_measured_ceiling_slow_growth.py run
    python scripts/hatteras_ms/experiments/HAT_measured_ceiling_slow_growth.py score
    python scripts/hatteras_ms/experiments/HAT_measured_ceiling_slow_growth.py ladder   # 12.1 m, slower rates, 1996; own table
    python scripts/hatteras_ms/experiments/HAT_measured_ceiling_slow_growth.py advisor  # Dmax 10.5 m above the berm; own table

Dmaxel at the 1996 DEM's 99th-percentile crest (9.03 m NAVD88) and at its tallest
crest (12.1 m), both with rmin-rmax 0.03-0.07, groin on, on 1996-2009 and
2009-2025, in the low growth-rate study's folder through the grid study's
launcher. Scored exactly as HAT_low_growth_test_period.py scores the slow-growth
candidate (5.0 m, 0.03-0.07) and the control, so the rows sit beside its table.
Asked by Hannah 2026-10-20: set the ceiling to the measured maximum and let the
growth rate keep the dunes from growing too quickly. "ladder" adds slower rates
at 12.1 m from the 1996 start only (LADDER_12), asked the same day: is there a rate
at which the 12.1 m ceiling's 2010 dunes match the 2009 lidar? "advisor" sets Dmax
(Barrier3D's Dmaxel - BermEl) to 10.5 m above the berm, as Hannah's advisor asked, with
the rate that matched the 2009 median (0-0.01), on both windows.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_low_growth_test_period as T  # noqa: E402
L, G, U = T.L, T.G, T.U

# --- CONFIG ------------------------------------------------------------------
SETTINGS = {"0.03-0.07, 9.03 m (DEM P99)": (9.03, (0.03, 0.07)),
            "0.03-0.07, 12.1 m (DEM max)": (12.1, (0.03, 0.07))}
LADDER_12 = {f"{lo:.2f}-{hi:.2f}, 12.1 m (DEM max)": (12.1, (lo, hi))
             for lo, hi in ((0.0, 0.0), (0.0, 0.01), (0.0, 0.02), (0.01, 0.01), (0.01, 0.02), (0.01, 0.03),
                            (0.02, 0.04))}   # 0-0.01 rungs added when 0.01-0.02 came closest
BERM_NAVD88 = 1.7            # m; Barrier3D BermEl, the same in every domain
DMAX_ABOVE_BERM = 10.5       # m; the advisor's Dmax
ADVISOR = {"0.00-0.01, Dmax 10.5 m above the berm": (round(BERM_NAVD88 + DMAX_ABOVE_BERM, 2), (0.0, 0.01))}
OUT_ADVISOR = L.EXP_DIR / "tables" / "dmax_10p5_above_berm.csv"
OUT = L.EXP_DIR / "tables" / "measured_ceiling_slow_growth.csv"
OUT_LADDER = L.EXP_DIR / "tables" / "measured_ceiling_12m_ladder_1996.csv"
# -----------------------------------------------------------------------------


def jobs():
    return [(s, G.name(c, r, groin=True)) for s in (1996, 2009) for c, r in SETTINGS.values()]


def ladder_jobs():
    return [(1996, G.name(c, r, groin=True)) for c, r in LADDER_12.values()]


def advisor_jobs():
    return [(s, G.name(c, r, groin=True)) for s in (1996, 2009) for c, r in ADVISOR.values()]


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "score"
    if mode == "run":
        G.run_all(jobs())
    if mode == "ladder":
        G.run_all(ladder_jobs())
    if mode == "advisor":
        G.run_all(advisor_jobs())
    pairs = {"ladder": lambda: list(zip(ladder_jobs(), LADDER_12)),
             "advisor": lambda: list(zip(advisor_jobs(), list(ADVISOR) * 2))}.get(
        mode, lambda: list(zip(jobs(), list(SETTINGS) * 2)))()
    rows = []
    for (start, variant), label in pairs:
        if G.run_dir(start, variant) is None:
            print(f"missing {start} {variant}")
            continue
        m, tr, d = G.metrics(start, variant), L.trajectory(start, variant), T.cells(start, variant)
        row = dict(window=f"{start}_{U.PERIODS[start][0]}", setting=label, variant=variant,
                   low_dune_p10_rmse=round(tr["p10_rmse"], 3), low_dune_p10_bias=round(tr["p10_bias"], 3),
                   share_below_4p3_rmse=round(tr["below_rmse"], 3), domain_r=round(tr["domain_r_mean"], 3),
                   overall_skill=m["skill"], hits=m["hits"], false_alarms=m["false_alarms"],
                   shoreline_rmse=round(m["rmse_interior_m_yr"], 4))
        for image in T.MAJOR[start]:
            row[f"image_{image}"] = T.major(d, image)
        rows.append(row)
    t = pd.DataFrame(rows)
    out = {"ladder": OUT_LADDER, "advisor": OUT_ADVISOR}.get(mode, OUT)
    out.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(out, index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 30):
        for w, g in t.groupby("window"):
            print(g.dropna(axis=1, how="all").drop(columns=["variant"]).to_string(index=False), "\n")


if __name__ == "__main__":
    main()
