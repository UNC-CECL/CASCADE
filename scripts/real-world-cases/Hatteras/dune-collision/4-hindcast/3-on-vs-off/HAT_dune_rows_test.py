"""
A back dune row at the window median, with the collision module eroding the front row (2026-11-11, test): crests
through time on GIS 9-15 against the lidar, and the fit of each run. The position-change figure is
HAT_collision_on_off_position.py maxmed (same runs, MAXMED_RUNS there).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_dune_rows_test.py [lateral]

"lateral" draws the same for the median back row with lateral erosion (MAXMED_LAT_RUNS there; 2026-11-11), into
dune_rows_maxmed_lateral_crests_over_time_GIS9-15.png and supporting/dune_rows_maxmed_lateral_metrics.csv.

Runs (HAT_dune_collision.py, 1996-2009, module-only setup): cal_step2_growth_noscrape and _maxmed (module off, both
rows at the window max / back row at the window median), cal_toeberm_cs2 (both rows at max, both eroded in proportion,
the current setting), cal_toeberm_cs2_front (both at max, front row eroded) and cal_toeberm_cs2_maxmed (back row at
median, front row eroded).

Writes, in the workspace's dune-collision/4-hindcast/3-on-vs-off/6-other-tests/dune-rows-maxmed/ (D:):
- dune_rows_maxmed_crests_over_time_GIS9-15.png: the domain-median crest (the higher row) through time, and the
  front and back rows of the max/median run
- supporting/dune_rows_maxmed_metrics.csv: crest RMSE, 2009 offset, overwash, shoreline skill per run

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-11
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
import HAT_collision_row_mode_sweep as R  # noqa: E402

H = R.H

# --- CONFIG ------------------------------------------------------------------
OUT = P.OUT / "6-other-tests" / "dune-rows-maxmed"
TS_GIS = range(9, 16)
EVENTS = ((1998.65, "Bonnie"), (1999.67, "Dennis"), (2003.71, "Isabel"))
# -----------------------------------------------------------------------------


def main():
    lateral = len(sys.argv) > 1 and sys.argv[1] == "lateral"
    runs = P.MAXMED_LAT_RUNS if lateral else P.MAXMED_RUNS
    tag = "dune_rows_maxmed_lateral" if lateral else "dune_rows_maxmed"
    back_of = "cal_toeberm_lat_cs10_maxmed" if lateral else "cal_toeberm_cs2_maxmed"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer import hat_figure_style as S
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    S.apply_style()
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    m = pd.DataFrame([R.metrics(r[0], 0.0) for r in runs]).drop(columns="cs")
    m.round(4).to_csv(OUT / "supporting" / f"{tag}_metrics.csv", index=False)
    pd.set_option("display.width", 220)
    print(m.round(3).to_string(index=False))

    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(TS_GIS)].copy()
    t = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = t.dt.year + (t.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), yr=("yr", "first")).reset_index()
    fig, axs = plt.subplots(len(TS_GIS), 1, figsize=(S.FIG_W_DOUBLE * 0.75, 10.6), sharex=True, sharey=True,
                            layout="constrained")
    for v, lab, colour, lw, *ls in runs:
        rd = H.run_dir(v)
        c = load_saved_cascade(rd / f"{rd.name}.npz")
        d = DOM.for_total(len(c.barrier3d))
        for ax, g in zip(axs, TS_GIS):
            b = c.barrier3d[d.gis_to_pad(g)]
            rows = (np.asarray(b.DuneDomain, float) + b._BermEl) * 10 + H.MHW_NAVD88_M     # (years, cells, rows)
            yrs = H.RUN_START + 1 + np.arange(rows.shape[0])
            ax.plot(yrs, np.median(rows.max(axis=2), axis=1), color=colour, lw=lw - 0.3, ls=ls[0] if ls else "-",
                    label=lab)
            if v == back_of:
                ax.plot(yrs, np.median(rows[:, :, 1], axis=1), color=colour, lw=0.9, ls=":",
                        label="  its back row")
    for ax, g in zip(axs, TS_GIS):
        L = lid[lid.gis == g].sort_values("yr")
        ax.plot(L.yr, L.z, "o", color="#1A1A1A", ms=4.5, label="Lidar (median)", zorder=5)
        ax.axhline(H.TOE_AT_BERM_M, color="0.55", lw=0.8, label="Berm, 1.7 m (the toe)")
        for x, _ in EVENTS:
            ax.axvline(x, color="0.7", lw=0.6, ls=":")
        ax.set_title(f"GIS {g}", loc="left", fontsize=9, fontweight="bold")
        ax.set_ylabel("Crest (m)", fontsize=8)
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    for x, n in EVENTS:
        axs[0].text(x, 1.03, n, transform=axs[0].get_xaxis_transform(), ha="center", fontsize=7, color="0.4")
    axs[-1].set_xlabel("Year")
    axs[-1].set_xlim(1996, 2010)
    axs[0].set_ylim(0, 9)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=2, fontsize=7.5, frameon=False)
    out = OUT / f"{tag}_crests_over_time_GIS9-15.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Domain-median dune crest (m NAVD88; the higher of Barrier3D's two dune rows) through 1996-2009, GIS 9-15, "
        "against the lidar domain median at each survey. Grey: module off; purple: both rows at the window maximum, "
        "proportional erosion, Cs 2e-6 (current); blue: both rows at the maximum, lateral erosion, Cs 1e-5; oranges "
        "to red (dashed): the back row at the window median, lateral erosion, Cs 3e-6, 1e-5, 2e-5, the dotted line "
        "the back row of the 1e-5 run. Toe at the berm. Model state at 1 January each year; dotted verticals: "
        "Bonnie, Dennis, Isabel.") if lateral else (
        "Domain-median dune crest (m NAVD88; the higher of Barrier3D's two dune rows) through 1996-2009, GIS 9-15, "
        "against the lidar domain median at each survey: both rows at the highest cell of each profile's hand-picked "
        "window (v3) or the landward row at the median of that window (dashed); the dune collision module off "
        "(grey) or on (toe at the berm, Cs 2e-6), eroding both rows in proportion (purple, current) or the front "
        "row only (blue, orange). Orange dotted: the back row of the max/median run. Model state at 1 January each "
        "year; dotted verticals: Bonnie, Dennis, Isabel."))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
