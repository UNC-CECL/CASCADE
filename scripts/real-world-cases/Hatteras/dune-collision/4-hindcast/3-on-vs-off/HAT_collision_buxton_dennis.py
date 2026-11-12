"""
What the collision module does at Buxton (GIS 11-14), where the base model missed the dune lowering after Bonnie
(Aug 1998) and Dennis (Sep 1999): domain-median dune crests through 1996-2009, lidar against the model with the
module off, with Larson alone (the default, Cs 1.9e-5) and with Larson x the beach factor (Cs 1.02e-4, the default
until 2026-11-04).

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_buxton_dennis.py

Runs (HAT_dune_collision.py, nothing rerun): cal_step2_growth (module off), cal_step3_module_la (Larson alone),
cal_step3_module_dt (beach factor); all with growth 0.02-0.04, one toe per domain, no road rule, Rodanthe 40%, so
the module is the only difference. Crests: each run's dune_crest_vs_lidar.csv (domain medians at each survey).

Writes, in the workspace's dune-collision/4-hindcast/7-storm-examples/ (D:): collision_buxton_dennis.png, the PDF and
supporting/collision_buxton_dennis.csv.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
OUT = P.OUT / "7-storm-examples"   # this script's figures (topic subfolder)

sys.path.insert(0, str(P._HAT_DIR / "dune-collision" / "4-hindcast" / "1-calibration-1996-2009"))

# --- CONFIG ------------------------------------------------------------------
GIS = (11, 12, 13, 14)
RUNS = (("cal_step2_growth", "Module off", "#7F7F7F", "-"),
        ("cal_step3_module_la", "Collision module (Larson, Cs 1.9e-5)", "#6A51A3", "-"))
STORMS = {"Bonnie": 1998.65, "Dennis": 1999.67, "Isabel": 2003.71}
# -----------------------------------------------------------------------------


def decimal_year(state):
    return pd.to_datetime(state).dt.year + (pd.to_datetime(state).dt.dayofyear - 1) / 365.25


# Isabel (Sep 2003): top row, the four domains the module lowered most across it; bottom row, the four the lidar
# shows dropping most (HAT_collision_isabel.py, 2001 to 2004)
ISABEL_ROWS = (((23, 31, 58, 59), "Module lowered most at Isabel"), ((7, 9, 51, 52), "Lidar dropped most at Isabel"))


def main(rows_gis=((GIS, None),), out_name="collision_buxton_dennis.png", caption_lead=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import HAT_dune_collision as C
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    from site_layer import hat_figure_style as S

    S.apply_style()
    every = [g for gs, _ in rows_gis for g in gs]
    rows = []
    for v, name, _, _ in RUNS:
        t = pd.read_csv(C.run_dir(v) / TABLES / TABLE)
        t = t[t.gis.isin(every)].assign(run=v, label=name)
        rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    # the survey's own date (first month flown) for the lidar; the model state the survey is read against
    t["survey_year"] = t.flown.str.slice(0, 7).pipe(lambda s: pd.to_datetime(s, format="%Y-%m")).pipe(
        lambda d: d.dt.year + (d.dt.month - 0.5) / 12)
    nr = len(rows_gis)
    fig, ax_grid = plt.subplots(nr, 4, figsize=(S.FIG_W_DOUBLE, 2.5 * nr + 0.3), sharey=True, sharex=True,
                                layout="constrained", squeeze=False)
    axes = ax_grid.ravel()
    for (gs, title), row in zip(rows_gis, ax_grid):
        if title:
            row[0].annotate(title, xy=(0, 1.13), xycoords="axes fraction", ha="left", va="bottom", fontsize=9,
                            fontweight="bold")
    for i, (ax, g) in enumerate(zip(axes, every)):
        d = t[t.gis == g]
        lid = d.drop_duplicates("survey").sort_values("survey_year")
        for _, x in STORMS.items():
            ax.axvline(x, color=S.INK_MUTED, lw=0.6, ls=":")
        for v, name, c, ls in RUNS:
            r = d[d.run == v].sort_values("survey_year")
            ax.plot(r.survey_year, r.model_median_m, color=c, ls=ls, lw=1.6, label=name)
        ax.plot(lid.survey_year, lid.survey_median_m, "o", color="#2166AC", ms=4, label="Lidar", zorder=5)
        ax.set_title(f"GIS {g}", loc="center", fontsize=9)
        ax.text(0.03, 0.97, f"({chr(97 + i)})", transform=ax.transAxes, ha="left", va="top", fontweight="bold")
        ax.set_xlim(1996, 2010)
        ax.set_xticks([1998, 2002, 2006, 2010])
        ax.grid(axis="y", color=S.GRID_C, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    for name, x in STORMS.items():
        axes[0].text(x - 0.15, 2.15, name, rotation=90, ha="right", va="bottom", fontsize=6.5, color=S.INK_MUTED)
    for row in ax_grid:
        row[0].set_ylabel("Dune crest, domain median\n(m NAVD88)")
    axes[0].set_ylim(2, 8.5)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="outside lower center", ncol=4, frameon=False, fontsize=7.5)
    out = OUT / out_name
    S.save(fig, out, bbox_inches="tight", close=True)
    t.to_csv(OUT / "supporting" / out_name.replace(".png", ".csv"), index=False)
    S.record_caption(out, (
        (caption_lead or "Dune crests at Buxton (GIS 11-14), 1996-2009")
        + ": lidar (dots, domain median) against the model with the dune "
        "collision module off (grey) and on (purple; Larson et al. 2004 per storm, Cs 1.9e-5). Both runs share every "
        "other setting (dune growth 0.02-0.04, one toe per domain, no road rule). Dotted lines: Bonnie (Aug 1998), "
        "Dennis (Sep 1999), Isabel (Sep 2003)."))
    print(f"wrote {out}")
    p = t[t.survey.astype(str).isin(["1997", "1998", "1999", "2000"])].pivot_table(
        index=["gis", "survey"], columns="run", values="model_median_m")
    p["lidar"] = t.groupby(["gis", "survey"]).survey_median_m.first()
    print(p.round(2).to_string())


if __name__ == "__main__":
    import argparse
    _p = argparse.ArgumentParser()
    _p.add_argument("--isabel", action="store_true",
                    help="the Isabel version: the domains the module and the lidar lowered most across it")
    if _p.parse_args().isabel:
        main(ISABEL_ROWS, "collision_isabel_domains.png",
             "Dune crests, 1996-2009, at the four domains the collision module lowered most across Hurricane Isabel "
             "(top: GIS 23, 31, 58, 59) and the four the lidar shows dropping most from 2001 to 2004 (bottom: GIS 7, "
             "9, 51, 52)")
    else:
        main()
