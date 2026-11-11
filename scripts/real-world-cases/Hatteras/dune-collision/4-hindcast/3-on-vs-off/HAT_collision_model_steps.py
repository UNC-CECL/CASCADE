"""
How much each model change improves the fit, calibration (1996-2009) and test (2009-2025): the changes added one at
a time, each run on the one before.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_model_steps.py

The steps (2026-11-04), all on the full-management edgeBE hindcast with the Buxton groin (0.16), the nourishment
fills and storms from the first year:
    base               the dune collision module off, no 198 m road rule, every village 40% of overwash to the
                       shoreface
    + dune growth      calibration only: the growth bounds paired with the module (0.02-0.04, the 1996 default
                       since 2026-10-29) in place of the 0.01-0.03 before it (the base); the test runs keep the
                       2009 rate (0.02-0.04, fitted to the test surveys) throughout
    + collision module Larson et al. (2004) alone, no beach factor (the default since 2026-11-04), Cs 1.9e-5
                       (the offline fit to the lidar crest changes, checked in the model by
                       HAT_collision_cs_sweep.py), off at Cape Point GIS 1-5; one toe per domain, the mean of its
                       screened lidar toes in the surveys nearest the period's start (1996 and 1997; 2004, 2005
                       and 2009, nothing after the start)
    + 198 m road rule  the road manager rebuilds dunes where NC 12 is within 198 m (650 ft) of the shoreline (the
                       FEIS's monitoring warning sign, not an NCDOT dune maintenance rule)
    + Rodanthe 90%     Rodanthe village (GIS 77-83) returns 90% of its overwash to the shoreface (the default since
                       2026-11-04)
Runs: HAT_dune_collision.py cal_step1_base, cal_step2_growth, cal_step3_module_la, cal_step4_road_la,
cal_step6_toeavg_la; HAT_collision_test_2009.py step1_base, step2_module_la, step3_road_la,
nocape_rd198_dry_toeavg_la. The same steps with the domain beach factor (W0 25 m, Cs 1.02e-4; cal_step3_module_dt,
cal_step4_road_dt, cal_step6_toeavg; step2_module_dt, step3_road_dt, nocape_rd198_dry_toeavg) were the figure until
the switch to Larson alone; their figures are in the workspace's beach-factor-steps-until-2026-11-04/. (The
per-cell-toe steps, cal_step3_module ... cal_step6_starttoes and step2_module, step3_road, nocape_rd198_dry,
nocape_rd198_dry_ts, are kept as run on 2026-11-04.)

Writes, in the workspace's dune-collision/4-hindcast/3-on-vs-off/: model_steps_position_<start>_<end>.png (shoreline
position change, each step and CoastSat), model_steps_fit.png (shoreline and dune crest error per step, island and
zones) and supporting/model_steps_fit.csv.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-04 (Larson alone)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
OUT = P.OUT / "2-model-steps"   # this script's figures (topic subfolder)

# --- CONFIG ------------------------------------------------------------------
# One gradient, deep indigo -> teal -> gold (viridis), as HAT_collision_toe_cs_figures.py: the base grey, the dune
# growth a muted lavender leading into it
STEP_STYLE = {"base": ("Base (module off)", "#8C8C8C", 2.4),
              "growth": ("+ dune growth 0.02-0.04", "#9C8FC4", 1.8),
              "module": ("+ collision module", "#3B1F6B", 2.0),
              "road": ("+ 198 m road rule", "#21918C", 2.0),
              "rodanthe": ("+ Rodanthe 90% overwash to the beach", "#E3B219", 2.2)}
# Every module run on one toe per domain with Larson alone (Cs 1.9e-5), the defaults since 2026-11-04
PERIODS = {
    (1996, 2009): ("HAT_dune_collision", "Calibration 1996–2009",
                   (("base", "cal_step1_base"), ("growth", "cal_step2_growth"), ("module", "cal_step3_module_la"),
                    ("road", "cal_step4_road_la"), ("rodanthe", "cal_step6_toeavg_la"))),
    (2009, 2025): ("HAT_collision_test_2009", "Test 2009–2025",
                   (("base", "step1_base"), ("module", "step2_module_la"), ("road", "step3_road_la"),
                    ("rodanthe", "nocape_rd198_dry_toeavg_la"))),
}
ZONES = (("island", "Island", 6, 90), ("buxton_avon_road", "Buxton–Avon road", 8, 20), ("avon", "Avon", 21, 31),
         ("avon_salvo_road", "Avon–Salvo road", 32, 67), ("tri_village", "Tri-Village", 68, 83),
         ("rodanthe_north", "North of Rodanthe", 84, 90))
# -----------------------------------------------------------------------------


# Shoreline: RMSE of the modelled shoreline change rate against CoastSat's (m/yr), island and per zone, as the
# period's own scoring; crest: RMSE of the domain-median crest against the lidar (m) over every survey in the run's
# crest table, GIS 6-90 (the module is off on Cape Point)
def fit(module, steps, start):
    import importlib
    from cascade_pipeline.coastsat_lowess import build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    sys.path.insert(0, str(P._HAT_DIR / "figure_making" / "model_output"))
    from rerender_run_figures import COASTSAT_DATASETS, LOWESS_CONFIG
    from run_figures import TARGET_WINDOW
    m = importlib.import_module(module)
    cs = build_coastsat_series(COASTSAT_DATASETS, active_period_start=start, lowess_config=LOWESS_CONFIG, domains=DOM)
    target = build_target_table(next(c for c in cs if c["active"]), LOWESS_CONFIG, DOM, TARGET_WINDOW)
    target = target.set_index("gis_domain")["target_lrr_m_yr"]
    rows = []
    for step, v in steps:
        rd = m.run_dir(v)
        rate = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain").lrr_m_yr
        err = (rate - target)
        crest = pd.read_csv(rd / TABLES / TABLE)
        for z, _, lo, hi in ZONES:
            e = err.loc[max(lo, 2):min(hi, 89)].dropna()
            c = crest[crest.gis.between(lo, hi)].model_minus_survey_m.dropna()
            rows.append(dict(start=start, step=step, run=v, zone=z, shoreline_rmse_m_yr=np.sqrt((e ** 2).mean()),
                             crest_rmse_m=np.sqrt((c ** 2).mean()), crest_bias_m=c.mean()))
    return pd.DataFrame(rows)


def fit_figure(t):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import record_caption
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharey="row", constrained_layout=True)
    grey = ("#2B2B2B", "#8C8C8C")
    zone_colors = dict(zip([z for z, *_ in ZONES[1:]], ("#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#E31A1C")))   # no brown
    for j, ((start, end), (_, title, steps)) in enumerate(PERIODS.items()):
        d = t[t.start == start]
        names = [s for s, _ in steps]
        x = np.arange(len(names))
        for i, (col, ylab) in enumerate((("shoreline_rmse_m_yr", "Shoreline change rate RMSE (m/yr)"),
                                         ("crest_rmse_m", "Dune crest RMSE against lidar (m)"))):
            ax = axes[i, j]
            for z, label, *_ in ZONES:
                y = d[d.zone == z].set_index("step")[col].reindex(names).to_numpy()
                if z == "island":
                    ax.plot(x, y, "-o", color=grey[0], lw=2.6, ms=7, zorder=5, label="Island (GIS 6-90)")
                    for xi, yi in zip(x, y):
                        ax.annotate(f"{yi:.2f}", (xi, yi), textcoords="offset points", xytext=(0, 8),
                                    ha="center", fontsize=8.5, color=grey[0])
                else:
                    ax.plot(x, y, "-o", color=zone_colors[z], lw=1.2, ms=4, alpha=0.8, label=label)
            ax.set_xticks(x, [STEP_STYLE[s][0].replace(" overwash to the beach", "").replace(" (module off)", "")
                              for s in names], rotation=30, ha="right", fontsize=9)
            ax.grid(axis="y", color="#E0E0E0", lw=0.8)
            ax.spines[["top", "right"]].set_visible(False)
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(title, fontsize=12, fontweight="bold")
    axes[0, 1].legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=9)
    fig.suptitle("Model fit as each change is added (lower is better)", fontsize=13, fontweight="bold")
    out = OUT / "model_steps_fit.png"
    fig.savefig(out, dpi=300, facecolor="white")
    plt.close(fig)
    record_caption(out, (
        "Model fit as each change is added, each step on the one before; lower is better. Top: RMSE of the modelled "
        "shoreline change rate against CoastSat's (7-domain LOWESS), per zone and over the island. Bottom: RMSE of "
        "the modelled domain-median dune crest against the lidar's, over every lidar survey in the period, GIS 6-90 "
        "(Cape Point, where the module is off, left out). Steps: the base (dune collision module off, no 198 m road "
        "rule, villages 40% of overwash to the shoreface); calibration only, the dune growth bounds paired with the "
        "module (0.02-0.04); the module (Larson et al. 2004 alone, no beach factor, Cs 1.9e-5, one toe per domain from the lidar nearest the "
        "period's start, 1996-1997 or 2004-2009); the 198 m road-dune rule; Rodanthe village (GIS 77-83) 90% of "
        "overwash to the shoreface."))
    return out


def main():
    tables = []
    for (start, end), (module, title, steps) in PERIODS.items():
        runs = tuple((v, *STEP_STYLE[s]) for s, v in steps)
        out, gis, ch = P.figure(
            start, end, module, runs,
            f"Cumulative effect of the model changes, {start}–{end}",
            f"2-model-steps/model_steps_position_{start}_{end}.png",
            "each line adding one change to the one before: the base (dune collision module off, no 198 m road rule, "
            "villages 40% of overwash to the shoreface)"
            + ("; the dune growth bounds paired with the module (0.02-0.04)" if start == 1996 else "")
            + "; the module (Larson et al. 2004 alone, no beach factor, Cs 1.9e-5, off at Cape Point GIS 1-5, one toe per domain: the mean "
            "of its lidar toes in the surveys nearest the period's start, 1996-1997 or 2004-2009); the 198 m "
            "road-dune rule; Rodanthe village (GIS 77-83) 90% of overwash to the shoreface.")
        print(f"wrote {out}")
        tables.append(fit(module, steps, start))
    t = pd.concat(tables, ignore_index=True)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(4).to_csv(OUT / "supporting" / "model_steps_fit.csv", index=False)
    print(f"wrote {fit_figure(t)}")
    pd.set_option("display.width", 200)
    for col in ("shoreline_rmse_m_yr", "crest_rmse_m"):
        print(f"\n{col}")
        print(t.pivot_table(index=["start", "zone"], columns="step", values=col, sort=False).round(3).to_string())


if __name__ == "__main__":
    main()
