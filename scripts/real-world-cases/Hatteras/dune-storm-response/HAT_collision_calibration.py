"""
Calibrating the model with the collision module on the 2026-10-29 defaults: the dune growth rate (island-wide and
by zone), the wave height and the groin strength.

    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_collision_calibration.py growth
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_collision_calibration.py runs [run ...]
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_collision_calibration.py groin [run ...]
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_collision_calibration.py hs
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_collision_calibration.py gsweep
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_collision_calibration.py buxton

The runs are HAT_dune_collision.py's cal_* variants (managed edgeBE 1996-2009, groin, StormStart 1, village dune cap
8.0 m NAVD88; the module at Cs 1.24e-4, W0 45 m), filed under its runs/. This script only scores them.

growth  every cal_r* run against the lidar crests, zone by zone: per survey, the zone's median modelled crest
        (same cells as the lidar) less the lidar's, and its RMSE over the surveys after the start; the best
        rate per zone and island-wide. Zones follow the management and the beach: the cape (Cape Point to
        the groin), the Buxton-Avon road, Avon, the Avon-Salvo road, the Tri-Village and the Rodanthe hotspot.
        Writes tables/growth_by_zone.csv, tables/growth_best.csv and figures/growth_by_zone.png
runs    shoreline (CoastSat LRR) and dune scores of the named runs, by zone and island-wide.
        Writes tables/run_scores.csv
groin   the groin's own fit for the named runs: the annual GIS 5 - 6 gap against CoastSat and the net change.
        Writes tables/groin_scores.csv
hs      the wave-height sweep (HS_RUNS) from run_scores.csv and groin_scores.csv (score both first): shoreline,
        crest and groin scores against Hs, with the end rates each Hs was solved to.
        Writes tables/hs_sweep.csv and figures/hs_sweep.png
gsweep  the groin strength sweep (GROIN_RUNS) the same way: the groin's fit and the shoreline by zone against
        the strength. Writes tables/groin_sweep.csv and figures/groin_sweep.png
buxton  is the CoastSat erosion at GIS 10-15 (Buxton-Avon) steady or in steps at the storms? Per window
        (BUXTON_WINDOWS: 1996-2009, and 2009-2025 with the 2017-18 and 2022 fills), per transect the annual
        median chainage (seaward +) less its median over the window's base years, then the median per GIS;
        the change across each storm or fill, the trends between them, and for 1996-2009 the model's own
        shoreline (cal_control, cal_r03). Writes tables/buxton_coastsat_annual_<window>.csv,
        tables/buxton_steps_<window>.csv, tables/buxton_trends_<window>.csv and figures/buxton_coastsat_<window>.png

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-29
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_collision as H  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "2026-10-29-collision-calibration"
OUT = H.PROJECT_ROOT / "output" / "raw_runs" / "dune-storm-response" / TAG
# (name, first GIS, last GIS)
ZONES = (("cape", 1, 7), ("buxton_avon_road", 8, 20), ("avon", 21, 31), ("avon_salvo_road", 32, 67),
         ("tri_village", 68, 83), ("rodanthe_north", 84, 90))
GROWTH_RUNS = {"cal_r02": 0.02, "cal_r03": 0.03, "cal_r04": 0.04, "cal_r05": 0.05, "cal_r06": 0.06}
# The wave-height sweep: growth 0.02-0.04, each on the end rates solved at its own Hs (cal_r03 is Hs 3.0)
HS_RUNS = {"cal_hs2": 2.0, "cal_hs2p5": 2.5, "cal_r03": 3.0, "cal_hs3p5": 3.5, "cal_hs4": 4.0}
# The groin strength sweep at Hs 3.0 (cal_r03 is the 0.10 default)
GROIN_RUNS = {"cal_r03": 0.10, "cal_g13": 0.13, "cal_g16": 0.16, "cal_g19": 0.19, "cal_g22": 0.22}
# The Buxton-Avon erosion check, per window: the years shown, each transect's zero (its median over base), the
# trend extended (trend), the events marked (decimal years; fills as start, end), the changes across events
# (label, years before, years after) and the trends between them; the model is drawn for 1996-2009 only
BUXTON_GIS = range(10, 16)
BUXTON_MIN_OBS = 3                    # images a transect needs in a year for that year's median
BUXTON_WINDOWS = {
    "1996_2009": dict(
        years=(1988, 2012), base=(1994, 1997), trend=(1988, 1997),
        storms={"Bonnie": 1998.65, "Dennis / Floyd": 1999.70, "Isabel": 2003.71},
        fills={},
        steps=[("Bonnie, Dennis / Floyd", (1996, 1997), (1999, 2000)), ("Isabel", (2002, 2003), (2005, 2006))],
        trends=[(1988, 1997), (2000, 2003), (2006, 2010)],
        model=("cal_control", "cal_r03")),
    "2009_2025": dict(
        years=(2005, 2025), base=(2008, 2010), trend=(2005, 2010),
        storms={"Irene": 2011.65, "Sandy": 2012.83, "Matthew": 2016.77, "Dorian": 2019.68, "Erin": 2025.63},
        # Buxton beach nourishment 2017-06-21 to 2018-02-27; Buxton shore protection 2022-06-30 to 2022-08-16 (GIS 6-16)
        fills={"2017-18 fill": (2017.47, 2018.16), "2022 fill": (2022.49, 2022.62)},
        steps=[("Irene, Sandy", (2009, 2010), (2012, 2013)), ("2017-18 fill", (2015, 2016), (2018, 2018)),
               ("Dorian", (2018, 2018), (2020, 2020)), ("2022 fill", (2021, 2021), (2023, 2023)),
               ("Erin (2025 partial)", (2024, 2024), (2025, 2025))],
        trends=[(2005, 2010), (2013, 2016), (2018, 2021), (2023, 2025)],
        model=()),
}
# -----------------------------------------------------------------------------


def zone_of(g):
    return next(z for z, lo, hi in ZONES if lo <= g <= hi)


# Per survey and GIS, the median modelled and lidar crest over the lidar's cells (m NAVD88)
def crests(variant):
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rd = H.run_dir(variant)
    if rd is None:
        raise SystemExit(f"no run for {variant}")
    c = load_saved_cascade(rd / f"{rd.name}.npz")
    dom = DOM.for_total(len(c.barrier3d))
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[sv.status == "ok"].copy()
    sv["row"] = pd.to_datetime(sv.model_state).dt.year - (H.RUN_START + 1)
    model = {}
    for g in range(1, 91):
        b = c.barrier3d[dom.gis_to_pad(g)]
        model[g] = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + H.MHW_NAVD88_M
    sv = sv[(sv.row >= 0) & (sv.row < len(model[1]))]
    sv["model"] = [model[g][r][cell] for g, r, cell in zip(sv.gis, sv.row, sv.cell)]
    out = sv.groupby(["survey", "gis"]).agg(model=("model", "median"), lidar=("crest_m_navd88", "median")).reset_index()
    out["zone"] = out.gis.map(zone_of)
    return out


# Zone (and island) median crest per survey, and the RMSE of model - lidar over the surveys after the start
def zone_scores(cr):
    first = cr.survey.min()
    rows = []
    for zone, x in [*cr.groupby("zone"), ("island", cr)]:
        per = x.groupby("survey")[["model", "lidar"]].median()
        later = per.loc[per.index != first]
        d = later.model - later.lidar
        last = per.index.max()
        rows.append(dict(zone=zone, crest_rmse=float(np.sqrt((d ** 2).mean())), crest_bias=float(d.mean()),
                         end_minus_lidar=float(per.loc[last, "model"] - per.loc[last, "lidar"]),
                         r_domains=float(np.corrcoef(*x[x.survey != first][["model", "lidar"]].to_numpy().T)[0, 1])))
    return pd.DataFrame(rows)


def cmd_growth(_):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import apply_style, record_caption
    apply_style()
    tabs = []
    for v, r in GROWTH_RUNS.items():
        if H.run_dir(v) is None:
            print(f"missing {v}")
            continue
        tabs.append(zone_scores(crests(v)).assign(run=v, growth_mean=r))
    t = pd.concat(tabs, ignore_index=True)
    best = t.loc[t.groupby("zone").crest_rmse.idxmin(), ["zone", "run", "growth_mean", "crest_rmse", "crest_bias",
                                                          "end_minus_lidar"]]
    order = [z for z, *_ in ZONES] + ["island"]
    best = best.set_index("zone").loc[order].reset_index()
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    t.round(4).to_csv(OUT / "tables" / "growth_by_zone.csv", index=False)
    best.round(4).to_csv(OUT / "tables" / "growth_best.csv", index=False)
    pd.set_option("display.width", 200)
    print(t.pivot(index="zone", columns="growth_mean", values="crest_rmse").loc[order].round(3))
    print(t.pivot(index="zone", columns="growth_mean", values="crest_bias").loc[order].round(2))
    print(best.round(3).to_string(index=False))

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    shades = ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#333333"]
    for z, col in zip(order, shades):
        x = t[t.zone == z].sort_values("growth_mean")
        lw = 2.4 if z == "island" else 1.6
        axs[0].plot(x.growth_mean, x.crest_rmse, "o-", color=col, lw=lw, label=z.replace("_", " "))
        axs[1].plot(x.growth_mean, x.crest_bias, "o-", color=col, lw=lw)
    axs[1].axhline(0, color="0.6", lw=0.8)
    axs[0].set_ylabel("Zone median crest, RMSE against lidar (m)")
    axs[1].set_ylabel("Mean model − lidar (m)")
    for ax in axs:
        ax.set_xlabel("Dune growth rate, mean of the range (± 0.01)")
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    out = OUT / "figures" / "growth_by_zone.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Dune growth rate against the lidar crests by zone, with the collision module (Cs 1.24e-4, W0 45 m) on the "
        "2026-10-29 defaults (village dune cap 8.0 m NAVD88, Hs 3.0): per survey after 1996, the zone's median "
        "modelled crest over the lidar's cells less the lidar's. Left, RMSE over the surveys; right, the mean "
        f"difference. Experiment {TAG}."))
    print(f"wrote {out}")


def cmd_runs(a):
    import json
    from cascade_pipeline.coastsat_lowess import build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    sys.path.insert(0, str(H._HAT_DIR / "figure_making" / "model_output"))
    from rerender_run_figures import COASTSAT_DATASETS, LOWESS_CONFIG
    from run_figures import TARGET_WINDOW
    cs = build_coastsat_series(COASTSAT_DATASETS, active_period_start=H.RUN_START, lowess_config=LOWESS_CONFIG,
                               domains=DOM)
    target = build_target_table(next(c for c in cs if c["active"]), LOWESS_CONFIG, DOM, TARGET_WINDOW)
    target = target.set_index("gis_domain")["target_lrr_m_yr"]
    rows = []
    for v in a.runs:
        rd = H.run_dir(v)
        rate = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain").lrr_m_yr
        meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
        err = (rate - target).loc[2:89]
        row = dict(run=v, rmse=float(meta["skill"]["rmse_interior_m_yr"]),
                   bias=float(meta["skill"]["mean_bias_interior_m_yr"]))
        for z, lo, hi in ZONES:
            e = err.loc[max(lo, 2):min(hi, 89)]
            row[f"lrr_rmse_{z}"] = float(np.sqrt((e ** 2).mean()))
        zs = zone_scores(crests(v)).set_index("zone")
        for z in zs.index:
            row[f"crest_rmse_{z}"] = zs.loc[z, "crest_rmse"]
        rows.append(row)
    t = pd.DataFrame(rows)
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    path = OUT / "tables" / "run_scores.csv"
    if path.exists():
        old = pd.read_csv(path)
        t = pd.concat([old[~old.run.isin(t.run)], t], ignore_index=True)
    t.round(4).to_csv(path, index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(t.set_index("run").round(3).T.to_string())


# The groin's own fit (hard-structures/groin/.../2026-10-20-spread-relaxation): the annual GIS 5 - GIS 6 gap against
# CoastSat, and the end-of-window net change against the CoastSat target over GIS 2-89, 2-12 and 4-9 (m)
def cmd_groin(a):
    sys.path.insert(0, str(H._HAT_DIR / "hard-structures" / "groin" / "4-hindcast" / "3-wave-driven-1967-2025"
                           / "2026-10-20-spread-relaxation"))
    import spread_relaxation as sr
    tf = sr.tf
    import corrected_groin as cgm
    import site_layer.hat_scoring as common
    from matrix_vs_observed import model_end_change
    from site_layer.hat_observed_rates import net_change_domain_csv
    # schedule_refit puts the workspace twin of the CoastSat gap folder (D:) on the path; the script is in the repo
    sys.path.insert(0, str(H._HAT_DIR / "hard-structures" / "groin" / "1-observations" / "analysis" / "2-coastsat_gap"))
    period = H.RUN_START
    gap_obs = cgm.observed_change(period)
    net_obs = pd.read_csv(net_change_domain_csv(*tf.NET_KEY[period])).set_index("gis_domain")["net_change_lowess7_m"]
    rows = []
    for v in a.runs:
        rd = H.run_dir(v)
        raw = model_end_change(rd)
        net = common.smooth_like_target(raw)
        path = next(rd.glob("*_shoreline_matrix.npy"))
        # schedule_refit's GIS 5 / 6 columns assume 15 buffer domains a side; today's runs have 30
        pad = (np.load(path).shape[1] - 90) // 2
        cgm.sr.UP, cgm.sr.DOWN = pad + 6 - 1, pad + 5 - 1
        sc = cgm.score_annual(path, period, gap_obs)
        row = dict(run=v, gap_rmse=sc["rmse_annual"], gap_bias=sc["bias_annual"])
        for name, (lo, hi) in (("net_rmse_2_89", (2, 89)), ("net_rmse_2_12", (2, 12)), ("net_rmse_4_9", (4, 9))):
            r = (net - net_obs).loc[lo:hi].dropna()
            row[name] = float(np.sqrt(np.mean(r ** 2)))
        row.update({f"gis{g}": float(raw.loc[g]) for g in range(3, 12)})
        rows.append(row)
    t = pd.DataFrame(rows)
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    path = OUT / "tables" / "groin_scores.csv"
    if path.exists():
        old = pd.read_csv(path)
        t = pd.concat([old[~old.run.isin(t.run)], t], ignore_index=True)
    t.round(3).to_csv(path, index=False)
    pd.set_option("display.width", 250)
    print(t.round(1).to_string(index=False))


def cmd_hs(_):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import apply_style, record_caption
    apply_style()
    runs = pd.read_csv(OUT / "tables" / "run_scores.csv").set_index("run")
    groin = pd.read_csv(OUT / "tables" / "groin_scores.csv").set_index("run")
    have = [v for v in HS_RUNS if v in runs.index and v in groin.index]
    t = runs.loc[have].join(groin.loc[have, ["gap_rmse", "gap_bias", "net_rmse_2_89", "net_rmse_4_9"]])
    t.insert(0, "hs_m", [HS_RUNS[v] for v in have])
    t["end_gis1"] = [H.ENDS_BY_HS[HS_RUNS[v]][1] for v in have]
    t["end_gis90"] = [H.ENDS_BY_HS[HS_RUNS[v]][90] for v in have]
    t.round(4).to_csv(OUT / "tables" / "hs_sweep.csv")
    pd.set_option("display.width", 250)
    print(t.set_index("hs_m").round(3).T.to_string())

    order = [z for z, *_ in ZONES]
    shades = ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"]
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.6))
    ax = axs[0, 0]
    for z, col in zip(order, shades):
        ax.plot(t.hs_m, t[f"lrr_rmse_{z}"], "o-", color=col, lw=1.3, label=z.replace("_", " "))
    ax.plot(t.hs_m, t.rmse, "o-", color="#333333", lw=2.6, label="island (GIS 2-89)")
    ax.set_ylabel("Shoreline LRR RMSE against CoastSat (m/yr)")
    ax.legend(fontsize=7, ncol=2)
    ax = axs[0, 1]
    for z, col in zip(order, shades):
        ax.plot(t.hs_m, t[f"crest_rmse_{z}"], "o-", color=col, lw=1.3)
    ax.plot(t.hs_m, t.crest_rmse_island, "o-", color="#333333", lw=2.6)
    ax.set_ylabel("Median dune crest RMSE against lidar (m)")
    ax = axs[1, 0]
    ax.plot(t.hs_m, t.gap_rmse, "o-", color="#333333", lw=2.2, label="annual GIS 5 - 6 gap")
    ax.plot(t.hs_m, t.net_rmse_4_9, "s--", color="#D95F02", lw=1.6, label="net change, GIS 4-9")
    ax.plot(t.hs_m, t.net_rmse_2_89, "^:", color="#7570B3", lw=1.6, label="net change, GIS 2-89")
    ax.set_ylabel("Groin fit, RMSE against CoastSat (m)")
    ax.legend(fontsize=7)
    ax = axs[1, 1]
    ax.plot(t.hs_m, t.end_gis1, "o-", color="#1B9E77", lw=2, label="GIS 1 (Cape Point)")
    ax.plot(t.hs_m, t.end_gis90, "o-", color="#E6AB02", lw=2, label="GIS 90 (Rodanthe)")
    ax.set_ylabel("End rate solved to the net change (m/yr)")
    ax.legend(fontsize=7)
    for ax, k in zip(axs.flat, "abcd"):
        ax.set_xlabel("Wave height Hs (m)")
        ax.set_xticks(t.hs_m)
        ax.set_title(f"({k})", loc="left", fontsize=9)
    fig.tight_layout()
    out = OUT / "figures" / "hs_sweep.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Wave height sweep with the collision module (Cs 1.24e-4, W0 45 m, growth 0.02-0.04) on the 2026-10-29 "
        "defaults (village dune cap 8.0 m NAVD88, groin 0.10), each Hs on the GIS 1 / 90 end rates solved at that "
        "Hs (d). (a) Shoreline LRR RMSE against CoastSat by zone and island-wide; (b) zone median crest RMSE "
        "against the lidar; (c) the groin's fit: the annual GIS 5 - 6 gap and the end-of-window net change. "
        f"Experiment {TAG}."))
    print(f"wrote {out}")


def cmd_gsweep(_):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import apply_style, record_caption
    apply_style()
    runs = pd.read_csv(OUT / "tables" / "run_scores.csv").set_index("run")
    groin = pd.read_csv(OUT / "tables" / "groin_scores.csv").set_index("run")
    have = [v for v in GROIN_RUNS if v in runs.index and v in groin.index]
    t = runs.loc[have].join(groin.loc[have].drop(columns=[c for c in groin.columns if c in runs.columns]))
    t.insert(0, "strength", [GROIN_RUNS[v] for v in have])
    t.round(4).to_csv(OUT / "tables" / "groin_sweep.csv")
    pd.set_option("display.width", 250)
    print(t.set_index("strength").round(3).T.to_string())

    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
    ax = axs[0]
    ax.plot(t.strength, t.gap_rmse, "o-", color="#333333", lw=2.2, label="RMSE")
    ax.plot(t.strength, t.gap_bias, "s--", color="#D95F02", lw=1.6, label="bias (model - CoastSat)")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_ylabel("Annual GIS 5 - 6 gap against CoastSat (m)")
    ax.legend(fontsize=7)
    ax = axs[1]
    ax.plot(t.strength, t.net_rmse_4_9, "o-", color="#D95F02", lw=1.8, label="GIS 4-9")
    ax.plot(t.strength, t.net_rmse_2_12, "s--", color="#1B9E77", lw=1.6, label="GIS 2-12")
    ax.plot(t.strength, t.net_rmse_2_89, "^:", color="#7570B3", lw=1.6, label="GIS 2-89")
    ax.set_ylabel("Net change 1996-2009, RMSE against CoastSat (m)")
    ax.legend(fontsize=7)
    ax = axs[2]
    shades = ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"]
    for (z, *_), col in zip(ZONES, shades):
        ax.plot(t.strength, t[f"lrr_rmse_{z}"], "o-", color=col, lw=1.3, label=z.replace("_", " "))
    ax.plot(t.strength, t.rmse, "o-", color="#333333", lw=2.6, label="island (GIS 2-89)")
    ax.set_ylabel("Shoreline LRR RMSE against CoastSat (m/yr)")
    ax.legend(fontsize=7, ncol=2)
    for ax, k in zip(axs, "abc"):
        ax.set_xlabel("Groin strength")
        ax.set_xticks(t.strength)
        ax.set_title(f"({k})", loc="left", fontsize=9)
    fig.tight_layout()
    out = OUT / "figures" / "groin_sweep.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Groin strength sweep with the collision module (Cs 1.24e-4, W0 45 m, growth 0.02-0.04) at Hs 3.0 on its "
        "own ends, on the 2026-10-29 defaults (village dune cap 8.0 m NAVD88). (a) The annual GIS 5 - 6 shoreline "
        "gap across the groin against CoastSat; (b) the end-of-window net change; (c) shoreline LRR RMSE by zone. "
        f"Experiment {TAG}."))
    print(f"wrote {out}")


def cmd_buxton(_):
    import glob
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import apply_style, record_caption
    from site_layer.hat_observed_rates import COASTSAT_TIMESERIES, transect_lookup
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    apply_style()

    # Every transect's annual median chainage (m, seaward +), years with BUXTON_MIN_OBS images or more
    lookup = pd.read_csv(transect_lookup()).dropna(subset=["domain_number"])
    lookup = lookup[lookup.domain_number.astype(int).isin(BUXTON_GIS)]
    annual = {}
    for tid, g in zip(lookup.transect_id, lookup.domain_number.astype(int)):
        hits = glob.glob(str(COASTSAT_TIMESERIES / "*" / f"{tid}.csv"))
        if not hits:
            continue
        ts = pd.read_csv(hits[0])
        ts.columns = ["date", "chainage_m"]
        ts["year"] = pd.to_datetime(ts.date, utc=True).dt.year
        ann = ts.groupby("year").chainage_m.agg(["median", "size"])
        annual[tid] = (g, ann[ann["size"] >= BUXTON_MIN_OBS]["median"])

    # A straight-line fit over the years given: (slope m/yr, intercept)
    def fit(sr, a, b):
        x = sr.loc[a:b].dropna()
        return np.polyfit(x.index, x.to_numpy(), 1)

    for window, w in BUXTON_WINDOWS.items():
        rows = []
        for tid, (g, ann) in annual.items():
            base = ann.loc[w["base"][0]:w["base"][1]]
            if len(base) < 2:
                continue
            rel = ann.loc[w["years"][0]:w["years"][1]] - base.median()
            rows += [dict(transect=tid, gis=g, year=y, position_m=v) for y, v in rel.items()]
        t = pd.DataFrame(rows)
        per_gis = t.groupby(["gis", "year"]).position_m.median().unstack("gis")
        n_transects = t.groupby("gis").transect.nunique()
        stretch = t.groupby("year").position_m.median()
        (OUT / "tables").mkdir(parents=True, exist_ok=True)
        per_gis.assign(gis_10_15=stretch).round(2).to_csv(OUT / "tables" / f"buxton_coastsat_annual_{window}.csv")

        # The change across each storm or fill, and what the window's trend alone gives over the same years
        slope0 = fit(stretch, *w["trend"])[0]
        reaches = [*((f"GIS {g}", per_gis[g]) for g in per_gis), ("GIS 10-15", stretch)]
        steps, trends = [], []
        for name, sr in reaches:
            for label, before, after in w["steps"]:
                change = float(sr.loc[after[0]:after[1]].mean() - sr.loc[before[0]:before[1]].mean())
                years = np.mean(after) - np.mean(before)
                steps.append(dict(reach=name, event=label, change_m=change, trend_expected_m=slope0 * years,
                                  beyond_trend_m=change - slope0 * years))
            trends.append(dict(reach=name, **{f"trend_{a}_{b}_m_yr": fit(sr, a, b)[0] for a, b in w["trends"]}))
        steps, trends = pd.DataFrame(steps), pd.DataFrame(trends)
        steps.round(2).to_csv(OUT / "tables" / f"buxton_steps_{window}.csv", index=False)
        trends.round(2).to_csv(OUT / "tables" / f"buxton_trends_{window}.csv", index=False)
        pd.set_option("display.width", 220)
        print(f"\n{window}  transects per GIS: {n_transects.to_dict()}")
        print(steps.pivot(index="reach", columns="event", values="change_m")
              .loc[[r for r, _ in reaches], [e for e, *_ in w["steps"]]].round(1).to_string())
        print(trends.set_index("reach").round(2).to_string())

        # The model's shoreline over the same GIS, seaward + from its start (matrix rows are the start year, ...)
        model = {}
        for v in w["model"]:
            rd = H.run_dir(v)
            m = np.load(next(rd.glob("*_shoreline_matrix.npy")))
            dom = DOM.for_total(m.shape[1])
            pads = [dom.gis_to_pad(g) for g in BUXTON_GIS]
            model[v] = pd.Series(np.median(-(m[:, pads] - m[0, pads]), axis=1),
                                 index=H.RUN_START + np.arange(len(m)))

        fig, axs = plt.subplots(2, 1, figsize=(10, 8.2), sharex=True)
        shades = ["#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C"]
        ax = axs[0]
        for g, col in zip(per_gis, shades):
            ax.plot(per_gis.index, per_gis[g], "o-", color=col, lw=1.4, ms=3.5, label=f"GIS {g} ({n_transects[g]})")
        ax.set_ylabel(f"CoastSat shoreline position (m)\nseaward +, from the {w['base'][0]}-{w['base'][1] % 100:02d} "
                      "median")
        ax.legend(fontsize=7, ncol=3, loc="lower left", title="transects per GIS", title_fontsize=7)
        ax = axs[1]
        slope, icpt = fit(stretch, *w["trend"])
        yrs = np.arange(w["trend"][0], w["years"][1] + 1)
        ax.plot(yrs, np.polyval([slope, icpt], yrs), "--", color="0.45", lw=1.2,
                label=f"{w['trend'][0]}-{w['trend'][1] % 100:02d} trend extended ({slope:+.2f} m/yr)")
        ax.plot(stretch.index, stretch, "o-", color="#08519C", lw=2.2, ms=4, label="CoastSat, GIS 10-15 median")
        styles = {"cal_control": ("s-", "#8C8C8C", "model, no collision"),
                  "cal_r03": ("^-", "#D95F02", "model, collision module")}
        # the model starts at the CoastSat value in its start year so the two share a zero
        off = float(stretch.get(H.RUN_START, 0.0))
        for v, sr in model.items():
            mk, col, label = styles[v]
            ax.plot(sr.index, sr + off, mk, color=col, lw=1.7, ms=3.5, label=label)
        ax.set_ylabel("Shoreline position (m), seaward +")
        ax.legend(fontsize=7, loc="lower left")
        for ax, k in zip(axs, "ab"):
            ax.axhline(0, color="0.75", lw=0.8)
            for name, (x0, x1) in w["fills"].items():
                ax.axvspan(x0, x1, color="#E6AB02", alpha=0.25, lw=0)
                if k == "a":
                    ax.text(x0, 0.03, name, transform=ax.get_xaxis_transform(), fontsize=7, color="#8C6D00")
            for name, x in w["storms"].items():
                ax.axvline(x, color="#B2182B", lw=0.9, ls=":")
                if k == "a":
                    ax.text(x + 0.1, 0.97, name, transform=ax.get_xaxis_transform(), fontsize=7, color="#B2182B",
                            va="top")
            ax.set_title(f"({k})", loc="left", fontsize=9)
        axs[1].set_xlabel("Year")
        fig.tight_layout()
        out = OUT / "figures" / f"buxton_coastsat_{window}.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=200, facecolor="white")
        plt.close(fig)
        fills = (" Shaded: the Buxton nourishment placements (GIS 6-16), " + ", ".join(w["fills"]) + ".") if w["fills"] else ""
        record_caption(out, (
            f"Is the Buxton-Avon (GIS 10-15) erosion steady or in steps at the storms, {window.replace('_', '-')}? "
            "(a) CoastSat shoreline position per GIS: per transect the annual median chainage (years with at least "
            f"{BUXTON_MIN_OBS} images) less its {w['base'][0]}-{w['base'][1]} median, then the median over the "
            f"GIS's transects. (b) The GIS 10-15 median with its {w['trend'][0]}-{w['trend'][1]} trend extended"
            + (", and the model's shoreline over the same GIS (median of the domains, seaward +) without and with "
               "the collision module, started at the CoastSat value in its start year" if model else "")
            + f". Dotted lines: hurricanes {', '.join(w['storms'])}.{fills} Experiment {TAG}."))
        print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("action", choices=["growth", "runs", "groin", "hs", "gsweep", "buxton"])
    ap.add_argument("runs", nargs="*")
    a = ap.parse_args()
    {"growth": cmd_growth, "runs": cmd_runs, "groin": cmd_groin, "hs": cmd_hs, "gsweep": cmd_gsweep, "buxton": cmd_buxton}[a.action](a)


if __name__ == "__main__":
    main()
