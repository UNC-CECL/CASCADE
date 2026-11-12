"""
The collision module's Cs sweep with the toe at the berm, for two ways of taking the eroded volume from the dune
(cascade/dune_collision.py, rows): "proportional", the crest lowered and both rows in proportion (the module's
default), and "lateral", the face retreating through the front row, which keeps its height until its sand above the
toe is gone (2026-11-11). Which Cs fits each, and does lateral erosion catch the sudden dune losses?

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_row_mode_sweep.py

Runs (HAT_dune_collision.py, 1996-2009, all with the road fix, dune growth 0.02-0.04, groin 0.16, no road rule,
Rodanthe 40%, toe at the berm, 1.7 m NAVD88; only the module differs from the off run): cal_step2_growth_noscrape
(module off), cal_toeberm_cs* (proportional) and cal_toeberm_lat_cs* (lateral). Runs not yet made are skipped.

Writes, in the workspace's dune-collision/4-hindcast/3-collision-design/ (D:):
- lateral_vs_proportional_sensitivity.png and supporting/lateral_vs_proportional_sensitivity.csv: crest RMSE (GIS
  6-90 and 9-15), 2009 mean offset, share within 0.5 m, overwash hits and false alarms, shoreline RMSE, against Cs
- lateral_collapse_<variant>.png: collision volume and front-row cells collapsed by domain and year, for the best
  lateral Cs and the strongest one
- lateral_crests_over_time_GIS10-14.png: the domain-median crest through time, off, best proportional, best lateral

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-11
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402

sys.path.insert(0, str(P._HAT_DIR / "dune-collision" / "4-hindcast" / "1-calibration-1996-2009"))
import HAT_dune_collision as H  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = P.OUT / "3-collision-design"
OFF = "cal_step2_growth_noscrape"
PROPORTIONAL = {"cal_toeberm_cs0p1": 1e-7, "cal_toeberm_cs0p3": 3e-7, **{f"cal_toeberm_cs{c}": c * 1e-6
                for c in (1, 2, 5, 10, 30, 100)}}
LATERAL = {f"cal_toeberm_lat_cs{c}": c * 1e-6 for c in (1, 3, 10, 20, 30, 50, 70, 100, 150, 200, 250)}
MODES = (("Proportional (crest lowered)", PROPORTIONAL, "#5E3C99"), ("Lateral (face retreat)", LATERAL, "#E66101"))
GIS_LO, GIS_HI = 6, 90                              # where the module runs
TS_GIS = range(10, 15)                              # the crest time series
COLLAPSE_MIN_M = 0.5                                # a collapse counts from at least this far above the toe (as
                                                    # DuneCollisionManager.COLLAPSE_MIN_HEIGHT_DAM)
# -----------------------------------------------------------------------------


def _done(v):
    rd = H.run_dir(v)
    return rd is not None and (rd / "tables").exists()


def metrics(v, cs):
    rd = H.run_dir(v)
    from cascade_pipeline.plotting.dune_crest_comparison import TABLE, TABLES
    t = pd.read_csv(rd / TABLES / TABLE)
    t = t[t.gis.between(GIS_LO, GIS_HI) & (t.survey > 1996)]       # 1996 is the shared start
    f = lambda x: float(np.sqrt((x ** 2).mean()))
    sk = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))["skill"]
    log = (H.EXP_DIR / "logs" / f"{v}.log").read_text(encoding="utf-8", errors="ignore")
    ow = re.findall(r"overwash vs observed\s+hit ([\d.]+)\s+false alarm ([\d.]+)", log)
    hit, fa = (100 * float(ow[-1][0]), 100 * float(ow[-1][1])) if ow else (np.nan, np.nan)
    return dict(variant=v, cs=cs, crest_rmse_m=f(t.model_minus_survey_m),
                crest_rmse_gis9_15_m=f(t[t.gis.between(9, 15)].model_minus_survey_m),
                bias_2009_m=float(t[t.survey == 2009].model_minus_survey_m.mean()),
                within_0p5_pct=float(100 * (t.model_minus_survey_m.abs() <= 0.5).mean()),
                overwash_hit_pct=hit, overwash_false_alarm_pct=fa, shoreline_rmse_m_yr=float(sk["rmse_interior_m_yr"]))


def sensitivity(S, plt):
    off = metrics(OFF, 0.0)
    rows = [dict(mode="off", **off)]
    for name, runs, _ in MODES:
        rows += [dict(mode=name, **metrics(v, cs)) for v, cs in runs.items() if _done(v)]
    t = pd.DataFrame(rows)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(4).to_csv(OUT / "supporting" / "lateral_vs_proportional_sensitivity.csv", index=False)
    panels = (("crest_rmse_m", "Crest RMSE, GIS 6–90 (m)"), ("crest_rmse_gis9_15_m", "Crest RMSE, GIS 9–15 (m)"),
              ("bias_2009_m", "Mean crest offset, 2009 (m)"), ("within_0p5_pct", "Domain-surveys within ±0.5 m (%)"),
              ("overwash_false_alarm_pct", "Overwash false alarms (%)"), ("shoreline_rmse_m_yr", "Shoreline RMSE (m/yr)"))
    fig, axs = plt.subplots(3, 2, figsize=(S.FIG_W_DOUBLE, 7.8), sharex=True, layout="constrained")
    for k, (ax, (col, lab)) in enumerate(zip(axs.ravel(), panels)):
        ax.axhline(off[col], color="0.45", lw=1, ls=(0, (4, 2)), label="Module off" if k == 0 else None)
        if col == "bias_2009_m":
            ax.axhline(0, color="0.2", lw=0.7)
        for name, _, colour in MODES:
            m = t[t["mode"] == name].sort_values("cs")
            ax.plot(m.cs, m[col], "-o", color=colour, ms=4, lw=1.8, label=name if k == 0 else None)
        ax.set_xscale("log")
        ax.set_ylabel(lab, fontsize=8.5)
        ax.grid(color="0.92", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.13, 1.02, f"({'abcdef'[k]})", transform=ax.transAxes, fontsize=10, fontweight="bold")
    for ax in axs[-1]:
        ax.set_xlabel("Cs (log scale)")
    fig.legend(loc="outside lower center", ncol=3, fontsize=8, frameon=False)
    out = OUT / "lateral_vs_proportional_sensitivity.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, (
        "Sensitivity of the 1996-2009 hindcast to the collision module's Cs (Larson et al. 2004 alone, toe at the berm, "
        "1.7 m NAVD88; road fix, growth 0.02-0.04, groin 0.16, no road rule, Rodanthe 40%; only the module differs), "
        "for two ways of taking the eroded volume from the dune: proportional (the crest lowered, both rows in "
        "proportion; purple) and lateral (the face retreats through the front row, which keeps its height until its "
        "sand above the toe is gone; orange). Dashed: module off. (a, b) Dune crest RMSE against the lidar domain "
        "medians over the surveys after 1996, GIS 6-90 and 9-15; (c) mean model minus lidar crest in 2009; (d) share "
        "of domain-surveys within 0.5 m; (e) overwash false-alarm rate against the imagery; (f) shoreline RMSE "
        "against CoastSat."))
    print(f"wrote {out}")
    return t


def recount_collapses(v):
    """Front-row collapses per domain and year, counted again from the saved run: a cell-year counts once when the
    module recorded a front-row collapse and the cell started that year at least COLLAPSE_MIN_M above its toe. The
    runs of 2026-11-11 were made before the module's own threshold, and counted a cell already flattened, and topped
    up by Barrier3D's alongshore diffusion, again at every storm."""
    from cascade_pipeline.saved_runs import load_saved_cascade
    rd = H.run_dir(v)
    c = load_saved_cascade(rd / f"{rd.name}.npz")
    rows = []
    for pad, (b, m) in enumerate(zip(c.barrier3d, c.dune_collision)):
        if m is None or getattr(m, "dune_row_loss_TS", None) is None:
            continue
        gis = pad - (len(c.barrier3d) - 90) // 2 + 1
        front = np.asarray(b.DuneDomain, float)[:, :, 0]                      # dam above the berm
        above_toe = (front + b._BermEl - np.asarray(m.toe_elevation)[None, :]) * 10   # m
        lost = m.dune_row_loss_TS[:, :, 0] > 0
        for t in range(1, len(m.crest_lowering_TS) + 1):
            rows.append(dict(year=H.RUN_START + t, gis=gis,
                             front_rows_collapsed=int((lost[t] & (above_toe[t - 1] >= COLLAPSE_MIN_M)).sum())))
    return pd.DataFrame(rows)


def collapse(variants):
    (OUT / "logs").mkdir(exist_ok=True)
    for v in variants:
        src = H.EXP_DIR / "logs" / f"collision_cells_{v}.csv"
        if src.exists():
            t = pd.read_csv(src).drop(columns="front_rows_collapsed").merge(recount_collapses(v), on=["year", "gis"])
            t.to_csv(OUT / "logs" / src.name, index=False)
            print(f"{v}: {t.front_rows_collapsed.sum()} front-row collapses after the recount "
                  f"(the run's own count: {pd.read_csv(src).front_rows_collapsed.sum()})")
    H.EXP_DIR, saved = OUT, H.EXP_DIR        # plot_collapse reads logs/ and writes beside it: point it here
    H.plot_collapse(variants)
    for v in variants:                       # the house names: lateral_collapse_<variant>.png
        p = OUT / f"collision_collapse_{v}.png"
        if p.exists():
            p.replace(OUT / f"lateral_collapse_{v}.png")
    H.EXP_DIR = saved


def crests_over_time(S, plt, runs):
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    sv = pd.read_csv(H.SURVEY_CRESTS)
    sv = sv[(sv.status == "ok") & sv.gis.isin(TS_GIS)].copy()
    t = pd.to_datetime(sv.flown.astype(str).str[:7], errors="coerce")
    sv["yr"] = t.dt.year + (t.dt.month - 0.5) / 12
    lid = sv.groupby(["gis", "survey"]).agg(z=("crest_m_navd88", "median"), yr=("yr", "first")).reset_index()
    fig, axs = plt.subplots(len(TS_GIS), 1, figsize=(S.FIG_W_DOUBLE * 0.75, 8.6), sharex=True, sharey=True,
                            layout="constrained")
    for v, lab, colour in runs:
        c = load_saved_cascade(H.run_dir(v) / f"{H.run_dir(v).name}.npz")
        d = DOM.for_total(len(c.barrier3d))
        for ax, g in zip(axs, TS_GIS):
            b = c.barrier3d[d.gis_to_pad(g)]
            z = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + H.MHW_NAVD88_M
            ax.plot(H.RUN_START + 1 + np.arange(z.shape[0]), np.median(z, axis=1), color=colour, lw=2, label=lab)
    for ax, g in zip(axs, TS_GIS):
        L = lid[lid.gis == g].sort_values("yr")
        ax.plot(L.yr, L.z, "o", color="#1A1A1A", ms=5, label="Lidar (median)", zorder=5)
        ax.axhline(H.TOE_AT_BERM_M, color="#7B3F99", lw=1, ls=(0, (4, 2)), label="Berm (the toe)")
        for x in (1998.65, 1999.67, 2003.71):
            ax.axvline(x, color="0.6", lw=0.6, ls=":")
        ax.set_title(f"GIS {g}", loc="left", fontsize=9, fontweight="bold")
        ax.set_ylabel("Crest (m NAVD88)", fontsize=8)
        ax.grid(axis="y", color="0.9", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    for x, n in ((1998.65, "Bonnie"), (1999.67, "Dennis"), (2003.71, "Isabel")):
        axs[0].text(x, 1.02, n, transform=axs[0].get_xaxis_transform(), ha="center", fontsize=7, color="0.4")
    axs[-1].set_xlabel("Year")
    axs[-1].set_xlim(1996, 2010)
    axs[0].set_ylim(0, 9)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=3, fontsize=7.5, frameon=False)
    out = OUT / "lateral_crests_over_time_GIS10-14.png"
    S.save(fig, out, bbox_inches="tight", close=True)
    S.record_caption(out, ("Domain-median dune crest through 1996-2009, GIS 10-14: model with the collision module off "
                           "and on at the best-fitting Cs of each way of taking the eroded volume (proportional, "
                           "lateral), model state at 1 January each year, against the lidar domain median at each "
                           "survey; the berm, the module's toe (dashed). Dotted: Bonnie, Dennis, Isabel."))
    print(f"wrote {out}")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as S
    S.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    t = sensitivity(S, plt)
    best = {name: t[t["mode"] == name].sort_values("crest_rmse_m").iloc[0] for name, _, _ in MODES
            if (t["mode"] == name).any()}
    for name, r in best.items():
        print(f"best {name}: Cs {r.cs:.1e} ({r.variant}), crest RMSE {r.crest_rmse_m:.3f} m, "
              f"GIS 9-15 {r.crest_rmse_gis9_15_m:.3f} m, 2009 offset {r.bias_2009_m:+.2f} m")
    lat = [v for v in LATERAL if _done(v)]
    if lat:
        strongest = max(lat, key=LATERAL.get)
        collapse(sorted({best["Lateral (face retreat)"].variant, strongest}))
    runs = [(OFF, "Module off", "#8C8C8C")] + [(r.variant, f"{name}, Cs {r.cs:.0e}", colour)
                                               for (name, _, colour) in MODES if name in best
                                               for r in [best[name]]]
    crests_over_time(S, plt, runs)


if __name__ == "__main__":
    main()
