"""
Does the shoreline offset's averaging window (v1 calendar, v2 DEM-centred) change what the model does?

    python scripts/analyze_output/compare_runs/offset_source/shoreline_v1_vs_v2_comparison.py

Reads the full-management shoreline_v1 / shoreline_v2 pairs of the 2026-09-29
adopted-setup offset study and their offsets; writes figures and tables to
output/comparisons/offset_source/shoreline_v1_vs_v2/.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))

import offset_source_comparison as osc  # noqa: E402  (paths, style, offsets, target)
from offset_source_comparison import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, INTERIOR, PERIODS, YEARS, _REPO, _title,
    apply_style, figsize, open_frame, plt, record_caption, save, town_bands,
)
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from site_layer.hat_figure_style import C  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
STUDY = (_REPO / "output" / "raw_runs" / "experiments" / "island-offset"
         / "2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup")
OUT = osc.OUT / "shoreline_v1_vs_v2"
ARMS = ("v1", "v2")
COL = {"v1": C["EARLY"], "v2": C["LATE"]}   # house vintage pair: earlier window red, later blue
# -----------------------------------------------------------------------------


# The one run folder of a version and period
def run_dir(v, start, end):
    root = STUDY / "runs" / f"shoreline_{v}_full_management" / f"{start}_{end}" / "zeroBE"
    hits = [d for d in root.glob("*") if d.is_dir()]
    if len(hits) != 1:
        raise SystemExit(f"{root}: expected one run, found {len(hits)}")
    return hits[0]


# Per-domain table and per-period summary for both versions
def build():
    per, rows = [], []
    for start, end in PERIODS:
        cols = {}
        for v in ARMS:
            d = run_dir(v, start, end)
            got = osc.offset_version(d, "shoreline")
            if got != v:
                raise ValueError(f"{d}: ran on shoreline/{got}, filed as {v}")
            s = osc.seaward("shoreline", start, v)
            th, tu = osc.orientation(s)
            rt = pd.read_csv(d / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
            cols.update({f"offset_{v}_m": s, f"orientation_{v}_deg": th,
                         f"turning_{v}_deg_km": tu, f"model_change_{v}_m": rt.lrr_m_yr * YEARS})
        t = pd.DataFrame(cols)
        t.index.name = "gis_domain"
        for q, unit in (("offset", "m"), ("orientation", "deg"), ("turning", "deg_km"),
                        ("model_change", "m")):
            t[f"d_{q}_{unit}"] = t[f"{q}_v2_{unit}"] - t[f"{q}_v1_{unit}"]
        t.insert(0, "period", f"{start}_{end}")
        per.append(t.reset_index())

        i = t.loc[INTERIOR[0]:INTERIOR[1]]
        dm = i.d_model_change_m
        rows.append(dict(
            period=f"{start}-{end}",
            model_sd_v1_m=i.model_change_v1_m.std(), model_sd_v2_m=i.model_change_v2_m.std(),
            r_between_runs=np.corrcoef(i.model_change_v1_m, i.model_change_v2_m)[0, 1],
            mean_diff_m=dm.mean(), mean_abs_diff_m=dm.abs().mean(), sd_diff_m=dm.std(),
            max_abs_diff_m=dm.abs().max(), max_at_gis=int(dm.abs().idxmax()),
            diff_sd_over_model_sd=dm.std() / i.model_change_v1_m.std(),
            sd_d_offset_m=i.d_offset_m.std(), max_abs_d_offset_m=i.d_offset_m.abs().max(),
            r_d_offset_vs_d_model=np.corrcoef(i.d_offset_m, dm)[0, 1],
            r_d_orientation_vs_d_model=np.corrcoef(i.d_orientation_deg, dm)[0, 1],
            r_d_turning_vs_d_model=np.corrcoef(i.d_turning_deg_km, dm)[0, 1],
            v1_run=str(run_dir("v1", start, end).relative_to(_REPO)).replace("\\", "/"),
            v2_run=str(run_dir("v2", start, end).relative_to(_REPO)).replace("\\", "/"),
        ))
    return pd.concat(per, ignore_index=True), pd.DataFrame(rows)


# Each run against the projected target, interior: bias, RMS residual, explained, r
def vs_target(t, obs):
    rows = []
    for start, end in PERIODS:
        p = t[t.period == f"{start}_{end}"].set_index("gis_domain").loc[INTERIOR[0]:INTERIOR[1]]
        o = obs.reindex(p.index)
        for v in ARMS:
            m = p[f"model_change_{v}_m"]
            res = m - o
            rows.append(dict(period=f"{start}-{end}", offset=f"shoreline/{v}", bias_m=res.mean(),
                             rms_residual_m=np.sqrt((res ** 2).mean()),
                             variance_explained=1 - (res ** 2).sum() / ((o - o.mean()) ** 2).sum(),
                             r=np.corrcoef(m, o)[0, 1]))
    return pd.DataFrame(rows)


# Caption text shared by every figure
COMMON = (" Full management, option A waves (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle "
          "0.5), adopted setup (per-cell dune ceilings, beach/dune cap on added sand only, "
          "v3_split12_trim24 storms), no source/sink correction at the ends (zeroBE), "
          "relocations and groins off; within each period the two runs differ ONLY in the "
          "shoreline offset's averaging window: v1 the CoastSat mean over calendar 1995-1997 / "
          "2009-2011, v2 over +/-1 yr of the start DEM's lidar flights (1995-10-12 to "
          "1997-10-12 / 2008-08-17 to 2010-08-17). Domain 1 is Cape Point, 90 Pea Island; 500 m "
          "per domain. Runs from raw_runs/experiments/island-offset/2026-09-29-shoreline-offset-"
          "v1-vs-v2-adopted-setup. Statistics in tables/summary.csv, interior GIS 2-89.")


# Profiles: the two runs' change with the target, offset difference, change difference
def fig_profiles(t, obs):
    f, axes = plt.subplots(3, 2, figsize=figsize("double", height=6.6), sharex=True,
                           sharey="row", constrained_layout=True)
    for j, (start, end) in enumerate(PERIODS):
        p = t[t.period == f"{start}_{end}"].set_index("gis_domain")
        a0, a1, a2 = axes[:, j]
        a0.plot(obs.index, obs.values, color=INK, lw=1.6, zorder=5)
        for v in ARMS:
            a0.plot(p.index, p[f"model_change_{v}_m"], color=COL[v], lw=1.1, zorder=4)
        a1.bar(p.index, p.d_offset_m, width=0.85, color=INK_MUTED)
        a2.bar(p.index, p.d_model_change_m, width=0.85,
               color=[COL["v2"] if x > 0 else COL["v1"] for x in p.d_model_change_m])
        for k, ax in enumerate((a0, a1, a2)):
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(0.5, 90.5)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(k == 0))
            _title(ax, 2 * k + j, f"{start}–{end}")
        a2.set_xlabel(DOMAIN_AXIS_LABEL)
        if j == 0:
            a0.set_ylabel("Shoreline change (m)")
            a1.set_ylabel("Offset difference,\nv2 − v1 (m, seaward +)")
            a2.set_ylabel("Change difference,\nv2 − v1 (m)")
    f.legend(handles=[Line2D([], [], color=INK, lw=1.6, label="Projected CoastSat position"),
                      Line2D([], [], color=COL["v1"], lw=1.1, label="Shoreline v1 (calendar window)"),
                      Line2D([], [], color=COL["v2"], lw=1.1, label="Shoreline v2 (DEM-centred window)"),
                      Patch(color=COL["v2"], label="v2 start more accretional"),
                      Patch(color=COL["v1"], label="v1 start more accretional")],
             loc="outside lower center", ncol=3, frameon=False)
    png = OUT / "shoreline_v1_vs_v2_difference_full_management.png"
    save(f, png, close=True)
    s = t.attrs["summary"].set_index("period")
    record_caption(png, (
        "Model outcome with the shoreline island offset averaged over the calendar window (v1) "
        "or the DEM-centred window (v2). (a, b) Modelled total shoreline change, each run's LRR "
        "x 14 yr, seaward positive: red v1, blue v2; black, projected shoreline change (CoastSat "
        "LRR 1996-2024, LOWESS 7 domains, southern 10 raw, x 14 yr), the same profile in both "
        "panels; the model is unsmoothed. (c, d) The difference in the offset itself, v2 minus "
        "v1, seaward position with each build's mean removed (m). (e, f) The difference in "
        "modelled change, v2 start minus v1 start; blue where the v2 start is more accretional. "
        "Interior: r between the runs {:.2f} / {:.2f}, mean |difference| {:.1f} / {:.1f} m, "
        "largest {:.1f} m at GIS {} / {:.1f} m at GIS {} (1996-2010 / 2010-2024).".format(
            s.loc["1996-2010", "r_between_runs"], s.loc["2010-2024", "r_between_runs"],
            s.loc["1996-2010", "mean_abs_diff_m"], s.loc["2010-2024", "mean_abs_diff_m"],
            s.loc["1996-2010", "max_abs_diff_m"], s.loc["1996-2010", "max_at_gis"],
            s.loc["2010-2024", "max_abs_diff_m"], s.loc["2010-2024", "max_at_gis"]) + COMMON))
    return png


# Change difference against offset and turning difference, per domain
def fig_scatter(t, summary):
    f, axes = plt.subplots(2, 2, figsize=figsize("double", height=5.4), sharey=True,
                           sharex="col", constrained_layout=True)
    for j, (start, end) in enumerate(PERIODS):
        p = t[t.period == f"{start}_{end}"].set_index("gis_domain").loc[INTERIOR[0]:INTERIOR[1]]
        c = C["EARLY"] if j == 0 else C["LATE"]
        for k, (col, lab) in enumerate((("d_offset_m", "Offset difference (m, seaward +)"),
                                        ("d_turning_deg_km", "Turning difference (° per km)"))):
            ax = axes[j, k]
            x, y = p[col].to_numpy(), p.d_model_change_m.to_numpy()
            ax.scatter(x, y, s=10, color=c, alpha=0.8, lw=0)
            b = np.polyfit(x, y, 1)
            xx = np.linspace(x.min(), x.max(), 2)
            ax.plot(xx, np.polyval(b, xx), color=INK, lw=0.8)
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.axvline(0, color=INK_MUTED, lw=0.6)
            ax.grid(axis="y")
            open_frame(ax)
            _title(ax, 2 * j + k, f"{start}–{end}")
            if j == 1:
                ax.set_xlabel(lab + ", v2 − v1")
            if k == 0:
                ax.set_ylabel("Change difference (m)")
    png = OUT / "shoreline_v1_vs_v2_offset_vs_model_full_management.png"
    save(f, png, close=True)
    s = summary.set_index("period")
    record_caption(png, (
        "Which property of the offset change drives the difference between the v1 and v2 runs. "
        "Each point is one interior domain (GIS 2-89); y is the modelled change with the v2 "
        "start minus that with the v1 start (m). Left (a, c): against the change in the "
        "offset's seaward position (mean removed). Right (b, d): against the change in "
        "turning, the alongshore rate of change of the offset's angle (degrees per km; positive "
        "at an embayment). Line: least squares. r, offset: {:+.2f}, {:+.2f}; r, turning: "
        "{:+.2f}, {:+.2f} (1996-2010, 2010-2024).".format(
            s.loc["1996-2010", "r_d_offset_vs_d_model"], s.loc["2010-2024", "r_d_offset_vs_d_model"],
            s.loc["1996-2010", "r_d_turning_vs_d_model"],
            s.loc["2010-2024", "r_d_turning_vs_d_model"]) + COMMON))
    return png


# Run: build the tables, score against the target, draw two figures
def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    apply_style()
    t, summary = build()
    t.attrs["summary"] = summary
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / "tables" / "per_domain.csv", index=False)
    summary.to_csv(OUT / "tables" / "summary.csv", index=False)
    obs = osc.projected_target()
    scores = vs_target(t, obs)
    scores.to_csv(OUT / "tables" / "vs_projected.csv", index=False)
    for png in (fig_profiles(t, obs), fig_scatter(t, summary)):
        print(png.relative_to(_REPO))
    print(summary.drop(columns=["v1_run", "v2_run"]).round(2).T.to_string())
    print(scores.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
