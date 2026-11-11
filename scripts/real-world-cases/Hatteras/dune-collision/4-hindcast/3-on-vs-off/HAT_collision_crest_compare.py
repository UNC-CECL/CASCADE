"""
Modelled against surveyed dune crests, 1996-2009, with the collision module off and on: does the module bring the
dunes closer to the lidar?

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_crest_compare.py

(a, b) Model minus lidar domain-median dune crest, each GIS domain (columns) and lidar survey (rows), on the run
overviews' bands (cascade_pipeline.plotting.dune_crest_comparison), for the base and the module on.
(c) The island crest error at each survey (RMSE over GIS 6-90, the module is off at Cape Point), and (d) the mean
offset (bias), the same runs.

Runs (HAT_dune_collision.py; nothing rerun): cal_step1_base (the base, before any change) and
cal_step3_module_la_alltoes (Larson alone, Cs 1.9e-5, dune growth 0.02-0.04, one toe per domain averaged over every
1996-2009 survey, 1998 left out). Crests: each run's
tables/dune_crest_vs_lidar.csv.

Writes, in the workspace's dune-collision/4-hindcast/3-on-vs-off/ (D:): crest_offset_base_vs_module_1996_2009.png
and supporting/crest_offset_base_vs_module_1996_2009.csv.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_collision_on_off_position as P  # noqa: E402
OUT = P.OUT / "1-module-on-off"   # this script's figures (topic subfolder)

sys.path.insert(0, str(P._HAT_DIR / "dune-collision" / "4-hindcast" / "1-calibration-1996-2009"))

# --- CONFIG ------------------------------------------------------------------
RUNS = (("cal_step1_base", "Base (before any changes)", "#8C8C8C", "-"),
        ("cal_step3_module_la_alltoes", "Collision module, toes from all 1996–2009 surveys", "#652381", "-"))
GIS_LO, GIS_HI = 6, 90
OUT_NAME = "crest_offset_base_vs_module_1996_2009.png"
SURVEY_PANELS = False   # the error and bias at each survey under the maps (dropped 2026-11-05; the CSV keeps them)
# -----------------------------------------------------------------------------


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    import HAT_dune_collision as C
    from cascade_pipeline.plotting.dune_crest_comparison import NO_SURVEY_COLOR, OFFSET_BINS_M, TABLE, TABLES
    from site_layer import hat_figure_style as S

    S.apply_style()
    tabs = {}
    for v, *_ in RUNS:
        t = pd.read_csv(C.run_dir(v) / TABLES / TABLE)
        tabs[v] = t[t.gis.between(GIS_LO, GIS_HI)]
    surveys = sorted(tabs[RUNS[0][0]].survey.unique())
    flown = tabs[RUNS[0][0]].drop_duplicates("survey").set_index("survey").flown

    half = (len(OFFSET_BINS_M) - 2) // 2
    blues = plt.get_cmap("Blues_r")(np.linspace(0.0, 0.75, half + 1))
    reds = plt.get_cmap("Reds")(np.linspace(0.25, 1.0, half + 1))
    cmap = ListedColormap(list(blues[1:]) + ["white"] + list(reds[:-1]))
    cmap.set_under(blues[0])
    cmap.set_over(reds[-1])
    cmap.set_bad(NO_SURVEY_COLOR)
    norm = BoundaryNorm(OFFSET_BINS_M, cmap.N)

    n = len(RUNS)
    if SURVEY_PANELS:
        fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 3.6 + 1.8 * n))
        gs = fig.add_gridspec(n + 3, 2, height_ratios=(1,) * n + (0.08, 0.12, 1.15), width_ratios=(1, 1),
                              hspace=0.42, wspace=0.28, top=0.91, bottom=0.09, left=0.08, right=0.98)
    else:
        fig = plt.figure(figsize=(S.FIG_W_DOUBLE, 0.9 + 2.0 * n))
        gs = fig.add_gridspec(n + 1, 2, height_ratios=(1,) * n + (0.09,), hspace=0.42, top=1 - 0.6 / (0.9 + 2.0 * n), bottom=0.10,
                              left=0.08, right=0.98)
    rows = []
    for k, (v, label, colour, ls) in enumerate(RUNS):
        ax = fig.add_subplot(gs[k, :])
        grid = tabs[v].pivot_table(index="survey", columns="gis", values="model_minus_survey_m").reindex(
            index=surveys, columns=range(GIS_LO, GIS_HI + 1))
        im = ax.imshow(np.ma.masked_invalid(grid.values), aspect="auto", cmap=cmap, norm=norm, origin="lower",
                       extent=(GIS_LO - 0.5, GIS_HI + 0.5, -0.5, len(surveys) - 0.5), interpolation="nearest")
        ax.set_yticks(range(len(surveys)), [str(s) for s in surveys], fontsize=7)
        ax.set_xticks(range(10, 91, 10))
        ax.tick_params(axis="x", labelbottom=(k == len(RUNS) - 1), labelsize=7.5)
        if k == len(RUNS) - 1:
            ax.set_xlabel(S.DOMAIN_AXIS_LABEL)
        ax.set_ylabel("Lidar survey", fontsize=8)
        ax.set_title(f"({chr(97 + k)}) {label}", loc="left", fontsize=9, fontweight="bold")
        d = tabs[v].model_minus_survey_m
        ax.set_title(f"RMSE {np.sqrt((d ** 2).mean()):.2f} m, within ±0.5 m {100 * (d.abs() <= 0.5).mean():.0f}%",
                     loc="right", fontsize=8, color=S.INK_MUTED)
        g = tabs[v].assign(err2=lambda x: x.model_minus_survey_m ** 2).groupby("survey")
        rows.append(pd.DataFrame({"run": v, "survey": surveys,
                                  "rmse_m": np.sqrt(g.err2.mean()).reindex(surveys).to_numpy(),
                                  "bias_m": g.model_minus_survey_m.mean().reindex(surveys).to_numpy(),
                                  "within_0p5_pct": (100 * g.model_minus_survey_m.apply(
                                      lambda x: (x.abs() <= 0.5).mean())).reindex(surveys).to_numpy()}))
    cax = fig.add_subplot(gs[n, :])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both", ticks=OFFSET_BINS_M)
    cb.set_label("Model − lidar, domain-median dune crest (m)   ← model lower   |   model higher →", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)

    t = pd.concat(rows, ignore_index=True)
    year = {s: int(str(f)[:4]) + (int(str(f)[5:7]) - 0.5) / 12 for s, f in flown.items()}
    t["year"] = t.survey.map(year)
    for j, (col, ylab, ref) in enumerate((("rmse_m", "Crest RMSE (m)", None), ("bias_m", "Mean offset (m)", 0.0))
                                         if SURVEY_PANELS else ()):
        ax = fig.add_subplot(gs[n + 2, j])
        if ref is not None:
            ax.axhline(ref, color=S.INK, lw=0.7)
        for x, name in ((1998.65, "Bonnie"), (1999.67, "Dennis"), (2003.71, "Isabel")):
            ax.axvline(x, color=S.INK_MUTED, lw=0.6, ls=":")
            if j == 0:
                ax.text(x - 0.12, 0.97, name, transform=ax.get_xaxis_transform(), rotation=90, ha="right",
                        va="top", fontsize=6.5, color=S.INK_MUTED)
        for v, label, colour, ls in RUNS:
            r = t[t.run == v].sort_values("year")
            ax.plot(r.year, r[col], color=colour, ls=ls, lw=1.8, marker="o", ms=3.5, label=label)
        ax.set_ylabel(ylab)
        ax.set_xlim(1996, 2010)
        ax.set_xticks([1996, 1999, 2002, 2005, 2008])
        ax.grid(axis="y", color=S.GRID_C, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(f"({chr(97 + n + j)}) {'Error' if j == 0 else 'Bias'} at each survey", loc="left", fontsize=9,
                     fontweight="bold")
    if SURVEY_PANELS:
        h, lab = ax.get_legend_handles_labels()
        fig.legend(h, lab, loc="upper center", ncol=n, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 0.035))
    fig.suptitle("Modelled against surveyed dune crests with and without the collision module, 1996–2009",
                 fontsize=11, fontweight="bold", y=0.975 if SURVEY_PANELS else 0.985)
    out = OUT / OUT_NAME
    S.save(fig, out, bbox_inches="tight", close=True)
    (OUT / "supporting").mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(OUT / "supporting" / OUT_NAME.replace(".png", ".csv"), index=False)
    S.record_caption(out, (
        "Modelled against surveyed dune crests, 1996-2009, GIS 6-90 (the module is off at Cape Point). (a, b) Model "
        "minus lidar domain-median dune crest per GIS domain and lidar survey (red, model too high; blue, too low; "
        "white within ±0.5 m; grey, no survey), for the base (no collision module, dune growth 0.01-0.03), and the "
        "collision module (Larson et al. 2004 alone, Cs 1.9e-5, dune growth 0.02-0.04) with one toe per domain "
        "averaged over every 1996-2009 lidar survey (1998 left out). Right of each title: the crest RMSE and the "
        "share of domain-surveys within ±0.5 m over the whole run. No road rule, villages 40%. The 2001 survey is unclassified lidar, so its crests may include buildings and "
        "vegetation."))
    print(f"wrote {out}")
    print(t.pivot_table(index="survey", columns="run", values=["rmse_m", "bias_m"], sort=False).round(2).to_string())


if __name__ == "__main__":
    main()
