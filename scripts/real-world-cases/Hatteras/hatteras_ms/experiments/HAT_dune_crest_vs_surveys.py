"""
How well do the modelled dunes follow the surveyed ones through the calibration and test windows?

    python scripts/hatteras_ms/experiments/HAT_dune_crest_vs_surveys.py

Compares the dune crests of the runs from the ceiling x growth grid (the
winner, 5.0 m and r 0.25-0.45, and the control, 5.5 m and r 0.55-0.95, both
with the groin on) against the surveys flown inside each window, taken on that
start's dune windows by
scripts/input_prep/1-barrier3d-domains/2-survey-crests/HAT_survey_dune_crests.py:
1996-2009 against the 1997-2005 lidar (coarse 3-9 m grids, drawn apart) and the
2009 lidar; 2009-2025 against the 2014-2019 lidar.
Each survey is set against the model state on the 1 January its metadata
dates point to. A survey cell whose 2009 dune position is now beach counts as
no dune (the berm, 1.7 m NAVD88). Writes tables and a figure. Read-only.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
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
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_ceiling_growth_grid as G  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-19-dune-crest-vs-surveys"
EXP_DIR = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
# start year -> (topography product, first saved state: DuneDomain[0] dated 1 January of)
PERIODS = {1996: ("1996-start", 1997), 2009: ("2009-start", 2010)}
RUNS = {"Winner: 5.0 m, r 0.25-0.45": G.name(5.0, (0.25, 0.45), groin=True),
        "Control: 5.5 m, r 0.55-0.95": G.name(*G.CONTROL, groin=True)}
RUN_COLOURS = {"Winner: 5.0 m, r 0.25-0.45": fs.C["ACCENT"], "Control: 5.5 m, r 0.55-0.95": fs.INK_MUTED}
BERM_NAVD88 = 1.7
DESIGN_NAVD88 = 4.3
# -----------------------------------------------------------------------------


# Model crest per saved state, GIS 1-90 x 50 cells, m NAVD88: (states, 4500)
def model_crests(start, variant):
    c = U.load(G.run_dir(start, variant))
    cc = U.cell_crests(c, DOM)                                     # m above berm
    return np.stack([x for x in cc], axis=1).reshape(cc[0].shape[0], -1) + BERM_NAVD88


def surveys(start):
    d = tv.product_dir(PERIODS[start][0]) / "survey-crests"
    t = pd.read_csv(d / "survey_dune_crests.csv")
    t["crest"] = t.crest_m_navd88.where(t.status != "beach_landward_of_window", BERM_NAVD88)
    kinds = pd.read_csv(d / "survey_sources.csv").set_index("survey").kind if (d / "survey_sources.csv").exists() else {}
    t["kind"] = t.survey.map(lambda y: kinds.get(y, "mosaic"))
    return t


# 2009 start crests, the same cells, m NAVD88 (model year 0 equals them)
def start_crests(model):
    return model[0]


def compare(obs, mod):
    ok = np.isfinite(obs)
    o, m = obs[ok], mod[ok]
    o2, m2 = obs.reshape(90, 50), mod.reshape(90, 50)
    full = np.isfinite(o2).sum(1) >= 25                      # domains surveyed in at least half their cells
    od, md = np.nanmedian(o2[full], 1), np.median(m2[full], 1)
    return dict(survey_median=np.median(o), model_median=np.median(m),
                survey_p10=np.percentile(o, 10), model_p10=np.percentile(m, 10),
                survey_share_below_design=float((o < DESIGN_NAVD88).mean()),
                model_share_below_design=float((m < DESIGN_NAVD88).mean()),
                median_offset_m=float(np.median(m - o)), mean_abs_error_m=float(np.mean(np.abs(m - o))),
                domain_median_r=float(np.corrcoef(od, md)[0, 1]),
                domain_min_r=float(np.corrcoef(np.nanmin(o2[full], 1), m2[full].min(1))[0, 1]))


def main():
    for start, (prod, first) in PERIODS.items():
        s = surveys(start)
        models = {label: model_crests(start, v) for label, v in RUNS.items()}
        rows = []
        for label, mod in models.items():
            rows.append(dict(run=label, survey=f"{start} DEM (start)", kind="mosaic",
                             model_state=f"{first}-01-01", **compare(mod[0], mod[0])))
            for year, g in s.groupby("survey"):
                g = g.sort_values(["gis", "cell"])
                k = int(g.model_state.iloc[0][:4]) - first
                if k == 0:
                    continue                                       # the start itself: the bias check, not a target
                rows.append(dict(run=label, survey=f"{year} lidar", kind=g.kind.iloc[0],
                                 model_state=g.model_state.iloc[0], **compare(g.crest.to_numpy(), mod[k])))
        t = pd.DataFrame(rows).round(3)
        (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
        t.to_csv(EXP_DIR / "tables" / f"model_vs_survey_crests_{start}_start.csv", index=False)
        with pd.option_context("display.width", 250, "display.max_columns", 30):
            print(t.to_string(index=False))
        figure(start, first, s, models)


# Crest statistics through a window: model lines every year, surveys as points (coarse grids as squares)
def figure(start, first, s, models):
    fs.apply_style()
    fig, axes = plt.subplots(1, 3, figsize=fs.figsize("double", height=3.3), layout="constrained")
    yrs = {label: first + np.arange(m.shape[0]) for label, m in models.items()}
    sv = s[s.model_state.str[:4].astype(int) > first].groupby("survey")
    pts = pd.DataFrame({"year": [int(g.model_state.iloc[0][:4]) for _, g in sv],
                        "kind": [g.kind.iloc[0] for _, g in sv],
                        "median": [g.crest.median() for _, g in sv],
                        "p10": [g.crest.quantile(.1) for _, g in sv],
                        "below": [(g.crest < DESIGN_NAVD88).mean() * 100 for _, g in sv]})
    s0 = next(iter(models.values()))[0]
    pts = pd.concat([pd.DataFrame(dict(year=[first], kind=["mosaic"], median=[np.median(s0)],
                                       p10=[np.percentile(s0, 10)], below=[(s0 < DESIGN_NAVD88).mean() * 100])), pts])
    for ax, key, f, ylab in ((axes[0], "median", lambda a, axis: np.median(a, axis=axis), "Median crest (m NAVD88)"),
                             (axes[1], "p10", lambda a, axis: np.percentile(a, 10, axis=axis),
                              "Low dunes, 10th percentile (m NAVD88)"),
                             (axes[2], "below", lambda a, axis: (a < DESIGN_NAVD88).mean(axis=axis) * 100,
                              f"Cells below {DESIGN_NAVD88} m (%)")):
        for label, m in models.items():
            ax.plot(yrs[label], f(m, axis=1), color=RUN_COLOURS[label], lw=1.4, label=label)
        for kind, marker, lab in (("mosaic", "o", "Surveyed, 1 m DEM"), ("raster", "s", "Surveyed, coarse 3-9 m grid")):
            p = pts[pts.kind == kind]
            if len(p):
                ax.plot(p.year, p[key], ls="none", marker=marker, ms=5, mfc="white", mec=fs.INK, mew=1.2,
                        label=lab, zorder=5)
        ax.set_ylabel(ylab)
        ax.set_xlabel("1 January of")
        ax.grid(axis="y", color=fs.GRID_C, lw=0.5)
        fs.open_frame(ax)
    for i, ax in enumerate(axes):
        fs._letter_inside(ax, i)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=4, frameon=False, fontsize=7)
    end = U.PERIODS[start][0]
    out = EXP_DIR / "figures" / f"dune_crest_model_vs_surveys_{start}_{end}.png"
    fs.save(fig, out, close=True)
    fs.record_caption(out, (
        f"Dune crests through the {start}-{end} window, all 4,500 cells of 10 m along GIS 1-90: the model each "
        "1 January (lines; winner 5.0 m ceiling and growth 0.25-0.45, purple; control 5.5 m and 0.55-0.95, grey; "
        f"both with the groin) and the surveys flown inside the window, at the 1 January their flight dates point "
        f"to, all on the {start} dune windows (circles: 1 m DEMs; squares: the 1997-2005 lidar as downloaded, "
        "3-9 m grids, which read individual crests less precisely but match the 1 m 1996 survey in median, "
        "10th percentile and domain pattern, r 0.99). (a) median crest; (b) 10th percentile, the low dunes; "
        "(c) share of cells below the NC-12 design height. Survey cells whose starting dune position is now beach "
        "count as no dune (1.7 m)."))
    print(out)


if __name__ == "__main__":
    main()
