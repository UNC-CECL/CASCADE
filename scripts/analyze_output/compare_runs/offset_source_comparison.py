"""
==============================================================================
offset_source_comparison.py -- how much does the island offset's source
(dune line or CoastSat shoreline) change what the model does?
==============================================================================
Asked by Hannah on 2026-09-28: a simple comparison of the model output started
from the dune-line offset against the model output started from the shoreline
offset, for 1996-2010 and 2010-2024. The main question is how much the island's
ORIENTATION in the offset affects the outcome.

RUNS  (no new runs; the full-management pair of each period from)
    output/raw_runs/experiments/island-offset/
        2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/runs/
            {duneline,shoreline}_full_management/<period>/zeroBE/<run>/
    option A waves, zeroBE ends, relocations and groins off. The only thing
    that differs within a pair is the island offset.

OFFSETS  2-brie-offset/<start>/{duneline,shoreline}/<version>/*_unpadded.csv,
    at the version each run's metadata records. BRIE adds the offset to x_s,
    so larger = more LANDWARD; everything below uses the SEAWARD position
    s = -offset, with the mean removed (a uniform shift does not change what
    BRIE does, and the builds are not on a common datum).

QUANTITIES  per GIS domain (500 m), shoreline start minus dune-line start
    orientation   theta = atan(ds/dx), degrees
    turning       d(theta)/dx, degrees per km; positive where the line
                  bends landward (an embayment), negative at a bulge
    model change  each run's LRR x 14 yr (m, seaward positive)
    Statistics on the interior, GIS 2-89.

OUTPUT   output/comparisons/offset_source/
    offset_source_model_change_full_management.png  the two runs, both periods
    offset_source_model_change_vs_projected_full_management.png
                    the same, with the projected target on top (2026-09-29)
    offset_source_difference_full_management.png   profiles
    offset_source_orientation_vs_model_full_management.png   scatter
    tables/summary.csv, tables/per_domain.csv, tables/vs_projected.csv

TARGET  (second version of the change figure, Hannah 2026-09-29) projected
    shoreline change: the CoastSat LRR fitted on 1996-2024, LOWESS over 7
    domains (southern 10 raw), x 14 yr -- one profile, the same in both
    periods. The model stays unsmoothed.

USAGE
    python scripts/analyze_output/compare_runs/offset_source_comparison.py
==============================================================================
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_REPO = next(p for p in _HERE.parents if (p / "scripts").is_dir() and (p / "data").is_dir())
sys.path.insert(0, str(_REPO / "scripts"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,
    figsize, open_frame, record_caption, save, structures, town_bands,
)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402
from site_layer.hat_topo_version import offset_file  # noqa: E402

STUDY = (_REPO / "output" / "raw_runs" / "experiments" / "island-offset"
         / "2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a")
OUT = COMPARISONS_ROOT / "offset_source"
PERIODS = ((1996, 2010), (2010, 2024))
SOURCES = ("duneline", "shoreline")
COL = {"duneline": "#1b7f6b", "shoreline": "#6a3d9a"}   # as in the study's figures
DX_M = 500.0              # BarrierLength: one GIS domain
YEARS = 14
INTERIOR = (2, 89)
LOWESS_DOMAINS = 7         # the group's smoothing range
SKIP_SOUTHERN = 10


def projected_target():
    """Projected shoreline change (m): the 1996-2024 CoastSat LRR target, built
    as the runner builds it at LOWESS_DOMAINS, x YEARS."""
    from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.coastsat_lowess import (CoastSatDataset, LowessConfig,
                                                 build_coastsat_series)
    ds = CoastSatDataset(label="CoastSat LRR (1996-2024)", period_start=1996,
                         csv_path=str(COASTSAT_LRR_ROOT / "1996_2024" / "transect_lrr_full.csv"))
    cfg = LowessConfig(window_domains=(LOWESS_DOMAINS,), skip_southern_domains=SKIP_SOUTHERN)
    cs = build_coastsat_series([ds], active_period_start=1996, lowess_config=cfg,
                               domains=HATTERAS_DOMAINS)[0]
    return build_target_table(cs, cfg, HATTERAS_DOMAINS, LOWESS_DOMAINS).set_index(
        "gis_domain")["target_lrr_m_yr"] * YEARS


def vs_target(t, obs):
    """Each run against the projected target, interior GIS 2-89: bias and RMS
    residual are the numbers to read; explained and r beside them."""
    rows = []
    for start, end in PERIODS:
        p = t[t.period == f"{start}_{end}"].set_index("gis_domain").loc[INTERIOR[0]:INTERIOR[1]]
        o = obs.reindex(p.index)
        for src in SOURCES:
            m = p[f"model_change_{src}_m"]
            res = m - o
            rows.append(dict(period=f"{start}-{end}", offset=src, bias_m=res.mean(),
                             rms_residual_m=np.sqrt((res ** 2).mean()),
                             variance_explained=1 - (res ** 2).sum() / ((o - o.mean()) ** 2).sum(),
                             r=np.corrcoef(m, o)[0, 1]))
    return pd.DataFrame(rows)


def run_dir(src, start, end):
    root = STUDY / "runs" / f"{src}_full_management" / f"{start}_{end}" / "zeroBE"
    hits = [d for d in root.glob("*") if d.is_dir()]
    if len(hits) != 1:
        raise SystemExit(f"{root}: expected one run, found {len(hits)}")
    return hits[0]


def offset_version(d, src):
    md = json.loads(next(d.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    got = str(md["identity"]["island_offset_version"])       # e.g. "shoreline/v1"
    if not got.startswith(src + "/"):
        raise ValueError(f"{d}: ran on offset {got!r}, filed as {src}")
    return got.split("/", 1)[1]


def seaward(src, start, version):
    t = pd.read_csv(offset_file(start, "unpadded", version=version, source=src))
    s = -t.set_index(t.columns[0])[t.columns[1]].astype(float)
    s.index = s.index.astype(int)
    return s - s.mean()


def orientation(s):
    theta = np.degrees(np.arctan(np.gradient(s.to_numpy(), DX_M)))
    turning = np.gradient(theta, DX_M / 1000.0)
    return pd.Series(theta, s.index), pd.Series(turning, s.index)


def build():
    per, rows = [], []
    for start, end in PERIODS:
        cols = {}
        for src in SOURCES:
            d = run_dir(src, start, end)
            s = seaward(src, start, offset_version(d, src))
            th, tu = orientation(s)
            rt = pd.read_csv(d / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
            cols.update({f"offset_{src}_m": s, f"orientation_{src}_deg": th,
                         f"turning_{src}_deg_km": tu, f"model_change_{src}_m": rt.lrr_m_yr * YEARS})
        t = pd.DataFrame(cols)
        t.index.name = "gis_domain"
        for q in ("offset", "orientation", "turning", "model_change"):
            unit = {"offset": "m", "orientation": "deg", "turning": "deg_km",
                    "model_change": "m"}[q]
            t[f"d_{q}_{unit}"] = t[f"{q}_shoreline_{unit}"] - t[f"{q}_duneline_{unit}"]
        t.insert(0, "period", f"{start}_{end}")
        per.append(t.reset_index())

        i = t.loc[INTERIOR[0]:INTERIOR[1]]
        dm = i.d_model_change_m
        rows.append(dict(
            period=f"{start}-{end}",
            model_sd_duneline_m=i.model_change_duneline_m.std(),
            model_sd_shoreline_m=i.model_change_shoreline_m.std(),
            r_between_runs=np.corrcoef(i.model_change_duneline_m, i.model_change_shoreline_m)[0, 1],
            mean_diff_m=dm.mean(), mean_abs_diff_m=dm.abs().mean(), sd_diff_m=dm.std(),
            max_abs_diff_m=dm.abs().max(), max_at_gis=int(dm.abs().idxmax()),
            diff_sd_over_model_sd=dm.std() / i.model_change_duneline_m.std(),
            sd_d_offset_m=i.d_offset_m.std(),
            sd_d_orientation_deg=i.d_orientation_deg.std(),
            max_abs_d_orientation_deg=i.d_orientation_deg.abs().max(),
            r_d_orientation_vs_d_model=np.corrcoef(i.d_orientation_deg, dm)[0, 1],
            r_d_turning_vs_d_model=np.corrcoef(i.d_turning_deg_km, dm)[0, 1],
            r_d_offset_vs_d_model=np.corrcoef(i.d_offset_m, dm)[0, 1],
            duneline_run=str(run_dir("duneline", start, end).relative_to(_REPO)).replace("\\", "/"),
            shoreline_run=str(run_dir("shoreline", start, end).relative_to(_REPO)).replace("\\", "/"),
        ))
    return pd.concat(per, ignore_index=True), pd.DataFrame(rows)


COMMON = (" Full management, option A waves (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle "
          "0.5), no source/sink correction at the ends (zeroBE), relocations and groins off; "
          "within each period the two runs differ ONLY in the island offset (dune-line or "
          "mean CoastSat shoreline build, v1 of each start: 1997 dune line / 1995-1997 "
          "shoreline for 1996, 2009 dune line / 2009-2011 shoreline for 2010). Domain 1 is "
          "Cape Point, 90 Pea Island; 500 m per domain. Runs from "
          "raw_runs/experiments/island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-"
          "waves-option-a. Statistics in tables/summary.csv, interior GIS 2-89.")


def fig_profiles(t):
    f, axes = plt.subplots(3, 2, figsize=figsize("double", height=6.6), sharex=True,
                           sharey="row", constrained_layout=True)
    for j, (start, end) in enumerate(PERIODS):
        p = t[t.period == f"{start}_{end}"].set_index("gis_domain")
        a0, a1, a2 = axes[:, j]
        for src in SOURCES:
            a0.plot(p.index, p[f"model_change_{src}_m"], color=COL[src], lw=1.2)
        a1.bar(p.index, p.d_orientation_deg, width=0.85, color=INK_MUTED)
        a2.bar(p.index, p.d_model_change_m, width=0.85,
               color=[COL["shoreline"] if v > 0 else COL["duneline"] for v in p.d_model_change_m])
        for k, ax in enumerate((a0, a1, a2)):
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(0.5, 90.5)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(k == 0))
            _title(ax, 2 * k + j, f"{start}–{end}")
        a2.set_xlabel(DOMAIN_AXIS_LABEL)
        if j == 0:
            a0.set_ylabel("Modelled change (m)")
            a1.set_ylabel("Orientation difference,\nshoreline − dune line (°)")
            a2.set_ylabel("Change difference,\nshoreline − dune line (m)")
    f.legend(handles=[Line2D([], [], color=COL["duneline"], lw=1.2, label="Started from the dune line"),
                      Line2D([], [], color=COL["shoreline"], lw=1.2, label="Started from the shoreline"),
                      Patch(color=COL["shoreline"], label="Shoreline start more accretional"),
                      Patch(color=COL["duneline"], label="Dune-line start more accretional")],
             loc="outside lower center", ncol=2, frameon=False)
    png = OUT / "offset_source_difference_full_management.png"
    save(f, png, close=True)
    record_caption(png, (
        "Model outcome with the island offset taken from the dune line or from the shoreline. "
        "(a, b) Modelled total shoreline change, each run's LRR x 14 yr, seaward positive: "
        "green started from the dune line, purple from the shoreline. (c, d) The difference in "
        "the offset's orientation, shoreline minus dune line, atan of the alongshore slope of "
        "the seaward position, degrees. (e, f) The difference in modelled change, shoreline "
        "start minus dune-line start; purple where the shoreline start is more accretional." + COMMON))
    return png


def fig_change_only(t, obs=None, scores=None):
    """Panels (a, b) of the profile figure on their own: the two runs' change.
    With `obs`, the second version: the projected target drawn on top of them.
    Every label stays out of the data (Hannah, 2026-09-28): the villages and
    shoals are named in a strip above the highest line, and the groin and
    piers are drawn below that strip and named in the legend, not on the lines."""
    # Draw order, bottom to top: village and shoal bands (0-0.5), grid (below
    # everything, set_axisbelow), zero line (2), groin and piers (3), the two
    # runs (4), labels (7). The structure lines stop at DATA_TOP, under the
    # label strip, and the grid has no ticks inside it.
    DATA_TOP = 0.78       # axes fraction the highest line may reach
    f, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.1), sharey=True,
                           constrained_layout=True)
    vals = t[[f"model_change_{s}_m" for s in SOURCES]].to_numpy().ravel()
    if obs is not None:
        vals = np.concatenate([vals, obs.to_numpy()])
    lo, hi = np.nanmin(vals), np.nanmax(vals)
    y0 = lo - 0.04 * (hi - lo)
    ylim = (y0, y0 + (hi - y0) / DATA_TOP)
    for j, (start, end) in enumerate(PERIODS):
        p = t[t.period == f"{start}_{end}"].set_index("gis_domain")
        ax = axes[j]
        for src in SOURCES:
            ax.plot(p.index, p[f"model_change_{src}_m"], color=COL[src], lw=1.2, zorder=4)
        if obs is not None:
            ax.plot(obs.index, obs.values, color=INK, lw=1.8, zorder=5)
        ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
        ax.set_xlim(0.5, 90.5)
        ax.set_ylim(*ylim)
        ax.set_axisbelow(True)                            # grid under everything drawn
        ax.grid(axis="y")
        # no gridlines through the label strip
        ax.set_yticks([v for v in ax.get_yticks() if ylim[0] <= v <= hi])
        open_frame(ax)
        town_bands(ax)                                    # names at the top row
        for nm, (a_, b_) in HATTERAS_ANNOTATIONS.shoal_zones.items():
            ax.axvspan(a_ - 0.5, b_ + 0.5, color=C["ADDED"], alpha=0.12, lw=0, zorder=0.5)
            ax.text((a_ + b_) / 2, 0.905, nm, transform=ax.get_xaxis_transform(),
                    ha="center", va="top", fontsize=6, color="#8a620e", zorder=7)
        for pos in HATTERAS_ANNOTATIONS.groins.values():
            ax.axvline(pos, ymax=DATA_TOP, color=INK, lw=0.7, zorder=3)
        for pos, _frac in HATTERAS_ANNOTATIONS.piers.values():
            ax.axvline(pos, ymax=DATA_TOP, color=INK_MUTED, lw=0.6,
                       ls=(0, (1.5, 1.5)), zorder=3)
        _title(ax, j, f"{start}–{end}")
        ax.set_xlabel(DOMAIN_AXIS_LABEL)
    axes[0].set_ylabel("Modelled change (m)" if obs is None else "Shoreline change (m)")
    target = ([Line2D([], [], color=INK, lw=1.8,
                      label="Projected CoastSat position")]
              if obs is not None else [])
    f.legend(handles=target + [Line2D([], [], color=COL["duneline"], lw=1.2, label="Started from the dune line"),
                      Line2D([], [], color=COL["shoreline"], lw=1.2, label="Started from the shoreline"),
                      Line2D([], [], color=INK, lw=0.7, label="Buxton groin"),
                      Line2D([], [], color=INK_MUTED, lw=0.6, ls=(0, (1.5, 1.5)),
                             label="Avon and Rodanthe piers"),
                      Patch(color=C["ADDED"], alpha=0.25, label="Shoals"),
                      Patch(color="0.94", label="Villages")],
             loc="outside lower center", ncol=3, frameon=False)
    if obs is not None:
        png = OUT / "offset_source_model_change_vs_projected_full_management.png"
        save(f, png, close=True)
        sc = "; ".join(
            f"{r.period} {'dune line' if r.offset == 'duneline' else 'shoreline'} bias "
            f"{r.bias_m:+.1f} m, RMS residual {r.rms_residual_m:.1f} m"
            for r in scores.itertuples())
        record_caption(png, (
            "Modelled total shoreline change, each run's LRR x 14 yr, seaward positive, with the "
            "island offset taken from the dune line (green) or from the shoreline (purple), "
            "against PROJECTED shoreline change (black): the CoastSat LRR fitted on 1996-2024, "
            "LOWESS over 7 domains (southern 10 raw), x 14 yr -- the same profile in both "
            "panels, carried onto each period rather than fitted on it; the model is "
            "unsmoothed. (a) 1996-2010, (b) 2010-2024. The version of "
            "offset_source_model_change_full_management.png with the target on top. Interior "
            "GIS 2-89, model minus target: " + sc + " (tables/vs_projected.csv). Amber: Avon "
            "and Wimble Shoals; solid line: Buxton groin; dotted lines: the Avon (GIS 26) and "
            "Rodanthe (GIS 79) piers; grey bands: villages." + COMMON))
        return png
    png = OUT / "offset_source_model_change_full_management.png"
    save(f, png, close=True)
    record_caption(png, (
        "Modelled total shoreline change, each run's LRR x 14 yr, seaward positive, with the "
        "island offset taken from the dune line (green) or from the shoreline (purple): "
        "(a) 1996-2010, (b) 2010-2024. Amber: Avon and Wimble Shoals; solid line: Buxton "
        "groin; dotted lines: the Avon (GIS 26) and Rodanthe (GIS 79) piers; grey bands: villages. Panels (a, b) of "
        "offset_source_difference_full_management.png on their own." + COMMON))
    return png


def fig_scatter(t, summary):
    f, axes = plt.subplots(2, 2, figsize=figsize("double", height=5.4), sharey=True,
                           sharex="col", constrained_layout=True)
    for j, (start, end) in enumerate(PERIODS):
        p = t[t.period == f"{start}_{end}"].set_index("gis_domain").loc[INTERIOR[0]:INTERIOR[1]]
        c = "#b2182b" if j == 0 else "#2166ac"   # house vintages: earlier red, later blue
        for k, (col, lab) in enumerate((("d_orientation_deg", "Orientation difference (°)"),
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
                ax.set_xlabel(lab + ", shoreline − dune line")
            if k == 0:
                ax.set_ylabel("Change difference (m)")
    png = OUT / "offset_source_orientation_vs_model_full_management.png"
    save(f, png, close=True)
    s = summary.set_index("period")
    record_caption(png, (
        "Which property of the offset drives the difference between the two runs. Each point "
        "is one interior domain (GIS 2-89); y is the modelled change with the shoreline start "
        "minus that with the dune-line start (m). Left (a, c): against the difference in "
        "orientation (the angle itself). Right (b, d): against the difference in turning, the "
        "alongshore rate of change of that angle (degrees per km; positive where the line "
        "bends landward, an embayment, negative at a bulge). Line: least squares. The angle "
        "difference predicts almost nothing (r {:+.2f}, {:+.2f}); the turning difference "
        "predicts most of it (r {:+.2f}, {:+.2f}): where one offset makes a local embayment "
        "the other does not, that run fills it in, where it makes a bulge, that run cuts it "
        "back, as BRIE's alongshore transport diffuses the planform.".format(
            s.loc["1996-2010", "r_d_orientation_vs_d_model"],
            s.loc["2010-2024", "r_d_orientation_vs_d_model"],
            s.loc["1996-2010", "r_d_turning_vs_d_model"],
            s.loc["2010-2024", "r_d_turning_vs_d_model"]) + COMMON))
    return png


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    apply_style()
    t, summary = build()
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / "tables" / "per_domain.csv", index=False)
    summary.to_csv(OUT / "tables" / "summary.csv", index=False)
    obs = projected_target()
    scores = vs_target(t, obs)
    scores.to_csv(OUT / "tables" / "vs_projected.csv", index=False)
    for png in (fig_profiles(t), fig_change_only(t), fig_change_only(t, obs, scores),
                fig_scatter(t, summary)):
        print(png.relative_to(_REPO))
    print(summary.drop(columns=["duneline_run", "shoreline_run"]).round(2).T.to_string())
    print(scores.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
