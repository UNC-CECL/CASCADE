#!/usr/bin/env python3
"""
The hindcast result figure: both periods' full-management runs against the LOWESS reference curve.

    python scripts/figure_making/model_output/hindcast_final_figure_lowess.py [--preset edgeBE]

Scores over D11-D89 and the canonical D2-D89 on the figure; reads the nogroin
matrix runs of the 1996 -> 2010 -> 2024 chain. Writes
output/comparisons/hindcast_calibrated/ and output/figures/5-results/hindcast_<preset>.png.
Details: scripts/figure_making/model_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import pathlib
import sys

import numpy as np
import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
import sys as _hssys
from pathlib import Path as _HSP
_hssys.path.insert(0, str(next(_q for _q in _HSP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_figure_style as _hs  # noqa: E402
OUT_DIR = _hs.COMPARISONS_ROOT / "hindcast_calibrated"
# The manuscript copy goes to output/figures/5-results/
FIGURES_DIR = _hs.figure_dir("results")
# On since the layout was redrawn for the house-style column (2026-09-18)
PUBLISH = True

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

# House style (site_layer/hat_figure_style.py), applied at import
from site_layer.hat_figure_style import (apply_style, figsize,  # noqa: E402
                              FIG_W_DOUBLE, record_caption, save,
                              C_1984, C_1997, INK_MUTED, DOMAIN_AXIS_LABEL)
apply_style()
from cascade_pipeline.run_layout import resolve                 # noqa: E402
from cascade_pipeline.run_registry import find_run_dir          # noqa: E402
from site_layer.hat_observed_rates import lrr_csv                # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_PERIODS     # noqa: E402

# --- CONFIG ------------------------------------------------------------------
LOWESS_PATH = (PROJECT_BASE_DIR / "scripts" / "input_prep" / "7-source-sink"
              / "2-calibrate" / "be_zone_residual_fit.py")

# The canonical chain 1996 -> 2010 -> 2024; nogroin, as the matrix is
GROIN = "nogroin"
_CHAIN = ((1996, "road_bdm", "a", C_1984), (2010, "road_bdm_nourish", "b", C_1997))
PERIODS = {}
for _st, _scen, _panel, _colour in _CHAIN:
    _end = HATTERAS_PERIODS[_st]["end_year"]
    PERIODS[f"{_st}_{_end}"] = dict(panel=_panel, label=f"{_st}–{_end}",
                                    start=_st, end=_end, scenario=_scen,
                                    colour=_colour)
LOCKED = (1, 90)
SKIP_SOUTH = 10               # matches LOWESS_CONFIG.skip_southern_domains
SCORE_DOMAINS = range(SKIP_SOUTH + 1, 90)
RESERVED_COLOUR = "#FF8C00"

# Non-overlapping place-label spans (D9-D10 go to Cape Point, as assign_physical_zone does)
PLACE_LABELS = [
    (1, 10, "Cape Point"), (11, 20, "Buxton"), (21, 31, "Avon"),
    (32, 59, "Mid-island"), (60, 74, "Wimble Shoals"),
    (75, 83, "Rodanthe"), (84, 90, "Pea Island"),
]
# -----------------------------------------------------------------------------


# Load be_zone_residual_fit.py by path, for its LOWESS and observation loaders
def analysis_module():
    spec = importlib.util.spec_from_file_location("_lowess", LOWESS_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# One period's LOWESS curve, transect SD by domain, and the transects
def lowess_and_spread(module, start, csv_path):
    from cascade_pipeline.coastsat_lowess import (CoastSatDataset,
                                                 build_coastsat_series)
    series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {start}", period_start=start,
                         csv_path=csv_path)], start, module.LOWESS_CONFIG)
    cs = series[0]
    match = [w for w in cs["windows"] if w["window"] == module.TARGET_WINDOW]
    if not match:
        raise ValueError(f"window {module.TARGET_WINDOW} not computed")
    lowess = dict(zip(np.asarray(match[0]["gis_x"], dtype=int),
                     np.asarray(match[0]["smoothed"], dtype=float)))
    frame = pd.read_csv(csv_path, usecols=["transect_id", "domain_number",
                                           "lrr_m_yr"])
    frame = frame.dropna(subset=["domain_number", "lrr_m_yr"])
    frame["domain"] = frame["domain_number"].astype(int)
    frame["lrr"] = frame["lrr_m_yr"].astype(float)

    # Spread each domain's transects across it in transect order, not stacked at the integer
    order = frame["transect_id"].astype(str).str.split("_").str[-1]
    frame["t_index"] = pd.to_numeric(order, errors="coerce")
    frame = frame.sort_values(["domain", "t_index"])
    rank = frame.groupby("domain").cumcount()
    count = frame.groupby("domain")["lrr"].transform("size")
    frame["x"] = frame["domain"] - 0.5 + (rank + 0.5) / count

    sd = frame.groupby("domain")["lrr"].std().to_dict()
    return lowess, {int(k): float(v) for k, v in sd.items()}, frame


# One run's per-domain rates and its name, resolved through run_registry
def run_rates(period, preset, scenario):
    # Matrix runs carry the offset token since the metres fix (2026-09-24)
    name = f"HAT_{period}_{preset}_offsetmetres_{scenario}_{GROIN}"
    # Resolved through run_registry, never joined by hand
    try:
        run_dir = find_run_dir(RAW_RUNS, name, period, preset)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"{exc}  Run HAT_run_all.py first.") from None
    path = resolve(run_dir, "rate_csv", name)
    return pd.read_csv(path).set_index("gis_domain")["lrr_m_yr"], name


# The preset names the file, title and footnote, so presets never overwrite each other
PRESETS = {
    "calibBE": dict(
        stem="hindcast_calibrated",
        published="hindcast_calibrated.png",
        title="Calibrated CASCADE hindcast against the CoastSat LOWESS reference, Cape Hatteras",
        series="CASCADE, calibrated",
        config="calibrated source/sink field, full management (roadway + beach/dune; "
               "nourishment in 2004–2024), groin active at M = 60 m/yr, f = 0.6.",
        zone_note="Source/sink corrections were applied only within a zone set fixed "
                  "before calibration; D5–D7 were reserved for the groin, so the misfit "
                  "there is the groin module's and is not absorbed by the sediment budget."),
    "edgeBE": dict(
        stem="hindcast_edgeBE",
        published="hindcast_edgeBE.png",
        title="CASCADE hindcast with the source/sink calibration removed, against the CoastSat LOWESS reference, Cape Hatteras",
        series="CASCADE edgeBE",
        config="edgeBE source/sink: the D1 and D90 edge values only, no interior "
               "correction. Full management (roadway + beach/dune; nourishment "
               "in 2010–2024); the Buxton groin is off, the nogroin matrix arm.",
        zone_note="NO interior source/sink correction was applied anywhere, so the frozen "
                  "zone set does not divide these domains; the edge values are kept because "
                  "they are solved by buffer-cell reproduction rather than fitted to a "
                  "residual. The misfit is the job a source/sink field would do."),
    "zeroBE": dict(
        stem="hindcast_zeroBE",
        published="hindcast_zeroBE.png",
        title="CASCADE hindcast with no source/sink field, against the CoastSat LOWESS reference, Cape Hatteras",
        series="CASCADE zeroBE",
        config="zeroBE: no background-erosion field anywhere, the D1/D90 edges included. "
               "Full management (roadway + beach/dune; nourishment in 2010–2024); the "
               "Buxton groin is off, the nogroin matrix arm.",
        zone_note="No source/sink field was applied anywhere, edges included, so the edges "
                  "carry no absorber and are free to run away from the target."),
}


# Run: load both periods, score both windows, draw, save both copies
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    # calibBE is not solved on the 1996/2010 chain, so edgeBE is the default.
    parser.add_argument("--preset", default="edgeBE", choices=sorted(PRESETS))
    args = parser.parse_args()
    vocab = PRESETS[args.preset]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    import dataclasses

    from cascade_pipeline.annotations import (add_geographic_annotations,
                                              annotation_legend_handles)
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS

    # Label heights tuned for this figure's y range and layout (README)
    annotations = dataclasses.replace(
        HATTERAS_ANNOTATIONS,
        groin_label_y=0.78,
        piers={"Avon Pier": (26, 0.85), "Rodanthe Pier": (79, 0.72)})

    module = analysis_module()
    reserved = list(module.GROIN_RESERVED_DOMAINS)

    panels = {}
    for key, meta in PERIODS.items():
        csv = str(lrr_csv(meta["start"], meta["end"]))
        lowess, sd, transects = lowess_and_spread(module, meta["start"], csv)
        model, run_name = run_rates(key, args.preset, meta["scenario"])
        # The D2-D89 score is computed here, never quoted, so it cannot go stale
        spliced = module.load_observed(meta["start"], csv)[1]
        wide = [g for g in range(2, 90)
                if g in model.index and not np.isnan(spliced.get(g, np.nan))]
        wm = np.array([model[g] for g in wide])
        wo = np.array([spliced[g] for g in wide])
        companion = dict(rmse=float(np.sqrt(np.mean((wm - wo) ** 2))),
                         bias=float(np.mean(wm - wo)),
                         corr=float(np.corrcoef(wm, wo)[0, 1]),
                         n=len(wide))
        panels[key] = dict(lowess=lowess, sd=sd, transects=transects,
                           model=model, run=run_name, companion=companion,
                           **meta)

    gis = np.arange(1, 91)
    lo, hi = np.inf, -np.inf
    for d in panels.values():
        for source in (d["lowess"], d["model"]):
            values = np.array([source.get(g, np.nan) for g in gis], dtype=float)
            lo = min(lo, np.nanmin(values))
            hi = max(hi, np.nanmax(values))
    # Widen the axis for the D1-D10 scatter too, without extra label headroom
    s_lo, s_hi = np.inf, -np.inf
    for d in panels.values():
        south_t = d["transects"][d["transects"]["domain"] <= SKIP_SOUTH]["lrr"]
        if len(south_t):
            s_lo = min(s_lo, float(south_t.min()))
            s_hi = max(s_hi, float(south_t.max()))
    pad = 0.10 * (hi - lo)
    # Asymmetric padding: only the top needs room for place names
    ylim = (min(lo - pad * 0.45, s_lo - pad * 0.1),
            max(hi + pad * 1.5, s_hi + pad * 0.1))

    # House-style layout: titles at the panels, place names once, title and note in the caption
    figure, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.0),
                                sharex=True, sharey=True)

    summary = []
    for i, (axis, (key, d)) in enumerate(zip(axes, panels.items())):
        obs = np.array([d["lowess"].get(g, np.nan) for g in gis], dtype=float)
        mod = np.array([d["model"].get(g, np.nan) for g in gis], dtype=float)
        spread = np.array([d["sd"].get(g, np.nan) for g in gis], dtype=float)

        # Place layer from the shared annotations, named on (a) only
        add_geographic_annotations(axis, annotations, label=(i == 0))

        # Reserved and locked spans are stated in the caption, not drawn
        axis.axvspan(0.5, SKIP_SOUTH + 0.5, facecolor="none",
                     edgecolor="#999999", hatch="\\\\\\", linewidth=0.0,
                     alpha=0.30, zorder=0)

        axis.fill_between(gis, obs - spread, obs + spread, color="#999999",
                          alpha=0.22, zorder=2, linewidth=0,
                          label="observed spread (±1 SD of transects)")

        # D1-D10 transects drawn; the LOWESS is dashed there, where the smoother is unreliable
        south_t = d["transects"][d["transects"]["domain"] <= SKIP_SOUTH]
        axis.plot(south_t["x"], south_t["lrr"], linestyle="none",
                  marker="o", markersize=1.8, color=INK_MUTED, alpha=0.6,
                  zorder=5, label="CoastSat transect LRR (D1–D10)")
        axis.fill_between(gis, obs, mod, where=~(np.isnan(obs) | np.isnan(mod)),
                          color=d["colour"], alpha=0.18, zorder=3, linewidth=0,
                          label=f"misfit, {d['label']}")

        south = gis <= SKIP_SOUTH
        axis.plot(gis[~south], obs[~south], color="#1A1A1A", linewidth=1.8,
                  zorder=7, label="CoastSat LOWESS (7-domain)")
        axis.plot(gis[south], obs[south], color="#1A1A1A", linewidth=1.3,
                  linestyle=(0, (4, 2)), zorder=7,
                  label="LOWESS, D1–D10 (excluded by convention)")
        axis.plot(gis, mod, color=d["colour"], linewidth=1.6, zorder=6,
                  label=f"{vocab['series']}, {d['label']}")
        axis.axhline(0.0, color="#AAAAAA", linewidth=0.6, zorder=4)

        shared = [g for g in SCORE_DOMAINS
                  if g in d["model"].index and not np.isnan(d["lowess"].get(g, np.nan))]
        mm = np.array([d["model"][g] for g in shared])
        oo = np.array([d["lowess"][g] for g in shared])
        rmse = float(np.sqrt(np.mean((mm - oo) ** 2)))
        bias = float(np.mean(mm - oo))
        corr = float(np.corrcoef(mm, oo)[0, 1])
        summary.append((d["label"], d["run"], rmse, bias, corr, len(shared)))

        # Both scoring windows on the figure: D11-D89 (the LOWESS target) and D2-D89 (canonical)
        c = d["companion"]
        axis.set_title(f"({d['panel']})  {d['label']}", loc="left",
                       fontweight="bold")
        axis.set_title(
            f"D{SKIP_SOUTH + 1}–D89 (n = {len(shared)}): RMSE {rmse:.2f} m/yr, "
            f"bias {bias:+.2f}, r = {corr:.2f}\n"
            f"D2–D89 (n = {c['n']}): RMSE {c['rmse']:.2f} m/yr, "
            f"bias {c['bias']:+.2f}, r = {c['corr']:.2f}",
            loc="right", fontsize=7.5, color=INK_MUTED, linespacing=1.3)
        axis.set_ylim(ylim)

    figure.supylabel("Shoreline change rate (m/yr, + seaward)", fontsize=9)
    axes[1].set_xlabel(DOMAIN_AXIS_LABEL)
    axes[1].set_xlim(0, 91)
    axes[1].set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])

    # One key for both panels, in reading order
    found = {}
    for axis in axes:
        for h, l in zip(*axis.get_legend_handles_labels()):
            found.setdefault(l, h)
    order = ["CoastSat LOWESS (7-domain)",
             "LOWESS, D1–D10 (excluded by convention)",
             "CoastSat transect LRR (D1–D10)",
             "observed spread (±1 SD of transects)"]
    for d in panels.values():
        order += [f"{vocab['series']}, {d['label']}", f"misfit, {d['label']}"]
    labels = [l for l in order if l in found]
    handles = [found[l] for l in labels]
    extra = annotation_legend_handles(annotations)
    handles += extra
    labels += [h.get_label() for h in extra]
    # Legend at figure level, below both panels
    figure.legend(handles, labels, loc="lower center", ncol=3, fontsize=7.5,
                  frameon=False, bbox_to_anchor=(0.5, 0.0))
    _leg_rows = math.ceil(len(handles) / 3)
    figure.tight_layout(rect=(0.01, 0.026 * _leg_rows + 0.01, 1, 1))

    caption_text = (
        vocab["title"] + ". Configuration: " + vocab["config"]
        + " The observed curve is the 7-domain LOWESS of CoastSat transect "
        "rates; over D1–D10 (hatched) it is dashed because the project "
        "excludes the LOWESS there (the smoother is poorly constrained at the "
        "end of its range, and Cape Point's attachment-detachment cycle is "
        "short-wavelength signal a 3.5 km smoother removes), and the individual "
        "transect rates are drawn instead. The grey band is ±1 SD of the "
        "transect rates within each domain, and the tinted band between the "
        "curves is the misfit. Both scoring windows are printed above each "
        "panel. D11–D89 is the span where the LOWESS curve and the grading "
        "target are identical numbers, so that statistic describes the line "
        "actually drawn. D2–D89 is the project's canonical skill window "
        "(rmse_interior_m_yr in run_index.csv), which also takes in D2–D10, "
        "where the target is the raw spliced domain mean and the model is at "
        "its worst; the gap between the two numbers is that strip alone, not "
        "the smoothing. D1 and D90 are in neither: they are locked boundary "
        "absorbers, solved by buffer-cell reproduction, and would dominate. "
        + vocab["zone_note"])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"{vocab['stem']}_lowess_reference.png"
    figure.savefig(path)
    print(f"wrote {path}")
    if PUBLISH:
        published = FIGURES_DIR / vocab["published"]
        save(figure, published)
        record_caption(published, caption_text)
        print(f"wrote {published}")
    plt.close(figure)

    print(f"  shared y limits: {ylim[0]:.2f} to {ylim[1]:.2f} m/yr")
    for label, run, rmse, bias, corr, n in summary:
        print(f"  {label}  {run}")
        print(f"      RMSE {rmse:.4f}   bias {bias:+.4f}   r {corr:.3f}   n {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
