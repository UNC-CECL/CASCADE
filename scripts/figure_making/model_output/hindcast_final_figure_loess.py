#!/usr/bin/env python3
"""The calibrated hindcast against the LOESS reference curve — presentation figure.

THE hindcast result figure. It had a companion, HAT_hindcast_final_figure.py,
which drew the same two runs scored over D2-D89; on 2026-09-14 this figure took
on both scoring windows and the companion was retired to
model_output/superseded_20260914/. The three features below were what distinguished
the two, and are now simply what this figure has.

    SHARED Y AXIS
        Both periods on identical limits, so the panels can be read against
        each other rather than only against their own observations. That
        comparison is the point: period 2 sits almost entirely above period 1,
        which is the post-Isabel recovery and the nourishment era showing up as
        a whole-island shift in the rate, not as a local feature.

    THE LOESS CURVE ONLY, NOT THE SPLICED TARGET
        The calibration target is not one curve: GIS 1-10 are raw per-domain
        means and D11 north is the 10-domain LOESS. That splice is right for
        calibrating -- the raw means keep the short-wavelength signal the
        source/sink field has to answer for -- but it makes an awkward figure,
        because the eye reads a change of estimator as a change of coast.
        Here the LOESS is drawn throughout.

        D1-D10 IS DRAWN DASHED. The project excludes the LOESS there by
        convention (`skip_southern_domains = 10`), because the smoother is
        poorly constrained at the end of its range and Cape Point's
        attachment-detachment cycle is exactly the short-wavelength signal a
        5 km smoother destroys. Dashing it shows the data without implying it
        carries the same weight.

    BOTH SCORING WINDOWS, PRINTED PER PANEL
        D11-D89 is where the LOESS curve and the calibration target are the
        SAME numbers, so that statistic is the one for the line actually drawn.
        D2-D89 is the project's canonical skill window -- rmse_interior_m_yr in
        run_index.csv, and what the groin fit and the source/sink convergence
        were ranked on -- and additionally takes in D2-D10, where the target is
        the raw spliced domain mean and the model is at its worst (RMSE ~1.1
        against ~0.45 north of it). Both are correct for what they measure and
        mixing them is the error, which is why both are on the panel and
        labelled. They are COMPUTED here, both of them: the D2-D89 pair used to
        be a literal quoted from the companion figure, and it sat twelve days
        stale after the 1984 run was remade on topography v2.

ALSO ADDED FOR PRESENTATION
    An observational spread band (+/- 1 SD of transect rates within each
    domain, from `unc_m_yr`'s parent transect table). A domain is a 500 m
    average over ~10 transects, and showing that spread makes clear how much
    real alongshore variability the domain mean hides -- and therefore which
    model-observation gaps are meaningful and which sit inside the noise.

    Place names along the top, from PHYSICAL_ZONES, so an audience can locate
    features without a separate map.

Usage:
    python hindcast_final_figure_loess.py [--preset edgeBE|zeroBE]

ON THE 1996 -> 2010 -> 2024 CHAIN since 2026-09-18, nogroin arm, edgeBE by
default. calibBE is kept in PRESETS but is not solved on this chain.

Writes output/comparisons/hindcast_calibrated/hindcast_<preset>_loess_reference.png
and, with PUBLISH, output/figures/shoreline/hindcast_<preset>.png with its
caption in supporting/CAPTIONS.md. Both are the same house-style figure; the
title and note that were drawn on the canvas are the caption (2026-09-18).

Author: Hannah A. Henry, UNC CECL
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
# The manuscript copy goes with the other figures, by subject
# (output/figures/README.md); the presentation version stays in OUT_DIR.
FIGURES_DIR = _hs.figure_dir("shoreline")
# ON since the layout was redrawn for the house-style column (2026-09-18).
# Before that the labels, title and legend overlapped, and it was held off so
# a broken figure could not land in output/figures/.
PUBLISH = True

sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

# HOUSE STYLE: one typeface and one palette across every figure in this
# project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. Applied at
# MODULE level even though this file imports matplotlib inside a function --
# apply_style() sets rcParams, so it only has to run before the figure is
# built, and it was never called at all until 2026-09-17.
from site_layer.hat_figure_style import (apply_style, figsize,  # noqa: E402
                              FIG_W_DOUBLE, record_caption, save,
                              C_1984, C_1997, INK_MUTED, DOMAIN_AXIS_LABEL)
apply_style()
from cascade_pipeline.run_layout import resolve                 # noqa: E402
from cascade_pipeline.run_registry import find_run_dir          # noqa: E402
from site_layer.hat_observed_rates import lrr_csv                # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_PERIODS     # noqa: E402

LOESS_PATH = (PROJECT_BASE_DIR / "scripts" / "input_prep" / "7-source-sink"
              / "2-calibrate" / "be_zone_residual_fit.py")

# THE CANONICAL CHAIN, 1996 -> 2010 -> 2024 (moved from 1984/2004 on
# 2026-09-18). Ends come from HATTERAS_PERIODS. The full-management arm is
# road_bdm in 1996-2010, which schedules no fill, and road_bdm_nourish after.
# GROIN: the matrix on this chain is nogroin only -- the groin is still being
# fitted in the sweep -- so the figure draws the nogroin arm and says so.
GROIN = "nogroin"
_CHAIN = ((1996, "road_bdm", "a", C_1984), (2010, "road_bdm_nourish", "b", C_1997))
PERIODS = {}
for _st, _scen, _panel, _colour in _CHAIN:
    _end = HATTERAS_PERIODS[_st]["end_year"]
    PERIODS[f"{_st}_{_end}"] = dict(panel=_panel, label=f"{_st}–{_end}",
                                    start=_st, end=_end, scenario=_scen,
                                    colour=_colour)
LOCKED = (1, 90)
SKIP_SOUTH = 10               # matches LOESS_CONFIG.skip_southern_domains
SCORE_DOMAINS = range(SKIP_SOUTH + 1, 90)
RESERVED_COLOUR = "#FF8C00"

# Non-overlapping display spans. PHYSICAL_ZONES overlaps at D9-D10 (Cape Point
# and Buxton-Avon both claim them) and assign_physical_zone resolves that by
# first match; a label strip has to pick one, so it picks the same one.
PLACE_LABELS = [
    (1, 10, "Cape Point"), (11, 20, "Buxton"), (21, 31, "Avon"),
    (32, 59, "Mid-island"), (60, 74, "Wimble Shoals"),
    (75, 83, "Rodanthe"), (84, 90, "Pea Island"),
]


def analysis_module():
    spec = importlib.util.spec_from_file_location("_loess", LOESS_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def loess_and_spread(module, start, csv_path):
    """(loess_by_domain, sd_by_domain) for one period.

    The LOESS is the widest window's smoothed curve across ALL domains --
    build_target_table would splice raw means over D1-D10, which is what this
    figure is deliberately not doing.
    """
    from cascade_pipeline.coastsat_loess import (CoastSatDataset,
                                                 build_coastsat_series)
    series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {start}", period_start=start,
                         csv_path=csv_path)], start, module.LOESS_CONFIG)
    cs = series[0]
    match = [w for w in cs["windows"] if w["window"] == module.TARGET_WINDOW]
    if not match:
        raise ValueError(f"window {module.TARGET_WINDOW} not computed")
    loess = dict(zip(np.asarray(match[0]["gis_x"], dtype=int),
                     np.asarray(match[0]["smoothed"], dtype=float)))
    frame = pd.read_csv(csv_path, usecols=["transect_id", "domain_number",
                                           "lrr_m_yr"])
    frame = frame.dropna(subset=["domain_number", "lrr_m_yr"])
    frame["domain"] = frame["domain_number"].astype(int)
    frame["lrr"] = frame["lrr_m_yr"].astype(float)

    # A domain is 500 m holding ~10 transects. Plotted at the domain integer
    # they stack into a vertical column, which reads as one uncertain value
    # rather than as an alongshore gradient. Spread them across the domain in
    # transect order instead. The index is the LAST numeric field of the id --
    # the first is the CoastSat site number and grabbing it sorts by site.
    order = frame["transect_id"].astype(str).str.split("_").str[-1]
    frame["t_index"] = pd.to_numeric(order, errors="coerce")
    frame = frame.sort_values(["domain", "t_index"])
    rank = frame.groupby("domain").cumcount()
    count = frame.groupby("domain")["lrr"].transform("size")
    frame["x"] = frame["domain"] - 0.5 + (rank + 0.5) / count

    sd = frame.groupby("domain")["lrr"].std().to_dict()
    return loess, {int(k): float(v) for k, v in sd.items()}, frame


def run_rates(period, preset, scenario):
    name = f"HAT_{period}_{preset}_{scenario}_{GROIN}"
    # Resolved rather than joined by hand: a hand-built path has no slot for
    # the arm component and reads an arm-scoped run as missing. The retired
    # companion carried a twin of this function, in
    # model_output/superseded_20260914/HAT_hindcast_final_figure.py.
    try:
        run_dir = find_run_dir(RAW_RUNS, name, period, preset)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"{exc}  Run HAT_run_all.py first.") from None
    path = resolve(run_dir, "rate_csv", name)
    return pd.read_csv(path).set_index("gis_domain")["lrr_m_yr"], name


# THE PRESET OWNS EVERY WORD THAT NAMES THE CONFIGURATION. Added 2026-09-14
# alongside the edgeBE companion. --preset already existed here, but the
# filename, the title and the footnote were hardcoded to calibBE, so any other
# preset overwrote the calibrated figure with a page still calling itself
# calibrated. `stem` must stay in step with the same table in
# the retired companion's table (model_output/superseded_20260914/), which named
# the same stems -- kept aligned so its output and this figure's still sort
# together in the folder.
PRESETS = {
    "calibBE": dict(
        stem="hindcast_calibrated",
        published="hindcast_calibrated.png",
        title="Calibrated CASCADE hindcast against the CoastSat LOESS reference, Cape Hatteras",
        series="CASCADE, calibrated",
        config="calibrated source/sink field, full management (roadway + beach/dune; "
               "nourishment in 2004–2024), groin active at M = 60 m/yr, f = 0.6.",
        zone_note="Source/sink corrections were applied only within a zone set fixed "
                  "before calibration; D5–D7 were reserved for the groin, so the misfit "
                  "there is the groin module's and is not absorbed by the sediment budget."),
    "edgeBE": dict(
        stem="hindcast_edgeBE",
        published="hindcast_edgeBE.png",
        title="CASCADE hindcast with the source/sink calibration removed, against the CoastSat LOESS reference, Cape Hatteras",
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
        title="CASCADE hindcast with no source/sink field, against the CoastSat LOESS reference, Cape Hatteras",
        series="CASCADE zeroBE",
        config="zeroBE: no background-erosion field anywhere, the D1/D90 edges included. "
               "Full management (roadway + beach/dune; nourishment in 2010–2024); the "
               "Buxton groin is off, the nogroin matrix arm.",
        zone_note="No source/sink field was applied anywhere, edges included, so the edges "
                  "carry no absorber and are free to run away from the target."),
}


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

    # Label heights are tuned PER FIGURE, not in the shared config, because
    # they depend on this figure's y range and on what is drawn where. The
    # defaults (groin 0.68, piers 0.76) put "Buxton Groin" through the D1-D10
    # transect scatter and "Avon Pier" through both shoreline curves.
    #
    #   Buxton Groin  0.78  below the "Buxton" town label (0.92 put the two on
    #                       top of each other at the 09-18 column width) and
    #                       still above the scatter in (a); the
    #                       "Buxton" town label is centred at D7.5 so its box
    #                       clears D5.5 horizontally.
    #   Avon Pier     0.85  the "Avon" town span is ALSO centred on D26, and it
    #                       sits at 0.90, so the pier has to hang below it.
    #   Rodanthe Pier 0.72  the "Rodanthe" village line is at D80 with its label
    #                       at 0.84; dropping the pier clears that, and D79 is
    #                       deeply erosional in both periods so the mid-panel is
    #                       empty there.
    annotations = dataclasses.replace(
        HATTERAS_ANNOTATIONS,
        groin_label_y=0.78,
        piers={"Avon Pier": (26, 0.85), "Rodanthe Pier": (79, 0.72)})

    module = analysis_module()
    reserved = list(module.GROIN_RESERVED_DOMAINS)

    panels = {}
    for key, meta in PERIODS.items():
        csv = str(lrr_csv(meta["start"], meta["end"]))
        loess, sd, transects = loess_and_spread(module, meta["start"], csv)
        model, run_name = run_rates(key, args.preset, meta["scenario"])
        # WHAT THE COMPANION FIGURE REPORTS, COMPUTED, NOT REMEMBERED. The
        # footnote quotes the D2-D89 score so a reader can see why it differs
        # from this figure's D11-D89 one. It used to be a literal, and on
        # 2026-09-14 it was found to be twelve days stale: written 09-02, it
        # still held the pre-v2 numbers after the 1984 run was remade on
        # topography v2 on 09-07. Same target, same model, same window as
        # the runner's own interior metric -- so it cannot drift from it.
        spliced = module.load_observed(meta["start"], csv)[1]
        wide = [g for g in range(2, 90)
                if g in model.index and not np.isnan(spliced.get(g, np.nan))]
        wm = np.array([model[g] for g in wide])
        wo = np.array([spliced[g] for g in wide])
        companion = dict(rmse=float(np.sqrt(np.mean((wm - wo) ** 2))),
                         bias=float(np.mean(wm - wo)),
                         corr=float(np.corrcoef(wm, wo)[0, 1]),
                         n=len(wide))
        panels[key] = dict(loess=loess, sd=sd, transects=transects,
                           model=model, run=run_name, companion=companion,
                           **meta)

    gis = np.arange(1, 91)
    lo, hi = np.inf, -np.inf
    for d in panels.values():
        for source in (d["loess"], d["model"]):
            values = np.array([source.get(g, np.nan) for g in gis], dtype=float)
            lo = min(lo, np.nanmin(values))
            hi = max(hi, np.nanmax(values))
    # The D1-D10 transect scatter is drawn too, and at Cape Point it runs past
    # both curves; a limit set on the curves alone clipped it (09-18). It only
    # widens the axis as far as it reaches, without the label headroom the
    # curves get -- the place names sit well north of D1-D10.
    s_lo, s_hi = np.inf, -np.inf
    for d in panels.values():
        south_t = d["transects"][d["transects"]["domain"] <= SKIP_SOUTH]["lrr"]
        if len(south_t):
            s_lo = min(s_lo, float(south_t.min()))
            s_hi = max(s_hi, float(south_t.max()))
    pad = 0.10 * (hi - lo)
    # Asymmetric: the top needs room for the place-name strip, the bottom does
    # not, and a symmetric pad on a shared axis wastes a band of panel (b).
    ylim = (min(lo - pad * 0.45, s_lo - pad * 0.1),
            max(hi + pad * 1.5, s_hi + pad * 0.1))

    # HOUSE STYLE (2026-09-18). This block drew for a 15 in presentation page
    # -- 11-15 pt type, a suptitle and a paragraph of note on the canvas --
    # and when figsize() put it on the 190 mm column on 09-17 the ylabels
    # clipped and the note ran under the legend. Now: house type sizes, the
    # letter and period at the left above each panel, the two scoring windows
    # at the right above it (they used to sit in a box over the curves), the
    # place names once on (a), and the title and note in CAPTIONS.md.
    figure, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.0),
                                sharex=True, sharey=True)

    summary = []
    for i, (axis, (key, d)) in enumerate(zip(axes, panels.items())):
        obs = np.array([d["loess"].get(g, np.nan) for g in gis], dtype=float)
        mod = np.array([d["model"].get(g, np.nan) for g in gis], dtype=float)
        spread = np.array([d["sd"].get(g, np.nan) for g in gis], dtype=float)

        # Communities, village centres, piers, groins and shoal zones, from the
        # shared layer every other cascade_pipeline figure uses -- so this
        # figure cannot disagree with the annotated run plots about where the
        # Buxton groin or the Avon pier is. Named on (a) only; (b) has the
        # same bands, and the same eight names twice is clutter.
        add_geographic_annotations(axis, annotations, label=(i == 0))

        # The D5-D7 (groin-reserved) and D1/D90 (locked) spans were dropped
        # from this figure: with the geographic annotation layer in place the
        # Buxton area carried four overlapping fills and read as clutter. Both
        # facts are stated in the caption instead.
        axis.axvspan(0.5, SKIP_SOUTH + 0.5, facecolor="none",
                     edgecolor="#999999", hatch="\\\\\\", linewidth=0.0,
                     alpha=0.30, zorder=0)

        axis.fill_between(gis, obs - spread, obs + spread, color="#999999",
                          alpha=0.22, zorder=2, linewidth=0,
                          label="observed spread (±1 SD of transects)")

        # Individual transects over D1-D10. The LOESS is dashed there because
        # the smoother is unreliable at the end of its range; the scatter is
        # the actual evidence, and it shows the Cape Point spread the smooth
        # curve cannot represent.
        south_t = d["transects"][d["transects"]["domain"] <= SKIP_SOUTH]
        axis.plot(south_t["x"], south_t["lrr"], linestyle="none",
                  marker="o", markersize=1.8, color=INK_MUTED, alpha=0.6,
                  zorder=5, label="CoastSat transect LRR (D1–D10)")
        axis.fill_between(gis, obs, mod, where=~(np.isnan(obs) | np.isnan(mod)),
                          color=d["colour"], alpha=0.18, zorder=3, linewidth=0,
                          label=f"misfit, {d['label']}")

        south = gis <= SKIP_SOUTH
        axis.plot(gis[~south], obs[~south], color="#1A1A1A", linewidth=1.8,
                  zorder=7, label="CoastSat LOESS (10-domain)")
        axis.plot(gis[south], obs[south], color="#1A1A1A", linewidth=1.3,
                  linestyle=(0, (4, 2)), zorder=7,
                  label="LOESS, D1–D10 (excluded by convention)")
        axis.plot(gis, mod, color=d["colour"], linewidth=1.6, zorder=6,
                  label=f"{vocab['series']}, {d['label']}")
        axis.axhline(0.0, color="#AAAAAA", linewidth=0.6, zorder=4)

        shared = [g for g in SCORE_DOMAINS
                  if g in d["model"].index and not np.isnan(d["loess"].get(g, np.nan))]
        mm = np.array([d["model"][g] for g in shared])
        oo = np.array([d["loess"][g] for g in shared])
        rmse = float(np.sqrt(np.mean((mm - oo) ** 2)))
        bias = float(np.mean(mm - oo))
        corr = float(np.corrcoef(mm, oo)[0, 1])
        summary.append((d["label"], d["run"], rmse, bias, corr, len(shared)))

        # BOTH SCORING WINDOWS, ON THE FIGURE. Added 2026-09-14. This figure
        # scores D11-D89, the span where the LOESS curve IS the calibration
        # target; the project's canonical skill column (rmse_interior_m_yr,
        # recorded for every run in run_index.csv) is D2-D89, which also takes
        # in D2-D10, where the target is the raw spliced mean and the model is
        # at its worst. Printing only the narrow one reads as a better model
        # rather than a shorter ruler, and printing it in a separate figure is
        # what let the two drift apart.
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

    # One key for both panels, in reading order: the observations, then the
    # model and its misfit per period, then the place layer.
    found = {}
    for axis in axes:
        for h, l in zip(*axis.get_legend_handles_labels()):
            found.setdefault(l, h)
    order = ["CoastSat LOESS (10-domain)",
             "LOESS, D1–D10 (excluded by convention)",
             "CoastSat transect LRR (D1–D10)",
             "observed spread (±1 SD of transects)"]
    for d in panels.values():
        order += [f"{vocab['series']}, {d['label']}", f"misfit, {d['label']}"]
    labels = [l for l in order if l in found]
    handles = [found[l] for l in labels]
    extra = annotation_legend_handles(annotations)
    handles += extra
    labels += [h.get_label() for h in extra]
    # Figure-level, below both panels. Too many entries to sit inside an axis
    # without covering something.
    figure.legend(handles, labels, loc="lower center", ncol=3, fontsize=7.5,
                  frameon=False, bbox_to_anchor=(0.5, 0.0))
    _leg_rows = math.ceil(len(handles) / 3)
    figure.tight_layout(rect=(0.01, 0.026 * _leg_rows + 0.01, 1, 1))

    caption_text = (
        vocab["title"] + ". Configuration: " + vocab["config"]
        + " The observed curve is the 10-domain LOESS of CoastSat transect "
        "rates; over D1–D10 (hatched) it is dashed because the project "
        "excludes the LOESS there (the smoother is poorly constrained at the "
        "end of its range, and Cape Point's attachment-detachment cycle is "
        "short-wavelength signal a 5 km smoother removes), and the individual "
        "transect rates are drawn instead. The grey band is ±1 SD of the "
        "transect rates within each domain, and the tinted band between the "
        "curves is the misfit. Both scoring windows are printed above each "
        "panel. D11–D89 is the span where the LOESS curve and the grading "
        "target are identical numbers, so that statistic describes the line "
        "actually drawn. D2–D89 is the project's canonical skill window "
        "(rmse_interior_m_yr in run_index.csv), which also takes in D2–D10, "
        "where the target is the raw spliced domain mean and the model is at "
        "its worst; the gap between the two numbers is that strip alone, not "
        "the smoothing. D1 and D90 are in neither: they are locked boundary "
        "absorbers, solved by buffer-cell reproduction, and would dominate. "
        + vocab["zone_note"])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"{vocab['stem']}_loess_reference.png"
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
