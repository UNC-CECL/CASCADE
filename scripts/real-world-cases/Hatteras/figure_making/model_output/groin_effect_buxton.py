#!/usr/bin/env python3
"""
Shoreline position change at Buxton (GIS 1-10), the full-management run with and without the groin.

    python groin_effect_buxton.py [--window 1996_2009]

Draws the run's shoreline_position_change_with_buffers figure with the same
plotting call, cuts it to GIS 1-10 at +/-55 m and adds the no-groin run. Reads the
edgeBE pair in PAIRS and the CoastSat net change for the window.
Writes output/comparisons/groin_effect/<window>/.
Details: scripts/figure_making/model_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import numpy as np

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in Path(__file__).resolve().parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
from cascade_pipeline.annotations import annotation_legend_handles  # noqa: E402
from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    CoastSatDataset, LowessConfig, build_coastsat_series)
from cascade_pipeline.plotting.net_change_comparison import smooth_model  # noqa: E402
from cascade_pipeline.plotting.rate_comparison import (  # noqa: E402
    DEFAULT_RATE_COMPARISON, plot_annotated_rate_comparison)
from cascade_pipeline.run_info import RunInfo  # noqa: E402
from cascade_pipeline.run_layout import resolve  # noqa: E402
from cascade_pipeline.shoreline import compute_change_rate  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, COMPARISONS_ROOT, mark_offaxis, offaxis_clause, record_caption, save)
from site_layer.hat_observed_rates import (  # noqa: E402
    NET_CHANGE_WINDOWS, WINDOW_ROLE, net_change_dir)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
MATRIX = REPO / "output" / "raw_runs" / "matrix"
# Full management, edgeBE, groin on vs off; nothing else differs within a pair
PAIRS = {
    "1996_2009": ("HAT_1996_2009_edgeBE_offsetmetres_road_bdm_groin",
                  "HAT_1996_2009_edgeBE_offsetmetres_road_bdm_nogroin"),
    "2009_2025": ("HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_groin",
                  "HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_nogroin"),
}
GIS_RANGE = (1, 10)
Y_HALF = 55.0
FLIP_SIGN_MODEL = True    # x_s grows landward; positive is seaward on the figure
LOWESS_CONFIG = LowessConfig(window_domains=(7,), skip_southern_domains=10)
C_NOGROIN = C["BASE"]
# Legend entries from the run figure that have something to point at in GIS 1-10
KEEP_FEATURES = ("Community", "Groin")
OUT = COMPARISONS_ROOT / "groin_effect"
# -----------------------------------------------------------------------------


# End minus start (m, seaward positive) over every padded domain of one run
def position_change(run_name: str, window: str) -> np.ndarray:
    run_dir = MATRIX / window / "edgeBE" / run_name
    shoreline_m = np.load(resolve(run_dir, "matrix", run_name))
    return compute_change_rate(shoreline_m, span_years=1, flip_sign=FLIP_SIGN_MODEL)


# The CoastSat net-change series for the window, as the run figure builds it
def observed(start: int, end: int) -> list:
    tag = f"{start}_{end}"
    dataset = CoastSatDataset(
        label=f"CoastSat net change ({start}–{end})", period_start=start,
        csv_path=str(net_change_dir(start, end) / f"transect_net_change_{tag}.csv"),
        rate_col="net_change_m")
    return build_coastsat_series([dataset], active_period_start=start,
                                 lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)


# Run: draw the run figure for the groin run, zoom it, add the no-groin run, caption, save
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", default="1996_2009", choices=sorted(PAIRS))
    args = ap.parse_args()
    start, end = (int(y) for y in args.window.split("_"))
    groin_run, nogroin_run = PAIRS[args.window]
    with_groin = smooth_model(position_change(groin_run, args.window))
    without = smooth_model(position_change(nogroin_run, args.window))
    cs_series = observed(start, end)

    # The run figure's own call and settings (net_change_comparison.py), y axis at +/-Y_HALF
    role = WINDOW_ROLE.get((start, end), "").replace(" (held out)", "")
    config = dataclasses.replace(
        DEFAULT_RATE_COMPARISON, quantity="position", publication_text=True,
        show_features=True, plot_domain_means=True, raw_lrr_southern_only=True,
        observed_label="CoastSat (7-domain LOWESS)", model_label="CASCADE, groin",
        title=f"Modelled and observed shoreline position change, {start}–{end}",
        ylim=(-Y_HALF, Y_HALF),
        # The observed domain mean is the only CoastSat line here; drawn at the model's weight
        domain_mean_lw=2.0, domain_mean_alpha=1.0,
        subtitle=" · ".join(p for p in (role, "Buxton, GIS 1–10",
                                        "with and without the groin") if p))
    run = RunInfo(run_name=groin_run, run_dir=str(MATRIX / args.window / "edgeBE" / groin_run),
                  start_year=start, end_year=end, Hs=2.0, flip_sign_model=FLIP_SIGN_MODEL)
    fig, ax = plot_annotated_rate_comparison(
        with_groin, cs_series, run, estimator=None, domains=HATTERAS_DOMAINS,
        annotations=HATTERAS_ANNOTATIONS, lowess_config=LOWESS_CONFIG, config=config)

    # Zoom: one tick per domain; hide the labels of features off to the north
    ax.set_xlim(GIS_RANGE[0] - 0.5, GIS_RANGE[1] + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    for txt in ax.texts:
        x_is_data = txt.get_transform() not in (ax.transAxes, fig.transFigure)
        if x_is_data and not GIS_RANGE[0] - 0.5 <= txt.get_position()[0] <= GIS_RANGE[1] + 0.5:
            txt.set_visible(False)
    # Domain d's transects run d to d+1 in along-coast/500 + 1; centre them on d like the model
    for coll in ax.collections:
        if "transect" in str(coll.get_label()):
            coll.set_offsets(coll.get_offsets() - [0.5, 0.0])

    gis = np.arange(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1)
    real = slice(HATTERAS_DOMAINS.start_real_index, HATTERAS_DOMAINS.end_real_index)
    ax.plot(gis, without[real], color=C_NOGROIN, lw=2.0, ls="--", zorder=5.8)
    # Triangles at the edge for model values past +/-Y_HALF; named in the caption
    in_view = (gis >= GIS_RANGE[0]) & (gis <= GIS_RANGE[1])
    offaxis = offaxis_clause(
        [(name, n) for name, n in (
            ("CASCADE, groin", mark_offaxis(ax, gis[in_view], with_groin[real][in_view], Y_HALF,
                                            color=HATTERAS_ANNOTATIONS.model_color)),
            ("CASCADE, no groin", mark_offaxis(ax, gis[in_view], without[real][in_view], Y_HALF,
                                               color=C_NOGROIN))) if n], Y_HALF)

    # The legend again, for what this view shows: two runs, domain means, transects
    for leg in list(fig.legends):
        leg.remove()
    blue = config.window_colors.get(max(LOWESS_CONFIG.window_domains),
                                    config.window_color_default)
    handles = [
        Line2D([0], [0], color=HATTERAS_ANNOTATIONS.model_color, lw=2.0, label="CASCADE, groin"),
        Line2D([0], [0], color=C_NOGROIN, lw=2.0, ls="--", label="CASCADE, no groin"),
        Line2D([0], [0], color=blue, lw=config.domain_mean_lw, alpha=config.domain_mean_alpha,
               label="CoastSat domain mean"),
        Line2D([0], [0], color=config.raw_color, marker="o", ms=3, ls="none",
               alpha=config.raw_scatter_alpha, label="Individual transects"),
    ] + [h for h in annotation_legend_handles(HATTERAS_ANNOTATIONS)
         if h.get_label() in KEEP_FEATURES]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False,
               handlelength=2.2, columnspacing=1.8)

    path = OUT / args.window / f"groin_effect_buxton_{args.window}.png"
    save(fig, path, bbox_inches="tight", facecolor="white", close=True)

    (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[(start, end)]
    keep = (gis >= GIS_RANGE[0]) & (gis <= GIS_RANGE[1])
    x, diff = gis[keep], (with_groin - without)[real][keep]
    i = int(np.argmax(np.abs(diff)))
    nourish = "nourishment, " if "nourish" in groin_run.split("_") else ""
    record_caption(path, (
        f"Shoreline position change at the cape end of Hatteras Island, GIS domains "
        f"{GIS_RANGE[0]}–{GIS_RANGE[1]}, {start}–{end} ({role.lower()}); positive is seaward. "
        f"The layout of the run's shoreline_position_change_with_buffers figure, cut to these "
        f"domains with the y axis at ±{Y_HALF:g} m. Modelled (CASCADE): the shoreline at the end "
        f"of the run minus the start, per 500 m domain (unsmoothed here, as south of domain 10 in "
        f"every figure), for the full-management run (road maintenance, beach and dune "
        f"management, {nourish}source/sink at the end domains only) with the Buxton groin "
        f"(orange, solid) and without it (grey, dashed); everything else is the same. The groin's "
        f"effect is largest at GIS {x[i]}, {diff[i]:+.1f} m. Both runs use the same edgeBE end "
        f"value at GIS 1, solved with the groin on, so near GIS 1 the no-groin line is the groin "
        f"removed with the ends held, not a separately fitted no-groin model. Observed "
        f"(CoastSat, blue): per transect, the mean position over {e0} to {e1} minus the mean "
        f"over {s0} to {s1} (dots, each centred on its own domain), and the mean of each "
        f"domain's transects (line). Shaded band: Buxton; dotted: the groin field, between "
        f"domains 5 and 6; dashed at 10.5: where the observed LOWESS begins.{offaxis} "
        f"Runs {groin_run} and {nogroin_run}."))
    print(f"Saved {path}")
    for g, d in zip(x, diff):
        print(f"  GIS {g:2d}: groin minus no groin {d:+6.1f} m")


if __name__ == "__main__":
    main()
