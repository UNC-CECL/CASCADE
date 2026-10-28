#!/usr/bin/env python3
"""
A finished hindcast run's figures and animations, drawn from its run folder alone.

    python scripts/real-world-cases/Hatteras/figure_making/model_output/run_figures.py <run_dir> [--no-gifs]

    from run_figures import make_run_figures      # what HAT_hindcast.py calls at the end of a run

Reads the run's metadata, shoreline matrix and, when the run wrote them, its overwash and road
relocation arrays; never the model. Moved out of HAT_hindcast.py on 2026-10-23 so the runner only
runs the model. Details: scripts/real-world-cases/Hatteras/figure_making/model_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# The packages live in Hatteras/ (the folder holding hatteras_paths.py), which is not installed
_HAT_DIR = next(_p for _p in Path(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from cascade_pipeline import reports  # noqa: E402
from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    CoastSatDataset, LowessConfig, build_coastsat_series)
from cascade_pipeline.hindcast import build_coastsat_target  # noqa: E402
from cascade_pipeline.plotting.dune_crest_comparison import plot_dune_crest_comparison  # noqa: E402
from cascade_pipeline.plotting.net_change_comparison import plot_net_change_comparison, smooth_model  # noqa: E402
from cascade_pipeline.plotting.rate_figure import plot_rate_figure  # noqa: E402
from cascade_pipeline.plotting.shoreline_gif import GifConfig, make_all_shoreline_gifs  # noqa: E402
from cascade_pipeline.run_info import RunInfo  # noqa: E402
from cascade_pipeline.run_layout import resolve, write_path  # noqa: E402
from cascade_pipeline.shoreline import compute_change_rate, compute_lrr  # noqa: E402
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_BASE_DIR  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS, HATTERAS_GEOMETRY_EXTENDED, HATTERAS_PERIODS)

# --- CONFIG ------------------------------------------------------------------
HATTERAS_DATA_BASE = WORKSPACE / "data" / "hatteras_init"

# 7-domain smoothing, the group's range (2026-09-28); the runner scores against the same
LOWESS_CONFIG = LowessConfig(window_domains=(7,), skip_southern_domains=10)
TARGET_WINDOW = 7               # the runner's TARGET_WINDOW: the GIFs' CoastSat target line
GIF_FPS = 2                     # frames (model years) per second; was 3 until 2026-10-26

# The estimator the figures draw: "lrr" matches the CoastSat target; "endpoint" only redraws old figures
RATE_ESTIMATOR = "lrr"
FLIP_SIGN_MODEL = True          # x_s_TS increases landward; flip so up = seaward
PLOT_REAL_DOMAINS_ONLY = True   # GIS 1-90 axis; False adds the buffer domains

# range: "real" | "all" | "groin" | "groin_span" | (lo, hi); mode: position | displacement | difference
GIF_JOBS = [
    dict(range="real", mode="displacement"),
    dict(range="real", mode="position"),
    dict(range="groin", mode="position", pad=9),
    dict(range="groin", mode="difference", pad=9),
]

# The arrays the runner saves for the figures, beside the shoreline matrix
OVERWASH_NPY = "{run}_overwash.npy"                  # Qow, years x padded domains
RELOCATIONS_NPY = "{run}_road_relocations.npy"       # True in the year a domain's road moved
# -----------------------------------------------------------------------------


# Every CoastSat window, the run's own active; an extended geometry swaps in its own rows
def coastsat_series(start_year):
    datasets = [
        CoastSatDataset(label="CoastSat LRR (1984-2004)", period_start=1984,
                        csv_path=str(COASTSAT_BASE_DIR / "1984_2004" / "transect_lrr_full.csv")),
        CoastSatDataset(label="CoastSat LRR (2004-2024)", period_start=2004,
                        csv_path=str(COASTSAT_BASE_DIR / "2004_2024" / "transect_lrr_full.csv")),
        *[CoastSatDataset(
            label=f"CoastSat LRR ({y}-{HATTERAS_PERIODS[y]['end_year']})", period_start=y,
            csv_path=str(COASTSAT_BASE_DIR / f"{y}_{HATTERAS_PERIODS[y]['end_year']}"
                         / "transect_lrr_full.csv"))
          for y in (1996, 2009)],
    ]
    if HATTERAS_GEOMETRY_EXTENDED:
        datasets = [
            CoastSatDataset(label=ds.label + " + extension", period_start=ds.period_start,
                            csv_path=str(Path(ds.csv_path).parent / "ext" / "transect_lrr_with_base.csv"))
            if ds.period_start == start_year else ds
            for ds in datasets]
    return build_coastsat_series(datasets, active_period_start=start_year,
                                 lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)


# The GIFs' target line; short, so the three-column legend fits the frame
def target_label(start_year, end_year):
    return f"CoastSat {end_year} shoreline"


# The run's own saved array of one kind, or None if it did not write one
def _optional_npy(run_dir, template, run_name):
    path = Path(run_dir) / template.format(run=run_name)
    return np.load(path) if path.is_file() else None


def make_run_figures(run_dir, gifs=True, show=False, cs_series=None):
    """Draws one finished run's figures and animations into its run folder.

    Args:
        run_dir: The run folder, holding <run>_run_metadata.json and <run>_shoreline_matrix.npy.
        gifs: Also draw the shoreline animations (slow).
        show: Show each figure as well as saving it.
        cs_series: The CoastSat series for the run's period, if the caller already built it.

    Returns:
        The list of saved animation paths.
    """
    # An empty AnnotationConfig would strip every figure's geography with no error: fail here
    if not any([HATTERAS_ANNOTATIONS.town_spans, HATTERAS_ANNOTATIONS.village_lines,
                HATTERAS_ANNOTATIONS.piers, HATTERAS_ANNOTATIONS.groins,
                HATTERAS_ANNOTATIONS.shoal_zones]):
        raise RuntimeError("HATTERAS_ANNOTATIONS is empty -- every figure would render with no "
                           "geographic layer. Check hatteras_site_config.")

    run_dir = Path(run_dir)
    run_name = run_dir.name
    meta = json.loads((run_dir / f"{run_name}_run_metadata.json").read_text(encoding="utf-8"))
    shoreline_m = np.load(resolve(run_dir, "matrix", run_name))
    # Every figure below draws the site config's buffers; a run made with others (15 a side before 2026-10-21)
    # would be drawn domains off
    if shoreline_m.shape[1] != HATTERAS_DOMAINS.total_domains:
        raise SystemExit(f"{run_name}: {shoreline_m.shape[1]} domains, the site config draws "
                         f"{HATTERAS_DOMAINS.total_domains}")
    period = meta["period"]
    start_year, end_year = int(period["start_year"]), int(period["end_year"])
    run_years = int(period["run_years"])
    wave_climate = meta["scenario"]["wave climate"]

    run = RunInfo(
        run_name=run_name, run_dir=str(run_dir),
        start_year=start_year, end_year=end_year,
        Hs=float(meta["wave climate"]["wave_height_m"]),
        first_model_year=int(period["first_model_year"]),
        flip_sign_model=FLIP_SIGN_MODEL,
        background_erosion_on=bool(meta["source/sink"]["background_erosion_on"]),
    )
    if cs_series is None:
        cs_series = coastsat_series(start_year)

    # The model drawn smoothed as CoastSat is: 7-domain LOWESS, GIS 1-10 and 90 raw (2026-10-20)
    if RATE_ESTIMATOR == "lrr":
        plotted_rate, _ = compute_lrr(shoreline_m, span_years=run_years, flip_sign=FLIP_SIGN_MODEL)
    else:
        plotted_rate = compute_change_rate(shoreline_m, span_years=run_years, flip_sign=FLIP_SIGN_MODEL)
    plotted_smoothed = smooth_model(plotted_rate)

    # One rate figure since 2026-10-26, styled as the position figure; the plain real-domains one is no longer drawn
    plot_rate_figure(
        plotted_smoothed, cs_series, run, write_path(run_dir, "figure_rate_buffers", run_name),
        lowess_config=LOWESS_CONFIG, estimator=RATE_ESTIMATOR, show=show)

    # The end-minus-start shoreline against the CoastSat net change; skipped for a window with no target
    plot_net_change_comparison(
        shoreline_m, run, write_path(run_dir, "figure_net_change", run_name),
        wave_climate=wave_climate, flip_sign=FLIP_SIGN_MODEL, show=show)

    # The modelled dune crests against every lidar survey in the window, ten domains a figure; needs the saved model
    plot_dune_crest_comparison(run_dir)

    # Where the run overwashed against the washover in the imagery; the map needs the GIS drive, so a failure only skips it
    qow = _optional_npy(run_dir, OVERWASH_NPY, run_name)
    if qow is None:
        print("  overwash vs observed  skipped: the run saved no overwash array")
    else:
        sys.path.insert(0, str(_HAT_DIR / "input_prep" / "8-overwash-analysis" / "4-vs-model"))
        try:
            from overwash_vs_model import run_figure
            run_figure(qow, run_dir, run_name, start_year,
                       HATTERAS_DATA_BASE / meta["identity"]["storm_file"])
        except (Exception, SystemExit) as _err:
            print(f"  overwash vs observed  skipped: {_err}")
        # And how much it overwashed over the whole run, where and when
        try:
            from overwash_total import run_figure as overwash_total_figure
            overwash_total_figure(qow, run_dir, run_name, start_year)
        except (Exception, SystemExit) as _err:
            print(f"  overwash total        skipped: {_err}")

    # The CoastSat end-year shoreline, the run's scoring target x the span (was the surveyed dune line
    # until 2026-10-26; the run starts from shorelines and is graded on CoastSat); buffers are NaN
    target_m, observed_change_m = build_coastsat_target(
        shoreline_m[0], cs_series, LOWESS_CONFIG, HATTERAS_DOMAINS,
        TARGET_WINDOW, end_year - start_year)
    reports.target_misfit_report(
        target_m=target_m, observed_change_m=observed_change_m,
        shoreline_m=shoreline_m, label=target_label(start_year, end_year),
        geometry=HATTERAS_DOMAINS)

    gif_paths = []
    if gifs:
        # A groin run differences against its paired no-groin run, named in its metadata
        baseline_name = meta.get("groin", {}).get("baseline_run")
        baseline_npy = (resolve(run_dir.parent / baseline_name, "matrix", baseline_name)
                        if baseline_name else None)
        relocations = _optional_npy(run_dir, RELOCATIONS_NPY, run_name)
        if relocations is not None:
            print(f"  roadway relocations   {int(relocations.sum())} events "
                  f"across {int(relocations.any(axis=0).sum())} domains")
        gif_paths = make_all_shoreline_gifs(
            shoreline_m, run, GIF_JOBS,
            baseline_npy=str(baseline_npy) if baseline_npy is not None and baseline_npy.exists() else None,
            target_m=target_m, relocations=relocations,
            domains=HATTERAS_DOMAINS, annotations=HATTERAS_ANNOTATIONS,
            gif_config=GifConfig(
                fps=GIF_FPS, year_stride=1, annotate=True, auto_open=False, keep_frames=False,
                save_matrix=False,         # the runner saved the matrix
                ocean_at_bottom=True,      # Hatteras' real cross-shore layout
                baseline_label="no-groin baseline",
                target_label=target_label(start_year, end_year)))

    plt.close("all")
    return gif_paths


# Run: one run folder from the command line
def main():
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("run_dir", type=Path, help="the run folder")
    ap.add_argument("--no-gifs", action="store_true", help="skip the animations")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    make_run_figures(args.run_dir, gifs=not args.no_gifs)


if __name__ == "__main__":
    main()
