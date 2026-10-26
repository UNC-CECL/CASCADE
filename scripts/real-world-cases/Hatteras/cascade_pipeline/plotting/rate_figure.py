"""
A run's shoreline change rate against the CoastSat LRR, in the style of the position-change figure.

    from cascade_pipeline.plotting.rate_figure import plot_rate_figure

The one rate figure a run makes since 2026-10-26 (the plain real-domains figure is no longer
drawn by run_figures.py). Drawn by plot_annotated_rate_comparison and finished by
comparison_style.restyle, so it matches the position-change figure; the title and the y axis
name the quantity (rate, LRR, m/yr) so the two cannot be mistaken for each other.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-26
"""
from __future__ import annotations

import dataclasses
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from cascade_pipeline.plotting.comparison_style import restyle, run_subtitle
from cascade_pipeline.plotting.rate_comparison import (
    DEFAULT_RATE_COMPARISON, plot_annotated_rate_comparison)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS

# --- CONFIG ------------------------------------------------------------------
# One symmetric y axis per window (m/yr), so every run of a window compares (2026-10-26). Fitted to the
# 45 matrix runs and CoastSat: 1996-2009 spans -6.0 to +3.6, 2009-2025 -12.5 to +12.1 (the Cape Point
# dots). A window not listed scales to its data
Y_HALF_BY_START = {1996: 8.0, 2009: 16.0}
Y_TICKS_PER_HALF = 4
# -----------------------------------------------------------------------------


def plot_rate_figure(change_rate, cs_series, run, save_path, lowess_config,
                     estimator="lrr", show=False):
    """Draw the rate figure for one run; returns the path.

    Args:
        change_rate: The model's rate per padded domain (m/yr), smoothed as CoastSat is.
        cs_series: build_coastsat_series output, the run's window active.
        run: RunInfo for this run.
        save_path: Where the PNG goes (run_layout kind "figure_rate_buffers").
        lowess_config: The LOWESS the model and CoastSat were smoothed with.
        estimator: "lrr" or "endpoint", as the rate was computed.
    """
    window = (int(run.start_year), int(run.end_year))
    rate_name = "LRR" if estimator == "lrr" else "end-point rate"
    y_half = Y_HALF_BY_START.get(window[0])
    config = dataclasses.replace(
        DEFAULT_RATE_COMPARISON, quantity="rate", publication_text=True,
        show_features=True, plot_domain_means=False, raw_lrr_southern_only=True,
        smoothing_marker=False,   # read as a village line; the caption says where the LOWESS starts
        observed_label="CoastSat LRR (7-domain LOWESS)", model_label="CASCADE (7-domain LOWESS)",
        title=f"Modelled and observed shoreline change rate ({rate_name}, m/yr), "
              f"{window[0]}–{window[1]}",
        ylim=(-y_half, y_half) if y_half else None,
        subtitle=run_subtitle(window, run.run_name))
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plot_annotated_rate_comparison(
        change_rate, cs_series, run, estimator=estimator, domains=HATTERAS_DOMAINS,
        annotations=HATTERAS_ANNOTATIONS, lowess_config=lowess_config, config=config,
        save_path=str(save_path), show=show)
    if y_half:
        ax.yaxis.set_major_locator(MultipleLocator(y_half / Y_TICKS_PER_HALF))
    restyle(fig, ax, config.subtitle)
    fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close("all")
    print(f"  Saved rate plot: {save_path}")
    return save_path
