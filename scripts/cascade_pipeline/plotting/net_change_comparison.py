"""
A run's end-minus-start shoreline against the CoastSat net change, in the layout of the annotated rate figure.

The observed side is the DEM-to-DEM target (coastsat/net_change/<window>): end-window
mean minus start-window mean per transect, its domain means, and its 7-domain LOWESS.
The model is smoothed the same way (7-domain LOWESS, GIS 1-10 raw, as every skill score
smooths it) and drawn over its unsmoothed domain values. A window with no net-change
target draws nothing. Called by the runner after the rate figures and by
rerender_run_figures.py.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-06
"""
from __future__ import annotations

import dataclasses
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

from cascade_pipeline.coastsat_lowess import CoastSatDataset, LowessConfig, build_coastsat_series
from cascade_pipeline.plotting.rate_comparison import (
    DEFAULT_RATE_COMPARISON, coastsat_domain_mean, plot_annotated_rate_comparison)
from cascade_pipeline.shoreline import compute_change_rate
from site_layer.hat_figure_style import mark_offaxis, offaxis_clause, record_caption
from site_layer.hat_observed_rates import (
    NET_CHANGE_CENTRES, NET_CHANGE_WINDOWS, WINDOW_ROLE, net_change_dir)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS

# --- CONFIG ------------------------------------------------------------------
LOWESS_CONFIG = LowessConfig(window_domains=(7,), skip_southern_domains=10)
COASTSAT_END = "2026-01-13"
# One symmetric y axis for every run, so figures compare (the project's position-change range)
Y_HALF = 100.0
PRESET_TEXT = {"zeroBE": "no source/sink",
               "edgeBE": "source/sink at the end domains only (GIS 1 and 90)",
               "calibBE": "zone source/sink (calibBE)",
               "domainBE": "per-domain source/sink (set 1)"}
MANAGEMENT_TOKENS = (("road", "road maintenance"), ("bdm", "beach and dune management"),
                     ("nourish", "nourishment"), ("reloc", "road relocation"))
# -----------------------------------------------------------------------------


def has_target(start_year, end_year):
    return (int(start_year), int(end_year)) in NET_CHANGE_WINDOWS


# The model smoothed as every score smooths it: 7-domain LOWESS, GIS 1-10 raw; the last domain raw too,
# because the end solve matches that domain's own value (Hannah, 2026-10-06)
def smooth_model(position_m, domains=HATTERAS_DOMAINS, lowess_config=LOWESS_CONFIG):
    real = np.asarray(position_m, dtype=float)[domains.start_real_index:domains.end_real_index]
    gis = np.arange(domains.first_gis_id, domains.last_gis_id + 1, dtype=float)
    sm = lowess(real, gis, frac=max(lowess_config.window_domains) / len(gis),
                return_sorted=False)
    raw = (gis <= lowess_config.skip_southern_domains) | (gis == domains.last_gis_id)
    sm[raw] = real[raw]
    out = np.array(position_m, dtype=float)
    out[domains.start_real_index:domains.end_real_index] = sm
    return out


# Source/sink, groin and management in plain words, from the run name's tokens
def setup_text(run_name):
    tokens = run_name.split("_")
    preset = tokens[3] if len(tokens) > 3 else ""
    groin = ("blocking groin" if "groinblock" in tokens
             else "groin" if "groin" in tokens else "no groin")
    management = [text for token, text in MANAGEMENT_TOKENS if token in tokens]
    return PRESET_TEXT.get(preset, preset), groin, management


# Triangles at the edge for everything beyond +/-Y_HALF; returns the caption clause naming it
def _mark_offaxis(ax, smoothed, raw, cs_series, domains=HATTERAS_DOMAINS):
    D, model_c = domains, HATTERAS_ANNOTATIONS.model_color
    gis = np.arange(D.first_gis_id, D.last_gis_id + 1)
    real = slice(D.start_real_index, D.end_real_index)
    named = [("CASCADE (smoothed)", mark_offaxis(ax, gis, np.asarray(smoothed)[real], Y_HALF, color=model_c)),
             ("CASCADE (unsmoothed)", mark_offaxis(ax, gis, np.asarray(raw)[real], Y_HALF, color=model_c))]
    clause, dots = "", ""
    for cs in cs_series:
        if not cs["active"]:
            continue
        mx, my = coastsat_domain_mean(cs)
        named.append(("CoastSat domain mean", mark_offaxis(ax, mx, my, Y_HALF, color="#6BAED6")))
        for win in cs["windows"]:
            named.append(("CoastSat LOWESS", mark_offaxis(ax, win["gis_x"], win["smoothed"], Y_HALF,
                                                          color="#6BAED6")))
        south = cs["transect_domains"] <= LOWESS_CONFIG.skip_southern_domains
        y = np.asarray(cs["transect_rates"])[south]
        x = cs["transect_along_coast"][south] / D.domain_spacing_m + D.first_gis_id
        out = np.isfinite(y) & (np.abs(y) > Y_HALF)
        if out.any():
            ax.scatter(x[out], np.sign(y[out]) * Y_HALF, s=8, marker="^", c="#5BA3C9",
                       zorder=14, clip_on=False)
            g = np.asarray(cs["transect_domains"])[south][out]
            dots = (f" {int(out.sum())} transects at GIS {g.min()}–{g.max()} are off the axis too "
                    f"(up to {y[out][np.argmax(np.abs(y[out]))]:+.0f} m), marked at the edge.")
    clause = offaxis_clause([(k, v) for k, v in named if v], Y_HALF)
    return clause + dots


def plot_net_change_comparison(shoreline_m, run, save_path, wave_climate=None,
                               flip_sign=True, show=False):
    """Draw the figure for one run; returns the path, or None when the window has no target.

    Args:
        shoreline_m: The run's annual shoreline matrix (states x padded domains).
        run: RunInfo; its start_year/end_year pick the target window.
        save_path: Where the PNG goes (run_layout kind "figure_net_change").
        wave_climate: The run's wave settings as one line, for the caption.
        flip_sign: The runner's FLIP_SIGN_MODEL (x_s grows landward).
    """
    window = (int(run.start_year), int(run.end_year))
    if not has_target(*window):
        print(f"  net-change figure: no target for {window[0]}-{window[1]}; skipped")
        return None
    tag = f"{window[0]}_{window[1]}"
    dataset = CoastSatDataset(
        label=f"CoastSat net change ({window[0]}–{window[1]})", period_start=window[0],
        csv_path=str(net_change_dir(*window) / f"transect_net_change_{tag}.csv"),
        rate_col="net_change_m")
    cs_series = build_coastsat_series([dataset], active_period_start=window[0],
                                      lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)
    position_m = compute_change_rate(shoreline_m, span_years=1, flip_sign=flip_sign)

    source_sink, groin, management = setup_text(run.run_name)
    role = WINDOW_ROLE.get(window, "").replace(" (held out)", "")
    config = dataclasses.replace(
        DEFAULT_RATE_COMPARISON, quantity="position", publication_text=True,
        show_features=True, plot_domain_means=True, raw_lrr_southern_only=True,
        observed_label="CoastSat (7-domain LOWESS)", model_label="CASCADE (7-domain LOWESS)",
        title=f"Modelled and observed shoreline position change, {window[0]}–{window[1]}",
        ylim=(-Y_HALF, Y_HALF),
        subtitle=" · ".join(p for p in (role, source_sink, groin) if p))
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    smoothed = smooth_model(position_m)
    fig, ax = plot_annotated_rate_comparison(
        smoothed, cs_series, run, estimator=None, domains=HATTERAS_DOMAINS,
        annotations=HATTERAS_ANNOTATIONS, lowess_config=LOWESS_CONFIG, config=config,
        save_path=str(save_path), show=show, change_rate_raw=position_m)
    offaxis = _mark_offaxis(ax, smoothed, position_m, cs_series)
    fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close("all")

    (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[window]
    c0, c1 = NET_CHANGE_CENTRES[window]
    years = np.asarray(shoreline_m).shape[0] - 1
    end_note = (f" CoastSat ends on {COASTSAT_END}, so that mean covers {e0} to {COASTSAT_END}."
                if e1 > COASTSAT_END else "")
    managed = ", ".join(management) if management else "no management"
    record_caption(save_path, (
        f"Shoreline position change along Hatteras Island, {window[0]}–{window[1]} "
        f"({role.lower() or 'model period'}), from GIS domain 1 (Cape Point) to 90 (Pea Island); "
        f"positive is seaward. Modelled (CASCADE, orange): the shoreline at the end of the "
        f"{years}-year run (1 Jan {window[0] + years}) minus the start, with {source_sink}, "
        f"{groin} and {managed}; thick: smoothed with a 7-domain LOWESS north of domain 10 "
        f"(unsmoothed south of it and at domain 90, the end value the source/sink solve matches); thin: each domain's unsmoothed value. "
        f"The start shoreline is the CoastSat mean over {s0} to {s1}. Observed (CoastSat, blue): "
        f"per transect, the mean shoreline position over {e0} to {e1} (centred on {c1}) minus "
        f"the mean over {s0} to {s1} (centred on {c0}).{end_note} Thin line: the mean of each "
        f"500 m domain's transects. Thick line: the transect values smoothed with a 7-domain "
        f"LOWESS north of domain 10; south of it the individual transects are shown as dots. "
        f"Blue shading: the smoothed observed change against zero. Shaded bands: communities; "
        f"dashed lines: village centres; hatched: shoal zones; dash-dot: piers; dotted: Buxton "
        f"groin. The y axis is fixed at ±{Y_HALF:g} m on every run so figures compare."
        + offaxis + (f" Wave climate: {wave_climate}." if wave_climate else "")
        + f" Run {run.run_name}."))
    print(f"  Saved net-change plot: {save_path}")
    return save_path
