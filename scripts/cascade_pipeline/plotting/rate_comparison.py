"""Modeled-vs-CoastSat shoreline-change-rate figures.

Two entry points:
  plot_rate_comparison           the working REAL-domains/ALL-domains QC
                                  figure (toggle via real_domains_only)
  plot_annotated_rate_comparison the publication/poster figure with the
                                  full geographic annotation layer

Both consume cs_series from cascade_pipeline.coastsat_loess.build_coastsat_series
-- this module only renders, it doesn't load or smooth CoastSat data itself.

STYLE. Drawn under the house standard (scripts/hat_figure_style.py): printed
width, 8-9 pt type, one alongshore axis label, village bands. The one place
these figures depart from it is the TITLE. Every other figure in the project
moves its title sentence into a CAPTIONS.md beside the image; these two are
written automatically into a run folder and opened months later with nothing
around them, so the run's identity -- period, scope, background erosion, run
name -- stays on the canvas. It is set small and muted beside the subject
rather than as a 14 pt banner, but it is not moved off the image and there is
no captions file in a run folder. The wave height and the SLR rate came OFF
that line on 2026-09-10 (Hannah): they crowded it, and both are in the run's
metadata JSON and TXT beside the PNG.
"""

import dataclasses

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.lines import Line2D

from cascade_pipeline.annotations import (
    DEFAULT_ANNOTATIONS,
    add_geographic_annotations,
    annotation_legend_handles,
)
from cascade_pipeline.coastsat_loess import DEFAULT_LOESS, splice_loess_with_raw_south
from cascade_pipeline.domains import DEFAULT_DOMAINS

# `scripts/` is on sys.path already -- cascade_pipeline lives inside it, so
# importing this package at all means the style module is importable too.
from hat_figure_style import (
    C,
    DOMAIN_AXIS_LABEL,
    INK_MUTED,
    apply_style,
    figsize,
    open_frame,
    town_bands,
)

apply_style()


# COLOUR HERE IS NOT THE HOUSE SEMANTIC PAIR, ON PURPOSE (Hannah, 2026-09-10).
# This is the figure read most often in the project and it has one reading:
# the ORANGE model curve against the BLUE observation. Orange is the site
# config's -- HATTERAS_ANNOTATIONS.model_color -- and this module reads it off
# the `annotations` object every call site already passes, so the hex is
# defined in exactly one place and stays there. The blues below are the
# observed layer: the darker for the widest LOESS window, the lighter for a
# narrower one, and a muted blue for the per-transect cloud underneath.
# The red/blue VINTAGE pair is a different case -- two periods drawn together
# -- and deliberately does not appear here.


@dataclasses.dataclass(frozen=True)
class RateComparisonConfig:
    """Styling for the modeled-vs-CoastSat shoreline-change-rate figures.

    Attributes:
        window_colors: {window_domains: color} for each LOESS window.
        window_color_default: Fallback color for an unlisted window size.
        window_styles: [(linewidth, linestyle, alpha_factor), ...] matched
            to loess_config.window_domains by position; extra windows fall
            back to (1.5, "-", 0.80).
        raw_color: Color for the individual-transect scatter.
            The modelled curve's colour is NOT here: it is the site config's
            `AnnotationConfig.model_color` (HATTERAS_ANNOTATIONS defines the
            orange), which is the single authority for it.
        plot_raw_lrr: Show the transect scatter at all.
        raw_lrr_southern_only: True -> scatter only for the domains where
            LOESS is suppressed (D1-loess_config.skip_southern_domains).
            False -> scatter for every real domain. Only read by
            plot_annotated_rate_comparison.
        plot_reference_period: Also show the CoastSat period that doesn't
            match the run's start year (faded).
        raw_scatter_size: Marker area, points^2.
        raw_scatter_alpha: Opacity for the active period; the reference
            period (if shown) uses 0.35x this.
        domain_tick_step: X-axis tick spacing, in GIS domains.
    """

    window_colors: dict = dataclasses.field(
        default_factory=lambda: {7: "#6BAED6", 10: "#08519C"}
    )
    window_color_default: str = "#4A7C8E"
    window_styles: tuple = ((1.6, "-", 1.00), (1.8, "-", 1.00))
    raw_color: str = "#5BA3C9"
    plot_raw_lrr: bool = True
    raw_lrr_southern_only: bool = True
    plot_reference_period: bool = False
    raw_scatter_size: float = 6
    raw_scatter_alpha: float = 0.60
    domain_tick_step: int = 5


DEFAULT_RATE_COMPARISON = RateComparisonConfig()


def plot_coastsat_overlay(ax, cs_series, loess_config, config, x_transform, gis_x_transform=None):
    """Draw CoastSat raw-transect scatter + LOESS lines for one axis.

    Public because the sensitivity figures draw many model curves on one axis
    and need the SAME observed layer underneath them; a second implementation
    there would let the two drift apart on styling and on the southern splice.

    Shared by both branches of plot_rate_comparison (REAL vs ALL domains);
    the annotated figure has extra features (fill_between, legend handle
    tracking) and builds its overlay separately.

    Args:
        ax: Axes to draw on.
        cs_series: Output of build_coastsat_series.
        loess_config: LoessConfig (only window_domains is read here, to
            find the widest window).
        config: RateComparisonConfig.
        x_transform: along_coast_m array -> x-axis coordinate array, for
            the raw scatter.
        gis_x_transform: gis-domain-ID array -> x-axis coordinate array,
            for the LOESS lines. Defaults to identity (GIS-ID x-axis).
    """
    if gis_x_transform is None:
        gis_x_transform = lambda gis_x: gis_x

    widest_win = max(loess_config.window_domains)
    for cs in cs_series:
        is_active = cs["active"]
        if not is_active and not config.plot_reference_period:
            continue
        if config.plot_raw_lrr:
            x = x_transform(cs["transect_along_coast"])
            raw_alpha = config.raw_scatter_alpha if is_active else config.raw_scatter_alpha * 0.35
            ax.scatter(x, cs["transect_rates"], color=config.raw_color,
                       s=config.raw_scatter_size, alpha=raw_alpha, zorder=1, linewidths=0)
        for idx, win in enumerate(cs["windows"]):
            cs_color = config.window_colors.get(win["window"], config.window_color_default)
            lw_base, ls, alpha_factor = (
                config.window_styles[idx] if idx < len(config.window_styles) else (1.5, "-", 0.80)
            )
            is_widest = (win["window"] == widest_win)
            plot_gis_x, plot_y = splice_loess_with_raw_south(
                win["gis_x"], win["smoothed"], cs["transect_domains"], cs["transect_rates"],
                skip_n=loess_config.skip_southern_domains, is_widest_window=is_widest,
            )
            plot_x = gis_x_transform(plot_gis_x)
            lbl = f"{cs['label']} — {win['window']}-domain LOESS"
            if is_active:
                ax.plot(plot_x, plot_y, color=cs_color, lw=lw_base, ls=ls,
                        alpha=alpha_factor, zorder=4, label=lbl)
            else:
                ax.plot(plot_x, plot_y, color=cs_color, lw=lw_base * 0.85, ls=ls,
                        alpha=0.40 * alpha_factor, zorder=3, label=lbl + " (ref)")


ESTIMATOR_LABELS = {
    "lrr": "LRR",
    "endpoint": "endpoint difference",
}


def _rate_axis_label(estimator, title_case=False):
    """Axis label naming the estimator the modelled curve was built with.

    The observed curve is always an LRR -- CoastSat supplies a per-transect
    OLS slope -- so a figure that does not say which estimator the MODEL
    used is inviting the reader to assume they match. They only match when
    estimator is "lrr".

    Args:
        estimator: "lrr", "endpoint", or None to leave the estimator
            unnamed (the pre-2026-08-22 label, kept so an old figure can
            be redrawn unchanged).
        title_case: Match the annotated figure's title-case axis labels.

    Returns:
        The y-axis label string.
    """
    stem = ("Shoreline Change Rate" if title_case
            else "Shoreline change rate")
    if estimator is None:
        return f"{stem} (m/yr)"
    if estimator not in ESTIMATOR_LABELS:
        raise ValueError(
            f"unknown estimator {estimator!r}; expected one of "
            f"{sorted(ESTIMATOR_LABELS)} or None")
    return f"{stem}, {ESTIMATOR_LABELS[estimator]} (m/yr)"


def _run_parameters(run, scope, domains, endpoints=True):
    """The run's identity, as two short lines for the provenance title.

    What a reader needs in order to know WHICH run they are looking at, months
    later, out of a run folder: the scope of the figure, whether background
    erosion was on, and the run name. Plus the alongshore endpoints, which the
    axis label no longer carries.

    NOT the wave height and NOT the SLR rate (Hannah, 2026-09-10): both were
    crowding the line, and both are in the run's metadata JSON and TXT beside
    the figure, so nothing is lost by leaving them off the canvas. The period
    belongs to the subject line, not here.

    `endpoints` off for a figure that already names them at the axis corners.
    """
    bits = [scope, f"background erosion "
                   f"{'on' if run.background_erosion_on else 'off'}"]
    tail = f"run {run.run_name}"
    if endpoints:
        tail = (f"domain {domains.first_gis_id} Cape Point, "
                f"domain {domains.last_gis_id} Pea Island  ·  {tail}")
    return "  ·  ".join(bits) + "\n" + tail


PROVENANCE_SIZE = 7.5


def _provenance(ax, subject, params, extra_pad=0.0):
    """Subject on the title line, the run's identity stacked under it.

    The house rule is that nothing on the canvas belongs in a caption; the
    documented exception is a per-run artefact, which has no caption file to
    be read beside (see the module docstring). Small and muted, not a banner.

    STACKED, not set to the right of the subject: the run name alone is 36
    characters, so a right-aligned block ran straight through the title at the
    printed width (seen on the first render, 2026-09-10). The title pad is
    computed from the number of lines so constrained_layout reserves the room
    and the two never touch.
    """
    lines = str(params).count("\n") + 1
    pad = 4.0 + PROVENANCE_SIZE * 1.45 * lines + float(extra_pad)
    ax.set_title(subject, loc="left", pad=pad)
    ax.annotate(params, xy=(0.0, 1.0), xycoords="axes fraction",
                xytext=(0, 3 + extra_pad), textcoords="offset points",
                ha="left", va="bottom", fontsize=PROVENANCE_SIZE,
                color=INK_MUTED, linespacing=1.45, annotation_clip=False)


def _tick_step(span, step, max_ticks=20):
    """`step`, widened until `span` carries no more than `max_ticks` labels.

    At the printed width a 120-domain axis ticked every 5 is a grey smear;
    the config's step is the floor, not the answer, on the widest layout.
    """
    step = max(int(step), 1)
    while span / step > max_ticks:
        step *= 2
    return step


def plot_rate_comparison(change_rate, cs_series, run, real_domains_only=True,
                          estimator=None,
                          domains=DEFAULT_DOMAINS,
                          annotations=DEFAULT_ANNOTATIONS,
                          loess_config=DEFAULT_LOESS,
                          config=DEFAULT_RATE_COMPARISON,
                          sea_level_rise_rate_m_yr=None,
                          save_path=None, show=False):
    """Modeled vs. observed shoreline-change-rate figure (REAL or ALL domains).

    Two layouts, chosen by real_domains_only:
      True  -> x-axis is GIS domain IDs (domains.first_gis_id to
               domains.last_gis_id) only. Community spans are drawn from
               annotations.town_spans, same source of truth as the
               annotated figure -- label text is whatever's stored there.
      False -> x-axis is all domains.total_domains padded indices, buffers
               shaded red, GIS IDs on a secondary top axis.

    Args:
        change_rate: 1-D array, length domains.total_domains, m/yr (from
            cascade_pipeline.shoreline.compute_change_rate).
        cs_series: Output of build_coastsat_series.
        run: RunInfo for this run.
        real_domains_only: Selects the layout (see above).
        estimator: Which estimator built change_rate -- "lrr", "endpoint",
            or None to leave it unnamed. Names it on the y axis, since
            the observed curve is always an LRR and a figure that does
            not say invites the reader to assume the two match.
        sea_level_rise_rate_m_yr: Accepted and ignored by the figure
            since 2026-09-10 -- the SLR rate crowded the provenance
            line and is in the run's metadata beside the PNG. Kept in
            the signature so no caller has to change.
        save_path: If given, fig.savefig(save_path, dpi=300, bbox_inches="tight").
        show: Call plt.show() before returning.

    Returns:
        (fig, ax, fig_suffix): fig_suffix is "REAL_DOMAINS_ONLY" or
        "ALL_DOMAINS_WITH_BUFFERS", handy for building the output filename.
    """
    subject = (f"Modelled against {annotations.obs_source_name} shoreline "
               f"change rate, {run.start_year}–{run.end_year}")
    model_label = (f"{run.model_name}" if run.Hs is None
                   else f"{run.model_name}, Hs {run.Hs} m")

    if real_domains_only:
        gis_ids = np.arange(domains.first_gis_id, domains.last_gis_id + 1)
        fig, ax = plt.subplots(figsize=figsize("double", aspect=0.48),
                               constrained_layout=True)

        real_rate = change_rate[domains.start_real_index:domains.end_real_index]
        ax.plot(gis_ids, real_rate, color=annotations.model_color,
                linewidth=1.8, label=model_label, zorder=6)

        plot_coastsat_overlay(
            ax, cs_series, loess_config, config,
            x_transform=lambda along_m: along_m / domains.domain_spacing_m + domains.first_gis_id,
        )

        ax.axhline(0.0, linestyle=(0, (4, 3)), linewidth=0.8, color=INK_MUTED,
                   zorder=2)

        ax.set_xlim(domains.first_gis_id - 0.5, domains.last_gis_id + 0.5)
        step = _tick_step(domains.num_real_domains, config.domain_tick_step)
        ax.set_xticks(np.arange(domains.first_gis_id,
                                domains.last_gis_id + 1, step))
        # After the limits, so a span outside the view is skipped and a label
        # is clamped to the visible part of its span.
        town_bands(ax, spans=annotations.town_spans)

        ax.set_xlabel(DOMAIN_AXIS_LABEL)
        ax.set_ylabel(_rate_axis_label(estimator))
        _provenance(ax, subject,
                    _run_parameters(run, "real domains only", domains))
        ax.grid(axis="y")
        open_frame(ax)
        ax.set_axisbelow(True)
        fig.legend(loc="outside lower center", ncol=3, frameon=False)

        fig_suffix = "REAL_DOMAINS_ONLY"

    else:
        domain_numbers = np.arange(domains.total_domains)
        fig, ax = plt.subplots(figsize=figsize("double", aspect=0.46),
                               constrained_layout=True)

        ax.axvspan(0, domains.start_real_index - 0.5, color=C["BASE_FILL"],
                   zorder=0, lw=0, label="buffer domains")
        ax.axvspan(domains.end_real_index - 0.5, domains.total_domains - 1,
                   color=C["BASE_FILL"], zorder=0, lw=0)

        ax.plot(domain_numbers, change_rate, color=annotations.model_color,
                linewidth=1.8, label=model_label, zorder=6)

        plot_coastsat_overlay(
            ax, cs_series, loess_config, config,
            x_transform=lambda along_m: along_m / domains.domain_spacing_m + domains.start_real_index,
            gis_x_transform=lambda gis_x: gis_x - domains.first_gis_id + domains.start_real_index,
        )

        for edge in (domains.start_real_index, domains.end_real_index - 1):
            ax.axvline(edge, linestyle=(0, (4, 3)), linewidth=0.8,
                       color=INK_MUTED, zorder=2)
        ax.axhline(0.0, linestyle=(0, (4, 3)), linewidth=0.8, color=INK_MUTED,
                   zorder=2)

        ax.set_xlim(-0.5, domains.total_domains - 0.5)
        step = _tick_step(domains.total_domains, config.domain_tick_step)
        ax.set_xticks(np.arange(0, domains.total_domains, step))
        town_bands(ax, spans={
            name: (domains.gis_to_pad(gis_lo), domains.gis_to_pad(gis_hi))
            for name, (gis_lo, gis_hi) in annotations.town_spans.items()})
        ax.set_xlabel(f"{run.model_name} domain index, buffers included "
                      f"(0\u2013{domains.total_domains - 1})")

        top_ax = ax.secondary_xaxis("top")
        top_positions, top_labels = [], []
        for gis_id in range(domains.first_gis_id, domains.last_gis_id + 1, step):
            top_positions.append(domains.start_real_index + (gis_id - domains.first_gis_id))
            top_labels.append(str(gis_id))
        top_ax.set_xticks(top_positions)
        top_ax.set_xticklabels(top_labels)
        top_ax.set_xlabel(DOMAIN_AXIS_LABEL)

        ax.set_ylabel(_rate_axis_label(estimator))
        # Extra pad: the title row sits above the secondary GIS axis, which
        # the parent axes' title placement knows nothing about.
        _provenance(ax, subject,
                    _run_parameters(run, "all domains, buffers included",
                                    domains),
                    extra_pad=26)
        ax.grid(axis="y")
        open_frame(ax)
        ax.set_axisbelow(True)
        fig.legend(loc="outside lower center", ncol=3, frameon=False)

        fig_suffix = "ALL_DOMAINS_WITH_BUFFERS"

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  Saved plot: {save_path}")
    if show:
        plt.show()

    return fig, ax, fig_suffix


def plot_annotated_rate_comparison(change_rate, cs_series, run,
                                    estimator=None,
                                    domains=DEFAULT_DOMAINS,
                                    annotations=DEFAULT_ANNOTATIONS,
                                    loess_config=DEFAULT_LOESS,
                                    config=DEFAULT_RATE_COMPARISON,
                                    sea_level_rise_rate_m_yr=None,
                                    save_path=None, show=False):
    """Publication/poster figure: modeled rate + full geographic annotation layer.

    Always uses the real-domains-only (GIS first_gis_id-last_gis_id) x-axis,
    regardless of the REAL/ALL toggle used for plot_rate_comparison -- this
    figure is for sharing, not QC.

    Args:
        change_rate: 1-D array, length domains.total_domains, m/yr.
        cs_series: Output of build_coastsat_series.
        run: RunInfo for this run.
        estimator: Which estimator built change_rate -- "lrr",
            "endpoint", or None to leave it unnamed. See
            plot_rate_comparison.
        sea_level_rise_rate_m_yr: Accepted and ignored by the figure
            since 2026-09-10 -- the SLR rate crowded the provenance
            line and is in the run's metadata beside the PNG. Kept in
            the signature so no caller has to change.
        save_path: If given, fig.savefig(save_path, dpi=300,
            bbox_inches="tight", facecolor="white").
        show: Call plt.show() before returning.

    Returns:
        (fig, ax)
    """
    gis_ids = np.arange(domains.first_gis_id, domains.last_gis_id + 1)
    real_rate = change_rate[domains.start_real_index:domains.end_real_index]

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.62),
                           constrained_layout=True)

    # Geographic annotations drawn first so data renders on top.
    add_geographic_annotations(ax, annotations)

    data_handles = []
    widest_win = max(loess_config.window_domains)
    for cs in cs_series:
        is_active = cs["active"]
        if not is_active and not config.plot_reference_period:
            continue

        scatter_x = cs["transect_along_coast"] / domains.domain_spacing_m + domains.first_gis_id
        if config.plot_raw_lrr:
            if config.raw_lrr_southern_only:
                south_mask = cs["transect_domains"] <= loess_config.skip_southern_domains
                scatter_x_plot = scatter_x[south_mask]
                scatter_y_plot = cs["transect_rates"][south_mask]
                raw_lbl = (f"{cs['label']} — transect LRR "
                           f"(D{domains.first_gis_id}-{loess_config.skip_southern_domains})"
                           if is_active else None)
            else:
                scatter_x_plot = scatter_x
                scatter_y_plot = cs["transect_rates"]
                raw_lbl = f"{cs['label']} — transect LRR" if is_active else None
            raw_alpha = config.raw_scatter_alpha if is_active else config.raw_scatter_alpha * 0.35
            ax.scatter(scatter_x_plot, scatter_y_plot, color=config.raw_color,
                       s=config.raw_scatter_size, alpha=raw_alpha, zorder=1,
                       linewidths=0, label=raw_lbl)
            if is_active:
                data_handles.append(
                    Line2D([0], [0], color=config.raw_color, marker=".", ms=5,
                           ls="none", alpha=config.raw_scatter_alpha, label=raw_lbl)
                )

        for idx, win in enumerate(cs["windows"]):
            cs_color = config.window_colors.get(win["window"], config.window_color_default)
            lw_base, ls, alpha_factor = (
                config.window_styles[idx] if idx < len(config.window_styles) else (1.5, "-", 0.80)
            )
            is_widest = (win["window"] == widest_win)
            gis_x, rate = splice_loess_with_raw_south(
                win["gis_x"], win["smoothed"], cs["transect_domains"], cs["transect_rates"],
                skip_n=loess_config.skip_southern_domains, is_widest_window=is_widest,
            )
            lbl = f"{cs['label']} — {win['window']}-domain LOESS"
            if is_active:
                if is_widest:
                    ax.fill_between(gis_x, rate, 0, where=(rate < 0), alpha=0.14,
                                    color=cs_color, interpolate=True)
                    ax.fill_between(gis_x, rate, 0, where=(rate >= 0), alpha=0.10,
                                    color=cs_color, interpolate=True)
                ax.plot(gis_x, rate, color=cs_color, lw=lw_base, ls=ls,
                        alpha=alpha_factor, zorder=5, label=lbl)
                data_handles.append(
                    Line2D([0], [0], color=cs_color, lw=lw_base, ls=ls,
                           alpha=alpha_factor, label=lbl)
                )
            else:
                ax.plot(gis_x, rate, color=cs_color, lw=lw_base * 0.85, ls=ls,
                        alpha=0.40 * alpha_factor, zorder=4)
                data_handles.append(
                    Line2D([0], [0], color=cs_color, lw=lw_base * 0.85, ls=ls,
                           alpha=0.40 * alpha_factor, label=lbl + " (ref)")
                )

    model_label = (f"{run.model_name}" if run.Hs is None
                   else f"{run.model_name}, Hs {run.Hs} m")
    ax.plot(gis_ids, real_rate, color=annotations.model_color, linewidth=2.0,
            zorder=6, label=model_label)
    data_handles.append(
        Line2D([0], [0], color=annotations.model_color, lw=2.0,
               label=model_label)
    )

    ax.axhline(0, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 3)),
               zorder=3)

    # Scatter/LOESS transition marker -- only when southern domains show
    # raw scatter only; marks where dots end and LOESS lines begin.
    if config.plot_raw_lrr and config.raw_lrr_southern_only and loess_config.skip_southern_domains > 0:
        ax.axvline(loess_config.skip_southern_domains + 0.5, color=INK_MUTED,
                   lw=0.8, ls=(0, (4, 4)), zorder=2)

    ax.set_xlim(domains.first_gis_id - 0.5, domains.last_gis_id + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="both", which="minor", length=1.8)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)

    # Lock ylim to data range, then place accretion/erosion side labels.
    all_vals = np.concatenate(
        [real_rate] + [w["smoothed"][np.isfinite(w["smoothed"])]
                       for cs in cs_series for w in cs["windows"]]
    )
    ymin, ymax = all_vals.min(), all_vals.max()
    ypad = (ymax - ymin) * 0.06
    ax.set_ylim(ymin - ypad, ymax + ypad)

    ybot, ytop = ax.get_ylim()
    zero_frac = (0 - ybot) / (ytop - ybot)
    acc_y = (annotations.label_accretion_y if annotations.label_accretion_y is not None
             else zero_frac + (1 - zero_frac) / 2)
    ero_y = (annotations.label_erosion_y if annotations.label_erosion_y is not None
             else zero_frac / 2)
    # A white backing: the erosion label sits at the right edge, which is where
    # the observed curve runs on this period.
    for _y, _txt in ((acc_y, "accretion \u25b2"), (ero_y, "erosion \u25bc")):
        ax.text(1.0, _y, _txt, transform=ax.transAxes, fontsize=8,
                color=INK_MUTED, ha="right", va="center", zorder=8,
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none",
                          boxstyle="square,pad=0.2"))

    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(_rate_axis_label(estimator, title_case=True))
    ax.text(0.0, 1.005, f"\u2190 {annotations.low_end_label}",
            transform=ax.transAxes, fontsize=7.5, color=INK_MUTED, ha="left",
            va="bottom", clip_on=False)
    ax.text(1.0, 1.005, f"{annotations.high_end_label} \u2192",
            transform=ax.transAxes, fontsize=7.5, color=INK_MUTED, ha="right",
            va="bottom", clip_on=False)

    # The run's identity, kept on the canvas -- see the module docstring. The
    # separate italic "Model | Obs | SLR | Run" footnote this figure used to
    # carry beneath the legend said the same things twice; it is folded in
    # here, where the eye already is.
    _provenance(
        ax,
        f"Modelled against {annotations.obs_source_name} shoreline change, "
        f"{run.start_year}–{run.end_year}",
        _run_parameters(run,
                        f"{annotations.obs_source_name} LRR per "
                        f"{int(domains.domain_spacing_m)} m domain",
                        domains, endpoints=False),
        extra_pad=12)

    fig.legend(
        handles=data_handles + annotation_legend_handles(annotations),
        loc="outside lower center", ncol=4, frameon=False,
    )

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved annotated plot: {save_path}")
    if show:
        plt.show()

    return fig, ax
