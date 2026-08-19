"""Modeled-vs-CoastSat shoreline-change-rate figures.

Two entry points:
  plot_rate_comparison           the working REAL-domains/ALL-domains QC
                                  figure (toggle via real_domains_only)
  plot_annotated_rate_comparison the publication/poster figure with the
                                  full geographic annotation layer

Both consume cs_series from cascade_pipeline.coastsat_loess.build_coastsat_series
-- this module only renders, it doesn't load or smooth CoastSat data itself.
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
    window_styles: tuple = ((1.8, "-", 1.00), (2.0, "-", 1.00))
    raw_color: str = "#5BA3C9"
    plot_raw_lrr: bool = True
    raw_lrr_southern_only: bool = True
    plot_reference_period: bool = False
    raw_scatter_size: float = 6
    raw_scatter_alpha: float = 0.60
    domain_tick_step: int = 5


DEFAULT_RATE_COMPARISON = RateComparisonConfig()


def _plot_coastsat_overlay(ax, cs_series, loess_config, config, x_transform, gis_x_transform=None):
    """Draw CoastSat raw-transect scatter + LOESS lines for one axis.

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
            lbl = f"{cs['label']} — LOESS {win['window']}-dom"
            if is_active:
                ax.plot(plot_x, plot_y, color=cs_color, lw=lw_base, ls=ls,
                        alpha=alpha_factor, zorder=4, label=lbl)
            else:
                ax.plot(plot_x, plot_y, color=cs_color, lw=lw_base * 0.85, ls=ls,
                        alpha=0.40 * alpha_factor, zorder=3, label=lbl + " (ref)")


def plot_rate_comparison(change_rate, cs_series, run, real_domains_only=True,
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
        sea_level_rise_rate_m_yr: Shown in the title if given.
        save_path: If given, fig.savefig(save_path, dpi=300, bbox_inches="tight").
        show: Call plt.show() before returning.

    Returns:
        (fig, ax, fig_suffix): fig_suffix is "REAL_DOMAINS_ONLY" or
        "ALL_DOMAINS_WITH_BUFFERS", handy for building the output filename.
    """
    be_label = "on" if run.background_erosion_on else "off"
    slr_txt = (f"SLR={sea_level_rise_rate_m_yr * 1000:.1f} mm/yr | "
               if sea_level_rise_rate_m_yr is not None else "")

    if real_domains_only:
        gis_ids = np.arange(domains.first_gis_id, domains.last_gis_id + 1)
        fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)

        real_rate = change_rate[domains.start_real_index:domains.end_real_index]
        ax.plot(gis_ids, real_rate, color=annotations.model_color, linewidth=2,
                label=f"Model Hs={run.Hs} m", zorder=6)

        _plot_coastsat_overlay(
            ax, cs_series, loess_config, config,
            x_transform=lambda along_m: along_m / domains.domain_spacing_m + domains.first_gis_id,
        )

        for span_label, (gis_lo, gis_hi) in annotations.town_spans.items():
            ax.axvspan(gis_lo - 0.5, gis_hi + 0.5, alpha=0.10, color="steelblue", zorder=0)
            ax.text((gis_lo + gis_hi) / 2, ax.get_ylim()[0], span_label,
                    ha="center", va="bottom", fontsize=7, style="italic", color="steelblue")

        ax.axhline(0.0, linestyle="--", linewidth=1, color="gray", alpha=0.7)

        xticks = np.arange(domains.first_gis_id, domains.last_gis_id + 1, config.domain_tick_step)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

        ax.set_xlabel(f"GIS Domain ID ({domains.first_gis_id}\u2013{domains.last_gis_id})")
        ax.set_ylabel("Shoreline change rate (m/yr)")
        ax.set_title(
            f"Modeled vs {annotations.obs_source_name} Shoreline Change Rate \u2013 "
            f"{annotations.region_name} | "
            f"Real domains only | {slr_txt}"
            f"{run.start_year}\u2013{run.end_year} | Hs={run.Hs} m | BE={be_label}"
        )
        ax.grid(alpha=0.3)
        ax.legend()

        fig_suffix = "REAL_DOMAINS_ONLY"

    else:
        domain_numbers = np.arange(domains.total_domains)
        fig, ax = plt.subplots(figsize=(22, 6), constrained_layout=True)

        ax.axvspan(0, domains.start_real_index - 0.5, alpha=0.12, color="red", label="Buffer")
        ax.axvspan(domains.end_real_index - 0.5, domains.total_domains - 1, alpha=0.12, color="red")

        for span_label, (gis_lo, gis_hi) in annotations.town_spans.items():
            pad_lo, pad_hi = domains.gis_to_pad(gis_lo), domains.gis_to_pad(gis_hi)
            ax.axvspan(pad_lo - 0.5, pad_hi + 0.5, alpha=0.10, color="steelblue", zorder=0)

        ax.plot(domain_numbers, change_rate, color=annotations.model_color,
                linewidth=2, label=f"Model Hs={run.Hs} m", zorder=6)

        _plot_coastsat_overlay(
            ax, cs_series, loess_config, config,
            x_transform=lambda along_m: along_m / domains.domain_spacing_m + domains.start_real_index,
            gis_x_transform=lambda gis_x: gis_x - domains.first_gis_id + domains.start_real_index,
        )

        ax.axvline(domains.start_real_index, linestyle="--", linewidth=1, color="k", alpha=0.5)
        ax.axvline(domains.end_real_index - 1, linestyle="--", linewidth=1, color="k", alpha=0.5)
        ax.axhline(0.0, linestyle="--", linewidth=1, color="gray", alpha=0.7)

        ax.text((0 + domains.start_real_index) / 2, 0, "Left\nbuffer",
                ha="center", va="center", fontsize=8, style="italic", alpha=0.6)
        ax.text((domains.start_real_index + domains.end_real_index) / 2, 0,
                f"Real island (GIS {domains.first_gis_id}\u2013{domains.last_gis_id})",
                ha="center", va="center", fontsize=9, fontweight="bold", alpha=0.5)
        ax.text((domains.end_real_index + domains.total_domains) / 2, 0, "Right\nbuffer",
                ha="center", va="center", fontsize=8, style="italic", alpha=0.6)

        xticks = np.arange(0, domains.total_domains, config.domain_tick_step)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)
        ax.set_xlabel(f"{run.model_name} domain index (including buffers, 0\u2013{domains.total_domains - 1})")

        top_ax = ax.secondary_xaxis("top")
        top_positions, top_labels = [], []
        for gis_id in range(domains.first_gis_id, domains.last_gis_id + 1, config.domain_tick_step):
            top_positions.append(domains.start_real_index + (gis_id - domains.first_gis_id))
            top_labels.append(str(gis_id))
        top_ax.set_xticks(top_positions)
        top_ax.set_xticklabels(top_labels, fontsize=9)
        top_ax.set_xlabel(f"GIS Domain ID ({domains.first_gis_id}\u2013{domains.last_gis_id})")

        ax.set_ylabel("Shoreline change rate (m/yr)")
        ax.set_title(
            f"Modeled vs {annotations.obs_source_name} Shoreline Change Rate \u2013 "
            f"{annotations.region_name} | "
            f"All domains (buffers included) | {slr_txt}"
            f"{run.start_year}\u2013{run.end_year} | Hs={run.Hs} m | BE={be_label}"
        )
        ax.grid(alpha=0.3)
        ax.legend()

        fig_suffix = "ALL_DOMAINS_WITH_BUFFERS"

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  Saved plot: {save_path}")
    if show:
        plt.show()

    return fig, ax, fig_suffix


def plot_annotated_rate_comparison(change_rate, cs_series, run,
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
        sea_level_rise_rate_m_yr: Shown in the caption if given.
        save_path: If given, fig.savefig(save_path, dpi=300,
            bbox_inches="tight", facecolor="white").
        show: Call plt.show() before returning.

    Returns:
        (fig, ax)
    """
    gis_ids = np.arange(domains.first_gis_id, domains.last_gis_id + 1)
    real_rate = change_rate[domains.start_real_index:domains.end_real_index]

    fig, ax = plt.subplots(figsize=(13, 7.0), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

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
            lbl = f"{cs['label']} — LOESS {win['window']}-dom"
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

    ax.plot(gis_ids, real_rate, color=annotations.model_color, linewidth=2.4,
            zorder=6, label=f"Model Hs={run.Hs} m")
    data_handles.append(
        Line2D([0], [0], color=annotations.model_color, lw=2.4, label=f"Model Hs={run.Hs} m")
    )

    ax.axhline(0, color="#2c2c2c", linewidth=1.2, linestyle="--", alpha=0.65, zorder=3)

    # Scatter/LOESS transition marker -- only when southern domains show
    # raw scatter only; marks where dots end and LOESS lines begin.
    if config.plot_raw_lrr and config.raw_lrr_southern_only and loess_config.skip_southern_domains > 0:
        ax.axvline(loess_config.skip_southern_domains + 0.5, color="0.60", lw=0.8,
                   ls=(0, (4, 4)), alpha=0.55, zorder=2)

    ax.set_xlim(domains.first_gis_id - 0.5, domains.last_gis_id + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="both", which="major", labelsize=11, direction="in", length=5)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)
    ax.grid(True, which="major", linestyle=":", linewidth=0.6, alpha=0.4, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.2)

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
    ax.text(1.0, acc_y, "Accretion \u25b2", transform=ax.transAxes,
            fontsize=9.5, color="#555555", ha="right", va="center", style="italic")
    ax.text(1.0, ero_y, "Erosion \u25bc", transform=ax.transAxes,
            fontsize=9.5, color="#555555", ha="right", va="center", style="italic")

    ax.set_xlabel(f"{run.model_name} Model Domain ({int(domains.domain_spacing_m)} m alongshore)",
                  fontsize=12, fontweight="bold", labelpad=4)
    ax.set_ylabel("Shoreline Change Rate (m/yr)", fontsize=12, fontweight="bold", labelpad=8)
    ax.text(0.0, 1.01, f"\u2190 {annotations.low_end_label}", transform=ax.transAxes, fontsize=9,
            color="#444444", ha="left", va="bottom", style="italic", clip_on=False)
    ax.text(1.0, 1.01, f"{annotations.high_end_label} \u2192", transform=ax.transAxes, fontsize=9,
            color="#444444", ha="right", va="bottom", style="italic", clip_on=False)

    be_label = "on" if run.background_erosion_on else "off"
    ax.set_title(
        f"Modeled vs Observed Shoreline Change \u2014 {annotations.region_name}  |  "
        f"{run.start_year}\u2013{run.end_year}  |  Hs={run.Hs} m  |  BE={be_label}",
        fontsize=12, fontweight="bold", pad=12, color="#1a2a3a"
    )

    ax.legend(
        handles=data_handles + annotation_legend_handles(annotations),
        loc="upper center", bbox_to_anchor=(0.5, -0.10), bbox_transform=ax.transAxes,
        fontsize=9.0, framealpha=0.95, edgecolor="#cccccc", frameon=True, ncol=5,
    )

    slr_txt = (f"SLR={sea_level_rise_rate_m_yr * 1000:.1f} mm/yr  |  "
               if sea_level_rise_rate_m_yr is not None else "")
    ax.annotate(
        f"Model: {run.model_name}  |  Obs: {annotations.obs_source_name} LRR per "
        f"{int(domains.domain_spacing_m)}-m domain  |  "
        f"{slr_txt}Run: {run.run_name}",
        xy=(0, 0), xycoords="axes fraction", xytext=(0, -0.22), textcoords="axes fraction",
        fontsize=7.5, color="#666666", ha="left", va="top", style="italic", annotation_clip=False,
    )

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved annotated plot: {save_path}")
    if show:
        plt.show()

    return fig, ax
