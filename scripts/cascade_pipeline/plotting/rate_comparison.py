"""
Modelled against CoastSat shoreline-change-rate figures.

    from cascade_pipeline.plotting.rate_comparison import plot_rate_comparison, plot_annotated_rate_comparison

The QC figure and the annotated publication figure; renders only, CoastSat comes from
coastsat_lowess. Details: scripts/cascade_pipeline/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import dataclasses
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.lines import Line2D

from cascade_pipeline.annotations import (
    DEFAULT_ANNOTATIONS,
    add_geographic_annotations,
    annotation_legend_handles,
)
from cascade_pipeline.coastsat_lowess import (
    DEFAULT_LOWESS,
    compute_domain_means,
    splice_lowess_with_raw_south,
)
from cascade_pipeline.domains import DEFAULT_DOMAINS

# scripts/ is already on sys.path, so the style module imports with this package
from site_layer.hat_figure_style import (
    C,
    DOMAIN_AXIS_LABEL,
    INK,
    INK_MUTED,
    record_caption,
    apply_style,
    figsize,
    open_frame,
    town_bands,
)

apply_style()


# Colour: the orange model (annotations.model_color) against the blue observation, not the house pair


@dataclasses.dataclass(frozen=True)
class RateComparisonConfig:
    """Styling for the modeled-vs-CoastSat shoreline-change-rate figures.

    Attributes:
        window_colors: {window_domains: color} for each LOWESS window.
        window_color_default: Fallback color for an unlisted window size.
        window_styles: [(linewidth, linestyle, alpha_factor), ...] matched
            to lowess_config.window_domains by position; extra windows fall
            back to (1.5, "-", 0.80).
        raw_color: Color for the individual-transect scatter.
            The modelled curve's colour is NOT here: it is the site config's
            `AnnotationConfig.model_color` (HATTERAS_ANNOTATIONS defines the
            orange), which is the single authority for it.
        plot_domain_means: Draw the UNSMOOTHED per-domain mean of the
            transect rates under the LOWESS curve, the whole island (Hannah,
            2026-09-22). Without it the only observed curve on a model figure
            is the smoothed one, and a reader cannot tell a real alongshore
            feature from one the 5 km window flattened -- the observation
            figures beside these (rates_figures.lrr_figures) draw the same
            two curves, so the two products no longer disagree.
        domain_mean_lw, domain_mean_alpha: Its weight; thinner and fainter
            than the LOWESS it underlies, in the same colour as that window.
        plot_raw_lrr: Show the transect scatter at all.
        raw_lrr_southern_only: True -> scatter only for the domains where
            LOWESS is suppressed (D1-lowess_config.skip_southern_domains).
            False -> scatter for every real domain. Read by
            plot_annotated_rate_comparison always, and by
            plot_coastsat_overlay only when overlay_raw_southern_only is set.
        overlay_raw_southern_only: Make plot_coastsat_overlay (the real-
            domains and with-buffers run figures) honour
            raw_lrr_southern_only too. Off by default so the matrix figures
            keep the whole-island scatter; the sensitivity cells turn it on
            with plot_domain_means=False (Hannah, 2026-09-29: "only keeping
            the 7 domain LOWESS, just use the dots for domains 1-10").
        south_mean_line: Draw the unsmoothed per-domain mean over
            D1..skip_southern_domains as a DASHED line, so the observed curve
            runs the whole island but the stretch that is NOT a LOWESS looks
            different from the stretch that is (Hannah, 2026-09-29: "add the
            line through D1-10 ... but ensure it isnt LOWESS").
        show_features: Draw the shoal zones, piers and groin on the
            real-domains figure (the annotated figure always has them).
        model_label_plain: Legend the model as just run.model_name
            ("CASCADE"); the wave settings go in the subtitle instead.
        show_wave_climate: Add run.wave_climate as a subtitle line.
        quantity: "rate" (m/yr, the default) or "position" -- the caller
            passes the model's end-minus-start change (m) as `change_rate`
            and a CoastSat series already scaled to metres
            (coastsat_lowess.scale_coastsat_series). Only the publication
            labels and caption read it.
        observed_label: Legend text for the observed LOWESS curve in the
            publication style; None keeps "CoastSat LRR (7-domain LOWESS)".
        title, subtitle: With publication_text, a centred title (and a
            smaller centred line under it) in place of the left-aligned
            period tag. None keeps the tag (Hannah, 2026-10-06: "academic
            and informative").
        model_label: Publication legend text for the model curve; None
            keeps "CASCADE".
        observed_description: Caption phrase for the observed quantity in
            position mode, e.g. "total change: the per-transect LRR fitted on
            1996-2010 x 14 yr".
        publication_text: Replace the title and provenance lines with ONE
            short tag above the plot (period and wave settings), write the
            rest to CAPTIONS.md beside the PNG, and legend in the house
            wording (Hannah, 2026-09-29: "more academic and professional").
            Both the real-domains and the annotated (with-buffers) figure;
            the annotated one also drops its S/N compass text, which the
            x label already says and which sat where the tag goes.
        These four are off by default; rerender_run_figures.py --lowess-only
        turns them on together (the sensitivity style).
        plot_reference_period: Also show the CoastSat period that doesn't
            match the run's start year (faded).
        raw_scatter_size: Marker area, points^2.
        raw_scatter_alpha: Opacity for the active period; the reference
            period (if shown) uses 0.35x this.
        domain_tick_step: X-axis tick spacing, in GIS domains.
        ylim: (low, high) in m/yr to fix the y axis, or None (the default) to
            fit it to the data. Set when a SET of runs has to be read on one
            scale (Hannah, 2026-09-27: every figure in raw_runs/matrix on one
            y axis); rerender_run_figures.py --ylim passes it. Applies to
            the with-buffers figures.
        ylim_real: the same for the real-domains-only figure, which has no
            buffer swings to hold and can be tighter (Hannah, 2026-09-27:
            "so we are efficient with space"). None falls back to ylim.
            rerender_run_figures.py --ylim-real passes it.
    """

    window_colors: dict = dataclasses.field(
        default_factory=lambda: {7: "#6BAED6", 10: "#08519C"}
    )
    window_color_default: str = "#4A7C8E"
    window_styles: tuple = ((1.6, "-", 1.00), (1.8, "-", 1.00))
    raw_color: str = "#5BA3C9"
    plot_domain_means: bool = True
    domain_mean_lw: float = 0.7
    domain_mean_alpha: float = 0.55
    plot_raw_lrr: bool = True
    raw_lrr_southern_only: bool = True
    overlay_raw_southern_only: bool = False
    south_mean_line: bool = False
    show_features: bool = False
    model_label_plain: bool = False
    show_wave_climate: bool = False
    publication_text: bool = False
    quantity: str = "rate"
    observed_label: str = None
    observed_description: str = None
    title: str = None
    subtitle: str = None
    model_label: str = None
    plot_reference_period: bool = False
    raw_scatter_size: float = 6
    raw_scatter_alpha: float = 0.60
    domain_tick_step: int = 5
    ylim: tuple = None
    ylim_real: tuple = None


DEFAULT_RATE_COMPARISON = RateComparisonConfig()


def coastsat_domain_mean(cs):
    """(gis ids, unsmoothed per-domain mean rate) for one CoastSat series.

    The plain mean of the transect LRRs in each 500 m domain -- what
    domain_lrr_summary.csv holds -- recomputed from the same transect arrays
    the LOWESS is fitted to so the two curves on a figure are the same data
    two ways, not two files.
    """
    dom = cs["transect_domains"]
    return compute_domain_means(dom, cs["transect_rates"],
                                int(np.min(dom)), int(np.max(dom)))


def south_domain_mean_line(ax, cs, lowess_config, config, gis_x_transform=None,
                           zorder=4):
    """The raw domain means over D1..skip, dashed, for one CoastSat series.

    Not a LOWESS: south of skip the target IS the unsmoothed domain mean, and
    the dashes say so against the solid LOWESS line north of it. Returns the
    label drawn, or None when there is nothing to draw.
    """
    skip = lowess_config.skip_southern_domains
    if skip <= 0:
        return None
    if gis_x_transform is None:
        gis_x_transform = lambda gis_x: gis_x
    mean_x, mean_y = coastsat_domain_mean(cs)
    keep = np.asarray(mean_x) <= skip
    if not keep.any():
        return None
    widest = max(lowess_config.window_domains)
    lbl = f"{cs['label']} — domain mean, D1–{skip} (not smoothed)"
    ax.plot(gis_x_transform(np.asarray(mean_x)[keep]), np.asarray(mean_y)[keep],
            color=config.window_colors.get(widest, config.window_color_default),
            lw=config.window_styles[0][0], ls=SOUTH_MEAN_DASH, zorder=zorder,
            label=lbl)
    return lbl


SOUTH_MEAN_DASH = (0, (3, 1.5))


def _model_label(run, config):
    """Legend text for the model curve."""
    if config.model_label_plain or run.Hs is None:
        return f"{run.model_name}"
    return f"{run.model_name}, Hs {run.Hs} m"


def _features_only(annotations):
    """The annotation layer minus the town spans and village lines, for a
    figure whose towns are already drawn by town_bands. The groin label goes
    beside the line at the bottom, south of it: D1-5 run positive in both
    windows, so the lower-left corner is the one empty place near the groin."""
    return dataclasses.replace(annotations, town_spans={}, village_lines={},
                               groin_label_y=0.03, groin_label_side="left")


def plot_coastsat_overlay(ax, cs_series, lowess_config, config, x_transform, gis_x_transform=None):
    """Draw CoastSat raw-transect scatter + LOWESS lines for one axis.

    Public because the sensitivity figures draw many model curves on one axis
    and need the SAME observed layer underneath them; a second implementation
    there would let the two drift apart on styling and on the southern splice.

    Shared by both branches of plot_rate_comparison (REAL vs ALL domains);
    the annotated figure has extra features (fill_between, legend handle
    tracking) and builds its overlay separately.

    Args:
        ax: Axes to draw on.
        cs_series: Output of build_coastsat_series.
        lowess_config: LowessConfig (only window_domains is read here, to
            find the widest window).
        config: RateComparisonConfig.
        x_transform: along_coast_m array -> x-axis coordinate array, for
            the raw scatter.
        gis_x_transform: gis-domain-ID array -> x-axis coordinate array,
            for the LOWESS lines. Defaults to identity (GIS-ID x-axis).
    """
    if gis_x_transform is None:
        gis_x_transform = lambda gis_x: gis_x

    widest_win = max(lowess_config.window_domains)
    for cs in cs_series:
        is_active = cs["active"]
        if not is_active and not config.plot_reference_period:
            continue
        if config.plot_raw_lrr:
            x = x_transform(cs["transect_along_coast"])
            y = cs["transect_rates"]
            if config.overlay_raw_southern_only and config.raw_lrr_southern_only:
                south = np.asarray(cs["transect_domains"]) <= lowess_config.skip_southern_domains
                x, y = np.asarray(x)[south], np.asarray(y)[south]
            raw_alpha = config.raw_scatter_alpha if is_active else config.raw_scatter_alpha * 0.35
            ax.scatter(x, y, color=config.raw_color,
                       s=config.raw_scatter_size, alpha=raw_alpha, zorder=1, linewidths=0)
        if config.plot_domain_means:
            mean_x, mean_y = coastsat_domain_mean(cs)
            if len(mean_x):
                ax.plot(gis_x_transform(mean_x), mean_y,
                        color=config.window_colors.get(widest_win,
                                                       config.window_color_default),
                        lw=config.domain_mean_lw, zorder=2,
                        alpha=config.domain_mean_alpha * (1.0 if is_active else 0.40))
        if config.south_mean_line and is_active:
            south_domain_mean_line(ax, cs, lowess_config, config,
                                   gis_x_transform=gis_x_transform)
        for idx, win in enumerate(cs["windows"]):
            cs_color = config.window_colors.get(win["window"], config.window_color_default)
            lw_base, ls, alpha_factor = (
                config.window_styles[idx] if idx < len(config.window_styles) else (1.5, "-", 0.80)
            )
            is_widest = (win["window"] == widest_win)
            plot_gis_x, plot_y = splice_lowess_with_raw_south(
                win["gis_x"], win["smoothed"], cs["transect_domains"], cs["transect_rates"],
                skip_n=lowess_config.skip_southern_domains, is_widest_window=is_widest,
            )
            plot_x = gis_x_transform(plot_gis_x)
            lbl = f"{cs['label']} — {win['window']}-domain LOWESS"
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


def _run_parameters(run, scope, domains, endpoints=True, wave_climate=None):
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
    head = "  ·  ".join(bits)
    # The sensitivity style names the wave settings here; its model legend is just "CASCADE"
    if wave_climate:
        head = f"waves: {wave_climate}\n{head}"
    return head + "\n" + tail


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


def _wave_tag(wave_climate):
    """"Hs 2.0 m, Tp 7.5 s, asym 0.6, high-angle 0.5" in figure wording."""
    if not wave_climate:
        return None
    text = " · ".join(part.strip() for part in str(wave_climate).split(","))
    return (text.replace("asym ", "asymmetry ")
                .replace("high-angle ", "high-angle fraction "))


def _publication_axes(ax, run, config=DEFAULT_RATE_COMPARISON):
    """Axis labels and the one-line tag (period and wave settings) that
    replace the title and provenance lines."""
    ax.set_xlabel("Alongshore position (GIS domain, south → north)")
    ax.set_ylabel("Shoreline position change (m)"
                  if config.quantity == "position"
                  else "Shoreline change rate, LRR (m/yr)")
    if config.title:
        # Title above, the one-line context under it, both centred on the plot
        ax.set_title(config.title, loc="center", fontsize=10, color=INK,
                     pad=20 if config.subtitle else 8)
        if config.subtitle:
            ax.text(0.5, 1.02, config.subtitle, transform=ax.transAxes,
                    ha="center", va="bottom", fontsize=8, color=INK_MUTED)
        return
    tag = f"{run.start_year}–{run.end_year}"
    wave = _wave_tag(run.wave_climate)
    if wave:
        tag += f"  ·  {wave}"
    ax.set_title(tag, loc="left", fontsize=9, color=INK, pad=6)


def _publication_legend(fig, config, annotations, lowess_config, extra=(),
                        ncol=None, extra_model_raw=False):
    """The legend in house wording (STYLE.md, 2026-09-19): short noun
    phrases, no period or dataset repeats -- those are in the caption."""
    skip = lowess_config.skip_southern_domains
    widest = max(lowess_config.window_domains)
    blue = config.window_colors.get(widest, config.window_color_default)
    lw = config.window_styles[0][0]
    handles = [
        Line2D([0], [0], color=annotations.model_color, lw=1.8,
               label=config.model_label or "CASCADE"),
        Line2D([0], [0], color=blue, lw=lw,
               label=config.observed_label
               or f"CoastSat LRR ({widest}-domain LOWESS)"),
    ]
    if config.plot_domain_means:
        handles.append(Line2D([0], [0], color=blue, lw=config.domain_mean_lw,
                              alpha=config.domain_mean_alpha,
                              label="Domain mean (unsmoothed)"))
    if extra_model_raw:
        handles.insert(1, Line2D([0], [0], color=annotations.model_color,
                                 lw=config.domain_mean_lw, alpha=config.domain_mean_alpha,
                                 label="CASCADE domain value (unsmoothed)"))
    if config.south_mean_line and skip > 0:
        handles.append(Line2D([0], [0], color=blue, lw=lw, ls=SOUTH_MEAN_DASH,
                              label=f"Domain mean, D1–{skip} (unsmoothed)"))
    if config.plot_raw_lrr:
        handles.append(Line2D([0], [0], color=config.raw_color, marker="o",
                              ms=3, ls="none", alpha=config.raw_scatter_alpha,
                              label="Individual transects"))
    handles += list(extra)
    fig.legend(handles=handles, loc="outside lower center",
               ncol=ncol or len(handles), frameon=False, handlelength=2.2,
               columnspacing=1.8)


def _publication_caption(run, config, annotations, lowess_config, domains,
                         annotated=False):
    """What the title and the three provenance lines used to say, as a
    caption for CAPTIONS.md (house rule: nothing on the canvas that belongs
    in a caption)."""
    skip = lowess_config.skip_southern_domains
    widest = max(lowess_config.window_domains)
    if config.quantity == "position":
        parts = [
            f"Shoreline position change along {annotations.region_name}, "
            f"{run.start_year}–{run.end_year}, GIS domain "
            f"{domains.first_gis_id} (Cape Point) to {domains.last_gis_id} "
            f"(Pea Island); positive is seaward (accretion). Modelled "
            f"(CASCADE, orange): shoreline position at the end of the run "
            f"minus the start. Observed ({annotations.obs_source_name}, "
            f"blue): {config.observed_description}.",
            f"The observed line is smoothed with a {widest}-domain LOWESS "
            f"north of domain {skip}; south of it LOWESS is not applied and "
            f"the line is the unsmoothed domain mean (dashed), with the "
            f"individual transects shown as dots.",
        ]
        return _caption_tail(parts, run, config, annotated, skip, widest)
    parts = [
        f"Modelled (CASCADE, orange) and observed ({annotations.obs_source_name}, "
        f"blue) shoreline change rate along {annotations.region_name}, "
        f"{run.start_year}–{run.end_year}, GIS domain {domains.first_gis_id} "
        f"(Cape Point) to {domains.last_gis_id} (Pea Island); positive is "
        f"seaward (accretion).",
        f"Observed rate is the linear regression rate (LRR) per transect, "
        f"smoothed with a {widest}-domain LOWESS north of domain {skip}; south "
        f"of it LOWESS is not applied and the line is the unsmoothed domain "
        f"mean (dashed), with the individual transects shown as dots.",
    ]
    return _caption_tail(parts, run, config, annotated, skip, widest)


def _caption_tail(parts, run, config, annotated, skip, widest):
    """The legend key, wave climate and run identity every caption ends with."""
    what = "change" if config.quantity == "position" else "rate"
    if annotated:
        parts.append(f"Blue shading: the {widest}-domain LOWESS {what} against "
                     f"zero. Grey dashed line at domain {skip}.5: where the "
                     f"LOWESS begins. Shaded bands: communities; grey dashed "
                     f"lines: village centres. Hatched: shoal zones. "
                     f"Dash-dot: piers. Dotted: Buxton groin.")
    elif config.show_features:
        parts.append("Grey bands: communities. Hatched: shoal zones. "
                     "Dash-dot: piers. Dotted: Buxton groin.")
    wave = _wave_tag(run.wave_climate)
    if wave:
        parts.append(f"Wave climate: {wave.replace(' · ', ', ')}.")
    parts.append(f"Background erosion "
                 f"{'on' if run.background_erosion_on else 'off'}. "
                 f"Run {run.run_name}.")
    return " ".join(parts)


def plot_rate_comparison(change_rate, cs_series, run, real_domains_only=True,
                          estimator=None,
                          domains=DEFAULT_DOMAINS,
                          annotations=DEFAULT_ANNOTATIONS,
                          lowess_config=DEFAULT_LOWESS,
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
    model_label = _model_label(run, config)
    wave_line = run.wave_climate if config.show_wave_climate else None

    if real_domains_only:
        gis_ids = np.arange(domains.first_gis_id, domains.last_gis_id + 1)
        fig, ax = plt.subplots(figsize=figsize("double", aspect=0.48),
                               constrained_layout=True)

        real_rate = change_rate[domains.start_real_index:domains.end_real_index]
        ax.plot(gis_ids, real_rate, color=annotations.model_color,
                linewidth=1.8, label=model_label, zorder=6)

        plot_coastsat_overlay(
            ax, cs_series, lowess_config, config,
            x_transform=lambda along_m: along_m / domains.domain_spacing_m + domains.first_gis_id,
        )

        ax.axhline(0.0, linestyle=(0, (4, 3)), linewidth=0.8, color=INK_MUTED,
                   zorder=2)

        ax.set_xlim(domains.first_gis_id - 0.5, domains.last_gis_id + 0.5)
        step = _tick_step(domains.num_real_domains, config.domain_tick_step)
        ax.set_xticks(np.arange(domains.first_gis_id,
                                domains.last_gis_id + 1, step))
        # After the limits: spans outside the view are skipped, labels clamped to the visible part
        town_bands(ax, spans=annotations.town_spans)
        if config.show_features:
            add_geographic_annotations(ax, _features_only(annotations))

        if config.publication_text:
            _publication_axes(ax, run, config)
        else:
            ax.set_xlabel(DOMAIN_AXIS_LABEL)
            ax.set_ylabel(_rate_axis_label(estimator))
            _provenance(ax, subject,
                        _run_parameters(run, "real domains only", domains,
                                        wave_climate=wave_line))
        ax.grid(axis="y")
        open_frame(ax)
        ax.set_axisbelow(True)
        if config.publication_text:
            _publication_legend(fig, config, annotations, lowess_config)
        else:
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
            ax, cs_series, lowess_config, config,
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
        # Extra pad: the title row sits above the secondary GIS axis
        _provenance(ax, subject,
                    _run_parameters(run, "all domains, buffers included",
                                    domains),
                    extra_pad=26)
        ax.grid(axis="y")
        open_frame(ax)
        ax.set_axisbelow(True)
        fig.legend(loc="outside lower center", ncol=3, frameon=False)

        fig_suffix = "ALL_DOMAINS_WITH_BUFFERS"

    fixed = (config.ylim_real if real_domains_only and config.ylim_real is not None
             else config.ylim)
    if fixed is not None:
        ax.set_ylim(*fixed)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"  Saved plot: {save_path}")
        if config.publication_text and real_domains_only:
            record_caption(Path(save_path), _publication_caption(
                run, config, annotations, lowess_config, domains))
    if show:
        plt.show()

    return fig, ax, fig_suffix


def plot_annotated_rate_comparison(change_rate, cs_series, run,
                                    estimator=None,
                                    domains=DEFAULT_DOMAINS,
                                    annotations=DEFAULT_ANNOTATIONS,
                                    lowess_config=DEFAULT_LOWESS,
                                    config=DEFAULT_RATE_COMPARISON,
                                    sea_level_rise_rate_m_yr=None,
                                    save_path=None, show=False,
                                    change_rate_raw=None):
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
        change_rate_raw: Optional, same shape as change_rate: the model
            before smoothing, drawn thin and faint under it (when
            change_rate is the smoothed model, as the observation is).

    Returns:
        (fig, ax)
    """
    gis_ids = np.arange(domains.first_gis_id, domains.last_gis_id + 1)
    real_rate = change_rate[domains.start_real_index:domains.end_real_index]
    real_raw = (None if change_rate_raw is None
                else change_rate_raw[domains.start_real_index:domains.end_real_index])

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.62),
                           constrained_layout=True)

    # Geographic annotations drawn first so data renders on top.
    add_geographic_annotations(
        ax, dataclasses.replace(annotations, groin_label_y=0.03,
                                groin_label_side="left")
        if config.publication_text else annotations)

    data_handles = []
    widest_win = max(lowess_config.window_domains)
    # What one transect value is: a rate, or a position change in metres
    per_transect = "transect net change" if config.quantity == "position" else "transect LRR"
    for cs in cs_series:
        is_active = cs["active"]
        if not is_active and not config.plot_reference_period:
            continue

        scatter_x = cs["transect_along_coast"] / domains.domain_spacing_m + domains.first_gis_id
        if config.plot_raw_lrr:
            if config.raw_lrr_southern_only:
                south_mask = cs["transect_domains"] <= lowess_config.skip_southern_domains
                scatter_x_plot = scatter_x[south_mask]
                scatter_y_plot = cs["transect_rates"][south_mask]
                raw_lbl = (f"{cs['label']} — {per_transect} "
                           f"(D{domains.first_gis_id}-{lowess_config.skip_southern_domains})"
                           if is_active else None)
            else:
                scatter_x_plot = scatter_x
                scatter_y_plot = cs["transect_rates"]
                raw_lbl = f"{cs['label']} — {per_transect}" if is_active else None
            raw_alpha = config.raw_scatter_alpha if is_active else config.raw_scatter_alpha * 0.35
            ax.scatter(scatter_x_plot, scatter_y_plot, color=config.raw_color,
                       s=config.raw_scatter_size, alpha=raw_alpha, zorder=1,
                       linewidths=0, label=raw_lbl)
            if is_active:
                data_handles.append(
                    Line2D([0], [0], color=config.raw_color, marker=".", ms=5,
                           ls="none", alpha=config.raw_scatter_alpha, label=raw_lbl)
                )

        if config.plot_domain_means:
            mean_x, mean_y = coastsat_domain_mean(cs)
            if len(mean_x):
                mean_c = config.window_colors.get(widest_win,
                                                  config.window_color_default)
                mean_a = config.domain_mean_alpha * (1.0 if is_active else 0.40)
                mean_lbl = f"{cs['label']} — domain mean (unsmoothed)"
                ax.plot(mean_x, mean_y, color=mean_c, lw=config.domain_mean_lw,
                        alpha=mean_a, zorder=3.5,
                        label=mean_lbl if is_active else None)
                if is_active:
                    data_handles.append(
                        Line2D([0], [0], color=mean_c, lw=config.domain_mean_lw,
                               alpha=mean_a, label=mean_lbl)
                    )
        for idx, win in enumerate(cs["windows"]):
            cs_color = config.window_colors.get(win["window"], config.window_color_default)
            lw_base, ls, alpha_factor = (
                config.window_styles[idx] if idx < len(config.window_styles) else (1.5, "-", 0.80)
            )
            is_widest = (win["window"] == widest_win)
            gis_x, rate = splice_lowess_with_raw_south(
                win["gis_x"], win["smoothed"], cs["transect_domains"], cs["transect_rates"],
                skip_n=lowess_config.skip_southern_domains, is_widest_window=is_widest,
            )
            lbl = f"{cs['label']} — {win['window']}-domain LOWESS"
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

    if config.south_mean_line:
        for cs in cs_series:
            if not cs["active"]:
                continue
            lbl = south_domain_mean_line(ax, cs, lowess_config, config, zorder=5)
            if lbl:
                data_handles.append(Line2D(
                    [0], [0], color=config.window_colors.get(
                        widest_win, config.window_color_default),
                    lw=config.window_styles[0][0], ls=SOUTH_MEAN_DASH,
                    label=lbl))

    model_label = _model_label(run, config)
    if real_raw is not None:
        ax.plot(gis_ids, real_raw, color=annotations.model_color,
                lw=config.domain_mean_lw, alpha=config.domain_mean_alpha, zorder=5.5)
    ax.plot(gis_ids, real_rate, color=annotations.model_color, linewidth=2.0,
            zorder=6, label=model_label)
    data_handles.append(
        Line2D([0], [0], color=annotations.model_color, lw=2.0,
               label=model_label)
    )

    ax.axhline(0, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 3)),
               zorder=3)

    # Scatter/LOWESS transition marker, only when the southern domains show raw scatter
    if config.plot_raw_lrr and config.raw_lrr_southern_only and lowess_config.skip_southern_domains > 0:
        ax.axvline(lowess_config.skip_southern_domains + 0.5, color=INK_MUTED,
                   lw=0.8, ls=(0, (4, 4)), zorder=2)

    ax.set_xlim(domains.first_gis_id - 0.5, domains.last_gis_id + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    # 1 m/yr steps suit a rate; a position change spans +/-50 m.
    ax.yaxis.set_major_locator(
        ticker.MaxNLocator(nbins=8, steps=[1, 2, 2.5, 5, 10])
        if config.quantity == "position" else ticker.MultipleLocator(1))
    ax.tick_params(axis="both", which="minor", length=1.8)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)

    # Lock ylim to data range, then place accretion/erosion side labels.
    all_vals = np.concatenate(
        [real_rate] + ([real_raw] if real_raw is not None else []) + [w["smoothed"][np.isfinite(w["smoothed"])]
                       for cs in cs_series for w in cs["windows"]]
        # The unsmoothed means swing wider than the LOWESS, so they set the bound too
        + ([coastsat_domain_mean(cs)[1] for cs in cs_series]
           if config.plot_domain_means else [])
        # So do the southern dots, which run past every curve at Cape Point
        + ([cs["transect_rates"][cs["transect_domains"]
                                 <= lowess_config.skip_southern_domains]
            for cs in cs_series if cs["active"]]
           if config.plot_raw_lrr and config.raw_lrr_southern_only else [])
    )
    ymin, ymax = all_vals.min(), all_vals.max()
    ypad = (ymax - ymin) * 0.06
    ax.set_ylim(*(config.ylim if config.ylim is not None
                  else (ymin - ypad, ymax + ypad)))

    ybot, ytop = ax.get_ylim()
    zero_frac = (0 - ybot) / (ytop - ybot)
    acc_y = (annotations.label_accretion_y if annotations.label_accretion_y is not None
             else zero_frac + (1 - zero_frac) / 2)
    ero_y = (annotations.label_erosion_y if annotations.label_erosion_y is not None
             else zero_frac / 2)
    # A white backing: the erosion label can sit on the observed curve
    for _y, _txt in ((acc_y, "accretion \u25b2"), (ero_y, "erosion \u25bc")):
        ax.text(1.0, _y, _txt, transform=ax.transAxes, fontsize=8,
                color=INK_MUTED, ha="right", va="center", zorder=8,
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none",
                          boxstyle="square,pad=0.2"))

    if config.publication_text:
        _publication_axes(ax, run, config)
        _publication_legend(fig, config, annotations, lowess_config,
                            extra=annotation_legend_handles(annotations), ncol=5,
                            extra_model_raw=real_raw is not None)
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight",
                        facecolor="white")
            record_caption(Path(save_path), _publication_caption(
                run, config, annotations, lowess_config, domains,
                annotated=True))
            print(f"  Saved annotated plot: {save_path}")
        if show:
            plt.show()
        return fig, ax

    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Shoreline Position Change (m)" if config.quantity == "position"
                  else _rate_axis_label(estimator, title_case=True))
    ax.text(0.0, 1.005, f"\u2190 {annotations.low_end_label}",
            transform=ax.transAxes, fontsize=7.5, color=INK_MUTED, ha="left",
            va="bottom", clip_on=False)
    ax.text(1.0, 1.005, f"{annotations.high_end_label} \u2192",
            transform=ax.transAxes, fontsize=7.5, color=INK_MUTED, ha="right",
            va="bottom", clip_on=False)

    # The run's identity, kept on the canvas (see the module docstring)
    _provenance(
        ax,
        f"Modelled against {annotations.obs_source_name} shoreline change, "
        f"{run.start_year}–{run.end_year}",
        _run_parameters(run,
                        f"{annotations.obs_source_name} "
                        f"{'net change' if config.quantity == 'position' else 'LRR'} per "
                        f"{int(domains.domain_spacing_m)} m domain",
                        domains, endpoints=False,
                        wave_climate=(run.wave_climate
                                      if config.show_wave_climate else None)),
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
