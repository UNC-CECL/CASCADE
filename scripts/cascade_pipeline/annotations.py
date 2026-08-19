"""Geographic reference annotations shared by every cascade_pipeline figure.

Community spans, village centers, piers, groins, and shoal-influence zones
are drawn identically on the shoreline GIF and both rate-comparison
figures, so the styling lives in one place instead of three.

Ships with NO site content: AnnotationConfig()'s dict fields default to
empty and its text fields default to generic placeholders. Build a
populated instance for your own site (see e.g. hatteras_site_config.py)
and pass it explicitly -- a reusable library shouldn't silently draw
somebody else's place names.
"""

import dataclasses

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory


@dataclasses.dataclass(frozen=True)
class AnnotationConfig:
    """Geographic annotation layer + site labels, keyed by GIS domain ID.

    Attributes:
        town_spans: {label: (gis_lo, gis_hi)} shaded community zones.
        village_lines: {label: gis_id} dashed village-center lines.
        piers: {label: (gis_id, label_y)} pier lines; label_y is an axes
            fraction (0=bottom, 1=top) for the rotated label.
        groins: {label: gis_id} groin lines; gis_id may be fractional
            (e.g. 5.5 = the boundary between domains 5 and 6).
        shoal_zones: {label: (gis_lo, gis_hi)} shoal/shoreline-position
            influence zones -- any number of them, not a fixed pair.
        pier_label_y: Default label_y for piers (informational; the actual
            y used is whatever is stored per-pier in `piers`).
        groin_label_y: Label y (axes fraction) for groin lines.
        label_accretion_y: Fixed axes-fraction y for the "Accretion" side
            label on the annotated figure, or None to auto-place at the
            midpoint between the zero line and the top.
        label_erosion_y: Same as label_accretion_y, for "Erosion" (midpoint
            between the bottom and the zero line).
        region_name: Site name used in figure titles, e.g. "Hatteras Island".
        low_end_label: Text for the low-domain-ID end of the compass
            annotation, e.g. "S" or "S | South Landmark".
        high_end_label: Text for the high-domain-ID end, e.g. "N" or
            "North Landmark | N".
        obs_source_name: Name of the observational dataset shown against
            the model, e.g. "CoastSat".
        color_*: Colors for each annotation layer.
        model_color: Line color for the modeled shoreline-change-rate curve.
    """

    town_spans: dict = dataclasses.field(default_factory=dict)
    village_lines: dict = dataclasses.field(default_factory=dict)
    piers: dict = dataclasses.field(default_factory=dict)
    groins: dict = dataclasses.field(default_factory=dict)
    shoal_zones: dict = dataclasses.field(default_factory=dict)

    pier_label_y: float = 0.76
    groin_label_y: float = 0.68
    label_accretion_y: float = None
    label_erosion_y: float = None

    region_name: str = "Study Area"
    low_end_label: str = "S"
    high_end_label: str = "N"
    obs_source_name: str = "Observed"

    color_town_span: str = "#90AFC5"
    color_shoal: str = "#E0A800"
    color_village_line: str = "0.40"
    color_pier: str = "#1565C0"
    color_groin: str = "#B71C1C"
    model_color: str = "#FF8C00"


DEFAULT_ANNOTATIONS = AnnotationConfig()


def add_geographic_annotations(ax, config=DEFAULT_ANNOTATIONS):
    """Add every geographic reference annotation to an axis (bottom -> top).

    Draw order: shoal zones, community spans, village center lines, pier
    lines, groin lines. Label y-positions use blended axes-fraction coords,
    so they hold position regardless of the y-data range. X-axis must be in
    GIS domain IDs. Any category left empty in `config` simply draws
    nothing -- there's no minimum set of annotations required.

    Args:
        ax: Matplotlib Axes to annotate.
        config: AnnotationConfig; defaults to an empty (no-op) layout.
    """
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    for zone_label, (z_lo, z_hi) in config.shoal_zones.items():
        ax.axvspan(z_lo - 0.5, z_hi + 0.5,
                   color=config.color_shoal, alpha=0.10, zorder=0,
                   hatch="///", edgecolor=config.color_shoal, linewidth=0)
        ax.text((z_lo + z_hi) / 2.0, 0.04,
                f"{zone_label}\nPosition", transform=trans,
                ha="center", va="bottom", fontsize=7, color="#7A5800", style="italic",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    for span_label, (d_lo, d_hi) in config.town_spans.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=config.color_town_span, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    for vname, dom in config.village_lines.items():
        ax.axvline(dom, color=config.color_village_line, lw=0.9, ls="--", alpha=0.65, zorder=1)
        ax.text(dom, 0.84, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    for pname, (dom, lbl_y) in config.piers.items():
        ax.axvline(dom, color=config.color_pier, lw=1.0, ls="-.", alpha=0.80, zorder=2)
        ax.text(dom, lbl_y, pname, transform=trans,
                ha="center", va="top", fontsize=7, color=config.color_pier, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    for gname, dom in config.groins.items():
        ax.axvline(dom, color=config.color_groin, lw=1.1, ls=":", alpha=0.85, zorder=2)
        ax.text(dom, config.groin_label_y, gname, transform=trans,
                ha="center", va="top", fontsize=7, color=config.color_groin, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))


def annotation_legend_handles(config=DEFAULT_ANNOTATIONS):
    """Proxy legend artists matching add_geographic_annotations' layers.

    Only returns a handle for categories that are actually populated in
    `config`, so an empty site config produces an empty legend rather than
    five swatches for things that weren't drawn.

    Args:
        config: AnnotationConfig; should match what was passed to
            add_geographic_annotations, or the legend will describe
            something different from what's drawn.

    Returns:
        List of matplotlib legend handles (Patch/Line2D).
    """
    handles = []
    if config.town_spans:
        handles.append(Patch(fc=config.color_town_span, alpha=0.30, label="Community"))
    if config.shoal_zones:
        handles.append(Patch(fc=config.color_shoal, alpha=0.25, hatch="///",
                              edgecolor=config.color_shoal, linewidth=0,
                              label="Shoal influence zone"))
    if config.village_lines:
        handles.append(Line2D([0], [0], color=config.color_village_line, lw=0.9, ls="--",
                               label="Village center"))
    if config.piers:
        handles.append(Line2D([0], [0], color=config.color_pier, lw=1.0, ls="-.", label="Pier"))
    if config.groins:
        handles.append(Line2D([0], [0], color=config.color_groin, lw=1.1, ls=":", label="Groin"))
    return handles
