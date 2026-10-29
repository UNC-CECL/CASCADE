"""
The house style shared by a run's two comparison figures, the change rate and the position change.

    from cascade_pipeline.plotting.comparison_style import run_subtitle, restyle

Both figures are drawn by plot_annotated_rate_comparison with publication_text on; restyle()
then gives them the same finish: a bold title, the scenario in the subtitle, a km axis on top,
the legend in two headed groups, and the accretion/erosion labels outside the right edge.
Moved out of net_change_comparison.py on 2026-10-26 so the rate figure could share it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-26
"""
from __future__ import annotations

from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation
import numpy as np

from cascade_pipeline import nourishment
from cascade_pipeline.annotations import annotation_legend_handles
from cascade_pipeline.run_names import parse_run_name
from site_layer.hat_observed_rates import WINDOW_ROLE
from site_layer.hatteras_site_config import (
    HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS, HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_PERIODS, first_model_year, last_model_year)

# --- CONFIG ------------------------------------------------------------------
# The top km axis: tick step (km), and the room the title and subtitle make for it (points)
KM_TICK = 5
KM_AXIS_PT = 30
PRESET_TEXT = {"zeroBE": "no source/sink",
               "edgeBE": "source/sink at the end domains only (GIS 1 and 90)",
               "calibBE": "zone source/sink (calibBE)",
               "domainBE": "per-domain source/sink (set 1)"}
MANAGEMENT_TOKENS = (("road", "road maintenance"), ("bdm", "beach and dune management"),
                     ("nourish", "nourishment"), ("reloc", "road relocation"))
# The accretion/erosion labels outside the right edge, turned to read bottom to top; the arrows,
# turned with the text, point up and down
SIGN_LABELS = {"accretion ▲": "accretion →", "erosion ▼": "← erosion"}
# -----------------------------------------------------------------------------


# Whether fills were placed: the old spelling says; the new leaves it to the period's schedule
def _fills_applied(parsed):
    if parsed.fills is not None:
        return parsed.fills
    if not parsed.beach_dune or parsed.fills_withheld:
        return False
    if parsed.start_year not in HATTERAS_PERIODS:
        return False
    return bool(nourishment.build_schedule(
        HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS,
        first_model_year(parsed.start_year), last_model_year(parsed.start_year)).projects)


# Source/sink, groin and management in plain words, from the run name (either spelling, run_names)
def setup_text(run_name):
    parsed = parse_run_name(run_name)
    if parsed is None:
        # Pre-metres and blocking-groin names: read token by token, as before 2026-10-29
        tokens = run_name.split("_")
        preset = tokens[3] if len(tokens) > 3 else ""
        groin = ("blocking groin" if "groinblock" in tokens
                 else "groin" if "groin" in tokens else "no groin")
        management = [text for token, text in MANAGEMENT_TOKENS if token in tokens]
        return PRESET_TEXT.get(preset, preset), groin, management
    on = {"road": parsed.roadway, "bdm": parsed.beach_dune,
          "nourish": _fills_applied(parsed),
          "reloc": parsed.relocations}
    management = [text for token, text in MANAGEMENT_TOKENS if on[token]]
    return (PRESET_TEXT.get(parsed.preset, parsed.preset),
            "groin" if parsed.groin else "no groin", management)


# The scenario (HAT_hindcast_setup.SCENARIOS) from the run name, for the subtitle
def scenario_text(run_name):
    parsed = parse_run_name(run_name)
    if parsed is None:
        tokens = run_name.split("_")
        road, bdm = "road" in tokens, "bdm" in tokens
        withheld, reloc = "nonourish" in tokens, "reloc" in tokens
    else:
        road, bdm = parsed.roadway, parsed.beach_dune
        withheld, reloc = parsed.fills_withheld, parsed.relocations
    scenario = ("full management" if road and bdm else "roadway only" if road
                else "beach/dune only" if bdm else "natural")
    if road and bdm and withheld:
        scenario += " without fills"
    if reloc:
        scenario += " + road relocation"
    groin = setup_text(run_name)[1]
    return scenario + (", no groin" if groin == "no groin" else f" + {groin}")


# The subtitle both figures carry: the window's role, the scenario with the groin, the source/sink
def run_subtitle(window, run_name):
    role = WINDOW_ROLE.get(tuple(int(y) for y in window), "").replace(" (held out)", "")
    return " · ".join(p for p in (role, scenario_text(run_name), setup_text(run_name)[0]) if p)


# The shared one-row legend redrawn as one centred legend whose columns are two headed groups: data, then map
def group_legend(fig, ax):
    old = fig.legends.pop()
    entries = list(zip(old.legend_handles, (t.get_text() for t in old.get_texts())))
    map_labels = {h.get_label() for h in annotation_legend_handles(HATTERAS_ANNOTATIONS)}
    data = [e for e in entries if e[1] not in map_labels]
    features = [e for e in entries if e[1] in map_labels]
    rows = max(1, len(data))

    # Filled column by column: each column is a heading row (blank after a group's first) and `rows` entries
    blank = (Line2D([], [], ls="none"), "")
    columns, headings = [], set()
    for title, group in (("Data", data), ("Map features", features)):
        for i in range(0, len(group), rows):
            body = group[i:i + rows]
            columns.append([(blank[0], title if i == 0 else "")] + body + [blank] * (rows - len(body)))
            if i == 0:
                headings.add(len(columns) - 1)
    cells = [cell for column in columns for cell in column]
    legend = fig.legend(*zip(*cells), loc="outside lower center", ncol=len(columns),
                        frameon=False, handlelength=2.2, columnspacing=2.5)
    for n, text in enumerate(legend.get_texts()):
        if n % (rows + 1) == 0 and n // (rows + 1) in headings:
            text.set_fontweight("bold")
            text.set_x(-legend.handlelength * legend.prop.get_size_in_points())


# Kilometres from the south end of GIS 1 along the top (domain g is centred at (g - 0.5) x spacing);
# the title (bold) and subtitle move up by KM_AXIS_PT to clear it
def km_axis(fig, ax, subtitle, domains=HATTERAS_DOMAINS):
    km = domains.domain_spacing_m / 1000.0
    first = domains.first_gis_id
    top = ax.secondary_xaxis("top", functions=(lambda g: (g - first + 0.5) * km,
                                               lambda d: d / km + first - 0.5))
    top.set_xticks(np.arange(0, (domains.last_gis_id - first + 1) * km + 1e-9, KM_TICK))
    top.set_xlabel("Alongshore distance from the south end of GIS 1 (km)")

    shift = ScaledTranslation(0, KM_AXIS_PT / 72, fig.dpi_scale_trans)
    for text in ax.texts:
        if text.get_text() == subtitle:
            text.set_transform(text.get_transform() + shift)
    title = ax.title
    ax.set_title(title.get_text(), loc="center", fontsize=title.get_fontsize(),
                 color=title.get_color(), fontweight="bold", pad=20 + KM_AXIS_PT)


# The accretion/erosion labels moved from inside the right edge to just outside it, where no curve can run
# under them, each centred on its half of the axis: 0 to the top, the bottom to 0 (axes fractions)
def sign_labels_outside(ax):
    lo, hi = ax.get_ylim()
    zero = -lo / (hi - lo)
    middle = {"accretion ▲": (zero + 1) / 2, "erosion ▼": zero / 2}
    for text in ax.texts:
        if text.get_text() in SIGN_LABELS:
            text.set_y(middle[text.get_text()])
            text.set_text(SIGN_LABELS[text.get_text()])
            text.set_x(1.012)
            # Rotated in anchor mode: ha runs along the text (centred), va across it (top = the side facing the plot)
            text.set_ha("center")
            text.set_va("top")
            text.set_rotation(90)
            text.set_rotation_mode("anchor")
            text.set_bbox(None)
            text.set_clip_on(False)


# The whole finish, on a figure plot_annotated_rate_comparison drew with publication_text on
def restyle(fig, ax, subtitle):
    group_legend(fig, ax)
    km_axis(fig, ax, subtitle)
    sign_labels_outside(ax)
