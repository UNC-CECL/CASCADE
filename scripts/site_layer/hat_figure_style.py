"""
One typographic and colour standard for every Hatteras figure.

    from site_layer.hat_figure_style import apply_style, figsize, save
    python scripts/site_layer/hat_figure_style.py   # writes STYLE.md and the style sheet

Typeface, palette, column widths, panel letters, captions to CAPTIONS.md and
figure folders; the rules are written out in scripts/figure_making/STYLE.md. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""

from __future__ import annotations

import datetime as _dt
import functools
import re
from pathlib import Path

import matplotlib as mpl
import numpy as np
import matplotlib.patheffects as pe
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.transforms import offset_copy

# Run as a file, scripts/ is not on sys.path: add it so the site_layer imports resolve
if __package__ in (None, ""):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
from site_layer.hat_map_layers import STYLE_DOC, STYLE_SHEET_DIR  # noqa: E402,F401

# Where figures go: a subject outside FIGURE_SUBJECTS raises rather than starting a new folder
OUTPUT_ROOT = PROJECT_ROOT / "output"
FIGURES_ROOT = OUTPUT_ROOT / "figures"
COMPARISONS_ROOT = OUTPUT_ROOT / "comparisons"       # cross-run figures
OBSERVATIONS_OUT = OUTPUT_ROOT / "observations"      # the observed record itself
# The numbered layout, in the paper's order; the retired subjects raise
FIGURE_SUBJECTS = {
    "site": "1-site",                  # the reach, the 90 domains, one domain
    "observations": "2-observations",  # CoastSat shoreline, dune lines, imagery
    "inputs": "3-model-inputs",        # how each model input is built, by step
    "mechanics": "4-model-mechanics",  # how Barrier3D, BRIE and CASCADE work
    "results": "5-results",            # hindcasts and scenarios
    "style": "style",
    "talk": "talk",                    # projector versions, mirroring the above
}
# the model-input steps under 3-model-inputs/, numbered as data/hatteras_init/
INPUT_STEPS = ("0-elevation", "1-domains", "2-brie-offset", "3-forcing",
               "4-management", "5-observed-target", "7-source-sink")


def figure_dir(subject: str, *parts: str) -> Path:
    """output/figures/<numbered subject>[/<parts>...]; subject must be a known one."""
    if subject not in FIGURE_SUBJECTS:
        raise ValueError(f"unknown figure subject {subject!r}; one of "
                         f"{', '.join(FIGURE_SUBJECTS)} (output/figures/README.md)")
    if subject == "inputs" and parts and parts[0] not in INPUT_STEPS:
        raise ValueError(f"unknown model-input step {parts[0]!r}; one of {INPUT_STEPS}")
    return FIGURES_ROOT.joinpath(FIGURE_SUBJECTS[subject], *parts)

# Type, ink, the vintage pair

FONT_STACK = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
INK = "0.15"            # text, axes, baselines
INK_MUTED = "0.42"      # secondary labels, rulers, guide lines
GRID_C = "0.88"         # hairline grid

# The ColorBrewer RdBu poles: red the earlier vintage, blue the later, a fill the band between
C_1984 = "#b2182b"      # RdBu, dark red
C_1997 = "#2166ac"      # RdBu, dark blue
C_1984_FILL = "#f4a582"  # RdBu, light red  - the band where 1984 lies seaward
C_1997_FILL = "#92c5de"  # RdBu, light blue - the band where 1984 lies landward

# Colour-blind safe: grey base against a purple accent, distinct from the vintage pair
C = {
    "BASE": "#7f7f7f",
    "BASE_FILL": "#d9d9d9",
    "ACCENT": "#7b3294",
    "ACCENT_FILL": "#c2a5cf",
    "ROAD": "#1a1a1a",
    "ADDED": "#c8880f",
    "ADDED_FILL": "#f6e3b8",
    "WATER": "#a8c8e0",
    "REF": "#2c6e49",
    "GRID": GRID_C,
    "INK": INK,
    "INK_MUTED": INK_MUTED,
    "EARLY": C_1984,
    "LATE": C_1997,
    "EARLY_FILL": C_1984_FILL,
    "LATE_FILL": C_1997_FILL,
    "LOCATOR": "#1b7f8c",     # the study area on a locator inset: teal, never the groin red
    "GROIN": "#d7191c",       # the Buxton groin field marker on maps
}

# An ordered variable (the smoothing width) gets a light-to-dark blue ramp, not the vintage pair
SMOOTH_RAMP = ("#9ecae1", "#6baed6", C_1997, "#08306b")

CELL_M = 10.0           # the Barrier3D cell; scale bars under 1 km say it

# Column widths: a figure is drawn at its printed width, so its type prints at 8-9 pt
FIG_W_SINGLE = 3.54     # in, 90 mm
FIG_W_DOUBLE = 7.48     # in, 190 mm
FIG_H_MAX = 9.4         # in, a page less its caption


def figsize(width="double", aspect=0.5, height=None):
    """(w, h) in inches: `width` is "single", "double" or a number of inches;
    the height is `height` or width * aspect, capped at a page."""
    w = {"single": FIG_W_SINGLE, "double": FIG_W_DOUBLE}.get(width, width)
    w = float(w)
    h = float(height) if height else w * float(aspect)
    return (w, min(h, FIG_H_MAX))


# One label for the alongshore axis; the endpoints belong in the caption
DOMAIN_AXIS_LABEL = "GIS domain (south → north)"


def town_bands(ax, where="top", label=True, shade="0.94", spans=None,
               fontsize=7, strip=None):
    """Village spans as light bands behind an alongshore axis, named once.
    `spans` is {name: (first_gis, last_gis)}; default the site config's.
    Call it AFTER the axis limits are set: spans outside the view are skipped
    and a label is clamped to the visible part of its span.

    `strip` draws the spans as a band of that fraction of the axes height
    against the `where` edge, instead of a full-height wash. Use it when the
    panel already shades something else: two full-height greys on one panel
    cannot be told apart (the road alongshore figure, 2026-09-10)."""
    if spans is None:
        try:
            from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS
            spans = HATTERAS_ANNOTATIONS.town_spans
        except ImportError:
            return
    # Only the spans this panel shows, labels clamped to the visible part
    x_lo, x_hi = sorted(ax.get_xlim())
    for name, (lo, hi) in spans.items():
        if hi + 0.5 < x_lo or lo - 0.5 > x_hi:
            continue
        if strip:
            y0 = 1.0 - float(strip) if where == "top" else 0.0
            ax.add_patch(mpl.patches.Rectangle(
                (lo - 0.5, y0), (hi + 0.5) - (lo - 0.5), float(strip),
                transform=ax.get_xaxis_transform(), facecolor=shade,
                edgecolor="none", zorder=0, clip_on=True))
        else:
            ax.axvspan(lo - 0.5, hi + 0.5, color=shade, lw=0, zorder=0)
        if label:
            mid = (max(lo - 0.5, x_lo) + min(hi + 0.5, x_hi)) / 2
            if strip:
                y = 1.0 - float(strip) / 2 if where == "top" else float(strip) / 2
                va = "center"
            else:
                y = 0.985 if where == "top" else 0.015
                va = "top" if where == "top" else "bottom"
            ax.text(mid, y, name, transform=ax.get_xaxis_transform(),
                    ha="center", va=va, fontsize=fontsize, color=INK_MUTED,
                    zorder=1, clip_on=True)


STRUCTURE_LABEL_PT = 6.5

# Where a structure label may sit, tried in order
_LABEL_SLOTS = (("bottom", "right"), ("bottom", "left"),
                ("top", "right"), ("top", "left"))


def _points_under(ax, txt) -> int:
    """How many plotted points fall inside `txt`'s box (a little padded), over
    every Line2D and LineCollection drawn in data coordinates."""
    renderer = ax.figure.canvas.get_renderer()
    bb = txt.get_window_extent(renderer).expanded(1.3, 1.08)
    (x0, y0), (x1, y1) = ax.transData.inverted().transform(
        [[bb.x0, bb.y0], [bb.x1, bb.y1]])
    def inside(x, y):
        # Test the line's segments, not its vertices: 20 points along each
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        if x.size < 2:
            return 0
        t = np.linspace(0.0, 1.0, 20)
        xd = (x[:-1, None] + (x[1:] - x[:-1])[:, None] * t).ravel()
        yd = (y[:-1, None] + (y[1:] - y[:-1])[:, None] * t).ravel()
        return int(((xd >= x0) & (xd <= x1) & (yd >= y0) & (yd <= y1)).sum())

    n = 0
    for line in ax.get_lines():
        if line.get_transform() is not ax.transData or not line.get_visible():
            continue
        n += inside(line.get_xdata(), line.get_ydata())
    for coll in ax.collections:
        if not hasattr(coll, "get_segments"):
            continue
        for seg in coll.get_segments():
            seg = np.asarray(seg, dtype=float)
            n += inside(seg[:, 0], seg[:, 1])
    return n


def _place_label(ax, pos, name, label_pt):
    """Put `name` along the line at `pos` in the first slot that covers no
    plotted point, else the slot that covers fewest (Hannah, 2026-09-15:
    labels were sitting on the data at Rodanthe Pier)."""
    best, best_n = None, None
    for where, ha in _LABEL_SLOTS:
        y, va = (0.03, "bottom") if where == "bottom" else (0.86, "top")
        # A 2 pt gap between the text and its line, so the halo does not hide the line
        tr = offset_copy(ax.get_xaxis_transform(), fig=ax.figure,
                         x=(-2.0 if ha == "right" else 2.0), units="points")
        txt = ax.text(pos, y, name, transform=tr,
                      rotation=90, ha=ha, va=va, fontsize=label_pt,
                      color=INK_MUTED, zorder=7, path_effects=_halo(2.0))
        n = _points_under(ax, txt)
        if n == 0:
            if best is not None:
                best.remove()
            return txt
        if best_n is None or n < best_n:
            if best is not None:
                best.remove()
            best, best_n = txt, n
        else:
            txt.remove()
    return best


def structures(ax, label=True, label_pt=STRUCTURE_LABEL_PT, spans=None):
    """The Buxton groin (solid hairline) and the two piers (dotted hairlines)
    on an alongshore axis, each named once along its own line, reading upward
    (Hannah, 2026-09-15). The lines stop short of the top so the village
    labels there stay clear of them. A label goes at the bottom of its line
    unless data is drawn there, in which case it moves to the other side of
    the line or to the top -- so call this AFTER the data is plotted AND
    after anything that resizes the axes at draw time (an outside legend, a
    colourbar): the test is made in the layout as it stands.
    `spans` is a site AnnotationConfig; default the Hatteras one. Lived in
    coastsat_lrr_windows.py until 2026-09-15, when a third alongshore figure
    wanted it."""
    if spans is None:
        try:
            from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS
            spans = HATTERAS_ANNOTATIONS
        except ImportError:
            return
    if label:
        ax.figure.canvas.draw()   # settle the layout so text extents are real

    for name, pos in spans.groins.items():
        ax.axvline(pos, ymax=0.88, color=INK, lw=0.7, zorder=6)
        if label:
            _place_label(ax, pos, name, label_pt)
    for name, (pos, _frac) in spans.piers.items():
        ax.axvline(pos, ymax=0.88, color=INK_MUTED, lw=0.6,
                   ls=(0, (1.5, 1.5)), zorder=6)
        if label:
            _place_label(ax, pos, name, label_pt)


# A figure folder shows figures: everything else a script writes goes under supporting/
SUPPORT_DIR = "supporting"


def support_dir(folder) -> Path:
    """`<folder>/supporting/`, created. Where a figure script puts everything
    that is not a PNG."""
    d = Path(folder) / SUPPORT_DIR
    d.mkdir(parents=True, exist_ok=True)
    return d


def save(fig, path, vector=True, close=False, **kwargs):
    """PNG at 300 dpi in the folder and, for `vector`, a PDF with the same
    stem under `supporting/`, so a line or bar figure stays sharp in a
    manuscript without the folder showing two files per figure. Returns the
    paths."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = [path.with_suffix(".png")]
    fig.savefig(out[0], **kwargs)
    if vector:
        out.append(support_dir(path.parent) / path.with_suffix(".pdf").name)
        fig.savefig(out[-1], **kwargs)
    if close:
        import matplotlib.pyplot as plt
        plt.close(fig)
    return out

STYLE_RC = {
    "font.family": "sans-serif", "font.sans-serif": FONT_STACK,
    "font.size": 9,
    "axes.titlesize": 10, "axes.titleweight": "normal", "axes.titlepad": 6,
    "axes.titlelocation": "left",
    "axes.labelsize": 9, "axes.labelcolor": INK,
    "axes.edgecolor": INK, "axes.linewidth": 0.6,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "grid.color": GRID_C, "grid.linewidth": 0.5, "grid.linestyle": "-",
    "legend.fontsize": 8, "legend.frameon": True, "legend.framealpha": 0.95,
    "legend.edgecolor": "none", "legend.fancybox": False,
    "legend.handlelength": 1.8, "legend.borderpad": 0.5,
    "legend.labelspacing": 0.4, "legend.columnspacing": 1.4,
    "lines.solid_capstyle": "butt",
    "text.color": INK,
    "figure.dpi": 130,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.dpi": 300,
}


def apply_style() -> None:
    """Idempotent. Called at the top of every figure so the style holds no
    matter which entry point drew it, including an import from elsewhere."""
    mpl.rcParams.update(STYLE_RC)


# Elevation classes

# Elevation, m MHW: the terrain ramp over land, one colour for water (Hannah, 2026-10-01)
# The tick edges of the old classes; callers use them for colorbar ticks and legends
ELEV_BOUNDS = [-99.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 99.0]
ELEV_MAX_M = 4.0                 # the top of the land ramp; higher ground takes its colour
ELEV_WATER = "#333399"           # plt.cm.terrain(0.0), terrain's own navy
_TERRAIN_SEA_POS = 0.35          # where terrain's land colours begin on its 0-1 scale
_ELEV_STEPS = 256


def elevation_cmap():
    """(cmap, norm, bounds) for domain topography: land on the terrain ramp from 0
    to ELEV_MAX_M m MHW, anything higher the top colour, and water below 0 m
    terrain's navy, as the colormap's under colour. A MASKED cell stays
    transparent, as before: a panel that masks water sets its background to
    ELEV_WATER. `bounds` are colorbar tick positions, not bins.

    Until 2026-10-01 this returned seven brown classes; Hannah asked for terrain
    everywhere domain topography is drawn, as on the inputs overview. A plain
    Normalize, not a fine BoundaryNorm: the colorbar drew a tick at every step."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    land = plt.cm.terrain(_np_linspace(_TERRAIN_SEA_POS, 1.0, _ELEV_STEPS))
    cmap = ListedColormap(land, name="hat_elevation")
    cmap.set_under(ELEV_WATER)
    cmap.set_over(land[-1])
    return cmap, Normalize(vmin=0.0, vmax=ELEV_MAX_M), ELEV_BOUNDS


# Error and cost surfaces

# Error surfaces are drawn without hue, truncated at both ends, so the marks on them stand out
_ERROR_LO, _ERROR_HI = 0.08, 0.86


def error_cmap(reverse: bool = False):
    """The greyscale ramp for a scalar error or cost surface; dark is worse.

    `reverse=True` for a surface where HIGH is better (a score, a share
    explained), so that dark still means the outcome you do not want."""
    import matplotlib.pyplot as plt
    base = plt.get_cmap("Greys")
    lo, hi = (_ERROR_HI, _ERROR_LO) if reverse else (_ERROR_LO, _ERROR_HI)
    return ListedColormap(base(_np_linspace(lo, hi, 256)),
                          name="hat_error" + ("_r" if reverse else ""))


def _np_linspace(a, b, n):
    """numpy.linspace without importing numpy at module scope: this module is
    imported by scripts that have not yet chosen a backend, and it has stayed
    free of the numeric stack."""
    step = (b - a) / (n - 1)
    return [a + step * i for i in range(n)]


# Panels, maps, frames

def _letter(i: int) -> str:
    return f"({chr(ord('a') + i)})"


def _title(ax, i: int, text: str) -> None:
    """Panel letter at the left, title centred - both above the axes."""
    ax.set_title(_letter(i), loc="left", fontweight="bold")
    ax.set_title(text, loc="center")


def _letter_inside(ax, i: int) -> None:
    """The letter inside the top-left corner, for a panel whose title is wide
    enough to run under a letter placed beside it."""
    ax.text(0.03, 0.985, _letter(i), transform=ax.transAxes, ha="left",
            va="top", fontsize=10, fontweight="bold", color=INK, zorder=20,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                      boxstyle="square,pad=0.15"))


def panel_title(letter: str, text: str) -> str:
    """"(a) text" as one string, for a script that sets the title itself.
    `_title()` is preferred: it puts the letter in bold at the left and the
    title centred, which is what the house figures do."""
    return "({}) {}".format(letter, text)


def _halo(lw: float = 2.5):
    return [pe.withStroke(linewidth=lw, foreground="white")]


def _north_arrow(ax, x=0.90, y=0.10, length=0.045) -> None:
    """A north arrow in axes fraction. Only on maps WITHOUT coordinate ticks -
    a labelled UTM frame is already north-up by construction."""
    ax.annotate("", xy=(x, y + length), xytext=(x, y),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0,
                                shrinkA=0, shrinkB=0,
                                mutation_scale=11), zorder=20)
    ax.text(x, y + length + 0.008, "N", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=8.5, fontweight="bold", color=INK,
            zorder=20,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none",
                      boxstyle="square,pad=0.1"))


def _scalebar(ax, length_m: float = 5 * CELL_M, cell_m: float = CELL_M,
              show_cells: bool | None = None) -> None:
    """A bar, because the coordinate ticks are gone. Drawn in data units, so
    it scales with the panel and cannot disagree with it. White halo rather
    than a box, so it sits on relief without blanking it. Under 1 km the label
    also says how many Barrier3D cells that is (pass show_cells=False to
    suppress it on a map that is not a model grid)."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    bx = x0 + 0.06 * (x1 - x0)
    by = y0 + 0.028 * (y1 - y0)
    tick = 0.004 * (y1 - y0)
    ax.plot([bx, bx + length_m], [by, by], color=INK, lw=2.2,
            solid_capstyle="butt", zorder=12, path_effects=_halo(4.2))
    for e in (bx, bx + length_m):
        ax.plot([e, e], [by - tick, by + tick], color=INK, lw=1.2,
                zorder=12, path_effects=_halo(3.0))
    if show_cells is None:
        show_cells = length_m < 1000
    if length_m >= 1000:
        label = f"{length_m / 1000:g} km"
    elif show_cells:
        label = f"{length_m:.0f} m  ({length_m / cell_m:.0f} cells)"
    else:
        label = f"{length_m:.0f} m"
    ax.text(bx + length_m / 2, by + 1.6 * tick, label,
            ha="center", va="bottom", fontsize=8, color=INK, zorder=12,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none",
                      boxstyle="square,pad=0.15"))


def spines_for_image(ax) -> None:
    """All four spines, closed: an imshow or map panel needs its frame."""
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(STYLE_RC["axes.linewidth"])
        s.set_edgecolor(INK)


def open_frame(ax) -> None:
    """Top and right spines off: a chart, not an image."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


# public names without the underscore, for scripts that prefer them
title = _title
letter = _letter
letter_inside = _letter_inside
north_arrow = _north_arrow
scalebar = _scalebar
halo = _halo


# Captions: off the canvas, into CAPTIONS.md beside the figure

# A figure's method goes to CAPTIONS.md beside the PNG, never onto the canvas

_CAPTION_ATTR = "_hat_caption"


def caption(fig, text: str, y: float = 0.005, size: float = 7.8) -> None:
    """Register `text` as the figure's caption. Written to CAPTIONS.md beside
    the image by the figure's next savefig; nothing is drawn. `y` and `size`
    are accepted for old callers and ignored."""
    setattr(fig, _CAPTION_ATTR, text)
    if not getattr(fig, "_hat_savefig_wrapped", False):
        original = fig.savefig

        @functools.wraps(original)
        def _savefig(fname, *args, **kwargs):
            out = original(fname, *args, **kwargs)
            cap = getattr(fig, _CAPTION_ATTR, None)
            # Only the PNG gets an entry (the PDF beside it shares the stem)
            if cap and isinstance(fname, (str, Path)) and Path(fname).suffix.lower() == ".png":
                record_caption(Path(fname), cap)
            return out

        fig.savefig = _savefig
        fig._hat_savefig_wrapped = True


def _prune_captions(body: str, fig_dir: Path) -> str:
    """Drop entries whose figure is no longer beside the CAPTIONS.md.

    `record_caption` replaces an entry by name but never removed one, so a
    renamed or deleted figure left its caption behind forever (68 such
    entries had accumulated by 2026-09-21, after the filename rename).
    Pruning on every write keeps the file honest without anyone having to
    remember.

    Only an entry whose named file is MISSING is dropped, so a run that
    redraws one figure never touches the captions of the others.
    """
    pattern = re.compile(r"^\*\*`([\w.\-]+)`\.\*\*.*?(?=\n\*\*`|\Z)", re.S | re.M)
    kept = [m.group(0).strip() for m in pattern.finditer(body)
            if (fig_dir / m.group(1)).exists()]
    head = body[:m.start()] if (m := pattern.search(body)) else body
    return head.rstrip("\n") + "\n\n" + "\n\n".join(kept) + "\n"


def mark_offaxis(ax, x, y, half, color=None, size=13.0):
    """Mark values beyond +/-`half` at the axis edge, and say which they were.

    A LINE that leaves the axis is simply cut by `ylim`, with nothing on the
    canvas to say so, so a domain at +136 m reads as +100 (found 2026-09-22,
    when the metre figures were fixed to +/-100 m). Scatter points already get
    this treatment; this is the same for a line, and the caption clause below
    names the values so nothing is lost silently.

    Args:
        ax: the axes, already at its final ylim.
        x, y: the series as drawn. y beyond +/-half is what gets marked.
        half: the axis half-range.
        color: marker colour; the axes' ink by default.
        size: marker area in points squared.

    Returns:
        [(x, y), ...] for the off-axis points, largest |y| first, for
        `offaxis_clause()`.
    """
    import numpy as _np
    x = _np.asarray(x, dtype=float)
    y = _np.asarray(y, dtype=float)
    out = _np.isfinite(y) & (_np.abs(y) > half)
    if not out.any():
        return []
    ax.scatter(x[out], _np.where(y[out] > 0, half, -half), s=size,
               marker="^", zorder=14, clip_on=False,
               c=[color or INK],
               transform=ax.transData)
    pairs = sorted(zip(x[out], y[out]), key=lambda t: -abs(t[1]))
    return [(float(a), float(b)) for a, b in pairs]


def offaxis_clause(named, half, unit="m"):
    """" Beyond +/-100 m, off the axis: the observed change +136 m at GIS 1."

    Args:
        named: [(label, [(x, y), ...]), ...] as returned by `mark_offaxis`.
        half: the axis half-range, for the sentence.
        unit: the y unit.
    """
    hits = [f"{label} {v:+.0f} {unit} at GIS {int(g)}"
            for label, pts in named for g, v in pts]
    if not hits:
        return ""
    return (f" Beyond \u00b1{half:g} {unit}, off the axis and marked with a "
            "triangle at the edge: " + "; ".join(hits) + ".")


def compare_header(fig, lines, size=8.5):
    """The 'what is being compared' line(s) above the panels.

    Hannah, 2026-09-22: on a figure that puts two DIFFERENT measurements side
    by side, which is which -- and over what dates -- has to be on the canvas,
    not only in the caption. This is the one thing the style lets above the
    panel titles, and it is a NAMING line, not a result: what each side is and
    what interval it spans, never the numbers that came out. The summary stays
    in `supporting/CAPTIONS.md`.

    `target_comparison` carried this idea first (the source/sink line); it
    is here so the dune-line comparisons place it identically.

    Args:
        fig: the figure.
        lines: one string, or a sequence joined with newlines.
        size: point size; 8.5 sits just under the 10 pt panel titles.
    """
    text = lines if isinstance(lines, str) else chr(10).join(lines)
    fig.suptitle(text, fontsize=size, color=INK)


def record_caption(png_path: Path, text: str) -> Path:
    """Write or replace the entry for `png_path.name` in the CAPTIONS.md under
    `supporting/` beside it. Entries are '**`<file>`.** text' paragraphs;
    other content is kept."""
    png_path = Path(png_path)
    md = support_dir(png_path.parent) / "CAPTIONS.md"
    text = " ".join(str(text).split())
    entry = f"**`{png_path.name}`.** {text}\n"
    if md.is_file():
        body = md.read_text(encoding="utf-8")
        pattern = re.compile(r"^\*\*`" + re.escape(png_path.name) + r"`\.\*\*.*?(?=\n\*\*`|\Z)",
                             re.S | re.M)
        if pattern.search(body):
            body = pattern.sub(lambda _m: entry.rstrip("\n"), body)
        else:
            body = body.rstrip("\n") + "\n\n" + entry
        # One blank line back before every entry, which a replaced entry's match consumes
        body = re.sub(r"(?<!\n)\n(\*\*`)", r"\n\n\1", body)
        body = _prune_captions(body, png_path.parent)
    else:
        body = (f"# Captions — {png_path.parent.name}\n\n"
                f"Written by the figure scripts through `hat_figure_style.caption()`; "
                f"the images carry no titles or footnotes, this file does.\n\n" + entry)
    md.write_text(body, encoding="utf-8")
    return md


# Maps: the house map elements (from the study-area map, 2026-10-02)

# One label rule on maps: white type, a thin dark halo, 8 pt
MAP_HALO = [pe.withStroke(linewidth=1.2, foreground="0.12")]
MAP_TEXT = dict(color="white", fontsize=8, zorder=8, path_effects=MAP_HALO)
# The same rule for a light map (an outline on a pale canvas): ink type, a thin white halo
MAP_TEXT_DARK = dict(color=INK, fontsize=8, zorder=8,
                     path_effects=[pe.withStroke(linewidth=1.2, foreground="white")])
# The groin field: a red bar with a white edge
GROIN_MARKER = dict(marker="|", ms=11, mew=2.4, color=C["GROIN"], ls="none",
                    path_effects=[pe.withStroke(linewidth=4.4, foreground="white")])


def spaced_caps(text: str) -> str:
    """Water-body names in upright letter-spaced capitals: 'P A M L I C O   S O U N D'."""
    return "   ".join(" ".join(w.upper()) for w in text.split())


def water_label(ax, x, y, text, text_kw=None, **kw):
    """A water body (sound, ocean, inlet): upright letter-spaced capitals, centred."""
    return ax.text(x, y, spaced_caps(text), **{"ha": "center", "va": "center",
                                               **(text_kw or MAP_TEXT), **kw})


def place_label(ax, x, y, text, text_kw=None, **kw):
    """A village or place name (Buxton, Cape Point, Pea Island): italic."""
    return ax.text(x, y, text, **{"ha": "center", "va": "center", "fontstyle": "italic",
                                  **(text_kw or MAP_TEXT), **kw})


def map_label(ax, x, y, text, text_kw=None, **kw):
    """Anything else on a map (domain numbers, scale figures): plain."""
    return ax.text(x, y, text, **{"ha": "center", "va": "center", **(text_kw or MAP_TEXT), **kw})


def scale_bar_km(ax, length_m: float = 10_000, segments: int = 2, x: float = 0.035,
                 y: float = 0.085, text_kw=None, unit: str = "km") -> None:
    """A cartographic scale bar: alternating black and white segments with a white keyline,
    labelled 0, 5, 10 km below with a clear gap. `x`, `y` place its left end in axes
    fraction; drawn in data units, so the map must be equal-aspect in metres. `unit="m"`
    labels the ticks in metres (for a bar under 1 km)."""
    from matplotlib.patches import Rectangle
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    h = 0.013 * (y1 - y0)
    bx, by = x0 + x * (x1 - x0), y0 + y * (y1 - y0)
    seg = length_m / segments
    # a white keyline so the black segment reads on dark water
    ax.add_patch(Rectangle((bx, by), length_m, h, facecolor="none", edgecolor="white", lw=2.2,
                           zorder=11))
    for i in range(segments):
        ax.add_patch(Rectangle((bx + i * seg, by), seg, h, facecolor=INK if i % 2 == 0 else "white",
                               edgecolor=INK, lw=0.6, zorder=12))
    div = 1000.0 if unit == "km" else 1.0
    for i in range(segments + 1):
        v = i * seg / div
        lab = f"{v:g} {unit}" if i == segments else f"{v:g}"
        ax.text(bx + i * seg, by - 1.5 * h, lab, ha="center", va="top",
                **{**(text_kw or MAP_TEXT), "zorder": 12})


def north_dart(ax, c, north=(0.0, 1.0), arrow_m: float = 2000.0, text_kw=None) -> None:
    """A split-dart north arrow, half black and half white with a white keyline, its "N"
    beyond the tip. `c` is its centre in data units; `north` is the unit vector of true
    north in the map's own frame (a quarter-turned map passes its rotation's north)."""
    from matplotlib.patches import Polygon
    n = np.asarray(north, float)
    perp = np.array([-n[1], n[0]])
    c = np.asarray(c, float)
    tip, back = c + n * arrow_m / 2, c - n * arrow_m / 2
    notch = c - n * arrow_m * 0.22
    w = 0.28 * arrow_m
    ax.add_patch(Polygon([tip, back + perp * w, notch, back - perp * w], closed=True, facecolor="none",
                         edgecolor="white", lw=2.2, zorder=11))
    ax.add_patch(Polygon([tip, back + perp * w, notch], closed=True, facecolor=INK, edgecolor=INK,
                         lw=0.6, zorder=12))
    ax.add_patch(Polygon([tip, notch, back - perp * w], closed=True, facecolor="white", edgecolor=INK,
                         lw=0.6, zorder=12))
    lab = tip + n * 0.40 * arrow_m
    ax.text(lab[0], lab[1], "N", ha="center", va="center", fontweight="bold",
            **{**(text_kw or MAP_TEXT), "fontsize": 8.5, "zorder": 12})


def groin_marker(ax, x, y, **kw):
    """The groin field at (x, y): the red, white-edged bar."""
    return ax.plot(x, y, **{"zorder": 9, **GROIN_MARKER, **kw})


def groin_handle(label: str = "Buxton groins"):
    """The legend entry for `groin_marker`, a little smaller than on the map."""
    from matplotlib.lines import Line2D
    return Line2D([], [], label=label, **{**GROIN_MARKER, "ms": 9, "mew": 2.0})


# The house map legend: translucent white, no frame, 8 pt, tight
MAP_LEGEND = dict(frameon=True, framealpha=0.72, edgecolor="none", facecolor="white", fontsize=8,
                  handlelength=1.4, labelspacing=0.3, borderpad=0.4, borderaxespad=0.0,
                  handletextpad=0.5)


def letter_at(ax, i: int, x: float, y: float, ha: str = "left") -> None:
    """The bold panel letter in its white box, top-left corner at (x, y) in axes fraction."""
    ax.text(x, y, f"({chr(ord('a') + i)})", transform=ax.transAxes, ha=ha, va="top",
            fontsize=10, fontweight="bold", color=INK, zorder=20,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.15"))


def top_row(ax, ax_in, handles, letter: int | None = 0, letter_x: float = 0.012, gap: float = 0.008,
            north=None, arrow_gap_m: float = 3200.0, arrow_m: float = 2000.0, legend_kw=None,
            text_kw=None):
    """The top row of a map with an inset in its upper right: the legend top-aligned
    immediately left of the inset, the panel letter at the left edge on the same top line,
    and (given `north`) the north arrow just left of the legend, centred on it. `ax_in`
    is pinned to its upper-right corner first, so the measured top is the drawn top.
    With `handles=None` there is no legend and the arrow sits just left of the inset,
    level with its upper part. Returns the legend (or None)."""
    fig = ax.figure
    ax_in.set_anchor("NE")
    ax_in.apply_aspect()
    to_ax = ax.transAxes.inverted()
    in_left = to_ax.transform(fig.transFigure.transform(ax_in.get_position().p0))[0]
    in_top = to_ax.transform(fig.transFigure.transform(ax_in.get_position().p1))[1]
    leg = None
    if handles is not None:
        leg = ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(in_left - gap, in_top),
                        ncol=1, **{**MAP_LEGEND, **(legend_kw or {})})
    if letter is not None:
        # the letter box's own padding, so its edge and the inset's share one line
        pad = 1.5 / 72 / (ax.get_position().height * fig.get_size_inches()[1])
        letter_at(ax, letter, letter_x, in_top - pad)
    if north is not None and leg is not None:
        fig.canvas.draw()
        lb = leg.get_window_extent().transformed(ax.transData.inverted())
        north_dart(ax, (lb.x0 - arrow_gap_m, (lb.y0 + lb.y1) / 2), north=north, arrow_m=arrow_m,
                   text_kw=text_kw)
    elif north is not None:
        x_d, y_d = ax.transData.inverted().transform(ax.transAxes.transform((in_left, in_top)))
        y_lo = ax.transData.inverted().transform(ax.transAxes.transform((in_left, in_top - 0.25)))[1]
        north_dart(ax, (x_d - arrow_gap_m, (y_d + y_lo) / 2), north=north, arrow_m=arrow_m,
                   text_kw=text_kw)
    return leg


def axes_point(ax, x: float, y: float):
    """An axes-fraction point in data units, for placing `north_dart` by eye."""
    return ax.transData.inverted().transform(ax.transAxes.transform((x, y)))


# The style sheet: output/figures/style/, and STYLE.md beside the figure code

def write_style_sheet(out_dir: Path | None = None) -> tuple[Path, Path]:
    """A swatch figure and a STYLE.md stating the rules, so the standard can be
    seen and read without opening this file. Re-run after changing anything
    above; both files say when they were written."""
    import matplotlib.pyplot as plt
    import numpy as np

    # A given out_dir takes both files (a trial); by default the sheet and STYLE.md go apart
    md = (Path(out_dir) / "STYLE.md") if out_dir else STYLE_DOC
    out_dir = Path(out_dir) if out_dir else STYLE_SHEET_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    apply_style()

    fig = plt.figure(figsize=figsize("double", aspect=0.86), constrained_layout=True)
    # Three rows: the error ramp needs its own strip
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 0.30, 1.15])

    # (a) the palette
    ax = fig.add_subplot(gs[0, 0])
    swatches = [("C_1984  earlier vintage", C_1984), ("C_1984_FILL", C_1984_FILL),
                ("C_1997  later vintage", C_1997), ("C_1997_FILL", C_1997_FILL),
                ("C['BASE']  unmodified input", C["BASE"]), ("C['ACCENT']  modification", C["ACCENT"]),
                ("C['ADDED']  fabricated ground", C["ADDED"]), ("C['WATER']", C["WATER"]),
                ("C['REF']  reference value", C["REF"]), ("C['ROAD']  NC-12", C["ROAD"]),
                ("C['GROIN']  groin field on maps", C["GROIN"]),
                ("C['LOCATOR']  study area on an inset", C["LOCATOR"]),
                ("INK", INK), ("INK_MUTED", INK_MUTED)]
    for k, (name, col) in enumerate(swatches):
        y = len(swatches) - 1 - k
        ax.add_patch(plt.Rectangle((0, y + 0.15), 1.2, 0.7, facecolor=col, edgecolor="none"))
        ax.text(1.45, y + 0.5, f"{name}   {col}", va="center", fontsize=8, color=INK)
    ax.set_xlim(0, 8)
    ax.set_ylim(0, len(swatches))
    ax.set_xticks([])
    ax.set_yticks([])
    _title(ax, 0, "colours, and what each one means")

    # (b) the elevation classes
    ax = fig.add_subplot(gs[0, 1])
    cmap, norm, bounds = elevation_cmap()
    demo = np.linspace(-0.5, 4.5, 400)[None, :].repeat(20, axis=0)
    ax.imshow(demo, cmap=cmap, norm=norm, aspect="auto", extent=(-0.5, 4.5, 0, 1))
    ax.set_yticks([])
    ax.set_xlabel("elevation (m above MHW)")
    spines_for_image(ax)
    _title(ax, 1, "elevation: terrain")

    # (c) the error ramp, in its own strip
    axe = fig.add_subplot(gs[1, 1])
    ramp = np.linspace(0, 1, 400)[None, :].repeat(20, axis=0)
    axe.imshow(ramp, cmap=error_cmap(), aspect="auto", extent=(0, 1, 0, 1))
    axe.plot([0.28], [0.5], marker="*", ms=11, color=C["ACCENT"],
             markeredgecolor="white", markeredgewidth=0.7)
    axe.plot([0.62], [0.5], marker="o", ms=7, color="none",
             markeredgecolor=C["REF"], markeredgewidth=1.6)
    axe.set_xticks([])
    axe.set_yticks([])
    axe.set_xlabel("no hue, dark is worse; the marks carry the colour")
    spines_for_image(axe)
    _title(axe, 2, "an error surface")

    # (d) a chart in the style
    ax = fig.add_subplot(gs[2, 0])
    x = np.arange(1, 13)
    rng = np.random.default_rng(4)
    early = 20 + 6 * np.sin(x / 2) + rng.normal(0, 1, x.size)
    late = early - 4 - 1.5 * np.cos(x / 3)
    ax.fill_between(x, early, late, where=early >= late, color=C_1984_FILL, alpha=0.8, lw=0)
    ax.plot(x, early, color=C_1984, lw=1.6, label="1984 (earlier vintage)")
    ax.plot(x, late, color=C_1997, lw=1.6, label="1997 (later vintage)")
    ax.axhline(20, color=C["REF"], lw=1.0, ls=(0, (4, 3)), label="reference value")
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax, spans={"a village": (3, 5), "another": (9, 11)})
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("quantity (m)")
    ax.legend(loc="lower left")
    _title(ax, 3, "a chart: open frame, hairline grid")

    # (e) a map panel in the style
    ax = fig.add_subplot(gs[2, 1])
    yy, xx = np.mgrid[0:60, 0:80]
    relief = 2.5 * np.exp(-((yy - 30) / 12.0) ** 2) - 0.3 + 0.2 * np.sin(xx / 7.0)
    ax.imshow(relief, cmap=cmap, norm=norm, extent=(0, 800, 0, 600), origin="lower")
    ax.plot([0, 800], [330, 300], color=C_1984, lw=2.0)
    ax.plot([0, 800], [350, 335], color=C_1997, lw=2.0)
    ax.plot([0, 800], [150, 140], color=C["ROAD"], lw=2.6)
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)
    ax.set_xlim(0, 800)
    ax.set_ylim(0, 600)
    scale_bar_km(ax, 200.0, unit="m", y=0.10)
    north_dart(ax, (60, 470), arrow_m=60.0)
    groin_marker(ax, 690, 60)
    place_label(ax, 560, 400, "a village")
    water_label(ax, 400, 40, "Ocean")
    _title(ax, 4, "a map: scale bar, arrow, labels")

    swatch = save(fig, out_dir / "HAT_figure_style_sheet.png", close=True)[0]

    stamp = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    md.write_text(f"""# Hatteras figure style

Written {stamp} by `scripts/site_layer/hat_figure_style.py` (`write_style_sheet()`); the
module is the source, this page is its rendering. `HAT_figure_style_sheet.png`
beside it shows every colour, the elevation classes, a chart and a map drawn
under the rules.

## How a script uses it

```python
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import apply_style, C, C_1984, C_1997, INK, INK_MUTED, _title, _scalebar, _north_arrow, caption
apply_style()                 # before any figure is made
```

A map also takes the map elements: `scale_bar_km`, `north_dart`, `water_label`,
`place_label`, `map_label`, `groin_marker`, `groin_handle`, `top_row`, `MAP_TEXT`.
`figure_making/island/study_area_figures.py` (`fig_study_area`) is the worked example.

Every figure script under `scripts/` that draws for this project calls
`apply_style()` first. `0-elevation/3-figures/HAT_plot_duneline_offset.py` also
re-exports these names, so `import HAT_plot_duneline_offset as off` still gives
`off.INK`, `off._title()` and so on to the scripts that take their map loaders
from it.

## The rules

| | |
|---|---|
| size | drawn at the printed width: `figsize("single")` = {FIG_W_SINGLE:.2f} in (90 mm), `figsize("double")` = {FIG_W_DOUBLE:.2f} in (190 mm), height from `aspect`; never wider, so the type below is the type on the page |
| typeface | {", ".join(FONT_STACK)}: the first one installed; 9 pt body, 10 pt panel titles, 8 pt ticks and legends |
| alongshore axis | one label, `DOMAIN_AXIS_LABEL` = "{DOMAIN_AXIS_LABEL}"; villages as light bands named once by `town_bands(ax)`; the endpoints (1 at Cape Point, 90 at Pea Island) go in the caption |
| ink | text and axes `{INK}`, secondary text and rulers `{INK_MUTED}`, grid `{GRID_C}`; axes 0.6 pt |
| panels | a bold letter at the left of the title, the title centred (`_title(ax, i, text)`); inside the corner when the title is wide (`_letter_inside`) |
| maps | closed frame (`spines_for_image`), no coordinate ticks. Since 2026-10-02 the house map is the study-area map: a segmented scale bar, `scale_bar_km(ax)` (black and white segments with a white keyline, 0 / 5 / 10 km below it with a clear gap; `unit="m"` under 1 km), and the split-dart north arrow `north_dart(ax, c, north=...)`, half black and half white, pointing to true north in the map's own frame (a rotated map passes its rotation's north). `_scalebar` (says the cell count under 1 km) and `_north_arrow` remain for the model-grid panels. A labelled UTM frame needs neither |
| map labels | one rule, `MAP_TEXT`: white 8 pt type with a thin dark halo (1.2 pt). Water bodies upright in letter-spaced capitals (`water_label`, `spaced_caps`), villages and places italic (`place_label`), domain numbers and scale figures plain (`map_label`). On a pale map (an outline on a light canvas) the same rule in ink with a thin white halo, `MAP_TEXT_DARK`. A name never sits on land: a short leader takes it out over water |
| map corners | with an inset in the upper right, `top_row(ax, ax_in, handles, north=...)` puts the legend (`MAP_LEGEND`: translucent white, no frame, 8 pt, tight) top-aligned immediately left of the inset, the panel letter at the left edge on the same top line, and the north arrow just left of the legend; the scale bar sits alone in the lower left |
| map symbols | the Buxton groin field is `groin_marker` (red `{C["GROIN"]}` bar with a white edge; `groin_handle()` for the legend, labelled "Buxton groins"); the study area on a locator inset is teal `{C["LOCATOR"]}`, never red, so it cannot be read as the groins; NC-12 on imagery is road yellow `#ffd23f`; a locator inset carries a light graticule without labels, its spacing in the caption |
| charts | top and right spines off (`open_frame`), hairline grid on the value axis only when it helps |
| legends | frameless (a faint white backing when inside), outside the axes where the layout allows: `fig.legend(handles, loc="outside lower center", ncol=n, frameon=False)` under `constrained_layout` |
| vintages | the earlier line or surface is red `{C_1984}`, the later blue `{C_1997}`, everywhere the two are drawn together; the light fills `{C_1984_FILL}` / `{C_1997_FILL}` are the band between them |
| semantic colours | `C["BASE"]` {C["BASE"]} unmodified input · `C["ACCENT"]` {C["ACCENT"]} the modification under test · `C["ROAD"]` {C["ROAD"]} NC-12 · `C["ADDED"]` {C["ADDED"]} fabricated ground · `C["WATER"]` {C["WATER"]} · `C["REF"]` {C["REF"]} a reference value |
| elevation | terrain, everywhere domain topography is drawn (since 2026-10-01; classes before): `elevation_cmap()` gives water below 0 m MHW one colour, terrain's navy `{ELEV_WATER}`, and land the terrain ramp from 0 to {ELEV_MAX_M:g} m, higher ground its top colour; colorbar ticks at {", ".join(f"{b:g}" for b in ELEV_BOUNDS[1:-1])} m |
| error surfaces | greyscale, no hue: `error_cmap()` (dark is worse; `reverse=True` where high is better). A scalar error or cost over a parameter grid is BACKGROUND, and all colour is reserved for what is marked on top of it -- the best cell, the chosen pair, a constraint, an iso-product curve |
| the canvas | no title sentences, statistics lines or footnote paragraphs on the image. That text goes in `supporting/CAPTIONS.md` beside the figure. ONE exception since 2026-09-22: a figure comparing two different measurements may carry a `compare_header()` line saying WHAT each side is and over what dates - a naming line, never a result. `caption(fig, text)` writes it there on the figure's next `savefig`; scripts with their own captions file (dune-line offset, footprint, road relocation) write it themselves |
| the folder | a figure folder shows figures: PNGs at the top level and nothing else. The PDFs, `CAPTIONS.md`, any table or `PROVENANCE.md` a figure script writes go under `{SUPPORT_DIR}/` (`save()` and `record_caption()` do this; a script's own files use `support_dir(folder)`). Since 2026-09-15 |
| legend wording | no working vocabulary: not "today's setback", "v2"/"v3", "as placed", "blank". Say what the thing is: "setback measured on the 1996 surface", "1984 setback (model input)", "rows inserted landward of NC-12", "centreline unchanged between surveys" |
| output | `save(fig, path)`: a 300 dpi PNG and, under `{SUPPORT_DIR}/`, a PDF with the same stem for anything drawn with lines and bars (`vector=False` for image-only panels); white background; `bbox_inches="tight"` only when nothing is positioned absolutely |
| semantic accent | `C["ACCENT"]` is purple since 2026-09-10; it was a red indistinguishable from the 1984 vintage red, so "the change under test" and "1984" read as one colour |

## Where it came from

The 2026-09-04 restyle of the dune-line figures (Hannah: "more academic /
professionally styled so they are informative and look good to present") set
these rules; the older `hat_figure_style.py` of the row-insert work carried the
colour semantics and the elevation classes. The two were merged here on
2026-09-10, with the 09-04 rules winning wherever they disagreed (typeface
order, legend frames, captions on the canvas), so that one style can be applied
across every figure.
""", encoding="utf-8")
    return swatch, md


if __name__ == "__main__":
    s, m = write_style_sheet()
    print(f"wrote {s}\n      {m}")
