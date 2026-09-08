"""
hat_figure_style.py
==============================================================================
One typographic and colour standard for the Hatteras figures.

WHY THIS EXISTS
    Each plotting script had been choosing its own font sizes, its own greys and
    reds, its own elevation ramp and its own way of labelling panels. Put two of
    them side by side in a document and they read as coming from different
    papers. Worse, the same quantity was drawn in different colours in different
    figures, so "grey" meant "v1" in one and "not applied" in another.

    This module fixes the vocabulary: colours carry meaning, and the meaning is
    the same everywhere.

COLOUR SEMANTICS -- do not reassign these locally
    BASE      the unmodified input (v1 as extracted)
    ACCENT    the modification under test (v2, the insert)
    ROAD      NC-12
    ADDED     ground that was fabricated
    WATER     cells at or below sea level
    REF       a reference value: a median, a target, an observation

ELEVATION IS DRAWN IN CLASSES, NOT A RAMP
    A continuous ramp is the wrong tool here. The back-barrier sits a few
    decimetres below MHW and the dune is five metres above it, so a linear ramp
    renders the entire island as one flat tone and hides the only distinction
    that matters -- which cells are land. `elevation_cmap()` returns a discrete
    scale with a hard break at 0 m.

USAGE
    from hat_figure_style import apply_style, C, panel_label, elevation_cmap
    apply_style()
==============================================================================
"""

from __future__ import annotations

import matplotlib as mpl
from matplotlib.colors import BoundaryNorm, ListedColormap

# Colour-blind safe: the base/accent pair is grey against a dark red, which
# separates on luminance as well as hue, so it survives greyscale printing.
C = {
    "BASE": "#7f7f7f",
    "BASE_FILL": "#d9d9d9",
    "ACCENT": "#9e2a2b",
    "ACCENT_FILL": "#e8c4c4",
    "ROAD": "#1a1a1a",
    "ADDED": "#c8880f",
    "ADDED_FILL": "#f6e3b8",
    "WATER": "#a8c8e0",
    "REF": "#2c6e49",
    "GRID": "#cccccc",
}

# Elevation classes, m MHW. The first edge is the water break.
ELEV_BOUNDS = [-99.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 99.0]
_ELEV_COLOURS = [
    C["WATER"],
    "#fbf7e3", "#f0e6bf", "#e0cf96",
    "#cdb26f", "#b3904f", "#8f6b36", "#5f4520",
]


def apply_style() -> None:
    """Set rcParams. Call once, before creating any figure."""
    mpl.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 300,
        # NOT "tight": it re-crops after layout, which drags
        # absolutely-positioned colourbars and captions on top of
        # the panels. Figures here set their own margins.
        "savefig.bbox": "standard",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 9,
        "axes.titlesize": 9.5,
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.titlepad": 6,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#333333",
        # Top and right spines carry no information on any of these panels.
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": C["GRID"],
        "grid.linewidth": 0.6,
        "grid.alpha": 0.7,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "legend.handlelength": 1.8,
        "lines.linewidth": 1.4,
        "lines.solid_capstyle": "butt",
    })


def elevation_cmap():
    """(cmap, norm, bounds) for the discrete elevation scale."""
    cmap = ListedColormap(_ELEV_COLOURS)
    return cmap, BoundaryNorm(ELEV_BOUNDS, cmap.N), ELEV_BOUNDS


def panel_title(letter: str, text: str) -> str:
    """"(a) text" -- the panel letter carried IN the title.

    Preferred over a floating label in axes coordinates, which has to be nudged
    per panel depending on whether that panel has a y-axis label, and collides
    with the title the moment a figure is resized.
    """
    return "({}) {}".format(letter, text)


def caption(fig, text: str, y: float = 0.005, size: float = 7.8) -> None:
    """A figure caption, left-aligned under everything.

    Figures in this project get read months later out of a folder, detached from
    whatever conversation produced them, so each one states its own method.
    """
    fig.text(0.055, y, text, fontsize=size, va="bottom", ha="left",
             color="#333333", wrap=True)


def spines_for_image(ax) -> None:
    """Restore all four spines. An imshow panel needs its frame closed."""
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.8)
        s.set_edgecolor("#333333")
