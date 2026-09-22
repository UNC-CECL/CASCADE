"""
initialization_figures.py
==============================================================================
The island as the model starts it: one page figure per orientation year in
YEARS, showing the t=0 elevation surface every domain is initialised from.

    python scripts/figure_making/island/initialization_figures.py

Writes to output/figures/initialization/<year>/<scheme>/ (ORGANIZATION.md
rule 1), PNG at the top of the scheme folder and PDF + CAPTIONS.md under its
supporting/, through the house `save()`. Three figures per folder:

    island_<year>.png                  DETRENDED, three panels: (a) the BRIE
                                       shoreline offset, (b) the 90 real
                                       domains, (c) the same with the 15
                                       buffer domains at each end
    island_<year>_absolute.png         ABSOLUTE, one panel: every domain at
                                       its true offset, 90 real domains
    island_<year>_absolute_buffers.png ABSOLUTE, one panel: the same with the
                                       buffer domains

...for each of the four hindcast starts (1984, 1996, 2004, 2010) in each of
the two elevation treatments (classes, terrain): 24 figures from two
topography loads. The views are the same data placed different ways, and each
caption points at its companions.

ONE FOLDER PER START, THEN PER TREATMENT
    Two dozen figures in one folder read as two dozen unrelated images; the
    pairing that matters is the views of ONE start, in ONE treatment. So
    <year>/<scheme>/, each with its own CAPTIONS.md, and the file names are
    identical across schemes -- only the path says which treatment it is
    (Hannah, 2026-09-17). The names keep their year, because a figure pulled
    out of its folder must still say which start it is.

WHY THE ABSOLUTE VIEW IS TWO FILES, NOT TWO PANELS
    The real-domain map and the with-buffers map were panels (a) and (b) of
    one figure until 2026-09-17. Each is a map of the whole reach in its own
    right and each is used on its own, so as a pair they cost the page two
    half-height panels to say one thing twice. Split, each gets the full
    column width and its own cross-shore window -- the real-domain map no
    longer reserves the extra kilometre the buffer offsets need. The price is
    that the two no longer share metres per inch, so each caption states its
    span and its exaggeration.

WHICH SURFACE EACH YEAR GETS
    Not one. 1984 and 1996 read the 1984-start extraction (2009+2014 with the
    1996 ALACE graft); 2004 and 2010 read 2004-start (2009+2014 alone).
    YEAR_PRODUCT in hat_topo_version pairs them and product_for_year() is the
    accessor. Until 2026-09-17 this script resolved topo_dirs() ONCE at import,
    so it drew the 2004-start island under a 1984 label -- the exact failure
    product_for_year()'s docstring was written to prevent.

WHY ONE FIGURE DRAWS THE ISLAND LEVEL, NOT AS A DIAGONAL
    Each domain sits at its own BRIE shoreline offset, and those offsets span
    about 6 km across the reach. Drawn in absolute cross-shore space the 2 km
    island became a thin diagonal ribbon crossing an 8 km canvas that was
    three-quarters empty water, and no cross-shore detail survived at page
    width. So the map panels are DETRENDED: every domain is drawn against its
    own frame origin, which is what its elevation array actually holds, and
    the offset that would have displaced it is drawn as panel (a) above, on
    the same alongshore axis. Nothing is hidden -- the offset becomes a number
    you can read off an axis instead of a slope you have to estimate
    (Hannah, 2026-09-17).

    The absolute view is kept alongside it, restyled, because a single canvas
    IS the initial condition and the detrended panels are a rearrangement of
    it (Hannah, 2026-09-17). It carries about 9 km of cross-shore, so it takes
    its own exaggeration -- see VERTICAL_EXAGGERATION_ABSOLUTE.

    Panels (b) and (c) share one alongshore scale in metres per inch, so the
    real span in (b) sits directly above its own position in (c) and the
    buffer domains are visibly the part that sticks out.

THE CROSS-SHORE IS EXAGGERATED, AND SAYS SO
    A 60 km reach and a 2 km island cannot share a scale on a 190 mm page:
    at 1:1 the map panels would be 0.22 in tall. VERTICAL_EXAGGERATION and
    VERTICAL_EXAGGERATION_ABSOLUTE set the factor for their figure, are the
    only place a panel height comes from, and both captions state them. The
    old poster stretched the cross-shore by an unstated `fig_w * aspect * 1.8`
    capped at 7.5 in.

HOUSE STYLE
    Elevation is drawn in CLASSES, not a ramp (`elevation_cmap()`): the back
    barrier sits a few decimetres below MHW and the dune is metres above it,
    so a continuous ramp renders the whole island as one flat tone. Until
    2026-09-17 this script used `plt.cm.terrain` under a `FuncNorm`, which
    also gave the figure TWO colours for water -- the -3 m sentinel fill came
    out terrain navy, the uncovered canvas came out house blue -- so the
    padded back-barrier frame read as deep ocean. One class, one colour now.

    The terrain treatment kept beside it obeys the SAME water rule: every cell
    below 0 m MHW is masked and painted C["WATER"], and the ramp is truncated
    to the 0-4 m land range it actually draws from, so its colourbar cannot
    advertise a blue no cell uses. The two schemes differ only in how they
    colour LAND (Hannah, 2026-09-17).

    The figure is sized by `figsize()` at the width it will be printed, its
    caption lives in CAPTIONS.md rather than on the canvas, and the local
    rcParams block that overrode the house ink with '#1a1a2e' is gone.

Author: Hannah Henry (UNC Chapel Hill)
==============================================================================
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

import numpy as np

# HOUSE STYLE: one typeface and one palette across every figure in this
# project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
# found by searching upward (ORGANIZATION.md rule 5), so this block is
# independent of whatever this script calls its own repository variable.
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(next(_q for _q in _P(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED, GRID_C,
                              DOMAIN_AXIS_LABEL, FIG_W_DOUBLE, elevation_cmap,
                              figsize, record_caption, save, spines_for_image,
                              open_frame, _title)
apply_style()
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap, Normalize

from cascade_pipeline.domains import DomainGeometry
from cascade_pipeline.plotting import init_planview

# =============================================================================
# CONFIGURATION
# =============================================================================

# Derived from this file's location (scripts/figure_making/island/) so the
# script runs on any checkout without editing a hardcoded path.
PROJECT_BASE_DIR   = str(next(
    _p for _p in Path(__file__).resolve().parents
    if (_p / "pyproject.toml").exists()))
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_DIR         = str(next(_q for _q in Path(__file__).resolve().parents
                              if (_q / "pyproject.toml").exists()
                              ) / "output" / "figures" / "initialization")   # rule 1

# Orientation years to render: every hindcast period start. Offsets come from
# 2-brie-offset/<year>/, the same PADDED_120 files the hindcast run script
# initializes from, and the topography from that year's product (below).
YEARS = (1984, 1996, 2004, 2010)

# WHICH EXTRACTION -- resolved, not pinned, the same way the hindcast runner and
# the groin sweep worker resolve it. This was hardcoded '2009_v2', a directory
# that has since been moved into 2009-dune-topo/incorrect/, so the poster either
# failed to load or drew arrays the run does not use. topo_dirs() reads VERSION
# out of HAT_dune_topo_extractor.py, so the figure always shows the surface the
# model is actually initialised from.
from site_layer.hat_topo_version import (topo_dirs, array_name,  # scripts/, on sys.path above
                              product_for_year)

# Paths come from what topo_dirs() RETURNED, not re-joined from parts - the
# tree went period-first on 2026-08-25 and re-joining would have kept pointing
# at a folder that no longer exists. topo_dirs() with no product resolves
# 2004-start, which is the surface this poster showed before the restructure.
from site_layer.hat_topo_version import BUFFER_DIR as _BUFFER_DIR   # noqa: E402

# array_name() is the single definition of these filenames - the same one
# the extractor writes with. Nothing here spells a name.

# ...AND THE PRODUCT IS RESOLVED PER YEAR, NOT ONCE. The four starts do not
# share one surface: 1984 and 1996 read 1984-start (the 2009+2014 mosaic with
# the 1996 ALACE graft), 2004 and 2010 read 2004-start (2009+2014 alone).
# YEAR_PRODUCT in hat_topo_version is the single definition of that pairing,
# and product_for_year() its accessor -- whose docstring names precisely the
# bug this script had until 2026-09-17: "a loop body that resolves topo_dirs()
# once, outside the loop, and silently gives every year the same interiors."
# It did, so the 1984 figure was drawing the 2004-start island.
from site_layer.hat_topo_version import DOMAIN_ROOT as _B3D_ROOT  # noqa: E402
BARRIER3D_DIR       = str(_B3D_ROOT)
BUFFER_DIR          = str(_BUFFER_DIR)

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120

START_REAL_INDEX   = NUM_BUFFER_DOMAINS        # 15  (buffer)
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS  # 105

FIRST_FILE_NUMBER  = 1

ELEV_MIN_M  = -1.0
ELEV_MAX_M  =  4.0
SEA_LEVEL   =  0.0
DAM_TO_M    = 10.0   # Barrier3D stores elevation in decameters
CELL_SIZE_M = 10.0   # ...on a 10 m grid, so offsets in metres convert to rows by the same factor
DOMAIN_M    = 500.0  # alongshore extent of one domain

# The extractor works in a 200-row (2000 m) cross-shore frame but writes the
# topography trimmed to the interior, so the back-barrier water rows are absent
# from the .npy files. Refill them with the extractor's water sentinel so every
# domain spans the full frame. Both values come from RUN_MANIFEST.txt
# (TOPO_ROWS, SENTINEL_WATER_M / WATER_CLAMP_M) in the version folder.
TOPO_ROWS        = 200
SENTINEL_WATER_M = -3.0   # metres MHW

# HOW FAR THE CROSS-SHORE IS STRETCHED. The alongshore scale is fixed by the
# page (120 domains across a 190 mm column); this is the only other free
# number in the layout, and the caption states it. 8 puts the map panels at
# about 1.75 in, enough to read a dune line against a back barrier.
VERTICAL_EXAGGERATION = 8.0

# The absolute-placement figure carries about 9 km of cross-shore against the
# detrended figure's 2, so it cannot use the same factor and still fit a page:
# 8x would make one of its panels 7.9 in tall. 3x puts them near 3 in.
VERTICAL_EXAGGERATION_ABSOLUTE = 3.0


def dune_offset_file(year):
    """Padded BRIE dune-offset CSV for one orientation year, THROUGH THE
    VERSION THE RUN READS.

    Every start's offsets were put under `<year>/v<n>/` with a CURRENT marker
    on 2026-09-15, so the flat `<year>/Island_Dune_Offsets_...csv` this used to
    read no longer exists and the poster could not be redrawn. The version is
    resolved the way the runner resolves it, so a re-version moves this with
    it rather than breaking it again.

    Since 2026-09-18 that is literally true: hat_topo_version.offset_file is
    the function the runner's config calls. The copy that stood here ignored
    HAT_OFFSET_VERSION_<year> and fell back to the newest v<n> where the
    runner raises."""
    from site_layer.hat_topo_version import offset_file
    path = offset_file(year, "padded", TOTAL_DOMAINS)
    if not path.is_file():
        raise FileNotFoundError(f"no padded offset file for {year}: {path}")
    return str(path)


# =============================================================================
# ELEVATION ARRAYS, ONE SET PER TOPOGRAPHY PRODUCT
# =============================================================================

# Compositing (padding, unit conversion, canvas assembly) is shared with the QC
# notebook via cascade_pipeline.plotting.init_planview; only the page-figure
# styling below is local to this script.
GEOMETRY = DomainGeometry(num_real_domains=NUM_REAL_DOMAINS,
                          num_buffer_domains=NUM_BUFFER_DOMAINS,
                          first_gis_id=FIRST_FILE_NUMBER,
                          domain_spacing_m=DOMAIN_M)
PLAN_VIEW = init_planview.PlanViewConfig(
    topo_rows=TOPO_ROWS, sentinel_water_m=SENTINEL_WATER_M,
    cell_size_m=CELL_SIZE_M, dam_to_m=DAM_TO_M,
    elev_min_m=ELEV_MIN_M, elev_max_m=ELEV_MAX_M, sea_level_m=SEA_LEVEL)

_PRODUCT_CACHE = {}


def elevation_file_paths(topo_dir):
    """The 120 padded-order array paths for one product's topography dir.

    The 15 domains at each end are the SAME sampled array repeated, which is
    why the buffers read as a regular comb in the absolute figure.
    """
    buffer_path = os.path.join(BUFFER_DIR, 'sample_1_topography.npy')
    paths = [buffer_path] * START_REAL_INDEX
    paths += [os.path.join(str(topo_dir), array_name('topography', n))
              for n in range(FIRST_FILE_NUMBER,
                             FIRST_FILE_NUMBER + NUM_REAL_DOMAINS)]
    paths += [buffer_path] * (TOTAL_DOMAINS - END_REAL_INDEX)
    return paths


def product_for(year):
    """(grids, product, version, interior row range) for one start year.

    Loaded once per PRODUCT, not once per year: 1984 and 1996 share a surface,
    as do 2004 and 2010, so four figures cost two loads.
    """
    product = product_for_year(year)
    if product not in _PRODUCT_CACHE:
        topo_dir, _dune_dir, version = topo_dirs(product)
        paths = elevation_file_paths(topo_dir)
        missing = [q for q in paths if not os.path.exists(q)]
        if missing:
            raise FileNotFoundError(
                f"{len(missing)} of {TOTAL_DOMAINS} elevation files missing for "
                f"{product} {version} under {topo_dir}\n"
                f"  first missing: {missing[0]}")
        print(f"Loading {TOTAL_DOMAINS} domain arrays from {product} {version}...")
        grids = init_planview.load_domain_grids(paths, PLAN_VIEW)
        rows = [np.load(q).shape[0] for q in paths[START_REAL_INDEX:END_REAL_INDEX]]
        print(f"  Interior rows {min(rows)}-{max(rows)}, padded to {TOPO_ROWS} "
              f"({TOPO_ROWS * int(CELL_SIZE_M)} m) with water at {SENTINEL_WATER_M} m")
        _PRODUCT_CACHE[product] = (grids, product, version, (min(rows), max(rows)))
    return _PRODUCT_CACHE[product]


def load_offsets_m(year):
    """The padded BRIE offsets for one year, in metres, one per domain."""
    offsets_m = np.loadtxt(dune_offset_file(year), skiprows=1, delimiter=',')
    if offsets_m.ndim != 1 or offsets_m.size != TOTAL_DOMAINS:
        raise ValueError(f"{year} offsets: expected {TOTAL_DOMAINS} values, "
                         f"got shape {offsets_m.shape}")
    return offsets_m


def detrended_canvas(year, include_buffers):
    """The composited surface with every domain on its OWN frame origin.

    Passing zero offsets to the shared compositor is exactly the detrending:
    `build_canvas` places domain i at row `offset_cells[i]`, so zeros put each
    one where its elevation array starts. The offsets themselves are drawn in
    panel (a) instead of being spent on 6 km of empty canvas.
    """
    grids = product_for(year)[0]
    zero = np.zeros(TOTAL_DOMAINS, dtype=int)
    canvas, _, _, _ = init_planview.build_canvas(
        grids, zero, GEOMETRY, include_buffers, PLAN_VIEW)
    # build_canvas leaves max(offset) + topo_rows + 5 rows; with no offsets
    # that is five rows of NaN above the frame.
    return canvas[:TOPO_ROWS]


def absolute_canvas(year, include_buffers):
    """The composited surface with every domain at its TRUE BRIE offset.

    The initial condition as one map, which is what the figure showed before
    2026-09-17 and what `island_<year>_absolute.png` shows again: the reach
    bends seaward by about 6 km from Pea Island to Cape Point, and here that
    bend is geometry rather than a curve on a separate axis.
    """
    grids = product_for(year)[0]
    offset_cells = np.round(load_offsets_m(year) / CELL_SIZE_M).astype(int)
    canvas, _, _, _ = init_planview.build_canvas(
        grids, offset_cells, GEOMETRY, include_buffers, PLAN_VIEW)
    return canvas


# =============================================================================
# THE PAGE FIGURE
# =============================================================================
# One alongshore scale, in metres per inch, shared by both map panels and the
# offset panel. Everything else in the layout is derived from it, so the panels
# cannot disagree about where a domain is.
_MARGIN_L, _MARGIN_R = 0.82, 0.26           # in
_AXW_FULL = FIG_W_DOUBLE - _MARGIN_L - _MARGIN_R          # the 120-domain width
_M_PER_IN = TOTAL_DOMAINS * DOMAIN_M / _AXW_FULL
_AXW_REAL = NUM_REAL_DOMAINS * DOMAIN_M / _M_PER_IN       # the 90-domain width
_MAP_H = (TOPO_ROWS * CELL_SIZE_M / _M_PER_IN) * VERTICAL_EXAGGERATION

# GIS-domain coordinates of the padded array: index 0 is the 15th buffer south
# of GIS 1, so domain i sits at i - NUM_BUFFER_DOMAINS + FIRST_FILE_NUMBER.
_GIS = np.arange(TOTAL_DOMAINS) - NUM_BUFFER_DOMAINS + FIRST_FILE_NUMBER
_GIS_LO, _GIS_HI = _GIS[0] - 0.5, _GIS[-1] + 0.5
_REAL_LO, _REAL_HI = FIRST_FILE_NUMBER - 0.5, FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 0.5
_XTICKS = [1, 15, 30, 45, 60, 75, 90]


def _map_panel(ax, canvas, xlim, cmap, norm, water, top_km=None,
               yticks=(0, 1, 2)):
    """One elevation strip, seaward at the bottom.

    The image extent comes from the canvas's OWN row count, so a cell always
    lands at its true cross-shore distance; `top_km` then sets the window.
    Two panels drawn from canvases of different heights -- which is what the
    real-only and with-buffers canvases are, the buffer offsets running further
    seaward than any real domain -- therefore still share one y scale.
    """
    rows_km = canvas.shape[0] * CELL_SIZE_M / 1000.0
    # WATER IS MASKED, NOT COLOURED. Every cell below 0 m MHW -- real
    # back-barrier, and the -3 m sentinel the extractor's trimmed rows are
    # refilled with -- goes to `set_bad`, which both schemes set to
    # C["WATER"]. Masking rather than letting each colormap render its own
    # low end is what makes the two treatments agree about water; it also
    # matches the class scale, whose first bound is [-99, 0), so nothing
    # changes in the class figures. Strictly BELOW: a cell at exactly 0.0 m
    # falls in the [0, 0.5) land class, and it still does here.
    canvas = np.ma.masked_less(np.ma.masked_invalid(canvas), SEA_LEVEL)
    im = ax.imshow(canvas, cmap=cmap, norm=norm, origin="lower",
                   extent=(xlim[0], xlim[1], 0.0, rows_km),
                   interpolation="nearest", aspect="auto", zorder=1)
    ax.set_facecolor(water)               # any panel edge beyond the canvas
    ax.set_xlim(*xlim)
    ax.set_ylim(0.0, top_km if top_km is not None else rows_km)
    ax.set_yticks(list(yticks))
    ax.set_ylabel("cross-shore\n(km)")
    ax.set_xticks(_XTICKS)
    spines_for_image(ax)                  # a map panel: the frame is the data edge
    return im


# =============================================================================
# THE TWO ELEVATION TREATMENTS
# =============================================================================
# The house rule is CLASSES (hat_figure_style: a hard break at 0 m, one colour
# for water, so the only distinction that matters -- which cells are land -- is
# the one you see first). The terrain ramp the figure used until 2026-09-17 is
# kept beside it, at Hannah's request, because a continuous ramp shows the
# smooth cross-shore gradient a class boundary hides, and because it is what
# the extractor's own QC view still draws.
#
# They are NOT mixed in a folder: a reader who sees both treatments of one
# island in one place has to work out that they are the same data. Each gets
# its own <year>/<scheme>/ with its own supporting/ and CAPTIONS.md, so the
# file names are identical and only the path says which is which
# (Hannah, 2026-09-17).
SCHEMES = ("classes", "terrain")

# ONE WATER RULE, BOTH SCHEMES (Hannah, 2026-09-17). Water is every cell below
# 0 m MHW and it is C["WATER"], full stop -- the class rule, applied to the
# ramp as well. `_map_panel` masks those cells and `set_bad` paints them, so
# both treatments take the identical path and cannot drift apart.
#
# Terrain used to draw its own bottom (navy) for everything under its -1 m
# floor, which on this canvas is the -3 m sentinel that refills each domain's
# trimmed rows. That put TWO blues in one image -- navy for the padded
# back-barrier frame inside a domain, pale blue for canvas no domain covers --
# and the navy read as deep ocean. Now there is one.
#
# Because water no longer comes off the ramp, the ramp must not claim it: the
# terrain colormap is TRUNCATED to the land part it actually uses, its old
# sea-level position upward, under a plain 0-4 m Normalize. Land colours are
# unchanged -- the old FuncNorm sent elevation e >= 0 to colormap position
# SEA_LEVEL_POS + (1 - SEA_LEVEL_POS) * e / ELEV_MAX_M, which is exactly what
# the truncated map under a linear norm does -- but the colourbar no longer
# shows a blue band no cell is drawn from.
SEA_LEVEL_POS = PLAN_VIEW.sea_level_pos          # 0.35, the old sea-level pin

# THE TERRAIN SCHEME KEEPS TERRAIN'S OWN NAVY (Hannah, 2026-09-17). It is
# plt.cm.terrain(0.0), the floor the ramp used to render the -3 m sentinel at,
# and it is what the figure looked like before this rewrite. The pale house
# water was tried first and sits too close in value to the saturated green the
# ramp starts at for the shoreline to read.
#
# What is NOT restored is the old two-blue split. Back then navy came off the
# ramp (the padded back-barrier frame inside a domain) while the uncovered
# canvas came off set_bad in house blue, so one image had two colours for
# water and the navy read as deep ocean. Water is masked under one rule now,
# so this navy is every cell below 0 m MHW and nothing else. The class figures
# keep C["WATER"] exactly; the schemes share the rule, not the shade.
WATER_TERRAIN = "#333399"        # plt.cm.terrain(0.0)

SCHEME_WATER = {"classes": C["WATER"], "terrain": WATER_TERRAIN}

# WHAT IS DRAWN ON TOP OF THE WATER FOLLOWS IT. The dashed lines bracketing the
# real span cross open water for most of their length, so their colour is a
# property of the water they sit on, not of the figure: near-black INK on the
# pale class blue, white on the terrain navy. Picking one ink for both would
# lose the line in whichever scheme disagreed (Hannah, 2026-09-17).
SCHEME_MARK = {"classes": INK, "terrain": "white"}

SCHEME_NOTE = {
    "classes": ("in the elevation classes of the house style, m MHW, water "
                "(below 0 m) one colour"),
    "terrain": ("on a continuous terrain ramp over the land range 0-"
                f"{ELEV_MAX_M:g} m MHW, with water (below 0 m) a single deeper "
                "blue under the same rule the class scheme uses"),
}


def scheme_colours(scheme):
    """(cmap, norm, colourbar ticks, water colour) for one elevation treatment.

    `_map_panel` masks every cell below sea level and `set_bad` paints it, so
    both schemes take the identical path to water; only the shade differs.
    """
    water = SCHEME_WATER[scheme]
    if scheme == "classes":
        cmap, norm, bounds = elevation_cmap()
        cmap = cmap.copy()
        cmap.set_bad(water)
        # bounds[0] is the water class, which no cell reaches any more; the
        # ticks start at the 0 m break either way.
        return cmap, norm, bounds[1:-1], water
    if scheme == "terrain":
        land = plt.cm.terrain(np.linspace(SEA_LEVEL_POS, 1.0, 256))
        cmap = ListedColormap(land, name="hat_terrain_land")
        cmap.set_bad(water)
        return cmap, Normalize(vmin=SEA_LEVEL, vmax=ELEV_MAX_M), [0, 1, 2, 3, 4], water
    raise ValueError(f"unknown scheme {scheme!r}; expected one of {SCHEMES}")


def year_dir(year, scheme):
    """`output/figures/initialization/<year>/<scheme>/`. `save()` creates it,
    and the `supporting/` inside it, so nothing here makes a directory."""
    return Path(OUTPUT_DIR) / str(year) / scheme


def figure_for_year(year, scheme):
    """Draw and save the detrended figure for one start year and treatment."""
    offsets_m = load_offsets_m(year)
    cmap, norm, cb_ticks, water = scheme_colours(scheme)
    mark = SCHEME_MARK[scheme]

    _, product, version, rows = product_for(year)
    real = detrended_canvas(year, include_buffers=False)
    buff = detrended_canvas(year, include_buffers=True)

    strip_h, cbar_h = 0.60, 0.11
    top, gap_a, gap_b, xlab, cbar_gap, bottom = 0.22, 0.32, 0.30, 0.42, 0.16, 0.46
    fh = (top + strip_h + gap_a + _MAP_H + gap_b + _MAP_H
          + xlab + cbar_gap + cbar_h + bottom)
    fig = plt.figure(figsize=figsize("double", height=fh))

    def add_axes(y_top_in, w_in, h_in, left_in=_MARGIN_L):
        return fig.add_axes([left_in / FIG_W_DOUBLE, (fh - y_top_in - h_in) / fh,
                             w_in / FIG_W_DOUBLE, h_in / fh])

    # ---- (a) the BRIE shoreline offset the map panels were detrended by ----
    y = top
    ax_o = add_axes(y, _AXW_FULL, strip_h)
    for lo, hi in ((_GIS_LO, _REAL_LO), (_REAL_HI, _GIS_HI)):
        ax_o.axvspan(lo, hi, color="0.94", lw=0, zorder=0)
    ax_o.plot(_GIS[:START_REAL_INDEX + 1], offsets_m[:START_REAL_INDEX + 1],
              color=INK_MUTED, lw=1.0, ls=(0, (3, 2)), zorder=2)
    ax_o.plot(_GIS[END_REAL_INDEX - 1:], offsets_m[END_REAL_INDEX - 1:],
              color=INK_MUTED, lw=1.0, ls=(0, (3, 2)), zorder=2)
    ax_o.plot(_GIS[START_REAL_INDEX:END_REAL_INDEX],
              offsets_m[START_REAL_INDEX:END_REAL_INDEX],
              color=INK, lw=1.2, zorder=3)
    ax_o.set_xlim(_GIS_LO, _GIS_HI)
    ax_o.set_xticks(_XTICKS)
    ax_o.tick_params(labelbottom=False)
    ax_o.set_yticks([0, 2000, 4000, 6000])
    ax_o.set_ylabel("BRIE offset\n(m)")
    ax_o.grid(axis="y", color=GRID_C, lw=0.5)
    ax_o.set_axisbelow(True)
    open_frame(ax_o)                      # a chart, not an image
    _title(ax_o, 0, "")

    # ---- (b) the 90 real domains, on the same metres-per-inch as (c) ----
    y += strip_h + gap_a
    ax_r = add_axes(y, _AXW_REAL, _MAP_H,
                    left_in=_MARGIN_L + (_AXW_FULL - _AXW_REAL) / 2)
    _map_panel(ax_r, real, (_REAL_LO, _REAL_HI), cmap, norm, water)
    ax_r.tick_params(labelbottom=False)
    _title(ax_r, 1, "")

    # ---- (c) the same with the interpolated buffer domains ----
    y += _MAP_H + gap_b
    ax_b = add_axes(y, _AXW_FULL, _MAP_H)
    im = _map_panel(ax_b, buff, (_GIS_LO, _GIS_HI), cmap, norm, water)
    for x in (_REAL_LO, _REAL_HI):
        ax_b.axvline(x, color=mark, lw=0.9, ls=(0, (4, 3)), zorder=6)
    for x_mid in ((_GIS_LO + _REAL_LO) / 2, (_REAL_HI + _GIS_HI) / 2):
        ax_b.text(x_mid, 0.93, f"buffer\n({NUM_BUFFER_DOMAINS} domains)",
                  transform=ax_b.get_xaxis_transform(), ha="center", va="top",
                  fontsize=7.5, color=INK, zorder=7, linespacing=1.15,
                  bbox=dict(boxstyle="square,pad=0.25", facecolor="white",
                            edgecolor="none", alpha=0.85))
    ax_b.set_xlabel(DOMAIN_AXIS_LABEL)
    _title(ax_b, 2, "")

    # ---- the shared elevation scale, in classes ----
    y += _MAP_H + xlab + cbar_gap
    cax = fig.add_axes([(_MARGIN_L + 0.30 * _AXW_FULL) / FIG_W_DOUBLE,
                        (fh - y - cbar_h) / fh,
                        0.40 * _AXW_FULL / FIG_W_DOUBLE, cbar_h / fh])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=cb_ticks)
    cb.set_label("elevation (m MHW)")
    cb.outline.set_linewidth(0.5)

    out = save(fig, year_dir(year, scheme) / f"island_{year}.png",
               vector=True, dpi=300)
    record_caption(out[0],
        f"The island as the model starts it in {year}, from the {product} {version} "
        f"extraction: every domain's elevation array at the model's 10 m resolution, "
        f"{SCHEME_NOTE[scheme]}. "
        f"(a) The BRIE shoreline offset each domain is initialised at, m, solid over "
        f"the 90 real domains and dashed over the padded buffers. (b) The 90 real "
        f"domains and (c) the same with the {NUM_BUFFER_DOMAINS} buffer domains at each "
        f"end, whose topography is one sampled array repeated and carries no survey -- "
        f"the regular comb at both ends of (c) is that repeat. Both map "
        f"panels are DETRENDED: each domain is drawn against its own frame origin, so "
        f"the offsets of (a) are not spent on canvas, and the ragged landward edge is "
        f"the island's own width, {rows[0] * int(CELL_SIZE_M)}-"
        f"{rows[1] * int(CELL_SIZE_M)} m over the real domains, against the "
        f"{TOPO_ROWS * int(CELL_SIZE_M)} m extraction frame refilled with water at "
        f"{SENTINEL_WATER_M:g} m. Seaward at the bottom; the two map panels share one "
        f"alongshore scale, and the cross-shore is exaggerated "
        f"{VERTICAL_EXAGGERATION:g}x against it.")
    plt.close(fig)
    return out[0]


# =============================================================================
# THE ABSOLUTE-PLACEMENT FIGURE
# =============================================================================
# The reach as one map: every domain at its true BRIE cross-shore offset, so
# the seaward bend from Pea Island to Cape Point is geometry rather than a
# curve on a separate axis. This is the view the figure had before 2026-09-17,
# kept because a single canvas IS the initial condition and a detrended panel
# is a rearrangement of it -- but drawn to the same rules as everything else.


def figure_absolute(year, scheme, include_buffers):
    """Draw and save ONE absolute-placement map: a single panel, its own file.

    The real-domain map and the with-buffers map were panels (a) and (b) of one
    figure until 2026-09-17. They are separate figures now (Hannah): each is a
    map of the whole reach in its own right, each is used on its own, and as a
    pair they cost the page two half-height panels to say one thing twice.

    Standing alone, each also gets the full column width and its OWN
    cross-shore window, so neither carries the other's dead space -- the
    real-domain map no longer reserves the extra kilometre the buffer offsets
    need. The price is that the two no longer share metres per inch, so the
    caption states the span and the exaggeration of each.
    """
    _, product, version, _rows = product_for(year)
    cmap, norm, cb_ticks, water = scheme_colours(scheme)
    mark = SCHEME_MARK[scheme]

    canvas = absolute_canvas(year, include_buffers)
    xlim = (_GIS_LO, _GIS_HI) if include_buffers else (_REAL_LO, _REAL_HI)
    span_domains = TOTAL_DOMAINS if include_buffers else NUM_REAL_DOMAINS
    m_per_in = span_domains * DOMAIN_M / _AXW_FULL
    top_km = canvas.shape[0] * CELL_SIZE_M / 1000.0
    map_h = (top_km * 1000.0 / m_per_in) * VERTICAL_EXAGGERATION_ABSOLUTE
    yticks = list(range(0, int(top_km) + 1, 2))

    cbar_h = 0.11
    top, xlab, cbar_gap, bottom = 0.16, 0.42, 0.16, 0.46
    fh = top + map_h + xlab + cbar_gap + cbar_h + bottom
    fig = plt.figure(figsize=figsize("double", height=fh))

    ax = fig.add_axes([_MARGIN_L / FIG_W_DOUBLE, (fh - top - map_h) / fh,
                       _AXW_FULL / FIG_W_DOUBLE, map_h / fh])
    im = _map_panel(ax, canvas, xlim, cmap, norm, water, top_km=top_km,
                    yticks=yticks)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    # No panel letter: a one-panel figure has nothing to letter against.

    if include_buffers:
        for x in (_REAL_LO, _REAL_HI):
            ax.axvline(x, color=mark, lw=0.9, ls=(0, (4, 3)), zorder=6)
        # The two buffer ramps sit at opposite ends of the cross-shore window
        # -- the south one high, the north one low -- so each tag goes to the
        # corner its own ramp leaves empty rather than onto the island.
        for x_mid, y_ax, va in (((_GIS_LO + _REAL_LO) / 2, 0.06, "bottom"),
                                ((_REAL_HI + _GIS_HI) / 2, 0.94, "top")):
            ax.text(x_mid, y_ax, f"buffer\n({NUM_BUFFER_DOMAINS} domains)",
                    transform=ax.get_xaxis_transform(), ha="center", va=va,
                    fontsize=7.5, color=INK, zorder=7, linespacing=1.15,
                    bbox=dict(boxstyle="square,pad=0.25", facecolor="white",
                              edgecolor="none", alpha=0.85))

    cax = fig.add_axes([(_MARGIN_L + 0.30 * _AXW_FULL) / FIG_W_DOUBLE,
                        (fh - top - map_h - xlab - cbar_gap - cbar_h) / fh,
                        0.40 * _AXW_FULL / FIG_W_DOUBLE, cbar_h / fh])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=cb_ticks)
    cb.set_label("elevation (m MHW)")
    cb.outline.set_linewidth(0.5)

    real_off = load_offsets_m(year)[START_REAL_INDEX:END_REAL_INDEX]
    stem = f"island_{year}_absolute" + ("_buffers" if include_buffers else "")
    out = save(fig, year_dir(year, scheme) / f"{stem}.png", vector=True, dpi=300)
    which = (f"the 90 real domains and the {NUM_BUFFER_DOMAINS} buffer domains at "
             f"each end, whose topography is one sampled array repeated and carries "
             f"no survey"
             if include_buffers else "the 90 real domains only")
    companion = (f"island_{year}_absolute.png drops the buffers"
                 if include_buffers
                 else f"island_{year}_absolute_buffers.png adds the buffer domains")
    record_caption(out[0],
        f"The island as the model starts it in {year}, from the {product} {version} "
        f"extraction, with every domain at its TRUE cross-shore position: each "
        f"domain's elevation array placed at its own BRIE shoreline offset, so this "
        f"canvas is the initial condition the run begins from, {which}. Elevation "
        f"{SCHEME_NOTE[scheme]}. Cross-shore 0 is the frame origin of the most "
        f"seaward domain; the offsets run {real_off.min():.0f}-{real_off.max():.0f} m "
        f"over the real domains, which is the seaward bend from Pea Island to Cape "
        f"Point and why the island crosses the canvas rather than running level. "
        f"Seaward at the bottom; {span_domains} domains across the page at "
        f"{m_per_in:,.0f} m per inch, with the cross-shore exaggerated "
        f"{VERTICAL_EXAGGERATION_ABSOLUTE:g}x against the alongshore. Companions: "
        f"{companion}, and island_{year}.png detrends this onto each domain's own "
        f"frame and carries the offsets as a panel of their own.")
    plt.close(fig)
    return out[0]


# =============================================================================
# RENDER
# =============================================================================

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Alongshore scale {_M_PER_IN:,.0f} m/in; map panels {_MAP_H:.2f} in "
          f"at {VERTICAL_EXAGGERATION:g}x vertical exaggeration")
    for _year in YEARS:
        print(f"Building {_year} ({product_for_year(_year)})...")
        _offsets = load_offsets_m(_year)
        _real = _offsets[START_REAL_INDEX:END_REAL_INDEX]
        print(f"  Offsets (m): min={_real.min():.1f}, max={_real.max():.1f}")
        for _scheme in SCHEMES:
            print(f"  [{_scheme}]")
            print(f"    {figure_for_year(_year, _scheme)}")
            for _buffers in (False, True):
                print(f"    {figure_absolute(_year, _scheme, _buffers)}")
