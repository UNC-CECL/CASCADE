"""
The island as the model starts it: the t=0 elevation surface of every domain, one figure per start year.

    python scripts/figure_making/island/initialization_figures.py

Detrended and absolute-placement maps, each in two elevation treatments
(classes, terrain); topography and offsets resolved per start year. Writes to
output/figures/3-model-inputs/1-domains/initial_island/<year>/<scheme>/. Details: scripts/figure_making/island/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

import numpy as np

# House style (site_layer/hat_figure_style.py), applied at import
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


# Derived from this file's location, so any checkout works
# --- CONFIG ------------------------------------------------------------------
PROJECT_BASE_DIR   = str(next(
    _p for _p in Path(__file__).resolve().parents
    if (_p / "pyproject.toml").exists()))
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
from site_layer.hat_figure_style import figure_dir as _figure_dir  # noqa: E402
OUTPUT_DIR         = str(_figure_dir("inputs", "1-domains", "initial_island"))   # rule 1

# Orientation years: every hindcast period start
YEARS = (1984, 1996, 2004, 2010)

# Which extraction: resolved like the runner resolves it, never pinned
from site_layer.hat_topo_version import (topo_dirs, array_name,  # scripts/, on sys.path above
                              product_for_year)

# Paths come from what topo_dirs() returned, not re-joined from parts
from site_layer.hat_topo_version import BUFFER_DIR as _BUFFER_DIR   # noqa: E402

# array_name() is the single definition of these filenames, shared with the extractor

# The topography product is resolved per year: 1984/1996 and 2004/2010 differ (README)
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

# Refill the trimmed back-barrier rows with the extractor's water sentinel
TOPO_ROWS        = 200
SENTINEL_WATER_M = -3.0   # metres MHW

# Cross-shore stretch of the detrended maps (the caption states it)
VERTICAL_EXAGGERATION = 8.0

# The absolute-placement maps need less stretch to fit a page
VERTICAL_EXAGGERATION_ABSOLUTE = 3.0
# -----------------------------------------------------------------------------


# Padded BRIE dune-offset CSV for one year, at the version the run reads
def dune_offset_file(year):
    from site_layer.hat_topo_version import offset_file
    path = offset_file(year, "padded", TOTAL_DOMAINS)
    if not path.is_file():
        raise FileNotFoundError(f"no padded offset file for {year}: {path}")
    return str(path)


# Elevation arrays, one set per topography product

# Compositing is shared with the QC notebook (init_planview); styling is local
GEOMETRY = DomainGeometry(num_real_domains=NUM_REAL_DOMAINS,
                          num_buffer_domains=NUM_BUFFER_DOMAINS,
                          first_gis_id=FIRST_FILE_NUMBER,
                          domain_spacing_m=DOMAIN_M)
PLAN_VIEW = init_planview.PlanViewConfig(
    topo_rows=TOPO_ROWS, sentinel_water_m=SENTINEL_WATER_M,
    cell_size_m=CELL_SIZE_M, dam_to_m=DAM_TO_M,
    elev_min_m=ELEV_MIN_M, elev_max_m=ELEV_MAX_M, sea_level_m=SEA_LEVEL)

_PRODUCT_CACHE = {}


# The 120 padded-order array paths for one product's topography dir
def elevation_file_paths(topo_dir):
    buffer_path = os.path.join(BUFFER_DIR, 'sample_1_topography.npy')
    paths = [buffer_path] * START_REAL_INDEX
    paths += [os.path.join(str(topo_dir), array_name('topography', n))
              for n in range(FIRST_FILE_NUMBER,
                             FIRST_FILE_NUMBER + NUM_REAL_DOMAINS)]
    paths += [buffer_path] * (TOTAL_DOMAINS - END_REAL_INDEX)
    return paths


# (grids, product, version, interior row range) for one start year
def product_for(year):
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


# The padded BRIE offsets for one year, in metres, one per domain
def load_offsets_m(year):
    offsets_m = np.loadtxt(dune_offset_file(year), skiprows=1, delimiter=',')
    if offsets_m.ndim != 1 or offsets_m.size != TOTAL_DOMAINS:
        raise ValueError(f"{year} offsets: expected {TOTAL_DOMAINS} values, "
                         f"got shape {offsets_m.shape}")
    return offsets_m


# The composited surface with every domain on its OWN frame origin
def detrended_canvas(year, include_buffers):
    grids = product_for(year)[0]
    zero = np.zeros(TOTAL_DOMAINS, dtype=int)
    canvas, _, _, _ = init_planview.build_canvas(
        grids, zero, GEOMETRY, include_buffers, PLAN_VIEW)
    # build_canvas leaves five NaN rows above the frame when offsets are zero
    return canvas[:TOPO_ROWS]


# The composited surface with every domain at its TRUE BRIE offset
def absolute_canvas(year, include_buffers):
    grids = product_for(year)[0]
    offset_cells = np.round(load_offsets_m(year) / CELL_SIZE_M).astype(int)
    canvas, _, _, _ = init_planview.build_canvas(
        grids, offset_cells, GEOMETRY, include_buffers, PLAN_VIEW)
    return canvas


# The page figure

# One alongshore scale (m per inch) shared by every panel
_MARGIN_L, _MARGIN_R = 0.82, 0.26           # in
_AXW_FULL = FIG_W_DOUBLE - _MARGIN_L - _MARGIN_R          # the 120-domain width
_M_PER_IN = TOTAL_DOMAINS * DOMAIN_M / _AXW_FULL
_AXW_REAL = NUM_REAL_DOMAINS * DOMAIN_M / _M_PER_IN       # the 90-domain width
_MAP_H = (TOPO_ROWS * CELL_SIZE_M / _M_PER_IN) * VERTICAL_EXAGGERATION

# GIS coordinates of the padded array: index 0 is the 15th buffer south of GIS 1
_GIS = np.arange(TOTAL_DOMAINS) - NUM_BUFFER_DOMAINS + FIRST_FILE_NUMBER
_GIS_LO, _GIS_HI = _GIS[0] - 0.5, _GIS[-1] + 0.5
_REAL_LO, _REAL_HI = FIRST_FILE_NUMBER - 0.5, FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 0.5
_XTICKS = [1, 15, 30, 45, 60, 75, 90]


# One elevation strip, seaward at the bottom
def _map_panel(ax, canvas, xlim, cmap, norm, water, top_km=None,
               yticks=(0, 1, 2)):
    rows_km = canvas.shape[0] * CELL_SIZE_M / 1000.0
    # Water is masked, not coloured: every cell below 0 m MHW (README)
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


# Two treatments: classes (house rule) and the terrain ramp, in separate folders (README)
SCHEMES = ("classes", "terrain")

# One water rule for both schemes; the terrain ramp is truncated to land (README)
SEA_LEVEL_POS = PLAN_VIEW.sea_level_pos          # 0.35, the old sea-level pin

# The terrain scheme keeps terrain's own navy for water (README)
WATER_TERRAIN = "#333399"        # plt.cm.terrain(0.0)

SCHEME_WATER = {"classes": C["WATER"], "terrain": WATER_TERRAIN}

# Lines over water take the water's contrast colour in each scheme
SCHEME_MARK = {"classes": INK, "terrain": "white"}

SCHEME_NOTE = {
    "classes": ("in the elevation classes of the house style, m MHW, water "
                "(below 0 m) one colour"),
    "terrain": ("on a continuous terrain ramp over the land range 0-"
                f"{ELEV_MAX_M:g} m MHW, with water (below 0 m) a single deeper "
                "blue under the same rule the class scheme uses"),
}


# (cmap, norm, colourbar ticks, water colour) for one elevation treatment
def scheme_colours(scheme):
    water = SCHEME_WATER[scheme]
    if scheme == "classes":
        cmap, norm, bounds = elevation_cmap()
        cmap = cmap.copy()
        cmap.set_bad(water)
        # Ticks start at the 0 m break; the water class is no longer reached
        return cmap, norm, bounds[1:-1], water
    if scheme == "terrain":
        land = plt.cm.terrain(np.linspace(SEA_LEVEL_POS, 1.0, 256))
        cmap = ListedColormap(land, name="hat_terrain_land")
        cmap.set_bad(water)
        return cmap, Normalize(vmin=SEA_LEVEL, vmax=ELEV_MAX_M), [0, 1, 2, 3, 4], water
    raise ValueError(f"unknown scheme {scheme!r}; expected one of {SCHEMES}")


# Output folder for one year and scheme (initial_island/<year>/<scheme>/)
def year_dir(year, scheme):
    return Path(OUTPUT_DIR) / str(year) / scheme


# Draw and save the detrended figure for one start year and treatment
def figure_for_year(year, scheme):
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

    # (a) the BRIE shoreline offset the map panels were detrended by
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

    # (b) the 90 real domains, on the same metres-per-inch as (c)
    y += strip_h + gap_a
    ax_r = add_axes(y, _AXW_REAL, _MAP_H,
                    left_in=_MARGIN_L + (_AXW_FULL - _AXW_REAL) / 2)
    _map_panel(ax_r, real, (_REAL_LO, _REAL_HI), cmap, norm, water)
    ax_r.tick_params(labelbottom=False)
    _title(ax_r, 1, "")

    # (c) the same with the interpolated buffer domains
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

    # the shared elevation scale, in classes
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


# The reach as one map: every domain at its true BRIE offset


# Draw and save one absolute-placement map
def figure_absolute(year, scheme, include_buffers):
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
        # Each buffer tag goes to the corner its own ramp leaves empty
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


# Run: every start year, both treatments, detrended and absolute maps
def main():
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


if __name__ == "__main__":
    main()
