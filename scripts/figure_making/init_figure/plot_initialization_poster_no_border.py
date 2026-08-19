"""
plot_initialization_poster.py
------------------------------
Poster-quality plan-view plot of CASCADE initialization (t=0) for Hatteras Island.
Shows REAL DOMAINS ONLY (no buffer domains) with clean scientific styling.

Renders one figure per orientation year in YEARS. The 2009 topography is the
same initial surface for both periods; only the BRIE dune offsets differ, so
the elevation arrays are loaded once and re-composited per year.

Author: Hannah Henry (UNC Chapel Hill)
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # scripts/

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.colors import TwoSlopeNorm, FuncNorm
from matplotlib.gridspec import GridSpec

from cascade_pipeline.domains import DomainGeometry
from cascade_pipeline.plotting import init_planview

# =============================================================================
# CONFIGURATION
# =============================================================================

# Derived from this file's location (scripts/figure_making/init_figure/) so the
# script runs on any checkout without editing a hardcoded path.
PROJECT_BASE_DIR   = str(Path(__file__).resolve().parents[3])
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_DIR         = str(Path(__file__).resolve().parent)

# Orientation years to render. Offsets come from 2-brie-offset/hindcast_<year>/,
# the same PADDED_120 files the hindcast run script initializes from.
YEARS = (1984, 2004)

# Layout: 1-barrier3d-domains/2009-dune-topo/<version>/topography/domain_<n>_topography_<year>.npy
TOPO_DUNE_INIT_YEAR = '2009'      # year label in the filename
TOPO_DUNE_VERSION   = '2009_v2'   # version folder; matches the hindcast run script
BARRIER3D_DIR       = os.path.join(HATTERAS_DATA_BASE, '1-barrier3d-domains')
DUNE_TOPO_DIR       = os.path.join(BARRIER3D_DIR, '2009-dune-topo', TOPO_DUNE_VERSION)
BUFFER_DIR          = os.path.join(BARRIER3D_DIR, '2009-buffer')

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120

START_REAL_INDEX   = NUM_BUFFER_DOMAINS        # 15  (buffer)
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS  # 105

FIRST_FILE_NUMBER  = 1

# Road domains (1-indexed real domain numbering)
# FIRST_ROAD_DOMAIN  = 9
# LAST_ROAD_DOMAIN   = 90
# Convert to global domain indices
# START_ROAD_INDEX   = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS
# END_ROAD_INDEX     = (LAST_ROAD_DOMAIN  - 1) + NUM_BUFFER_DOMAINS

ELEV_MIN_M  = -1.0
ELEV_MAX_M  =  4.0
SEA_LEVEL   =  0.0
DAM_TO_M    = 10.0   # Barrier3D stores elevation in decameters
CELL_SIZE_M = 10.0   # ...on a 10 m grid, so offsets in metres convert to rows by the same factor

# The extractor works in a 200-row (2000 m) cross-shore frame but writes the
# topography trimmed to the interior, so the back-barrier water rows are absent
# from the .npy files. Refill them with the extractor's water sentinel so every
# domain spans the full frame. Both values come from RUN_MANIFEST.txt
# (TOPO_ROWS, SENTINEL_WATER_M / WATER_CLAMP_M) in the version folder.
TOPO_ROWS        = 200
SENTINEL_WATER_M = -3.0   # metres MHW


def dune_offset_file(year):
    """Padded BRIE dune-offset CSV for one orientation year."""
    return os.path.join(HATTERAS_DATA_BASE, '2-brie-offset', f'hindcast_{year}',
                        f'Island_Dune_Offsets_{year}_PADDED_{TOTAL_DOMAINS}.csv')


# =============================================================================
# FILE PATHS  (build full list including buffers for raw_offset loading)
# =============================================================================

ELEVATION_FILE_PATHS = []

for _ in range(START_REAL_INDEX):
    ELEVATION_FILE_PATHS.append(os.path.join(BUFFER_DIR, 'sample_1_topography.npy'))

for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
    ELEVATION_FILE_PATHS.append(os.path.join(
        DUNE_TOPO_DIR, 'topography',
        f'domain_{file_num}_topography_{TOPO_DUNE_INIT_YEAR}.npy'))

for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
    ELEVATION_FILE_PATHS.append(os.path.join(BUFFER_DIR, 'sample_1_topography.npy'))

_missing = [p for p in ELEVATION_FILE_PATHS if not os.path.exists(p)]
if _missing:
    raise FileNotFoundError(
        f"{len(_missing)} of {TOTAL_DOMAINS} elevation files missing - check "
        f"TOPO_DUNE_VERSION ({TOPO_DUNE_VERSION!r}) and DUNE_TOPO_DIR.\n"
        f"  first missing: {_missing[0]}"
    )

# =============================================================================
# LOAD ELEVATION ARRAYS (all domains, for consistent canvas height)
# =============================================================================

# Compositing (padding, unit conversion, canvas assembly) is shared with the QC
# notebook via cascade_pipeline.plotting.init_planview; only the poster-specific
# styling below is local to this script.
GEOMETRY = DomainGeometry(num_real_domains=NUM_REAL_DOMAINS,
                          num_buffer_domains=NUM_BUFFER_DOMAINS,
                          first_gis_id=FIRST_FILE_NUMBER,
                          domain_spacing_m=500.0)
PLAN_VIEW = init_planview.PlanViewConfig(
    topo_rows=TOPO_ROWS, sentinel_water_m=SENTINEL_WATER_M,
    cell_size_m=CELL_SIZE_M, dam_to_m=DAM_TO_M,
    elev_min_m=ELEV_MIN_M, elev_max_m=ELEV_MAX_M, sea_level_m=SEA_LEVEL)

print(f"Loading {TOTAL_DOMAINS} domain elevation arrays from {TOPO_DUNE_VERSION}...")

domain_grids = init_planview.load_domain_grids(ELEVATION_FILE_PATHS, PLAN_VIEW)

_interior_rows = [np.load(p).shape[0] for p in ELEVATION_FILE_PATHS[START_REAL_INDEX:END_REAL_INDEX]]
print(f"  Interior rows {min(_interior_rows)}-{max(_interior_rows)}, "
      f"padded to {TOPO_ROWS} ({TOPO_ROWS * int(CELL_SIZE_M)} m) with water at {SENTINEL_WATER_M} m")

# =============================================================================
# BUILD COMPOSITE CANVAS
# =============================================================================

def build_canvas(year, include_buffers=False):
    """Composites domains onto one canvas for a given orientation year.

    Thin wrapper over cascade_pipeline.plotting.init_planview.build_canvas that
    adds the year's offset file loading and validation.

    Args:
        year: Orientation year, keying into 2-brie-offset/hindcast_<year>/.
        include_buffers: Whether to include the 15-domain buffers on each end.

    Returns:
        A (canvas, domain_col_starts, cells_per_domain, first_real_idx) tuple.
    """
    offsets_m = np.loadtxt(dune_offset_file(year), skiprows=1, delimiter=',')
    if offsets_m.ndim != 1 or offsets_m.size != TOTAL_DOMAINS:
        raise ValueError(f"{year} offsets: expected {TOTAL_DOMAINS} values, got shape {offsets_m.shape}")

    offset_cells = np.round(offsets_m / CELL_SIZE_M).astype(int)
    canvas, col_starts, cells, first_real = init_planview.build_canvas(
        domain_grids, offset_cells, GEOMETRY, include_buffers, PLAN_VIEW)

    plotted = offset_cells if include_buffers else offset_cells[
        GEOMETRY.start_real_index:GEOMETRY.end_real_index]
    print(f"  Offsets (m): min={min(plotted) * CELL_SIZE_M:.1f}, "
          f"max={max(plotted) * CELL_SIZE_M:.1f}  ({len(cells)} domains)")
    return canvas, col_starts, cells, first_real


# Road column extents within the real-domain canvas
# Road starts at real domain index (FIRST_ROAD_DOMAIN - 1) within real_grids
# road_local_start = FIRST_ROAD_DOMAIN - 1   # 0-indexed within real domains
# road_local_end   = LAST_ROAD_DOMAIN - 1

# road_col_start = domain_col_starts[road_local_start]
# road_col_end   = domain_col_starts[road_local_end] + cells_per_domain[road_local_end] - 1

# =============================================================================
# POSTER-QUALITY PLOT
# =============================================================================

plt.rcParams.update({
    'font.family':      'DejaVu Sans',
    'font.size':        11,
    'axes.linewidth':   1.0,
    'xtick.direction':  'out',
    'ytick.direction':  'out',
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
    'text.color':       '#1a1a2e',
    'axes.labelcolor':  '#1a1a2e',
    'xtick.color':      '#1a1a2e',
    'ytick.color':      '#1a1a2e',
})

# ---- Colormap ----
# Keep the full 'terrain' colormap (blue → teal → green → tan → white) but shift
# sea level to map at position ~0.35 instead of 0.5, so green falls on low
# above-sea-level terrain rather than being buried in the sub-sea-level range.
SEA_LEVEL_POS = 0.35   # colormap position for sea level (0.35 → sits in the green band)


def _forward(x):
    result = np.where(
        x < SEA_LEVEL,
        SEA_LEVEL_POS * (x - ELEV_MIN_M) / (SEA_LEVEL - ELEV_MIN_M),
        SEA_LEVEL_POS + (1.0 - SEA_LEVEL_POS) * x / ELEV_MAX_M,
    )
    result = np.where(np.isnan(x), np.nan, result)   # preserve NaN so set_bad fires
    return result


def _inverse(x):
    return np.where(
        x < SEA_LEVEL_POS,
        ELEV_MIN_M + (x / SEA_LEVEL_POS) * (SEA_LEVEL - ELEV_MIN_M),
        (x - SEA_LEVEL_POS) / (1.0 - SEA_LEVEL_POS) * ELEV_MAX_M,
    )


def plot_poster(year, canvas, domain_col_starts, cells_per_domain, first_real_idx,
                include_buffers=False):
    """Render and save the poster figure for one orientation year."""
    # Aspect: keep cross-shore visually readable — not too squashed
    n_cs   = canvas.shape[0]
    n_al   = canvas.shape[1]
    aspect = n_cs / n_al   # natural cell aspect
    fig_w  = 20
    fig_h  = max(4.5, fig_w * aspect * 1.8)   # slightly exaggerated cross-shore for readability
    fig_h  = min(fig_h, 7.5)                  # cap so it fits on a poster

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')
    ax  = fig.add_axes([0.06, 0.18, 0.88, 0.68])   # [left, bottom, width, height]
    ax.set_facecolor('#b0cfe8')   # match set_bad ocean color for any remaining NaN bleed-through

    cmap = plt.cm.terrain.copy()
    cmap.set_bad(color='#b0cfe8')   # light blue for NaN/ocean cells

    norm = FuncNorm((_forward, _inverse), vmin=ELEV_MIN_M, vmax=ELEV_MAX_M)

    im = ax.pcolormesh(np.ma.masked_invalid(canvas), cmap=cmap, norm=norm, shading='auto', rasterized=True)

    # ---- Colorbar ----
    cax  = fig.add_axes([0.955, 0.18, 0.013, 0.68])
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label('Elevation (m MSL)', fontsize=12, color='#1a1a2e', labelpad=10, rotation=270)
    cbar.ax.yaxis.set_tick_params(color='#1a1a2e', labelcolor='#1a1a2e')
    cbar.outline.set_edgecolor('#cccccc')
    cbar.set_ticks([-1, 0, 1, 2, 3, 4])

    # ---- Road extent ----
    # road_lw = 1.8
    # road_color = '#c0392b'   # deep red — clear on white, matches scientific figure convention
    # ax.axvline(road_col_start - 0.5, color=road_color, lw=road_lw, ls='--', zorder=5, alpha=0.85)
    # ax.axvline(road_col_end   + 0.5, color=road_color, lw=road_lw, ls='--', zorder=5, alpha=0.85)

    # Subtle road span shading
    # ax.axvspan(road_col_start - 0.5, road_col_end + 0.5,
               # color=road_color, alpha=0.04, zorder=3)

    # (NC-12 label shown in legend)

    # ---- X axis: real domain labels (1-indexed) every 5 ----
    tick_real_indices = list(range(0, NUM_REAL_DOMAINS, 5))   # 0, 5, 10, ... within the real span
    tick_cols  = [domain_col_starts[first_real_idx + i] + cells_per_domain[first_real_idx + i] // 2
                  for i in tick_real_indices]
    tick_labels = [str(i + 1) for i in tick_real_indices]     # 1-indexed display

    ax.set_xticks(tick_cols)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.set_xlabel('Domain (S → N,  Cape Hatteras to Rodanthe)', fontsize=12, labelpad=8)

    # ---- Buffer zones: bracket the real span and label the interpolated ends ----
    if include_buffers:
        last_real_idx = first_real_idx + NUM_REAL_DOMAINS - 1
        real_col_start = domain_col_starts[first_real_idx]
        real_col_end   = domain_col_starts[last_real_idx] + cells_per_domain[last_real_idx]

        for x in (real_col_start, real_col_end):
            ax.axvline(x, color='#555555', lw=1.2, ls='--', alpha=0.8, zorder=6)

        # S label sits low (its buffer domains ride high on the canvas), N label high.
        for x_mid, side, y, va in ((real_col_start / 2, 'S', 0.03, 'bottom'),
                                   ((real_col_end + canvas.shape[1]) / 2, 'N', 0.97, 'top')):
            ax.text(x_mid, y, f'{side} buffer\n({NUM_BUFFER_DOMAINS} domains, interpolated)',
                    transform=ax.get_xaxis_transform(), ha='center', va=va,
                    fontsize=9, color='#333333', zorder=7,
                    bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                              edgecolor='#cccccc', alpha=0.9))

    # Clamp axes limits to canvas edges to eliminate white margin bars
    ax.set_xlim(0, canvas.shape[1])
    ax.set_ylim(0, canvas.shape[0])

    # ---- Y axis: qualitative cross-shore labels ----
    # n_rows = canvas.shape[0]
    # ax.set_yticks([0, n_rows * 0.15, n_rows * 0.55, n_rows * 0.90])
    # ax.set_yticklabels(['Ocean', 'Shoreline /\ndune toe', 'Dune crest /\ninterior', 'Back-\nbarrier'], fontsize=8.5)
    # ax.tick_params(axis='y', pad=4)

    # ---- Spine styling ----
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_color('#999999')

    # ---- Legend ----
    # legend_handles = [
        # mpatches.Patch(facecolor=road_color, alpha=0.7,
                       # label=f'NC-12 corridor  (domains {FIRST_ROAD_DOMAIN}–{LAST_ROAD_DOMAIN})'),
    # ]
    # leg = ax.legend(handles=legend_handles, loc='upper right', fontsize=10,
                    # framealpha=0.9, facecolor='white', edgecolor='#cccccc',
                    # labelcolor='#1a1a2e')

    # ---- Title ----
    subtitle = (f'  |  with {NUM_BUFFER_DOMAINS}-domain buffers' if include_buffers else '')
    ax.set_title(
        f'Hatteras Island — CASCADE Initialization  |  {year}{subtitle}',
        fontsize=14, fontweight='bold', color='#1a1a2e', pad=12,
    )

    # ---- Subtle domain grid lines every 10 domains ----
    for i in range(0, NUM_REAL_DOMAINS, 10):
        x = domain_col_starts[first_real_idx + i] - 0.5
        ax.axvline(x, color='#aaaaaa', lw=0.4, alpha=0.5, zorder=2)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    variant  = 'with_buffers' if include_buffers else 'concise'
    out_path = os.path.join(OUTPUT_DIR, f'initialization_poster_{year}_{variant}.png')
    fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return out_path


# =============================================================================
# RENDER
# =============================================================================

for _year in YEARS:
    for _with_buffers in (False, True):
        _label = 'with buffers' if _with_buffers else 'real domains only'
        print(f"Building {_year} canvas ({_label})...")
        _canvas, _col_starts, _cells, _first_real = build_canvas(_year, include_buffers=_with_buffers)
        print(f"  Canvas shape: {_canvas.shape}")
        print(f"  Saved: {plot_poster(_year, _canvas, _col_starts, _cells, _first_real, _with_buffers)}")
