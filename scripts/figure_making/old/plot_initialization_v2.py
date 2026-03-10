"""
plot_initialization.py
----------------------
Static plan-view plot of CASCADE initialization (t=0) for Hatteras Island.
Mirrors the domain assembly logic of FIX_gif_plot_script_improved.py exactly:
  - np.fliplr applied to each domain
  - Domains placed at canvas row = dune_offset (equivalent to x_s)
  - pcolormesh so row 0 is at bottom (ocean at bottom, inland at top)

Author: Hannah Henry (UNC Chapel Hill)
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")  # change to "TkAgg" for interactive window
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_BASE_DIR   = r'/'
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_DIR         = os.path.join(PROJECT_BASE_DIR, 'scripts', 'analyze_outputs', 'plots')
OUTPUT_NAME        = 'initialization_overview_1978_interpolated.png'

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120

START_REAL_INDEX   = NUM_BUFFER_DOMAINS       # 15
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS  # 105
FIRST_FILE_NUMBER  = 1

YEAR_COLUMN_INDEX  = 0   # 0 = 1978, 1 = 1997
YEAR_LABEL         = 1978

FIRST_ROAD_DOMAIN  = 9
LAST_ROAD_DOMAIN   = 90
START_ROAD_INDEX   = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS
END_ROAD_INDEX     = (LAST_ROAD_DOMAIN  - 1) + NUM_BUFFER_DOMAINS

DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE, 'island_offset', 'hindcast_1978_1997',
    f'Island_Dune_Offsets_1978_1997_PADDED_{TOTAL_DOMAINS}.csv'
)

ELEV_MIN_M = -1.0
ELEV_MAX_M =  4.0
SEA_LEVEL  =  0.0
DAM_TO_M   = 10.0

# =============================================================================
# FILE PATHS
# =============================================================================

ELEVATION_FILE_PATHS = []

for _ in range(START_REAL_INDEX):
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy'))

for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
    ELEVATION_FILE_PATHS.append(os.path.join(
        HATTERAS_DATA_BASE, 'topography', '2009',
        f'domain_{file_num}_topography_2009.npy'))

for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy'))

# =============================================================================
# LOAD DUNE OFFSETS  (equivalent to x_s in gif_plot_script)
# =============================================================================

print("Loading dune offsets...")
dune_offset_all   = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=',')
dune_offset_m     = dune_offset_all[:, YEAR_COLUMN_INDEX]           # meters
dune_offset_cells = np.round(dune_offset_m / DAM_TO_M).astype(int)  # cells (1 cell = 1 dam)

print(f"  Offsets (m): min={dune_offset_m.min():.1f}, max={dune_offset_m.max():.1f}")
print(f"  Offset cells: min={dune_offset_cells.min()}, max={dune_offset_cells.max()}")

# =============================================================================
# LOAD ELEVATION ARRAYS
# =============================================================================

print(f"Loading {TOTAL_DOMAINS} domain elevation arrays...")

domain_grids   = []
domain_offsets = []

for i, elev_path in enumerate(ELEVATION_FILE_PATHS):
    try:
        arr = np.load(elev_path).astype(float)
    except FileNotFoundError:
        print(f"  WARNING domain {i}: not found — {os.path.basename(elev_path)}")
        arr = np.zeros((30, 1))
    except Exception as e:
        print(f"  WARNING domain {i}: {e}")
        arr = np.zeros((30, 1))

    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)

    domain_grids.append(arr * DAM_TO_M)   # dam → m
    domain_offsets.append(int(dune_offset_cells[i]) if i < len(dune_offset_cells) else 0)

print(f"  Cross-shore rows: {min(g.shape[0] for g in domain_grids)}–{max(g.shape[0] for g in domain_grids)}")

# =============================================================================
# BUILD COMPOSITE CANVAS  — mirrors FIX_gif_plot_script_improved.py exactly
#
# gif_plot_script key lines:
#   Domain = np.fliplr(Domain)
#   OriginTstart = int(x_s)          ← dune_offset is our x_s
#   AnimateDomain[OriginTstart:OriginTstop, xOrigin:xOrigin+BarrierLength] = Domain
#   pcolormesh(AnimateDomain) → row 0 at bottom = ocean
# =============================================================================

max_cs_raw  = max(g.shape[0] for g in domain_grids)
max_offset  = max(domain_offsets)
canvas_rows = max_offset + max_cs_raw + 5   # same padding logic as AniDomainWidth

total_cols  = sum(g.shape[1] for g in domain_grids)
canvas      = np.full((canvas_rows, total_cols), np.nan)

col_cursor        = 0
domain_col_starts = []
cells_per_domain  = []

for i, grid in enumerate(domain_grids):
    n_rows, n_cols = grid.shape
    domain_col_starts.append(col_cursor)
    cells_per_domain.append(n_cols)

    # Match gif_plot_script: np.fliplr(Domain)
    grid = np.fliplr(grid)

    # Place at row = offset, matching OriginTstart = int(x_s) in gif script
    origin   = domain_offsets[i]
    row_end  = origin + n_rows
    if row_end <= canvas_rows:
        canvas[origin:row_end, col_cursor:col_cursor + n_cols] = grid
    else:
        canvas[origin:canvas_rows, col_cursor:col_cursor + n_cols] = grid[:canvas_rows - origin, :]

    col_cursor += n_cols

domain_col_starts = np.array(domain_col_starts)

print(f"Canvas shape: {canvas.shape}  (cross-shore × alongshore cells)")
print(f"  Row 0 = ocean side (bottom of plot), row {canvas_rows-1} = inland")

# =============================================================================
# PLOT — pcolormesh so row 0 is at bottom (ocean at bottom, matching gif script)
# =============================================================================

fig_width  = min(30, max(14, TOTAL_DOMAINS * 0.22))
fig_height = 7
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

cmap = plt.cm.terrain.copy()
cmap.set_bad(color='#1a1a5e')  # NaN = deep blue (ocean)

norm = TwoSlopeNorm(vmin=ELEV_MIN_M, vcenter=SEA_LEVEL, vmax=ELEV_MAX_M)

# pcolormesh: row 0 naturally at bottom → ocean at bottom
im = ax.pcolormesh(canvas, cmap=cmap, norm=norm, shading='auto')

cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cbar.set_label('Elevation (m MSL)', fontsize=11)

n_rows = canvas.shape[0]

# --- Buffer shading ---
buf_color       = 'steelblue'
left_end_col    = domain_col_starts[START_REAL_INDEX]
right_start_col = domain_col_starts[END_REAL_INDEX]

ax.axvspan(-0.5,                   left_end_col - 0.5,     color=buf_color, alpha=0.15, zorder=3)
ax.axvspan(right_start_col - 0.5,  canvas.shape[1] - 0.5,  color=buf_color, alpha=0.15, zorder=3)
ax.axvline(left_end_col - 0.5,     color=buf_color, lw=2.0, ls='--', zorder=4)
ax.axvline(right_start_col - 0.5,  color=buf_color, lw=2.0, ls='--', zorder=4)

ax.text(left_end_col / 2, n_rows * 0.95,
        f'S Buffer\n(0–{START_REAL_INDEX - 1})',
        ha='center', va='top', fontsize=8, color='white', style='italic')
ax.text((right_start_col + canvas.shape[1]) / 2, n_rows * 0.95,
        f'N Buffer\n({END_REAL_INDEX}–{TOTAL_DOMAINS - 1})',
        ha='center', va='top', fontsize=8, color='white', style='italic')

# --- Road extent ---
road_col_start = domain_col_starts[START_ROAD_INDEX]
road_col_end   = domain_col_starts[END_ROAD_INDEX] + cells_per_domain[END_ROAD_INDEX] - 1
ax.axvline(road_col_start - 0.5, color='saddlebrown', lw=1.5, ls=':', zorder=4)
ax.axvline(road_col_end   + 0.5, color='saddlebrown', lw=1.5, ls=':', zorder=4)

# --- X axis: domain ticks every 5 ---
tick_domains = list(range(0, TOTAL_DOMAINS, 5))
tick_cols    = [domain_col_starts[i] + cells_per_domain[i] // 2 for i in tick_domains]
ax.set_xticks(tick_cols)
ax.set_xticklabels([str(i) for i in tick_domains], fontsize=7, rotation=45)
ax.set_xlabel(f'Domain index  (S → N)  |  real domains {START_REAL_INDEX}–{END_REAL_INDEX - 1}', fontsize=11)

# --- Y axis ---
ax.set_yticks([0, n_rows // 3, 2 * n_rows // 3, n_rows - 1])
ax.set_yticklabels(['Ocean', 'Dune/beach', 'Interior', 'Back-barrier'], fontsize=9)

# --- Legend ---
legend_handles = [
    mpatches.Patch(facecolor=buf_color, alpha=0.35, label=f'Buffer domains (±{NUM_BUFFER_DOMAINS})'),
    mpatches.Patch(facecolor='saddlebrown', alpha=0.5,
                   label=f'Road extent (domains {FIRST_ROAD_DOMAIN}–{LAST_ROAD_DOMAIN})'),
]
ax.legend(handles=legend_handles, loc='upper right', fontsize=9, framealpha=0.75)

ax.set_title(
    f'CASCADE Initialization — Hatteras Island  |  {YEAR_LABEL}  |  '
    f'{TOTAL_DOMAINS} total domains  ({NUM_BUFFER_DOMAINS} buf + {NUM_REAL_DOMAINS} real + {NUM_BUFFER_DOMAINS} buf)\n'
    f'Elevation (m MSL)  |  ocean at bottom, back-barrier at top  |  domain assembly matches gif_plot_script',
    fontsize=12,
)

plt.tight_layout()

# =============================================================================
# SAVE
# =============================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)
out_path = os.path.join(OUTPUT_DIR, OUTPUT_NAME)
fig.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\n✓ Saved to:\n  {out_path}")
plt.close(fig)
