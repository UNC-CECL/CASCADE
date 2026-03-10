"""
plot_initialization_poster.py
------------------------------
Poster-quality plan-view plot of CASCADE initialization (t=0) for Hatteras Island.
Shows REAL DOMAINS ONLY (no buffer domains) with clean scientific styling.

Author: Hannah Henry (UNC Chapel Hill)
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_BASE_DIR   = r'C:\Users\hanna\PycharmProjects\CASCADE'
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_DIR         = os.path.join(PROJECT_BASE_DIR, 'scripts', 'figure_making', 'plots')
OUTPUT_NAME        = 'initialization_poster_1978_concise.png'

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120

START_REAL_INDEX   = NUM_BUFFER_DOMAINS        # 15  (buffer)
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS  # 105

FIRST_FILE_NUMBER  = 1

YEAR_COLUMN_INDEX  = 0   # 0 = 1978, 1 = 1997
YEAR_LABEL         = 1978

# Road domains (1-indexed real domain numbering)
# FIRST_ROAD_DOMAIN  = 9
# LAST_ROAD_DOMAIN   = 90
# Convert to global domain indices
# START_ROAD_INDEX   = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS
# END_ROAD_INDEX     = (LAST_ROAD_DOMAIN  - 1) + NUM_BUFFER_DOMAINS

DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE, 'island_offset', 'hindcast_1978_1997',
    f'Island_Dune_Offsets_1978_1997_PADDED_{TOTAL_DOMAINS}.csv'
)

ELEV_MIN_M = -1.0
ELEV_MAX_M =  4.0
SEA_LEVEL  =  0.0
DAM_TO_M   = 10.0

# =============================================================================
# FILE PATHS  (build full list including buffers for offset loading)
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
# LOAD DUNE OFFSETS
# =============================================================================

print("Loading dune offsets...")
dune_offset_all   = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=',')
dune_offset_m     = dune_offset_all[:, YEAR_COLUMN_INDEX]
dune_offset_cells = np.round(dune_offset_m / DAM_TO_M).astype(int)

print(f"  Offsets (m): min={dune_offset_m.min():.1f}, max={dune_offset_m.max():.1f}")

# =============================================================================
# LOAD ELEVATION ARRAYS (all domains, for consistent canvas height)
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

    domain_grids.append(arr * DAM_TO_M)
    domain_offsets.append(int(dune_offset_cells[i]) if i < len(dune_offset_cells) else 0)

# =============================================================================
# BUILD COMPOSITE CANVAS — REAL DOMAINS ONLY
# =============================================================================

# Subset to real domains only
real_grids   = domain_grids[START_REAL_INDEX:END_REAL_INDEX]
real_offsets = domain_offsets[START_REAL_INDEX:END_REAL_INDEX]

max_cs_raw  = max(g.shape[0] for g in real_grids)
max_offset  = max(real_offsets)
canvas_rows = max_offset + max_cs_raw + 5

total_cols  = sum(g.shape[1] for g in real_grids)
canvas      = np.full((canvas_rows, total_cols), np.nan)

col_cursor        = 0
domain_col_starts = []   # column start for each REAL domain (0-indexed within canvas)
cells_per_domain  = []

for i, grid in enumerate(real_grids):
    n_rows, n_cols = grid.shape
    domain_col_starts.append(col_cursor)
    cells_per_domain.append(n_cols)

    grid_lr = np.fliplr(grid)

    origin  = real_offsets[i]
    row_end = origin + n_rows
    if row_end <= canvas_rows:
        canvas[origin:row_end, col_cursor:col_cursor + n_cols] = grid_lr
    else:
        canvas[origin:canvas_rows, col_cursor:col_cursor + n_cols] = grid_lr[:canvas_rows - origin, :]

    col_cursor += n_cols

domain_col_starts = np.array(domain_col_starts)

print(f"Canvas shape (real domains only): {canvas.shape}")

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

# Aspect: keep cross-shore visually readable — not too squashed
n_cs   = canvas.shape[0]
n_al   = canvas.shape[1]
aspect = n_cs / n_al   # natural cell aspect
fig_w  = 20
fig_h  = max(4.5, fig_w * aspect * 1.8)   # slightly exaggerated cross-shore for readability
fig_h  = min(fig_h, 7.5)                  # cap so it fits on a poster

fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')
ax  = fig.add_axes([0.06, 0.18, 0.88, 0.68])   # [left, bottom, width, height]

# ---- Colormap ----
# matplotlib 'terrain': deep blue (below MSL) → teal → green → tan/sand → white peaks
# Classic topographic map style; high end stays light so dune crests are visible
cmap = plt.cm.terrain.copy()
cmap.set_bad(color='#b0cfe8')   # light blue for NaN/ocean cells

norm = TwoSlopeNorm(vmin=ELEV_MIN_M, vcenter=SEA_LEVEL, vmax=ELEV_MAX_M)

im = ax.pcolormesh(canvas, cmap=cmap, norm=norm, shading='auto', rasterized=True)

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
tick_real_indices = list(range(0, NUM_REAL_DOMAINS, 5))   # 0, 5, 10, ... within real_grids
tick_cols  = [domain_col_starts[i] + cells_per_domain[i] // 2 for i in tick_real_indices]
tick_labels = [str(i + 1) for i in tick_real_indices]     # 1-indexed display

ax.set_xticks(tick_cols)
ax.set_xticklabels(tick_labels, fontsize=9)
ax.set_xlabel('Domain (S → N,  Cape Hatteras to Rodanthe)', fontsize=12, labelpad=8)

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
ax.set_title(
    f'Hatteras Island — CASCADE Initialization  |  {YEAR_LABEL}',
    fontsize=14, fontweight='bold', color='#1a1a2e', pad=12,
)

# ---- Subtle domain grid lines every 10 domains ----
for i in range(0, NUM_REAL_DOMAINS, 10):
    x = domain_col_starts[i] - 0.5
    ax.axvline(x, color='#aaaaaa', lw=0.4, alpha=0.5, zorder=2)

# =============================================================================
# SAVE
# =============================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)
out_path = os.path.join(OUTPUT_DIR, OUTPUT_NAME)
fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
print(f"\n✓ Saved to:\n  {out_path}")
plt.close(fig)
