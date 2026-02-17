"""
plot_initialization.py
----------------------
Static plan-view plot of CASCADE initialization (t=0) for Hatteras Island,
showing all 120 domains (real + buffers) with their shoreline offsets applied.

Loads data the same way HAT_1978_1997_status.py does — no CASCADE run needed.
Output: a single PNG saved to OUTPUT_DIR.

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
# CONFIGURATION — mirrors HAT_1978_1997_status.py exactly
# =============================================================================

PROJECT_BASE_DIR   = r'C:\Users\hanna\PycharmProjects\CASCADE'
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_DIR         = os.path.join(PROJECT_BASE_DIR, 'scripts', 'analyze_outputs', 'plots')
OUTPUT_NAME        = 'initialization_overview_1978_v5.png'

# Domain layout
NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120

START_REAL_INDEX   = NUM_BUFFER_DOMAINS          # 15
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS  # 105
FIRST_FILE_NUMBER  = 1

# Which year column to use from the offset CSV (0 = 1978, 1 = 1997)
YEAR_COLUMN_INDEX  = 0
YEAR_LABEL         = 1978

# Road extent (domain numbers, 1-indexed within real domains)
FIRST_ROAD_DOMAIN  = 9
LAST_ROAD_DOMAIN   = 90
START_ROAD_INDEX   = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS   # global index
END_ROAD_INDEX     = (LAST_ROAD_DOMAIN  - 1) + NUM_BUFFER_DOMAINS   # global index (inclusive)

# Dune offset file
DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE, 'island_offset', 'hindcast_1978_1997_constant',
    f'Island_Dune_Offsets_1978_1997_constant_PADDED_{TOTAL_DOMAINS}.csv'
)

# Elevation display range (meters)
ELEV_MIN_M = -1.0
ELEV_MAX_M =  4.0
SEA_LEVEL  =  0.0
DAM_TO_M   = 10.0   # CASCADE internal: decameters → meters

# =============================================================================
# BUILD FILE PATH LISTS (identical logic to init script)
# =============================================================================

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS      = []

# South buffer
for _ in range(START_REAL_INDEX):
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy'))
    DUNE_FILE_PATHS.append(     os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy'))

# Real domains
for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'topography', '2009',
                                              f'domain_{file_num}_topography_2009.npy'))
    DUNE_FILE_PATHS.append(     os.path.join(HATTERAS_DATA_BASE, 'dunes',      '2009',
                                              f'domain_{file_num}_dune_2009.npy'))

# North buffer
for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy'))
    DUNE_FILE_PATHS.append(     os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy'))

# =============================================================================
# LOAD DUNE OFFSETS
# =============================================================================

print(f"Loading dune offsets from:\n  {DUNE_OFFSET_FILE}")
dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=',')
dune_offset_dam = dune_offset_all[:, YEAR_COLUMN_INDEX] / 10.0   # m → dam
dune_offset_cells = np.round(dune_offset_dam).astype(int)         # integer cell shifts

print(f"  Offsets (dam): min={dune_offset_dam.min():.1f}, max={dune_offset_dam.max():.1f}")
print(f"  Offset range in cells: {dune_offset_cells.min()} to {dune_offset_cells.max()}")

# =============================================================================
# LOAD ELEVATION ARRAYS AND APPLY OFFSETS
# =============================================================================

print(f"\nLoading {TOTAL_DOMAINS} domain elevation arrays...")

domain_grids   = []   # raw interior grids (cross-shore × 1 or few cols each)
domain_offsets = []

for i, elev_path in enumerate(ELEVATION_FILE_PATHS):
    try:
        arr = np.load(elev_path).astype(float)   # shape varies: (cross_shore, alongshore_cols)
    except FileNotFoundError:
        print(f"  ⚠  Domain {i}: file not found — {os.path.basename(elev_path)}")
        arr = np.zeros((30, 1))   # fallback blank column
    except Exception as e:
        print(f"  ⚠  Domain {i}: load error — {e}")
        arr = np.zeros((30, 1))

    # Ensure 2D
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)

    domain_grids.append(arr * DAM_TO_M)   # convert dam → m
    domain_offsets.append(int(dune_offset_cells[i]) if i < len(dune_offset_cells) else 0)

print(f"  Shapes: cross-shore rows in [" +
      f"{min(g.shape[0] for g in domain_grids)}, "
      f"{max(g.shape[0] for g in domain_grids)}]")

# =============================================================================
# BUILD COMPOSITE PLAN-VIEW ARRAY WITH OFFSETS
# =============================================================================

# Each domain occupies its own alongshore column(s).
# The offset shifts the domain oceanward/shoreward in cross-shore.
# We place each domain in a canvas large enough to fit all offsets.

max_cs_raw = max(g.shape[0] for g in domain_grids)
max_offset  = max(abs(o) for o in domain_offsets)
canvas_rows = max_cs_raw + max_offset + 5   # extra padding

total_cols  = sum(g.shape[1] for g in domain_grids)
canvas      = np.full((canvas_rows, total_cols), np.nan)

col_cursor  = 0
domain_col_starts = []

for i, (grid, offset) in enumerate(zip(domain_grids, domain_offsets)):
    n_rows, n_cols = grid.shape
    domain_col_starts.append(col_cursor)

    # Positive offset = domain shifted in positive cross-shore direction (landward in Barrier3D)
    # We place ocean side at row 0 → offset shifts start row
    row_start = max(0, offset)
    row_end   = row_start + n_rows

    if row_end <= canvas_rows:
        canvas[row_start:row_end, col_cursor:col_cursor + n_cols] = grid
    else:
        # Clip if somehow oversized
        clip = canvas_rows - row_start
        canvas[row_start:canvas_rows, col_cursor:col_cursor + n_cols] = grid[:clip, :]

    col_cursor += n_cols

domain_col_starts = np.array(domain_col_starts)
cells_per_domain  = [g.shape[1] for g in domain_grids]

print(f"\nComposite canvas shape: {canvas.shape}  (cross-shore × alongshore cells)")

# =============================================================================
# PLOT
# =============================================================================

fig_width  = min(30, max(14, TOTAL_DOMAINS * 0.22))
fig_height = 7
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

cmap = plt.cm.terrain.copy()
cmap.set_bad(color='#d0d0d0')   # NaN (empty/gap) = light gray

norm = TwoSlopeNorm(vmin=ELEV_MIN_M, vcenter=SEA_LEVEL, vmax=ELEV_MAX_M)

im = ax.imshow(
    canvas,
    aspect='auto',
    cmap=cmap,
    norm=norm,
    origin='upper',          # row 0 = ocean side
    interpolation='nearest',
)

cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cbar.set_label('Elevation (m MSL)', fontsize=11)

# --- Buffer domain shading ---
buf_color = 'steelblue'

left_end_col  = domain_col_starts[START_REAL_INDEX]
right_start_col = domain_col_starts[END_REAL_INDEX]

ax.axvspan(-0.5,                  left_end_col - 0.5,    color=buf_color, alpha=0.15, zorder=3)
ax.axvspan(right_start_col - 0.5, canvas.shape[1] - 0.5, color=buf_color, alpha=0.15, zorder=3)
ax.axvline(left_end_col - 0.5,    color=buf_color, lw=2.0, ls='--', zorder=4)
ax.axvline(right_start_col - 0.5, color=buf_color, lw=2.0, ls='--', zorder=4)

# --- Buffer domain labels ---
ax.text(left_end_col / 2, 2, f'S Buffer\n(0–{START_REAL_INDEX - 1})',
        ha='center', va='top', fontsize=8, color='navy', style='italic')
ax.text((right_start_col + canvas.shape[1]) / 2, 2, f'N Buffer\n({END_REAL_INDEX}–{TOTAL_DOMAINS - 1})',
        ha='center', va='top', fontsize=8, color='navy', style='italic')

# --- Road extent ---
road_col_start = domain_col_starts[START_ROAD_INDEX]
road_col_end   = domain_col_starts[END_ROAD_INDEX] + cells_per_domain[END_ROAD_INDEX] - 1
ax.axvline(road_col_start - 0.5, color='saddlebrown', lw=1.5, ls=':', zorder=4)
ax.axvline(road_col_end   + 0.5, color='saddlebrown', lw=1.5, ls=':', zorder=4)

# --- Domain index ticks every 5 domains ---
tick_every = 5
tick_domains = list(range(0, TOTAL_DOMAINS, tick_every))
tick_cols    = [domain_col_starts[i] + cells_per_domain[i] // 2 for i in tick_domains]
ax.set_xticks(tick_cols)
ax.set_xticklabels([str(i) for i in tick_domains], fontsize=7, rotation=45)
ax.set_xlabel(f'Domain index  (alongshore, S → N)  —  real domains {START_REAL_INDEX}–{END_REAL_INDEX - 1}',
              fontsize=11)

# --- Offset profile annotation (thin line showing offset curve along top) ---
offset_y = [domain_offsets[i] for i in range(TOTAL_DOMAINS)]
offset_x = [domain_col_starts[i] + cells_per_domain[i] // 2 for i in range(TOTAL_DOMAINS)]
ax_twin = ax.twinx()
ax_twin.plot(offset_x, [o * DAM_TO_M for o in offset_y], color='white', lw=1.2,
             alpha=0.85, label='Shoreline offset (m)', zorder=5)
ax_twin.plot(offset_x, [o * DAM_TO_M for o in offset_y], color='black', lw=0.5,
             alpha=0.5, zorder=5)
ax_twin.set_ylabel('Shoreline offset (m)', fontsize=9, color='white')
ax_twin.tick_params(axis='y', colors='gray', labelsize=8)
ax_twin.set_ylim(bottom=0)

# --- Y-axis cross-shore labels ---
n_rows = canvas.shape[0]
ax.set_yticks([0, n_rows // 3, 2 * n_rows // 3, n_rows - 1])
ax.set_yticklabels(['Ocean\n(row 0)', 'Dune front', 'Interior', 'Back-barrier'], fontsize=8)

# --- Legend ---
legend_handles = [
    mpatches.Patch(facecolor=buf_color, alpha=0.35,
                   label=f'Buffer domains (±{NUM_BUFFER_DOMAINS})'),
    mpatches.Patch(facecolor='saddlebrown', alpha=0.5,
                   label=f'Road extent (domains {FIRST_ROAD_DOMAIN}–{LAST_ROAD_DOMAIN})'),
    plt.Line2D([0], [0], color='white', lw=2, label='Shoreline offset profile'),
]
ax.legend(handles=legend_handles, loc='lower right', fontsize=9, framealpha=0.75)

ax.set_title(
    f'CASCADE Initialization — Hatteras Island  |  {YEAR_LABEL}  |  '
    f'{TOTAL_DOMAINS} total domains  ({NUM_BUFFER_DOMAINS} buf + {NUM_REAL_DOMAINS} real + {NUM_BUFFER_DOMAINS} buf)\n'
    f'Elevation (m MSL), offsets applied, ocean at top',
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
