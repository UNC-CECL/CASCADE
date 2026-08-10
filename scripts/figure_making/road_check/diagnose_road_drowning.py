"""
diagnose_road_drowning.py
=========================
Visualises the INITIAL cross-shore elevation profile for any set of CASCADE
domains and marks exactly where the road sits and which cells the drowning
check inspects at t=0.

Loads topography from the same .npy files that CASCADE reads at startup —
NOT from the saved NPZ (which reflects end-of-run state) — so the check
cells reflect exactly what RoadwayManager sees on the very first time step.

Elevation units: Barrier3D stores interior topography in decametres (dam)
relative to MHW.  Multiply by 10 to convert to metres MHW.
Cell spacing: 10 m cross-shore, 10 m alongshore.

Usage
-----
Edit the CONFIG block, then run directly in PyCharm.

Placement
---------
    CASCADE/scripts/hatteras_ms/diagnose_road_drowning.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# =============================================================================
# CONFIG — match your run script
# =============================================================================

PROJECT_BASE_DIR   = r"C:\Users\hanna\PycharmProjects\CASCADE"
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
NPZ_PATH = (
    r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
    r"\HAT_1984_2004_base\HAT_1984_2004_base.npz"
)

# Topo/dune file settings — must match run script
TOPO_DUNE_INIT_YEAR = "2009"
TOPO_DUNE_SUBFOLDER = r"2009\2009_v2"

# Domain index constants — must match run script
NUM_BUFFER_DOMAINS = 15
NUM_REAL_DOMAINS   = 90
FIRST_FILE_NUMBER  = 1
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS   # = 105
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS * 2 + NUM_REAL_DOMAINS  # = 120

# Grid spacing
CELL_SIZE_M = 10.0   # metres per cross-shore cell (and alongshore cell)
DAM_TO_M    = 10.0   # 1 decametre = 10 metres

# Drowning check parameters (roadway_manager.py defaults)
DROWN_THRESHOLD_M = 0.0    # metres MHW — cells at or below this are "water"
PERCENT_THRESHOLD = 0.20   # 20% of border cells must be water to trigger drown

# Domains to plot — padded indices
# Drowning domains identified by diagnostic: GIS 12, 14, 15, 50, 51, 86
PLOT_DOMAINS    = [26, 28, 29, 64, 65, 100]
# Healthy domains for comparison
HEALTHY_DOMAINS = [30, 50, 70]

# =============================================================================
# HELPERS
# =============================================================================

def pad_to_gis(pad):
    return FIRST_FILE_NUMBER + (pad - START_REAL_INDEX)

def road_indices(setback_m, road_width_m):
    """Return (road_start, road_end) grid indices.
    Mirrors roadway_manager.py: road_start = int(setback_m / dy) where dy=10 m.
    """
    road_start = int(setback_m / CELL_SIZE_M)
    road_end   = road_start + int(road_width_m / CELL_SIZE_M)
    return road_start, road_end

def load_initial_interior(pad):
    """Load the initial interior topography .npy for a padded domain index.

    Returns array in dam MHW, shape (cross_shore_cells, alongshore_cells).
    Buffer domains use the shared buffer file.
    Real domains use domain_{N}_topography_{year}.npy.
    """
    if pad < START_REAL_INDEX or pad >= END_REAL_INDEX:
        fpath = os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_topography.npy")
    else:
        file_num = FIRST_FILE_NUMBER + (pad - START_REAL_INDEX)
        fpath = os.path.join(
            HATTERAS_DATA_BASE, "topography", TOPO_DUNE_SUBFOLDER,
            f"domain_{file_num}_topography_{TOPO_DUNE_INIT_YEAR}.npy"
        )
    arr = np.load(fpath, allow_pickle=True)
    # Ensure cross-shore is axis 0 (200 cells) not axis 1 (50 cells)
    if arr.ndim == 2 and arr.shape[0] < arr.shape[1]:
        arr = arr.T
    return arr   # units: dam MHW

def water_fraction(interior_dam, row_idx):
    """Fraction of alongshore cells in a cross-shore row at or below MHW.

    interior_dam : array in dam MHW, shape (cross_shore, alongshore)
    row_idx      : cross-shore row to check
    Returns float in [0, 1].
    """
    if 0 <= row_idx < interior_dam.shape[0]:
        row_m = interior_dam[row_idx, :] * DAM_TO_M   # dam → m MHW
        return float(np.mean(row_m <= DROWN_THRESHOLD_M))
    return np.nan

def mean_profile_m(interior_dam):
    """Alongshore-mean cross-shore elevation profile in metres MHW."""
    return np.mean(interior_dam, axis=1) * DAM_TO_M   # dam → m MHW

# =============================================================================
# LOAD NPZ — for road manager state only (setbacks, drown_break)
# =============================================================================

print(f"Loading NPZ for road manager state: {NPZ_PATH}")
if not os.path.exists(NPZ_PATH):
    raise FileNotFoundError(f"NPZ not found:\n  {NPZ_PATH}")

data     = np.load(NPZ_PATH, allow_pickle=True)
cascade  = data["cascade.npy"].item()
roadways = cascade.roadways
print("  Loaded. Road manager setbacks and drown_break flags available.\n")

# =============================================================================
# PLOT
# =============================================================================

all_domains = PLOT_DOMAINS + HEALTHY_DOMAINS
n     = len(all_domains)
ncols = 3
nrows = int(np.ceil(n / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for ax_idx, pad in enumerate(all_domains):
    ax  = axes[ax_idx]
    rm  = roadways[pad] if roadways is not None else None
    gis = pad_to_gis(pad) if START_REAL_INDEX <= pad < END_REAL_INDEX else f"buf{pad}"

    # ── Load initial topography ───────────────────────────────────────────────
    try:
        interior_dam = load_initial_interior(pad)
    except FileNotFoundError as e:
        ax.set_title(f"GIS {gis} (pad {pad})\nFILE NOT FOUND", fontsize=8)
        ax.text(0.5, 0.5, str(e), transform=ax.transAxes,
                ha="center", va="center", fontsize=6, color="red")
        continue

    profile_m = mean_profile_m(interior_dam)           # m MHW, length = n_cross_shore
    n_cells   = len(profile_m)
    x_m       = np.arange(n_cells) * CELL_SIZE_M      # metres from seaward edge

    # ── Road position ─────────────────────────────────────────────────────────
    # Use the first entry of _road_setback_TS for the true initial setback.
    # Fall back to _road_setback if TS is unavailable or all-zero.
    if rm is not None and hasattr(rm, '_road_setback_TS'):
        ts = np.asarray(rm._road_setback_TS)
        nonzero = ts[ts > 0]
        setback_m = float(nonzero[0]) if len(nonzero) > 0 else float(rm._road_setback)
    elif rm is not None and hasattr(rm, '_road_setback'):
        setback_m = float(rm._road_setback)
    else:
        setback_m = 30.0   # fallback

    road_width_m = float(rm._road_width) if rm and hasattr(rm, '_road_width') else 20.0
    drown_break  = int(getattr(rm, '_drown_break', 0)) if rm else 0

    rs, re = road_indices(setback_m, road_width_m)
    rs = min(rs, n_cells - 2)
    re = min(re, n_cells - 1)

    sea_row = rs - 1
    bay_row = re + 1
    sf = water_fraction(interior_dam, sea_row)
    bf = water_fraction(interior_dam, bay_row)

    init_drowned = (sf > PERCENT_THRESHOLD) or (bf > PERCENT_THRESHOLD)
    was_drowned  = drown_break == 1
    suffix       = " ⬅ DROWNED" if was_drowned else " (healthy)"

    # ── Plot profile ──────────────────────────────────────────────────────────
    line_color = "#d32f2f" if was_drowned else "#1976d2"
    ax.plot(x_m, profile_m, color=line_color, lw=1.5, label=f"GIS {gis}{suffix}")
    ax.axhline(DROWN_THRESHOLD_M, color="navy", lw=0.9, ls="--", alpha=0.6,
               label="MHW (0 m)")

    # Road span
    ax.axvspan(rs * CELL_SIZE_M, re * CELL_SIZE_M,
               color="#FFA726", alpha=0.45, label="Road width")

    # Seaward check cell
    if 0 <= sea_row < n_cells:
        ec = "#d32f2f" if sf > PERCENT_THRESHOLD else "#43a047"
        ax.axvspan(sea_row * CELL_SIZE_M, (sea_row + 1) * CELL_SIZE_M,
                   color=ec, alpha=0.35,
                   label=f"Sea check — {sf*100:.0f}% ≤ MHW")

    # Bayward check cell
    if 0 <= bay_row < n_cells:
        ec = "#d32f2f" if bf > PERCENT_THRESHOLD else "#43a047"
        ax.axvspan(bay_row * CELL_SIZE_M, (bay_row + 1) * CELL_SIZE_M,
                   color=ec, alpha=0.35,
                   label=f"Bay check — {bf*100:.0f}% ≤ MHW")

    status_str = "⚠ DROWNED (groin_init topo)" if init_drowned else "✓ OK (groin_init topo)"
    ax.set_title(
        f"GIS {gis} (pad {pad}){suffix}\n"
        f"setback={setback_m:.0f} m  |  road rows {rs}–{re}  |  "
        f"sea={sf*100:.0f}%  bay={bf*100:.0f}%  →  {status_str}",
        fontsize=7.5
    )
    ax.set_ylabel("Elevation (m MHW)", fontsize=8)
    ax.set_xlabel("Cross-shore distance from seaward edge (m)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_ylim(-1.5, 6.0)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=6, loc="upper right")

    # ── Zoom x-axis to road neighbourhood ──────────────────────────────────────────────
    # Show from x=0 to 100 m beyond the bayward check cell so the road
    # and both check cells are clearly visible without the full bay.
    x_zoom = max((bay_row + 3) * CELL_SIZE_M, setback_m + road_width_m + 150)
    x_zoom = min(x_zoom, n_cells * CELL_SIZE_M)
    ax.set_xlim(0, x_zoom)

for ax in axes[n:]:
    ax.set_visible(False)

legend_patches = [
    mpatches.Patch(color="#FFA726", alpha=0.5,  label="Road width"),
    mpatches.Patch(color="#d32f2f", alpha=0.35, label="Check cell > 20% at/below MHW"),
    mpatches.Patch(color="#43a047", alpha=0.35, label="Check cell ≤ 20% at/below MHW"),
    plt.Line2D([0], [0], color="navy",    ls="--", lw=1,  label="MHW (0 m)"),
    plt.Line2D([0], [0], color="#d32f2f", lw=1.5,         label="Drowned domain (drown_break=1)"),
    plt.Line2D([0], [0], color="#1976d2", lw=1.5,         label="Healthy domain"),
]
fig.legend(handles=legend_patches, loc="lower center",
           ncol=3, fontsize=8, bbox_to_anchor=(0.5, -0.02))

fig.suptitle(
    "Cross-shore Road Position Diagnostic — INITIAL TOPOGRAPHY (2009 DEM)\n"
    "Orange = road width  |  Red check cell = drowning trigger  |  Navy dashed = MHW (0 m)",
    fontsize=10, fontweight="bold"
)
plt.tight_layout(rect=[0, 0.05, 1, 1])

out_path = os.path.join(os.path.dirname(NPZ_PATH), "road_drowning_diagnostic_initial.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"✓ Saved: {out_path}")
plt.show()

# =============================================================================
# NUMERICAL SUMMARY
# =============================================================================

print("\nNumerical summary — INITIAL topography:")
print(f"{'GIS':>4}  {'pad':>4}  {'setback_m':>10}  {'road rows':>10}  "
      f"{'sea wet%':>9}  {'bay wet%':>9}  {'groin_init status':>22}  {'run status':>10}")
print("-" * 90)

for pad in all_domains:
    gis = pad_to_gis(pad) if START_REAL_INDEX <= pad < END_REAL_INDEX else "buf"
    rm  = roadways[pad] if roadways is not None else None

    try:
        interior_dam = load_initial_interior(pad)
    except FileNotFoundError:
        print(f"{str(gis):>4}  {pad:>4}  FILE NOT FOUND")
        continue

    if rm and hasattr(rm, '_road_setback_TS'):
        ts = np.asarray(rm._road_setback_TS)
        nonzero = ts[ts > 0]
        setback_m = float(nonzero[0]) if len(nonzero) > 0 else float(rm._road_setback)
    elif rm and hasattr(rm, '_road_setback'):
        setback_m = float(rm._road_setback)
    else:
        setback_m = float('nan')

    road_width_m = float(rm._road_width) if rm and hasattr(rm, '_road_width') else 20.0
    drown_break  = int(getattr(rm, '_drown_break', 0)) if rm else 0

    if not np.isnan(setback_m):
        rs, re = road_indices(setback_m, road_width_m)
        rs = min(rs, interior_dam.shape[0] - 2)
        re = min(re, interior_dam.shape[0] - 1)
    else:
        rs, re = 0, 0

    sf = water_fraction(interior_dam, rs - 1)
    bf = water_fraction(interior_dam, re + 1)

    init_ok  = (sf <= PERCENT_THRESHOLD) and (bf <= PERCENT_THRESHOLD)
    run_ok   = drown_break != 1
    init_str = "✓ OK" if init_ok  else "⚠ DROWNED (groin_init topo)"
    run_str  = "✓ OK" if run_ok   else "⚠ drowned"

    print(f"{str(gis):>4}  {pad:>4}  {setback_m:>10.1f}  "
          f"{str(rs)+'-'+str(re):>10}  {sf*100:>8.0f}%  {bf*100:>8.0f}%  "
          f"{init_str:>22}  {run_str:>10}")
