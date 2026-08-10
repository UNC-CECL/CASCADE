# ==============================================================================
# Hatteras CASCADE — Dune & Topography Handoff Extractor
#
# PURPOSE
# -------
# Extracts the final-timestep dune heights and interior topography from a
# completed CASCADE calibration run (.npz) and writes them as .npy groin_init files
# for the next period's run (e.g. 2004–2024 validation).
#
# The comparison files match the naming convention and array shapes expected by
# HAT_hindcast_1984_2024_updated.py (Section 8), so no changes to the run
# script are needed — just update TOPO_DUNE_INIT_YEAR and TOPO_DUNE_SUBFOLDER
# to point to the new post_run folders.
#
# INPUTS
# ------
#   - Calibration run .npz file (CASCADE comparison)
#     Contains cascade object with cascade.barrier3d[i].DuneDomain and
#     cascade.barrier3d[i].DomainTS for all 120 padded domains.
#
# OUTPUTS
# -------
#   domain_{N}_dune_2004.npy       — shape (50,)       dune height above berm (dam)
#   domain_{N}_topography_2004.npy — shape (200, 50)   island interior (dam)
#   Written for GIS domain IDs 1–90 (real domains only; buffers use existing files).
#
# UNITS
#   DuneDomain : decameters above berm  (Barrier3D native; no conversion needed)
#   DomainTS   : decameters             (Barrier3D native; no conversion needed)
#   Sentinel   : -0.3 dam  (= -3.0 m MHW-relative, matching dune_topo_extractor_V1.py)
#
# TO USE IN THE VALIDATION RUN
# ----------------------------
# In HAT_hindcast_1984_2024_updated.py, update Section 3:
#
#   TOPO_DUNE_INIT_YEAR = "2004"
#   TOPO_DUNE_SUBFOLDER = "post_run"
#
# ==============================================================================

import os
import sys
import numpy as np
from pathlib import Path

# =============================================================================
# SECTION 1: CONFIGURATION — edit these before running
# =============================================================================

# Path to the calibration run .npz file
CALIB_NPZ_PATH = (
    r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1984_2004_BASE_Hs2.5\HAT_1984_2004_BASE_Hs2.5.npz"
)

# Final timestep to extract (0-indexed).
# For a 20-year run (1984–2004), the last annual timestep is index 20.
# The GIF script uses range(TMAX_SIM), so t=20 is the post-final-year state.
# Set to None to auto-detect as the last available timestep.
T_FINAL = None   # <- set to an int to override, e.g. T_FINAL = 20

# Output directories (will be created if they don't exist)
TOPO_SAVE_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\topography\post_run"
DUNE_SAVE_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\dunes\post_run"

# Year label used in comparison filenames — must match TOPO_DUNE_INIT_YEAR in run script
INIT_YEAR_LABEL = "2004"

# Domain numbering — must match HAT_hindcast_1984_2024_updated.py Section 1
NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 1    # GIS domain IDs: 1..90

# Expected array dimensions (must match dune_topo_extractor_V1.py)
TOPO_ROWS  = 200   # cross-shore rows in topography array
ALONG_COLS = 50    # alongshore columns in topography array

# Sentinel value for water / missing cells (dam), matching V1 extractor
SENTINEL_DAM = -0.3

# =============================================================================
# SECTION 2: DIAGNOSTIC SETTINGS
# =============================================================================

# Print per-domain summary stats (min/max dune, mean interior) for spot-checking
PRINT_DOMAIN_STATS = True

# Number of domains to print stats for (set to NUM_REAL_DOMAINS to print all)
STATS_PRINT_N = 10

# =============================================================================
# SECTION 3: LOAD CALIBRATION RUN
# =============================================================================

print("=" * 70)
print("Hatteras CASCADE — Dune & Topography Handoff Extractor")
print("=" * 70)
print(f"\nLoading calibration run:")
print(f"  {CALIB_NPZ_PATH}")

if not os.path.exists(CALIB_NPZ_PATH):
    print(f"\n❌ ERROR: File not found: {CALIB_NPZ_PATH}")
    sys.exit(1)

output  = np.load(CALIB_NPZ_PATH, allow_pickle=True)
cascade = output["cascade"][0]
b3d     = cascade.barrier3d

total_domains = len(b3d)
start_real    = NUM_BUFFER_DOMAINS                          # padded index of first real domain
end_real      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS       # padded index one past last real domain

print(f"✓ Loaded. Total padded domains: {total_domains}")
print(f"  Real domain padded indices: {start_real}–{end_real - 1}")
print(f"  GIS domain IDs:             {FIRST_FILE_NUMBER}–{FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1}")

if total_domains != NUM_BUFFER_DOMAINS * 2 + NUM_REAL_DOMAINS:
    print(
        f"\n⚠️  WARNING: Expected {NUM_BUFFER_DOMAINS * 2 + NUM_REAL_DOMAINS} total domains, "
        f"got {total_domains}. Check NUM_BUFFER_DOMAINS / NUM_REAL_DOMAINS."
    )

# =============================================================================
# SECTION 4: DETERMINE FINAL TIMESTEP
# =============================================================================

# DomainTS is a list of arrays, one per timestep.
# Use the first real domain to find how many timesteps were recorded.
sample_domain_ts = b3d[start_real].DomainTS
n_timesteps = len(sample_domain_ts)

if T_FINAL is None:
    T_FINAL = n_timesteps - 1
    print(f"\nAuto-detected T_FINAL = {T_FINAL}  ({n_timesteps} total timesteps in DomainTS)")
else:
    print(f"\nUsing configured T_FINAL = {T_FINAL}  ({n_timesteps} total timesteps in DomainTS)")

if T_FINAL >= n_timesteps:
    print(
        f"❌ ERROR: T_FINAL={T_FINAL} exceeds available timesteps (0–{n_timesteps - 1})."
    )
    sys.exit(1)

# =============================================================================
# SECTION 5: DIAGNOSTIC — DuneDomain axis layout
# =============================================================================
# Actual DuneDomain shape in Barrier3D: (time, alongshore, dune_rows)
# i.e. axis 1 = alongshore cells (50), axis 2 = dune rows (2)
#
# Dune row convention: index 0 = SEAWARD crest (the active dune face that
# controls overwash), index -1 = landward row.
#
# Correct extraction for the seaward crest across all 50 alongshore cells:
#   DuneDomain[T_FINAL, :, 0]   ← all alongshore, seaward row

sample_dune  = b3d[start_real].DuneDomain   # shape: (time, alongshore, dune_rows)
n_along_dune = sample_dune.shape[1]
n_dune_rows  = sample_dune.shape[2]

print(f"\nDuneDomain dimensions (sample domain {start_real}):")
print(f"  Shape: {sample_dune.shape}  →  (time={sample_dune.shape[0]}, "
      f"alongshore={n_along_dune}, dune_rows={n_dune_rows})")
print(f"  At T_FINAL={T_FINAL}:")
print(f"    Col  0 (seaward row,  all alongshore): "
      f"min={sample_dune[T_FINAL, :, 0].min():.4f}  "
      f"max={sample_dune[T_FINAL, :, 0].max():.4f} dam")
if n_dune_rows > 1:
    print(f"    Col -1 (landward row, all alongshore): "
          f"min={sample_dune[T_FINAL, :, -1].min():.4f}  "
          f"max={sample_dune[T_FINAL, :, -1].max():.4f} dam")
print("  → DuneDomain[T_FINAL, :, 0] (seaward crest) will be used.")

# =============================================================================
# SECTION 6: EXTRACT AND SAVE
# =============================================================================

os.makedirs(TOPO_SAVE_PATH, exist_ok=True)
os.makedirs(DUNE_SAVE_PATH, exist_ok=True)

print(f"\nExtracting T_FINAL={T_FINAL} for {NUM_REAL_DOMAINS} real domains...")
print(f"  Topo → {TOPO_SAVE_PATH}")
print(f"  Dune → {DUNE_SAVE_PATH}")
print()

n_written      = 0
n_shape_warn   = 0
n_sentinel_fix = 0

for pad_idx in range(start_real, end_real):

    gis_id   = FIRST_FILE_NUMBER + (pad_idx - start_real)   # 1..90
    domain   = b3d[pad_idx]

    # ------------------------------------------------------------------
    # DUNE: shape (50,) — seaward dune crest, decameters above berm
    # DuneDomain axes: (time, alongshore, dune_rows)
    # [:, 0] = seaward row across all 50 alongshore cells
    # ------------------------------------------------------------------
    dune_raw = domain.DuneDomain[T_FINAL, :, 0]   # (alongshore=50,), seaward crest

    # Allocate fixed-size comparison (ALONG_COLS = 50) with sentinel fill
    dune_dm = np.full(ALONG_COLS, fill_value=SENTINEL_DAM, dtype=float)

    n_copy = min(len(dune_raw), ALONG_COLS)
    dune_dm[:n_copy] = dune_raw[:n_copy]

    # Clamp any negative dune heights to a small positive value (0.01 dam = 0.1 m)
    # so CASCADE doesn't interpret them as subaqueous cells.
    # This mirrors the "dune_h_m < 0 → 0.1 m" guard in dune_topo_extractor_V1.py.
    dune_dm = np.where(dune_dm < 0.0, 0.01, dune_dm)

    # ------------------------------------------------------------------
    # TOPOGRAPHY: shape (200, 50) — island interior, decameters
    # ------------------------------------------------------------------
    topo_raw = domain.DomainTS[T_FINAL]   # shape: (cross_shore, alongshore)

    # Allocate fixed-size comparison with sentinel fill
    topo_dm = np.full((TOPO_ROWS, ALONG_COLS), fill_value=SENTINEL_DAM, dtype=float)

    # DomainTS may not be exactly (200, 50) — copy what fits, warn if different
    raw_rows, raw_cols = topo_raw.shape
    copy_rows = min(raw_rows, TOPO_ROWS)
    copy_cols = min(raw_cols, ALONG_COLS)

    # Note: varying cross-shore rows are expected — CASCADE adjusts island width
    # dynamically as overwash and erosion occur. The alongshore dimension (50)
    # must stay fixed; the cross-shore dimension is allowed to vary.
    # A mismatch in the *alongshore* dimension would be a real problem.
    if raw_cols != ALONG_COLS:
        n_shape_warn += 1
        print(
            f"  ❌  GIS {gis_id} (pad {pad_idx}): DomainTS alongshore cols {raw_cols} "
            f"≠ expected {ALONG_COLS}. This is a real error — check domain config."
        )

    topo_dm[:copy_rows, :copy_cols] = topo_raw[:copy_rows, :copy_cols]

    # ------------------------------------------------------------------
    # SAVE
    # ------------------------------------------------------------------
    topo_filename = f"domain_{gis_id}_topography_{INIT_YEAR_LABEL}.npy"
    dune_filename = f"domain_{gis_id}_dune_{INIT_YEAR_LABEL}.npy"

    np.save(os.path.join(TOPO_SAVE_PATH, topo_filename), topo_dm)
    np.save(os.path.join(DUNE_SAVE_PATH, dune_filename), dune_dm)
    n_written += 1

    # ------------------------------------------------------------------
    # OPTIONAL DIAGNOSTIC STATS (first N domains)
    # ------------------------------------------------------------------
    if PRINT_DOMAIN_STATS and gis_id <= FIRST_FILE_NUMBER + STATS_PRINT_N - 1:
        interior_valid = topo_dm[topo_dm > SENTINEL_DAM]
        print(
            f"  GIS {gis_id:2d} (pad {pad_idx:3d}): "
            f"dune min={dune_dm.min():.4f} max={dune_dm.max():.4f} mean={dune_dm.mean():.4f} dam  |  "
            f"topo valid cells={interior_valid.size}  "
            f"mean={interior_valid.mean():.4f} dam"
            if interior_valid.size > 0 else
            f"  GIS {gis_id:2d} (pad {pad_idx:3d}): "
            f"dune min={dune_dm.min():.4f} max={dune_dm.max():.4f} dam  |  "
            f"topo: ALL SENTINEL (check DomainTS)"
        )

print()
print("=" * 70)
print(f"✓ Done. {n_written} domain pairs written.")
if n_shape_warn > 0:
    print(f"❌  {n_shape_warn} domain(s) had unexpected alongshore column count — see errors above.")
print()
print("Next steps:")
print("  1. Visually spot-check a few comparison files with np.load() and print shapes/values.")
print("  2. In HAT_hindcast_1984_2024_updated.py Section 3, set:")
print(f"       TOPO_DUNE_INIT_YEAR = \"{INIT_YEAR_LABEL}\"")
print(f"       TOPO_DUNE_SUBFOLDER = \"post_run\"")
print("  3. Run the 2004–2024 validation period.")
print("=" * 70)
