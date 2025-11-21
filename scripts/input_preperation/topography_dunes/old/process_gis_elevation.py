# ==============================================================================
# Hatteras CASCADE Dune & Topography Extractor — Benton-Parity Version
# ------------------------------------------------------------------------------
# Matches Benton Franklin's logic exactly, adapted for your input orientation:
#   - INPUT:  rows = alongshore, cols = cross-shore, OCEAN on the RIGHT
#   - OUTPUT: Topography (100 x 50), Dune vector (50), in decameters (dam)
#
# Benton’s behaviors preserved:
#   * Subtract MHW (0.26 m), then clamp any value < -1 m to -3.0 m (MHW-rel)
#   * For each alongshore profile:
#       - Flip so OCEAN is at index 0
#       - Find first index where z > 0.5 m  (MHW-relative)
#       - Dune elevation = max over next 8 pixels
#       - Dune location = first index in the WHOLE flipped profile equal to that max
#         (Yes, this can jump outside the 8-pixel window — we keep this to be exact.)
#       - Interior = prof[(dune_loc + 1) : -1]
#       - Write interior down rows of a 100x50 matrix (unused rows stay -3.0 m)
#       - Dune height above berm = max(dune_elev - 1.7 m, 0.1 m)
#   * Convert both outputs to decameters (× 0.1): sentinel becomes -0.3 dam,
#     min dune becomes 0.01 dam
#
# Notes:
#   * This script assumes there are at least 50 alongshore profiles per domain.
#     If more, extra profiles are ignored; if fewer, the trailing columns remain -3.
# ==============================================================================

import os
import copy
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --- load one test domain ---
f = Path(r"/data/hatteras_init/elevations/2009_pea_hatteras/domain_7_resampled.npy")
arr = np.load(f)

print("Array shape:", arr.shape)  # rows = alongshore, cols = cross-shore (ocean on RIGHT)

# 1) Quick image check
plt.figure(figsize=(7,5))
plt.imshow(arr, cmap="terrain", origin="upper", aspect="auto")
plt.title("Raw elevations (NAVD88); ocean should be on the RIGHT")
plt.colorbar(label="Elevation (m NAVD88)")
plt.show()

# 2) Pick a few alongshore rows, flip to ocean→land
for i in [0, arr.shape[0]//2, arr.shape[0]-1]:
    prof_lr = arr[i, :]           # land→ocean
    prof = np.flip(prof_lr)       # ocean→land
    plt.figure(figsize=(7,3))
    plt.plot(prof, lw=1)
    plt.axhline(0.26, ls='--', label='MHW (0.26 m NAVD88)')
    plt.title(f"Row {i}: flipped profile (ocean index 0)")
    plt.legend()
    plt.tight_layout()
    plt.show()

# 3) Percentile check: ocean side should be lower
edge = 5
left_q  = np.nanpercentile(arr[:, :edge],  25)
right_q = np.nanpercentile(arr[:, -edge:], 25)
print(f"25th percentile near edges: left={left_q:.2f}, right={right_q:.2f} (right should be lower if ocean is on RIGHT)")



# --- USER PATHS --------------------------------------------------------------
LOAD_PATH = r"/data/hatteras_init/elevations/2009_pea_hatteras"
TOPO_SAVE_PATH = r"/data/hatteras_init/topography_dunes/2009"
DUNE_SAVE_PATH = r"/data/hatteras_init/dunes/2009"

os.makedirs(TOPO_SAVE_PATH, exist_ok=True)
os.makedirs(DUNE_SAVE_PATH, exist_ok=True)

# --- CONSTANTS (Benton’s) ----------------------------------------------------
MHW_M             = 0.26   # meters (NAVD88)
BERM_ELEV_NAVD_M  = 1.70   # meters (NAVD88)
BEACH_START_THR_M = 0.50   # meters (MHW-relative), strict '>' comparison
DUNE_WINDOW_PX    = 8      # pixels (not meters)
SENTINEL_WATER_M  = -3.0   # meters (MHW-relative)
TOPO_ROWS         = 100    # fixed
ALONG_COLS        = 50     # fixed

def process_domain_file(in_path: Path, topo_out_dir: Path, dune_out_dir: Path) -> None:
    arr = np.load(in_path).astype(float, copy=False)

    # Subtract MHW and clamp <-1 m to -3.0 m (MHW-relative), BEFORE any slicing
    arr = arr - MHW_M
    arr[arr < -1.0] = SENTINEL_WATER_M

    # Input is (alongshore, cross_shore) with OCEAN on RIGHT.
    n_along, n_cross = arr.shape

    # --- NEW: quick sanity checks -------------------------------------------
    if n_along < ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} < {ALONG_COLS}; "
              f"trailing output columns will remain sentinel (-0.3 dam).")
    elif n_along > ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} > {ALONG_COLS}; "
              f"only the first {ALONG_COLS} profiles will be processed.")

    edge = min(5, max(1, n_cross // 20))
    left_q  = np.nanpercentile(arr[:, :edge],  25)
    right_q = np.nanpercentile(arr[:, -edge:], 25)
    print(f"[check] {in_path.name}: 25th pct left={left_q:.2f} m, right={right_q:.2f} m "
          f"(right should be lower if ocean is on RIGHT)")

    # Allocate fixed-size outputs
    topo_m = np.full((TOPO_ROWS, ALONG_COLS), fill_value=SENTINEL_WATER_M, dtype=float)
    dune_m = np.full((ALONG_COLS,), fill_value=SENTINEL_WATER_M, dtype=float)
    ...

    # Process up to 50 alongshore profiles to match Benton’s fixed output width
    n_cols_to_fill = min(ALONG_COLS, n_along)
    for i in range(n_cols_to_fill):
        # Cross-shore profile for alongshore row i (land ... ocean)
        prof_lr = arr[i, :]            # land -> ocean (RIGHT)
        prof = np.flip(prof_lr)        # OCEAN index 0  ... LAND end

        # 1) Find beginning of beach: first index where z > 0.5 m (MHW-relative)
        #    (strict '>', not '>=')
        idx = np.where(prof > BEACH_START_THR_M)[0]
        if idx.size == 0:
            continue
        start_beach = int(idx[0])

        # 2) Define 8-pixel window landward of that point
        end_beach = start_beach + DUNE_WINDOW_PX
        end_beach = min(end_beach, prof.size)
        if end_beach <= start_beach:
            continue

        # 3) Dune elevation = max in the window
        window = prof[start_beach:end_beach]
        if window.size == 0:
            continue
        dune_elev = float(np.max(window))

        # 4) Dune location = FIRST index in the ENTIRE profile equal to dune_elev
        #    (This replicates Benton exactly, even though it can select outside the window.)
        matches = np.where(prof == dune_elev)[0]
        if matches.size == 0:
            continue
        dune_loc = int(matches[0])

        # 5) Island interior starts immediately landward of the dune
        start_island = dune_loc + 1
        #   Note: Benton slices to -1 (drops the very last element)
        use_elev = prof[start_island:-1] if start_island < (prof.size - 1) else np.array([], dtype=float)

        # 6) Write interior into the column of topo matrix (top-down), meter units still
        if use_elev.size > 0:
            rows_to_copy = min(TOPO_ROWS, use_elev.size)
            topo_m[0:rows_to_copy, i] = use_elev[:rows_to_copy]

        # 7) Dune height above berm, with minimum 0.1 m
        dune_h_m = dune_elev - (BERM_ELEV_NAVD_M - MHW_M)  # both MHW-relative
        if dune_h_m < 0.0:
            dune_h_m = 0.1
        dune_m[i] = dune_h_m

    # Convert to decameters (×0.1), like Benton
    topo_dm = topo_m * 0.1
    dune_dm = dune_m * 0.1

    # Save with Benton’s naming convention
    stem = in_path.stem  # e.g., "domain_1_resampled"
    topo_out = Path(topo_out_dir) / f"{stem}_topography_2009.npy"
    dune_out = Path(dune_out_dir) / f"{stem}_dune_2009.npy"

    np.save(topo_out, topo_dm)
    np.save(dune_out, dune_dm)

def main():
    load_dir = Path(LOAD_PATH)
    topo_dir = Path(TOPO_SAVE_PATH)
    dune_dir = Path(DUNE_SAVE_PATH)
    topo_dir.mkdir(parents=True, exist_ok=True)
    dune_dir.mkdir(parents=True, exist_ok=True)

    npy_files = [f for f in os.listdir(load_dir) if f.endswith(".npy")]
    npy_files.sort()
    print(f"[info] Found {len(npy_files)} domain files in {load_dir}")

    for name in npy_files:
        in_path = load_dir / name
        process_domain_file(in_path, topo_dir, dune_dir)
        print(f"Processed: {name}")

if __name__ == "__main__":
    main()
