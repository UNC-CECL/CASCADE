# ==============================================================================
# Hatteras CASCADE Dune & Topography Extractor — Benton-Parity (200-row version)
# INPUT  : numpy array shape = (alongshore_rows, cross_shore_cols), OCEAN on RIGHT
# OUTPUT : topography (200 x 50) and dune (50), in decameters (dam)
# Notes  : Works for any cross-shore length; writes up to 200 rows inland.
# ==============================================================================

import os
from pathlib import Path
import numpy as np

# --- USER PATHS --------------------------------------------------------------
LOAD_PATH       = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\elevations\2009_pea_hatteras"   # folder with *resampled*.npy
TOPO_SAVE_PATH  = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\topography\2009"
DUNE_SAVE_PATH  = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\dunes\2009"

os.makedirs(TOPO_SAVE_PATH, exist_ok=True)
os.makedirs(DUNE_SAVE_PATH, exist_ok=True)

# --- CONSTANTS (Benton’s) ----------------------------------------------------
MHW_M             = 0.26   # meters (NAVD88)
BERM_ELEV_NAVD_M  = 1.70   # meters (NAVD88)
BEACH_START_THR_M = 0.50   # meters (MHW-relative), strict '>' comparison
DUNE_WINDOW_PX    = 8      # pixels (not meters)
SENTINEL_WATER_M  = -3.0   # meters (MHW-relative)
TOPO_ROWS         = 200    # <-- updated: write 200 inland rows
ALONG_COLS        = 50     # still 50 alongshore profiles
AUTO_ORIENT       = True   # auto-detect if ocean is on right and flip if needed

def process_domain_file(in_path: Path, topo_out_dir: Path, dune_out_dir: Path) -> None:
    arr = np.load(in_path).astype(float, copy=False)
    if arr.ndim != 2:
        print(f"[skip] {in_path.name}: expected 2D array, got {arr.ndim}D")
        return

    # Subtract MHW and clamp <-1 m to -3.0 m (MHW-relative) BEFORE any slicing
    arr = arr - MHW_M
    arr[arr < -1.0] = SENTINEL_WATER_M

    n_along, n_cross = arr.shape

    # Sanity messages (informational)
    if n_along < ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} < {ALONG_COLS}; trailing output cols remain sentinel.")
    elif n_along > ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} > {ALONG_COLS}; only first {ALONG_COLS} profiles used.")

    # Auto-detect which edge is ocean (lower elevations)
    ocean_on_right = True
    if AUTO_ORIENT:
        edge = max(1, min(5, n_cross // 20))
        left_q  = np.nanpercentile(arr[:, :edge],  25)
        right_q = np.nanpercentile(arr[:, -edge:], 25)
        ocean_on_right = right_q <= left_q
        print(f"[check] {in_path.name}: edge p25 left={left_q:.2f}, right={right_q:.2f} -> ocean_on_right={ocean_on_right}")

    # Allocate fixed-size outputs
    topo_m = np.full((TOPO_ROWS, ALONG_COLS), fill_value=SENTINEL_WATER_M, dtype=float)
    dune_m = np.full((ALONG_COLS,),             fill_value=SENTINEL_WATER_M, dtype=float)

    # Process up to 50 alongshore profiles
    n_cols_to_fill = min(ALONG_COLS, n_along)
    for i in range(n_cols_to_fill):
        prof_lr = arr[i, :]                              # as stored (land ... ocean-right)
        prof = np.flip(prof_lr) if ocean_on_right else prof_lr  # ensure index 0 == ocean

        # 1) First index where z > 0.5 m (MHW-relative)
        idx = np.where(prof > BEACH_START_THR_M)[0]
        if idx.size == 0:
            continue
        start_beach = int(idx[0])

        # 2) 8-pixel window landward of that point
        end_beach = min(start_beach + DUNE_WINDOW_PX, prof.size)
        if end_beach <= start_beach:
            continue

        # 3) Dune elevation = max in the window
        window = prof[start_beach:end_beach]
        if window.size == 0:
            continue
        dune_elev = float(np.max(window))

        # 4) Dune location = first index in ENTIRE profile equal to dune_elev
        matches = np.where(prof == dune_elev)[0]
        if matches.size == 0:
            continue
        dune_loc = int(matches[0])

        # 5) Island interior starts immediately landward of the dune
        start_island = dune_loc + 1
        use_elev = prof[start_island:-1] if start_island < (prof.size - 1) else np.array([], dtype=float)

        # 6) Write interior into rows (top-down); truncate/pad to TOPO_ROWS (200)
        if use_elev.size > 0:
            rows_to_copy = min(TOPO_ROWS, use_elev.size)
            topo_m[0:rows_to_copy, i] = use_elev[:rows_to_copy]
            # remaining rows stay as sentinel

        # 7) Dune height above berm, min 0.1 m
        dune_h_m = dune_elev - (BERM_ELEV_NAVD_M - MHW_M)  # both MHW-relative
        if dune_h_m < 0.0:
            dune_h_m = 0.1
        dune_m[i] = dune_h_m

    # Convert to decameters (× 0.1), keeping sentinel at -0.3 dam
    topo_dm = topo_m * 0.1
    dune_dm = dune_m * 0.1

    # Save
    stem = in_path.stem  # e.g., "domain_7_resampled"
    topo_out = Path(topo_out_dir) / f"{stem}_topography_2009.npy"
    dune_out = Path(dune_out_dir) / f"{stem}_dune_2009.npy"
    np.save(topo_out, topo_dm)
    np.save(dune_out, dune_dm)
    print(f"[ok] wrote {topo_out.name} (shape {topo_dm.shape}), {dune_out.name} (len {dune_dm.size})")

def main():
    load_dir = Path(LOAD_PATH)
    topo_dir = Path(TOPO_SAVE_PATH)
    dune_dir = Path(DUNE_SAVE_PATH)
    topo_dir.mkdir(parents=True, exist_ok=True)
    dune_dir.mkdir(parents=True, exist_ok=True)

    # Only process resampled domain arrays
    names = sorted([n for n in os.listdir(load_dir) if n.endswith(".npy") and "resampled" in n])
    print(f"[info] Found {len(names)} resampled domain file(s) in {load_dir}")

    for name in names:
        process_domain_file(load_dir / name, topo_dir, dune_dir)

if __name__ == "__main__":
    main()
