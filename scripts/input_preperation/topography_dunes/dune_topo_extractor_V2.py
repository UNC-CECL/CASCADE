# ==============================================================================
# Hatteras CASCADE Dune & Topography Extractor (TRULY FIXED VERSION)
# INPUT  : numpy array shape = (alongshore_rows, cross_shore_cols), OCEAN on RIGHT
# OUTPUT : topography_dunes (200 x 50) and dune (50), in decameters (dam)
#
# KEY FIX: Uses LOCAL beach/berm elevation measured from each profile,
#          not a fixed value. Berm = minimum elevation between beach start and dune.
# ==============================================================================

import os
from pathlib import Path
import numpy as np

# --- USER PATHS --------------------------------------------------------------
LOAD_PATH = r"/data/hatteras_init/elevations/2009_pea_hatteras"
TOPO_SAVE_PATH = r"/data/hatteras_init/topography/2009_FIXED"
DUNE_SAVE_PATH = r"/data/hatteras_init/dunes/2009_FIXED"

os.makedirs(TOPO_SAVE_PATH, exist_ok=True)
os.makedirs(DUNE_SAVE_PATH, exist_ok=True)

# --- CONSTANTS ---------------------------------------------------------------
MHW_M = 0.26  # meters (NAVD88)
BEACH_START_THR_M = 0.50  # meters (MHW-relative), strict '>' comparison
DUNE_WINDOW_PX = 8  # pixels - search window for dune peak
BERM_WINDOW_PX = 5  # pixels - window to find berm before dune
SENTINEL_WATER_M = -3.0  # meters (MHW-relative)
TOPO_ROWS = 200  # number of inland rows to write
ALONG_COLS = 50  # number of alongshore profiles
AUTO_ORIENT = True  # auto-detect if ocean is on right

# NEW: Realistic dune height constraints
MIN_DUNE_HEIGHT_M = 0.1  # Minimum dune height (m)
MAX_DUNE_HEIGHT_M = 4.0  # Maximum realistic dune height (m) - reduced from 6m
EXPECTED_MEAN_HEIGHT_M = 2.0  # Expected typical dune height


def process_domain_file(in_path: Path, topo_out_dir: Path, dune_out_dir: Path) -> None:
    """Process a single domain elevation array and write topography and dune outputs."""
    arr = np.load(in_path).astype(float, copy=False)
    if arr.ndim != 2:
        print(f"[skip] {in_path.name}: expected 2D array, got {arr.ndim}D")
        return

    # Subtract MHW and clamp <-1 m to -3.0 m (MHW-relative)
    arr = arr - MHW_M
    arr[arr < -1.0] = SENTINEL_WATER_M

    n_along, n_cross = arr.shape

    # Sanity messages
    if n_along < ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} < {ALONG_COLS}")
    elif n_along > ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} > {ALONG_COLS}")

    # Auto-detect which edge is ocean
    ocean_on_right = True
    if AUTO_ORIENT:
        edge = max(1, min(5, n_cross // 20))
        left_q = np.nanpercentile(arr[:, :edge], 25)
        right_q = np.nanpercentile(arr[:, -edge:], 25)
        ocean_on_right = right_q <= left_q

    # Allocate fixed-size outputs
    topo_m = np.full((TOPO_ROWS, ALONG_COLS), fill_value=SENTINEL_WATER_M, dtype=float)
    dune_m = np.full((ALONG_COLS,), fill_value=SENTINEL_WATER_M, dtype=float)

    # Diagnostics
    dune_heights_found = []
    berm_elevations_found = []
    profiles_processed = 0
    profiles_with_tall_dunes = 0

    # Process up to ALONG_COLS alongshore profiles
    n_cols_to_fill = min(ALONG_COLS, n_along)
    for i in range(n_cols_to_fill):
        # Profile as stored
        prof_lr = arr[i, :]
        # Ensure index 0 == ocean
        prof = np.flip(prof_lr) if ocean_on_right else prof_lr

        # 1) First index where z > 0.5 m (MHW-relative) = beach start
        idx = np.where(prof > BEACH_START_THR_M)[0]
        if idx.size == 0:
            continue
        start_beach = int(idx[0])

        # 2) Define search window for dune peak
        end_dune_search = min(start_beach + DUNE_WINDOW_PX, prof.size)
        if end_dune_search <= start_beach:
            continue

        # 3) Find dune peak elevation in search window
        dune_search_window = prof[start_beach:end_dune_search]
        if dune_search_window.size == 0:
            continue
        dune_elev = float(np.max(dune_search_window))

        # 4) Find LOCAL berm elevation
        # The berm is the relatively flat area before the dune
        # Use a window before the dune peak to find minimum elevation
        # This represents the beach/berm surface
        berm_search_end = min(start_beach + BERM_WINDOW_PX, end_dune_search)
        berm_window = prof[start_beach:berm_search_end]

        # Berm elevation = minimum in the window (the flat beach before dune)
        # Could also use median or percentile, but min is safest
        beach_elev_local = float(np.min(berm_window))

        # Alternative: use median if you want to avoid outliers
        # beach_elev_local = float(np.median(berm_window))

        # 5) Calculate dune HEIGHT above LOCAL berm
        dune_h_m = dune_elev - beach_elev_local

        # Quality control checks
        if dune_h_m < MIN_DUNE_HEIGHT_M:
            dune_h_m = MIN_DUNE_HEIGHT_M
        elif dune_h_m > MAX_DUNE_HEIGHT_M:
            print(f"[warn] {in_path.name} profile {i}: dune height {dune_h_m:.2f}m "
                  f"exceeds max ({MAX_DUNE_HEIGHT_M}m). "
                  f"(dune_elev={dune_elev:.2f}m, berm={beach_elev_local:.2f}m) Capping.")
            dune_h_m = MAX_DUNE_HEIGHT_M
            profiles_with_tall_dunes += 1

        dune_m[i] = dune_h_m
        dune_heights_found.append(dune_h_m)
        berm_elevations_found.append(beach_elev_local)
        profiles_processed += 1

        # 6) Find dune location for interior extraction
        matches = np.where(prof == dune_elev)[0]
        if matches.size == 0:
            continue
        dune_loc = int(matches[0])

        # 7) Island interior starts immediately landward of the dune
        start_island = dune_loc + 1
        use_elev = (
            prof[start_island:-1]
            if start_island < (prof.size - 1)
            else np.array([], dtype=float)
        )

        # 8) Write interior into rows (top-down)
        if use_elev.size > 0:
            rows_to_copy = min(TOPO_ROWS, use_elev.size)
            topo_m[0:rows_to_copy, i] = use_elev[:rows_to_copy]

    # Diagnostics summary
    if dune_heights_found:
        mean_height = np.mean(dune_heights_found)
        min_height = np.min(dune_heights_found)
        max_height = np.max(dune_heights_found)
        mean_berm = np.mean(berm_elevations_found)
        min_berm = np.min(berm_elevations_found)
        max_berm = np.max(berm_elevations_found)

        print(f"[stats] {in_path.name}: processed {profiles_processed} profiles")
        print(f"        Dune heights: mean={mean_height:.2f}m, "
              f"min={min_height:.2f}m, max={max_height:.2f}m")
        print(f"        Berm elevations: mean={mean_berm:.2f}m, "
              f"min={min_berm:.2f}m, max={max_berm:.2f}m")

        if mean_height > 3.0:
            print(f"[WARN]  Mean dune height ({mean_height:.2f}m) is high but may be valid")
        elif mean_height < 0.5:
            print(f"[WARN]  Mean dune height ({mean_height:.2f}m) seems very low!")
        else:
            print(f"[GOOD]  Dune heights look reasonable!")

        if profiles_with_tall_dunes > 0:
            pct = 100 * profiles_with_tall_dunes / profiles_processed
            print(f"[WARN]  {profiles_with_tall_dunes} profiles ({pct:.0f}%) had dunes >{MAX_DUNE_HEIGHT_M}m (capped)")

    # Convert to decameters (× 0.1)
    topo_dm = topo_m * 0.1
    dune_dm = dune_m * 0.1

    # Save
    stem = in_path.stem
    topo_out = Path(topo_out_dir) / f"{stem}_topography_2009.npy"
    dune_out = Path(dune_out_dir) / f"{stem}_dune_2009.npy"
    np.save(topo_out, topo_dm)
    np.save(dune_out, dune_dm)
    print(f"[ok] wrote {topo_out.name} and {dune_out.name}\n")


def main():
    load_dir = Path(LOAD_PATH)
    topo_dir = Path(TOPO_SAVE_PATH)
    dune_dir = Path(DUNE_SAVE_PATH)
    topo_dir.mkdir(parents=True, exist_ok=True)
    dune_dir.mkdir(parents=True, exist_ok=True)

    # Process domain elevation arrays
    names = sorted([
        n for n in os.listdir(load_dir)
        if n.endswith(".npy") and n.startswith("domain_")
    ])
    print(f"[info] Found {len(names)} domain file(s) in {load_dir}\n")
    print("=" * 80)
    print("TRULY FIXED VERSION:")
    print("- Uses LOCAL berm elevation (minimum in window before dune)")
    print("- Not a fixed berm elevation value")
    print("=" * 80 + "\n")

    all_dune_heights = []
    all_berm_elevs = []

    for name in names:
        process_domain_file(load_dir / name, topo_dir, dune_dir)

    print("\n" + "=" * 80)
    print("PROCESSING COMPLETE")
    print("=" * 80)
    print(f"\nOutput directories:")
    print(f"  Topography: {topo_dir}")
    print(f"  Dunes:      {dune_dir}")
    print("\nREVIEW RESULTS:")
    print("- Expected mean dune heights: 1.5-3.0m")
    print("- Expected berm elevations: variable (0.5-9m depending on location)")
    print("- If most domains show good stats, proceed with CASCADE")
    print("=" * 80)


if __name__ == "__main__":
    main()