# ==============================================================================
# Hatteras CASCADE Dune & Topography Extractor
# INPUT  : numpy array shape = (alongshore_rows, cross_shore_cols)
# OUTPUT : topography_dunes and dune, in decameters (dam)
# ==============================================================================

import os, copy
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats


def remove_water_cols(domain_array, w_elev):
    n_cols = np.shape(domain_array)[1]
    c_list = []
    for c in range(n_cols):
        if not np.all(domain_array[:, c] == w_elev):
            c_list.append(c)
    min_c = min(c_list)
    max_c = max(c_list)
    domain = domain_array[:, min_c:max_c + 1]

    return domain

def remove_water_rows(domain_array, w_elev):
    n_rows = np.shape(domain_array)[0]
    r_list = []
    for r in range(n_rows):
        if not np.all(domain_array[r, :] == w_elev):
            r_list.append(r)
    min_r = min(r_list)
    max_r = max(r_list)
    domain = domain_array[min_r:max_r + 1, :]

    return domain


def process_domain_file(
        in_path: Path,
        topo_out_dir: Path,
        dune_out_dir: Path,
        MHW_M=0.421,
        BERM_ELEV_NAVD_M=1.95,
        BEACH_START_THR_M=0.75,
        DUNE_WINDOW_PX=15,
        SENTINEL_WATER_M=-3.0,
        TOPO_ROWS=200,
        ALONG_COLS=50,
        OCEAN_LOC="bottom",
        SHIFT_CELLS=False,
        use_const_interior=False,
        ) -> None:
    """Process a single domain elevation array and write topography and dune outputs."""
    arr = np.load(in_path).astype(float, copy=False)
    if arr.ndim != 2:
        print(f"[skip] {in_path.name}: expected 2D array, got {arr.ndim}D")
        return

    # flip domains so that the ocean is on the right since that is how the code was written
    if OCEAN_LOC == "left":
        arr = np.flip(arr)
    elif OCEAN_LOC == "bottom":
        arr = np.rot90(arr)
    elif OCEAN_LOC == "top":
        arr = np.flip(arr)  # now it is the same orientation as "bottom"
        arr = np.rot90(arr)

    # Subtract MHW and clamp <-1 m to -3.0 m (MHW-relative) before any slicing
    # note: lexi edited this since we are keeping the marsh cells which are between 0 and -3 m MHW
    arr = arr - MHW_M
    arr[arr < -3.0] = SENTINEL_WATER_M
    # remove columns that are all water (need to keep rows so that alongshore length is always 50)
    arr = remove_water_cols(arr, -3)

    n_along, n_cross = arr.shape
    TOPO_ROWS = n_cross

    # Messages
    if n_along < ALONG_COLS:
        print(
            f"[warn] {in_path.name}: alongshore={n_along} < {ALONG_COLS}; "
            f"trailing output cols remain sentinel."
            )
    elif n_along > ALONG_COLS:
        print(
            f"[warn] {in_path.name}: alongshore={n_along} > {ALONG_COLS}; "
            f"only first {ALONG_COLS} profiles used."
            )

    # Allocate fixed-size outputs
    topo_m = np.full((TOPO_ROWS, ALONG_COLS), fill_value=SENTINEL_WATER_M, dtype=float)
    dune_m = np.full((ALONG_COLS,), fill_value=SENTINEL_WATER_M, dtype=float)

    dune_loc_array = []

    # Process up to ALONG_COLS alongshore profiles
    n_cols_to_fill = min(ALONG_COLS, n_along)
    for i in range(n_cols_to_fill):
        # Profile as stored (land ... ocean-right)
        prof_lr = arr[i, :]  # this is now one column with marsh/back-barrier on top, ocean on bottom
        # Ensure index 0 == ocean
        prof = np.flip(prof_lr)  # ocean is now on top

        # 1) First index where z > 0.5 m (MHW-relative)
        idx = np.where(prof > BEACH_START_THR_M)[0]
        if idx.size == 0:
            # there are no cells above the beach threshold
            dune_loc_array.append(np.nan)  # just make the dune cell the last cell
        else:
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
            dune_loc_array.append(dune_loc)

            # 5) Island interior starts immediately landward of the dune
            start_island = dune_loc + 1
            use_elev = (
                prof[start_island:-1]
                if start_island < (prof.size - 1)
                else np.array([], dtype=float)
            )

            # 6) Write interior into rows (top-down); truncate/pad to TOPO_ROWS
            if use_elev.size > 0:
                rows_to_copy = min(TOPO_ROWS, use_elev.size)
                topo_m[0:rows_to_copy, i] = use_elev[:rows_to_copy]
                # remaining rows stay as sentinel

            # 7) Dune height above berm, min 0.1 m
            dune_h_m = dune_elev - (BERM_ELEV_NAVD_M - MHW_M)  # both MHW-relative
            if dune_h_m < 0.0:
                dune_h_m = 0.1
            dune_m[i] = dune_h_m

    if SHIFT_CELLS:
        # shift interior columns based on dune locations
        dune_loc_array = np.array(dune_loc_array)
        min_dune_loc = np.nanmin(dune_loc_array)  # this is the dune position that is closest to the ocean
        dune_shift = dune_loc_array - min_dune_loc
        avg_interior_per_row = np.mean(topo_m, axis=1)  # we are going to fill in cells with the average value of the ROW (alongshore direction)
        # add extra water cells to the end of the interior domain to make sure we dont lose elevation cells
        n_rows_to_add = np.nanmax(dune_shift)
        new_rows = np.ones([int(n_rows_to_add), np.shape(topo_m)[1]])*SENTINEL_WATER_M  # same number of cols as topo, new rows are just water
        topo_m = np.vstack((topo_m, new_rows))
        # add new cells to the top of the array, remove cells from the bottom (basically re-index)
        for c in range(np.shape(topo_m)[1]-1):
            current_col = topo_m[:,c]
            shift = dune_shift[c]
            if shift > 0 and not np.isnan(shift):
                new_cells = avg_interior_per_row[0:int(shift)]
                new_col = np.transpose(np.hstack((new_cells, current_col[0:-int(shift)])))
                topo_m[:,c] = new_col
    if use_const_interior:
        arr_ocean_top = np.rot90(arr)
        dune_loc_array = np.array(dune_loc_array)
        max_dune_loc = np.nanmax(dune_loc_array)
        min_dune_loc = np.nanmin(dune_loc_array)
        mode_dune_loc = stats.mode(dune_loc_array).mode
        start_island = int(max_dune_loc + 1)
        if "domain_3" in str(in_path):
            start_island = int(min_dune_loc + 1)
        elif "domain_19" in str(in_path):
            start_island = int(mode_dune_loc + 1)
        topo_m = arr_ocean_top[start_island:, :]

    # remove rows that are all water cells
    topo_m = remove_water_rows(topo_m, -3)

    topo_dm = topo_m * 0.1
    dune_dm = dune_m * 0.1

    # Save
    stem = in_path.stem  # e.g., "domain_7"
    topo_out = Path(topo_out_dir) / f"{stem}_interior_20192020.npy"
    dune_out = Path(dune_out_dir) / f"{stem}_dunes_20192020.npy"
    np.save(topo_out, topo_dm)
    np.save(dune_out, dune_dm)
    print(
        f"[ok] wrote {topo_out.name} (shape {topo_dm.shape}), "
        f"{dune_out.name} (len {dune_dm.size})"
        )

    return topo_m, dune_m


############################################################################################
# main script
############################################################################################
plt.rcParams["font.size"] = 14

# --- PATHS --------------------------------------------------------------
version = "v7"  # save version to append to folder name
LOAD_PATH = r"C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data\DEM\processed_domains\rotated_domains_npys"
TOPO_SAVE_PATH = r"C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data\DEM\processed_domains\cascade_domains\topography\rotated_domains_{0}_edited".format(version)
DUNE_SAVE_PATH = r"C:\Users\agfig\OneDrive - University of North Carolina at Chapel Hill\UNC\data\DEM\processed_domains\cascade_domains\dunes\rotated_domains_{0}_edited".format(version)

os.makedirs(TOPO_SAVE_PATH, exist_ok=True)
os.makedirs(DUNE_SAVE_PATH, exist_ok=True)

# --- CONSTANTS ----------------------------------------------------
MHW_M = 0.421              # meters (NAVD88)
BERM_ELEV_NAVD_M = 1.95    # meters (NAVD88)
BEACH_START_THR_M = 0.5    # meters (MHW-relative), strict '>' comparison
DUNE_WINDOW_PX = 5         # pixels
SENTINEL_WATER_M = -3.0    # meters (MHW-relative)
TOPO_ROWS = 200            # number of inland rows to write
ALONG_COLS = 50            # number of alongshore profiles
OCEAN_LOC = "bottom"       # "top", "bottom", "left", or "right"
shift_interior = False     # add cells to the beginning of the interior domain to keep it aligned
use_const_interior = True  # select a start row for the interior (most landward dune cell + 1)

load_dir = Path(LOAD_PATH)
topo_dir = Path(TOPO_SAVE_PATH)
dune_dir = Path(DUNE_SAVE_PATH)
topo_dir.mkdir(parents=True, exist_ok=True)
dune_dir.mkdir(parents=True, exist_ok=True)

# Process domain elevation arrays: domain_#.npy
names = sorted(
    [
        n for n in os.listdir(load_dir)
        if n.endswith(".npy") and n.startswith("domain_21")
        ]
    )
print(f"[info] Found {len(names)} domain file(s) in {load_dir}")
topo_domain = []
for name in names:
    topo_domain, dune_domain = process_domain_file(
        load_dir / name,
        topo_dir,
        dune_dir,
        MHW_M=MHW_M,
        BERM_ELEV_NAVD_M=BERM_ELEV_NAVD_M,
        BEACH_START_THR_M=BEACH_START_THR_M,
        DUNE_WINDOW_PX=DUNE_WINDOW_PX,
        SENTINEL_WATER_M=SENTINEL_WATER_M,
        TOPO_ROWS=TOPO_ROWS,
        ALONG_COLS=ALONG_COLS,
        OCEAN_LOC=OCEAN_LOC,
        SHIFT_CELLS=shift_interior,
        use_const_interior=use_const_interior
        )

domain = np.vstack((dune_domain+(BERM_ELEV_NAVD_M - MHW_M), topo_domain))
minz = -3
maxz = 5
xlabel = "alongshore distance (dam)"
ylabel = "cross-shore distance (dam)"
fig1 = plt.figure(figsize=[6,10])
ax1 = fig1.add_subplot(111)
mat1 = ax1.matshow(
    domain,
    cmap="terrain",
    vmin=minz,
    vmax=maxz,
    )
ax1.set_ylabel(ylabel)
ax1.set_title("dune window: {0}".format(DUNE_WINDOW_PX))
plt.gca().xaxis.tick_bottom()
plt.show()

