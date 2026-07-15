# ==============================================================================
# Hatteras CASCADE Dune & Topography Extractor
# INPUT  : numpy array shape = (alongshore_rows, cross_shore_cols)
# OUTPUT : topography_dunes and dune, in decameters (dam)
# ==============================================================================

import os, copy
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt


def plot_domain(domain, size, name, xlabel="alongshore distance (dam)", ylabel="cross-shore distance (dam)"):
    # plot domains

    minz = -3
    maxz = 5

    # initial
    fig = plt.figure(figsize=size)
    ax1 = fig.add_subplot(111)
    mat1 = ax1.matshow(
        domain,
        cmap="terrain",
        vmin=minz,
        vmax=maxz,
        )
    cbar = fig.colorbar(mat1)
    cbar.set_label('m MHW', rotation=270, labelpad=10)
    ax1.set_title(name)
    ax1.set_ylim([np.shape(domain)[0] - 1, 0])
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax1.xaxis.set_ticks_position("bottom")
    plt.show()


def plot_both_domains(gis_domain, casc_domain, size, xlabel="alongshore distance (dam)",
                      ylabel="cross-shore distance (dam)"):
    minz = -3
    maxz = 5
    # gis_domain = np.rot90(gis_domain)

    fig1 = plt.figure(figsize=size)
    ax1 = fig1.add_subplot(121)
    mat1 = ax1.matshow(
        gis_domain,
        cmap="terrain",
        vmin=minz,
        vmax=maxz,
        )
    ax1.set_title("GIS domain")
    ax1.set_ylabel(ylabel)

    plt.gca().xaxis.tick_bottom()
    ax2 = fig1.add_subplot(122)
    mat2 = ax2.matshow(
        casc_domain,
        cmap="terrain",
        vmin=minz,
        vmax=maxz,
        )
    cbar = fig1.colorbar(mat2)
    cbar.set_label('m MHW', rotation=270, labelpad=10)
    ax2.set_title("cascade domain")
    plt.gca().xaxis.tick_bottom()
    fig1.text(0.5, 0.01, xlabel, ha='center', va='center')
    fig1.tight_layout()
    plt.show()


def remove_water(domain_array, w_elev):
    n_cols = np.shape(domain_array)[1]
    c_list = []
    for c in range(n_cols):
        if not np.all(domain_array[:, c] == w_elev):
            c_list.append(c)
    min_c = min(c_list)
    max_c = max(c_list)
    domain = domain_array[:, min_c:max_c + 1]

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
    arr = remove_water(arr, -3)

    n_along, n_cross = arr.shape
    # note: lexi changed this bc she has variable barrier widths and wants to keep the full domain for now
    ALONG_COLS = n_along  # should always be 50 anyway
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
        dune_loc_array.append(dune_loc)

        # 5) Island interior starts immediately landward of the dune
        start_island = dune_loc + 1
        use_elev = (
            prof[start_island:-1]
            if start_island < (prof.size - 1)
            else np.array([], dtype=float)
        )

        start_channel = len(use_elev)

        # 6) Write interior into rows (top-down); truncate/pad to TOPO_ROWS
        if use_elev.size > 0:
            rows_to_copy = min(TOPO_ROWS, use_elev.size, start_channel)
            topo_m[0:rows_to_copy, i] = use_elev[:rows_to_copy]
            # remaining rows stay as sentinel

        # 7) Dune height above berm, min 0.1 m
        dune_h_m = dune_elev - (BERM_ELEV_NAVD_M - MHW_M)  # both MHW-relative
        if dune_h_m < 0.0:
            dune_h_m = 0.1
        dune_m[i] = dune_h_m

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


############################################################################################
# main script
############################################################################################
plt.rcParams["font.size"] = 14

# --- PATHS --------------------------------------------------------------
version = "v3"  # save version to append to folder name
LOAD_PATH = r"C:\Users\Lexi\Documents\UNC\data\DEM\processed_domains\rotated_domains_npys"
TOPO_SAVE_PATH = r"C:\Users\Lexi\Documents\UNC\data\DEM\processed_domains\cascade_domains\topography\rotated_domains_{0}".format(version)
DUNE_SAVE_PATH = r"C:\Users\Lexi\Documents\UNC\data\DEM\processed_domains\cascade_domains\dunes\rotated_domains_{0}".format(version)

os.makedirs(TOPO_SAVE_PATH, exist_ok=True)
os.makedirs(DUNE_SAVE_PATH, exist_ok=True)

# --- CONSTANTS ----------------------------------------------------
MHW_M = 0.421             # meters (NAVD88)
BERM_ELEV_NAVD_M = 1.95   # meters (NAVD88)
BEACH_START_THR_M = 0.50  # meters (MHW-relative), strict '>' comparison
DUNE_WINDOW_PX = 30       # pixels
SENTINEL_WATER_M = -3.0   # meters (MHW-relative)
TOPO_ROWS = 200           # number of inland rows to write
ALONG_COLS = 50           # number of alongshore profiles
OCEAN_LOC = "bottom"      # "top", "bottom", "left", or "right"

load_dir = Path(LOAD_PATH)
topo_dir = Path(TOPO_SAVE_PATH)
dune_dir = Path(DUNE_SAVE_PATH)
topo_dir.mkdir(parents=True, exist_ok=True)
dune_dir.mkdir(parents=True, exist_ok=True)

# Process domain elevation arrays: domain_#.npy
names = sorted(
    [
        n for n in os.listdir(load_dir)
        if n.endswith(".npy") and n.startswith("domain_")
        ]
    )
print(f"[info] Found {len(names)} domain file(s) in {load_dir}")

for name in names:
    process_domain_file(
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
        )


