# ==============================================================================
# PEA_2011_dune_topo_extractor_local_2row_dunes_v3.py
#
# Pea Island CASCADE dune & interior topography extractor for domains 80–119,
# with a freehand dune band, two real dune rows, local interior alignment,
# and the complete 200-row water-padded interior preserved.
#
# INPUT  : domain_#.npy, shape = (alongshore_rows, cross_shore_cols),
#          elevation in m NAVD88
# OUTPUT : interior topography and dune height arrays, in decameters (dam)
#          + a JSON of the per-domain dune search windows (re-runnable)
#          + a middle-domain cross-section figure for dune-location review
#          + an aligned binary road-mask overlay in the interactive and QC plots
#
# Adapted from Lexi's dune_topo_extractor_v3.py and the earlier
# dune_topo_extractor_from_GIS.py.
#
# WORKFLOW
#   1. MODE = "pick"          -> step through domains and DRAW a dune-search
#                                band on the map; saves after every domain
#   2. MODE = "run"           -> extract using the saved JSON bands
#   3. MODE = "pick_and_run"  -> both in one pass
#
# PICKER LAYOUT / KEYS
#   The picker draws the domain with the OCEAN AT THE BOTTOM: cross-shore runs
#   vertically (cell 0 = ocean, landward upward), alongshore runs left-right.
#   Left panel = elevation map, right panel = profile stack, shared y-axis.
#
#   left-drag on map : draw the ACTIVE boundary as a tilted/curved freehand line
#   1                : activate the seaward boundary
#   2                : activate the landward boundary
#   u                : undo the last boundary edit
#   enter            : accept the locally varying search band
#   r                : reset both boundaries to the default horizontal band
#   s                : skip this domain (use DEFAULT_WINDOW_PX)
#   q / esc          : quit picking, keep everything saved so far
# ==============================================================================

from __future__ import annotations

import csv
import json
import os
import re
from datetime import datetime
from pathlib import Path

import numpy as np

import matplotlib

# IMPORTANT: choose the GUI backend BEFORE importing pyplot. ``force=True``
# prevents PyCharm/Jupyter from silently replacing TkAgg with a noninteractive
# inline backend. The picker must open in its own desktop window.
try:
    import tkinter  # noqa: F401
    matplotlib.use("TkAgg", force=True)
except Exception as exc:
    raise RuntimeError(
        "The dune picker requires TkAgg, but tkinter/TkAgg could not be loaded. "
        "On macOS, run this script from the PyCharm Run button or Terminal and "
        "make sure your Python installation includes tkinter."
    ) from exc

import matplotlib.pyplot as plt
from matplotlib.colors import FuncNorm, ListedColormap
from matplotlib.transforms import blended_transform_factory

plt.rcParams["font.size"] = 11

# ==============================================================================
# CONFIG
# ==============================================================================

# --- PEA ISLAND DOMAIN CONFIGURATION ------------------------------------
FIRST_REAL_DOMAIN_ID = 80
LAST_REAL_DOMAIN_ID = 119
REAL_DOMAIN_IDS = list(range(FIRST_REAL_DOMAIN_ID, LAST_REAL_DOMAIN_ID + 1))
NUM_REAL_DOMAINS = len(REAL_DOMAIN_IDS)

# The input folder must contain these 40 files:
#     resampled_2011_domains_80.npy, ..., resampled_2011_domains_119.npy
# Each file must be a 2-D elevation array in meters NAVD88.
#
# The input filenames are preserved on disk, but the script standardizes all
# picker keys and CASCADE outputs to domain_80, ..., domain_119.

# --- MODE ---------------------------------------------------------------
MODE = "run"  # use saved bands when available; then review and extract
PICK_DOMAINS = REAL_DOMAIN_IDS

# False skips the old stand-alone dune-band picker when a saved band exists and
# opens the new combined plan-view/cross-section editor directly. Set True only
# when you also want to redraw the dune-search band before choosing the interior.
REPICK_EXISTING = False

# Stop immediately if Matplotlib is using a noninteractive backend rather than
# silently extracting dunes without giving you a drawing window.
REQUIRE_INTERACTIVE_PICKER = True
SAVE_QC_FIGS = True        # per-domain dune-detection QC figure
SAVE_COMPARISON_FIGS = True  # per-domain raw GIS vs processed CASCADE input figure
SAVE_CROSS_SECTION_FIGS = True   # save the final interactive partition QC
SHOW_CROSS_SECTION_FIGS = True     # display plan view and cross-section together
REQUIRE_CROSS_SECTION_APPROVAL = True
# Interactive partition picker controls:
#   PLAN VIEW left-drag = draw the interior-start line across the domain
#   PLAN VIEW right-click = select an alongshore profile
#   CROSS-SECTION click/drag = set interior start for the selected profile
#   Left/Right arrows = move through profiles
#   U = undo, R = reset to immediately landward of dune, Enter/A = accept
#   NOTE: the manually chosen interior start may be seaward or landward of the dune.
#   S = skip domain, Q/Esc = stop run
CROSS_SECTION_CONTEXT_PX = 8       # extra cells around the local search band
SAVE_SETTINGS_SHEET = True   # per-domain settings/results sheet (csv + xlsx)

# --- PATHS --------------------------------------------------------------
# Everything for one settings variant lands in ONE run folder, so comparing
# versions means comparing two directories:
#
#   data/PeaIsland_init/dune_topo_new/
#       picks/
#           PEA_dune_search_windows_2011_modified_RSLR.json  <- saved picks
#       2011_v1/
#           RUN_MANIFEST.txt                     <- every setting that made this folder
#           PEA_dune_topo_settings_2011_v1.xlsx  <- per-domain sheet (+ .csv)
#           PEA_dune_topo_summary_2011_v1.png    <- all domains on one page
#           topography\  domain_80_topography_2011.npy   <- CASCADE reads these two
#           dunes\       domain_80_dune_2011.npy
#           figures\
#               gis_vs_processed\  domain_080_gis_vs_processed.png
#               qc\                domain_080_qc.png
#       2011_v2/   ... same shape, nothing shared, nothing overwritten
#
# NOTE: topography\ and dunes\ moved. Point the hindcast runner at
#       RUN_DIR\topography and RUN_DIR\dunes for whichever version you're using.
VERSION = "Roya_1996_cross_ROAD_2row"             # local two-row dune + profile-following interior
DEM_YEAR = "2011"

# A simple folder name only -- do not put the full path in DEM_NAME.
DEM_NAME = "Hybrid_DEM"
RUN_NAME = f"{DEM_YEAR}_{VERSION}"

PROJECT_ROOT = Path(r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE")
INIT_ROOT = PROJECT_ROOT / "data" / "PeaIsland_init"

# Folder containing resampled_2011_domains_80.npy through
# resampled_2011_domains_119.npy.
LOAD_PATH = PROJECT_ROOT / "data" / DEM_NAME

# --- ROAD MASK OVERLAY --------------------------------------------------
# The road masks must have been created in ArcGIS with the SAME snap raster,
# extent, cell size, CRS, clipping polygons, and RasterToNumPyArray workflow as
# the elevation domains. Values may be 0/1 or any positive road value; this
# script converts every value > 0 to True.
SHOW_ROAD = True
REQUIRE_ROAD_MASKS = True
ROAD_MASK_FOLDER = PROJECT_ROOT / "data" /"PeaIsland_init" /"road_masks"
ROAD_MASK_PATTERNS = [
    "road_mask_domain_{domain}.npy",
    "road_mask_domain_{domain:03d}.npy",
    "domains_{domain}.npy",
    "resampled_road_domains_{domain}.npy",
]
ROAD_LABEL = "N.C. 12 road"
ROAD_COLOR = "#111111"
ROAD_EDGE_COLOR = "#FFFFFF"
ROAD_PLAN_ALPHA = 0.52
ROAD_CROSS_SECTION_ALPHA = 0.16
SHOW_ROAD_CENTER_LINE = True

DUNE_TOPO_ROOT = INIT_ROOT / "1996"/"dune_topo_new"
RUN_DIR = DUNE_TOPO_ROOT / RUN_NAME
TOPO_SAVE_PATH = RUN_DIR / "topography"
DUNE_SAVE_PATH = RUN_DIR / "dunes"
FIG_DIR_CMP = RUN_DIR / "figures" / "gis_vs_processed"
FIG_DIR_QC = RUN_DIR / "figures" / "qc"
FIG_DIR_XSEC = RUN_DIR / "figures" / "middle_cross_sections"
SHEET_SAVE_PATH = RUN_DIR / f"PEA_dune_topo_settings_{RUN_NAME}"   # .csv + .xlsx
SUMMARY_FIG_PATH = RUN_DIR / f"PEA_dune_topo_summary_{RUN_NAME}.png"
ISLAND_FIG_PATH = RUN_DIR / f"PEA_dune_topo_island_offsets_{RUN_NAME}.png"
PLAN_FIG_STEM = f"PEA_dune_topo_island_planview_{RUN_NAME}"  # + _{year}.png
MANIFEST_PATH = RUN_DIR / "RUN_MANIFEST.txt"

# Picks live OUTSIDE the run folder. They are the only artifact here you cannot
# regenerate, and they describe where the dune sits in the DEM -- not which
# settings variant you're testing -- so v2, v3... reuse them by default. Set
# PICK_SET = RUN_NAME instead if you want a version to carry its own picks.
PICK_SET = DEM_NAME
PICKS_DIR = DUNE_TOPO_ROOT / "picks"
WINDOW_JSON = Path(r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/1996/dune_topo_new/picks/PEA_dune_search_windows_Hybrid_DEM.json")

TAG = DEM_YEAR             # appended to the .npy filenames CASCADE reads

# --- OPTIONAL ISLAND OFFSETS ---------------------------------------------
# Offsets are NOT required to create the dune and interior .npy files.
# Keep the two SAVE flags False until a Pea Island offset CSV is available.
#
# Supported offset CSV layouts:
#   1) 40 real-domain rows for domains 80..119, or
#   2) 111 padded rows = 71 left buffers + 40 real domains + 0 right buffers.
#
# Add a file like this when available:
# OFFSET_FILES = {
#     2011: OFFSET_DIR / "Island_Dune_Offsets_2011_PADDED_71.csv",
# }
SAVE_ISLAND_FIG = False
OFFSET_DIR = INIT_ROOT / "island_offset"
OFFSET_FILES = {}

CELL_SIZE_M = 10.0
N_LEFT_BUFFER_DOMAINS = 71
N_RIGHT_BUFFER_DOMAINS = 0
OFFSET_COLUMN = 0
OFFSET_ROW_ORDER = "ascending"   # "ascending" = 80..119; "descending" = 119..80
OFFSET_SEAWARD_POSITIVE = True

# Plan-view canvas, reproducing plot_initialization_poster_no_border.py exactly:
#   offset_cells = round(offset_m / 10); each domain's topo row 0 (ocean side)
#   lands on canvas row = offset_cells; alongshore flipped with np.fliplr.
SAVE_ISLAND_PLAN_FIG = False
ISLAND_FLIP_ALONGSHORE = True  # the np.fliplr in the poster script
ISLAND_INCLUDE_DUNE = True     # write the dune crest into canvas row offset-1
ISLAND_ELEV_MIN_M = -1.0       # poster value. Raising this flattens land contrast:
ISLAND_ELEV_MAX_M = 4.0        #   land maps to ramp 0.35-1.0, so max=4 puts a 2 m
                               #   cell at tan, max=5 at pale yellow-green.
ISLAND_SEA_LEVEL_POS = 0.35    # colormap position for 0 m, as in the poster
ISLAND_SENTINEL_AS_OCEAN = False  # False = poster behaviour: sentinel (-3 m) water
                                  #   cells clip to the dark navy bottom of terrain,
                                  #   so the model's cross-shore extent stays visible
                                  #   and only outside-canvas NaN is light blue.
                                  # True = sentinel also renders light blue.
ISLAND_OCEAN_COLOR = "#b0cfe8"    # poster set_bad colour
# Display only, does not touch the .npy files. ONE FIGURE PER YEAR PER MODE, so a
# single run gives you both versions to compare:
#   "trimmed" = each domain only as tall as its own island, as stored
#   "padded"  = every domain given the same cross-shore extent (ISLAND_PAD_ROWS),
#               sentinel-filled landward, so the water behind each domain shows as
#               the navy wedge in the poster figure
# If TRIM_INTERIOR_ROWS = False the arrays are already padded and the two match.
ISLAND_CROSS_SHORE_MODES = ["trimmed", "padded"]
ISLAND_PAD_ROWS = 200               # cells of constant cross-shore extent when padded.
                                    # 200 = 2000 m, equal to TOPO_ROWS, so nothing is
                                    #       cropped (the poster look)
                                    # 100 = 1000 m, tighter; the script warns if that
                                    #       crops real land off the sound side
# --- OPTIONAL PEA ISLAND SECTIONS ----------------------------------------
# Leave empty unless you want named bands on summary figures. The extractor does
# not assume which Pea Island domain IDs correspond to specific road/habitat areas.
# Example: SECTIONS = [((80, 89), "Section A"), ((90, 99), "Section B")]
SECTIONS = []

# --- DATUMS / THRESHOLDS ------------------------------------------------
MHW_M = 0.36               # m NAVD88
BERM_ELEV_NAVD_M = 1.70    # m NAVD88 (matches the 1.7 m storm/collision threshold)
BEACH_START_THR_M = 0.50   # m MHW-relative, strict '>' comparison
WATER_CLAMP_M = -3.0       # m MHW-relative; below this -> sentinel.
                           #   -3.0 keeps back-barrier marsh cells (Lexi's v3 edit)
                           #   -1.0 was the old dune_topo_extractor_from_GIS behavior
SENTINEL_WATER_M = -3.0    # m MHW-relative
MIN_DUNE_H_M = 0.1         # m, floor on dune height above berm

# --- GEOMETRY -----------------------------------------------------------
TOPO_ROWS = 200            # max interior rows written
ALONG_COLS = 50            # alongshore profiles per domain (500 m / 10 m)
OCEAN_LOC = "right"        # verify for the 2011 Pea DEM; change if needed
ALONGSHORE_FLIP = False    # True reverses alongshore order after orienting

# --- DUNE SEARCH --------------------------------------------------------
DEFAULT_WINDOW_PX = 8      # fallback band width (px landward of beach start)
CLIP_WINDOW_TO_BEACH = True  # per profile, start search at max(i0_line, beach_start)
                             # so a drawn band cannot wander onto the wet beach
                             # where the shoreline curves within a domain
MIN_DRAWN_BAND_WIDTH_PX = 2  # minimum local distance between the two drawn lines
FREEHAND_SMOOTH_PX = 3       # odd moving-average width; 1 disables smoothing
AUTO_SWITCH_DRAW_BOUNDARY = True  # after a stroke, switch to the other boundary

# --- DUNE / INTERIOR PARTITION -----------------------------------------
DUNE_ROWS = 2              # Number of real 10 m cross-shore dune cells.
                           # Saved dune .npy shape = (ALONG_COLS, DUNE_ROWS).
                           # For Pea Island this is (50, 2): row 0 starts at
                           # the detected crest and row 1 is the next landward cell.
FOLLOW_LOCAL_DUNE_LINE = True
                           # Each alongshore profile is re-indexed independently:
                           #   dune row 0 = detected crest cell
                           #   interior row 0 = the next landward GIS cell
                           # This makes the rectangular output follow the curved dune.
FILL_MISSING_DUNE = True   # interpolate missing crest locations alongshore and fill
                           # invalid dune heights with MIN_DUNE_H_M.
TRIM_INTERIOR_ROWS = False # Preserve all TOPO_ROWS = 200 rows, including
                           # back-barrier/landward sentinel water at -3 m MHW.


# ==============================================================================
# ARRAY HELPERS
# ==============================================================================

def water_col_bounds(domain_array: np.ndarray, w_elev: float) -> tuple[int, int]:
    """First/last (inclusive) columns that are not entirely water."""
    keep = [c for c in range(domain_array.shape[1])
            if not np.all(domain_array[:, c] == w_elev)]
    if not keep:
        raise ValueError("domain is entirely water after clamping")
    return min(keep), max(keep)


def remove_water_rows(domain_array: np.ndarray, w_elev: float) -> np.ndarray:
    """Trim leading/trailing rows that are entirely water."""
    keep = [r for r in range(domain_array.shape[0])
            if not np.all(domain_array[r, :] == w_elev)]
    if not keep:
        return domain_array[:0, :]
    return domain_array[min(keep):max(keep) + 1, :]




def trim_landward_water_rows(domain_array: np.ndarray, w_elev: float) -> np.ndarray:
    """Remove only trailing landward all-water rows; preserve local row 0.

    With profile-following re-indexing, row 0 has a physical meaning: it is the
    first GIS cell immediately landward of the reserved dune cells. Therefore, unlike
    ``remove_water_rows``, this function never trims rows from the ocean-side
    beginning of the processed interior.
    """
    keep = [r for r in range(domain_array.shape[0])
            if not np.all(domain_array[r, :] == w_elev)]
    if not keep:
        return domain_array[:0, :]
    return domain_array[:max(keep) + 1, :]


def orient_ocean_right(arr: np.ndarray, ocean_loc: str) -> np.ndarray:
    """
    Return an (alongshore, cross_shore) array with the ocean in the LAST column.

    NOTE: v3's "top"/"left" branches used np.flip(arr), which flips BOTH axes and
    silently reverses the alongshore order. These branches do not.
    """
    if ocean_loc == "right":
        out = arr
    elif ocean_loc == "left":
        out = arr[:, ::-1]
    elif ocean_loc == "bottom":
        out = np.rot90(arr)        # last row -> last column
    elif ocean_loc == "top":
        out = np.rot90(arr, -1)    # first row -> last column
    else:
        raise ValueError(f"OCEAN_LOC must be right/left/top/bottom, got {ocean_loc!r}")
    if ALONGSHORE_FLIP:
        out = out[::-1, :]
    return np.ascontiguousarray(out)


def domain_number(name: str | Path) -> int | None:
    """
    Extract a Pea Island domain ID from either input or standardized names.

    Recognized examples:
        resampled_2011_domains_80.npy
        resampled_2011_domains_80
        domain_80.npy
        domain_80_topography_2011.npy
    """
    filename = Path(name).name
    m = re.search(r"domains?_(\d+)", filename, flags=re.IGNORECASE)
    return int(m.group(1)) if m else None


def canonical_domain_stem(name: str | Path) -> str:
    """Return the standardized CASCADE/picker name, e.g. domain_80."""
    n = domain_number(name)
    if n is None:
        raise ValueError(f"cannot determine domain ID from {name!r}")
    return f"domain_{n}"


def is_real_domain_file(name: str | Path) -> bool:
    """True for a .npy input whose parsed ID is within domains 80..119."""
    p = Path(name)
    n = domain_number(p.name)
    return p.suffix.lower() == ".npy" and n in REAL_DOMAIN_IDS


def natural_key(name: str | Path) -> tuple:
    n = domain_number(name)
    return (n if n is not None else 10**9, str(name))


def section_for(stem: str) -> str:
    """Optional Pea Island section label for a domain."""
    n = domain_number(stem)
    if n is None:
        return ""
    for (lo, hi), label in SECTIONS:
        if lo <= n <= hi:
            return label
    return ""


def fig_stem(stem: str) -> str:
    """Zero-padded figure name so 90+ domains sort correctly in Explorer."""
    n = domain_number(stem)
    return f"domain_{n:03d}" if n is not None else stem


def fig_title(stem: str) -> str:
    sec = section_for(stem)
    return f"{stem}  ({sec})" if sec else stem



# ==============================================================================
# ROAD-MASK HELPERS
# ==============================================================================

def road_mask_candidates(domain_id: int) -> list[Path]:
    """Return the configured candidate filenames for one domain."""
    return [
        Path(ROAD_MASK_FOLDER) / pattern.format(
            domain=domain_id,
            domain_id=domain_id,
        )
        for pattern in ROAD_MASK_PATTERNS
    ]


def find_road_mask_path(domain_id: int) -> Path | None:
    """Return the first existing road-mask file for one domain."""
    for candidate in road_mask_candidates(domain_id):
        if candidate.is_file():
            return candidate
    return None


def load_road_mask_for_domain(
        dem_path: Path,
        expected_shape: tuple[int, int],
        ) -> tuple[np.ndarray | None, Path | None]:
    """Load a binary road mask in the ORIGINAL GIS-array orientation.

    The mask must have exactly the same shape as the matching DEM .npy. No
    transpose or resizing is performed because that could hide a grid-alignment
    error. The DEM and road are oriented and trimmed together later.
    """
    if not SHOW_ROAD:
        return None, None

    domain_id = domain_number(dem_path.name)
    if domain_id is None:
        raise ValueError(f"cannot determine domain ID for road mask: {dem_path.name}")

    road_path = find_road_mask_path(domain_id)
    if road_path is None:
        tried = "\n".join(f"  - {p}" for p in road_mask_candidates(domain_id))
        message = (
            f"No road mask found for domain {domain_id}. Tried:\n{tried}"
        )
        if REQUIRE_ROAD_MASKS:
            raise FileNotFoundError(message)
        print(f"[road warning] {message}")
        return np.zeros(expected_shape, dtype=bool), None

    road = np.load(road_path)
    road = np.squeeze(road)
    if road.ndim != 2:
        raise ValueError(
            f"{road_path.name}: expected a 2-D road mask, got shape {road.shape}"
        )
    if tuple(road.shape) != tuple(expected_shape):
        raise ValueError(
            f"{road_path.name}: road shape {road.shape} does not match "
            f"DEM shape {expected_shape}. Re-export it in ArcGIS using the exact "
            "same snap raster, extent, cell size, and domain polygon."
        )

    binary = np.isfinite(road) & (road > 0)
    print(
        f"[road] {dem_path.name}: {road_path.name}, "
        f"{int(binary.sum())} road cell(s)"
    )
    return np.ascontiguousarray(binary), road_path


def validate_road_mask(
        road_mask: np.ndarray | None,
        expected_shape: tuple[int, int],
        ) -> np.ndarray | None:
    """Validate an already oriented/trimmed road mask before plotting."""
    if road_mask is None:
        return None
    road = np.asarray(road_mask, dtype=bool)
    if road.shape != expected_shape:
        raise ValueError(
            f"road mask shape {road.shape} does not match elevation shape "
            f"{expected_shape}"
        )
    return road


def road_profile_positions(
        road_mask: np.ndarray | None,
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return seaward edge, landward edge, and center cell for each profile."""
    if road_mask is None:
        return (np.array([], dtype=float),) * 3

    n_along = road_mask.shape[0]
    seaward = np.full(n_along, np.nan, dtype=float)
    landward = np.full(n_along, np.nan, dtype=float)
    center = np.full(n_along, np.nan, dtype=float)

    for i in range(n_along):
        cells = np.flatnonzero(road_mask[i])
        if cells.size:
            seaward[i] = float(cells.min())
            landward[i] = float(cells.max())
            center[i] = float(np.mean(cells))

    return seaward, landward, center


def contiguous_true_runs(values: np.ndarray) -> list[tuple[int, int]]:
    """Return inclusive index ranges for contiguous True cells."""
    idx = np.flatnonzero(np.asarray(values, dtype=bool))
    if idx.size == 0:
        return []
    split_at = np.where(np.diff(idx) > 1)[0] + 1
    groups = np.split(idx, split_at)
    return [(int(g[0]), int(g[-1])) for g in groups if g.size]


def add_road_plan_overlay(
        ax,
        road_mask: np.ndarray | None,
        *,
        label: str = ROAD_LABEL,
        zorder: float = 5.0,
        ) -> None:
    """Overlay exact binary road cells on an alongshore/cross-shore plan view."""
    if road_mask is None or not np.any(road_mask):
        return

    n_along, n_cross = road_mask.shape
    display = np.ma.masked_where(~road_mask.T, np.ones((n_cross, n_along)))
    ax.imshow(
        display,
        aspect="auto",
        origin="lower",
        extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
        interpolation="nearest",
        cmap=ListedColormap([ROAD_COLOR]),
        vmin=0.0,
        vmax=1.0,
        alpha=ROAD_PLAN_ALPHA,
        zorder=zorder,
    )

    # A contour makes the road visible on both dark water and light land.
    if np.any(~road_mask):
        ax.contour(
            np.arange(n_along),
            np.arange(n_cross),
            road_mask.T.astype(float),
            levels=[0.5],
            colors=ROAD_EDGE_COLOR,
            linewidths=0.8,
            zorder=zorder + 0.2,
        )

    seaward, landward, center = road_profile_positions(road_mask)
    valid = np.isfinite(center)
    if SHOW_ROAD_CENTER_LINE and valid.any():
        ax.plot(
            np.arange(n_along), center,
            color=ROAD_COLOR, lw=1.4, ls="--",
            label=f"{label} center", zorder=zorder + 0.4,
        )

    # Legend proxy for the filled road cells.
    ax.plot([], [], color=ROAD_COLOR, lw=6, alpha=ROAD_PLAN_ALPHA,
            label=label)


def add_road_profile_envelope(ax, road_mask: np.ndarray | None) -> None:
    """Show the road's overall cross-shore envelope on a profile-stack panel."""
    if road_mask is None or not np.any(road_mask):
        return
    seaward, landward, center = road_profile_positions(road_mask)
    valid = np.isfinite(center)
    if not valid.any():
        return
    ax.axhspan(
        float(np.nanmin(seaward[valid])) - 0.5,
        float(np.nanmax(landward[valid])) + 0.5,
        color=ROAD_COLOR,
        alpha=0.08,
        zorder=0,
        label="road cross-shore envelope",
    )
    ax.axhline(
        float(np.nanmedian(center[valid])),
        color=ROAD_COLOR,
        lw=1.3,
        ls="--",
        label="median road position",
    )


def road_profile_summary(
        road_mask: np.ndarray | None,
        profile_index: int,
        dune_cell: int | None = None,
        ) -> dict:
    """Return road cells and dune-to-road distance for one alongshore profile."""
    result = {
        "present": False,
        "cells": np.array([], dtype=int),
        "seaward_cell": np.nan,
        "landward_cell": np.nan,
        "center_cell": np.nan,
        "dune_to_road_m": np.nan,
        "note": "road not present in this profile",
    }
    if road_mask is None or not (0 <= profile_index < road_mask.shape[0]):
        return result

    cells = np.flatnonzero(road_mask[profile_index])
    if cells.size == 0:
        return result

    center = float(np.mean(cells))
    result.update({
        "present": True,
        "cells": cells,
        "seaward_cell": int(cells.min()),
        "landward_cell": int(cells.max()),
        "center_cell": center,
    })
    if dune_cell is not None and dune_cell >= 0:
        distance = (center - float(dune_cell)) * CELL_SIZE_M
        result["dune_to_road_m"] = distance
        result["note"] = (
            f"road center cell {center:.1f}; dune-to-road {distance:+.0f} m "
            "(+ = road landward)"
        )
    else:
        result["note"] = f"road center cell {center:.1f}"
    return result


def add_road_cross_section_overlay(
        ax,
        road_mask: np.ndarray | None,
        profile_index: int,
        profile_elev: np.ndarray,
        *,
        dune_cell: int | None = None,
        include_label: bool = True,
        ) -> dict:
    """Shade exact road cells and mark their elevations on one cross-section."""
    info = road_profile_summary(road_mask, profile_index, dune_cell)
    if not info["present"]:
        return info

    first = True
    for start, end in contiguous_true_runs(road_mask[profile_index]):
        left_m = max(0.0, (start - 0.5) * CELL_SIZE_M)
        right_m = (end + 0.5) * CELL_SIZE_M
        ax.axvspan(
            left_m, right_m,
            color=ROAD_COLOR,
            alpha=ROAD_CROSS_SECTION_ALPHA,
            zorder=1,
            label=ROAD_LABEL if include_label and first else None,
        )
        first = False

    cells = info["cells"]
    valid = (
        (cells >= 0)
        & (cells < len(profile_elev))
        & np.isfinite(profile_elev[cells])
        & (profile_elev[cells] > SENTINEL_WATER_M + 1e-9)
    )
    if valid.any():
        good_cells = cells[valid]
        ax.scatter(
            good_cells * CELL_SIZE_M,
            profile_elev[good_cells],
            marker="s", s=48,
            facecolor=ROAD_COLOR,
            edgecolor="white",
            linewidth=0.7,
            zorder=10,
            label="road-mask cell elevation" if include_label else None,
        )

    ax.axvline(
        info["center_cell"] * CELL_SIZE_M,
        color=ROAD_COLOR,
        lw=1.7,
        ls="--",
        zorder=9,
        label="road center" if include_label else None,
    )
    return info


def processed_interior_road_mask(
        road_mask: np.ndarray | None,
        interior_start_line: np.ndarray,
        n_rows: int,
        n_cols: int,
        ) -> np.ndarray:
    """Map source road cells into the profile-shifted CASCADE interior grid."""
    out = np.zeros((n_rows, n_cols), dtype=bool)
    if road_mask is None:
        return out
    n_fill = min(n_cols, road_mask.shape[0], len(interior_start_line))
    for i in range(n_fill):
        start = int(interior_start_line[i])
        if start < 0:
            continue
        source_cells = np.flatnonzero(road_mask[i])
        destination = source_cells - start
        keep = (destination >= 0) & (destination < n_rows)
        out[destination[keep], i] = True
    return out


def add_processed_road_overlay(ax, road_grid: np.ndarray) -> None:
    """Overlay a (cross-shore rows, alongshore cols) road grid on CASCADE topo."""
    if road_grid is None or not np.any(road_grid):
        return
    n_rows, n_cols = road_grid.shape
    display = np.ma.masked_where(~road_grid, np.ones_like(road_grid, dtype=float))
    ax.imshow(
        display,
        aspect="auto",
        origin="lower",
        extent=[-0.5, n_cols - 0.5, -0.5, n_rows - 0.5],
        interpolation="nearest",
        cmap=ListedColormap([ROAD_COLOR]),
        vmin=0.0,
        vmax=1.0,
        alpha=ROAD_PLAN_ALPHA,
        zorder=5,
    )
    ax.plot([], [], color=ROAD_COLOR, lw=6, alpha=ROAD_PLAN_ALPHA,
            label=ROAD_LABEL)


def attach_road_statistics(res: dict, dom: dict) -> None:
    """Add per-domain road coverage and dune-to-road distances to result metadata."""
    road_mask = dom.get("road_mask")
    res["road_mask_file"] = dom.get("road_mask_name", "")
    if road_mask is None or not np.any(road_mask):
        res.update({
            "road_profiles": 0,
            "mean_dune_to_road_m": np.nan,
            "min_dune_to_road_m": np.nan,
            "max_dune_to_road_m": np.nan,
        })
        return

    _, _, road_center = road_profile_positions(road_mask)
    dune_loc = np.asarray(res.get("dune_loc_used", []), dtype=float)
    n = min(len(road_center), len(dune_loc), ALONG_COLS)
    valid = np.isfinite(road_center[:n]) & np.isfinite(dune_loc[:n]) & (dune_loc[:n] >= 0)
    distances = (road_center[:n][valid] - dune_loc[:n][valid]) * CELL_SIZE_M
    res["road_profiles"] = int(np.isfinite(road_center[:n]).sum())
    if distances.size:
        res["mean_dune_to_road_m"] = float(np.mean(distances))
        res["min_dune_to_road_m"] = float(np.min(distances))
        res["max_dune_to_road_m"] = float(np.max(distances))
    else:
        res["mean_dune_to_road_m"] = np.nan
        res["min_dune_to_road_m"] = np.nan
        res["max_dune_to_road_m"] = np.nan


# ==============================================================================
# LOAD / PREP
# ==============================================================================

def load_profiles(in_path: Path) -> dict:
    """Load one elevation domain and its aligned binary road mask.

    Returns
    -------
    raw
        Original GIS elevations in m NAVD88, oriented OCEAN-FIRST and untrimmed.
    z
        Elevations relative to MHW, clamped and water-trimmed, OCEAN-FIRST.
    start_beach
        First cell in each profile where elevation exceeds BEACH_START_THR_M.
    c0
        Cross-shore trim offset between ``raw`` and ``z``.
    road_raw
        Binary road mask in the same orientation/shape as ``raw``.
    road_mask
        Binary road mask in the same orientation/shape as ``z``.
    """
    arr = np.load(in_path).astype(float, copy=False)
    if arr.ndim != 2:
        raise ValueError(f"expected 2D array, got {arr.ndim}D")

    road_original, road_path = load_road_mask_for_domain(in_path, arr.shape)

    # Apply exactly the same orientation operations to DEM and road mask.
    arr = orient_ocean_right(arr, OCEAN_LOC)
    road_oriented = (
        orient_ocean_right(road_original, OCEAN_LOC)
        if road_original is not None else None
    )

    # Sanity check that OCEAN_LOC is actually correct.
    edge = max(1, min(5, arr.shape[1] // 20))
    left_q = np.nanpercentile(arr[:, :edge], 25)
    right_q = np.nanpercentile(arr[:, -edge:], 25)
    if right_q > left_q:
        print(f"[warn] {in_path.name}: after orienting, the LEFT edge is lower "
              f"(p25 left={left_q:.2f}, right={right_q:.2f}). Check OCEAN_LOC.")

    # Convert to ocean-first coordinates: index 0 is the ocean-side cell.
    raw = np.ascontiguousarray(arr[:, ::-1])
    road_raw = (
        np.ascontiguousarray(road_oriented[:, ::-1], dtype=bool)
        if road_oriented is not None else None
    )

    z = raw - MHW_M
    z[z < WATER_CLAMP_M] = SENTINEL_WATER_M

    # Trim DEM and road with the same cross-shore indices.
    c0, c1 = water_col_bounds(z, SENTINEL_WATER_M)
    z = z[:, c0:c1 + 1]
    road_mask = (
        np.ascontiguousarray(road_raw[:, c0:c1 + 1], dtype=bool)
        if road_raw is not None else None
    )

    n_along = z.shape[0]
    if n_along < ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} < {ALONG_COLS}; "
              f"trailing output cols remain sentinel.")
    elif n_along > ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} > {ALONG_COLS}; "
              f"only first {ALONG_COLS} profiles used.")

    above = z > BEACH_START_THR_M
    start_beach = np.where(above.any(axis=1), above.argmax(axis=1), -1)

    if road_mask is not None:
        road_seaward, road_landward, road_center = road_profile_positions(road_mask)
    else:
        road_seaward = np.full(n_along, np.nan)
        road_landward = np.full(n_along, np.nan)
        road_center = np.full(n_along, np.nan)

    return {
        "raw": raw,
        "z": z,
        "start_beach": start_beach,
        "c0": int(c0),
        "name": in_path.name,
        "road_raw": road_raw,
        "road_mask": road_mask,
        "road_seaward": road_seaward,
        "road_landward": road_landward,
        "road_center": road_center,
        "road_mask_path": road_path,
        "road_mask_name": road_path.name if road_path is not None else "",
    }

def masked_profiles(prof_arr: np.ndarray) -> np.ndarray:
    """NaN out sentinel water cells for plotting/statistics."""
    return np.where(prof_arr <= SENTINEL_WATER_M + 1e-9, np.nan, prof_arr)


def default_window(prof_arr: np.ndarray, start_beach: np.ndarray) -> tuple[int, int]:
    valid = start_beach[start_beach >= 0]
    base = int(np.median(valid)) if valid.size else 0
    return base, min(base + DEFAULT_WINDOW_PX, prof_arr.shape[1])


# ==============================================================================
# INTERACTIVE PICKER
# ==============================================================================

def _resize_line(values, n_along: int) -> np.ndarray:
    """Resize a saved/drawn alongshore line to n_along by normalized interpolation."""
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0:
        raise ValueError("empty boundary line")
    if arr.size == 1:
        return np.full(n_along, float(arr[0]))
    if arr.size == n_along:
        return arr.copy()
    old_x = np.linspace(0.0, 1.0, arr.size)
    new_x = np.linspace(0.0, 1.0, n_along)
    return np.interp(new_x, old_x, arr)


def _smooth_drawn_line(line: np.ndarray, width: int) -> np.ndarray:
    """Small edge-preserving moving average for noisy freehand mouse strokes."""
    width = int(width)
    if width <= 1:
        return line.copy()
    if width % 2 == 0:
        width += 1
    if line.size < width:
        return line.copy()
    pad = width // 2
    padded = np.pad(line, (pad, pad), mode="edge")
    kernel = np.ones(width, dtype=float) / width
    return np.convolve(padded, kernel, mode="valid")


def _stroke_to_line(points, n_along: int, n_cross: int) -> np.ndarray:
    """
    Convert a freehand mouse stroke into one y-value for every alongshore profile.

    The user may drag in either direction or backtrack. Duplicate x positions are
    averaged, then the line is interpolated/extrapolated across the full domain.
    """
    if not points:
        raise ValueError("no points were drawn")
    pts = np.asarray(points, dtype=float)
    pts = pts[np.isfinite(pts).all(axis=1)]
    if pts.size == 0:
        raise ValueError("drawn stroke contains no finite points")

    x = np.clip(pts[:, 0], 0, n_along - 1)
    y = np.clip(pts[:, 1], 0, n_cross)

    # Bin the freehand stroke to integer alongshore cells. Averaging duplicates
    # handles backtracking and a slow mouse without creating interpolation spikes.
    xi = np.rint(x).astype(int)
    unique_x = np.unique(xi)
    unique_y = np.array([np.mean(y[xi == xx]) for xx in unique_x], dtype=float)

    if unique_x.size == 1:
        line = np.full(n_along, unique_y[0], dtype=float)
    else:
        line = np.interp(
            np.arange(n_along, dtype=float), unique_x.astype(float), unique_y,
            left=unique_y[0], right=unique_y[-1],
        )
    line = _smooth_drawn_line(line, FREEHAND_SMOOTH_PX)
    return np.clip(line, 0, n_cross)


def normalize_window_lines(i0_line, i1_line, n_along: int, n_cross: int,
                           min_width: int = MIN_DRAWN_BAND_WIDTH_PX
                           ) -> tuple[np.ndarray, np.ndarray]:
    """
    Return integer, per-profile [i0, i1) limits with a valid local band width.

    This also provides backward compatibility: scalar i0/i1 values become flat
    lines, and saved lines with a different alongshore length are interpolated.
    """
    a = _resize_line(i0_line, n_along)
    b = _resize_line(i1_line, n_along)
    lo = np.minimum(a, b)
    hi = np.maximum(a, b)

    min_width = max(1, int(min_width))
    lo = np.clip(np.rint(lo), 0, max(0, n_cross - 1)).astype(int)
    hi = np.clip(np.rint(hi), 1, n_cross).astype(int)

    too_narrow = hi - lo < min_width
    hi[too_narrow] = np.minimum(n_cross, lo[too_narrow] + min_width)
    still_bad = hi <= lo
    lo[still_bad] = np.maximum(0, hi[still_bad] - 1)
    return lo, hi


def describe_window_lines(i0_line, i1_line) -> str:
    """Compact printable description of a locally varying search band."""
    a = np.asarray(i0_line, dtype=float)
    b = np.asarray(i1_line, dtype=float)
    return (f"seaward {a.min():.0f}..{a.max():.0f}, "
            f"landward {b.min():.0f}..{b.max():.0f}, "
            f"mean width {np.mean(b - a):.1f} cells")


def pick_window(stem: str, prof_arr: np.ndarray, start_beach: np.ndarray,
                init, road_mask=None) -> tuple[str, np.ndarray, np.ndarray]:
    """
    Draw a profile-varying dune-search band on the map.

    The two freehand boundaries may be horizontal, tilted, or curved. They are
    interpolated to every alongshore profile and returned as integer [i0, i1)
    arrays in ocean-first cross-shore coordinates.
    """
    n_along, n_cross = prof_arr.shape
    zm = masked_profiles(prof_arr)
    road_mask = validate_road_mask(road_mask, prof_arr.shape)

    if len(init) != 2:
        raise ValueError("init must contain seaward and landward boundaries")
    init_i0, init_i1 = normalize_window_lines(
        init[0], init[1], n_along, n_cross, MIN_DRAWN_BAND_WIDTH_PX
    )

    state = {
        "i0_line": init_i0.astype(float),
        "i1_line": init_i1.astype(float),
        "active": "seaward",
        "action": None,
        "dragging": False,
        "points": [],
        "history": [],
    }

    fig, (ax_map, ax_prof) = plt.subplots(
        1, 2, figsize=(14, 9), sharey=True,
        gridspec_kw={"width_ratios": [1.45, 1.0], "wspace": 0.06},
    )
    try:
        fig.canvas.manager.set_window_title(f"draw dune search band - {stem}")
    except Exception:
        pass

    finite = zm[np.isfinite(zm)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0

    im = ax_map.imshow(
        np.ma.masked_invalid(zm.T), aspect="auto", origin="lower",
        extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
        cmap="terrain", vmin=-1.0, vmax=max(vmax, 2.0),
    )
    add_road_plan_overlay(ax_map, road_mask, zorder=3.0)
    sb = np.where(start_beach >= 0, start_beach, np.nan)
    ax_map.plot(np.arange(n_along), sb, color="k", lw=1.2, label="beach start")
    ax_map.set_xlabel("alongshore cell")
    ax_map.set_ylabel("cross-shore cell  (0 = ocean, landward up)")

    y = np.arange(n_cross)
    ax_prof.plot(zm.T, y, color="0.75", lw=0.6)
    add_road_profile_envelope(ax_prof, road_mask)
    med = np.nanmedian(zm, axis=0)
    ax_prof.plot(med, y, color="k", lw=2.0, label="median profile")
    ax_prof.axvline(BEACH_START_THR_M, color="#1565C0", ls="--", lw=1.0,
                    label=f"beach thr {BEACH_START_THR_M} m")
    ax_prof.axvline(BERM_ELEV_NAVD_M - MHW_M, color="#B71C1C", ls=":", lw=1.2,
                    label=f"berm {BERM_ELEV_NAVD_M} m NAVD88")
    ax_prof.set_xlabel("elev (m MHW)")
    ax_prof.set_xlim(-1.2, max(vmax, 2.0) + 0.5)
    ax_prof.set_ylim(-0.5, n_cross - 0.5)
    ax_prof.legend(loc="upper right", fontsize=9, framealpha=0.9)
    plt.setp(ax_prof.get_yticklabels(), visible=False)

    try:
        fig.colorbar(im, ax=[ax_map, ax_prof], location="right",
                     fraction=0.035, pad=0.02, label="elev (m MHW)")
    except (TypeError, ValueError):
        fig.colorbar(im, ax=ax_prof, label="elev (m MHW)")

    xcells = np.arange(n_along)
    artists = {"band": None, "seaward": None, "landward": None,
               "stroke": None, "envelope": None}

    def normalized_for_display():
        return normalize_window_lines(
            state["i0_line"], state["i1_line"], n_along, n_cross,
            MIN_DRAWN_BAND_WIDTH_PX,
        )

    def redraw():
        i0_line, i1_line = normalized_for_display()
        state["i0_line"] = i0_line.astype(float)
        state["i1_line"] = i1_line.astype(float)

        for key in ("band", "seaward", "landward", "envelope"):
            if artists[key] is not None:
                artists[key].remove()
                artists[key] = None

        artists["band"] = ax_map.fill_between(
            xcells, i0_line, i1_line, color="#FF8C00", alpha=0.24,
            label="drawn search band", zorder=2,
        )
        artists["seaward"], = ax_map.plot(
            xcells, i0_line, color="#1565C0", lw=2.0,
            label="seaward boundary", zorder=4,
        )
        artists["landward"], = ax_map.plot(
            xcells, i1_line, color="#B71C1C", lw=2.0,
            label="landward boundary", zorder=4,
        )

        # The profile panel cannot show a different band for every alongshore
        # profile, so display only the total cross-shore envelope there.
        artists["envelope"] = ax_prof.axhspan(
            int(i0_line.min()), int(i1_line.max()), color="#FF8C00", alpha=0.12,
            zorder=0,
        )
        ax_map.legend(loc="upper right", fontsize=8, framealpha=0.9)

        active_label = ("SEAWARD (blue)" if state["active"] == "seaward"
                        else "LANDWARD (red)")
        fig.suptitle(
            f"{stem}   active = {active_label}\n"
            f"{describe_window_lines(i0_line, i1_line)}\n"
            "left-drag on MAP = draw active boundary | 1/2 = choose boundary | "
            "u = undo | enter = accept | r = reset | s = skip | q = quit",
            fontsize=11,
        )
        fig.canvas.draw_idle()

    def save_history():
        state["history"].append((
            state["i0_line"].copy(), state["i1_line"].copy(), state["active"]
        ))
        if len(state["history"]) > 20:
            state["history"].pop(0)

    def on_press(event):
        if event.inaxes is not ax_map or event.button != 1:
            return
        if event.xdata is None or event.ydata is None:
            return
        save_history()
        state["dragging"] = True
        state["points"] = [(event.xdata, event.ydata)]
        if artists["stroke"] is not None:
            artists["stroke"].remove()
        artists["stroke"], = ax_map.plot(
            [event.xdata], [event.ydata], color="white", lw=2.5,
            ls="--", zorder=6,
        )
        fig.canvas.draw_idle()

    def on_motion(event):
        if not state["dragging"] or event.inaxes is not ax_map:
            return
        if event.xdata is None or event.ydata is None:
            return
        state["points"].append((event.xdata, event.ydata))
        pts = np.asarray(state["points"])
        artists["stroke"].set_data(pts[:, 0], pts[:, 1])
        fig.canvas.draw_idle()

    def on_release(event):
        if not state["dragging"]:
            return
        state["dragging"] = False
        if event.inaxes is ax_map and event.xdata is not None and event.ydata is not None:
            state["points"].append((event.xdata, event.ydata))
        try:
            line = _stroke_to_line(state["points"], n_along, n_cross)
            key = "i0_line" if state["active"] == "seaward" else "i1_line"
            state[key] = line
            if AUTO_SWITCH_DRAW_BOUNDARY:
                state["active"] = (
                    "landward" if state["active"] == "seaward" else "seaward"
                )
        except ValueError as exc:
            print(f"[picker] {stem}: {exc}")
        finally:
            state["points"] = []
            if artists["stroke"] is not None:
                artists["stroke"].remove()
                artists["stroke"] = None
            redraw()

    def on_key(event):
        if event.key in ("enter", "return"):
            state["action"] = "accept"
            plt.close(fig)
        elif event.key == "1":
            state["active"] = "seaward"
            redraw()
        elif event.key == "2":
            state["active"] = "landward"
            redraw()
        elif event.key == "u":
            if state["history"]:
                i0_old, i1_old, active_old = state["history"].pop()
                state["i0_line"] = i0_old
                state["i1_line"] = i1_old
                state["active"] = active_old
                redraw()
        elif event.key == "r":
            save_history()
            state["i0_line"] = init_i0.astype(float)
            state["i1_line"] = init_i1.astype(float)
            state["active"] = "seaward"
            redraw()
        elif event.key == "s":
            state["action"] = "skip"
            plt.close(fig)
        elif event.key in ("q", "escape"):
            state["action"] = "quit"
            plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", on_press)
    fig.canvas.mpl_connect("motion_notify_event", on_motion)
    fig.canvas.mpl_connect("button_release_event", on_release)
    fig.canvas.mpl_connect("key_press_event", on_key)

    redraw()
    fig.canvas.draw_idle()
    fig.show()
    plt.show(block=True)  # blocking desktop GIS window

    if state["action"] is None:  # window closed with X
        state["action"] = "accept"
    i0_line, i1_line = normalized_for_display()
    return state["action"], i0_line, i1_line


# ==============================================================================
# EXTRACTION
# ==============================================================================

def find_dunes(prof_arr, start_beach, i0_line, i1_line):
    """Dune elevation/location per profile inside a locally varying drawn band."""
    n_along, n_cross = prof_arr.shape
    i0_line, i1_line = normalize_window_lines(
        i0_line, i1_line, n_along, n_cross, MIN_DRAWN_BAND_WIDTH_PX
    )
    dune_elev = np.full(n_along, np.nan)
    dune_loc = np.full(n_along, -1, dtype=int)

    for i in range(n_along):
        a = int(i0_line[i])
        if CLIP_WINDOW_TO_BEACH and start_beach[i] >= 0:
            a = max(a, int(start_beach[i]))
        b = min(int(i1_line[i]), n_cross)
        if b <= a:
            continue
        w = prof_arr[i, a:b]
        valid = w > SENTINEL_WATER_M + 1e-9
        if not valid.any():
            continue
        # Argmax is restricted to this profile's locally drawn search interval.
        k = int(np.argmax(np.where(valid, w, -np.inf)))
        dune_elev[i] = float(w[k])
        dune_loc[i] = a + k

    return dune_elev, dune_loc, i0_line, i1_line


def fill_missing_dune_locations(dune_loc: np.ndarray, n_cross: int) -> np.ndarray:
    """Interpolate missing crest locations alongshore for geometric re-indexing.

    The detected locations are retained wherever they exist. Missing profiles are
    filled by linear interpolation between neighboring valid profiles (constant
    extrapolation at either end). Locations are clipped so DUNE_ROWS real dune cells can be separated from at least
    one landward interior cell whenever possible.
    """
    loc = np.asarray(dune_loc, dtype=float).copy()
    valid = np.isfinite(loc) & (loc >= 0)
    if not valid.any():
        raise ValueError("cannot build local interior: no valid dune locations")

    x = np.arange(loc.size)
    loc[~valid] = np.interp(x[~valid], x[valid], loc[valid])
    loc = np.rint(loc).astype(int)

    max_loc = max(0, n_cross - DUNE_ROWS - 1)
    return np.clip(loc, 0, max_loc)



# ==============================================================================
# INTERACTIVE PRE-EXTRACTION PARTITION PICKER
# ==============================================================================

def normalize_interior_start_line(
        values,
        dune_loc_used: np.ndarray,
        n_along: int,
        n_cross: int,
        ) -> np.ndarray:
    """Return an integer interior-start cell for every alongshore profile.

    The detected dune crest is only a visual/reference feature. A manually chosen
    interior-start line may be either seaward or landward of the detected crest.
    The only enforced limits are the cross-shore bounds of the DEM.

    When ``values`` is None, the initial/default line is placed immediately
    landward of the detected dune. Pressing R in the picker returns to this default,
    but the user remains free to move the line seaward afterward.
    """
    dune_loc_used = np.asarray(dune_loc_used, dtype=int)
    default_line = np.where(
        dune_loc_used >= 0,
        dune_loc_used + DUNE_ROWS,
        0,
    )

    if values is None:
        line = default_line.astype(float)
    else:
        line = _resize_line(values, n_along)

    line = np.rint(line).astype(int)

    # Do not constrain the line relative to the detected dune crest.
    maximum = max(0, n_cross - 1)
    return np.clip(line, 0, maximum)


def pre_extraction_cross_section_qc(
        stem: str,
        prof_arr: np.ndarray,
        start_beach: np.ndarray,
        i0_line,
        i1_line,
        fig_dir: Path,
        initial_interior_line=None,
        road_mask=None,
        ) -> tuple[str, dict, np.ndarray]:
    """Interactively choose the profile-varying interior-start line.

    The layout is deliberately spacious:
      * plan view on the left;
      * full active-profile cross-section on the upper right;
      * clean transition zoom on the lower right;
      * a dedicated colorbar column;
      * status/instructions in a reserved header above the axes.
    """
    n_along, n_cross = prof_arr.shape
    n_fill = min(ALONG_COLS, n_along)
    road_mask = validate_road_mask(road_mask, prof_arr.shape)
    if n_fill <= 0 or n_cross <= 0:
        return "skip", {}, np.array([], dtype=int)

    dune_elev, dune_loc, i0_norm, i1_norm = find_dunes(
        prof_arr, start_beach, i0_line, i1_line
    )
    if not np.any(dune_loc[:n_fill] >= 0):
        print(f"[partition] {stem}: no dune crest found in any profile")
        return "skip", {}, np.array([], dtype=int)

    if FILL_MISSING_DUNE:
        dune_loc_used = fill_missing_dune_locations(dune_loc, n_cross)
    else:
        dune_loc_used = np.asarray(dune_loc, dtype=int).copy()

    default_interior = normalize_interior_start_line(
        None, dune_loc_used, n_along, n_cross
    )
    interior_line = normalize_interior_start_line(
        initial_interior_line, dune_loc_used, n_along, n_cross
    )

    zm = masked_profiles(prof_arr)
    x_along = np.arange(n_along)
    berm_mhw = BERM_ELEV_NAVD_M - MHW_M

    state = {
        "action": None,
        "active": min(n_fill // 2, n_along - 1),
        "interior": interior_line.copy(),
        "history": [],
        "plan_dragging": False,
        "plan_points": [],
        "cross_dragging": False,
    }

    # ------------------------------------------------------------------
    # CLEAN LAYOUT
    # Reserve the top ~14% for status and controls. A narrow middle
    # column holds the colorbar so it cannot squeeze or overlap the map.
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(17.5, 10.5))
    gs = fig.add_gridspec(
        2, 3,
        width_ratios=[1.20, 0.040, 1.0],
        height_ratios=[1.0, 1.0],
        left=0.055,
        right=0.975,
        bottom=0.075,
        top=0.845,
        wspace=0.24,
        hspace=0.38,
    )
    ax_map = fig.add_subplot(gs[:, 0])
    cax = fig.add_subplot(gs[:, 1])
    ax_full = fig.add_subplot(gs[0, 2])
    ax_zoom = fig.add_subplot(gs[1, 2])

    try:
        fig.canvas.manager.set_window_title(
            f"INTERACTIVE DUNE / INTERIOR PARTITION - {stem}"
        )
    except Exception:
        pass

    # Three separate header lines prevent the old suptitle from colliding
    # with the axes titles.
    status_artist = fig.text(
        0.5, 0.972, "",
        ha="center", va="top",
        fontsize=12.5, fontweight="semibold",
    )
    road_status_artist = fig.text(
        0.5, 0.942, "",
        ha="center", va="top",
        fontsize=10.2,
    )
    fig.text(
        0.5, 0.912,
        "PLAN: left-drag boundary; right-click selects profile    |    "
        "PROFILE: click/drag sets interior start    |    "
        "←/→ profile   U undo   R reset   Enter accept   S skip   Q quit",
        ha="center", va="top",
        fontsize=9.2,
    )

    finite = zm[np.isfinite(zm)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0
    im = ax_map.imshow(
        np.ma.masked_invalid(zm.T),
        aspect="auto",
        origin="lower",
        extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
        cmap="terrain",
        vmin=-1.0,
        vmax=max(vmax, 2.0),
        interpolation="nearest",
    )
    add_road_plan_overlay(ax_map, road_mask, zorder=3.0)
    ax_map.fill_between(
        x_along, i0_norm, i1_norm,
        color="#FF8C00", alpha=0.18,
        label="dune-search band",
    )
    ax_map.plot(
        x_along,
        np.where(start_beach >= 0, start_beach, np.nan),
        color="k", lw=1.0,
        label="beach start",
    )
    ax_map.plot(
        x_along,
        np.where(dune_loc >= 0, dune_loc, np.nan),
        color="#6A1B9A", lw=1.5, marker=".", ms=4,
        label="detected dune crest",
    )

    interior_artist, = ax_map.plot(
        x_along, state["interior"],
        color="#D81B60", lw=2.3, marker=".", ms=4,
        label="chosen interior start",
        zorder=7,
    )
    active_artist = ax_map.axvline(
        state["active"], color="white", lw=2.0, ls="--",
        label="active profile", zorder=8,
    )
    stroke_artist, = ax_map.plot(
        [], [], color="white", lw=2.8, ls="--", zorder=10,
    )

    ax_map.set_xlabel("alongshore cell")
    ax_map.set_ylabel("cross-shore cell  (0 = ocean; landward upward)")
    ax_map.set_title("Plan view — drag the interior-start boundary",
                     fontsize=12, pad=8)
    ax_map.legend(
        loc="upper right",
        fontsize=7.5,
        framealpha=0.92,
        borderpad=0.55,
        labelspacing=0.35,
        handlelength=2.2,
    )

    cb = fig.colorbar(im, cax=cax)
    cb.set_label("elevation relative to MHW (m)", fontsize=10)
    cb.ax.tick_params(labelsize=9)

    def save_history():
        state["history"].append(state["interior"].copy())
        if len(state["history"]) > 30:
            state["history"].pop(0)

    def active_default_start(profile_index: int) -> int:
        """Default reference: first cell immediately landward of the dune."""
        loc = int(dune_loc_used[profile_index])
        return max(0, loc + DUNE_ROWS) if loc >= 0 else 0

    def set_active_profile(profile_index: int):
        state["active"] = int(np.clip(profile_index, 0, n_fill - 1))
        redraw()

    def set_active_interior_from_distance(distance_m: float):
        if distance_m is None or not np.isfinite(distance_m):
            return
        i = state["active"]
        cell = int(np.rint(distance_m / CELL_SIZE_M))
        cell = int(np.clip(cell, 0, n_cross - 1))
        state["interior"][i] = cell
        redraw()

    def compact_full_legend(ax):
        """Keep only the useful legend entries and remove duplicates."""
        handles, labels = ax.get_legend_handles_labels()
        keep_exact = {
            ROAD_LABEL,
            "road center",
            "MHW",
            "drawn search band",
        }
        keep_prefixes = (
            "original profile",
            "beach threshold",
            "berm ",
            "beach start cell",
            "detected crest cell",
            "filled crest cell",
            "CHOSEN interior start",
        )

        selected_handles = []
        selected_labels = []
        seen = set()
        for handle, label in zip(handles, labels):
            useful = label in keep_exact or label.startswith(keep_prefixes)
            if useful and label not in seen:
                selected_handles.append(handle)
                selected_labels.append(label)
                seen.add(label)

        if selected_handles:
            ax.legend(
                selected_handles,
                selected_labels,
                loc="upper right",
                fontsize=6.9,
                framealpha=0.90,
                ncol=2,
                borderpad=0.45,
                columnspacing=0.8,
                labelspacing=0.30,
                handlelength=1.8,
            )

    def draw_profile_axis(ax, zoom: bool):
        ax.clear()
        i = state["active"]
        profile = np.asarray(prof_arr[i], dtype=float)
        cells = np.arange(n_cross)
        distance_m = cells * CELL_SIZE_M

        # Do not draw the long -3 m sentinel-water plateau. It was the main
        # reason the full cross-section looked visually heavy.
        profile_for_plot = np.where(
            np.isfinite(profile)
            & (profile > SENTINEL_WATER_M + 1e-9),
            profile,
            np.nan,
        )

        drawn_start = int(i0_norm[i])
        drawn_end = int(i1_norm[i])
        beach = int(start_beach[i]) if int(start_beach[i]) >= 0 else -1
        effective_start = max(drawn_start, beach) if (
            CLIP_WINDOW_TO_BEACH and beach >= 0
        ) else drawn_start
        detected = int(dune_loc[i])
        used = int(dune_loc_used[i])
        chosen = int(state["interior"][i])

        ax.plot(
            distance_m, profile_for_plot,
            color="0.15", lw=1.9,
            marker="o", ms=3.0,
            label=f"original profile {i}",
        )
        road_info = add_road_cross_section_overlay(
            ax, road_mask, i, profile,
            dune_cell=used, include_label=not zoom,
        )
        ax.axhline(0.0, color="#1565C0", lw=1.2, label="MHW")
        ax.axhline(
            BEACH_START_THR_M, color="#1565C0", ls="--", lw=1.0,
            label=f"beach threshold {BEACH_START_THR_M:.2f} m",
        )
        ax.axhline(
            berm_mhw, color="#B71C1C", ls=":", lw=1.2,
            label=f"berm {berm_mhw:.2f} m MHW",
        )
        ax.axvspan(
            drawn_start * CELL_SIZE_M,
            drawn_end * CELL_SIZE_M,
            color="#FF8C00", alpha=0.14,
            label="drawn search band",
        )
        if effective_start != drawn_start:
            ax.axvspan(
                effective_start * CELL_SIZE_M,
                drawn_end * CELL_SIZE_M,
                color="#FF8C00", alpha=0.24,
                label="effective searched band",
            )
        if beach >= 0:
            ax.axvline(
                beach * CELL_SIZE_M,
                color="k", ls="--", lw=1.0,
                label=f"beach start cell {beach}",
            )
        if 0 <= detected < n_cross:
            ax.scatter(
                detected * CELL_SIZE_M,
                profile[detected],
                s=105, color="#6A1B9A", zorder=8,
                label=f"detected crest cell {detected}",
            )
        if used != detected and 0 <= used < n_cross:
            ax.scatter(
                used * CELL_SIZE_M,
                profile[used],
                s=90, marker="D", color="#00897B", zorder=9,
                label=f"filled crest cell {used}",
            )

        default_start = active_default_start(i)
        if chosen > default_start:
            ax.axvspan(
                default_start * CELL_SIZE_M,
                chosen * CELL_SIZE_M,
                color="#D81B60", alpha=0.08,
                label="additional landward cells excluded",
            )
        elif chosen < default_start:
            ax.axvspan(
                chosen * CELL_SIZE_M,
                default_start * CELL_SIZE_M,
                color="#D81B60", alpha=0.08,
                label="interior begins seaward of dune",
            )
        ax.axvline(
            chosen * CELL_SIZE_M,
            color="#D81B60", lw=2.3,
            label=f"CHOSEN interior start cell {chosen}",
        )

        ax.grid(alpha=0.18)
        ax.tick_params(labelsize=9)
        ax.margins(x=0.01)

        if zoom:
            focus = [drawn_start, drawn_end, effective_start, used, chosen]
            if road_info["present"]:
                focus.extend([
                    int(road_info["seaward_cell"]),
                    int(road_info["landward_cell"]),
                ])

            context = max(1, int(CROSS_SECTION_CONTEXT_PX))
            z0 = max(0, min(focus) - context)
            z1 = min(n_cross, max(focus) + context + 1)

            ax.set_xlim(
                (z0 - 0.5) * CELL_SIZE_M,
                (z1 - 0.5) * CELL_SIZE_M,
            )

            visible = profile[z0:z1]
            good = visible[
                np.isfinite(visible)
                & (visible > SENTINEL_WATER_M + 1e-9)
            ]
            if good.size:
                spread = float(np.ptp(good))
                pad = max(0.32, 0.12 * (spread + 1.0))
                ax.set_ylim(
                    float(np.min(good)) - pad,
                    float(np.max(good)) + pad,
                )

            # Cell numbers belong on the axis, not on top of every point.
            # This keeps exact cell selection possible without label clutter.
            tick_step = 1 if (z1 - z0) <= 16 else 2
            tick_cells = np.arange(z0, z1, tick_step, dtype=int)
            ax.set_xticks(tick_cells * CELL_SIZE_M)
            ax.set_xticklabels(
                [str(c) for c in tick_cells],
                rotation=0 if len(tick_cells) <= 12 else 45,
                ha="center" if len(tick_cells) <= 12 else "right",
                fontsize=8,
            )
            ax.set_xlabel(
                f"cross-shore cell index  (1 cell = {CELL_SIZE_M:.0f} m)"
            )
            ax.set_ylabel("")
            ax.set_title(
                "Zoom — click or drag to choose the interior-start cell",
                fontsize=11.5, pad=8,
            )
        else:
            ax.set_xlabel("cross-shore distance from ocean-side edge (m)")
            ax.set_ylabel("elevation relative to MHW (m)")
            ax.set_title(
                f"Cross-section — active alongshore profile {i}",
                fontsize=11.5, pad=8,
            )
            compact_full_legend(ax)

    def redraw():
        interior_artist.set_data(x_along, state["interior"])
        active_artist.set_xdata([state["active"], state["active"]])

        draw_profile_axis(ax_full, zoom=False)
        draw_profile_axis(ax_zoom, zoom=True)

        i = state["active"]
        active_road = road_profile_summary(
            road_mask, i, int(dune_loc_used[i])
        )
        status_artist.set_text(
            f"{fig_title(stem)}   |   active profile {i}   |   "
            f"dune cell {int(dune_loc_used[i])}   |   "
            f"interior cell {int(state['interior'][i])}"
        )
        road_status_artist.set_text(active_road["note"])
        fig.canvas.draw_idle()

    def on_press(event):
        if event.inaxes is ax_map:
            if event.button == 3 and event.xdata is not None:
                set_active_profile(int(np.rint(event.xdata)))
                return
            if (
                event.button == 1
                and event.xdata is not None
                and event.ydata is not None
            ):
                state["plan_dragging"] = True
                state["plan_points"] = [(event.xdata, event.ydata)]
                stroke_artist.set_data([event.xdata], [event.ydata])
                fig.canvas.draw_idle()
                return

        if event.inaxes in (ax_full, ax_zoom) and event.button == 1:
            save_history()
            state["cross_dragging"] = True
            set_active_interior_from_distance(event.xdata)

    def on_motion(event):
        if state["plan_dragging"] and event.inaxes is ax_map:
            if event.xdata is None or event.ydata is None:
                return
            state["plan_points"].append((event.xdata, event.ydata))
            pts = np.asarray(state["plan_points"])
            stroke_artist.set_data(pts[:, 0], pts[:, 1])
            fig.canvas.draw_idle()
            return

        if state["cross_dragging"] and event.inaxes in (ax_full, ax_zoom):
            set_active_interior_from_distance(event.xdata)

    def on_release(event):
        if state["cross_dragging"]:
            state["cross_dragging"] = False

        if not state["plan_dragging"]:
            return

        state["plan_dragging"] = False
        if (
            event.inaxes is ax_map
            and event.xdata is not None
            and event.ydata is not None
        ):
            state["plan_points"].append((event.xdata, event.ydata))

        pts = np.asarray(state["plan_points"], dtype=float)
        state["plan_points"] = []
        stroke_artist.set_data([], [])

        if pts.size == 0:
            redraw()
            return

        # A short click selects a profile; a real stroke redraws the line.
        x_span = float(np.ptp(pts[:, 0])) if len(pts) > 1 else 0.0
        y_span = float(np.ptp(pts[:, 1])) if len(pts) > 1 else 0.0
        if len(pts) < 3 or (x_span < 0.75 and y_span < 0.75):
            set_active_profile(int(np.rint(pts[-1, 0])))
            return

        save_history()
        try:
            freehand = _stroke_to_line(pts.tolist(), n_along, n_cross)
            state["interior"] = normalize_interior_start_line(
                freehand, dune_loc_used, n_along, n_cross
            )
            state["active"] = int(
                np.clip(np.rint(pts[-1, 0]), 0, n_fill - 1)
            )
        except ValueError as exc:
            print(f"[partition] {stem}: {exc}")
        redraw()

    def on_key(event):
        key = (event.key or "").lower()
        if key in ("enter", "return", "a"):
            state["action"] = "accept"
            plt.close(fig)
        elif key == "u":
            if state["history"]:
                state["interior"] = state["history"].pop()
                redraw()
        elif key == "r":
            save_history()
            state["interior"] = default_interior.copy()
            redraw()
        elif key in ("left", "["):
            set_active_profile(state["active"] - 1)
        elif key in ("right", "]"):
            set_active_profile(state["active"] + 1)
        elif key == "s":
            state["action"] = "skip"
            plt.close(fig)
        elif key in ("q", "escape"):
            state["action"] = "quit"
            plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", on_press)
    fig.canvas.mpl_connect("motion_notify_event", on_motion)
    fig.canvas.mpl_connect("button_release_event", on_release)
    fig.canvas.mpl_connect("key_press_event", on_key)

    redraw()
    fig.canvas.draw()
    fig.canvas.flush_events()

    if REQUIRE_CROSS_SECTION_APPROVAL:
        print(
            f"[partition] opening plan view + cross-section for {stem}; "
            f"left-drag plan line or click cross-section anywhere"
        )
        plt.show(block=True)
        if state["action"] is None:
            state["action"] = "accept"
    elif SHOW_CROSS_SECTION_FIGS:
        plt.show()
        state["action"] = "accept"
    else:
        state["action"] = "accept"

    final_line = normalize_interior_start_line(
        state["interior"], dune_loc_used, n_along, n_cross
    )
    active = int(np.clip(state["active"], 0, n_fill - 1))
    profile = np.asarray(prof_arr[active], dtype=float)
    used_loc = int(dune_loc_used[active])
    used_elev = (
        float(profile[used_loc]) if 0 <= used_loc < n_cross else np.nan
    )

    active_road = road_profile_summary(road_mask, active, used_loc)
    qc = {
        "middle_profile_index": active,
        "middle_road_present": bool(active_road["present"]),
        "middle_road_seaward_cell": active_road["seaward_cell"],
        "middle_road_landward_cell": active_road["landward_cell"],
        "middle_road_center_cell": active_road["center_cell"],
        "middle_dune_to_road_m": active_road["dune_to_road_m"],
        "middle_beach_start": (
            int(start_beach[active]) if int(start_beach[active]) >= 0 else -1
        ),
        "middle_window_start": int(i0_norm[active]),
        "middle_window_end": int(i1_norm[active]),
        "middle_effective_start": max(
            int(i0_norm[active]),
            int(start_beach[active]) if (
                CLIP_WINDOW_TO_BEACH and int(start_beach[active]) >= 0
            ) else int(i0_norm[active]),
        ),
        "middle_effective_end": int(i1_norm[active]),
        "middle_dune_loc_detected": int(dune_loc[active]),
        "middle_dune_loc_used": used_loc,
        "middle_dune_elev_m_mhw": used_elev,
        "middle_dune_elev_m_navd88": (
            used_elev + MHW_M if np.isfinite(used_elev) else np.nan
        ),
        "middle_dune_height_above_berm_m": (
            max(used_elev - berm_mhw, MIN_DUNE_H_M)
            if np.isfinite(used_elev) else np.nan
        ),
        "middle_interior_start": int(final_line[active]),
        "chosen_interior_start_mean": float(np.mean(final_line[:n_fill])),
        "chosen_interior_start_min": int(np.min(final_line[:n_fill])),
        "chosen_interior_start_max": int(np.max(final_line[:n_fill])),
    }

    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_path = fig_dir / f"{fig_stem(stem)}_interactive_partition.png"
    qc["middle_cross_section_file"] = out_path.name
    if SAVE_CROSS_SECTION_FIGS and state["action"] == "accept":
        fig.savefig(out_path, dpi=175, facecolor="white")

    plt.close(fig)
    print(
        f"[partition] {stem}: action={state['action']}, "
        f"interior range={final_line[:n_fill].min()}.."
        f"{final_line[:n_fill].max()}"
    )
    return state["action"], qc, final_line

def build_local_dunes_and_interior(
        prof_arr: np.ndarray,
        dune_loc: np.ndarray,
        interior_start_override=None,
        ):
    """Separate DUNE_ROWS real local dune cells and align the remaining interior.

    Returns
    -------
    dune_m : ndarray, shape (ALONG_COLS, DUNE_ROWS)
        Dune heights above the berm in meters. This orientation matches
        Barrier3D's ``DuneDomain[:, q]`` convention: alongshore x dune width.
    topo_m : ndarray, shape (TOPO_ROWS, ALONG_COLS)
        Interior elevations relative to MHW. Row 0 is the manually selected
        source-GIS cell for each profile. It may be seaward or landward of the
        detected dune crest.
    dune_loc_used : ndarray, shape (n_along,)
        Detected locations with missing values interpolated when requested.
    interior_start : ndarray, shape (n_along,)
        Source-GIS cross-shore index used as interior row 0 in each profile.
    """
    n_along, n_cross = prof_arr.shape
    n_fill = min(ALONG_COLS, n_along)

    if n_cross <= DUNE_ROWS:
        raise ValueError(
            f"cross-shore width {n_cross} is too small for {DUNE_ROWS} dune rows"
        )

    if FILL_MISSING_DUNE:
        dune_loc_used = fill_missing_dune_locations(dune_loc, n_cross)
    else:
        dune_loc_used = np.asarray(dune_loc, dtype=int).copy()

    dune_m = np.full(
        (ALONG_COLS, DUNE_ROWS), SENTINEL_WATER_M, dtype=float
    )
    topo_m = np.full(
        (TOPO_ROWS, ALONG_COLS), SENTINEL_WATER_M, dtype=float
    )
    interior_start = np.full(n_along, -1, dtype=int)
    berm_mhw = BERM_ELEV_NAVD_M - MHW_M

    chosen_interior_start = normalize_interior_start_line(
        interior_start_override,
        dune_loc_used,
        n_along,
        n_cross,
    )

    for i in range(n_fill):
        loc = int(dune_loc_used[i])
        if loc < 0:
            continue

        # Reserve DUNE_ROWS consecutive REAL GIS cells for the dune field.
        # q=0 is the detected crest cell; q=1 is the next landward cell, etc.
        # Each cell keeps its own DEM-derived elevation.
        for q in range(DUNE_ROWS):
            src = loc + q
            if src >= n_cross:
                h = MIN_DUNE_H_M
            else:
                z = float(prof_arr[i, src])
                if not np.isfinite(z) or z <= SENTINEL_WATER_M + 1e-9:
                    h = MIN_DUNE_H_M if FILL_MISSING_DUNE else SENTINEL_WATER_M
                else:
                    h = max(z - berm_mhw, MIN_DUNE_H_M)
            dune_m[i, q] = h

        # Use the manually selected profile-varying interior boundary exactly.
        # It may be seaward or landward of the detected/reserved dune cell.
        start = int(chosen_interior_start[i])
        interior_start[i] = start
        use = prof_arr[i, start:]
        if use.size:
            rows = min(TOPO_ROWS, use.size)
            topo_m[:rows, i] = use[:rows]
            if use.size > TOPO_ROWS:
                print(
                    f"       [warn] profile {i} interior truncated: "
                    f"{use.size} -> {TOPO_ROWS} rows"
                )

    return dune_m, topo_m, dune_loc_used, interior_start


def extract_domain(stem, prof_arr, start_beach, i0_line, i1_line,
                   topo_dir, dune_dir, interior_start_override=None):
    n_along, n_cross = prof_arr.shape
    dune_elev, dune_loc, i0_line, i1_line = find_dunes(
        prof_arr, start_beach, i0_line, i1_line
    )

    n_found = int(np.isfinite(dune_elev).sum())
    if n_found == 0:
        print(f"[skip] {stem}: no dune found in drawn profile-varying band")
        return None

    n_fill = min(ALONG_COLS, n_along)
    missing = dune_loc[:n_fill] < 0
    if missing.any():
        idx = np.where(missing)[0].tolist()
        action = "interpolated alongshore" if FILL_MISSING_DUNE else "left empty"
        print(f"       [warn] {stem}: no detected crest at profiles {idx} -> {action}")

    dune_m, topo_m, dune_loc_used, interior_start_line = \
        build_local_dunes_and_interior(
            prof_arr,
            dune_loc,
            interior_start_override=interior_start_override,
        )

    if TRIM_INTERIOR_ROWS:
        topo_m = trim_landward_water_rows(topo_m, SENTINEL_WATER_M)
    else:
        # Keep the complete fixed-size CASCADE interior, including sentinel water.
        expected_shape = (TOPO_ROWS, ALONG_COLS)
        if topo_m.shape != expected_shape:
            padded = np.full(expected_shape, SENTINEL_WATER_M, dtype=float)
            rows = min(TOPO_ROWS, topo_m.shape[0])
            cols = min(ALONG_COLS, topo_m.shape[1])
            padded[:rows, :cols] = topo_m[:rows, :cols]
            topo_m = padded
        print(
            f"       [water preserved] interior retained at {topo_m.shape}; "
            f"sentinel water = {SENTINEL_WATER_M:.1f} m MHW"
        )

    topo_dm = topo_m * 0.1

    # --------------------------------------------------------------------------
    # TWO REAL DUNE ROWS
    # --------------------------------------------------------------------------
    # Keep the physically extracted dune field as a 2-D array internally:
    #     (alongshore cells, cross-shore dune rows) = (50, 2)
    # This is the orientation used by the QC/plotting routines below.
    dune_dm_2d = dune_m * 0.1

    expected_internal_shape = (ALONG_COLS, DUNE_ROWS)
    if dune_dm_2d.shape != expected_internal_shape:
        raise ValueError(
            f"{stem}: internal dune array shape {dune_dm_2d.shape} does not "
            f"match expected {expected_internal_shape}"
        )

    # Barrier3D's native multiple-row input format is one packed 1-D array:
    #     first  ALONG_COLS values -> dune row 0
    #     second ALONG_COLS values -> dune row 1
    # and so on for any additional rows.
    # Transpose first so each cross-shore dune row is stored as one contiguous
    # alongshore block before flattening.
    dune_dm_save = dune_dm_2d.T.reshape(-1)

    expected_saved_shape = (ALONG_COLS * DUNE_ROWS,)
    if dune_dm_save.shape != expected_saved_shape:
        raise ValueError(
            f"{stem}: saved dune array shape {dune_dm_save.shape} does not "
            f"match expected {expected_saved_shape}"
        )

    # Safety check: unpack the saved representation and verify that it exactly
    # reconstructs the internal (alongshore, dune-row) array.
    dune_check = dune_dm_save.reshape(DUNE_ROWS, ALONG_COLS).T
    if not np.allclose(dune_check, dune_dm_2d, equal_nan=True):
        raise RuntimeError(f"{stem}: dune packing/unpacking validation failed")

    print(
        f"Two-row dune input: internal={dune_dm_2d.shape}, "
        f"saved={dune_dm_save.shape}"
    )

    # Row 0 begins at the detected crest cell; the following columns are
    # consecutive landward dune cells from the DEM.
    crest_h = dune_m[:n_fill, 0]
    valid_crest_h = crest_h[crest_h > SENTINEL_WATER_M + 1e-9]
    start_valid = interior_start_line[:n_fill]
    start_valid = start_valid[start_valid >= 0]
    start_island = float(np.mean(start_valid)) if start_valid.size else None

    topo_out = Path(topo_dir) / f"{stem}_topography_{TAG}.npy"
    dune_out = Path(dune_dir) / f"{stem}_dune_{TAG}.npy"
    np.save(topo_out, topo_dm)
    np.save(dune_out, dune_dm_save)

    band_desc = describe_window_lines(i0_line, i1_line)
    print(f"[ok] {stem}: {band_desc}, {n_found}/{n_along} detected crests, "
          f"dunes internal {dune_dm_2d.shape}, saved {dune_dm_save.shape}, "
          f"interior {topo_dm.shape}, mean crest h = "
          f"{np.mean(valid_crest_h):.2f} m")

    return {
        "stem": stem,
        # Mean scalar values are retained for old summary/report code, while the
        # complete profile-varying geometry is carried in i0_line/i1_line.
        "i0": float(np.mean(i0_line)),
        "i1": float(np.mean(i1_line)),
        "i0_line": i0_line,
        "i1_line": i1_line,
        "i0_min": int(np.min(i0_line)),
        "i0_max": int(np.max(i0_line)),
        "i1_min": int(np.min(i1_line)),
        "i1_max": int(np.max(i1_line)),
        "window_width_mean": float(np.mean(i1_line - i0_line)),
        "window_width_min": int(np.min(i1_line - i0_line)),
        "window_width_max": int(np.max(i1_line - i0_line)),
        "n_found": n_found, "n_along": n_along,
        "mean_dune_h_m": float(np.mean(valid_crest_h)),
        "min_dune_h_m": float(np.min(valid_crest_h)),
        "max_dune_h_m": float(np.max(valid_crest_h)),
        "dune_rows": DUNE_ROWS,
        "dune_shape": tuple(dune_dm_2d.shape),
        "dune_saved_shape": tuple(dune_dm_save.shape),
        "n_filled": int(missing.sum()),
        "interior_rows": int(topo_dm.shape[0]),
        "interior_cols": int(topo_dm.shape[1]),
        "mean_interior_elev_m": float(np.mean(topo_m[topo_m > SENTINEL_WATER_M + 1e-9]))
        if np.any(topo_m > SENTINEL_WATER_M + 1e-9) else np.nan,
        # start_island is retained as the mean for old summary code.
        "start_island": start_island,
        "interior_start_line": interior_start_line,
        "interior_start_min": int(np.min(start_valid)) if start_valid.size else -1,
        "interior_start_max": int(np.max(start_valid)) if start_valid.size else -1,
        "dune_elev": dune_elev,
        "dune_loc": dune_loc,
        "dune_loc_used": dune_loc_used,
        "dune_h": crest_h,
        "topo_dm": topo_dm, "dune_dm": dune_dm_2d,
        "topo_file": topo_out.name, "dune_file": dune_out.name,
    }


# ==============================================================================
# QC FIGURE
# ==============================================================================

def qc_figure(stem, prof_arr, start_beach, res, fig_dir: Path, road_mask=None):
    """QC figure showing the exact tilted/curved search band and detected crest."""
    n_along, n_cross = prof_arr.shape
    zm = masked_profiles(prof_arr)
    road_mask = validate_road_mask(road_mask, prof_arr.shape)
    i0_line = np.asarray(res["i0_line"])
    i1_line = np.asarray(res["i1_line"])
    dl = np.where(res["dune_loc"] >= 0, res["dune_loc"], np.nan)
    x = np.arange(n_along)

    fig = plt.figure(figsize=(12.5, 9.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1.0], height_ratios=[1.0, 0.35],
                          wspace=0.06, hspace=0.12)
    ax_map = fig.add_subplot(gs[0, 0])
    ax_prof = fig.add_subplot(gs[0, 1], sharey=ax_map)
    ax_h = fig.add_subplot(gs[1, 0], sharex=ax_map)

    finite = zm[np.isfinite(zm)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0

    im = ax_map.imshow(np.ma.masked_invalid(zm.T), aspect="auto", origin="lower",
                       extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
                       cmap="terrain", vmin=-1.0, vmax=max(vmax, 2.0))
    add_road_plan_overlay(ax_map, road_mask, zorder=3.0)
    ax_map.fill_between(x, i0_line, i1_line, color="#FF8C00", alpha=0.22,
                        label="drawn search band", zorder=2)
    ax_map.plot(x, i0_line, color="#1565C0", lw=1.5, label="seaward boundary")
    ax_map.plot(x, i1_line, color="#B71C1C", lw=1.5, label="landward boundary")
    ax_map.plot(x, np.where(start_beach >= 0, start_beach, np.nan),
                color="k", lw=1.0, label="beach start")
    ax_map.plot(x, dl, color="#6A1B9A", lw=1.1, marker=".", ms=4,
                label="detected dune crest")
    interior_start = np.asarray(res["interior_start_line"], dtype=float)
    interior_start[interior_start < 0] = np.nan
    ax_map.plot(x, interior_start, color="w", lw=1.4, ls="--",
                label=f"interior start after {DUNE_ROWS} dune cells")
    ax_map.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_map.legend(loc="upper right", fontsize=8, framealpha=0.9)
    plt.setp(ax_map.get_xticklabels(), visible=False)

    y = np.arange(n_cross)
    ax_prof.plot(zm.T, y, color="0.75", lw=0.5)
    add_road_profile_envelope(ax_prof, road_mask)
    ax_prof.plot(np.nanmedian(zm, axis=0), y, color="k", lw=2.0)
    ax_prof.axhspan(i0_line.min(), i1_line.max(), color="#FF8C00", alpha=0.12,
                    zorder=0, label="total band envelope")
    ax_prof.axvline(BERM_ELEV_NAVD_M - MHW_M, color="#B71C1C", ls=":", lw=1.2)
    ax_prof.plot(res["dune_elev"], dl, lw=0, marker=".", ms=4, color="#6A1B9A")
    ax_prof.set_xlabel("elev (m MHW)")
    ax_prof.set_xlim(-1.2, max(vmax, 2.0) + 0.5)
    ax_prof.set_ylim(-0.5, n_cross - 0.5)
    plt.setp(ax_prof.get_yticklabels(), visible=False)

    dune_x = x[:len(res["dune_h"])]
    ax_h.plot(dune_x, res["dune_h"], color="#FF8C00", lw=1.5,
              marker="o", ms=3, label="dune crest row")
    ax_h.axhline(MIN_DUNE_H_M, color="0.5", ls="--", lw=1.0)
    ax_h.set_xlabel("alongshore cell")
    ax_h.set_ylabel("dune height\nabove berm (m)")
    ax_h.set_xlim(-0.5, n_along - 0.5)
    ax_h.legend(loc="upper right", fontsize=8, framealpha=0.9)

    try:
        fig.colorbar(im, ax=[ax_map, ax_prof], location="right",
                     fraction=0.035, pad=0.02, label="elev (m MHW)")
    except (TypeError, ValueError):
        fig.colorbar(im, ax=ax_prof, label="elev (m MHW)")

    fig.suptitle(
        f"{fig_title(stem)}  |  {describe_window_lines(i0_line, i1_line)}  |  "
        f"{res['n_found']}/{n_along} dunes found", fontsize=12
    )

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"{fig_stem(stem)}_qc.png", dpi=130,
                bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# MIDDLE-DOMAIN CROSS-SECTION ANALYSIS
# ==============================================================================

def middle_cross_section_figure(stem, prof_arr, start_beach, res, fig_dir: Path, road_mask=None):
    """
    Plot the cross-shore profile from the middle alongshore row of a domain.

    The upper panel shows the complete ocean-to-land profile. The lower panel
    zooms around the locally drawn dune-search band so the selected dune crest
    can be checked against neighboring cells.

    Coordinates are ocean-first:
        cross-shore cell 0 = ocean side
        increasing cell number = landward

    The function also adds middle-profile values to ``res`` so they are written
    to the settings CSV/XLSX.
    """
    n_along, n_cross = prof_arr.shape
    road_mask = validate_road_mask(road_mask, prof_arr.shape)
    n_fill = min(ALONG_COLS, n_along)
    if n_fill <= 0 or n_cross <= 0:
        print(f"[cross section] {stem}: empty profile array; skipped")
        return

    # For 50 profiles this selects index 25, the center-right member of 24/25.
    mid = n_fill // 2
    profile = np.asarray(prof_arr[mid, :], dtype=float)
    x_cell = np.arange(n_cross, dtype=int)
    x_m = x_cell.astype(float) * CELL_SIZE_M

    i0_line = np.asarray(res["i0_line"], dtype=int)
    i1_line = np.asarray(res["i1_line"], dtype=int)
    dune_loc_detected_line = np.asarray(res["dune_loc"], dtype=int)
    dune_loc_used_line = np.asarray(res["dune_loc_used"], dtype=int)
    interior_start_line = np.asarray(res["interior_start_line"], dtype=int)

    drawn_start = int(np.clip(i0_line[mid], 0, max(0, n_cross - 1)))
    drawn_end = int(np.clip(i1_line[mid], drawn_start + 1, n_cross))
    beach = int(start_beach[mid]) if int(start_beach[mid]) >= 0 else -1

    effective_start = drawn_start
    if CLIP_WINDOW_TO_BEACH and beach >= 0:
        effective_start = max(effective_start, beach)
    effective_end = drawn_end

    detected_loc = int(dune_loc_detected_line[mid])
    used_loc = int(dune_loc_used_line[mid])
    interior_start = int(interior_start_line[mid])

    detected_elev = (
        float(profile[detected_loc])
        if 0 <= detected_loc < n_cross else np.nan
    )
    used_elev = (
        float(profile[used_loc])
        if 0 <= used_loc < n_cross else np.nan
    )
    berm_mhw = BERM_ELEV_NAVD_M - MHW_M

    # Save the middle-profile diagnostics into the result dictionary.
    res["middle_profile_index"] = mid
    res["middle_beach_start"] = beach
    res["middle_window_start"] = drawn_start
    res["middle_window_end"] = drawn_end
    res["middle_effective_start"] = effective_start
    res["middle_effective_end"] = effective_end
    res["middle_dune_loc_detected"] = detected_loc
    res["middle_dune_loc_used"] = used_loc
    res["middle_dune_elev_m_mhw"] = used_elev
    res["middle_dune_elev_m_navd88"] = (
        used_elev + MHW_M if np.isfinite(used_elev) else np.nan
    )
    res["middle_dune_height_above_berm_m"] = (
        max(used_elev - berm_mhw, MIN_DUNE_H_M)
        if np.isfinite(used_elev) else np.nan
    )
    res["middle_interior_start"] = interior_start

    fig, (ax_full, ax_zoom) = plt.subplots(
        2, 1, figsize=(14, 9),
        gridspec_kw={"height_ratios": [1.25, 1.0], "hspace": 0.30},
    )

    def add_reference_lines(ax, include_labels=True):
        ax.axhline(
            0.0, color="#1565C0", lw=1.4,
            label="MHW = 0 m" if include_labels else None,
        )
        ax.axhline(
            BEACH_START_THR_M, color="#1565C0", ls="--", lw=1.2,
            label=(
                f"beach threshold = {BEACH_START_THR_M:.2f} m MHW"
                if include_labels else None
            ),
        )
        ax.axhline(
            berm_mhw, color="#B71C1C", ls=":", lw=1.4,
            label=(
                f"berm = {berm_mhw:.2f} m MHW "
                f"({BERM_ELEV_NAVD_M:.2f} m NAVD88)"
                if include_labels else None
            ),
        )

    def add_geometry(ax, include_labels=True):
        add_road_cross_section_overlay(
            ax, road_mask, mid, profile, dune_cell=used_loc,
            include_label=include_labels,
        )
        # Drawn band and the actually searched portion after beach clipping.
        ax.axvspan(
            drawn_start * CELL_SIZE_M,
            drawn_end * CELL_SIZE_M,
            color="#FF8C00", alpha=0.18,
            label="drawn search band" if include_labels else None,
        )
        if effective_start != drawn_start:
            ax.axvspan(
                effective_start * CELL_SIZE_M,
                effective_end * CELL_SIZE_M,
                color="#FF8C00", alpha=0.30,
                label=(
                    "effective band after beach clipping"
                    if include_labels else None
                ),
            )

        if beach >= 0:
            ax.axvline(
                beach * CELL_SIZE_M, color="k", ls="--", lw=1.2,
                label=f"beach start: cell {beach}" if include_labels else None,
            )

        if 0 <= detected_loc < n_cross:
            ax.scatter(
                detected_loc * CELL_SIZE_M, detected_elev,
                s=90, color="#6A1B9A", zorder=8,
                label=(
                    f"detected crest: cell {detected_loc}, "
                    f"z={detected_elev:.2f} m MHW"
                    if include_labels else None
                ),
            )
            ax.axvline(
                detected_loc * CELL_SIZE_M, color="#6A1B9A",
                ls="-.", lw=1.3,
            )

        # If a crest was missing, the extraction may use an interpolated location.
        if used_loc != detected_loc and 0 <= used_loc < n_cross:
            ax.scatter(
                used_loc * CELL_SIZE_M, used_elev,
                s=90, marker="D", color="#00897B", zorder=9,
                label=(
                    f"interpolated/used crest: cell {used_loc}, "
                    f"z={used_elev:.2f} m MHW"
                    if include_labels else None
                ),
            )
            ax.axvline(
                used_loc * CELL_SIZE_M, color="#00897B",
                ls="-.", lw=1.3,
            )

        if 0 <= interior_start < n_cross:
            ax.axvline(
                interior_start * CELL_SIZE_M,
                color="#D81B60", ls=":", lw=2.0,
                label=(
                    f"interior starts: cell {interior_start}"
                    if include_labels else None
                ),
            )

    # Full ocean-to-land cross-section.
    ax_full.plot(
        x_m, profile, color="0.15", lw=2.0,
        marker="o", ms=2.5, label=f"profile {mid}",
    )
    add_reference_lines(ax_full, include_labels=True)
    add_geometry(ax_full, include_labels=True)
    ax_full.set_xlabel("cross-shore distance from ocean-side edge (m)")
    ax_full.set_ylabel("elevation relative to MHW (m)")
    ax_full.set_title(
        f"{fig_title(stem)} — middle alongshore cross-section (profile {mid})"
    )
    ax_full.grid(alpha=0.25)
    ax_full.legend(loc="best", fontsize=8.5, framealpha=0.92, ncol=2)

    # Zoomed dune-zone profile, with each 10 m GIS cell visible.
    context = max(0, int(CROSS_SECTION_CONTEXT_PX))
    useful_locations = [drawn_start, effective_start, drawn_end, effective_end]
    if used_loc >= 0:
        useful_locations.append(used_loc)
    zoom_start = max(0, min(useful_locations) - context)
    zoom_end = min(
        n_cross,
        max(useful_locations + [max(0, used_loc) + DUNE_ROWS + 1]) + context,
    )
    if zoom_end <= zoom_start:
        zoom_start, zoom_end = 0, n_cross

    ax_zoom.plot(
        x_m, profile, color="0.15", lw=2.0,
        marker="o", ms=5.0, label="GIS cell elevations",
    )
    add_reference_lines(ax_zoom, include_labels=False)
    add_geometry(ax_zoom, include_labels=False)
    ax_zoom.set_xlim(
        (zoom_start - 0.5) * CELL_SIZE_M,
        (zoom_end - 0.5) * CELL_SIZE_M,
    )

    zoom_values = profile[zoom_start:zoom_end]
    finite_zoom = zoom_values[np.isfinite(zoom_values)]
    finite_zoom = finite_zoom[finite_zoom > SENTINEL_WATER_M + 1e-9]
    if finite_zoom.size:
        spread = float(finite_zoom.max()) - float(finite_zoom.min())
        pad = max(0.35, 0.12 * (spread + 1.0))
        ax_zoom.set_ylim(
            float(finite_zoom.min()) - pad,
            float(finite_zoom.max()) + pad,
        )

    # Label the ocean-first cell number over each visible elevation point.
    for cell in range(zoom_start, zoom_end):
        z = profile[cell]
        if np.isfinite(z) and z > SENTINEL_WATER_M + 1e-9:
            ax_zoom.annotate(
                str(cell),
                xy=(cell * CELL_SIZE_M, z),
                xytext=(0, 7), textcoords="offset points",
                ha="center", va="bottom", fontsize=7, rotation=45,
            )

    ax_zoom.set_xlabel("cross-shore distance from ocean-side edge (m)")
    ax_zoom.set_ylabel("elevation relative to MHW (m)")
    ax_zoom.set_title(
        "Dune-zone zoom — labels are ocean-first cross-shore cell indices"
    )
    ax_zoom.grid(alpha=0.25)

    status = "detected"
    if detected_loc < 0 and used_loc >= 0:
        status = "interpolated from neighboring alongshore profiles"

    if np.isfinite(used_elev):
        elev_text = (
            f"{used_elev:.2f} m MHW / "
            f"{used_elev + MHW_M:.2f} m NAVD88"
        )
        dune_h_text = f"{res['middle_dune_height_above_berm_m']:.2f} m"
    else:
        elev_text = "not available"
        dune_h_text = "not available"

    middle_road = road_profile_summary(road_mask, mid, used_loc)
    analysis_text = (
        f"Profile: {mid}\n"
        f"Drawn band: [{drawn_start}, {drawn_end}) cells "
        f"= {drawn_start * CELL_SIZE_M:.0f}–{drawn_end * CELL_SIZE_M:.0f} m\n"
        f"Effective band: [{effective_start}, {effective_end})\n"
        f"Beach start: {beach if beach >= 0 else 'none'}\n"
        f"Dune crest used: cell {used_loc} ({status})\n"
        f"Crest elevation: {elev_text}\n"
        f"Dune height above berm: {dune_h_text}\n"
        f"Interior starts: cell {interior_start}\n"
        f"Road: {middle_road['note']}"
    )
    ax_zoom.text(
        0.995, 0.98, analysis_text,
        transform=ax_zoom.transAxes,
        ha="right", va="top", fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "white",
            "edgecolor": "0.5",
            "alpha": 0.92,
        },
    )

    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_path = fig_dir / f"{fig_stem(stem)}_middle_cross_section.png"
    res["middle_cross_section_file"] = out_path.name
    fig.savefig(out_path, dpi=160, bbox_inches="tight")

    if SHOW_CROSS_SECTION_FIGS:
        plt.show()
    plt.close(fig)

    print(
        f"[cross section] {stem}: profile {mid}, dune cell {used_loc}, "
        f"z={used_elev:.2f} m MHW, interior starts {interior_start} "
        f"-> {out_path.name}"
    )


# ==============================================================================
# GIS vs PROCESSED COMPARISON FIGURE
# ==============================================================================

def comparison_figure(stem, dom, res, fig_dir: Path):
    """Raw GIS versus saved CASCADE arrays, including the drawn local search band."""
    raw = dom["raw"]
    road_raw = dom.get("road_raw")
    road_mask = dom.get("road_mask")
    c0 = dom["c0"]
    n_along, n_raw = raw.shape
    x = np.arange(n_along)

    i0_raw = np.asarray(res["i0_line"]) + c0
    i1_raw = np.asarray(res["i1_line"]) + c0
    dl_raw = np.where(res["dune_loc"] >= 0, res["dune_loc"] + c0, np.nan)
    sb_raw = np.where(dom["start_beach"] >= 0, dom["start_beach"] + c0, np.nan)

    topo_m = res["topo_dm"] * 10.0
    dune_m = res["dune_dm"] * 10.0
    topo_disp = np.where(topo_m <= SENTINEL_WATER_M + 1e-9, np.nan, topo_m)
    n_int_rows, n_int_cols = topo_disp.shape

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[1.0, 0.14],
                          wspace=0.16, hspace=0.06)
    ax_raw = fig.add_subplot(gs[:, 0])
    ax_int = fig.add_subplot(gs[0, 1])
    ax_dun = fig.add_subplot(gs[1, 1])

    land = raw[raw > MHW_M]
    hi_navd = float(np.percentile(land, 99)) if land.size else MHW_M + 4.0
    vmin_r = MHW_M - 1.0
    vmax_r = max(hi_navd, MHW_M + 2.0)
    im_r = ax_raw.imshow(raw.T, aspect="auto", origin="lower",
                         extent=[-0.5, n_along - 0.5, -0.5, n_raw - 0.5],
                         cmap="terrain", vmin=vmin_r, vmax=vmax_r)
    try:
        fig.colorbar(im_r, ax=ax_raw, location="left", pad=0.02, fraction=0.04,
                     label="elev (m NAVD88)")
    except (TypeError, ValueError):
        fig.colorbar(im_r, ax=ax_raw, pad=0.01, fraction=0.04,
                     label="elev (m NAVD88)")
    ax_raw.contour(np.arange(n_along), np.arange(n_raw), raw.T, levels=[MHW_M],
                   colors="#1565C0", linewidths=1.0)
    add_road_plan_overlay(ax_raw, road_raw, zorder=3.0)
    ax_raw.fill_between(x, i0_raw, i1_raw, color="#FF8C00", alpha=0.22,
                        label="drawn dune-search band", zorder=2)
    ax_raw.plot(x, i0_raw, color="#1565C0", lw=1.5, label="seaward boundary")
    ax_raw.plot(x, i1_raw, color="#B71C1C", lw=1.5, label="landward boundary")
    ax_raw.plot(x, sb_raw, color="k", lw=1.0, label="beach start")
    ax_raw.plot(x, dl_raw, color="#6A1B9A", lw=1.1, marker=".", ms=4,
                label="detected dune crest")
    interior_start_raw = np.asarray(res["interior_start_line"], dtype=float) + c0
    interior_start_raw[np.asarray(res["interior_start_line"]) < 0] = np.nan
    ax_raw.plot(x, interior_start_raw, color="w", lw=1.4, ls="--",
                label=f"interior start after {DUNE_ROWS} dune cells")
    ax_raw.axhline(c0 - 0.5, color="0.3", lw=1.0, ls="-.", label="water trim edge")
    ax_raw.set_xlabel("alongshore cell")
    ax_raw.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_raw.set_title(f"RAW GIS input   {raw.shape}   m NAVD88\n{dom['name']}",
                     fontsize=11)
    ax_raw.legend(loc="upper right", fontsize=8, framealpha=0.9)

    im_i = ax_int.imshow(np.ma.masked_invalid(topo_disp), aspect="auto",
                         origin="lower",
                         extent=[-0.5, n_int_cols - 0.5, -0.5, n_int_rows - 0.5],
                         cmap="terrain", vmin=vmin_r - MHW_M, vmax=vmax_r - MHW_M)
    fig.colorbar(im_i, ax=ax_int, pad=0.01, fraction=0.04, label="elev (m MHW)")
    processed_road = processed_interior_road_mask(
        road_mask, np.asarray(res["interior_start_line"], dtype=int),
        n_int_rows, n_int_cols,
    )
    add_processed_road_overlay(ax_int, processed_road)
    if np.any(processed_road):
        ax_int.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax_int.set_ylabel("interior row\n(0 = ocean side)")
    ax_int.set_title(f"PROCESSED CASCADE input\ninterior {res['topo_file']}  "
                     f"{res['topo_dm'].shape}  (dam on disk)", fontsize=11)
    plt.setp(ax_int.get_xticklabels(), visible=False)

    dune_disp = np.where(
        dune_m <= SENTINEL_WATER_M + 1e-9, np.nan, dune_m
    )
    # Internal orientation is (alongshore, dune row). Ensure a 2-D display even
    # when only one dune row is used.
    dune_view = np.atleast_2d(dune_disp.T)
    im_d = ax_dun.imshow(np.ma.masked_invalid(dune_view), aspect="auto",
                         origin="lower",
                         extent=[-0.5, ALONG_COLS - 0.5, -0.5, DUNE_ROWS - 0.5],
                         cmap="YlOrBr", vmin=0.0,
                         vmax=max(float(np.nanmax(dune_view)), 0.5))
    fig.colorbar(im_d, ax=ax_dun, pad=0.01, fraction=0.04,
                 label="dune h above berm (m)")
    ax_dun.set_yticks(np.arange(DUNE_ROWS))
    ax_dun.set_yticklabels([f"dune {q + 1}" for q in range(DUNE_ROWS)])
    ax_dun.set_xlabel("alongshore cell")

    fig.suptitle(
        f"{fig_title(stem)}  —  {RUN_NAME}   |   "
        f"{describe_window_lines(res['i0_line'], res['i1_line'])}   |   "
        f"{res['n_found']}/{res['n_along']} dunes   |   "
        f"mean crest h {res['mean_dune_h_m']:.2f} m   |   "
        f"{DUNE_ROWS} dune rows   |   local interior starts "
        f"{res['interior_start_min']}-{res['interior_start_max']}",
        fontsize=12,
    )

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"{fig_stem(stem)}_gis_vs_processed.png", dpi=130,
                bbox_inches="tight")
    plt.close(fig)


# ==============================================================================
# SETTINGS SHEET & JSON I/O
# ==============================================================================

def global_config() -> dict:
    """Every global knob, for the config sheet / manifest / provenance."""
    return {
        "run name": RUN_NAME,
        "version": VERSION,
        "first real domain": FIRST_REAL_DOMAIN_ID,
        "last real domain": LAST_REAL_DOMAIN_ID,
        "number real domains": NUM_REAL_DOMAINS,
        "run_time": datetime.now().isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "root GIS domains": LOAD_PATH.name,
        "run folder": str(RUN_DIR),
        "load_path_full": str(LOAD_PATH),
        "topo_path_full": str(TOPO_SAVE_PATH),
        "dune_path_full": str(DUNE_SAVE_PATH),
        "window_json": str(WINDOW_JSON),
        "SHOW_ROAD": SHOW_ROAD,
        "REQUIRE_ROAD_MASKS": REQUIRE_ROAD_MASKS,
        "ROAD_MASK_FOLDER": str(ROAD_MASK_FOLDER),
        "ROAD_MASK_PATTERNS": ROAD_MASK_PATTERNS,
        "ROAD_LABEL": ROAD_LABEL,
        "ROAD_PLAN_ALPHA": ROAD_PLAN_ALPHA,
        "MHW_M": MHW_M,
        "BERM_ELEV_NAVD_M": BERM_ELEV_NAVD_M,
        "beach start": BEACH_START_THR_M,
        "WATER_CLAMP_M": WATER_CLAMP_M,
        "SENTINEL_WATER_M": SENTINEL_WATER_M,
        "MIN_DUNE_H_M": MIN_DUNE_H_M,
        "TOPO_ROWS": TOPO_ROWS,
        "ALONG_COLS": ALONG_COLS,
        "OCEAN_LOC": OCEAN_LOC,
        "ALONGSHORE_FLIP": ALONGSHORE_FLIP,
        "DEFAULT_WINDOW_PX": DEFAULT_WINDOW_PX,
        "CLIP_WINDOW_TO_BEACH": CLIP_WINDOW_TO_BEACH,
        "MIN_DRAWN_BAND_WIDTH_PX": MIN_DRAWN_BAND_WIDTH_PX,
        "FREEHAND_SMOOTH_PX": FREEHAND_SMOOTH_PX,
        "AUTO_SWITCH_DRAW_BOUNDARY": AUTO_SWITCH_DRAW_BOUNDARY,
        "SAVE_CROSS_SECTION_FIGS": SAVE_CROSS_SECTION_FIGS,
        "SHOW_CROSS_SECTION_FIGS": SHOW_CROSS_SECTION_FIGS,
        "REQUIRE_CROSS_SECTION_APPROVAL": REQUIRE_CROSS_SECTION_APPROVAL,
        "CROSS_SECTION_CONTEXT_PX": CROSS_SECTION_CONTEXT_PX,
        "cross section timing": "before dune/interior separation",
        "cross section profile": "middle alongshore row (n_fill // 2)",
        "window schema": "profile-varying freehand band v3",
        "DUNE_ROWS": DUNE_ROWS,
        "dune output orientation": "alongshore x dune rows",
        "follow local dune line": FOLLOW_LOCAL_DUNE_LINE,
        "shift interior": True,
        "constant interior row": False,
        "FILL_MISSING_DUNE": FILL_MISSING_DUNE,
    }


def settings_row(res: dict, w: dict | None) -> dict:
    """One sheet row per domain: settings used + profile-varying band statistics."""
    geometry = ("profile-varying" if
                (res["i0_min"] != res["i0_max"] or res["i1_min"] != res["i1_max"])
                else "horizontal")
    return {
        "domain": domain_number(res["stem"]),
        "stem": res["stem"],
        "section": section_for(res["stem"]),
        "root GIS domains": LOAD_PATH.name,
        "root dunes/topo": RUN_NAME,
        "figure name": f"{fig_stem(res['stem'])}_gis_vs_processed.png",
        "beach start": BEACH_START_THR_M,
        # Existing column names now contain alongshore means for compatibility.
        "dune window start": round(res["i0"], 2),
        "dune window end": round(res["i1"], 2),
        "dune window": round(res["window_width_mean"], 2),
        "window geometry": geometry,
        "seaward boundary min": res["i0_min"],
        "seaward boundary max": res["i0_max"],
        "landward boundary min": res["i1_min"],
        "landward boundary max": res["i1_max"],
        "band width min": res["window_width_min"],
        "band width max": res["window_width_max"],
        "window source": "picked" if w else "default",
        "window schema": (w or {}).get("schema", "default-horizontal"),
        "clip window to beach": CLIP_WINDOW_TO_BEACH,
        "shift interior": True,
        "constant interior row": False,
        "follow local dune line": FOLLOW_LOCAL_DUNE_LINE,
        "dune rows": DUNE_ROWS,
        "MHW (m NAVD88)": MHW_M,
        "berm (m NAVD88)": BERM_ELEV_NAVD_M,
        "water clamp (m MHW)": WATER_CLAMP_M,
        "interior start row": round(res["start_island"], 2)
        if res["start_island"] is not None else "",
        "interior start min": res["interior_start_min"],
        "interior start max": res["interior_start_max"],
        "dunes found": res["n_found"],
        "profiles": res["n_along"],
        "dunes filled": res["n_filled"],
        "mean dune h (m)": round(res["mean_dune_h_m"], 3),
        "min dune h (m)": round(res["min_dune_h_m"], 3),
        "max dune h (m)": round(res["max_dune_h_m"], 3),
        "dune array shape": str(res["dune_shape"]),
        "dune saved shape": str(res.get("dune_saved_shape", "")),
        "interior rows": res["interior_rows"],
        "interior cols": res["interior_cols"],
        "mean interior elev (m MHW)": round(res["mean_interior_elev_m"], 3),
        "middle profile index": res.get("middle_profile_index", ""),
        "middle beach start": res.get("middle_beach_start", ""),
        "middle drawn window start": res.get("middle_window_start", ""),
        "middle drawn window end": res.get("middle_window_end", ""),
        "middle effective window start": res.get("middle_effective_start", ""),
        "middle effective window end": res.get("middle_effective_end", ""),
        "middle detected dune cell": res.get("middle_dune_loc_detected", ""),
        "middle used dune cell": res.get("middle_dune_loc_used", ""),
        "middle dune elev (m MHW)": (
            round(res["middle_dune_elev_m_mhw"], 3)
            if np.isfinite(res.get("middle_dune_elev_m_mhw", np.nan)) else ""
        ),
        "middle dune elev (m NAVD88)": (
            round(res["middle_dune_elev_m_navd88"], 3)
            if np.isfinite(res.get("middle_dune_elev_m_navd88", np.nan)) else ""
        ),
        "middle dune h above berm (m)": (
            round(res["middle_dune_height_above_berm_m"], 3)
            if np.isfinite(res.get("middle_dune_height_above_berm_m", np.nan)) else ""
        ),
        "middle interior start": res.get("middle_interior_start", ""),
        "middle cross section figure": res.get("middle_cross_section_file", ""),
        "road mask file": res.get("road_mask_file", ""),
        "profiles containing road": res.get("road_profiles", 0),
        "mean dune-to-road (m)": (
            round(res["mean_dune_to_road_m"], 2)
            if np.isfinite(res.get("mean_dune_to_road_m", np.nan)) else ""
        ),
        "min dune-to-road (m)": (
            round(res["min_dune_to_road_m"], 2)
            if np.isfinite(res.get("min_dune_to_road_m", np.nan)) else ""
        ),
        "max dune-to-road (m)": (
            round(res["max_dune_to_road_m"], 2)
            if np.isfinite(res.get("max_dune_to_road_m", np.nan)) else ""
        ),
        "middle road present": res.get("middle_road_present", ""),
        "middle road center cell": res.get("middle_road_center_cell", ""),
        "middle dune-to-road (m)": (
            round(res["middle_dune_to_road_m"], 2)
            if np.isfinite(res.get("middle_dune_to_road_m", np.nan)) else ""
        ),
        "topo file": res["topo_file"],
        "dune file": res["dune_file"],
        "picked": (w or {}).get("picked", ""),
    }


def write_settings_sheet(rows: list, base_path: Path) -> None:
    """Write the per-domain settings sheet as CSV, plus XLSX if pandas is around."""
    if not rows:
        return
    base_path.parent.mkdir(parents=True, exist_ok=True)

    csv_path = base_path.with_suffix(".csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[sheet] {csv_path}")

    try:
        import pandas as pd
        cfg = global_config()
        cfg_df = pd.DataFrame({"setting": list(cfg.keys()),
                               "value": [str(v) for v in cfg.values()]})
        xlsx_path = base_path.with_suffix(".xlsx")
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
            pd.DataFrame(rows).to_excel(xw, sheet_name="domains", index=False)
            cfg_df.to_excel(xw, sheet_name="global_config", index=False)
        print(f"[sheet] {xlsx_path}")
    except ImportError:
        print("[info] pandas/openpyxl not available; wrote CSV only")


def write_manifest(rows: list, path: Path) -> None:
    """Plain-text record so a run folder found later explains itself."""
    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = global_config()
    width = max(len(k) for k in cfg)
    lines = [
        "=" * 78,
        f"Pea Island dune & topography extraction  --  run {RUN_NAME}",
        "=" * 78,
        "",
        "CONTENTS",
        "  topography\\            interior elevation arrays (dam), CASCADE input",
        "  dunes\\                 dune height above berm (dam), CASCADE input",
        "  figures\\gis_vs_processed\\   raw DEM vs what CASCADE receives, per domain",
        "  figures\\qc\\            dune detection QC, per domain",
        "  figures\\middle_cross_sections\\   interactive profile + road overlay",
        f"  {SHEET_SAVE_PATH.name}.xlsx / .csv   per-domain settings + results",
        f"  {SUMMARY_FIG_PATH.name}   all domains on one page",
        f"  {PLAN_FIG_STEM}_<year>_<mode>.png   plan view, per year per mode",
        "",
        "PICKS (not in this folder -- shared across versions, not regenerable)",
        f"  {WINDOW_JSON}",
        "",
        "SETTINGS",
    ]
    lines += [f"  {k:<{width}} : {v}" for k, v in cfg.items()]
    lines += ["", f"DOMAINS WRITTEN: {len(rows)}"]
    if rows:
        nums = [r["domain"] for r in rows if r["domain"] is not None]
        if nums:
            lines.append(f"  range: {min(nums)} - {max(nums)}")
        defaults = [r["stem"] for r in rows if r["window source"] == "default"]
        if defaults:
            lines.append(f"  !! fallback default window (never picked): "
                         f"{', '.join(defaults)}")
        filled = [r["stem"] for r in rows if r["dunes filled"] > 0]
        if filled:
            lines.append(f"  !! profiles with no dune found, filled with "
                         f"{MIN_DUNE_H_M} m: {', '.join(filled)}")
    lines.append("")
    path.write_text("\n".join(lines))
    print(f"[manifest] {path}")


def summary_figure(rows: list, path: Path) -> None:
    """Every domain on one page: windows, dune heights, interior extent."""
    if not rows:
        return
    rows = sorted([r for r in rows if r["domain"] is not None],
                  key=lambda r: r["domain"])
    if not rows:
        return
    d = np.array([r["domain"] for r in rows])

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(13, 9), sharex=True,
                                        gridspec_kw={"hspace": 0.12})

    # optional Pea Island section bands + labels
    for k, ((lo, hi), label) in enumerate(SECTIONS):
        if not ((d >= lo) & (d <= hi)).any():
            continue
        for ax in (ax0, ax1, ax2):
            ax.axvspan(lo - 0.5, hi + 0.5, color="#90AFC5",
                       alpha=0.18 if k % 2 == 0 else 0.08, lw=0, zorder=0)
        trans = blended_transform_factory(ax0.transData, ax0.transAxes)
        ax0.text((lo + hi) / 2, 0.96, label, transform=trans, ha="center",
                 va="top", fontsize=8, rotation=0,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                           alpha=0.85))

    # 1) search window and interior start
    ax0.fill_between(d, [r["dune window start"] for r in rows],
                     [r["dune window end"] for r in rows],
                     color="#FF8C00", alpha=0.45, step="mid",
                     label="mean dune search band")
    ax0.step(d, [r["interior start row"] for r in rows], where="mid",
             color="#B71C1C", lw=1.4, label="interior start row")
    ax0.set_ylabel("cross-shore cell\n(0 = ocean)")
    ax0.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax0.set_title(f"Pea Island dune & topography extraction  —  {RUN_NAME}  —  "
                  f"{len(rows)} domains", fontsize=12)

    # 2) dune height
    ax1.fill_between(d, [r["min dune h (m)"] for r in rows],
                     [r["max dune h (m)"] for r in rows],
                     color="#FF8C00", alpha=0.25, label="min-max")
    ax1.plot(d, [r["mean dune h (m)"] for r in rows], color="#FF8C00", lw=1.6,
             marker="o", ms=3, label="mean")
    ax1.axhline(MIN_DUNE_H_M, color="0.4", ls="--", lw=1.0,
                label=f"floor {MIN_DUNE_H_M} m")
    ax1.set_ylabel("dune height\nabove berm (m)")
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # 3) interior extent + data-quality flags
    ax2.plot(d, [r["interior rows"] for r in rows], color="#1565C0", lw=1.6,
             marker="o", ms=3, label="interior rows")
    ax2.set_ylabel("interior rows")
    ax2.set_xlabel("Pea Island CASCADE domain ID")

    flag = np.array([r["dunes filled"] > 0 or r["window source"] == "default"
                     for r in rows])
    if flag.any():
        trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)
        ax2.plot(d[flag], np.full(flag.sum(), 0.04), transform=trans2, lw=0,
                 marker="v", ms=6, color="#B71C1C", clip_on=False,
                 label="default window / filled dunes")
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax2.set_xlim(d.min() - 0.5, d.max() + 0.5)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[summary] {path}")


def resolve_offset_file(year: int, configured) -> Path | None:
    """
    Return the offset CSV for a year. If the configured path is wrong, search
    OFFSET_DIR rather than silently skipping -- the hindcast_* subfolder naming
    is a convention this script only partly knows.
    """
    configured = Path(configured)
    if configured.exists():
        return configured
    root = Path(OFFSET_DIR)
    if not root.is_dir():
        return None
    hits = sorted(p for p in root.rglob("Island_Dune_Offsets*.csv")
                  if str(year) in p.name)
    if len(hits) == 1:
        print(f"[path] offsets {year}: configured path not found, found instead")
        print(f"       {hits[0]}")
        print(f"       -> update OFFSET_FILES[{year}] to match")
        return hits[0]
    if len(hits) > 1:
        print(f"[path] offsets {year}: configured path not found; "
              f"{len(hits)} candidates, none picked (ambiguous):")
        for h in hits:
            print(f"       {h}")
    return None


def check_paths() -> bool:
    """Preflight every path before doing any work. False = don't run."""
    print("=" * 78)
    print(f"PATH CHECK  —  run {RUN_NAME}")
    print("=" * 78)
    ok = True

    def row(tag, label, path, note=""):
        print(f"  [{tag:^7}] {label:<20} {path}{note}")

    for label, pth in (("project root", PROJECT_ROOT), ("PeaIsland_init", INIT_ROOT)):
        if Path(pth).is_dir():
            row("ok", label, pth)
        else:
            row("MISSING", label, pth)
            ok = False

    if Path(LOAD_PATH).is_dir():
        all_npy = [f for f in os.listdir(LOAD_PATH)
                   if Path(f).suffix.lower() == ".npy"]
        files = [f for f in all_npy if is_real_domain_file(f)]

        found_ids = {domain_number(f) for f in files}
        found_ids.discard(None)
        expected = set(REAL_DOMAIN_IDS)
        missing = sorted(expected - found_ids)

        # Detect two files resolving to the same domain ID.
        by_id = {}
        for filename in files:
            by_id.setdefault(domain_number(filename), []).append(filename)
        duplicates = {did: vals for did, vals in by_id.items() if len(vals) > 1}

        if missing:
            row("MISSING", "DEM domains", LOAD_PATH,
                f"   missing {len(missing)} expected file(s): {missing}")
            ok = False
        elif duplicates:
            row("DUPLICATE", "DEM domains", LOAD_PATH,
                f"   multiple files resolve to the same domain: {duplicates}")
            ok = False
        else:
            row("ok", "DEM domains", LOAD_PATH,
                f"   ({len(files)} files; recognized IDs "
                f"{FIRST_REAL_DOMAIN_ID}..{LAST_REAL_DOMAIN_ID})")

        ignored_npy = sorted(set(all_npy) - set(files))
        if ignored_npy:
            preview = ignored_npy[:8]
            suffix = " ..." if len(ignored_npy) > len(preview) else ""
            row("IGNORE", "other NPY files", LOAD_PATH,
                f"   {preview}{suffix}")
    else:
        row("MISSING", "DEM domains", LOAD_PATH,
            f"   create this folder and add "
            f"resampled_1996_domains_{FIRST_REAL_DOMAIN_ID}.npy through "
            f"resampled_1996_domains_{LAST_REAL_DOMAIN_ID}.npy")
        ok = False

    if SHOW_ROAD:
        road_root = Path(ROAD_MASK_FOLDER)
        if not road_root.is_dir():
            row("MISSING", "road masks", road_root,
                "   create this folder and add aligned binary .npy masks")
            if REQUIRE_ROAD_MASKS:
                ok = False
        else:
            missing_road = [
                domain_id for domain_id in REAL_DOMAIN_IDS
                if find_road_mask_path(domain_id) is None
            ]
            found_count = NUM_REAL_DOMAINS - len(missing_road)
            if missing_road:
                tag = "MISSING" if REQUIRE_ROAD_MASKS else "WARN"
                row(tag, "road masks", road_root,
                    f"   found {found_count}/{NUM_REAL_DOMAINS}; "
                    f"missing domains {missing_road}")
                if REQUIRE_ROAD_MASKS:
                    ok = False
            else:
                row("ok", "road masks", road_root,
                    f"   ({found_count} aligned binary masks)")

    for year in sorted(OFFSET_FILES):
        found = resolve_offset_file(year, OFFSET_FILES[year])
        if found is None:
            row("MISSING", f"offsets {year}", OFFSET_FILES[year],
                "   (figures will skip this year)")
        elif Path(found) != Path(OFFSET_FILES[year]):
            row("FOUND", f"offsets {year}", found, "   (not where configured)")
        else:
            row("ok", f"offsets {year}", found)

    wj = Path(WINDOW_JSON)
    if wj.exists():
        try:
            n = len([k for k in json.loads(wj.read_text()) if not k.startswith("_")])
            row("ok", "picks", wj, f"   ({n} domains picked)")
        except Exception as e:
            row("BAD", "picks", wj, f"   unreadable: {e}")
            ok = False
    else:
        row("new", "picks", wj, "   (created on first pick)")

    for legacy in (Path(RUN_DIR) / f"{PLAN_FIG_STEM}.png",):
        if legacy.exists():
            row("STALE", "old plan view", legacy,
                "\n            ^ pre-dates the per-year/per-mode naming, nothing "
                "overwrites it — delete it")

    rd = Path(RUN_DIR)
    row("exists" if rd.exists() else "new", "run folder", rd,
        "   (will be overwritten)" if rd.exists() else "   (created on run)")

    print("=" * 78)
    if not ok:
        print("[stop] required inputs missing — fix the paths above and re-run.\n")
    return ok


def load_offsets() -> dict:
    """Return {year: (Pea_domain_ids, offset_m)} for configured offset CSVs."""
    out = {}
    for year, configured in OFFSET_FILES.items():
        path = resolve_offset_file(year, configured)
        if path is None:
            print(f"[warn] offset file not found, skipping {year}: {configured}")
            continue

        values = np.loadtxt(path, skiprows=1, delimiter=",", ndmin=2).astype(float)
        if values.shape[1] > 1:
            print(f"[offset] {year}: {values.shape[1]} columns, "
                  f"using column {OFFSET_COLUMN}")
            values = values[:, OFFSET_COLUMN]
        else:
            values = values[:, 0]

        padded = N_LEFT_BUFFER_DOMAINS + NUM_REAL_DOMAINS + N_RIGHT_BUFFER_DOMAINS
        if values.size == padded:
            i0 = N_LEFT_BUFFER_DOMAINS
            i1 = i0 + NUM_REAL_DOMAINS
            print(f"[offset] {year}: padded file ({padded} rows), stripping "
                  f"left={N_LEFT_BUFFER_DOMAINS}, right={N_RIGHT_BUFFER_DOMAINS}")
            values = values[i0:i1]

        if values.size != NUM_REAL_DOMAINS:
            print(f"[warn] offsets {year}: expected {NUM_REAL_DOMAINS} real-domain "
                  f"values (or {padded} padded rows), found {values.size}; skipped")
            continue

        if OFFSET_ROW_ORDER == "ascending":
            dom = np.asarray(REAL_DOMAIN_IDS, dtype=int)
        elif OFFSET_ROW_ORDER == "descending":
            dom = np.asarray(REAL_DOMAIN_IDS[::-1], dtype=int)
        else:
            raise ValueError("OFFSET_ROW_ORDER must be 'ascending' or 'descending'")

        order = np.argsort(dom)
        out[year] = (dom[order], values[order])
        print(f"[offset] {year}: domains {FIRST_REAL_DOMAIN_ID}..{LAST_REAL_DOMAIN_ID}, "
              f"{np.nanmin(values):.0f}-{np.nanmax(values):.0f} m ({path.name})")
    return out


def _island_norm():
    """Poster colormap: terrain with 0 m pinned to colormap position 0.35."""
    lo, hi, pos = ISLAND_ELEV_MIN_M, ISLAND_ELEV_MAX_M, ISLAND_SEA_LEVEL_POS

    def fwd(x):
        out = np.where(x < 0.0, pos * (x - lo) / (0.0 - lo),
                       pos + (1.0 - pos) * x / hi)
        return np.where(np.isnan(x), np.nan, out)

    def inv(x):
        return np.where(x < pos, lo + (x / pos) * (0.0 - lo),
                        (x - pos) / (1.0 - pos) * hi)

    cmap = plt.cm.terrain.copy()
    cmap.set_bad(color=ISLAND_OCEAN_COLOR)
    return cmap, FuncNorm((fwd, inv), vmin=lo, vmax=hi)


def _build_island_canvas(recs, offset_m_by_domain, mode):
    """Stitch processed domains onto one plan-view canvas at their dune offsets."""
    use = [(n, r) for n, r in recs if n in offset_m_by_domain]
    if not use:
        return None, None, None
    off_cells = [int(round(offset_m_by_domain[n] / CELL_SIZE_M)) for n, _ in use]

    grids, dunes, cropped = [], [], []
    for n_dom, r in use:
        g = r["topo_dm"] * CELL_SIZE_M                 # dam -> m MHW
        if mode == "padded":
            # pad landward, matching where dune_topo_extractor_from_GIS.py left its
            # sentinel: interior row 0 is the ocean side, rows increase landward
            n = g.shape[0]
            if n < ISLAND_PAD_ROWS:
                g = np.vstack([g, np.full((ISLAND_PAD_ROWS - n, g.shape[1]),
                                          SENTINEL_WATER_M)])
            elif n > ISLAND_PAD_ROWS:
                n_land = int(np.sum(g[ISLAND_PAD_ROWS:] > SENTINEL_WATER_M + 1e-9))
                if n_land:
                    cropped.append((n_dom, n, n_land))
                g = g[:ISLAND_PAD_ROWS]
        if ISLAND_SENTINEL_AS_OCEAN:
            g = np.where(g <= SENTINEL_WATER_M + 1e-9, np.nan, g)
        # Saved dunes are (alongshore, dune rows); transpose to
        # (dune rows, alongshore) for the stitched plan-view canvas.
        d = (r["dune_dm"] * CELL_SIZE_M
             + (BERM_ELEV_NAVD_M - MHW_M)).T   # -> m MHW
        if ISLAND_FLIP_ALONGSHORE:
            g = np.fliplr(g)
            d = np.fliplr(d)
        grids.append(g)
        dunes.append(d)

    if cropped:
        print(f"[planview] ISLAND_PAD_ROWS={ISLAND_PAD_ROWS} "
              f"({ISLAND_PAD_ROWS * CELL_SIZE_M:.0f} m) crops real land from "
              f"{len(cropped)} domain(s):")
        for n_dom, n_rows, n_land in cropped[:8]:
            print(f"           D{n_dom}: {n_rows} rows stored, {n_land} land cells lost")
        if len(cropped) > 8:
            print(f"           ... and {len(cropped) - 8} more")

    max_rows = max(g.shape[0] for g in grids)
    canvas_rows = max(off_cells) + max_rows + 5
    total_cols = sum(g.shape[1] for g in grids)
    canvas = np.full((canvas_rows, total_cols), np.nan)

    col = 0
    starts = []
    for k, (g, d) in enumerate(zip(grids, dunes)):
        n_rows, n_cols = g.shape
        starts.append(col)
        origin = off_cells[k]
        end = min(origin + n_rows, canvas_rows)
        canvas[origin:end, col:col + n_cols] = g[:end - origin, :]
        if ISLAND_INCLUDE_DUNE:
            n_dune_rows, n_dune_cols = d.shape
            dune_top = max(0, origin - n_dune_rows)
            src_top = n_dune_rows - (origin - dune_top)
            width = min(n_cols, n_dune_cols)
            if origin > dune_top and width > 0:
                canvas[dune_top:origin, col:col + width] = d[src_top:, :width]
        col += n_cols

    return canvas, np.array(starts), [n for n, _ in use]


def island_plan_figure(summary: list, offsets: dict, run_dir: Path) -> None:
    """
    Plan view of the processed dune + interior for Pea Island domains 80-119 at
    offsets. Writes ONE FIGURE PER OFFSET YEAR PER CROSS-SHORE MODE, styled to
    match plot_initialization_poster_no_border.py.
    """
    recs = sorted([(domain_number(r["stem"]), r) for r in summary
                   if domain_number(r["stem"]) is not None])
    if not recs or not offsets:
        return
    cmap, norm = _island_norm()

    for year in sorted(offsets):
        dom, v = offsets[year]
        omap = {int(a): float(b) for a, b in zip(dom, v)}
        for mode in ISLAND_CROSS_SHORE_MODES:
            canvas, starts, used = _build_island_canvas(recs, omap, mode)
            if canvas is None:
                print(f"[planview] {year}/{mode}: no domains overlap the offsets, skipped")
                continue

            n_cs, n_al = canvas.shape
            fig_w = 20.0
            fig_h = min(max(4.5, fig_w * (n_cs / n_al) * 1.8), 7.5)   # poster aspect
            fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
            ax = fig.add_axes([0.06, 0.18, 0.88, 0.68])
            ax.set_facecolor(ISLAND_OCEAN_COLOR)

            im = ax.pcolormesh(np.ma.masked_invalid(canvas), cmap=cmap, norm=norm,
                               shading="auto", rasterized=True)
            ax.set_xlim(0, n_al)
            ax.set_ylim(0, n_cs)

            cax = fig.add_axes([0.955, 0.18, 0.013, 0.68])
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label("Elevation (m MHW)", fontsize=12, color="#1a1a2e",
                           labelpad=10, rotation=270)
            cbar.ax.yaxis.set_tick_params(color="#1a1a2e", labelcolor="#1a1a2e")
            cbar.outline.set_edgecolor("#cccccc")
            cbar.set_ticks([-1, 0, 1, 2, 3, 4])

            ticks, labels = [], []
            for k, n in enumerate(used):
                if n % 5 == 0 or n == FIRST_REAL_DOMAIN_ID:
                    ticks.append(starts[k] + ALONG_COLS // 2)
                    labels.append(str(n))
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels, fontsize=9)
            ax.set_xlabel("Pea Island CASCADE domain ID", fontsize=12,
                          labelpad=8)
            ax.set_ylabel("Cross-shore cell (offset frame)", fontsize=12)
            for k, n in enumerate(used):
                if n % 10 == 0:
                    ax.axvline(starts[k] - 0.5, color="#aaaaaa", lw=0.4, alpha=0.5,
                               zorder=2)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("bottom", "left"):
                ax.spines[sp].set_color("#999999")

            what = "dune + interior" if ISLAND_INCLUDE_DUNE else "interior"
            if mode == "padded":
                extent = f"{ISLAND_PAD_ROWS} cells / {ISLAND_PAD_ROWS * CELL_SIZE_M:.0f} m"
                note = f"{what}, every domain padded to {extent} cross-shore"
            else:
                note = f"{what}, each domain trimmed to its own island"
            ax.set_title(f"Pea Island — CASCADE Initialization  |  {year} offsets  "
                         f"|  {DEM_YEAR} extracted {note}  ({len(used)} domains)",
                         fontsize=14, fontweight="bold", color="#1a1a2e", pad=12)

            path = Path(run_dir) / f"{PLAN_FIG_STEM}_{year}_{mode}.png"
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(f"[planview] {path}")


def island_figure(summary: list, offsets: dict, path: Path) -> None:
    """
    All domains together: measured dune offsets vs the crest this run extracted.

    The measured offsets and the extracted crest live in DIFFERENT frames -- the
    offsets are in the model's common cross-shore frame, the extracted crest is in
    the per-domain raw DEM array frame (m landward of that array's cell 0). They
    are NOT differenced here. They're plotted on separate panels, and the script
    reports the correlation between them so the frame relationship is testable
    rather than assumed.
    """
    recs = sorted([(domain_number(r["stem"]), r) for r in summary
                   if domain_number(r["stem"]) is not None])
    if not recs:
        return
    d = np.array([n for n, _ in recs])

    pos_mean, pos_lo, pos_hi = [], [], []
    for _, r in recs:
        loc = r["dune_loc"].astype(float)
        loc[loc < 0] = np.nan
        pm = (loc + r["c0"]) * CELL_SIZE_M
        pos_mean.append(np.nanmean(pm))
        pos_lo.append(np.nanmin(pm))
        pos_hi.append(np.nanmax(pm))
    pos_mean = np.array(pos_mean)

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(13, 10), sharex=True,
                                        gridspec_kw={"hspace": 0.12})

    for k, ((lo, hi), label) in enumerate(SECTIONS):
        for ax in (ax0, ax1, ax2):
            ax.axvspan(lo - 0.5, hi + 0.5, color="#90AFC5",
                       alpha=0.18 if k % 2 == 0 else 0.08, lw=0, zorder=0)
        trans = blended_transform_factory(ax0.transData, ax0.transAxes)
        ax0.text((lo + hi) / 2, 0.97, label, transform=trans, ha="center",
                 va="top", fontsize=8,
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                           alpha=0.85))

    # 1) measured offsets, all configured years
    for year in sorted(offsets):
        dom, v = offsets[year]
        ax0.plot(dom, v, lw=1.6, label=f"measured dune offset {year}")
    ax0.set_ylabel("measured offset\n(m, common frame)")
    if offsets:
        ax0.legend(loc="upper right", fontsize=8, framealpha=0.9)
    else:
        ax0.text(0.5, 0.5, "No offset CSV configured", ha="center", va="center",
                 transform=ax0.transAxes, fontsize=10, color="0.4")
    sign = "seaward +" if OFFSET_SEAWARD_POSITIVE else "landward +"
    ax0.set_title(f"Island dune offsets vs extracted crest  —  {RUN_NAME}  "
                  f"({len(recs)} domains, row0={OFFSET_ROW_ORDER}, {sign})",
                  fontsize=12)

    # 2) crest extracted from the DEM this run, raw-array frame
    ax1.fill_between(d, pos_lo, pos_hi, color="#FF8C00", alpha=0.25,
                     label="within-domain min-max")
    ax1.plot(d, pos_mean, color="#FF8C00", lw=1.6, marker="o", ms=3,
             label=f"extracted crest, {DEM_YEAR} DEM")
    ax1.set_ylabel("extracted crest\n(m from raw cell 0)")
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # correlation against each measured year -- diagnostic for the frame relation
    notes = []
    for year in sorted(offsets):
        dom, v = offsets[year]
        common = np.intersect1d(d, dom)
        if common.size < 3:
            continue
        a = pos_mean[np.isin(d, common)]
        b = v[np.isin(dom, common)]
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 3:
            continue
        r = float(np.corrcoef(a[ok], b[ok])[0, 1])
        notes.append(f"r(extracted, {year}) = {r:+.2f}  (n={ok.sum()})")
    if notes:
        trans = blended_transform_factory(ax1.transAxes, ax1.transAxes)
        ax1.text(0.01, 0.06, "   |   ".join(notes), transform=trans, fontsize=8,
                 va="bottom", ha="left",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                           alpha=0.9))

    # 3) measured change between earliest and latest configured offset years
    years = sorted(offsets)
    if len(years) >= 2:
        year_a, year_b = years[0], years[-1]
        dom_a, va = offsets[year_a]
        dom_b, vb = offsets[year_b]
        common = np.intersect1d(dom_a, dom_b)
        dt = float(year_b - year_a)
        ch = (vb[np.isin(dom_b, common)] - va[np.isin(dom_a, common)]) / dt
        ax2.axhline(0, color="0.5", lw=1.0)
        ax2.plot(common, ch, lw=1.4, marker="o", ms=3)
        ax2.set_ylabel(f"measured change\n{year_a}→{year_b} (m/yr)")
    else:
        ax2.text(0.5, 0.5, "Configure at least two offset years for change rates",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=10, color="0.4")
    ax2.set_xlabel("Pea Island CASCADE domain ID")
    ax2.set_xlim(min(d.min(), FIRST_REAL_DOMAIN_ID) - 0.5,
                  max(d.max(), LAST_REAL_DOMAIN_ID) + 0.5)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[island] {path}")
    for n in notes:
        print(f"         {n}")


def load_windows(path: Path) -> dict:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {"_meta": {}}


def saved_window_lines(w: dict, n_along: int, n_cross: int,
                       default_pair) -> tuple[np.ndarray, np.ndarray]:
    """Load new freehand bands or transparently upgrade old scalar i0/i1 picks."""
    if w is None:
        return normalize_window_lines(
            default_pair[0], default_pair[1], n_along, n_cross,
            MIN_DRAWN_BAND_WIDTH_PX,
        )

    if "i0_by_along" in w and "i1_by_along" in w:
        return normalize_window_lines(
            w["i0_by_along"], w["i1_by_along"], n_along, n_cross,
            MIN_DRAWN_BAND_WIDTH_PX,
        )

    # Backward compatibility with the original straight horizontal windows.
    if "i0" in w and "i1" in w:
        return normalize_window_lines(
            int(w["i0"]), int(w["i1"]), n_along, n_cross,
            MIN_DRAWN_BAND_WIDTH_PX,
        )

    print("[warn] saved window has no recognized boundaries; using default")
    return normalize_window_lines(
        default_pair[0], default_pair[1], n_along, n_cross,
        MIN_DRAWN_BAND_WIDTH_PX,
    )


def save_windows(path: Path, windows: dict) -> None:
    windows["_meta"] = {
        "updated": datetime.now().isoformat(timespec="seconds"),
        "schema": "profile-varying-freehand-band-v2",
        "load_path": str(LOAD_PATH),
        "mhw_m": MHW_M,
        "berm_elev_navd_m": BERM_ELEV_NAVD_M,
        "beach_start_thr_m": BEACH_START_THR_M,
        "water_clamp_m": WATER_CLAMP_M,
        "ocean_loc": OCEAN_LOC,
        "alongshore_flip": ALONGSHORE_FLIP,
        "min_drawn_band_width_px": MIN_DRAWN_BAND_WIDTH_PX,
        "freehand_smooth_px": FREEHAND_SMOOTH_PX,
        "note": "i0_by_along/i1_by_along are per-profile [start, stop) "
                "cross-shore indices in the ocean-first, water-trimmed array. "
                "Old scalar i0/i1 entries remain readable.",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(windows, f, indent=2, sort_keys=True)


# ==============================================================================
# MAIN
# ==============================================================================

def verify_interactive_backend() -> None:
    """Fail clearly when the picker cannot create a desktop GUI window."""
    backend = str(matplotlib.get_backend())
    print(f"[picker] Matplotlib backend: {backend}")

    noninteractive = {
        "agg", "pdf", "ps", "svg", "template", "cairo",
        "module://matplotlib_inline.backend_inline",
    }
    if REQUIRE_INTERACTIVE_PICKER and backend.lower() in noninteractive:
        raise RuntimeError(
            f"Matplotlib backend {backend!r} is noninteractive. The GIS picker "
            "cannot open. Run the .py file directly from PyCharm or Terminal, "
            "not inside a Jupyter notebook or Python Console."
        )


def main():
    needs_interactive = (
        MODE in ("pick", "pick_and_run")
        or (MODE in ("run", "pick_and_run") and REQUIRE_CROSS_SECTION_APPROVAL)
    )
    if needs_interactive:
        verify_interactive_backend()

    if MODE in ("pick", "pick_and_run"):
        print("[picker] Interactive GIS drawing is ENABLED")
        print(f"[picker] REPICK_EXISTING={REPICK_EXISTING}")
        print(f"[picker] Saved-pick file: {WINDOW_JSON}")
        if not REPICK_EXISTING:
            print("[picker] NOTE: domains already in the JSON will be skipped")

    if MODE in ("run", "pick_and_run") and REQUIRE_CROSS_SECTION_APPROVAL:
        print("[partition] Each domain opens plan view and cross-section together")
        print("[partition] drag plan boundary or click cross-section; Enter=accept")

    if not check_paths():
        return
    load_dir = Path(LOAD_PATH)
    topo_dir = Path(TOPO_SAVE_PATH)
    dune_dir = Path(DUNE_SAVE_PATH)
    topo_dir.mkdir(parents=True, exist_ok=True)
    dune_dir.mkdir(parents=True, exist_ok=True)

    names = sorted(
        [n for n in os.listdir(load_dir) if is_real_domain_file(n)],
        key=natural_key,
    )
    if PICK_DOMAINS is not None:
        wanted = set(PICK_DOMAINS)
        names = [n for n in names if domain_number(n) in wanted]
    print(f"[info] Found {len(names)} Pea Island domain file(s) in {load_dir}")

    windows = load_windows(Path(WINDOW_JSON))

    # ---------------- PICK-ONLY PASS ----------------
    # In pick_and_run mode, picking is handled inside the run loop so each
    # domain follows: picker -> cross-section QC -> extraction.
    if MODE == "pick":
        for name in names:
            stem = canonical_domain_stem(name)
            if not REPICK_EXISTING and stem in windows:
                w = windows[stem]
                if "i0_by_along" in w:
                    desc = describe_window_lines(w["i0_by_along"], w["i1_by_along"])
                else:
                    desc = f"horizontal [{w.get('i0')}, {w.get('i1')}]"
                print(f"[skip pick] {stem}: {desc} already saved; "
                      "set REPICK_EXISTING=True to reopen it")
                continue
            try:
                dom = load_profiles(load_dir / name)
            except Exception as e:
                print(f"[skip] {name}: {e}")
                continue

            prof_arr = dom["z"]
            n_along, n_cross = prof_arr.shape
            default_pair = default_window(prof_arr, dom["start_beach"])

            # REPICK_EXISTING=True opens the old drawing as the editable starting
            # geometry. Otherwise the picker begins from the default flat band.
            existing = windows.get(stem) if REPICK_EXISTING else None
            init_lines = saved_window_lines(
                existing, n_along, n_cross, default_pair
            )
            action, i0_line, i1_line = pick_window(
                stem, prof_arr, dom["start_beach"], init_lines,
                road_mask=dom.get("road_mask"),
            )

            if action == "quit":
                print("[info] quit picking; saved bands kept")
                break
            if action == "skip":
                print(f"[pick] {stem}: skipped -> default horizontal band {default_pair}")
                continue

            windows[stem] = {
                "schema": "profile-varying-freehand-band-v2",
                # Scalar means are retained so old inspection scripts do not fail.
                "i0": int(round(float(np.mean(i0_line)))),
                "i1": int(round(float(np.mean(i1_line)))),
                "i0_by_along": [int(v) for v in i0_line],
                "i1_by_along": [int(v) for v in i1_line],
                "n_cross_trimmed": int(n_cross),
                "n_along": int(n_along),
                "trim_offset_c0": int(dom["c0"]),
                "picked": datetime.now().isoformat(timespec="seconds"),
            }
            save_windows(Path(WINDOW_JSON), windows)
            print(f"[pick] {stem}: {describe_window_lines(i0_line, i1_line)} saved")

    # ---------------- RUN PASS ----------------
    if MODE in ("run", "pick_and_run"):
        summary, rows = [], []
        for name in names:
            stem = canonical_domain_stem(name)
            try:
                dom = load_profiles(load_dir / name)
            except Exception as e:
                print(f"[skip] {name}: {e}")
                continue
            prof_arr, start_beach = dom["z"], dom["start_beach"]
            n_along, n_cross = prof_arr.shape
            default_pair = default_window(prof_arr, start_beach)

            w = windows.get(stem)

            # --------------------------------------------------------------
            # DOMAIN-BY-DOMAIN PICKING FOR MODE = "pick_and_run"
            # --------------------------------------------------------------
            # This is intentionally inside the run loop. The order is now:
            #
            #     1. draw/review the dune-search band
            #     2. inspect the original middle cross-section
            #     3. accept the proposed dune/interior split
            #     4. create and save the CASCADE arrays
            #
            # Therefore the cross-section appears immediately after the picker
            # for the same domain, rather than after all 40 domains are picked.
            if MODE == "pick_and_run":
                should_open_picker = REPICK_EXISTING or w is None

                if should_open_picker:
                    init_lines = saved_window_lines(
                        w, n_along, n_cross, default_pair
                    )
                    picker_action, picked_i0, picked_i1 = pick_window(
                        stem, prof_arr, start_beach, init_lines,
                        road_mask=dom.get("road_mask"),
                    )

                    if picker_action == "quit":
                        print("[picker] run stopped; completed outputs are preserved")
                        break

                    if picker_action == "skip":
                        print(f"[picker] {stem}: skipped; no arrays written")
                        continue

                    i0_line, i1_line = picked_i0, picked_i1
                    windows[stem] = {
                        "schema": "profile-varying-freehand-band-v2",
                        "i0": int(round(float(np.mean(i0_line)))),
                        "i1": int(round(float(np.mean(i1_line)))),
                        "i0_by_along": [int(v) for v in i0_line],
                        "i1_by_along": [int(v) for v in i1_line],
                        "n_cross_trimmed": int(n_cross),
                        "n_along": int(n_along),
                        "trim_offset_c0": int(dom["c0"]),
                        "picked": datetime.now().isoformat(timespec="seconds"),
                    }
                    save_windows(Path(WINDOW_JSON), windows)
                    w = windows[stem]
                    print(
                        f"[pick] {stem}: "
                        f"{describe_window_lines(i0_line, i1_line)} saved"
                    )
                else:
                    i0_line, i1_line = saved_window_lines(
                        w, n_along, n_cross, default_pair
                    )
                    print(
                        f"[pick] {stem}: using saved band; "
                        "set REPICK_EXISTING=True to redraw it"
                    )

            else:
                # MODE = "run": use a saved band when available, otherwise the
                # fallback horizontal window.
                if w is None:
                    print(f"[warn] {stem}: no picked band, using default {default_pair}")
                else:
                    if w.get("n_cross_trimmed") != n_cross:
                        print(f"[warn] {stem}: trimmed width changed "
                              f"({w.get('n_cross_trimmed')} -> {n_cross}); "
                              "re-picking is best.")
                    if w.get("n_along") not in (None, n_along):
                        print(f"[warn] {stem}: alongshore count changed "
                              f"({w.get('n_along')} -> {n_along}); "
                              "saved lines interpolated.")

                i0_line, i1_line = saved_window_lines(
                    w, n_along, n_cross, default_pair
                )

            # --------------------------------------------------------------
            # PRE-EXTRACTION QC LOOP
            # --------------------------------------------------------------
            # The original middle cross-section is shown before
            # build_local_dunes_and_interior() is called. Press R to redraw
            # the freehand band and return to this QC view.
            skip_domain = False
            stop_run = False
            qc_info = {}
            interior_start_choice = (
                (w or {}).get("interior_start_by_along")
                if w is not None else None
            )

            while True:
                qc_action, qc_info, interior_start_choice = \
                    pre_extraction_cross_section_qc(
                        stem=stem,
                        prof_arr=prof_arr,
                        start_beach=start_beach,
                        i0_line=i0_line,
                        i1_line=i1_line,
                        fig_dir=Path(FIG_DIR_XSEC),
                        initial_interior_line=interior_start_choice,
                        road_mask=dom.get("road_mask"),
                    )

                if qc_action == "accept":
                    # Persist both the dune-search band and the manually chosen
                    # interior-start boundary so MODE="run" can reproduce it.
                    if stem not in windows:
                        windows[stem] = {
                            "schema": "profile-varying-freehand-band-v2",
                            "i0": int(round(float(np.mean(i0_line)))),
                            "i1": int(round(float(np.mean(i1_line)))),
                            "i0_by_along": [int(v) for v in i0_line],
                            "i1_by_along": [int(v) for v in i1_line],
                            "n_cross_trimmed": int(n_cross),
                            "n_along": int(n_along),
                            "trim_offset_c0": int(dom["c0"]),
                            "picked": datetime.now().isoformat(timespec="seconds"),
                        }
                    windows[stem]["interior_start_by_along"] = [
                        int(v) for v in interior_start_choice
                    ]
                    windows[stem]["interior_picked"] = datetime.now().isoformat(
                        timespec="seconds"
                    )
                    save_windows(Path(WINDOW_JSON), windows)
                    w = windows[stem]
                    break

                if qc_action == "skip":
                    print(f"[pre-QC] {stem}: skipped; no arrays written")
                    skip_domain = True
                    break

                if qc_action == "quit":
                    print("[pre-QC] run stopped; completed outputs are preserved")
                    stop_run = True
                    break

                if qc_action == "repick":
                    picker_action, new_i0, new_i1 = pick_window(
                        stem,
                        prof_arr,
                        start_beach,
                        (i0_line, i1_line),
                        road_mask=dom.get("road_mask"),
                    )
                    if picker_action == "quit":
                        stop_run = True
                        break
                    if picker_action == "skip":
                        print(f"[pre-QC] {stem}: picker skipped; domain not written")
                        skip_domain = True
                        break

                    i0_line, i1_line = new_i0, new_i1
                    windows[stem] = {
                        "schema": "profile-varying-freehand-band-v2",
                        "i0": int(round(float(np.mean(i0_line)))),
                        "i1": int(round(float(np.mean(i1_line)))),
                        "i0_by_along": [int(v) for v in i0_line],
                        "i1_by_along": [int(v) for v in i1_line],
                        "n_cross_trimmed": int(n_cross),
                        "n_along": int(n_along),
                        "trim_offset_c0": int(dom["c0"]),
                        "picked": datetime.now().isoformat(timespec="seconds"),
                    }
                    save_windows(Path(WINDOW_JSON), windows)
                    w = windows[stem]
                    print(
                        f"[pre-QC] {stem}: revised band saved; reopening QC"
                    )
                    continue

            if stop_run:
                break
            if skip_domain:
                continue

            # Only now are dune and interior arrays separated and written.
            res = extract_domain(
                stem,
                prof_arr,
                start_beach,
                i0_line,
                i1_line,
                topo_dir,
                dune_dir,
                interior_start_override=interior_start_choice,
            )
            if res is None:
                continue
            res.update(qc_info)
            res["c0"] = dom["c0"]
            attach_road_statistics(res, dom)
            if SAVE_QC_FIGS:
                qc_figure(
                    stem, prof_arr, start_beach, res, Path(FIG_DIR_QC),
                    road_mask=dom.get("road_mask"),
                )
            if SAVE_COMPARISON_FIGS:
                comparison_figure(stem, dom, res, Path(FIG_DIR_CMP))
            summary.append(res)
            rows.append(settings_row(res, w))

        if SAVE_SETTINGS_SHEET:
            write_settings_sheet(rows, Path(SHEET_SAVE_PATH))
        summary_figure(rows, Path(SUMMARY_FIG_PATH))
        if SAVE_ISLAND_FIG:
            _off = load_offsets()
            island_figure(summary, _off, Path(ISLAND_FIG_PATH))
            if SAVE_ISLAND_PLAN_FIG:
                island_plan_figure(summary, _off, Path(RUN_DIR))
        write_manifest(rows, Path(MANIFEST_PATH))

        print("\n" + "=" * 94)
        print(f"{'domain':<12}{'section':<24}{'mean band':>18}{'geometry':>14}"
              f"{'dunes':>8}{'mean h (m)':>12}{'interior':>10}")
        print("-" * 94)
        for r in rows:
            win = f"[{r['dune window start']:.1f},{r['dune window end']:.1f}]"
            print(f"{r['stem']:<12}{r['section'][:23]:<24}{win:>18}"
                  f"{r['window geometry']:>14}{r['dunes found']:>8}"
                  f"{r['mean dune h (m)']:>12.2f}{r['interior rows']:>10}")
        print("=" * 94)

        n_default = sum(1 for r in rows if r["window source"] == "default")
        n_filled = sum(1 for r in rows if r["dunes filled"] > 0)
        if n_default:
            print(f"[!] {n_default} domain(s) used the fallback default band "
                  f"— re-run MODE='pick' for those")
        if n_filled:
            print(f"[!] {n_filled} domain(s) had profiles with no dune found")
        if SHOW_ROAD:
            road_domains = sum(1 for r in rows if r.get("profiles containing road", 0) > 0)
            print(f"[road] road cells present in {road_domains}/{len(rows)} written domain(s)")
        print(f"\n[done] {len(summary)} domain(s). Everything for this run is in:"
              f"\n       {RUN_DIR}\n"
              f"       start with RUN_MANIFEST.txt and {SUMMARY_FIG_PATH.name}\n"
              f"       CASCADE inputs: {topo_dir}\n"
              f"                       {dune_dir}\n"
              f"       cross sections: {FIG_DIR_XSEC}")


if __name__ == "__main__":
    main()
