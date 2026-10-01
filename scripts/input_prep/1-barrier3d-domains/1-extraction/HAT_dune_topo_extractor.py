"""
Build each domain's Barrier3D dune and interior arrays from its DEM, with a hand-picked dune search window.

    python scripts/input_prep/1-barrier3d-domains/1-extraction/HAT_dune_topo_extractor.py

MODE (in CONFIG) picks the pass: "pick" drags a dune search window per domain
on the profile stack (NC-12 drawn for reference, saved after every domain),
"run" extracts with the saved windows, "pick_and_run" does both. Writes the
topography and dune arrays (dam), the picks JSON, a settings sheet, a manifest
and QC figures to the product's dune-topo/<version>/ folder.
Details: scripts/input_prep/1-barrier3d-domains/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""


from __future__ import annotations

import csv
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

import matplotlib
try:
    # Needed for a real, blocking, interactive picker window
    import tkinter  # noqa: F401
    matplotlib.use("TkAgg")
except Exception:
    print("[warn] TkAgg unavailable; the picker needs an interactive backend.")
import matplotlib.pyplot as plt
from matplotlib.colors import FuncNorm, ListedColormap
from matplotlib.transforms import blended_transform_factory
from matplotlib.widgets import SpanSelector

# Repo root, found by searching upward
_PATH_REPO = next(_p for _p in Path(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())

# Type sized for a projected slide
plt.rcParams.update({
    "font.size": 13,
    "axes.titlesize": 13.0,
    "axes.labelsize": 13.5,
    "xtick.labelsize": 12.5,
       "ytick.labelsize": 12.5,
    "legend.fontsize": 11.5,
    "axes.linewidth": 0.9,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

# --- CONFIG ------------------------------------------------------------------
# Mode

# v4 re-picks every window with the road on the picker, adjusting the saved v3 windows
MODE = "pick_and_run"      # "pick" | "run" | "pick_and_run"
# Domains to process, in both passes: the 90 modelled domains by default (options in README)
PICK_DOMAINS = list(range(1, 91))
# True for the v4 re-pick (seeded from v3); set False afterwards so pick_and_run resumes
REPICK_EXISTING = True     # False = skip domains already present in the JSON
SAVE_QC_FIGS = True        # per-domain dune-detection QC figure
SAVE_COMPARISON_FIGS = True  # per-domain raw GIS vs processed CASCADE input figure
SAVE_SETTINGS_SHEET = True   # per-domain settings/results sheet (csv + xlsx)

# Paths

# Run layout: one folder per settings version (tree in README); which product this run builds
TOPO_PRODUCT = "1984-start"

# scripts/ on the path, so array names come from the one resolver
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))   # scripts/ (1-extraction/ since 2026-09-09)
from site_layer.hat_topo_version import (array_name, year_for_product,  # noqa: E402
                              ROAD_LINE_VINTAGES)

VERSION = "v2"             # 2026-09-02 re-pick (was "v3" until the 2026-09-04
                           # renumber); what this writes: dune-topo/CURRENT decides what is read

# A label only, on two figures: no year in any filename since 2026-08-26
DEM_LABEL = TOPO_PRODUCT
# v5 reads the gap-filled arrays
DEM_NAME = "2009_pea_hatteras_filled"
# The product folder already says which run this is, so the run folder is the VERSION alone
RUN_NAME = VERSION

PROJECT_ROOT = Path(str(_PATH_REPO))
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
# Paths from hat_topo_version: the period-first layout since 2026-08-25
from site_layer.hat_topo_version import npy_dirs, product_dir  # noqa: E402
PRODUCT_DIR = product_dir(TOPO_PRODUCT)
LOAD_PATH = npy_dirs(TOPO_PRODUCT)[0]                            # the extraction half (2026-09-09)

DUNE_TOPO_ROOT = PRODUCT_DIR / "dune-topo"
RUN_DIR = DUNE_TOPO_ROOT / RUN_NAME
TOPO_SAVE_PATH = RUN_DIR / "topography"
DUNE_SAVE_PATH = RUN_DIR / "dunes"
FIG_DIR_CMP = RUN_DIR / "figures" / "gis_vs_processed"
FIG_DIR_QC = RUN_DIR / "figures" / "qc"
SHEET_SAVE_PATH = RUN_DIR / f"HAT_dune_topo_settings_{RUN_NAME}"   # .csv + .xlsx
SUMMARY_FIG_PATH = RUN_DIR / f"HAT_dune_topo_summary_{RUN_NAME}.png"
ISLAND_FIG_PATH = RUN_DIR / f"HAT_dune_topo_island_offsets_{RUN_NAME}.png"
PLAN_FIG_STEM = f"HAT_dune_topo_island_planview_{RUN_NAME}"  # + _{year}.png
MANIFEST_PATH = RUN_DIR / "RUN_MANIFEST.txt"

# Picks are frame-dependent and written back every domain, so v4 keeps its own file (seeded from v3)
PICK_SET = RUN_NAME
PICKS_DIR = PRODUCT_DIR / "1-extraction" / "picks"
WINDOW_JSON = PICKS_DIR / f"HAT_dune_search_windows_{PICK_SET}.json"

# No tag: names come from hat_topo_version.array_name()

# Island offsets

# Measured per-domain dune offsets, used to place domains in a common cross-shore frame
SAVE_ISLAND_FIG = True
# Each start's CURRENT build (2026-09-18)
from site_layer.hat_topo_version import BRIE_ROOT as OFFSET_DIR, offset_file  # noqa: E402
OFFSET_FILES = {
    1984: offset_file(1984, "input"),
    2004: offset_file(2004, "input"),
}
CELL_SIZE_M = 10.0             # DEM cell size, cross-shore and alongshore (DAM_TO_M)
NUM_REAL_DOMAINS = 90
N_BUFFER_DOMAINS = 15          # raw_offset files may be padded to 15+90+15 = 120 rows
OFFSET_COLUMN = 0              # for multi-year raw_offset files, which column to read

# No relevant year: the plan view falls back to every year it loaded
PRODUCT_YEAR = year_for_product(TOPO_PRODUCT, strict=False)

# Plan-view canvas, reproducing initialization_figures.py's absolute placement
SAVE_ISLAND_PLAN_FIG = True
# ISLAND_FLIP_ALONGSHORE was removed (2026-08-17)
ISLAND_INCLUDE_DUNE = True     # write the dune crest into canvas row raw_offset-1
ISLAND_ELEV_MIN_M = -1.0       # poster value. Raising this flattens land contrast:
ISLAND_ELEV_MAX_M = 4.0        #   land maps to ramp 0.35-1.0, so max=4 puts a 2 m
                               #   cell at tan, max=5 at pale yellow-green.
ISLAND_SEA_LEVEL_POS = 0.35    # colormap position for 0 m, as in the poster
ISLAND_SENTINEL_AS_OCEAN = False  # False = poster behaviour: sentinel (-3 m) water
                                  # cells clip to terrain's navy; True = sentinel also renders light blue
ISLAND_OCEAN_COLOR = "#b0cfe8"    # poster set_bad colour
# Display only, does not touch the .npy files
ISLAND_CROSS_SHORE_MODES = ["trimmed", "padded"]
ISLAND_PAD_ROWS = 200               # cells of constant cross-shore extent when padded.
                                    # 200 = 2000 m, nothing cropped; 100 = 1000 m, warns if land is cropped
OFFSET_ROW_ORDER = "D1_first"  # "D1_first" | "D90_first"
OFFSET_SEAWARD_POSITIVE = True

# Island sections

# D1 = Cape Point (south) to D90 = Pea Island (north); [] disables the labels
SECTIONS = [
    ((1, 6),   "Cape Point"),
    ((7, 8),   "Buxton"),
    ((9, 20),  "Buxton-Avon"),
    ((21, 31), "Avon"),
    ((32, 67), "Avon-Tri-Village / Wimble Shoals"),
    ((68, 83), "Tri-Village"),
    ((84, 90), "Pea Island / N. Rodanthe"),
]

# Datums / thresholds
MHW_M = 0.36               # m NAVD88
BERM_ELEV_NAVD_M = 1.70    # m NAVD88 (matches the 1.7 m storm/collision threshold)
BEACH_START_THR_M = 0.50   # m MHW-relative, strict '>' comparison
WATER_CLAMP_M = -3.0       # m MHW-relative; below this -> sentinel.
                           # -3.0 keeps back-barrier marsh (Lexi's v3 edit); -1.0 was the original
SENTINEL_WATER_M = -3.0    # m MHW-relative
MIN_DUNE_H_M = 0.1         # m, floor on dune height above berm

# No data is not water

# The -10 m no-return value is tracked apart from the water clamp
RAW_NODATA_MAX_NAVD = -9.0   # raw NoData is exactly -10.0 m NAVD88
NODATA_SENTINEL_M = -99.0    # internal only; never written to the topography
WRITE_NODATA_MASK = True

# Geometry
TOPO_ROWS = 200            # max interior rows written
ALONG_COLS = 50            # alongshore profiles per domain (500 m / 10 m)
OCEAN_LOC = "right"        # "right", "left", "top", "bottom" in the RAW array

# True reverses alongshore order after orienting
ALONGSHORE_FLIP = True

# Road overlay

# NC-12 drawn on the picker and every per-domain figure
SHOW_ROAD = True
ROAD_YEARS = list(ROAD_LINE_VINTAGES)
REQUIRE_ROAD_MASKS = True   # missing file or shape mismatch -> hard error.
                            # fires only if the raster tree moved or a re-export changed a grid
from site_layer.hat_topo_version import ROAD_RASTER_ROOT  # noqa: E402
ROAD_MASK_DIR_FMT = "{year}/masks"
ROAD_MASK_NAME_FMT = "domain_{domain}_road_{year}.npy"

# D1-D7 (Cape Point) have ZERO road cells in both vintages -- NC-12 does not reach the point
ROAD_COLORS = {1978: "#6A1B9A", 2008: "#111111"}   # purple 1978 line, near-black 2008 line
ROAD_EDGE_COLOR = "#FFFFFF"   # thin outline so the road reads on dark water AND
                              # light land, which a single fill colour cannot
ROAD_PLAN_ALPHA = 0.55        # filled road cells on a map panel
ROAD_ENVELOPE_ALPHA = 0.10    # cross-shore envelope band on a profile panel
SHOW_ROAD_CENTER_LINE = True  # dashed alongshore line through the road centre

# Dune search
DEFAULT_WINDOW_PX = 8      # fallback window length (px landward of beach start)
CLIP_WINDOW_TO_BEACH = True  # per profile, start the search at max(i0, beach_start)
                             # so a window cannot wander onto the wet beach where the shoreline curves

# Interior
USE_CONST_INTERIOR = False  # interior starts one cell landward of the MOST landward
                           # dune in the domain, keeping alignment; False = each profile behind its own dune
FILL_MISSING_DUNE = True   # profiles with no dune found get MIN_DUNE_H_M instead of
                           # the -3.0 m sentinel (a negative dune height is not a valid Barrier3D input)
TRIM_INTERIOR_ROWS = True  # CHANGES THE .npy CASCADE READS, not just the figures.
                           # True = Lexi's v3, drop all-water rows; False = fixed (TOPO_ROWS, ALONG_COLS), as 2009_v1

# Straighten

# Shear each profile so the shoreline runs horizontally before the window is picked (why: README)
STRAIGHTEN = True

# What to align on: 'beach', the first cell above BEACH_START_THR_M
STRAIGHTEN_REF = "beach"

# 'linear' shears by a straight-line fit to start_beach; 'raw' by start_beach itself
STRAIGHTEN_FIT = "linear"

# Below this many profiles with a beach, skip straightening and say so.
STRAIGHTEN_MIN_PROFILES = 10
# -----------------------------------------------------------------------------


# Array helpers

# First/last (inclusive) columns that are not entirely water
def water_col_bounds(domain_array: np.ndarray, w_elev: float) -> tuple[int, int]:
    keep = [c for c in range(domain_array.shape[1])
            if not np.all(domain_array[:, c] <= w_elev + 1e-9)]
    if not keep:
        raise ValueError("domain is entirely water after clamping")
    return min(keep), max(keep)


# Trim leading/trailing rows that are entirely water
def remove_water_rows(domain_array: np.ndarray, w_elev: float) -> np.ndarray:
    keep = [r for r in range(domain_array.shape[0])
            if not np.all(domain_array[r, :] <= w_elev + 1e-9)]
    if not keep:
        return domain_array[:0, :]
    return domain_array[min(keep):max(keep) + 1, :]


# Return an (alongshore, cross_shore) array with the ocean in the LAST column
def orient_ocean_right(arr: np.ndarray, ocean_loc: str) -> np.ndarray:
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


# Shear each alongshore profile so the shoreline is horizontal
def straighten_profiles(z: np.ndarray, start_beach: np.ndarray):
    n_along, n_cross = z.shape
    shear = np.zeros(n_along, dtype=int)

    if not STRAIGHTEN:
        return z, start_beach, shear, 0.0

    sb = np.asarray(start_beach)
    ok = sb >= 0
    if int(ok.sum()) < STRAIGHTEN_MIN_PROFILES:
        print(f"       [warn] only {int(ok.sum())} profiles have a beach; "
              f"not straightening this domain")
        return z, start_beach, shear, 0.0

    x = np.arange(n_along)
    slope, intercept = np.polyfit(x[ok], sb[ok].astype(float), 1)
    obliq = float(np.degrees(np.arctan(abs(slope))))

    if STRAIGHTEN_FIT == "linear":
        ref = slope * x + intercept
    elif STRAIGHTEN_FIT == "raw":
        ref = sb.astype(float).copy()
        if (~ok).any():
            ref[~ok] = np.interp(x[~ok], x[ok], sb[ok].astype(float))
    else:
        raise ValueError(f"STRAIGHTEN_FIT must be linear/raw, "
                         f"got {STRAIGHTEN_FIT!r}")

    shear = np.round(ref - np.nanmin(ref)).astype(int)
    shear = np.clip(shear, 0, n_cross - 1)

    # Cells shifted in from beyond the array are no-data, not water
    zs = np.full_like(z, NODATA_SENTINEL_M)
    for i in range(n_along):
        k = int(shear[i])
        if k:
            zs[i, :n_cross - k] = z[i, k:]
        else:
            zs[i, :] = z[i, :]

    sb_new = np.where(ok, np.maximum(sb - shear, 0), -1)
    return zs, sb_new, shear, obliq


# Apply an existing shear to another array on the same grid
def shear_like(arr: np.ndarray, shear: np.ndarray) -> np.ndarray:
    out = np.zeros_like(arr)
    n_along, n_cross = arr.shape
    for i in range(min(n_along, len(shear))):
        k = int(shear[i])
        if k:
            out[i, :n_cross - k] = arr[i, k:]
        else:
            out[i, :] = arr[i, :]
    return out


# Put a raw GIS mask into the frame the topography was saved in
def align_mask_to_topography(raw_mask: np.ndarray, dom: dict) -> np.ndarray:
    mask = np.squeeze(np.asarray(raw_mask))
    if mask.ndim != 2:
        raise ValueError(f"expected a 2-D mask, got shape {mask.shape}")

    binary = np.isfinite(mask) & (mask > 0)

    oriented = orient_ocean_right(binary, OCEAN_LOC)
    ocean_first = np.ascontiguousarray(oriented[:, ::-1])

    # shear_like fills with zeros, so the mask must stay bool
    sheared = shear_like(ocean_first, dom["shear"])

    c0 = int(dom["c0"])
    n_cross = dom["z"].shape[1]
    trimmed = sheared[:, c0:c0 + n_cross]

    if trimmed.shape != dom["z"].shape:
        raise ValueError(
            f"aligned mask {trimmed.shape} does not match topography "
            f"{dom['z'].shape}; the mask was not exported on the same grid"
        )
    return np.ascontiguousarray(trimmed, dtype=bool)


# Source cross-shore cell that becomes SAVED interior row 0, per profile
def interior_row0_line(prof_arr: np.ndarray,
                       dune_loc: np.ndarray) -> tuple[np.ndarray, int]:
    topo, start_island = build_interior(prof_arr, dune_loc)

    if TRIM_INTERIOR_ROWS:
        keep = [r for r in range(topo.shape[0])
                if not np.all(topo[r, :] <= SENTINEL_WATER_M + 1e-9)]
        lead_trim = int(min(keep)) if keep else 0
    else:
        lead_trim = 0

    if start_island is not None:
        # USE_CONST_INTERIOR: a horizontal cut, the same row 0 on every profile
        row0 = np.full(len(dune_loc), start_island + lead_trim, dtype=int)
    else:
        row0 = np.where(dune_loc >= 0, dune_loc + 1 + lead_trim, -1)
    return row0, lead_trim


# Sort key: the domain number in a name, then the name
def natural_key(name: str) -> tuple:
    m = re.search(r"domain_(\d+)", name)
    return (int(m.group(1)) if m else 10**9, name)


# The domain number in a file stem, or None
def domain_number(stem: str) -> int | None:
    m = re.search(r"domain_(\d+)", stem)
    return int(m.group(1)) if m else None


# Island section label for a domain, per the D1=Cape Point convention
def section_for(stem: str) -> str:
    n = domain_number(stem)
    if n is None:
        return ""
    for (lo, hi), label in SECTIONS:
        if lo <= n <= hi:
            return label
    return ""


# Zero-padded figure name so 90+ domains sort correctly in Explorer
def fig_stem(stem: str) -> str:
    n = domain_number(stem)
    return f"domain_{n:03d}" if n is not None else stem


# A domain's figure title, with its island section
def fig_title(stem: str) -> str:
    sec = section_for(stem)
    return f"{stem}  ({sec})" if sec else stem


# Road masks: read-only consumers of HAT_rasterize_road_to_domains.py's masks

# Where HAT_rasterize_road_to_domains.py put one domain's mask
def road_mask_path(year: int, domain_id: int) -> Path:
    return (Path(ROAD_RASTER_ROOT) / ROAD_MASK_DIR_FMT.format(year=year)
            / ROAD_MASK_NAME_FMT.format(domain=domain_id, year=year))


# Load every ROAD_YEARS mask for one domain, in both frames
def load_road_masks(dem_path: Path, dom: dict) -> tuple[dict, dict]:
    if not SHOW_ROAD:
        return {}, {}

    domain_id = domain_number(Path(dem_path).stem)
    if domain_id is None:
        raise ValueError(f"cannot determine domain ID for road mask: {dem_path.name}")

    raw_out, aligned_out = {}, {}
    for year in ROAD_YEARS:
        path = road_mask_path(year, domain_id)
        if not path.is_file():
            msg = f"no {year} road mask for domain {domain_id}: {path}"
            if REQUIRE_ROAD_MASKS:
                raise FileNotFoundError(msg)
            print(f"       [road warn] {msg}")
            continue

        mask = np.squeeze(np.load(path))
        if mask.ndim != 2:
            raise ValueError(f"{path.name}: expected a 2-D mask, got {mask.shape}")
        # The mask is checked against the raw DEM shape; no transpose or resize
        if tuple(mask.shape) != tuple(dom["raw_shape_unoriented"]):
            raise ValueError(
                f"{path.name}: road shape {mask.shape} does not match DEM shape "
                f"{tuple(dom['raw_shape_unoriented'])}. Re-run "
                f"HAT_rasterize_road_to_domains.py for {year}."
            )

        binary = np.isfinite(mask) & (mask > 0)
        oriented = orient_ocean_right(binary, OCEAN_LOC)
        raw_out[year] = np.ascontiguousarray(oriented[:, ::-1], dtype=bool)
        aligned_out[year] = align_mask_to_topography(mask, dom)

    return raw_out, aligned_out


# Seaward edge, landward edge and centre cell of the road, per profile
def road_profile_positions(mask: np.ndarray | None) -> tuple[np.ndarray, np.ndarray,
                                                             np.ndarray]:
    if mask is None:
        return (np.array([], dtype=float),) * 3
    n_along = mask.shape[0]
    seaward = np.full(n_along, np.nan)
    landward = np.full(n_along, np.nan)
    center = np.full(n_along, np.nan)
    for i in range(n_along):
        cells = np.flatnonzero(mask[i])
        if cells.size:
            seaward[i] = float(cells.min())
            landward[i] = float(cells.max())
            center[i] = float(np.mean(cells))
    return seaward, landward, center


# Re-index road cells into the SAVED interior grid
def processed_road_grid(mask: np.ndarray | None, row0_line: np.ndarray,
                        n_rows: int, n_cols: int) -> np.ndarray:
    out = np.zeros((n_rows, n_cols), dtype=bool)
    if mask is None:
        return out
    n_fill = min(n_cols, mask.shape[0], len(row0_line))
    for i in range(n_fill):
        start = int(row0_line[i])
        if start < 0:
            continue
        dest = np.flatnonzero(mask[i]) - start
        keep = (dest >= 0) & (dest < n_rows)
        out[dest[keep], i] = True
    return out


# Per-domain road geometry and setback from SAVED interior row 0
def road_offset_stats(mask: np.ndarray | None,
                      row0_line: np.ndarray) -> dict:
    empty = {
        "road_profiles": 0, "road_cells": 0,
        "road_span_cells": np.nan, "road_center_cell": np.nan,
        "road_width_cells": np.nan,
        "setback_median_m": np.nan, "setback_mean_m": np.nan,
        "setback_min_m": np.nan, "setback_max_m": np.nan,
        "center_median_m": np.nan, "n_seaward": 0,
    }
    if mask is None or not np.any(mask):
        return empty

    seaward, landward, center = road_profile_positions(mask)

    # Cap at ALONG_COLS, as the setback script does
    n = min(len(center), len(row0_line), ALONG_COLS)
    row0 = np.asarray(row0_line[:n], dtype=float)
    valid = np.isfinite(center[:n]) & (row0 >= 0)

    out = dict(empty)
    out["road_profiles"] = int(valid.sum())
    out["road_cells"] = int(np.count_nonzero(mask[:n]))
    if np.isfinite(center[:n]).any():
        out["road_span_cells"] = float(np.nanmax(landward[:n])
                                      - np.nanmin(seaward[:n]) + 1)
        out["road_center_cell"] = float(np.nanmedian(center[:n]))

    if valid.any():
        sb = (seaward[:n][valid] - row0[valid]) * CELL_SIZE_M
        ctr = (center[:n][valid] - row0[valid]) * CELL_SIZE_M
        width = (landward[:n][valid] - seaward[:n][valid] + 1)
        out.update({
            "setback_median_m": float(np.median(sb)),
            "setback_mean_m": float(np.mean(sb)),
            "setback_min_m": float(np.min(sb)),
            "setback_max_m": float(np.max(sb)),
            "center_median_m": float(np.median(ctr)),
            "road_width_cells": float(np.median(width)),
            # Counted on the seaward edge, the value floored before the model
            "n_seaward": int(np.count_nonzero(sb < 0)),
        })
    return out


# Draw exact road cells on an (alongshore x, cross-shore y) map panel
def add_road_plan_overlay(ax, mask: np.ndarray | None, year: int,
                          *, zorder: float = 5.0, label: bool = True) -> None:
    if mask is None or not np.any(mask):
        return
    color = ROAD_COLORS.get(year, "#111111")
    n_along, n_cross = mask.shape
    display = np.ma.masked_where(~mask.T, np.ones((n_cross, n_along)))
    ax.imshow(display, aspect="auto", origin="lower",
              extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
              interpolation="nearest", cmap=ListedColormap([color]),
              vmin=0.0, vmax=1.0, alpha=ROAD_PLAN_ALPHA, zorder=zorder)
    if np.any(~mask):
        ax.contour(np.arange(n_along), np.arange(n_cross), mask.T.astype(float),
                   levels=[0.5], colors=ROAD_EDGE_COLOR, linewidths=0.7,
                   zorder=zorder + 0.2)
    if SHOW_ROAD_CENTER_LINE:
        _, _, center = road_profile_positions(mask)
        if np.isfinite(center).any():
            ax.plot(np.arange(n_along), center, color=color, lw=1.2, ls="--",
                    zorder=zorder + 0.4)
    if label:
        ax.plot([], [], color=color, lw=5, alpha=ROAD_PLAN_ALPHA,
                label=f"NC-12 {year}")


# Road's cross-shore envelope on a profile-stack panel (elevation x, cell y)
def add_road_envelope(ax, mask: np.ndarray | None, year: int,
                      *, label: bool = True) -> None:
    if mask is None or not np.any(mask):
        return
    color = ROAD_COLORS.get(year, "#111111")
    seaward, landward, center = road_profile_positions(mask)
    if not np.isfinite(center).any():
        return
    ax.axhspan(float(np.nanmin(seaward)) - 0.5, float(np.nanmax(landward)) + 0.5,
               color=color, alpha=ROAD_ENVELOPE_ALPHA, zorder=0,
               label=f"NC-12 {year} envelope" if label else None)
    ax.axhline(float(np.nanmedian(center)), color=color, lw=1.1, ls="--",
               zorder=1)


# Load / prep

# Load one domain
def load_profiles(in_path: Path) -> dict:
    arr = np.load(in_path).astype(float, copy=False)
    if arr.ndim != 2:
        raise ValueError(f"expected 2D array, got {arr.ndim}D")

    raw_shape_unoriented = arr.shape   # the shape the road rasterizer snapped to
    arr = orient_ocean_right(arr, OCEAN_LOC)

    # Sanity check that OCEAN_LOC is right (a warning, not an auto-orient)
    edge = max(1, min(5, arr.shape[1] // 20))
    left_q = np.nanpercentile(arr[:, :edge], 25)
    right_q = np.nanpercentile(arr[:, -edge:], 25)
    if right_q > left_q:
        print(f"[warn] {in_path.name}: after orienting, the LEFT edge is lower "
              f"(p25 left={left_q:.2f}, right={right_q:.2f}). Check OCEAN_LOC.")

    raw = np.ascontiguousarray(arr[:, ::-1])  # NAVD88, ocean first, untrimmed

    # No-data is identified on the RAW array, before the clamp folds it into the water sentinel
    nodata = raw <= RAW_NODATA_MAX_NAVD

    z = raw - MHW_M
    z[z < WATER_CLAMP_M] = SENTINEL_WATER_M
    z[nodata] = NODATA_SENTINEL_M

    n_along = z.shape[0]
    if n_along < ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} < {ALONG_COLS}; "
              f"trailing output cols remain sentinel.")
    elif n_along > ALONG_COLS:
        print(f"[warn] {in_path.name}: alongshore={n_along} > {ALONG_COLS}; "
              f"only first {ALONG_COLS} profiles used.")

    # Order matters: start_beach, then straighten, then water trim
    above = z > BEACH_START_THR_M
    start_beach = np.where(above.any(axis=1), above.argmax(axis=1), -1)

    z, start_beach, shear, obliq = straighten_profiles(z, start_beach)

    c0, c1 = water_col_bounds(z, SENTINEL_WATER_M)
    z = z[:, c0:c1 + 1]
    start_beach = np.where(start_beach >= 0,
                           np.maximum(start_beach - c0, 0), -1)

    dom = {"raw": raw, "z": z, "start_beach": start_beach, "c0": int(c0),
           "shear": shear, "obliquity_deg": round(float(obliq), 2),
           "name": in_path.name,
           "raw_shape_unoriented": raw_shape_unoriented}

    # LAST, because align_mask_to_topography needs the finished z, c0 and shear
    dom["road_raw"], dom["road_masks"] = load_road_masks(in_path, dom)
    return dom


# NaN out sentinel water cells for plotting/statistics
def masked_profiles(prof_arr: np.ndarray) -> np.ndarray:
    return np.where(prof_arr <= SENTINEL_WATER_M + 1e-9, np.nan, prof_arr)


# The fallback window: DEFAULT_WINDOW_PX landward of the median beach start
def default_window(prof_arr: np.ndarray, start_beach: np.ndarray) -> tuple[int, int]:
    valid = start_beach[start_beach >= 0]
    base = int(np.median(valid)) if valid.size else 0
    return base, min(base + DEFAULT_WINDOW_PX, prof_arr.shape[1])


# Suggested window

# How far landward to look for a foredune, and the margin around the crest
SUGGEST_REACH_PX = 26      # cells landward of beach start to hunt the crest in
SUGGEST_PAD_SEAWARD = 3    # cells kept seaward of the crest
SUGGEST_PAD_LANDWARD = 2   # cells kept landward of the crest (INSIDE the window)


# A crest-aware starting window
def suggest_window(prof_arr: np.ndarray, start_beach: np.ndarray,
                   road_seaward: int | None = None) -> tuple[int, int, float]:
    n_along, n_cross = prof_arr.shape
    limit = n_cross

    locs = []
    for i in range(n_along):
        a = int(start_beach[i]) if start_beach[i] >= 0 else 0
        b = min(a + SUGGEST_REACH_PX, limit)
        if b <= a:
            continue
        w = prof_arr[i, a:b]
        valid = w > SENTINEL_WATER_M + 1e-9
        if not valid.any():
            continue
        locs.append(a + int(np.argmax(np.where(valid, w, -np.inf))))

    if not locs:
        i0, i1 = default_window(prof_arr, start_beach)
        return i0, i1, float("nan"), False

    crest = int(np.median(locs))
    i0 = max(0, crest - SUGGEST_PAD_SEAWARD)
    i1 = min(limit, crest + SUGGEST_PAD_LANDWARD + 1)
    if i1 <= i0:
        i1 = min(limit, i0 + 1)

    elev, loc = find_dunes(prof_arr, start_beach, i0, i1)
    ok = loc >= 0
    crest_el = float(np.median(elev[ok])) if ok.any() else float("nan")
    overlaps = road_seaward is not None and i1 > road_seaward
    return i0, i1, crest_el, overlaps


# (median crest, % of profiles pinned at i1-1, higher cell just outside the landward edge?)
def window_diagnostics(prof_arr, start_beach, i0, i1):
    elev, loc = find_dunes(prof_arr, start_beach, i0, i1)
    ok = loc >= 0
    if not ok.any():
        return float("nan"), 0.0, 0.0
    pinned = float(np.mean(loc[ok] == i1 - 1))
    if i1 < prof_arr.shape[1]:
        nxt = prof_arr[:, i1]
        higher = float(np.mean(nxt[ok] > elev[ok]))
    else:
        higher = 0.0
    return float(np.median(elev[ok])), pinned, higher


# Interactive picker

# SpanSelector with a props/rectprops fallback for matplotlib < 3.5
def _span_selector(ax, on_select):
    try:
        return SpanSelector(ax, on_select, "vertical", useblit=True,
                            props=dict(alpha=0.2, facecolor="#FF8C00"))
    except TypeError:
        return SpanSelector(ax, on_select, "vertical", useblit=True,
                            rectprops=dict(alpha=0.2, facecolor="#FF8C00"))


# Show the domain ocean-at-bottom and let the user drag a dune search window
def pick_window(stem: str, prof_arr: np.ndarray, start_beach: np.ndarray,
                init: tuple[int, int],
                road_masks: dict | None = None) -> tuple[str, int, int]:
    n_along, n_cross = prof_arr.shape
    zm = masked_profiles(prof_arr)
    state = {"i0": init[0], "i1": init[1], "action": None}

    fig, (ax_map, ax_prof) = plt.subplots(
        1, 2, figsize=(13, 9), sharey=True,
        gridspec_kw={"width_ratios": [1.3, 1.0], "wspace": 0.06},
    )
    try:
        fig.canvas.manager.set_window_title(f"dune search window - {stem}")
    except Exception:
        pass

    finite = zm[np.isfinite(zm)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0

    # Map panel: alongshore on x, cross-shore on y, ocean at the bottom
    im = ax_map.imshow(
        np.ma.masked_invalid(zm.T), aspect="auto", origin="lower",
        extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
        cmap="terrain", vmin=-1.0, vmax=max(vmax, 2.0),
    )
    sb = np.where(start_beach >= 0, start_beach, np.nan)
    ax_map.plot(np.arange(n_along), sb, color="k", lw=1.2, label="beach start")
    for _yr, _m in (road_masks or {}).items():
        add_road_plan_overlay(ax_map, _m, _yr)
    ax_map.set_xlabel("alongshore cell")
    ax_map.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_map.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # Profile panel: elevation on x, cross-shore on y
    y = np.arange(n_cross)
    ax_prof.plot(zm.T, y, color="0.75", lw=0.6)
    med = np.nanmedian(zm, axis=0)
    ax_prof.plot(med, y, color="k", lw=2.0, label="median profile")
    ax_prof.axvline(BEACH_START_THR_M, color="#1565C0", ls="--", lw=1.0,
                    label=f"beach thr {BEACH_START_THR_M} m")
    ax_prof.axvline(BERM_ELEV_NAVD_M - MHW_M, color="#B71C1C", ls=":", lw=1.2,
                    label=f"berm {BERM_ELEV_NAVD_M} m NAVD88")
    for _yr, _m in (road_masks or {}).items():
        add_road_envelope(ax_prof, _m, _yr)
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

    spans = [None, None]
    sug_lines = [None, None]

    # The road's seaward edge, so the suggestion can never propose a window that reaches NC-12
    _road_sea = None
    if road_masks:
        cols = [np.flatnonzero(m.any(axis=0)) for m in road_masks.values()
                if m is not None and m.any()]
        if cols:
            _road_sea = int(min(c.min() for c in cols))
    sug_i0, sug_i1, sug_crest, sug_on_road = suggest_window(
        prof_arr, start_beach, _road_sea)

    def redraw():
        for k, ax in enumerate((ax_map, ax_prof)):
            if spans[k] is not None:
                spans[k].remove()
            # HALF-CELL OFFSETS, so the band covers the cells actually SEARCHED
            spans[k] = ax.axhspan(state["i0"] - 0.5, state["i1"] - 0.5,
                                  color="#FF8C00", alpha=0.25, zorder=0)
            # The suggestion, as an outline so it reads as a proposal rather than a second selection
            if sug_lines[k] is not None:
                for ln in sug_lines[k]:
                    ln.remove()
            sug_lines[k] = [
                ax.axhline(sug_i0 - 0.5, color="#1b6ca8", lw=1.4, ls=(0, (5, 3))),
                ax.axhline(sug_i1 - 0.5, color="#1b6ca8", lw=1.4, ls=(0, (5, 3))),
            ]

        crest, pinned, higher = window_diagnostics(
            prof_arr, start_beach, state["i0"], state["i1"])
        warn = ""
        if pinned >= 0.5 and higher >= 0.5:
            warn = ("   ⚠ CLIPPED: argmax pins at the landward edge on "
                    f"{pinned:.0%} of profiles and cell {state['i1']} is higher "
                    f"on {higher:.0%} — widen landward")
        elif pinned >= 0.5:
            warn = (f"   · argmax sits on the last cell on {pinned:.0%} of "
                    "profiles (fine if the crest really is there)")
        fig.suptitle(
            f"{stem}   window = [{state['i0']}, {state['i1']}) = cells "
            f"{state['i0']}-{state['i1'] - 1}   crest {crest:.2f} m{warn}\n"
            f"suggested [{sug_i0}, {sug_i1}) crest {sug_crest:.2f} m"
            f"{'  [overlaps NC-12 - check it is the dune, not the embankment]' if sug_on_road else ''}  (blue "
            f"dashes)   |   a = adopt suggestion   enter = accept   r = reset   "
            f"s = skip   q = quit",
            fontsize=10.5,
        )
        fig.canvas.draw_idle()

    redraw()

    def on_select(ymin, ymax):
        i0 = int(np.clip(round(ymin), 0, n_cross - 2))
        i1 = int(np.clip(round(ymax), i0 + 1, n_cross))
        state["i0"], state["i1"] = i0, i1
        redraw()

    def on_key(event):
        if event.key in ("enter", "return"):
            state["action"] = "accept"
            plt.close(fig)
        elif event.key == "a":
            state["i0"], state["i1"] = sug_i0, sug_i1
            redraw()
        elif event.key == "r":
            state["i0"], state["i1"] = init
            redraw()
        elif event.key == "s":
            state["action"] = "skip"
            plt.close(fig)
        elif event.key in ("q", "escape"):
            state["action"] = "quit"
            plt.close(fig)

    selectors = [_span_selector(ax_map, on_select), _span_selector(ax_prof, on_select)]

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()  # blocks
    del selectors

    if state["action"] is None:  # window closed with the X
        state["action"] = "accept"
    return state["action"], state["i0"], state["i1"]


# Extraction

# Dune elevation and location per profile within the picked window
def find_dunes(prof_arr, start_beach, i0, i1):
    n_along, n_cross = prof_arr.shape
    dune_elev = np.full(n_along, np.nan)
    dune_loc = np.full(n_along, -1, dtype=int)

    for i in range(n_along):
        a = i0
        if CLIP_WINDOW_TO_BEACH and start_beach[i] >= 0:
            a = max(i0, int(start_beach[i]))
        b = min(i1, n_cross)
        if b <= a:
            continue
        w = prof_arr[i, a:b]
        valid = w > SENTINEL_WATER_M + 1e-9
        if not valid.any():
            continue
        # Argmax inside the window, not over the whole profile
        k = int(np.argmax(np.where(valid, w, -np.inf)))
        dune_elev[i] = float(w[k])
        dune_loc[i] = a + k

    return dune_elev, dune_loc


# Return interior topography, (cross_shore_rows, alongshore_cols), ocean-first
def build_interior(prof_arr, dune_loc):
    n_along, n_cross = prof_arr.shape
    topo = np.full((TOPO_ROWS, ALONG_COLS), SENTINEL_WATER_M, dtype=float)
    n_fill = min(ALONG_COLS, n_along)

    if USE_CONST_INTERIOR:
        valid = dune_loc[dune_loc >= 0]
        if valid.size == 0:
            return topo, None
        start_island = int(valid.max()) + 1
        block = prof_arr[:n_fill, start_island:].T  # (cross, along)
        rows = min(TOPO_ROWS, block.shape[0])
        topo[:rows, :n_fill] = block[:rows, :]
        if block.shape[0] > TOPO_ROWS:
            print(f"       [warn] interior truncated: {block.shape[0]} -> {TOPO_ROWS} rows")
        return topo, start_island

    for i in range(n_fill):
        if dune_loc[i] < 0:
            continue
        use = prof_arr[i, dune_loc[i] + 1:-1]
        if use.size:
            rows = min(TOPO_ROWS, use.size)
            topo[:rows, i] = use[:rows]
    return topo, None


# 'domain_11' -> '11'
def _gis_id(stem: str) -> str:
    if not stem.startswith("domain_"):
        raise ValueError(
            f"unexpected array stem {stem!r} - expected 'domain_<gis id>'")
    return stem[len("domain_"):]


# One domain: find the dunes, build the interior, save both arrays and the figures
def extract_domain(stem, prof_arr, start_beach, i0, i1, topo_dir, dune_dir,
                   shear=None, obliquity_deg=0.0, road_masks=None):
    n_along, n_cross = prof_arr.shape
    dune_elev, dune_loc = find_dunes(prof_arr, start_beach, i0, i1)

    n_found = int(np.isfinite(dune_elev).sum())
    if n_found == 0:
        print(f"[skip] {stem}: no dune found in window [{i0}, {i1}]")
        return None

    dune_h = dune_elev - (BERM_ELEV_NAVD_M - MHW_M)  # both MHW-relative
    dune_h[np.isfinite(dune_h) & (dune_h < 0.0)] = MIN_DUNE_H_M

    dune_m = np.full(ALONG_COLS, SENTINEL_WATER_M, dtype=float)
    n_fill = min(ALONG_COLS, n_along)
    dune_m[:n_fill] = dune_h[:n_fill]
    missing = ~np.isfinite(dune_m[:n_fill])
    if missing.any():
        idx = np.where(missing)[0].tolist()
        if FILL_MISSING_DUNE:
            dune_m[:n_fill][missing] = MIN_DUNE_H_M
            print(f"       [warn] {stem}: no dune at profiles {idx} -> "
                  f"filled with {MIN_DUNE_H_M} m")
        else:
            dune_m[:n_fill][missing] = SENTINEL_WATER_M
            print(f"       [warn] {stem}: no dune at profiles {idx} -> sentinel")

    topo_m, start_island = build_interior(prof_arr, dune_loc)
    if TRIM_INTERIOR_ROWS:
        topo_m = remove_water_rows(topo_m, SENTINEL_WATER_M)

    # Where SAVED interior row 0 sits on each profile, and the road measured against it
    row0_line, lead_trim = interior_row0_line(prof_arr, dune_loc)
    road_stats = {yr: road_offset_stats(m, row0_line)
                  for yr, m in (road_masks or {}).items()}

    # Split no-data back out; the topography written is unchanged
    topo_nodata = topo_m <= NODATA_SENTINEL_M + 1e-9
    dune_nodata = dune_m <= NODATA_SENTINEL_M + 1e-9
    topo_m = np.where(topo_nodata, SENTINEL_WATER_M, topo_m)
    dune_m = np.where(dune_nodata, SENTINEL_WATER_M, dune_m)

    topo_dm = topo_m * 0.1
    dune_dm = dune_m * 0.1

    gid = _gis_id(stem)
    topo_out = Path(topo_dir) / array_name("topography", gid)
    dune_out = Path(dune_dir) / array_name("dune", gid)
    np.save(topo_out, topo_dm)
    np.save(dune_out, dune_dm)
    if WRITE_NODATA_MASK:
        np.save(Path(topo_dir) / array_name("nodata", gid), topo_nodata)
        if topo_nodata.any():
            print(f"       [nodata] {stem}: {int(topo_nodata.sum())} of "
                  f"{topo_nodata.size} interior cells never surveyed "
                  f"({topo_nodata.mean():.1%})")
    print(f"[ok] {stem}: window [{i0}, {i1}], {n_found}/{n_along} dunes, "
          f"interior {topo_dm.shape}, mean dune h = "
          f"{np.nanmean(dune_h):.2f} m")

    return {
        "stem": stem, "i0": i0, "i1": i1, "n_found": n_found, "n_along": n_along,
        "mean_dune_h_m": float(np.nanmean(dune_h)),
        "min_dune_h_m": float(np.nanmin(dune_h)),
        "max_dune_h_m": float(np.nanmax(dune_h)),
        "n_filled": int(missing.sum()),
        "interior_rows": int(topo_dm.shape[0]),
        "interior_cols": int(topo_dm.shape[1]),
        "mean_interior_elev_m": float(np.mean(topo_m[topo_m > SENTINEL_WATER_M + 1e-9]))
        if np.any(topo_m > SENTINEL_WATER_M + 1e-9) else np.nan,
        "start_island": start_island,
        # The saved arrays are in the straightened frame
        "straightened": bool(STRAIGHTEN),
        "obliquity_deg": obliquity_deg,
        "shear_max_cells": int(np.max(shear)) if shear is not None else 0,
        "shear": shear,
        "row0_line": row0_line, "lead_trim_rows": lead_trim,
        "road_stats": road_stats,
        "dune_elev": dune_elev, "dune_loc": dune_loc, "dune_h": dune_h,
        "topo_dm": topo_dm, "dune_dm": dune_dm,
        "topo_file": topo_out.name, "dune_file": dune_out.name,
    }


# Qc figure

# QC figure in the same ocean-at-bottom orientation as the picker
def qc_figure(stem, prof_arr, start_beach, res, fig_dir: Path,
              road_masks: dict | None = None):
    n_along, n_cross = prof_arr.shape
    zm = masked_profiles(prof_arr)
    i0, i1 = res["i0"], res["i1"]
    dl = np.where(res["dune_loc"] >= 0, res["dune_loc"], np.nan)

    fig = plt.figure(figsize=(12.5, 9.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1.0], height_ratios=[1.0, 0.35],
                          wspace=0.06, hspace=0.12)
    ax_map = fig.add_subplot(gs[0, 0])
    ax_prof = fig.add_subplot(gs[0, 1], sharey=ax_map)
    ax_h = fig.add_subplot(gs[1, 0], sharex=ax_map)

    finite = zm[np.isfinite(zm)]
    vmax = float(np.percentile(finite, 99)) if finite.size else 3.0

    # Map: alongshore x, cross-shore y (ocean at bottom)
    im = ax_map.imshow(np.ma.masked_invalid(zm.T), aspect="auto", origin="lower",
                       extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
                       cmap="terrain", vmin=-1.0, vmax=max(vmax, 2.0))
    ax_map.axhspan(i0, i1, color="#FF8C00", alpha=0.20, zorder=0)
    ax_map.plot(np.arange(n_along), np.where(start_beach >= 0, start_beach, np.nan),
                color="k", lw=1.0, label="beach start")
    ax_map.plot(np.arange(n_along), dl, color="#B71C1C", lw=0, marker=".", ms=4,
                label="dune crest")
    if res["start_island"] is not None:
        ax_map.axhline(res["start_island"], color="w", lw=1.2, ls="--",
                       label="interior start")
    for _yr, _m in (road_masks or {}).items():
        add_road_plan_overlay(ax_map, _m, _yr)
    ax_map.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_map.legend(loc="upper right", fontsize=8, framealpha=0.9)
    plt.setp(ax_map.get_xticklabels(), visible=False)

    # Profiles: elevation x, cross-shore y
    y = np.arange(n_cross)
    ax_prof.plot(zm.T, y, color="0.75", lw=0.5)
    ax_prof.plot(np.nanmedian(zm, axis=0), y, color="k", lw=2.0)
    ax_prof.axhspan(i0, i1, color="#FF8C00", alpha=0.20, zorder=0)
    ax_prof.axvline(BERM_ELEV_NAVD_M - MHW_M, color="#B71C1C", ls=":", lw=1.2)
    ax_prof.plot(res["dune_elev"], dl, lw=0, marker=".", ms=4, color="#B71C1C")
    for _yr, _m in (road_masks or {}).items():
        add_road_envelope(ax_prof, _m, _yr, label=False)
    ax_prof.set_xlabel("elev (m MHW)")
    ax_prof.set_xlim(-1.2, max(vmax, 2.0) + 0.5)
    ax_prof.set_ylim(-0.5, n_cross - 0.5)
    plt.setp(ax_prof.get_yticklabels(), visible=False)

    # Dune height alongshore, aligned under the map
    ax_h.plot(np.arange(n_along), res["dune_h"], color="#FF8C00", lw=1.5,
              marker="o", ms=3)
    ax_h.axhline(MIN_DUNE_H_M, color="0.5", ls="--", lw=1.0)
    ax_h.set_xlabel("alongshore cell")
    ax_h.set_ylabel("dune height\nabove berm (m)")
    ax_h.set_xlim(-0.5, n_along - 0.5)

    try:
        fig.colorbar(im, ax=[ax_map, ax_prof], location="right",
                     fraction=0.035, pad=0.02, label="elev (m MHW)")
    except (TypeError, ValueError):
        fig.colorbar(im, ax=ax_prof, label="elev (m MHW)")

    road_note = ""
    for _yr in ROAD_YEARS:
        st = (res.get("road_stats") or {}).get(_yr)
        if st and np.isfinite(st["setback_median_m"]):
            road_note += (f"   |  NC-12 {_yr}: {st['setback_median_m']:+.0f} m "
                          f"from interior row 0")
        elif st is not None:
            road_note += f"   |  NC-12 {_yr}: not in domain"
    fig.suptitle(f"{fig_title(stem)}  |  search window [{i0}, {i1}]  |  "
                 f"{res['n_found']}/{n_along} dunes found{road_note}",
                 fontsize=12)

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"{fig_stem(stem)}_qc.png", dpi=130,
                bbox_inches="tight")
    plt.close(fig)


# GIS vs processed comparison figure

# The whole chain, left to right
def comparison_figure(stem, dom, res, fig_dir: Path):
    raw = dom["raw"]
    zs = dom["z"]
    c0 = dom["c0"]
    n_along, n_raw = raw.shape
    n_cross = zs.shape[1]
    i0, i1 = res["i0"], res["i1"]
    xs = np.arange(n_along)

    shear = np.asarray(dom.get("shear", np.zeros(n_along, dtype=int)))
    off = c0 + shear.astype(float)          # processed cell k -> raw column k+off
    obliq = dom.get("obliquity_deg", 0.0)

    dl = np.where(res["dune_loc"] >= 0, res["dune_loc"].astype(float), np.nan)
    sb = np.where(dom["start_beach"] >= 0, dom["start_beach"].astype(float),
                  np.nan)

    topo_m = res["topo_dm"] * 10.0          # back to m MHW
    dune_m = res["dune_dm"] * 10.0          # m above berm
    topo_disp = np.where(topo_m <= SENTINEL_WATER_M + 1e-9, np.nan, topo_m)
    n_int_rows, n_int_cols = topo_disp.shape

    fig = plt.figure(figsize=(19, 9.6), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.08, h_pad=0.05, wspace=0.04, hspace=0.02)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 1.0],
                          height_ratios=[1.0, 0.14])
    ax_raw = fig.add_subplot(gs[:, 0])
    ax_str = fig.add_subplot(gs[:, 1])
    ax_int = fig.add_subplot(gs[0, 2])
    ax_dun = fig.add_subplot(gs[1, 2])

    # Shared elevation span so the three panels are comparable by eye -- just labelled NAVD88 vs MHW
    land = raw[raw > MHW_M]
    hi_navd = float(np.percentile(land, 99)) if land.size else MHW_M + 4.0
    vmin_r = MHW_M - 1.0
    vmax_r = max(hi_navd, MHW_M + 2.0)

    # 1: raw GIS
    im_r = ax_raw.imshow(raw.T, aspect="auto", origin="lower",
                         extent=[-0.5, n_along - 0.5, -0.5, n_raw - 0.5],
                         cmap="terrain", vmin=vmin_r, vmax=vmax_r)
    cb_r = fig.colorbar(im_r, ax=ax_raw, location="right", pad=0.02,
                        fraction=0.045, aspect=32)
    cb_r.set_label("elev (m NAVD88)", fontsize=13)
    cb_r.ax.tick_params(labelsize=12)
    ax_raw.contour(xs, np.arange(n_raw), raw.T, levels=[MHW_M],
                   colors="#1565C0", linewidths=1.0)
    ax_raw.fill_between(xs, i0 + off, i1 + off, color="#FF8C00", alpha=0.22,
                        zorder=0, label="dune search window")
    ax_raw.plot(xs, sb + off, color="k", lw=1.0, label="beach start")
    ax_raw.plot(xs, dl + off, color="#B71C1C", lw=0, marker=".", ms=4,
                label="picked dune crest")
    if res["start_island"] is not None:
        ax_raw.plot(xs, res["start_island"] + off, color="w", lw=1.4, ls="--",
                    label="interior start")
    ax_raw.plot(xs, off - 0.5, color="0.3", lw=1.0, ls="-.",
                label="water trim edge")

    # The road in the raw frame: NC-12's real diagonal, before the shear
    for _yr, _m in (dom.get("road_raw") or {}).items():
        add_road_plan_overlay(ax_raw, _m, _yr)

    # Crop to the island so the diagonal is legible next to the straightened panel
    lo_r = max(float(np.nanmin(off)) - 4.0, -0.5)
    hi_r = min(float(np.nanmax(off)) + n_cross + 4.0, n_raw - 0.5)
    ax_raw.set_ylim(lo_r, hi_r)
    ax_raw.set_xlabel("alongshore cell")
    ax_raw.set_ylabel("cross-shore cell  (0 = ocean, landward up)")
    ax_raw.set_title(f"1. RAW GIS input   {raw.shape}   m NAVD88\n"
                     f"{dom['name']}", fontsize=13.5)
    ax_raw.legend(loc="upper right", fontsize=11, framealpha=0.92)

    # 2: straightened (the picking frame)
    zm = masked_profiles(zs)
    im_s = ax_str.imshow(np.ma.masked_invalid(zm.T), aspect="auto",
                         origin="lower",
                         extent=[-0.5, n_along - 0.5, -0.5, n_cross - 0.5],
                         cmap="terrain", vmin=vmin_r - MHW_M,
                         vmax=vmax_r - MHW_M)
    cb_s = fig.colorbar(im_s, ax=ax_str, pad=0.02, fraction=0.045, aspect=32)
    cb_s.set_label("elev (m MHW)", fontsize=13)
    cb_s.ax.tick_params(labelsize=12)
    ax_str.axhspan(i0, i1, color="#FF8C00", alpha=0.22, zorder=0,
                   label="dune search window")
    ax_str.plot(xs, sb, color="k", lw=1.0, label="beach start")
    ax_str.plot(xs, dl, color="#B71C1C", lw=0, marker=".", ms=4,
                label="picked dune crest")
    if res["start_island"] is not None:
        ax_str.axhline(res["start_island"], color="w", lw=1.4, ls="--",
                       label="interior start")
    for _yr, _m in (dom.get("road_masks") or {}).items():
        add_road_plan_overlay(ax_str, _m, _yr)
    ax_str.set_xlabel("alongshore cell")
    ax_str.set_ylabel("cross-shore cell  (straightened frame)")
    if STRAIGHTEN:
        t2 = (f"2. STRAIGHTENED — the picking frame   {zs.shape}\n"
              f"{STRAIGHTEN_FIT} fit on {STRAIGHTEN_REF} start  |  "
              f"obliquity {obliq:.1f}°  |  shear 0–{int(np.max(shear))} cells")
    else:
        t2 = (f"2. PROCESSED PROFILES — the picking frame   {zs.shape}\n"
              f"STRAIGHTEN = False, so this is the raw frame minus the trim")
    ax_str.set_title(t2, fontsize=13.5)
    ax_str.legend(loc="upper right", fontsize=11, framealpha=0.92)

    # 3: what CASCADE reads
    im_i = ax_int.imshow(np.ma.masked_invalid(topo_disp), aspect="auto",
                         origin="lower",
                         extent=[-0.5, n_int_cols - 0.5, -0.5, n_int_rows - 0.5],
                         cmap="terrain", vmin=vmin_r - MHW_M, vmax=vmax_r - MHW_M)
    cb_i = fig.colorbar(im_i, ax=ax_int, pad=0.02, fraction=0.045, aspect=32)
    cb_i.set_label("elev (m MHW)", fontsize=13)
    cb_i.ax.tick_params(labelsize=12)

    # The road re-indexed into the SAVED interior grid
    _row0 = np.asarray(res.get("row0_line", []), dtype=int)
    for _yr, _m in (dom.get("road_masks") or {}).items():
        if _row0.size:
            add_road_plan_overlay(
                ax_int,
                processed_road_grid(_m, _row0, n_int_rows, n_int_cols).T,
                _yr, label=False)
    ax_int.set_ylabel("interior row\n(0 = ocean side)")
    ax_int.set_title(f"3. PROCESSED CASCADE input\ninterior {res['topo_file']}  "
                     f"{res['topo_dm'].shape}  (dam on disk)", fontsize=13.5)
    plt.setp(ax_int.get_xticklabels(), visible=False)

    dune_disp = np.where(dune_m <= SENTINEL_WATER_M + 1e-9, np.nan, dune_m)
    im_d = ax_dun.imshow(np.ma.masked_invalid(dune_disp[None, :]), aspect="auto",
                         origin="lower", extent=[-0.5, ALONG_COLS - 0.5, -0.5, 0.5],
                         cmap="YlOrBr", vmin=0.0,
                         vmax=max(float(np.nanmax(dune_disp)), 0.5))
    # A taller colourbar for the short dune strip
    cb_d = fig.colorbar(im_d, ax=ax_dun, pad=0.02, fraction=0.045, aspect=5)
    cb_d.set_label("dune h above\nberm (m)", fontsize=13)
    cb_d.ax.tick_params(labelsize=12)
    ax_dun.set_yticks([0])
    ax_dun.set_yticklabels(["dune"])
    ax_dun.set_xlabel("alongshore cell")

    interior = ("sheared per profile" if res["start_island"] is None
                else f"const from row {res['start_island']}")
    road_note = ""
    for _yr in ROAD_YEARS:
        st = (res.get("road_stats") or {}).get(_yr)
        if st and np.isfinite(st["setback_median_m"]):
            road_note += f"   |   NC-12 {_yr} {st['setback_median_m']:+.0f} m"
    fig.suptitle(
        f"{fig_title(stem)}  —  {RUN_NAME}\n"
        f"window [{i0}, {i1}] ({i1 - i0} cells)   |   "
        f"{res['n_found']}/{res['n_along']} dunes   |   "
        f"mean dune h {res['mean_dune_h_m']:.2f} m   |   "
        f"obliquity {obliq:.1f}°   |   interior {interior}{road_note}",
        fontsize=14, linespacing=1.5,
    )

    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"{fig_stem(stem)}_gis_vs_processed.png", dpi=150,
                bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


# Settings sheet & json I/O

# Every global knob, for the config sheet / manifest / provenance
def global_config() -> dict:
    return {
        "run name": RUN_NAME,
        "version": VERSION,
        "run_time": datetime.now().isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "root GIS domains": LOAD_PATH.name,
        "run folder": str(RUN_DIR),
        "load_path_full": str(LOAD_PATH),
        "topo_path_full": str(TOPO_SAVE_PATH),
        "dune_path_full": str(DUNE_SAVE_PATH),
        "window_json": str(WINDOW_JSON),
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
        "shift interior": False,  # not implemented in this version (see header)
        "constant interior row": USE_CONST_INTERIOR,
        "FILL_MISSING_DUNE": FILL_MISSING_DUNE,
        "STRAIGHTEN": STRAIGHTEN,
        "STRAIGHTEN_REF": STRAIGHTEN_REF if STRAIGHTEN else "",
        "STRAIGHTEN_FIT": STRAIGHTEN_FIT if STRAIGHTEN else "",
        "TRIM_INTERIOR_ROWS": TRIM_INTERIOR_ROWS,
        # Road overlay settings, recorded; none affect the saved arrays
        "SHOW_ROAD": SHOW_ROAD,
        "ROAD_YEARS": str(ROAD_YEARS) if SHOW_ROAD else "",
        "ROAD_RASTER_ROOT": str(ROAD_RASTER_ROOT) if SHOW_ROAD else "",
        "REQUIRE_ROAD_MASKS": REQUIRE_ROAD_MASKS if SHOW_ROAD else "",
    }


# Per-year NC-12 columns for the settings sheet
def road_columns(res: dict) -> dict:
    def num(v, nd=1):
        return round(v, nd) if np.isfinite(v) else ""

    out = {}
    for year in ROAD_YEARS:
        st = (res.get("road_stats") or {}).get(year)
        if st is None:
            continue
        out[f"road profiles {year}"] = st["road_profiles"]
        out[f"road cells {year}"] = st["road_cells"]
        out[f"road span cells {year}"] = num(st["road_span_cells"])
        out[f"road width cells {year}"] = num(st["road_width_cells"])
        out[f"road center cell {year}"] = num(st["road_center_cell"])
        # the comparable one, named for what it is
        out[f"road setback {year} (m)"] = num(st["setback_median_m"])
        out[f"road setback mean {year} (m)"] = num(st["setback_mean_m"])
        out[f"road setback min {year} (m)"] = num(st["setback_min_m"])
        out[f"road setback max {year} (m)"] = num(st["setback_max_m"])
        # Centre-referenced, for continuity with the roya-style dune-to-road number
        out[f"road center offset {year} (m)"] = num(st["center_median_m"])
        # Profiles where the road is SEAWARD of interior row 0
        out[f"road seaward profiles {year}"] = st["n_seaward"]
    return out


# One sheet row per domain
def settings_row(res: dict, w: dict | None) -> dict:
    return {
        "domain": domain_number(res["stem"]),
        "stem": res["stem"],
        "section": section_for(res["stem"]),
        # Settings, mirroring the tracking sheet
        "root GIS domains": LOAD_PATH.name,
        "root dunes/topo": RUN_NAME,
        "figure name": f"{fig_stem(res['stem'])}_gis_vs_processed.png",
        "beach start": BEACH_START_THR_M,
        "dune window start": res["i0"],
        "dune window end": res["i1"],
        "dune window": res["i1"] - res["i0"],
        "window source": "picked" if w else "default",
        "clip window to beach": CLIP_WINDOW_TO_BEACH,
        "shift interior": False,
        "constant interior row": USE_CONST_INTERIOR,
        "straightened": res.get("straightened", ""),
        "straighten fit": STRAIGHTEN_FIT if STRAIGHTEN else "",
        "obliquity (deg)": res.get("obliquity_deg", ""),
        "shear max (cells)": res.get("shear_max_cells", ""),
        "shear max (m)": (round(res.get("shear_max_cells", 0) * CELL_SIZE_M, 1)
                          if res.get("shear_max_cells") else ""),
        "MHW (m NAVD88)": MHW_M,
        "berm (m NAVD88)": BERM_ELEV_NAVD_M,
        "water clamp (m MHW)": WATER_CLAMP_M,
        # Results
        "interior start row": res["start_island"],
        "lead trim rows": res.get("lead_trim_rows", ""),
        "dunes found": res["n_found"],
        "profiles": res["n_along"],
        "dunes filled": res["n_filled"],
        "mean dune h (m)": round(res["mean_dune_h_m"], 3),
        "min dune h (m)": round(res["min_dune_h_m"], 3),
        "max dune h (m)": round(res["max_dune_h_m"], 3),
        "interior rows": res["interior_rows"],
        "interior cols": res["interior_cols"],
        "mean interior elev (m MHW)": round(res["mean_interior_elev_m"], 3),
        "topo file": res["topo_file"],
        "dune file": res["dune_file"],
        "picked": (w or {}).get("picked", ""),
        **road_columns(res),
    }


# Write the per-domain settings sheet as CSV, plus XLSX if pandas is around
def write_settings_sheet(rows: list, base_path: Path) -> None:
    if not rows:
        return
    base_path.parent.mkdir(parents=True, exist_ok=True)

    # Fieldnames from every row, so one absent mask cannot crash the write
    fieldnames = list(dict.fromkeys(k for r in rows for k in r))
    csv_path = base_path.with_suffix(".csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, restval="")
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


# Plain-text record so a run folder found later explains itself
def write_manifest(rows: list, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = global_config()
    width = max(len(k) for k in cfg)
    lines = [
        "=" * 78,
        f"HAT dune & topography extraction  --  run {RUN_NAME}",
        "=" * 78,
        "",
        "CONTENTS",
        "  topography\\            interior elevation arrays (dam), CASCADE input",
        "  dunes\\                 dune height above berm (dam), CASCADE input",
        "  figures\\gis_vs_processed\\   raw DEM vs what CASCADE receives, per domain",
        "  figures\\qc\\            dune detection QC, per domain",
        f"  {SHEET_SAVE_PATH.name}.xlsx / .csv   per-domain settings + results",
        f"  {SUMMARY_FIG_PATH.name}   all domains on one page",
        f"  {PLAN_FIG_STEM}_<year>_<mode>.png   plan view, per year per mode",
        "",
        "PICKS (not in this folder -- not regenerable)",
        f"  {WINDOW_JSON}",
    ]
    if PICK_SET == RUN_NAME:
        lines += [
            "  ^ THIS VERSION CARRIES ITS OWN PICKS. It was seeded by copying the",
            "    shared straight set, then re-picked with NC-12 drawn on the",
            "    picker, so the shared set is untouched and earlier versions stay",
            "    reproducible.",
        ]
    else:
        lines += ["  ^ shared across versions -- re-picking here rewrites it for "
                  "every version that points at it"]
    if SHOW_ROAD:
        lines += [
            "",
            "ROAD OVERLAY (display + diagnostics only)",
            f"  masks read from : {ROAD_RASTER_ROOT}",
            f"  vintages        : {', '.join(str(y) for y in ROAD_YEARS)}",
            "  The road does NOT enter find_dunes, build_interior or",
            "  straighten_profiles. topography\\, dunes\\ and the nodata masks are",
            "  byte-identical with SHOW_ROAD either way; only the figures and the",
            "  road columns of the settings sheet change.",
            "  Road distances are metres from SAVED interior row 0, positive =",
            "  landward, matching setback_dunestart_m in RoadOffset_<year>_domains.csv.",
        ]
    lines += [
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


# Every domain on one page
def summary_figure(rows: list, path: Path) -> None:
    if not rows:
        return
    rows = sorted([r for r in rows if r["domain"] is not None],
                  key=lambda r: r["domain"])
    if not rows:
        return
    d = np.array([r["domain"] for r in rows])

    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(13, 9), sharex=True,
                                        gridspec_kw={"hspace": 0.12})

    # section bands + labels, matching the HAT_hindcast annotation style
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
                     label="dune search window")
    ax0.step(d, [r["interior start row"] for r in rows], where="mid",
             color="#B71C1C", lw=1.4, label="interior start row")

    # NC-12's median cross-shore cell per domain, on the same axis as the search window
    for _yr in ROAD_YEARS:
        key = f"road center cell {_yr}"
        if not any(isinstance(r.get(key), (int, float)) for r in rows):
            continue
        vals = np.array([r.get(key) if isinstance(r.get(key), (int, float))
                         else np.nan for r in rows], dtype=float)
        if np.isfinite(vals).any():
            ax0.plot(d, vals, color=ROAD_COLORS.get(_yr, "#111111"), lw=1.2,
                     ls="--", label=f"NC-12 {_yr} centre")
    ax0.set_ylabel("cross-shore cell\n(0 = ocean)")
    ax0.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax0.set_title(f"HAT dune & topography extraction  —  {RUN_NAME}  —  "
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
    ax2.set_xlabel("GIS domain (south → north)")

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


# Return the raw_offset CSV for a year
def resolve_offset_file(year: int, configured) -> Path | None:
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


# Preflight every path before doing any work
def check_paths() -> bool:
    print("=" * 78)
    print(f"PATH CHECK  —  run {RUN_NAME}")
    print("=" * 78)
    ok = True

    def row(tag, label, path, note=""):
        print(f"  [{tag:^7}] {label:<20} {path}{note}")

    for label, pth in (("project root", PROJECT_ROOT), ("hatteras_init", INIT_ROOT)):
        if Path(pth).is_dir():
            row("ok", label, pth)
        else:
            row("MISSING", label, pth)
            ok = False

    if Path(LOAD_PATH).is_dir():
        n = len([f for f in os.listdir(LOAD_PATH)
                 if f.startswith("domain_") and f.endswith(".npy")])
        if n:
            row("ok", "DEM domains", LOAD_PATH, f"   ({n} domain_*.npy)")
        else:
            row("EMPTY", "DEM domains", LOAD_PATH, "   no domain_*.npy here")
            ok = False
    else:
        row("MISSING", "DEM domains", LOAD_PATH)
        ok = False

    for year in sorted(OFFSET_FILES):
        found = resolve_offset_file(year, OFFSET_FILES[year])
        if found is None:
            row("MISSING", f"offsets {year}", OFFSET_FILES[year],
                "   (figures will skip this year)")
        elif Path(found) != Path(OFFSET_FILES[year]):
            row("FOUND", f"offsets {year}", found, "   (not where configured)")
        else:
            row("ok", f"offsets {year}", found)

    if SHOW_ROAD:
        for year in ROAD_YEARS:
            mdir = Path(ROAD_RASTER_ROOT) / ROAD_MASK_DIR_FMT.format(year=year)
            if not mdir.is_dir():
                row("MISSING", f"road masks {year}", mdir,
                    "   run HAT_rasterize_road_to_domains.py")
                if REQUIRE_ROAD_MASKS:
                    ok = False
                continue
            wanted = PICK_DOMAINS if PICK_DOMAINS is not None \
                else list(range(1, NUM_REAL_DOMAINS + 1))
            absent = [n for n in wanted if not road_mask_path(year, n).is_file()]
            if absent:
                row("MISSING" if REQUIRE_ROAD_MASKS else "WARN",
                    f"road masks {year}", mdir,
                    f"   {len(wanted) - len(absent)}/{len(wanted)} present, "
                    f"missing {absent[:8]}{' ...' if len(absent) > 8 else ''}")
                if REQUIRE_ROAD_MASKS:
                    ok = False
            else:
                row("ok", f"road masks {year}", mdir,
                    f"   ({len(wanted)} masks)")

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

    # Picking against a shared pick set overwrites other versions' picks: warn
    if MODE in ("pick", "pick_and_run") and PICK_SET != RUN_NAME:
        row("WARN", "pick set", wj,
            f"\n            ^ MODE = {MODE!r} will REWRITE this shared pick set "
            f"({PICK_SET!r}).\n"
            f"              Every version that reads it changes with it. To give "
            f"{RUN_NAME} its own,\n"
            f"              set PICK_SET = RUN_NAME and copy the old file to "
            f"HAT_dune_search_windows_{RUN_NAME}.json first.")

    if MODE in ("pick", "pick_and_run") and not REPICK_EXISTING and wj.exists():
        try:
            n_saved = len([k for k in json.loads(wj.read_text())
                           if not k.startswith("_")])
            n_want = (len(PICK_DOMAINS) if PICK_DOMAINS is not None else n_saved)
            if n_saved >= n_want:
                row("WARN", "pick pass", wj,
                    f"\n            ^ REPICK_EXISTING = False and all {n_want} "
                    f"requested domains are already saved,\n"
                    f"              so the pick pass will skip every one of them. "
                    f"Set REPICK_EXISTING = True to re-draw.")
        except Exception:
            pass

    for legacy in (Path(RUN_DIR) / f"{PLAN_FIG_STEM}.png",):
        if legacy.exists():
            row("STALE", "original plan view", legacy,
                "\n            ^ pre-dates the per-year/per-mode naming, nothing "
                "overwrites it — delete it")

    rd = Path(RUN_DIR)
    row("exists" if rd.exists() else "new", "run folder", rd,
        "   (will be overwritten)" if rd.exists() else "   (created on run)")

    print("=" * 78)
    if not ok:
        print("[stop] required inputs missing — fix the paths above and re-run.\n")
    return ok


# {year
def load_offsets() -> dict:
    out = {}
    for year, configured in OFFSET_FILES.items():
        path = resolve_offset_file(year, configured)
        if path is None:
            print(f"[warn] raw_offset file not found, skipping {year}: {configured}")
            continue
        v = np.loadtxt(path, skiprows=1, delimiter=",", ndmin=2).astype(float)
        if v.shape[1] > 1:
            print(f"[raw_offset] {year}: {v.shape[1]} columns, using column {OFFSET_COLUMN}")
            v = v[:, OFFSET_COLUMN]
        else:
            v = v[:, 0]
        padded = NUM_REAL_DOMAINS + 2 * N_BUFFER_DOMAINS
        if v.size == padded:
            print(f"[raw_offset] {year}: padded file ({padded} rows), stripping "
                  f"{N_BUFFER_DOMAINS} buffer domains from each end")
            v = v[N_BUFFER_DOMAINS:N_BUFFER_DOMAINS + NUM_REAL_DOMAINS]
        n = v.size
        if OFFSET_ROW_ORDER == "D1_first":
            dom = np.arange(n) + 1
        else:
            dom = n - np.arange(n)
        order = np.argsort(dom)
        out[year] = (dom[order], v[order])
        print(f"[raw_offset] {year}: {n} domains, {v.min():.0f}-{v.max():.0f} m "
              f"({path.name})")
    return out


# Terrain with 0 m pinned to colormap position 0.35
def _island_norm():
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


# Warn if the assembled alongshore axis is discontinuous at domain seams
def _assert_alongshore_continuity(grids, label: str, warn_ratio: float = 5.0):
    series = []
    for g in grids:
        land = (g > SENTINEL_WATER_M + 1e-9).sum(axis=0).astype(float)
        if land.size:
            series.append(land)
    if len(series) < 3:
        return float("nan")

    width = min(s.size for s in series)
    stack = np.array([s[:width] for s in series])
    seam = np.abs(stack[:-1, -1] - stack[1:, 0])
    inner = np.abs(np.diff(stack, axis=1))
    if not inner.size or inner.mean() <= 0:
        return float("nan")

    ratio = float(seam.mean() / inner.mean())
    if ratio > warn_ratio:
        print(f"       [ALONGSHORE WARNING] {label}: seam/inner discontinuity "
              f"ratio {ratio:.1f} (expected ~1-2).")
        print(f"       The alongshore axis is reversed WITHIN each domain "
              f"relative to the domain order, so every 500 m block is backwards.")
        print(f"       Either the arrays were extracted with ALONGSHORE_FLIP = "
              f"False, or a per-domain np.fliplr has been reintroduced here.")
    return ratio


# Stitch processed domains onto one plan-view canvas at their dune offsets
def _build_island_canvas(recs, offset_m_by_domain, mode):
    use = [(n, r) for n, r in recs if n in offset_m_by_domain]
    if not use:
        return None, None, None, None
    off_cells = [int(round(offset_m_by_domain[n] / CELL_SIZE_M)) for n, _ in use]

    grids, dunes, cropped = [], [], []
    roads = {yr: [] for yr in ROAD_YEARS}
    for n_dom, r in use:
        g = r["topo_dm"] * CELL_SIZE_M                 # dam -> m MHW

        # The road in the saved interior frame, padded and cropped with the grid
        _row0 = np.asarray(r.get("row0_line", []), dtype=int)
        for yr in ROAD_YEARS:
            m = (r.get("road_masks") or {}).get(yr)
            roads[yr].append(
                processed_road_grid(m, _row0, g.shape[0], g.shape[1])
                if (m is not None and _row0.size)
                else np.zeros(g.shape, dtype=bool))

        if mode == "padded":
            # Pad landward, matching where dune_topo_extractor_from_GIS.py left its sentinel
            n = g.shape[0]
            if n < ISLAND_PAD_ROWS:
                g = np.vstack([g, np.full((ISLAND_PAD_ROWS - n, g.shape[1]),
                                          SENTINEL_WATER_M)])
                for yr in ROAD_YEARS:
                    roads[yr][-1] = np.vstack([
                        roads[yr][-1],
                        np.zeros((ISLAND_PAD_ROWS - n, g.shape[1]), dtype=bool)])
            elif n > ISLAND_PAD_ROWS:
                n_land = int(np.sum(g[ISLAND_PAD_ROWS:] > SENTINEL_WATER_M + 1e-9))
                if n_land:
                    cropped.append((n_dom, n, n_land))
                g = g[:ISLAND_PAD_ROWS]
                # The road can be cropped off entirely here
                for yr in ROAD_YEARS:
                    roads[yr][-1] = roads[yr][-1][:ISLAND_PAD_ROWS]
        if ISLAND_SENTINEL_AS_OCEAN:
            g = np.where(g <= SENTINEL_WATER_M + 1e-9, np.nan, g)
        d = r["dune_dm"] * CELL_SIZE_M + (BERM_ELEV_NAVD_M - MHW_M)   # -> m MHW
        # No per-domain flip: the arrays already run south to north
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

    _assert_alongshore_continuity(grids, f"planview canvas ({mode})")

    max_rows = max(g.shape[0] for g in grids)
    canvas_rows = max(off_cells) + max_rows + 5
    total_cols = sum(g.shape[1] for g in grids)
    canvas = np.full((canvas_rows, total_cols), np.nan)

    road_canvas = {yr: np.zeros((canvas_rows, total_cols), dtype=bool)
                   for yr in ROAD_YEARS}

    col = 0
    starts = []
    for k, (g, d) in enumerate(zip(grids, dunes)):
        n_rows, n_cols = g.shape
        starts.append(col)
        origin = off_cells[k]
        end = min(origin + n_rows, canvas_rows)
        canvas[origin:end, col:col + n_cols] = g[:end - origin, :]
        if ISLAND_INCLUDE_DUNE and origin >= 1:
            canvas[origin - 1, col:col + min(n_cols, d.size)] = d[:n_cols]
        # The road placed by the topography's origin, columns and clip
        for yr in ROAD_YEARS:
            rg = roads[yr][k]
            if rg.shape[0] >= end - origin:
                road_canvas[yr][origin:end, col:col + n_cols] = rg[:end - origin, :]
        col += n_cols

    return canvas, np.array(starts), [n for n, _ in use], road_canvas


# Plan view of the processed dune and interior for domains 1-90 at the measured offsets, one per mode
def island_plan_figure(summary: list, offsets: dict, run_dir: Path) -> None:
    recs = sorted([(domain_number(r["stem"]), r) for r in summary
                   if domain_number(r["stem"]) is not None])
    if not recs or not offsets:
        return
    cmap, norm = _island_norm()

    if PRODUCT_YEAR is None:
        plot_years = sorted(offsets)
        print(f"[planview] {TOPO_PRODUCT}: no period year, plotting all "
              f"offset years {plot_years}")
    elif PRODUCT_YEAR in offsets:
        plot_years = [PRODUCT_YEAR]
        skipped = [y for y in sorted(offsets) if y != PRODUCT_YEAR]
        if skipped:
            print(f"[planview] {TOPO_PRODUCT}: plotting {PRODUCT_YEAR} offsets "
                  f"only; {skipped} loaded but not plotted (this topography is "
                  f"not that year's island)")
    else:
        print(f"[planview] {TOPO_PRODUCT}: {PRODUCT_YEAR} offsets not loaded "
              f"(have {sorted(offsets)}), no plan view written")
        return

    for year in plot_years:
        dom, v = offsets[year]
        omap = {int(a): float(b) for a, b in zip(dom, v)}
        for mode in ISLAND_CROSS_SHORE_MODES:
            canvas, starts, used, road_canvas = _build_island_canvas(
                recs, omap, mode)
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

            # NC-12 across the whole island, in the offset frame
            for _yr in ROAD_YEARS:
                rc = (road_canvas or {}).get(_yr)
                if rc is None or not rc.any():
                    continue
                ax.pcolormesh(np.ma.masked_where(~rc, rc.astype(float)),
                              cmap=ListedColormap([ROAD_COLORS.get(_yr, "#111111")]),
                              vmin=0.0, vmax=1.0, shading="auto",
                              rasterized=True, zorder=4)
                ax.plot([], [], color=ROAD_COLORS.get(_yr, "#111111"), lw=3,
                        label=f"NC-12 {_yr}")
            if any((road_canvas or {}).get(y) is not None
                   and (road_canvas or {})[y].any() for y in ROAD_YEARS):
                ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

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
                if n % 5 == 0 or n == 1:
                    ticks.append(starts[k] + ALONG_COLS // 2)
                    labels.append(str(n))
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels, fontsize=9)
            ax.set_xlabel("GIS domain (south → north)", fontsize=12,
                          labelpad=8)
            ax.set_ylabel("Cross-shore cell (raw_offset frame)", fontsize=12)
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
            ax.set_title(f"Hatteras Island — CASCADE Initialization  |  {year} offsets  "
                         f"|  {DEM_LABEL} extracted {note}  ({len(used)} domains)",
                         fontsize=14, fontweight="bold", color="#1a1a2e", pad=12)

            path = Path(run_dir) / f"{PLAN_FIG_STEM}_{year}_{mode}.png"
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
            plt.close(fig)
            print(f"[planview] {path}")


# All domains together
def island_figure(summary: list, offsets: dict, path: Path) -> None:
    recs = sorted([(domain_number(r["stem"]), r) for r in summary
                   if domain_number(r["stem"]) is not None])
    if not recs:
        return
    d = np.array([n for n, _ in recs])

    pos_mean, pos_lo, pos_hi = [], [], []
    for _, r in recs:
        loc = r["dune_loc"].astype(float)
        loc[loc < 0] = np.nan
        # Back to the raw cross-shore axis
        _sh = np.asarray(r.get("shear", 0))
        pm = (loc + r["c0"] + _sh) * CELL_SIZE_M
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

    # 1) Measured offsets, all domains: both years drawn, the product's year heavy
    colors = {1984: "#1565C0", 2004: "#B71C1C"}
    for year in sorted(offsets):
        dom, v = offsets[year]
        own = (PRODUCT_YEAR is None) or (year == PRODUCT_YEAR)
        ax0.plot(dom, v, color=colors.get(year, "0.4"),
                 lw=2.2 if own else 1.1,
                 ls="-" if own else (0, (5, 3)),
                 alpha=1.0 if own else 0.55,
                 zorder=3 if own else 2,
                 label=(f"measured dune raw_offset {year}" if own else
                        f"measured dune raw_offset {year}  (reference)"))
    if PRODUCT_YEAR is not None and PRODUCT_YEAR in offsets:
        ax0.text(0.005, 0.04, f"{TOPO_PRODUCT} initialises at {PRODUCT_YEAR}",
                 transform=ax0.transAxes, fontsize=8, va="bottom", ha="left",
                 color=colors.get(PRODUCT_YEAR, "0.2"), fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.75",
                           alpha=0.9))
    ax0.set_ylabel("measured raw_offset\n(m, common frame)")
    ax0.legend(loc="upper right", fontsize=8, framealpha=0.9)
    sign = "seaward +" if OFFSET_SEAWARD_POSITIVE else "landward +"
    ax0.set_title(f"Island dune offsets vs extracted crest  —  "
                  f"{TOPO_PRODUCT} {RUN_NAME}  "
                  f"({len(recs)} domains, row0={OFFSET_ROW_ORDER}, {sign})",
                  fontsize=12)

    # 2) crest extracted from the DEM this run, raw-array frame
    ax1.fill_between(d, pos_lo, pos_hi, color="#FF8C00", alpha=0.25,
                     label="within-domain min-max")
    ax1.plot(d, pos_mean, color="#FF8C00", lw=1.6, marker="o", ms=3,
             label=f"extracted crest, {DEM_LABEL}")
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
        tag = "" if year != PRODUCT_YEAR else "  <- this product"
        notes.append(f"r(extracted, {year}) = {r:+.2f}  (n={ok.sum()}){tag}")
    if notes:
        trans = blended_transform_factory(ax1.transAxes, ax1.transAxes)
        ax1.text(0.01, 0.06, "   |   ".join(notes), transform=trans, fontsize=8,
                 va="bottom", ha="left",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7",
                           alpha=0.9))

    # 3) measured change
    if 1984 in offsets and 2004 in offsets:
        dom_a, va = offsets[1984]
        dom_b, vb = offsets[2004]
        common = np.intersect1d(dom_a, dom_b)
        ch = (vb[np.isin(dom_b, common)] - va[np.isin(dom_a, common)]) / 20.0
        ax2.axhline(0, color="0.5", lw=1.0)
        ax2.plot(common, ch, color="#B71C1C", lw=1.4, marker="o", ms=3)
        ax2.set_ylabel("measured change\n1984→2004 (m/yr)")
    else:
        ax2.text(0.5, 0.5, "need both 1984 and 2004 offsets", ha="center",
                 va="center", transform=ax2.transAxes, fontsize=10, color="0.4")
    ax2.set_xlabel("GIS domain (south → north)")
    ax2.set_xlim(min(d.min(), 1) - 0.5, max(d.max(), 90) + 0.5)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[island] {path}")
    for n in notes:
        print(f"         {n}")


# The saved picks, or an empty set
def load_windows(path: Path) -> dict:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {"_meta": {}}


# Write the picks with the settings they were made under
def save_windows(path: Path, windows: dict) -> None:
    windows["_meta"] = {
        "updated": datetime.now().isoformat(timespec="seconds"),
        "load_path": str(LOAD_PATH),
        "mhw_m": MHW_M,
        "berm_elev_navd_m": BERM_ELEV_NAVD_M,
        "beach_start_thr_m": BEACH_START_THR_M,
        "water_clamp_m": WATER_CLAMP_M,
        "ocean_loc": OCEAN_LOC,
        "alongshore_flip": ALONGSHORE_FLIP,
        "straighten": STRAIGHTEN,
        "straighten_fit": STRAIGHTEN_FIT if STRAIGHTEN else "",
        "note": "i0/i1 are cross-shore indices in the ocean-first, "
                "water-trimmed profile array. If straighten is true they are "
                "in the STRAIGHTENED frame and do not apply to an "
                "unstraightened array -- the index range is still valid, it "
                "just points at different cells. alongshore_flip is recorded "
                "for provenance only: the flip reverses the alongshore axis "
                "but leaves every profile's cross-shore frame bit-identical "
                "(shear[i] -> shear[n-1-i], c0 unchanged), so these windows "
                "are valid under either setting. Verified for all 90 domains "
                "by HAT_alongshore_frame_check.py.",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(windows, f, indent=2, sort_keys=True)


# Run: preflight the paths, pick and/or extract every domain, then the sheet and figures
def main():
    print("=" * 86)
    print(f"HAT dune / topography extraction  |  {RUN_NAME}  |  MODE = {MODE}")
    print("=" * 86)
    if STRAIGHTEN:
        print(f"  STRAIGHTEN = True ({STRAIGHTEN_FIT} fit on "
              f"{STRAIGHTEN_REF} start)")
        print(f"    Each profile is sheared so the shoreline is horizontal")
        print(f"    BEFORE the window is picked. Windows should come out narrow")
        print(f"    (~4 cells) instead of spanning the diagonal (~12).")
        print(f"    Picks are frame-dependent: this WINDOW_JSON must have been")
        print(f"    picked with STRAIGHTEN=True or the run pass will refuse.")
        print(f"      {WINDOW_JSON}")
    else:
        print(f"  STRAIGHTEN = False -- the shoreline crosses each domain")
        print(f"    diagonally, so windows must span it and dune picks get")
        print(f"    loose in the high-obliquity domains.")
    if USE_CONST_INTERIOR:
        print(f"  USE_CONST_INTERIOR = True -- the interior is cut horizontally")
        print(f"    at max(dune_loc)+1, which throws away everything seaward of")
        print(f"    that on every other profile. Straightening shrinks the loss")
        print(f"    but does not remove it; residual dune variability is real.")
        print(f"    Consider False.")
    else:
        print(f"  USE_CONST_INTERIOR = False -- interior sheared per profile,")
        print(f"    no wedge lost. InteriorDomain[0] is each profile's dune.")
    if SHOW_ROAD:
        print(f"  SHOW_ROAD = True -- NC-12 "
              f"({', '.join(str(y) for y in ROAD_YEARS)}) drawn on the picker and")
        print(f"    every per-domain figure, in the SAME sheared/trimmed frame as")
        print(f"    the topography. Display and diagnostics only: the arrays")
        print(f"    CASCADE reads are byte-identical with SHOW_ROAD off.")
        print(f"    A search window sitting landward of the road is picking the")
        print(f"    road embankment, not a dune crest.")
    if MODE in ("pick", "pick_and_run"):
        print(f"  PICKS -> {WINDOW_JSON.name}")
        print(f"    rewritten after every domain. REPICK_EXISTING="
              f"{REPICK_EXISTING}.")
    print()
    if not check_paths():
        return
    load_dir = Path(LOAD_PATH)
    topo_dir = Path(TOPO_SAVE_PATH)
    dune_dir = Path(DUNE_SAVE_PATH)
    topo_dir.mkdir(parents=True, exist_ok=True)
    dune_dir.mkdir(parents=True, exist_ok=True)

    names = sorted(
        [n for n in os.listdir(load_dir)
         if n.endswith(".npy") and n.startswith("domain_")],
        key=natural_key,
    )
    if PICK_DOMAINS is not None:
        wanted = set(PICK_DOMAINS)
        names = [n for n in names if domain_number(Path(n).stem) in wanted]
    print(f"[info] Found {len(names)} domain file(s) in {load_dir}")

    windows = load_windows(Path(WINDOW_JSON))

    # Pick pass
    if MODE in ("pick", "pick_and_run"):
        for name in names:
            stem = Path(name).stem
            if not REPICK_EXISTING and stem in windows:
                w = windows[stem]
                if bool(w.get("straightened", False)) != bool(STRAIGHTEN):
                    print(f"[repick] {stem}: saved window is from the "
                          f"STRAIGHTEN={bool(w.get('straightened', False))} "
                          f"frame, this run is STRAIGHTEN={STRAIGHTEN} -- "
                          f"re-picking")
                else:
                    print(f"[skip pick] {stem}: window [{w['i0']}, {w['i1']}] "
                          f"already saved")
                    continue
            try:
                dom = load_profiles(load_dir / name)
            except Exception as e:
                print(f"[skip] {name}: {e}")
                continue

            # Open on the SAVED window when there is one, so a re-pick is an adjustment
            saved = windows.get(stem)
            if saved and bool(saved.get("straightened", False)) == bool(STRAIGHTEN):
                init = (int(saved["i0"]), int(saved["i1"]))
            else:
                init = default_window(dom["z"], dom["start_beach"])
            action, i0, i1 = pick_window(stem, dom["z"], dom["start_beach"], init,
                                         road_masks=dom.get("road_masks"))

            if action == "quit":
                print("[info] quit picking; saved windows kept")
                break
            if action == "skip":
                print(f"[pick] {stem}: skipped -> default window {init}")
                continue

            windows[stem] = {
                "i0": int(i0), "i1": int(i1),
                "n_cross_trimmed": int(dom["z"].shape[1]),
                "n_along": int(dom["z"].shape[0]),
                "trim_offset_c0": int(dom["c0"]),
                # The frame this window was picked in
                "straightened": bool(STRAIGHTEN),
                "obliquity_deg": dom["obliquity_deg"],
                "shear_max_cells": int(np.max(dom["shear"])),
                "picked": datetime.now().isoformat(timespec="seconds"),
                # What this window replaced, so the v3 -> v4 change is auditable
                "prev_i0": init[0], "prev_i1": init[1],
                "changed": bool((int(i0), int(i1)) != (init[0], init[1])),
            }
            save_windows(Path(WINDOW_JSON), windows)  # save every domain
            moved = "" if (int(i0), int(i1)) == (init[0], init[1]) else \
                f"  (was [{init[0]}, {init[1]}])"
            print(f"[pick] {stem}: window [{i0}, {i1}] saved{moved}")

    # Run pass
    if MODE in ("run", "pick_and_run"):
        summary, rows = [], []
        for name in names:
            stem = Path(name).stem
            try:
                dom = load_profiles(load_dir / name)
            except Exception as e:
                print(f"[skip] {name}: {e}")
                continue
            prof_arr, start_beach = dom["z"], dom["start_beach"]

            w = windows.get(stem)
            if w is None:
                i0, i1 = default_window(prof_arr, start_beach)
                print(f"[warn] {stem}: no picked window, using default [{i0}, {i1}]")
            else:
                # Not 'unknown, proceed': that would apply an unstraightened window to a straightened array
                w_str = bool(w.get("straightened", False))
                if w_str != bool(STRAIGHTEN):
                    print(f"[skip] {stem}: window was picked with "
                          f"STRAIGHTEN={w_str} but this run has "
                          f"STRAIGHTEN={STRAIGHTEN}. The frames differ by the "
                          f"shear, so the window points at different cells. "
                          f"Re-pick, or point WINDOW_JSON at the matching set.")
                    continue
                if w.get("n_cross_trimmed") != prof_arr.shape[1]:
                    print(f"[warn] {stem}: trimmed width changed "
                          f"({w.get('n_cross_trimmed')} -> {prof_arr.shape[1]}); "
                          f"the saved window may no longer line up. Re-pick.")
                i0, i1 = int(w["i0"]), int(w["i1"])

            res = extract_domain(stem, prof_arr, start_beach, i0, i1,
                                 topo_dir, dune_dir,
                                 shear=dom["shear"],
                                 obliquity_deg=dom["obliquity_deg"],
                                 road_masks=dom.get("road_masks"))
            if res is None:
                continue
            res["c0"] = dom["c0"]
            # Carried for island_plan_figure, which needs the road in the saved interior frame
            res["road_masks"] = dom.get("road_masks")
            if SAVE_QC_FIGS:
                qc_figure(stem, prof_arr, start_beach, res, Path(FIG_DIR_QC),
                          road_masks=dom.get("road_masks"))
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

        print("\n" + "=" * 86)
        print(f"{'domain':<12}{'section':<26}{'window':>12}{'dunes':>8}"
              f"{'mean h (m)':>12}{'interior rows':>16}")
        print("-" * 86)
        for r in rows:
            win = "[{},{}]".format(r["dune window start"], r["dune window end"])
            print(f"{r['stem']:<12}{r['section'][:25]:<26}{win:>12}"
                  f"{r['dunes found']:>8}{r['mean dune h (m)']:>12.2f}"
                  f"{r['interior rows']:>16}")
        print("=" * 86)

        n_default = sum(1 for r in rows if r["window source"] == "default")
        n_filled = sum(1 for r in rows if r["dunes filled"] > 0)
        if n_default:
            print(f"[!] {n_default} domain(s) used the fallback default window "
                  f"— re-run MODE='pick' for those")
        if n_filled:
            print(f"[!] {n_filled} domain(s) had profiles with no dune found")

        if SHOW_ROAD:
            for year in ROAD_YEARS:
                have = [r for r in rows if r.get(f"road profiles {year}", 0)]
                seaward = [r["stem"] for r in rows
                           if isinstance(r.get(f"road seaward profiles {year}"), int)
                           and r[f"road seaward profiles {year}"] > 0]
                print(f"[road {year}] NC-12 present in {len(have)}/{len(rows)} "
                      f"domain(s)")
                if seaward:
                    print(f"[road {year}] road SEAWARD of interior row 0 on some "
                          f"profiles in {len(seaward)} domain(s): "
                          f"{', '.join(seaward[:10])}"
                          f"{' ...' if len(seaward) > 10 else ''}")
                    print(f"            Barrier3D cannot place a road at a "
                          f"negative setback; see RoadOffset_dunestart_audit.md.")
        print(f"\n[done] {len(summary)} domain(s). Everything for this run is in:"
              f"\n       {RUN_DIR}\n"
              f"       start with RUN_MANIFEST.txt and "
              f"{SUMMARY_FIG_PATH.name}\n"
              f"       CASCADE inputs: {topo_dir}\n"
              f"                       {dune_dir}")


if __name__ == "__main__":
    main()
