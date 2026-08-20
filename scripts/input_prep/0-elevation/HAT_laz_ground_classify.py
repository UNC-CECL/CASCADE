"""
HAT_laz_ground_classify.py

Step 0 of 3: derives a bare-earth ground classification for the 2008 NOAA IOCM
point cloud, which ships 100% unclassified.

WHY THIS EXISTS
---------------
The delivered tiles are entirely Classification 1 (verified: 117,849,675 of
117,849,675 points in tile R0). There is no class 2, so any "ground returns"
filter against the raw tiles selects nothing. HAT_dem_gap_fill.py needs bare
earth, so it has to be produced here first.

The cloud is also topographic only - z bottoms out at -1.33 m NAVD88, with no
submerged returns. It cannot extend the surface into open water; it can only
recover the intertidal/marsh fringe between the 2008 and 2009 waterlines. That
is the whole point of the fill, and it is why the fill in step 1 is limited to
cells where this cloud actually has returns.

WHAT IT DOES
------------
TWO PHASES, AND WHY IT HAS TO BE TWO
-------------------------------------
PHASE 1 (once per tile)   crop to the domain footprint -> split into 2 km cells
                          with a 100 m overlap buffer -> write raw cell files.
                          No filters, so no spatial index is ever built.
PHASE 2 (once per cell)   elm -> outlier -> smrf -> keep class 2 -> ground LAZ.
                          A separate PDAL PROCESS per cell.

The obvious one-pipeline version of this does not work, and failed here before
the split into phases. PDAL holds every PointView for the life of the run and
accumulates each filtered output until the writer executes at the very end, so
peak memory is the whole cropped cloud PLUS all outputs PLUS the current index -
it only ever grows. Measured on tile R0: 7.9 GB at 4.5 min, 14.1 GB at 27 min,
killed at 15.7 GB at 29 min with 1.2 GB free, having written nothing.

Putting phase 2 in its own process per cell is what actually bounds it: the OS
reclaims everything when each process exits. It also makes the run restartable
(finished cells are skipped) and gives real progress, neither of which the
single-pipeline version could do since it wrote only at the end.

Phase 1 still holds the cropped cloud in memory, but with no index and no
filtered copy it is far smaller than the peak above.

The 100 m buffer gives SMRF real context across cell boundaries, at the price of
duplicated points in the overlaps. That is deliberate and harmless downstream:
HAT_dem_gap_fill.py grids these points by taking the MINIMUM z per 1 m cell, and
a duplicated point contributes the same value twice.

SMRF PARAMETERS - THE ONE THING TO SANITY CHECK
-----------------------------------------------
SMRF removes objects that are smaller than `window` and steeper than `slope`.
On a barrier island the risk is not vegetation, it is DUNES: they are genuine
ground, but they are also compact and steep, so an aggressive window will shave
their crests and quietly lower the thing the model cares most about.

Defaults below are set wide (window 16 m, slope 0.2) to protect dune crests at
the cost of leaving some buildings in. Verify with the printed dune-crest check
before trusting the output, and see the QC note at the bottom.

REQUIRES
--------
PDAL, installed in its own conda env so it does not disturb the project Python:
    mamba create -y -n pdal -c conda-forge pdal
This script shells out to the pdal executable; it does not import it.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import geopandas as gpd

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[3]
GIS_ROOT = Path(r"D:\Hatteras_GIS")

POINT_DIR = GIS_ROOT / "Elevation" / "Points" / "2008" / "2008_NOAA_IOCM_NC_VA_J1437738"
TILE_GLOB = "*.laz"

# Ground-only output. Kept beside the source on D: - these are hundreds of MB
# and do not belong in the repo.
GROUND_DIR = GIS_ROOT / "Elevation" / "Points" / "2008" / "ground_smrf"

# Phase 1 output: raw cropped cells, unclassified. Kept rather than deleted so
# a failed or interrupted phase 2 can resume without re-reading the tiles.
CELL_DIR = GROUND_DIR / "_cells"

DOMAIN_FILE = GIS_ROOT / "domains.geojson"

# Crop margin around the domain footprint. Wide enough that SMRF's window has
# real context at the domain edge rather than classifying against a cut edge.
CROP_BUFFER_M = 200.0

# --- SMRF ---
SMRF_CELL = 1.0        # matches the 1 m DEM grid
SMRF_SLOPE = 0.20      # rise/run; raise to keep steeper dune faces as ground
SMRF_WINDOW = 16.0     # m; objects wider than this are kept as ground
SMRF_THRESHOLD = 0.45  # m elevation threshold
SMRF_SCALAR = 1.25

# Statistical outlier removal, applied before SMRF.
#
# This is the most memory-hungry stage after the point data itself: it builds a
# KD-tree over every point in a view. On the first full run, peak RSS climbed
# ~2 GB/min through the filter phase and blew past 14 GB. If that happens,
# turn this OFF first - filters.elm already removes the LOW outliers that
# actually corrupt a ground surface, and SMRF's own threshold/scalar absorb
# much of the rest. Dropping it costs far less accuracy than banding costs time
# (banding re-reads the whole tile per band).
USE_OUTLIER_FILTER = True
OUTLIER_MEAN_K = 8
OUTLIER_MULTIPLIER = 3.0

# --- CHUNKING (see WHY THE SPLITTER above) ---
SPLIT_LENGTH_M = 2000.0   # edge of each processing cell
SPLIT_BUFFER_M = 100.0    # overlap, for SMRF context across cell edges

# Hard sanity gate. Anything outside this is noise, not Hatteras.
Z_MIN, Z_MAX = -5.0, 60.0

PDAL_ENV_NAME = "pdal"


# =============================================================================
# PDAL DISCOVERY
# =============================================================================

def find_pdal():
    """pdal on PATH, else the conda env we told the user to create."""
    exe = shutil.which("pdal")
    if exe:
        return exe
    for base in (Path.home() / "miniforge3", Path.home() / "miniconda3",
                 Path.home() / "anaconda3"):
        cand = base / "envs" / PDAL_ENV_NAME / "Library" / "bin" / "pdal.exe"
        if cand.exists():
            return str(cand)
        cand = base / "envs" / PDAL_ENV_NAME / "bin" / "pdal"
        if cand.exists():
            return str(cand)
    raise FileNotFoundError(
        "pdal executable not found. Create it with:\n"
        f"    mamba create -y -n {PDAL_ENV_NAME} -c conda-forge pdal"
    )


# =============================================================================
# PIPELINE
# =============================================================================

def domain_crop_bounds(domain_file, buffer_m):
    gdf = gpd.read_file(domain_file)
    minx, miny, maxx, maxy = gdf.total_bounds
    return (minx - buffer_m, maxx + buffer_m,
            miny - buffer_m, maxy + buffer_m), len(gdf), gdf.crs


def split_pipeline(in_laz, out_template, bounds):
    """PHASE 1: crop + split only. No filters, so no index is built."""
    minx, maxx, miny, maxy = bounds
    return {
        "pipeline": [
            {"type": "readers.las", "filename": str(in_laz)},
            {"type": "filters.crop",
             "bounds": f"([{minx},{maxx}],[{miny},{maxy}])"},
            {"type": "filters.range", "limits": f"Z[{Z_MIN}:{Z_MAX}]"},
            {"type": "filters.splitter",
             "length": SPLIT_LENGTH_M,
             "buffer": SPLIT_BUFFER_M},
            {"type": "writers.las",
             "filename": str(out_template),
             "compression": "laszip",
             "forward": "all"},
        ]
    }


def ground_pipeline(in_cell, out_cell):
    """PHASE 2: one cell, one process. Memory is reclaimed when it exits."""
    stages = [
        {"type": "readers.las", "filename": str(in_cell)},
        {"type": "filters.elm"},
    ]
    if USE_OUTLIER_FILTER:
        stages.append(
            {"type": "filters.outlier",
             "method": "statistical",
             "mean_k": OUTLIER_MEAN_K,
             "multiplier": OUTLIER_MULTIPLIER})
    return {
        "pipeline": stages + [
            {"type": "filters.smrf",
             "cell": SMRF_CELL,
             "slope": SMRF_SLOPE,
             "window": SMRF_WINDOW,
             "threshold": SMRF_THRESHOLD,
             "scalar": SMRF_SCALAR,
             # do not let flagged noise anchor the ground surface
             "ignore": "Classification[7:7]"},
            {"type": "filters.range", "limits": "Classification[2:2]"},
            {"type": "writers.las",
             "filename": str(out_cell),
             "compression": "laszip",
             "forward": "all"},
        ]
    }


def run_pipeline(pdal_exe, pipeline, tag, quiet=False):
    pipe_dir = GROUND_DIR / "_pipelines"
    pipe_dir.mkdir(parents=True, exist_ok=True)
    pipe_path = pipe_dir / f"{tag}.json"
    pipe_path.write_text(json.dumps(pipeline, indent=2))
    cmd = [pdal_exe, "pipeline", str(pipe_path)]
    if not quiet:
        print(f"  $ {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise RuntimeError(f"exit {proc.returncode}: "
                           f"{err[-1] if err else 'no output'}")
    return proc


def count_points(laz_path):
    import laspy
    with laspy.open(str(laz_path)) as fh:
        h = fh.header
        return h.point_count, (h.mins[2], h.maxs[2])


def count_many(paths):
    """Total points and overall z range across a tile's output cells."""
    total, zlo, zhi = 0, float("inf"), float("-inf")
    for p in paths:
        n, (lo, hi) = count_points(p)
        total += n
        zlo, zhi = min(zlo, lo), max(zhi, hi)
    if total == 0:
        return 0, (float("nan"), float("nan"))
    return total, (zlo, zhi)


# =============================================================================
# MAIN
# =============================================================================

def main():
    pdal_exe = find_pdal()
    print(f"PDAL: {pdal_exe}")
    ver = subprocess.run([pdal_exe, "--version"], capture_output=True, text=True)
    print(f"  {ver.stdout.strip().splitlines()[0] if ver.stdout else '?'}")

    if not POINT_DIR.exists():
        raise FileNotFoundError(f"point cloud folder not found: {POINT_DIR}")
    if not DOMAIN_FILE.exists():
        raise FileNotFoundError(f"domain file not found: {DOMAIN_FILE}")

    GROUND_DIR.mkdir(parents=True, exist_ok=True)

    bounds, n_dom, crs = domain_crop_bounds(DOMAIN_FILE, CROP_BUFFER_M)
    print(f"\nCrop to {n_dom} domains + {CROP_BUFFER_M:g} m ({crs}):")
    print(f"  x {bounds[0]:.1f} .. {bounds[1]:.1f}")
    print(f"  y {bounds[2]:.1f} .. {bounds[3]:.1f}")
    print(f"\nSMRF: cell={SMRF_CELL} slope={SMRF_SLOPE} window={SMRF_WINDOW} "
          f"threshold={SMRF_THRESHOLD} scalar={SMRF_SCALAR}")

    tiles = sorted(POINT_DIR.glob(TILE_GLOB))
    if not tiles:
        raise FileNotFoundError(f"no {TILE_GLOB} under {POINT_DIR}")
    print(f"\n{len(tiles)} tile(s) to classify")

    print(f"Splitter: {SPLIT_LENGTH_M:g} m cells, {SPLIT_BUFFER_M:g} m buffer")
    print(f"Outlier filter: {'on' if USE_OUTLIER_FILTER else 'OFF'}")
    CELL_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    for tile in tiles:
        tag = tile.stem
        print(f"\n{'=' * 70}\n{tile.name}", flush=True)
        n_in, z_in = count_points(tile)
        print(f"  in : {n_in:,} pts   z {z_in[0]:.2f} .. {z_in[1]:.2f}", flush=True)

        # ---- PHASE 1: crop + split ----
        cells = sorted(CELL_DIR.glob(f"{tag}_cell_*.laz"))
        if cells:
            print(f"  [1/2] SKIP split - {len(cells)} cell file(s) present")
        else:
            print(f"  [1/2] splitting into {SPLIT_LENGTH_M:g} m cells...", flush=True)
            run_pipeline(pdal_exe,
                         split_pipeline(tile, CELL_DIR / f"{tag}_cell_#.laz", bounds),
                         f"{tag}_split")
            cells = sorted(CELL_DIR.glob(f"{tag}_cell_*.laz"))
            print(f"        {len(cells)} cell(s) written")

        if not cells:
            print("  WARNING: crop does not intersect this tile - skipping.")
            results.append((tile.name, n_in, 0, 0, (float('nan'),) * 2))
            continue

        # ---- PHASE 2: classify each cell in its own process ----
        print(f"  [2/2] classifying {len(cells)} cell(s)", flush=True)
        for i, cell in enumerate(cells, 1):
            idx = cell.stem.rsplit("_", 1)[-1]
            out_cell = GROUND_DIR / f"{tag}_ground_{idx}.laz"
            if out_cell.exists():
                print(f"        [{i}/{len(cells)}] {cell.name} SKIP", flush=True)
                continue
            n_cell, _ = count_points(cell)
            try:
                run_pipeline(pdal_exe, ground_pipeline(cell, out_cell),
                             f"{tag}_ground_{idx}", quiet=True)
            except RuntimeError as e:
                # One bad cell must not lose the other 40
                print(f"        [{i}/{len(cells)}] {cell.name} FAILED: {e}",
                      flush=True)
                continue
            n_g = count_points(out_cell)[0] if out_cell.exists() else 0
            print(f"        [{i}/{len(cells)}] {cell.name}  {n_cell:>9,} -> "
                  f"{n_g:>9,} ground ({100 * n_g / max(n_cell, 1):4.1f}%)",
                  flush=True)

        existing = sorted(GROUND_DIR.glob(f"{tag}_ground_*.laz"))
        if existing:
            n_out, z_out = count_many(existing)
            pct = 100 * n_out / n_in if n_in else 0
            print(f"  out: {n_out:,} ground pts in {len(existing)} cell(s) "
                  f"({pct:.1f}% of tile, after crop)"
                  f"   z {z_out[0]:.2f} .. {z_out[1]:.2f}")
            results.append((tile.name, n_in, n_out, len(existing), z_out))
        else:
            print("  WARNING: no ground output written.")
            results.append((tile.name, n_in, 0, 0, (float('nan'),) * 2))

    print("\n" + "=" * 70)
    print("SUMMARY")
    for name, n_in, n_out, n_cells, z_out in results:
        print(f"  {name:55s} {n_out:>12,} ground in {n_cells:>3d} cell(s)")
    print(f"\nGround-only cells: {GROUND_DIR}")

    print("""
QC BEFORE YOU TRUST THIS
------------------------
1. DUNE CRESTS. Compare the classified max z against the 2009 DEM along a few
   domains. If ground z is systematically below the 2009 dune crest by more
   than the survey difference, SMRF_WINDOW is too small and is shaving crests -
   raise it and re-run.
2. BUILDINGS. A wide window keeps some structures as "ground". They sit on the
   developed side, away from the marsh fringe this fill targets, so they are
   mostly harmless - but check none land inside a domain's fill area.
3. WATER. elm + outlier should drop most water returns, but the 2008 waterline
   is where the fill stops anyway, so residual water points at the fringe are
   the main thing that could bias the filled marsh low.
Then run HAT_dem_gap_fill.py.""")


if __name__ == "__main__":
    main()
