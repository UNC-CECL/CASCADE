# ==============================================================================
# HAT_road_buffer_bias.py
#
# The setback is measured to the SEAWARD EDGE OF A BUFFERED MASK, not to the road.
# How far apart are those, and does it matter?
#
# WHY THIS EXISTS
#   HAT_rasterize_road_to_domains.py burns NC-12 with ROAD_BUFFER_M = 6.0 and
#   ALL_TOUCHED = True, giving a ~24 m mask for an ~8 m road.
#   HAT_road_offset_from_dune_start.py then takes that mask's seaward-most cell as
#   `road_start`, because `bulldoze` indexes [road_start : road_start + width].
#   So every setback inherits however far the buffer pushed that edge seaward, and
#   nothing else in this tree measures it -- the geojson is used to rasterize, to
#   check registration (HAT_check_geojson_vs_mask.py) and to draw maps, never to
#   produce a setback.
#
# WHAT IT MEASURES
#   bias = (mask seaward cell) - (geojson centerline position), in metres,
#   NEGATIVE meaning the mask edge sits SEAWARD of the true centerline and the
#   reported setback is therefore SMALLER than a centerline measurement.
#
#   Measured in the ORIGINAL raster frame. `setback = (seaward - row0) * cell` and
#   row0 is identical either way, so orient / flip / shear / trim all cancel out of
#   the difference -- no transform is needed, and the result is independent of the
#   extractor.
#
# WHAT IT FOUND (2026-08-19, both vintages)
#   median bias -6.8 m, p10 -10.8, p90 -2.7.
#
#   That is NOT a 6.8 m placement error, because of what it is compared against.
#   `road_width_m = 20.0` and bulldoze lays a 2-cell block LANDWARD from
#   road_start, so a block centred on the real road wants
#   road_start = centerline - 10 m. The buffer supplies centerline - 6.8 m. The
#   residual misplacement is 3.2 m -- about a third of a cell, and an order of
#   magnitude below the 20-40 m p90 placement error that the alongshore collapse
#   to one scalar already causes (HAT_method_comparison_figures.py).
#
#   DO NOT "CORRECT" IT. Per-profile setbacks are integer cell differences times
#   10 m, so domain medians land almost exactly on cell boundaries, and bulldoze
#   truncates with int(). Applying the 3.2 m adjustment moves road_start a FULL
#   cell (10 m) seaward in 83% of 1984 domains and 90% of 2004 domains, replacing
#   a 3 m error with a 7 m error in the other direction. The buffer is accidentally
#   doing very nearly the right job; the quantization is what would bite.
#
#   Cross-validated independently: the road_x/road_y columns
#   HAT_road_offset_from_dune_start.py writes into RoadOffset_<year>_profiles.csv
#   are produced by a completely different path (full index inversion to map
#   coordinates), and their shapely distance to the same geojson is 6.62-6.68 m
#   median -- the same number to ~0.2 m.
#
#   D8 is the one extreme outlier (-353 m in 1984). That is the Buxton bend, where
#   NC-12 runs parallel to the raster rows, so a "seaward-most cell" is not a
#   meaningful quantity -- the same reason D8 is already EXCLUDED_FROM_SPAN. Its
#   appearance here is a consistency check passing, not a new problem.
#
# USAGE
#     python HAT_road_buffer_bias.py
#
# REQUIREMENTS
#     numpy, geopandas, rasterio   (same set HAT_check_geojson_vs_mask.py needs)
# ==============================================================================

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import geopandas as gpd

# The centerline-to-raster-column projection is already solved next door, half-cell
# convention and all. Importing it rather than re-deriving it keeps ONE definition
# of "where is the centerline in this raster".
sys.path.insert(0, str(Path(__file__).resolve().parent))
from HAT_check_geojson_vs_mask import (  # noqa: E402
    geojson_cols_by_row, MASK_FMT, GEOJSON_FMT, TIF_FMT, DOMAINS, CELL_SIZE_M,
)

YEARS = [1984, 2004]

# bulldoze's modelled road block, from cascade_pipeline/roadway.py (road_width_m)
# and roadway_manager.bulldoze: road_end = road_start + int(road_width / dx).
ROAD_WIDTH_M = 20.0

# Domains excluded from the model-facing span, reported separately rather than
# dropped. Must match ROAD_SPAN in HAT_road_offset_from_dune_start.py.
EXCLUDED = [8]


def measure_year(year: int) -> tuple[np.ndarray, dict]:
    """Per-profile bias in metres, and the per-domain median of it."""
    geo = gpd.read_file(str(GEOJSON_FMT).format(year=year))
    per_profile: list[float] = []
    per_domain: dict[int, float] = {}

    for d in DOMAINS:
        tif = Path(str(TIF_FMT).format(d=d))
        mask_path = Path(str(MASK_FMT).format(year=year, d=d))
        if not tif.is_file() or not mask_path.is_file():
            continue

        mask = np.load(mask_path)
        mask = np.isfinite(mask) & (mask > 0)
        line_cols = geojson_cols_by_row(geo, tif)

        vals = []
        for r in range(mask.shape[0]):
            cells = np.flatnonzero(mask[r])
            if not cells.size or r not in line_cols:
                continue
            # Ocean is at the RIGHT in the raw frame (OCEAN_LOC = "right"), so
            # ocean-first index 0 is the LAST raw column and the SEAWARD-most mask
            # cell is the MAXIMUM raw column.
            seaward_raw = int(cells.max())
            # `c` is a fractional column from the inverse affine, where cell k
            # spans [k, k+1); a point therefore sits at cell-centre units c - 0.5.
            # Comparing against the raw integer injects half a cell (5 m).
            c = line_cols[r]
            vals.append(((c - 0.5) - seaward_raw) * CELL_SIZE_M)

        if vals:
            per_domain[d] = float(np.median(vals))
            per_profile.extend(vals)

    return np.asarray(per_profile), per_domain


def main() -> None:
    print("=" * 78)
    print("road buffer bias -- mask seaward edge vs geojson centerline")
    print("ORIGINAL raster frame; orient/flip/shear/trim cancel out of the diff")
    print("=" * 78)

    for year in YEARS:
        b, per_dom = measure_year(year)
        if not b.size:
            print(f"\n--- {year}: no data")
            continue

        keep = {d: v for d, v in per_dom.items() if d not in EXCLUDED}
        pd_ = np.asarray(list(keep.values()))

        print(f"\n--- {year}   {b.size} profiles, {len(per_dom)} domains")
        print(f"  per-profile bias  : median {np.median(b):+.1f} m   "
              f"p10 {np.percentile(b, 10):+.1f}   p90 {np.percentile(b, 90):+.1f}")
        print(f"  per-domain median : median {np.median(pd_):+.1f} m   "
              f"min {pd_.min():+.1f}   max {pd_.max():+.1f}   "
              f"(excluding {EXCLUDED})")
        print("  negative = mask edge SEAWARD of centerline, so the reported")
        print("             setback is SMALLER than a centerline measurement")

        # What the bias means once bulldoze's own geometry is accounted for.
        half = ROAD_WIDTH_M / 2.0
        residual = half + np.median(pd_)     # median bias is negative
        print(f"\n  bulldoze lays a {ROAD_WIDTH_M:.0f} m block LANDWARD from "
              f"road_start, so a block")
        print(f"  centred on the real road wants road_start = centerline "
              f"- {half:.0f} m.")
        print(f"  The buffer supplies centerline {np.median(pd_):+.1f} m, leaving a "
              f"residual of {residual:.1f} m")
        print(f"  ({residual / CELL_SIZE_M:.2f} of a cell) -- see the header for why "
              f"correcting it is worse.")

        for d in EXCLUDED:
            if d in per_dom:
                print(f"\n  [excluded] D{d}: {per_dom[d]:+.0f} m -- the Buxton bend, "
                      f"road parallel to the raster rows.")
                print(f"             Expected; D{d} is EXCLUDED_FROM_SPAN for the "
                      f"same reason.")

    print()


if __name__ == "__main__":
    main()
