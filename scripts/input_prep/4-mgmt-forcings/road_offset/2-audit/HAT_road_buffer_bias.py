"""
How far the setback's buffered-mask edge sits from the real road, and whether it matters.

    python scripts/input_prep/4-mgmt-forcings/road_offset/2-audit/HAT_road_buffer_bias.py

Mask seaward edge against the geojson centreline, per profile, in the
original raster frame; prints the bias. Details: scripts/input_prep/4-mgmt-forcings/road_offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-08-20
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import geopandas as gpd

# The centerline-to-raster-column projection is already solved next door, half-cell convention and all
sys.path.insert(0, str(Path(__file__).resolve().parent))
from HAT_check_geojson_vs_mask import (  # noqa: E402
    geojson_cols_by_row, MASK_FMT, GEOJSON_FMT, TIF_FMT, DOMAINS, CELL_SIZE_M,
)

# --- CONFIG ------------------------------------------------------------------
YEARS = [1984, 2004]

# Bulldoze's modelled road block, from cascade_pipeline/roadway.py (road_width_m) and roadway_manager.bulldoze
ROAD_WIDTH_M = 20.0

# Domains excluded from the model-facing span, reported separately rather than dropped
EXCLUDED = [8]
# -----------------------------------------------------------------------------


# Per-profile bias in metres, and the per-domain median of it
def measure_year(year: int) -> tuple[np.ndarray, dict]:
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
            # Ocean is at the RIGHT in the raw frame (OCEAN_LOC = "right")
            seaward_raw = int(cells.max())
            # Cell k spans [k, k+1), so a point sits at c - 0.5 in cell-centre units
            c = line_cols[r]
            vals.append(((c - 0.5) - seaward_raw) * CELL_SIZE_M)

        if vals:
            per_domain[d] = float(np.median(vals))
            per_profile.extend(vals)

    return np.asarray(per_profile), per_domain


# Run: both vintages, the bias per domain and overall
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
