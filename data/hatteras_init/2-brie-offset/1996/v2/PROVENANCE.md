# 1996 island offsets, v2

**Derived, not surveyed.** Built from `duneline_1997_v2.geojson`, the 1997
dune line with **local corrections** (Hannah, 2026-09-15; the rest of the line
is the same vertices as v1). No 1996 line exists, and 1997 is the nearest
island-wide survey. What follows from that one-year gap is in `../PROVENANCE.md`.

## How it was built

1. `scripts/input_prep/2-brie-offset/duneline_to_raw_offsets.py --duneline
   duneline_1997_v2.geojson --out 1997_v2_duneline_offset_raw.csv`
   intersects the line with the 100 m transects (`../../transects/`) in
   shapely and writes the station of each crossing from the offshore datum,
   one row per transect, 450 transects, five per domain, no transect missed
   and none crossed twice. This is the step ArcGIS did for every earlier raw
   file.
2. That file was copied to `raw_offsets/1996_duneline_offset_raw.csv`, the
   name the end-year target loader resolved for the 1996 start at the time.
   Later the same day the copy was removed: the file is
   `raw_offsets/1997_duneline_offset_raw.csv` (named for its vintage) and the
   1996 start finds it through `hat_topo_version.DUNE_LINE_FOR_YEAR`. A copy
   of it sits in this folder as the raw this build came from.
3. `island_offset_hybrid.py --year 1996 --version v2` averaged the five
   transects per domain, zeroed on the minimum (GIS 90 in both versions), and
   padded to 120.

## The intersection was validated before it was trusted

Run on the **v1** line and compared with the ArcGIS export of the same line
(`1997_duneline_offset_raw.csv`): all 450 transects match, mean −1.01 m, sd
0.32 m, worst 1.9 m. The constant metre is the export's station convention
(the landward-most 1 m point inside the 1.5 m buffer is listed first, and the
loaders take the first row), not geometry. It cancels in the per-year
zeroing here, and in any end-year difference of two files built by the same
script. **It does not cancel between a GIS-built and a shapely-built file**,
so a 1996→2010 target needs the 2010 line put through this same script, and a
1984→1996 difference (both raw files exist) carries the metre — a tenth of a
cell, below the digitising noise, but recorded.

## What changed from v1

`offset_1996_v1_vs_v2.csv` / `.png` / `.pdf` (caption in `CAPTIONS.md`),
from `HAT_compare_offset_versions.py`. The v1 side of the raw comparison is
the v1 line re-intersected by the same shapely script, so the metre above is
not in these numbers.

* **22 of 90 domains** differ by 0.5 m or more, and **every one moved
  landward** (a larger station from the datum): mean over all 90 domains
  +6.3 m.
* Where: GIS 16 (+20 m), 46–51 (+47, +4, +8, +17, +36, +13 m), 62–67
  (+31, +50, +56, +44, +63, +45 m), 79–81 (+53, +43, +10 m), 86 (+15 m), and
  under 5 m at 17, 30, 31, 35, 61.
* The per-year minimum moved +0.06 m, so the model-frame difference is the
  same picture: mean +6.4 m, range −0.2 to +63.4 m.
* Against the 1984 raw file, the island-wide 1984→1996 landward change goes
  from +3.5 m (v1) to +8.8 m (v2).

Sanity check that still holds: the profile sits between the 1984 and 2004
profiles.
