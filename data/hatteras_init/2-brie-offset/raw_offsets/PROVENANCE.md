# raw_offsets — which survey is behind each file

One file per dune-line survey, holding the per-transect distance from the
shared offshore baseline (`ORIG_LEN`). Everything in `hindcast_<year>/` is
built from these, and the end-year validation target differences two of them,
so the baseline must be the same layer in every file here.

## Files are named for the LINE'S VINTAGE, and a table pairs periods with them

`<vintage>_duneline_offset_raw.csv` is the current build of the line traced
from that year's imagery. A hindcast period finds its file through
`hat_topo_version.DUNE_LINE_FOR_YEAR` (start and end years alike; the
end-year target loader and `island_offset_hybrid.py` both go through it), the
same rule the NC-12 lines follow with `ROAD_LINE_FOR_YEAR`. Since 2026-09-15
(Hannah: "a table, no copies"); from 2026-09-11 to then the 1996 start read a
byte-identical copy of the 1997 file filed under the name `1996_...`, which
put two files on one survey.

| period year | line vintage | note |
|---|---|---|
| 1984 | 1984 | a 1984 survey, no stand-in |
| 1996 | **1997** | no 1996 line exists; the nearest island-wide survey |
| 2004 | 2004 | a 2004 survey, no stand-in |
| 2010 | **2009** | no 2010 aerial imagery; the 2009 line (Hannah, 2026-09-15) |
| 2024 | **2023** | the NOAA 2023 imagery; end year only, no padded build is made for it |

Consequences of the 1996 stand-in, stated rather than left to be discovered:

* the 1996 period starts from the island as surveyed a year later;
* a 1996-to-2010 shoreline change computed from these files spans thirteen
  years of survey while the run spans fourteen.

The 2010 and 2024 lines arrived on 2026-09-15 and went into the table under
the year of their imagery (2009-05-30 Google Earth; the NOAA 2023 set), which
is why the 2010 start reads `2009_...` and the 2024 end reads `2023_...`; the
raw file and the geojson carry the imagery vintage in their names. The same
stand-in caveat as 1996 applies: the 2010 period starts from the island as
surveyed a year earlier and ends a year early too, so a 2010-to-2024 change
from these files spans fourteen years of survey (2009 to 2023), and a
2004-to-2024 change spans nineteen (2004 to 2023) while the run spans twenty.

One vintage can have more than one digitisation (`duneline_1997.geojson`,
`duneline_1997_v2.geojson`). The raw file under the vintage's name is the
current one, and every build under `<year>/v<n>/` keeps a copy of the raw it
was made from, so an older build can always be reproduced with
`island_offset_hybrid.py --raw-file`.

## Which file was made how

| File | Line | Made by |
|---|---|---|
| `superseded_20260915_gis_exports/1997_v1_duneline_offset_raw.csv` | `duneline_1997.geojson` (v1, 2026-09-02) | ArcGIS: 1.5 m buffer, 1 m points along the transects, export. ~3 rows per transect. Was `1997_duneline_offset_raw.csv` until the v2 build took the vintage's name |
| `1997_duneline_offset_raw.csv` | `duneline_1997_v2.geojson` (local corrections, 2026-09-15) | `scripts/input_prep/2-brie-offset/duneline_to_raw_offsets.py`: shapely intersection with `../transects/transects_100m.geojson`, one row per transect. Was `1997_v2_...` for a few hours on 2026-09-15 |
| `1984_duneline_offset_raw.csv` | `duneline_1984.geojson` | the script, 2026-09-15; an ArcGIS export until then |
| `2004_duneline_offset_raw.csv` | `duneline_2004.geojson` | the script, 2026-09-15; an ArcGIS export until then |
| `2009_duneline_offset_raw.csv` | `duneline_2009.geojson` | the script, 2026-09-15 (by `build_island_offset.py --year 2010`) |
| `2023_duneline_offset_raw.csv` | `duneline_2023.geojson` | the script, 2026-09-15 (built as `2024_...` for an hour, before the line was renamed for its imagery) |

The shapely script was validated on the v1 line against the v1 export: all
450 transects match, mean −1.01 m, sd 0.32 m. **The metre is a convention,
not an error**: the export lists the landward-most 1 m station inside its
buffer first, and both `island_offset_hybrid.py` and
`hindcast.load_absolute_dune_distance` keep the first row per transect, so a
GIS file sits ~1 m landward of the exact crossing. It cancels within a year
(each build is zeroed on its own minimum) and between two shapely files. It
does NOT cancel between a GIS file and a shapely file, which is why the
**1984 and 2004 files were rebuilt by the script on 2026-09-15**: the ArcGIS
exports they replaced are under `superseded_20260915_gis_exports/` with the
per-transect comparison (1984: shapely −1.28 m from the export, sd 0.32;
2004: −1.01 m, sd 0.28; every transect crossed exactly once in both). So
1984, 1996 (= 1997 v2), 2004 and every line from here are one method, and a
change between any two of them carries no method term. The **2010 and 2024
lines go through the script** like the rest.

Every raw file from here is made by the script, not exported, so the column
set is the script's (`domain_id, LineID, ORIG_LEN, n_crossings, x, y` plus
the line's metadata). The three columns the loaders read are the same in
both kinds of file. No export is left at this level; the v1 1997 one, the
script's validation reference, is under `superseded_20260915_gis_exports/`
with the 1984 and 2004 ones.

## Building from a line

`scripts/input_prep/2-brie-offset/build_island_offset.py --duneline <geojson>
--year <start>` runs the intersection, writes the vintage's raw file here,
builds `<year>/v<n>/`, compares it with the previous build, sets `CURRENT`
and writes the provenance. It refuses a geojson whose vintage the table does
not pair with `--year`.

## Coverage

`2017_duneline_offset_raw.csv` is a **Buxton-only clip**, eleven domains, not
an island-wide survey. It cannot stand in for a hindcast start or an end-year
target. `1978` and `1997` carry domains past GIS 90; everything downstream
keeps 1 to 90 and drops the rest. The shapely-built files only ever hold 1 to 90.
