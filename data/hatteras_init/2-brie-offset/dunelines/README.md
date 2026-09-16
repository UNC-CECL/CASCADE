# dunelines

The digitised dune lines, in the repo rather than on `D:\Hatteras_GIS`. Here
since 2026-09-15 (Hannah: they are the island-offset input, so they belong
beside `../transects/` and `../raw_offsets/`, not under the topography tree
where they sat as `1-barrier3d-domains/raw-duneline-geojson/`).

**Naming:** `duneline_<vintage>[_v<n>].geojson`, one feature per file.
`<vintage>` is the year of the IMAGERY the line was traced from, never the
hindcast period it serves; a period finds its line through
`hat_topo_version.DUNE_LINE_FOR_YEAR`. `_v<n>` is a re-digitisation of the
same vintage. `hat_topo_version.DUNELINE_DIR` and `duneline_geojson(vintage,
version)` are how scripts spell the location.

**Properties every new line should carry:** `feature_type`, `year`,
`imagery_date`, `source_type`, `method`, `editor`, `edit_date`, `notes`.
`build_island_offset.py` copies them into the build's provenance and warns
when `imagery_date` is missing (the 1984 and 2004 dates had to be recovered
from a USGS metadata file and from memory).

| file | CRS | vertices | metadata | notes |
|---|---|---|---|---|
| `duneline_1967.geojson` | EPSG:26918 | 96 | none | oldest, unused so far |
| `duneline_1984.geojson` | EPSG:26918 | 495 | none | the 1984 initial condition; USGS photo of 1984-09-19 |
| `duneline_1997.geojson` | EPSG:3725 | 581 | full | digitized 2026-09-02 from the USGS photo of 1997-10-12; still what `HAT_measure_duneline_shift.py` reads |
| `duneline_1997_v2.geojson` | EPSG:3725 | 592 | full | 2026-09-15, local corrections to the line above; behind `2-brie-offset/1996/v2/` |
| `duneline_2004.geojson` | EPSG:3725 | 158 | partial | the 2004 initial condition; Google Earth capture of 2004-05-25 |
| `duneline_2009.geojson` | EPSG:26918 | 740 | none | added 2026-09-15; Google Earth capture of 2009-05-30; stands in for the **2010** start (no 2010 aerial imagery), `DUNE_LINE_FOR_YEAR[2010] == 2009` |
| `duneline_2023.geojson` | EPSG:26918 | 868 | none | added 2026-09-15 as `duneline_2024` and renamed for its imagery, the NOAA NGS 2023 set under `D:\Hatteras_GIS\Aerial3`; stands in for the **2024** end year, `DUNE_LINE_FOR_YEAR[2024] == 2023` |

`HAT_measure_duneline_shift.py` reads THIS directory first and falls back to
`D:\Hatteras_GIS\Dunelines` only for older invocations. The external drive is
not version-controlled, not present on another machine, and not something a run
can record the state of, so the curated copy is the source.

The loader reprojects from whatever CRS each file declares, so the mixed CRSs
above are handled rather than assumed away.

## Why 1997 exists

Measuring the 1984 line against the model's interior row 0 confounds two things:
how far the island moved (date) and the offset between a digitized line and the
extractor's row 0 (definition). Differencing 1984 against another digitized line
cancels the definitional term exactly, because the same feature is at both ends.

The split it produced, island-wide: total +18.9 m = feature **+16.2 m** + date
**+0.8 m**. So roughly 85% of the naive number was definitional. See
`../1984-start/duneline-shift/README.md`.

## The comparability caveat

`duneline_1997` carries `feature_type`, `method`, `editor` and `edit_date`; its
method is *"Digitized from light/dark elevation break (no DEM available)"* and
its `feature_type` reads *"Island orientation"*. **The 1984 and 1967 files carry
no metadata at all**, so "the same feature at both ends" is supported by the
numbers — the 1997-vs-row-0 offset is a near-constant +16.2 m, IQR +12.8 to
+21.0, which is what a definitional offset looks like and not what noise looks
like — but it is not documented in the files themselves.

Note also that 1997 is one year after the 1996 ALACE survey the DEM's beach comes
from. At the fastest measured rate (GIS 85, 4.4 m/yr) that is ~4 m, under half a
cell.

## 1997 v2

Local corrections only (the rest of the line is the same vertices). Where they
land, measured along the 100 m transects: 22 of 90 domains moved, all of them
landward, at GIS 16, 46-51, 62-67, 79-81 and 86, by up to 63 m (GIS 66); see
`2-brie-offset/1996/v2/PROVENANCE.md`. The 1996 hindcast start reads the v2
build. The dune-line shift measurement (`../1984-start/duneline-shift/`) was
made on v1 and has not been re-run; its +16.2 m definitional term was derived
on v1 and is quoted above as such.

