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
| `duneline_1997.geojson` | EPSG:3725 | 529 | full* | **re-digitized 2026-09-18** from the USGS photo of 1997-10-12; behind `2-brie-offset/1996/v1/` (CURRENT; built as v3, renumbered 2026-09-19). Also what `HAT_measure_duneline_shift.py` reads |
| `duneline_2004.geojson` | EPSG:3725 | 158 | partial | the 2004 initial condition; Google Earth capture of 2004-05-25 |
| `duneline_2009.geojson` | EPSG:26918 | 603 | none | **re-digitized 2026-09-18** (behind `2-brie-offset/2010/v1/`, CURRENT; built as v2, renumbered 2026-09-19); Google Earth capture of 2009-05-30; stands in for the **2010** start (no 2010 aerial imagery), `DUNE_LINE_FOR_YEAR[2010] == 2009` |
| `duneline_2023.geojson` | EPSG:26918 | 749 | none | **re-digitized 2026-09-18**; first added 2026-09-15 as `duneline_2024` and renamed for its imagery, the NOAA NGS 2023 set under `D:\Hatteras_GIS\Aerial3`; stands in for the **2024** end year, `DUNE_LINE_FOR_YEAR[2024] == 2023` |

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

## The 2026-09-18 re-digitization

Hannah re-digitized the 1997, 2009 and 2023 lines on 2026-09-18, and they
**replace** the earlier ones under the plain vintage name. That was her
decision: the earlier 1997, 1997_v2, 2009 and 2023 files survive only in git
history from before that date. Each old build folder under
`2-brie-offset/<year>/v<n>/` keeps the raw offsets it was made from, so it
can still be rebuilt.

\* The 1997 file's `edit_date` still reads 2025-11-04. The properties were
carried over from the earlier digitization, so that date does not describe
the 2026-09-18 geometry.

Measured per domain along the 100 m transects, every change went the same
way: the new line lies SEAWARD of the old.

| line | domains moved ≥ 0.5 m | largest | build |
|---|---|---|---|
| 1997 (vs 1997_v2) | 32 of 90 | 66.2 m, GIS 35 | `1996/v1/` (was v3), `offset_1996_v2_vs_v3.csv` |
| 2009 | 50 of 90 | 65.6 m, GIS 2 | `2010/v1/` (was v2), `offset_2010_v1_vs_v2.csv` |
| 2023 | 53 of 90 | 44.6 m, GIS 32 | no build (only an end year); `raw_offsets/2023_...` |

Rebuilt from the new lines the same day:
- the three raw offset files and the two Pea Island extension raws
- `1996/v3`, `2010/v2` and `1996/ext/n115` (since 2026-09-19 `1996/v1` and `2010/v1`: the numbering restarts at the re-digitized lines, the earlier builds are under `<start>/superseded_20260919_pre-redigitized/`)
- `5-scr/4-comparisons/shoreline_vs_duneline/net_change/` (was `duneline_vs_coastsat/`), for the three windows that
  use these lines

NOT rebuilt:
- `5-scr/3-rates/duneline_lrr/`, a fit through several lines. Retired later
  that day for `3-rates/duneline/endpoint/` (net change, built from these
  lines); the fit is in `5-scr/archive/duneline_lrr_retired_20260918/`
- the 1984 dune-line shift measurement
- every model run on the 1996 or 2010 inputs

## 1997 v2 (superseded 2026-09-18, kept as history)

Local corrections only (the rest of the line is the same vertices). Where they
land, measured along the 100 m transects: 22 of 90 domains moved, all of them
landward, at GIS 16, 46-51, 62-67, 79-81 and 86, by up to 63 m (GIS 66); see
`2-brie-offset/1996/v2/PROVENANCE.md`. The 1996 hindcast start reads the v2
build. The dune-line shift measurement (`../1984-start/duneline-shift/`) was
made on v1 and has not been re-run; its +16.2 m definitional term was derived
on v1 and is quoted above as such.

