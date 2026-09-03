# raw-duneline-geojson

Digitized NC-12-era dune lines, in the repo rather than on `D:\Hatteras_GIS`.

| file | CRS | vertices | metadata | notes |
|---|---|---|---|---|
| `duneline_1967.geojson` | EPSG:26918 | 96 | none | oldest, unused so far |
| `duneline_1984.geojson` | EPSG:26918 | 495 | none | the 1984 initial condition |
| `duneline_1997.geojson` | EPSG:3725 | 581 | full | digitized 2026-09-02 |
| `duneline_2004.geojson` | EPSG:3725 | 158 | partial | the 2004 initial condition |

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
