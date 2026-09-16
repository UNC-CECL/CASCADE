# superseded_20260915_gis_exports

The ArcGIS exports that were `raw_offsets/1984_duneline_offset_raw.csv` and
`2004_duneline_offset_raw.csv` until 2026-09-15, when both lines were re-run
through `scripts/input_prep/2-brie-offset/duneline_to_raw_offsets.py` so
that every raw file is built the same way (Hannah: "rerun the 1984 and 2004
lines through the shapely script now").

Also here since the same afternoon: `1997_v1_duneline_offset_raw.csv`, the
ArcGIS export of the v1 1997 line, which was `raw_offsets/1997_duneline_offset_raw.csv`
until the shapely build of the v2 line took that name (raw files are named
for the vintage, one current build each; see `../PROVENANCE.md`). It is the
file `duneline_to_raw_offsets.py --validate-against` was validated on.

Kept for the record, not for use. Nothing reads this folder.

| line | transects | shapely minus export, per transect | per domain mean |
|---|---|---|---|
| 1984 | 450 of 450, every one crossed once | mean −1.28 m, sd 0.32, range −2.25 .. −0.56 | −1.28, sd 0.18 |
| 2004 | 450 of 450, every one crossed once | mean −1.01 m, sd 0.28, range −1.58 .. −0.50 | −1.01, sd 0.13 |

The constant metre is the export convention (the landward-most 1 m station
inside a 1.5 m buffer, listed first), not geometry; see `../PROVENANCE.md`.
The 1984 export sits a quarter-metre further landward than the 2004 and 1997
ones did, so the old 1984→2004 change carried +0.27 m of method. The
per-transect comparisons are `1984_shapely_vs_gis_export.csv` and
`2004_shapely_vs_gis_export.csv`.
