# transects — the 100 m transect layer the offsets are measured along

`transects_100m.geojson` (EPSG:3725, 622 LineStrings, 0.8 MB, in LFS) is the
transect set every `raw_offsets/` file is built on. Each transect is 10 km
long, starts on the offshore datum line (x = 460198 m) and runs west across
the island; `ORIG_LEN` in a raw file is the station along one of these from
that start. `domain_id` is the ArcGIS spatial join onto the 500 m domain
polygons: 450 transects fall in GIS 1-90, exactly five per domain, LineID 12
through 463. The other 172 carry no domain and are ignored.

Copied here 2026-09-15 from
`hard-structures/groin/HAT-groin-gis-analysis/gis_data/transects_100m.geojson`
(byte-identical) so that `2-brie-offset/` holds everything a dune line needs
to become an offset file, and `duneline_to_raw_offsets.py` does not reach into
the groin study for an input. The ArcGIS exports in `raw_offsets/` were made
on this same layer (their `LineID`s and per-domain transect sets match it
exactly).

The column names are ArcGIS join names (`Transects_100m.LineID`,
`Transects_100m_AddSpatialJoin.domain_id`, ...). The script strips them to the
leaf name.
