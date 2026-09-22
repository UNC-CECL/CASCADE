# shoreline_inventory — edits to the source data

## 2026-09-22: the two Hatteras 1997 features restamped to 9/27/1997

`nc_shorelines.geojson` stamped every one of its 23 `SHR_YEAR = 1997`
features `1/1/1997`. That is a placeholder, not a survey date (Hannah,
2026-09-22); the Hatteras 1997 shoreline was flown **1997-09-27**, the date
also commented into `coastsat_domain_lrr_specific_dates.py`.

**Two features were changed**, both `SHR_DATE` and `DATE_`:

| LOCATION | was | now |
|---|---|---|
| Outer Banks - National Seashore | 1/1/1997 | 9/27/1997 |
| Outer Banks - North of Oregon Inlet | 1/1/1997 | 9/27/1997 |

**The other 21 were left alone.** They are elsewhere in North Carolina —
Carolina Beach, Masonboro, Topsail, Bogue Banks, Onslow Beach, Bald Head
Island and the rest — all `ANALYST = USGS`, all carrying the same `1/1/1997`
stamp. Those are different surveys on different days; 9/27/1997 is a Hatteras
flight date and asserting it statewide would replace one wrong date with
another. Which two features to change was decided spatially, by intersecting
each 1997 feature with `cascade_area.geojson`.

Note the file is **EPSG:3857** and `cascade_area.geojson` declares `EPSG:3725`,
which is not a real EPSG code — it is the ArcGIS export artefact
`coastsat_domain_mapping.fix_crs` already knows about, and was read as
EPSG:32618 (UTM 18N) for the intersection.

The file as it stood before the edit is in
`superseded_20260922_pre_restamp/`. These geojsons are gitignored, so that
copy is the only history they have.
