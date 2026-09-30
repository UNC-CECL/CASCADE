# 2-transect-frame — which transect belongs to which domain

One script, one job, and everything downstream depends on it. A CoastSat
transect is a line in space; a Barrier3D domain is a 500 m polygon. This builds
the lookup that joins them, and every rate fit in `../3-rates/` reads it.

```
coastsat_domain_mapping.py
    Spatially joins the CoastSat transect geometry to the 90 domain polygons:

        transect_id | domain_number | distance_to_domain_m | match_method

    Writes 2-transect-frame/transect_domains/transect_domain_lookup.csv.
```

**Nearest-feature snapping is off, deliberately.** Only point-in-polygon
matches are kept. A transect that falls outside every domain is left unmatched
rather than attached to whichever domain happens to be closest — a wrong
assignment is invisible downstream, where it silently pulls one domain's mean
rate toward its neighbour's.

This runs first and rarely: only when the transect layer or the domain
polygons change. It is not part of rebuilding a rate window.

The extension experiment numbers transects beyond GIS 90 by a different route
(`hat_extension_domains.join_origins`, used by
`../3-rates/coastsat/extension/coastsat_extension_lrr.py`), because the 90
polygons this lookup was built on simply stop there.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### coastsat_domain_mapping.py

Step 1 of 2: join CoastSat transects to the CASCADE domains, the lookup every rate script uses.

From the script's original header:

```text
CoastSat Transect → CASCADE Domain Spatial Mapping
Step 1 of 2 in the domain-level LRR workflow.

This script spatially joins CoastSat transect geometry to your CASCADE
domain geometry, producing a lookup table:

    transect_id  |  domain_number  |  distance_to_domain_m  |  match_method

That lookup table is then used by coastsat_domain_lrr.py to compute
domain-level LRR summaries.

Note: nearest-feature snapping is disabled. Only point-in-polygon matches
are kept. Transects outside all domain polygons are excluded.

Inputs
  - CoastSat transects GeoJSON (downloaded from coastsat.space)
  - CASCADE domain polygons GeoJSON or shapefile (exported from ArcGIS)

Outputs
  transect_domain_lookup.csv  –  saved to OUTPUT_DIR
  transect_domain_map.png     –  quick-look map to verify the join

Dependencies
  pip install geopandas pandas numpy matplotlib shapely pyproj
```

Notes that were in the code:

```text
--- CoastSat transect geometry ---
GeoJSON downloaded from coastsat.space → "transects" button
the script will automatically
clip it to your study area using the domain bounding box.
Anchored on this file 2026-09-12. The literals here were
drive-rooted and had never resolved; the data they name also
moved out of the scripts tree on that date.
Resolved through hat_observed_rates.py since 2026-09-18.
```

```text
--- CASCADE domain geometry ---
GeoJSON exported from ArcGIS
```

```text
Coordinate Reference System for distance calculations.
UTM Zone 18N is correct for the NC Outer Banks.
```

```text
Buffer (in degrees) added around the domain bounding box when
pre-filtering the global transect file. 0.1 deg ~ 10 km.
```

```text
Nearest-feature snapping is disabled — Hatteras Island's curvature
causes cross-curve mismatches. Only point-in-polygon matches are kept.
MAX_SNAP_DISTANCE_M is retained for reference but not used.
```

```text
Extract transect origin point from LineString and make it the
ONLY geometry column (avoids the 'active geometry not present' error)
```

```text
Fix invalid CRS — EPSG:3725 is not a real code; domains exported
from ArcGIS in UTM 18N should be EPSG:32618
```

<details><summary>Function notes (the original docstrings)</summary>

**`fix_crs()`**

```text
Validate the CRS of a GeoDataFrame. If the CRS is missing or
unrecognised (e.g. EPSG:3725 which is not a real EPSG code),
replace it with fallback_crs and warn the user.

This commonly happens when ArcGIS exports GeoJSON with a non-standard
CRS code stored internally that geopandas cannot resolve.
```

**`load_transects()`**

```text
Load CoastSat transect geometry.
If geometries are LineStrings, extracts the first point (transect origin)
and sets it as the active geometry for point-in-polygon joining.
The original LineString geometry column is dropped to avoid confusion.
```

**`clip_transects_to_study_area()`**

```text
Pre-filter the global transect dataset to only those within a
buffered bounding box around the domain extent.

Essential when working with the full global CoastSat GeoJSON
(230,000+ transects) — without this the spatial join is very slow.
```

**`spatial_join_polygon()`**

```text
Join transect origin points to domain polygons — point-in-polygon only.

Nearest-feature snapping is intentionally disabled: Hatteras Island's
shoreline curvature means that snapping unmatched transects to the nearest
domain frequently assigns them to the wrong domain across the curve.
Transects that do not fall within any domain polygon are left as NaN
and excluded from downstream LRR calculations.
```

**`make_verification_map()`**

```text
Quick-look map: transect origins coloured by match method,
domain polygons shown in light yellow with domain ID labels.
Inspect this carefully before using the lookup table.
```

</details>
