# example — step 1 on 35 real transects

Run it from this folder:

```
python ../transect_zone_join_template.py --out expected_output
```

## What is in it

```
data/
    transects.geojson    35 CoastSat transect lines (id = time-series filename), EPSG:4326
    zones.geojson        3 zone polygons, zone_id 83, 84, 85, EPSG:4326
expected_output/         what the command above writes
    transect_zones.csv            the lookup step 2 reads
    transect_zones_problems.csv   transects left out, and why
    transect_zones_map.png        the check map
```

The zones are three of the ~500 m CASCADE model domains on northern Hatteras
Island (83–85). The transects are every CoastSat transect whose origin falls in
them, `usa_NC_0036_0034` through `usa_NC_0036_0066`, plus one more beyond each
end (`_0033`, `_0067`).

## What you should see

- **33 matched, 2 no zone.** The two extra transects sit in domains 82 and 86,
  which are not in `zones.geojson`, so the join refuses them instead of snapping
  them into the nearest zone. They are the red crosses on the map and the two
  rows of the problems file. That is the behaviour to check for in your own data.
- The 33 matches agree one for one with the production join
  (`../../../2-transect-frame/coastsat_domain_mapping.py`).

If your run differs from `expected_output/`, the difference is in your
environment (geopandas/shapely version), not in the data.

Transects from CoastSat (coastsat.space); domain boxes from this project.
