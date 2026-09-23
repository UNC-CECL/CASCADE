# 1-observations — what was measured, not what we fitted

The shoreline record as it arrives: digitized lines, the NC Coastal Management
historical set, and the CoastSat satellite chainage. Nothing here fits a rate.
That is `../3-rates/`, and it cannot run until `../2-transect-frame/` has tied
these observations to domains.

```
shoreline_inventory/
    shoreline_inventory.py
        Cross-source inventory for the whole study area: what exists from each
        source, how the sources overlap in time, and where the gaps are. The
        study area comes from a spatial filter file, so the same script can be
        pointed at a narrower reach (the Buxton groin transects, say) without
        being edited.
        Writes 1-observations/shoreline_inventory/shoreline_position_output/.

mean_shoreline/
    coastsat_mean_shoreline.py
        One averaging window's MEAN satellite shoreline, as a line on the
        ground: each CoastSat transect's chainage averaged over the window,
        geolocated, and the mean points strung into a single polyline that
        2-brie-offset intersects with the 100 m transect frame exactly as it
        intersects a digitized dune line.

        It is the only script in 5-scr that needs a POSITION rather than a
        difference, and that is the whole job. Every other CoastSat product
        differences chainage, so each transect's arbitrary origin cancels;
        here it does not. Aggregated to the 90 domains, raw chainage spans
        124 m alongshore while the geolocated position spans 6222 m -- the
        origins follow the shore around the cape.
        Writes 1-observations/mean_shoreline/<start>_<end>/.

shoreline_patterns/
    shoreline_trajectory_classification.py
        Is a domain eroding steadily, stable, or reversing? Classifies each
        domain's trajectory from the CoastSat time series over 1984-2004,
        2004-2024 and the full record.
    shoreline_trajectory_map.py
        The same classification, and LRR magnitude, drawn on the island.
        USE_SATELLITE = True wants an Esri basemap and therefore the internet;
        False gives a plain ocean background and works offline.
```

Both trajectory scripts write under `1-observations/shoreline_patterns/`. Their
earlier outputs were deleted as stale on 2026-09-18 — the folder being empty
does not mean the scripts are broken, only that nobody has rerun them since.

Data paths are resolved through `scripts/site_layer/hat_observed_rates.py`;
none is typed here.
