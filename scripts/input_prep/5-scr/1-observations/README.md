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

## detrended_position

Why do the rate products disagree between windows? Not because the fits are
noisy: there is one signal the whole island shares, on top of every transect's
own trend.

```
coastsat_detrended_position.py
    Detrends every CoastSat transect against its OWN 1996-2024 fit, reduces it
    to an annual median and averages all 906, so anything local cancels. What
    survives is a landward sag through 2005-2020 and a +16.8 m step in the
    single year to 2021. That excursion alone biases the fitted rate by
    -0.37 m/yr over 1996-2010 and +0.88 m/yr over 2010-2024, which is the
    answer to why no short window recovers the long-term rate. Writes the
    detrended matrix the next script reads, so the detrending happens once.

coastsat_position_attribution.py
    What the step IS, in five tests, each stored with its verdict including
    the two that failed: per-transect noise (rejected), nourishment (rejected
    as the cause, confirmed as a visible signal), seasonal sampling (rejected),
    the independent dune line (corroborates), spatial structure (corroborates).
    The step is real, so the grading target's sensitivity to 2021-2024 is a
    question about which period the model should represent, not about which
    data to trust. Reads stored endpoint tables for the dune-line test and
    refits nothing.
```

Built 2026-09-23, out of `3-rates/coastsat/window_convergence/`.

