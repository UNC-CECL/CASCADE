# 1-observations — what was measured, not what we fitted

The shoreline record as it arrives: digitized lines, the NC Coastal Management
historical set, and the CoastSat satellite chainage. Nothing here fits a rate.
That is `../3-rates/`, and it cannot run until `../2-transect-frame/` has tied
these observations to domains.

```
shoreline_inventory/
    HAT_shoreline_inventory.py
        Cross-source inventory for the whole study area: what exists from each
        source, how the sources overlap in time, and where the gaps are. The
        study area comes from a spatial filter file, so the same script can be
        pointed at a narrower reach (the Buxton groin transects, say) without
        being edited.
        Writes 1-observations/shoreline_inventory/shoreline_position_output/.

shoreline_patterns/
    HAT_shoreline_trajectory_classification.py
        Is a domain eroding steadily, stable, or reversing? Classifies each
        domain's trajectory from the CoastSat time series over 1984-2004,
        2004-2024 and the full record.
    HAT_trajectory_map.py
        The same classification, and LRR magnitude, drawn on the island.
        USE_SATELLITE = True wants an Esri basemap and therefore the internet;
        False gives a plain ocean background and works offline.
```

Both trajectory scripts write under `1-observations/shoreline_patterns/`. Their
earlier outputs were deleted as stale on 2026-09-18 — the folder being empty
does not mean the scripts are broken, only that nobody has rerun them since.

Data paths are resolved through `scripts/site_layer/hat_observed_rates.py`;
none is typed here.
