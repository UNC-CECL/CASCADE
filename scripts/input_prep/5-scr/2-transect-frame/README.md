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
