# forecast - not built

Reserved for forecast-period domains, expected to come from a 2025 DEM.

Nothing here yet. To build it:

1. Add the DEM as a product in `0-elevation/` and register it in
   `scripts/hat_elevation_products.py`.
2. `HAT_export_to_numpy.py --product <dem> --target forecast`
3. Set `TOPO_PRODUCT = "forecast"` and a `VERSION` in
   `HAT_dune_topo_extractor.py`, then run a pick pass.
4. Add a `forecast` entry to `HATTERAS_PERIODS` carrying
   `"topo_product": "forecast"`.

`"forecast"` is already in `PRODUCTS` in `scripts/hat_topo_version.py`, so the
resolver will accept the name and fail with a useful message until the
directory exists.
