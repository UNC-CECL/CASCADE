# 1-observations — what the groin did to the shoreline

Measured or digitized. Nothing here is fitted to the model, and nothing here changes when the model does. Start with `structure_history.md`.

```
structure_history.md       the structure's dates and the fillet in numbers
gis_data/                  GIS layers: groins_hatteras, wet_dry_groin, transects_100m,
                           the offshore datum line, the CASCADE domains and study area
wetdry_photo_positions/    shoreline position per domain (GIS 2-12) from 24 dated wet/dry
                           surveys and the dune lines, as change since 1967; the
                           Change_from_wetdry table is the observed gap every fit reads
coastsat_shoreline/        CoastSat and the 100 m grid around the field: rates before,
                           during and after the groin's working life, profiles, GIFs
coastsat_groin_condition/  the GIS 5|6 gap from CoastSat, year by year, and a break-year
                           fit: the gap stopped widening in 1995
figures/                   groin_two_shorelines: the two shorelines and the gap
```

**Who reads these.**
- `wetdry_photo_positions/Change_from_wetdry_1967_D2_D12.csv`: every hindcast fit, the figures, `scripts/hatteras_ms/groin-sweep/HAT_groin_sweep_config.py`.
- `gis_data/groins_hatteras.geojson`: `scripts/input_prep/4-mgmt-forcings/beach_nourishment.py`.
- `gis_data/transects_100m.geojson`: the 2-brie-offset stage (`data/hatteras_init/2-brie-offset/transects/README.md`).
- `coastsat_shoreline/shoreline_output_coastsat/groin_analysis_chainage_all.csv`: the dipole sweep's target (`scripts/hatteras_ms/groin-sweep/`).
- `coastsat_groin_condition/`: the 2026-10-08 schedule refit (`../3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/`).

The photos and CoastSat disagree at one date, the 2004 survey (by 42 m); see `coastsat_groin_condition/README.md`, figure 7.
