# Beach nourishment

Where and when the beach was filled on Hatteras, what the model receives, and checks on where the sand went.
Nothing here is a model input. The runs read the fill list in `HATTERAS_NOURISHMENT_PROJECTS`
(`scripts/site_layer/hatteras_site_config.py`). Each folder keeps its PDFs and `CAPTIONS.md` under `supporting/`.

## Which file for what

| I need... | Look at |
|---|---|
| The one table of every fill the model uses (dates, extent, source, domains, volume) | `2-model-input/nourishment_summary_table.csv` (`.png` beside it) |
| A map of where the fills are across Hatteras | `3-maps/nourishment_domain_map.png` |
| The paper map figure | `3-maps/nourishment_maps_paper.png` |
| A close-up of one fill | `3-maps/per_fill/nourishment_map_<year>_<town>.png` |
| When and where, as a chart | `2-model-input/nourishment_when_where.png` |
| How much sand per metre each domain gets | `2-model-input/nourishment_volume_alongshore.png` |
| The original records | `1-sources/` (see its README) |
| Where the reported ends of each fill fall in the model | `4-extent-checks/reported_limits/` |
| Where CoastSat saw the shoreline move after each fill | `4-extent-checks/coastsat/` |

## Folders

```
1-sources/          the records: national database sheets, Hatteras_BN_data.xlsx,
                    national_beach_nourishment_database.csv, placement dates; README compares them
2-model-input/      what the runs receive: summary table, nourishment_projects.csv (every fill,
                    on and off the modelled reach), the when/where chart, m3/m per domain
3-maps/             the whole-island domain map and the paper map; per_fill/ one map per fill,
                    including the 2026 fills that are not in the model
4-extent-checks/
  reported_limits/  each reported limit geolocated and placed in a model domain
  coastsat/         the observed footprint after each fill; shorelines/ and animations/ per fill
```

## Producers

All in `scripts/input_prep/4-mgmt-forcings/`:

- `beach_nourishment.py`: `2-model-input/` and `3-maps/`
- `nourishment_reported_extent.py`: `4-extent-checks/reported_limits/`
- `nourishment_extent_coastsat.py`: `4-extent-checks/coastsat/`
- `nourishment_shorelines_coastsat.py`: `4-extent-checks/coastsat/shorelines/`
- `nourishment_animation_coastsat.py`, `nourishment_animation_map_coastsat.py`: `4-extent-checks/coastsat/animations/`

Reorganized 2026-10-06 from `datasets/`, `project_maps/`, `reported_extent/`, `coastsat_check/` and loose top-level files.
