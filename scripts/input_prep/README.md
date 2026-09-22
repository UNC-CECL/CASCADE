# input_prep - building what the model eats

One folder per stage, numbered to match `data/hatteras_init/`. A stage's code
is here; everything it reads or writes is in the data tree under the same
number. Inside a stage, the code is split into numbered STEPS in the order the
work runs, and where the data tree is grouped the same way the names match.

```
0-elevation/          the DEMs, from survey to 10 m domain rasters
    1-source-selection/  2-produce/  3-figures/
1-barrier3d-domains/  the extraction: interiors and dune arrays per domain
    1-extraction/  2-domain-reconstruction-1984/ (1-measurement .. 6-result)
2-brie-offset/        where each domain sits cross-shore at model year zero
    1-produce/        duneline_to_raw_offsets -> island_offset_hybrid, and
                      build_island_offset.py, the one command that runs both
    2-figures/        version comparison, the 1:1 offset profile
    old-linear-bridge/  the retired 1984 linear-bridge variant
3-env-forcings/       sea level, storms, waves  (data: 1-records/ 2-rslr/ 3-storms/)
    1-records/        the gauge downloader, the hurricane record figure
    2-rslr/           the sea-level fit
    3-storms/         historical_storm_creation_v3_HAT.py (the model input),
                      storm_check/ and storm_validation/ (both validators),
                      from_Hannah/ from_lexi/ from_roya/ (earlier generators,
                      kept for provenance; the from_Hannah figure scripts
                      name paths that no longer exist)
4-mgmt-forcings/      NC-12 and nourishment
    road_offset/ (1-produce .. 4-compare)  road_elevation/  road_relocation/
    road_relocation/from_roya/   Roya's Pea Island original of the
                                 relocation measurement (was roya_files/)
5-scr/                the observed shoreline record and its rate fits
    1-observations/  2-transect-frame/  3-rates/  4-comparisons/
    lib/              the shared OLS, and scr_paths.py (the module resolver)
    tools/            the WINDOWS.md generator, the pipeline verifier
6-scr-smooth/         smoothing the observed rates
7-source-sink/        the background-erosion calibration
    1-prepare/  2-calibrate/  3-figures/  4-export/
8-overwash-analysis/  the observed overwash record and its figures

HAT_units_datum_check.py   a cross-stage check: are the inputs in the units and
                           datum the model assumes (m, MHW vs NAVD88)?
superseded_20260902/       the parametric source/sink attempts; see its WHY.md
```

Note the plural: the code folder is `4-mgmt-forcings`, the data folder is
`4-mgmt-forcing`. That is a spelling accident, not a distinction; it is left
alone because too many paths spell it.

`5-scr/` mirrors its data tree (`1-observations/ .. 4-comparisons/`) since
2026-09-22. It had kept producer folders (`CoastSat/`, `duneline_endpoint/`,
...) because the `CoastSat/` scripts imported `coastsat_lrr_analysis.py` as a
sibling, so splitting them by job meant moving that module first. It now sits
in `5-scr/lib/coastsat_lrr.py`, and `5-scr/lib/scr_paths.py` is the one place
that knows where any shared 5-scr module lives. See `5-scr/README.md`.

## Naming: two conventions, and where the line is

Scripts are **bare** in the shoreline chain and **`HAT_`-prefixed** everywhere
else. That is a live split, not drift, and this is the boundary:

| bare | | `HAT_` prefixed | |
|---|---|---|---|
| `5-scr/` | 27 | `0-elevation/` | 10 |
| `7-source-sink/` | 10 | `1-barrier3d-domains/` | 36 |
| `8-overwash-analysis/` | 3 of 4 | `4-mgmt-forcings/` | 16 of 18 |
| `2-brie-offset/` | 3 of 5 | `3-env-forcings/` | 10 of 12 |
| `6-scr-smooth/` | 2 | | |

`5-scr` dropped the prefix on 2026-09-22 for a hard reason: `scr_paths.py` puts
five of its folders on `sys.path` at once, so its filenames are flat module
names, and `import HAT_rates_figures` reads badly where `import rates_figures`
does not. `6-scr-smooth` and `7-source-sink` followed the same day for a softer
one — they are read together with `5-scr`, which feeds them, and a reader
should not have to remember that one link in a chain spells its scripts
differently.

Nothing forces the rest. `1-barrier3d-domains` alone is 36 scripts, and
renaming them buys consistency at the price of a large sweep through code that
loads several of them by path. **If the split is ever closed, close it toward
bare** — that is the direction the four resolver modules in `site_layer/`,
`cascade_pipeline/` and `repo_tools/` already point, none of which is
prefixed.

Environment variables keep `HAT_` regardless (`HAT_BE_BASE_PRESET`,
`HAT_RUN_KIND`, `HAT_GEOMETRY`): they share one namespace with every other
program on the machine, which is exactly the case a prefix is for.

## Moving a script

Every script here finds the repository by searching upward for
`pyproject.toml`, and the data through `scripts/site_layer/`, so a script can
change depth without breaking. What does break is another script that loads it
BY PATH: `build_island_offset.py` runs its three steps that way, and
`be_zone_residual_fit.py`, `HAT_road_offset_from_dune_start.py`,
`HAT_road_placement_on_domains.py` and `HAT_dune_topo_extractor.py` are each
loaded by several others. Search for the file name before moving one.

The exception is `5-scr/`, where the six shared modules are registered in
`5-scr/lib/scr_paths.py`: moving one is editing its row there, and
`python 5-scr/lib/scr_paths.py` reports any row that no longer matches disk.

(2026-09-18: 2-brie-offset and 3-env-forcings were split into steps and
roya_files/ folded into road_relocation/. 2026-09-22: 5-scr was split to
mirror its data tree, and its by-path loads replaced by scr_paths.)

## Before running any of these

Most take the period as an argument and derive every path from it. Prefer that
to editing a literal: a folder name and a file name that are typed separately
are a pair that can disagree, and have.
