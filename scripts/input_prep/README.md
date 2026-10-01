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
    superseded_20260914/  the retired 1984 linear-bridge variant
                          (was old-linear-bridge/ until 2026-09-22)
3-env-forcings/       sea level, storms, waves  (data: 1-records/ 2-rslr/ 3-storms/)
    1-records/        the gauge downloader, the hurricane record figure
    2-rslr/           the sea-level fit
    3-storms/         historical_storm_creation_v3_HAT.py (the model input),
                      storm_validation/ (the validator; storm_check/ held
                      the two it replaced and was retired 2026-09-22),
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

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_units_datum_check.py

Units and vertical datum of every elevation-like quantity the hindcast hands Barrier3D.

From the script's original header:

```text
Verify the units and vertical datum of every elevation-like quantity the
Hatteras CASCADE hindcast feeds to Barrier3D, by tracing the source and by
checking the actual data files.

THE THING THIS EXISTS TO CATCH
"Everything in CASCADE is relative to MHW" is true of what the model COMPUTES
ON, and false of what you SUPPLY. There are two conventions side by side:

  CONVERTED FOR YOU -- supply in metres NAVD88, load_input.py converts:
      MHW          load_input.py:227    /10
      BermEl       load_input.py:241    /10 - MHW
      Dmaxel       load_input.py:304    /10 - MHW
      ShrubEl_*    load_input.py:346-7  /10 - MHW
      Dstart       load_input.py:240    /10

  ALREADY CONVERTED -- supply in the model's own frame, nothing touches it:
      elevation_file (.npy)   dam, MHW-relative
                              load_elevation() only does np.load; configuration.py
                              documents it as "[dam x dam x dam MHW]"
      dune_file (.npy)        dam, height ABOVE BERM
                              load_input.py:249 assigns DuneStart straight into
                              DuneDomain with no conversion
      road_ele                metres, MHW-relative
                              bulldoze() only does road_ele/dz

barrier3d.py:1197 pops MHW into self._MHW -- and _MHW appears NOWHERE ELSE in
the file. Barrier3D never applies it to a grid. Anything grid-shaped has to
arrive pre-converted.

That asymmetry is what made ROAD_ELEVATION = 1.45 ambiguous: it sits in the one
scalar family that takes pre-converted values while looking like the ones that
do not.

WHAT IS CHECKED
  1. The contract, as a table: supplied-as / converted-where / model-sees.
  2. The saved arrays, against the range each convention implies. A 10x unit
     error or a 0.36 m datum slip is far outside the plausible band, so this
     catches both.
  3. Berm elevation by two independent code paths (the extractor's and
     load_input's) -- they must agree.
  4. Road elevation against the interior elevation at the road -- ONCE PER
     PERIOD, since both halves are period-specific: the setback decides which
     rows to read, the topography decides what is in them.

     The two periods have different expected answers, and that is the whole
     content of the check. RoadElevation.csv is sampled on the 2009-2014
     baseline, which is the DEM behind 2004-start, so 2004 compares a surface
     against itself and must agree to ~0. 1984 runs 2009-2014-1996, where the
     1996 ALACE survey overwrites the corridor and its vertical offset is left
     uncorrected by design -- so 1984 is expected to sit HIGH by exactly that
     offset, which HAT_dem_1984_mosaic.py measured on the 2009/1996 overlap
     and writes to mosaic_1984_audit.csv every run.

     So the tolerance is applied to `gap - expected`, not to `gap`. That tests
     a claim -- "all of this gap is the 1996 survey offset" -- rather than
     widening a band until the number fits. Measured: 1984 gap +0.26 m against
     a recorded offset of +0.26 m, residual -0.001 m; 2004 residual -0.000 m.
     A missing or doubled MHW subtraction still shows up as a ~0.36 m residual
     in BOTH periods, because `expected` comes from a different file than the
     road elevation does.
  5. Whether each runner constant actually reaches the model, or is overridden.

REQUIREMENTS
  numpy (pandas optional, for the road-elevation cross-check)
```

Notes that were in the code:

```text
scripts/input_prep/HAT_units_datum_check.py -> repo root is parents[2].
This said parents[4], which resolved to C:\Users\hanna: every path below
pointed outside the repo and the script could not find one of its inputs.
```

```text
Resolved from the extractor rather than pinned, so this checks the units of
the arrays actually being run. It said "2009_v2", which has since been moved
to 2009-dune-topo/incorrect/.
```

```text
BOTH PERIODS, AND THE FILES THE RUNNER ACTUALLY SPENDS (2026-08-26).

Two things were wrong here, and both were silent.

(1) ONE PRODUCT. `topo_dirs()` with no argument resolves DEFAULT_PRODUCT,
i.e. 2004-start. Since 1-barrier3d-domains went period-first there are
two extractions and the 1984 period runs the other one - so this script
reported "the units of the arrays actually being run" having never
opened half of them. All 90 domains differ between the products and 65
differ in interior SHAPE, so it is not a formality: the road-vs-interior
gap below indexes a specific row of a specific array.

(2) THE LEGACY SETBACK FILE. ROAD_SETBACK_CSV pointed at
old_method_offset/2004/RoadSetback_2004.csv. The runner switched to the
dune-start method on 2026-08-18 (hatteras_site_config.py), and the two
put the road a median ~22 m apart in 2004. The most sensitive check in
this file - road elevation against the interior AT THE ROAD - was
therefore reading the interior at rows the model does not bulldoze.

Both are fixed by asking hatteras_site_config, which is what the runner
reads, instead of mirroring it. The header above says this CONFIG "must
mirror HAT_hindcast_1984_2024.py"; mirroring is how it drifted, so the
period-dependent forcings are now imported and cannot.
```

```text
One file, both periods -- not because there is one DEM (there are two) but
because they differ under the road only by the uncorrected 1996-vs-2009
survey offset. See HATTERAS_ROAD_ELEVATION_FILE in hatteras_site_config.py.
```

```text
WHERE RoadElevation.csv WAS SAMPLED. One file serves both periods and it is
built on the 2009-2014 baseline (HAT_road_elevation.py, FILL_SOURCE), which
is the DEM behind 2004-start. So for 2004 the road elevation and the interior
under the road are the SAME surface and must agree to ~0.
```

```text
The 1984 period does not run that surface. 2009-2014-1996 overwrites measured
ground wherever the 1996 ALACE survey has data, including through the road
corridor, and HAT_dem_1984_mosaic.py leaves the vertical offset UNCORRECTED
on purpose ("bias correction OFF, feathering OFF"). It writes the offset it
measured, per domain, to this file every run.
```

```text
Was "domain_*_topography_*.npy". The trailing _* required a year tag
that no longer exists, so the glob matched nothing after 2026-08-26.
```

```text
Expected: dam MHW-relative. Sentinel is SENTINEL_WATER_M/10 dam exactly.
Real barrier tops out a few metres above MHW -> a few tenths of a dam.
```

```text
Was "domain_*_dune_*.npy" - the trailing _* needed a year tag that no
longer exists, so this matched nothing after 2026-08-26.
```

```text
UNITS: dam above berm. If these were metres the median would be ~10x and
land outside this band; that is what the check discriminates. It says
nothing about whether the values are physically sensible.
```

```text
PLAUSIBILITY: separate question, and NOT a units problem. NC-12's
artificial dune ridge is typically 3-5 m NAVD88. The extractor's own
docstring warns the search window can be "wide enough to catch a
back-dune, a wooded ridge, or a house", which is what a high crest means.
```

```text
THE SENSITIVE ONE: road elevation vs the interior it is written into.

Per period, because both halves of it are period-specific: the setback
says WHICH ROWS to read and the topography is WHAT IS IN THEM. This used
to run once, against the legacy 2004 setback and the default product -
so it compared the 2004-start interior at rows the old method chose, for
a model that runs neither.
```

```text
WHAT GAP SHOULD THIS PERIOD SHOW?

For the period whose product IS the surface RoadElevation.csv was
sampled from, zero: same LiDAR, same corridor, so anything left is a
datum or unit error, which is what this check is for.

For 1984 it is NOT zero, and pretending otherwise would make this
check fail forever for a reason that is documented and deliberate.
The expected value is not invented here either - it is the offset
HAT_dem_1984_mosaic.py measured on the 2009/1996 overlap and wrote to
its own audit. Testing `gap - expected` therefore tests a specific
claim ("the whole gap is the 1996 survey offset") rather than
widening a tolerance until the number fits inside it. If the graft
is ever bias-corrected, or the road elevation is rebuilt on the 1984
product, expected goes to ~0 on its own and this keeps working.
```

<details><summary>Function notes (the original docstrings)</summary>

**`recorded_survey_offset()`**

```text
Median 1996-minus-2009 offset over `domains`, from the mosaic audit.

READ, NEVER ASSUMED. The point of returning it rather than hardcoding a
tolerance is that the road-vs-interior gap below can then be tested against
a number measured independently, by a different script, from the overlap of
the two surveys - instead of being written off as "about right".

Sign in the CSV is base - fill, i.e. 2009 - 1996, so it is negated here to
read as "1996 sits this much higher than 2009".
```

</details>
