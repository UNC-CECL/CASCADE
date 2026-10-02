# figure_making/pipeline — how each model input is built

One folder per model-input step, numbered as `data/hatteras_init/` is. Each
script draws how its input is made from the producer's own files or code, and
writes to `output/figures/3-model-inputs/<step>/`. All are redrawn by
`tools/regenerate_all_figures.py`.

```
0-elevation/dem_composition_figures.py            which survey supplies each part of the DEM
1-barrier3d-domains/domain_extraction_figures.py  how one domain's arrays are cut from the DEM
2-brie-offset/offset_build_figures.py             how the BRIE offset is built from the dune line
3-storms/storm_construction_figures.py            how the storm series is built (reproduced, checked)
4-mgmt-forcings/road_setback_figures.py           how the NC-12 setback is measured
5-scr/observed_target_figures.py                  how the observed shoreline target is built
7-source-sink/be_method_figures.py                how the source/sink end rates are solved
inputs_overview_figure.py                         every input the 1996 and 2010 runs read, one page each
```

What not to trust: nothing here regenerates an input. The extraction and
storm figures stop if their reproduction differs from the saved product, so a
figure that draws is a figure that matches what the model reads.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 0-elevation/dem_composition_figures.py

Which survey supplies each part of the model's topography, and how 1 m surveys become the 10 m grid.

From the script's original header:

```text
Which survey supplies each part of the model's topography, and how the 1 m
surveys become the 10 m grid Barrier3D reads.

    python scripts/figure_making/pipeline/0-elevation/dem_composition_figures.py

Writes to output/figures/3-model-inputs/0-elevation/:

    dem_sources_alongshore.png   per GIS domain, the share of measured cells
                                 each survey supplies, in both products:
                                 2009-2014 (the 2004/2010 starts) and
                                 2009-2014-1996 (the 1984/1996 starts)
    dem_resample_one_domain.png  one domain (GIS 45): the survey each 1 m cell
                                 comes from, the 1 m surface, and the 10 m
                                 resample the domain arrays are cut from

Read only: the products' clip_domain_<N>_{filled,survey}.tif (1 m) and
resampled_domain_<N>_{filled,survey}.tif (10 m), resolved through
scripts/site_layer/hat_elevation_products.py. Nothing is regenerated.

The tiles are north-up UTM boxes, 500 m alongshore by 2000 m cross-shore;
through the reach the ocean is to the east, so it sits at the right.
```

Notes that were in the code:

```text
survey code -> label, colour. The vintage pair: the earlier survey red, the
later blue; the 2009 base neutral.
```

### 1-barrier3d-domains/domain_extraction_figures.py

How one domain's Barrier3D arrays are cut from the 10 m DEM, step by step, with the extractor's own code.

From the script's original header:

```text
How one domain's Barrier3D arrays are cut from the 10 m DEM, step by step,
using the extractor's own functions on its own saved picks.

    python scripts/figure_making/pipeline/1-barrier3d-domains/domain_extraction_figures.py [--gis 45] [--product 2004-start]

Writes output/figures/3-model-inputs/1-domains/domain_extraction_gis<N>.png.

WHAT IT CALLS
    HAT_dune_topo_extractor.load_profiles (orient, MHW, clamp, beach start,
    straighten, trim), find_dunes (the crest inside the picked window) and
    build_interior (everything landward of the crest), with the window from
    the product's picks file. extract_domain() is NOT called: it writes the
    arrays. Instead the result is checked against the saved arrays of the
    product's CURRENT version, and the script stops if they differ, so the
    figure shows what the model reads. The road overlay is switched off.

Every cross-shore panel has the ocean at the RIGHT.
```

### 2-brie-offset/offset_build_figures.py

How the BRIE shoreline offset is built from a digitised dune line, step by step, for the 1996 and 2010 starts.

From the script's original header:

```text
How the BRIE shoreline offset is built from a digitised dune line, step by
step, for the builds the current runs read (1996 and 2010 starts).

    python scripts/figure_making/pipeline/2-brie-offset/offset_build_figures.py

Writes output/figures/3-model-inputs/2-brie-offset/offset_build_<year>.png.

THE STEPS DRAWN (the producers, not re-implemented)
    1. duneline_to_raw_offsets.py intersects the dune line with the 100 m
       transects. Each transect starts on the offshore datum line and runs
       west across the island; the STATION of a crossing is its distance
       along the transect from the datum. The raw file the build kept is read,
       not recomputed.
    2. island_offset_hybrid.py averages the ~5 transects of each 500 m domain
       and subtracts the smallest domain mean, so the offset is 0 at the most
       seaward domain and positive landward.
    3. cascade_pipeline.hindcast.pad_offset_ring closes BRIE's periodic line
       with 15 buffer domains per side (a cubic Hermite from GIS 90 back
       round to GIS 1). The padded file IS what the runner hands Cascade.

Every path resolves through site_layer.hat_topo_version: the dune-line vintage
from DUNE_LINE_FOR_YEAR, the build from offset_version (env > CURRENT > the
only v<n>). The ocean is on the right in the plan panels (easting across).
```

### 3-storms/storm_construction_figures.py

How the model's storm series is built: the generator's rule reproduced, checked, and drawn.

From the script's original header:

```text
How the storm series the model reads is built, drawn from the same records by
the same rules as the generator
(scripts/input_prep/3-env-forcings/3-storms/historical_storm_creation_v3_HAT.py).

    python scripts/figure_making/pipeline/3-storms/storm_construction_figures.py

Writes to output/figures/3-model-inputs/3-forcing/:

    storm_construction_steps.png   the chain on a worked stretch (Edouard and
                                   Fran, late August - early September 1996):
                                   Duck water level and WIS waves -> Stockdon
                                   run-up -> total water level against the berm
                                   -> hours above it grouped into events ->
                                   events split where the storm hours break ->
                                   each event trimmed to 24 h around its peak
    storm_events_by_duration.png   every event of 1996-2024 by its length above
                                   the berm and its Rhigh, kept whole or
                                   trimmed; and the storm hours each year
                                   against the hours the model receives

THE RULE (the series in use, v3_split12_trim24, adopted 2026-09-29)
    1 hours when total water level exceeds the berm are storm hours
    2 storm hours less than 24 h apart are grouped into one event
    3 an event is split wherever consecutive storm hours are >= 12 h apart; a
      piece shorter than 8 h is folded into the piece before it (after, for
      the first)
    4 an event of fewer than 8 storm hours is not a storm
    5 an event longer than 24 h is cut to the 24 storm hours centred on its
      peak total water level; Rhigh, Rlow and the period come from what is kept

THE REPRODUCTION
    The generator is a script with module-level execution, so its logic is
    re-implemented here (build_events) and CHECKED before any figure is drawn:
    the events must equal, row for row, the committed
    <window>_storms_v3_split12_trim24_summary.csv of 1996_2010 and 2010_2024.
    The run stops if they do not.

REWORKED 2026-09-29 (Hannah: "rework the construction figures"). Until then
this drew the v3_72 rule -- events over 72 h dropped whole -- which has not
been the model's input since 2026-09-28 (trim24) and did not have the split.
```

Notes that were in the code:

```text
the generator's inputs (historical_storm_creation_v3_HAT.py, "user inputs",
and the command line of the series in use)
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_merged()`**

```text
load_data(): an hourly index over the window, the Duck water level and
WIS Hs/Tp placed on it, rows with any gap dropped.
```

**`_split()`**

```text
The generator's split: cut where consecutive storm hours are >=
SPLIT_GAP_H apart, fold a piece under MIN_DUR_H into its neighbour.
```

**`build_events()`**

```text
create_storms() with the series-in-use options. Returns the hourly
table (with the grouped system id) and EVERY piece, kept or not, with the
hours the trim keeps. Units as the generator: Rhigh/Rlow in dam above MHW,
duration in hours above the berm.
```

</details>

### 4-mgmt-forcings/road_setback_figures.py

How the NC-12 road setback is measured on one domain, and what the 1996 and 2010 runs are handed.

From the script's original header:

```text
How the NC-12 road setback the model reads is measured, on one domain, and
what the 1996 and 2010 runs are handed along the reach.

    python scripts/figure_making/pipeline/4-mgmt-forcings/road_setback_figures.py

Writes output/figures/3-model-inputs/4-management/:
    road_setback_measurement.png   one domain: raw grid + rasterised road,
                                   the straightened profiles, the per-profile
                                   setbacks and the domain value
    road_setback_inputs.png        the setbacks the 1996 and 2010 runs read,
                                   the measurements they come from, and the
                                   corrections applied on the way

Everything is read from the producers' saved products (nothing re-measured):
    raster/<line vintage>/masks/domain_N_road_<vintage>.npy
        HAT_rasterize_road_to_domains.py: the NC-12 centreline buffered to
        the road width and burned onto each domain's raw 10 m grid
    dunestart_offset/measured/<year>/RoadOffset_<year>_profiles.csv, _domains.csv
        HAT_road_offset_from_dune_start.py: the mask sheared with the same
        per-profile shear as the topography, then per profile the distance
        from interior row 0 (one cell landward of the picked dune crest) to
        the road's seaward edge; the domain value is the median, then the
        negative floor (ocean side) and the drowning-road move (bay side)
    dunestart_offset/derived/<1996|2010>/RoadSetback_*_dunestart.csv
        HAT_road_setback_derived_vintages.py: 1996 = the 1984 measurement +
        the 1989 Pea Island relocation; 2010 = the 2004 measurement.
Paths resolve through site_layer.hat_topo_version. Ocean on the right.
```

Notes that were in the code:

```text
landward of the road's most seaward cell, in both frames (raw columns
grow toward the ocean; straightened cells grow away from it)
```

### 5-scr/observed_target_figures.py

How the observed shoreline target is built: CoastSat series, per-transect LRR, domain means, LOWESS.

From the script's original header:

```text
How the observed shoreline-change target the hindcast is scored against is
built from CoastSat, for the two current windows (1996-2010, 2010-2024).

    python scripts/figure_making/pipeline/5-scr/observed_target_figures.py

Writes output/figures/3-model-inputs/5-observed-target/observed_target_<window>.png.

THE STEPS DRAWN (the producers' own functions, not re-implemented)
    1. One transect's CoastSat shoreline positions (chainage, + seaward) in the
       calendar window, and the OLS slope through them: the LRR
       (5-scr/lib/coastsat_lrr.compute_lrr; the calendar-window filter of
       coastsat_domain_lrr.py, START y0-01-01, END y1-12-31). The slope drawn
       is checked against transect_lrr_full.csv.
    2. Transects grouped into their GIS domain (transect_domain_lookup.csv)
       and averaged: the raw domain mean.
    3. LOWESS at transect resolution over a 7-domain (3.5 km) window, averaged
       back to domains, with GIS 1-10 kept as raw domain means
       (cascade_pipeline.coastsat_lowess.spliced_lowess_series, the same two
       steps hindcast.build_target_table applies to the scoring target).
```

### 7-source-sink/be_method_figures.py

How the source/sink (BE) end rates are solved and what zone field the runs carry.

From the script's original header:

```text
The source/sink (background erosion, BE) method figures for the current
1996-2010 / 2010-2024 pair, which existed only for 1984/2004.

    python scripts/figure_making/pipeline/7-source-sink/be_method_figures.py

Writes output/figures/3-model-inputs/7-source-sink/:
    be_end_solve.png   how the two end values the current (edgeBE) runs carry
                       were solved: residual per Newton step at GIS 1 and 90,
                       and the direct probes that set 2010 GIS 1
    be_zone_field.png  the zone-by-zone calibration (calibBE) for the pair:
                       residuals, which domains were eligible, and the rates.
                       Made 2026-09-18, BEFORE the metres offset and wave
                       option A; the current runs do not use it

WHAT EXISTS AND WHAT DOES NOT
    The 1984/2004 convergence figure reads convergence_history.json, the
    per-pass RMSE of the iterated zone calibration. The 1996/2010 pair has no
    such file: its zone field (2-calibrate/1996_2010__2010_2024/) is a single
    pass. So no zone-iteration convergence figure can be drawn for it, and
    none is invented. What the current runs DO carry is the edge-only preset,
    solved on the adopted model (2026-09-28, secant steps plus direct probes
    for 2010 GIS 1) and, for 2010 GIS 90, re-solved after the dune-cap fix;
    be_end_solve.png draws both records (end-domain-boundaries/
    2026-09-28-ends-resolved-adopted/ and 2026-09-28-ends-resolved-dunecap/).
```

Notes that were in the code:

```text
The ends the runs carry (2026-09-29): the adopted-model solve for 1996-2010 and
2010 GIS 1, and the re-solve of 2010 GIS 90 after the dune-cap fix. Until
2026-09-29 this figure drew the option A solve of 2026-09-27
(end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/), which no run
carries any more. Since 2026-09-29 the runs carry the split12 re-solve (storms
v3_split12_trim24, GIS 1 moved +0.04 in each window), drawn after the steps above.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_gis1_target_2010()`**

```text
The 2010-2024 target at GIS 1: the raw domain mean of the CoastSat LRR
(the end values are solved against it; no LOWESS reaches GIS 1).
```

**`_direct_probes()`**

```text
The adopted solve's direct GIS 1 probes (2010, GIS 90 held at the secant
value): imposed rate and residual, read from the probe runs themselves.
```

</details>
