# cascade_pipeline - the library the runner imports

Not a stage and not a script: this is the shared code behind
`scripts/hatteras_ms/HAT_hindcast_1984_2024.py` and everything that reads its
output.

| Module | Holds |
|---|---|
| `hindcast.py` | building and running CASCADE, the run name, the shoreline target |
| `run_registry.py` | where a run lives, and the run index. **Address runs through this**, never by joining paths |
| `domains.py` | the padded and GIS domain geometry, and the conversions between them |
| `roadway.py`, `nourishment.py` | the management forcing a period carries |
| `coastsat_lowess.py` | the observed rate series and its smoothing |
| `annotations.py`, `plotting/` | the geography layer and the figure types |
| `reports.py` | the blocks the runner prints, so a run log states how it was driven |

Changing anything here changes every run made afterwards. The run index records
enough per run - preset values, digest, topography version, git commit - to
tell runs made before a change from runs made after.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### __init__.py

cascade_pipeline: post-run analysis and figures for CASCADE hindcasts.

From the script's original header:

```text
cascade_pipeline: post-run analysis and figures for CASCADE hindcasts.

Geometry, shoreline extraction, CoastSat-style transect/LOWESS processing,
and figure/GIF rendering for a completed CASCADE run. Consumes a run's
output; it doesn't drive the simulation itself (that stays in your own run
script). Ships with no site content -- see e.g. hatteras_site_config.py for
how one study site (Hatteras Island) supplies its own domain geometry and
place names on top of this package.

Import from submodules explicitly rather than from the package root, e.g.:

    from cascade_pipeline.domains import DomainGeometry
    from cascade_pipeline.run_info import RunInfo
    from cascade_pipeline.shoreline import build_shoreline_matrix, compute_change_rate, compute_lrr
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, LowessConfig, build_coastsat_series
    from cascade_pipeline.plotting.shoreline_gif import GifConfig, make_all_shoreline_gifs
    from cascade_pipeline.plotting.rate_comparison import plot_rate_comparison, plot_annotated_rate_comparison
```

### annotations.py

Geographic reference annotations shared by every cascade_pipeline figure.

From the script's original header:

```text
Geographic reference annotations shared by every cascade_pipeline figure.

Community spans, village centers, piers, groins, and shoal-influence zones
are drawn identically on the shoreline GIF and both rate-comparison
figures, so the styling lives in one place instead of three.

Ships with NO site content: AnnotationConfig()'s dict fields default to
empty and its text fields default to generic placeholders. Build a
populated instance for your own site (see e.g. hatteras_site_config.py)
and pass it explicitly -- a reusable library shouldn't silently draw
somebody else's place names.
```

<details><summary>Function notes (the original docstrings)</summary>

**`AnnotationConfig()`**

```text
Geographic annotation layer + site labels, keyed by GIS domain ID.

Attributes:
    town_spans: {label: (gis_lo, gis_hi)} shaded community zones.
    village_lines: {label: gis_id} dashed village-center lines.
    piers: {label: (gis_id, label_y)} pier lines; label_y is an axes
        fraction (0=bottom, 1=top) for the rotated label.
    groins: {label: gis_id} groin lines; gis_id may be fractional
        (e.g. 5.5 = the boundary between domains 5 and 6).
    shoal_zones: {label: (gis_lo, gis_hi)} shoal/shoreline-position
        influence zones -- any number of them, not a fixed pair.
    pier_label_y: Default label_y for piers (informational; the actual
        y used is whatever is stored per-pier in `piers`).
    groin_label_y: Label y (axes fraction) for groin lines.
    groin_label_side: "center" puts the label ON the line, top-anchored
        at groin_label_y; "left" puts it beside the line on its south
        side, bottom-anchored at groin_label_y. The Buxton groin sits in
        the densest data on the island (D1-10), and a label on the line
        covered the southern CoastSat curve (2026-09-29).
    label_accretion_y: Fixed axes-fraction y for the "Accretion" side
        label on the annotated figure, or None to auto-place at the
        midpoint between the zero line and the top.
    label_erosion_y: Same as label_accretion_y, for "Erosion" (midpoint
        between the bottom and the zero line).
    region_name: Site name used in figure titles, e.g. "Hatteras Island".
    low_end_label: Text for the low-domain-ID end of the compass
        annotation, e.g. "S" or "S | South Landmark".
    high_end_label: Text for the high-domain-ID end, e.g. "N" or
        "North Landmark | N".
    obs_source_name: Name of the observational dataset shown against
        the model, e.g. "CoastSat".
    color_*: Colors for each annotation layer.
    model_color: Line color for the modeled shoreline-change-rate curve.
```

**`add_geographic_annotations()`**

```text
Add every geographic reference annotation to an axis (bottom -> top).

Draw order: shoal zones, community spans, village center lines, pier
lines, groin lines. Label y-positions use blended axes-fraction coords,
so they hold position regardless of the y-data range. X-axis must be in
GIS domain IDs. Any category left empty in `config` simply draws
nothing -- there's no minimum set of annotations required.

Args:
    ax: Matplotlib Axes to annotate.
    config: AnnotationConfig; defaults to an empty (no-op) layout.
    label: Draw the name beside each band and line. Set False on a small
        multiple: the bands stay legible at any panel width but their
        labels do not, and on a six-panel grid the same eight names would
        be repeated six times (2026-09-17). Name them in the caption
        instead.
```

**`annotation_legend_handles()`**

```text
Proxy legend artists matching add_geographic_annotations' layers.

Only returns a handle for categories that are actually populated in
`config`, so an empty site config produces an empty legend rather than
five swatches for things that weren't drawn.

Args:
    config: AnnotationConfig; should match what was passed to
        add_geographic_annotations, or the legend will describe
        something different from what's drawn.

Returns:
    List of matplotlib legend handles (Patch/Line2D).
```

</details>

### coastsat_lowess.py

CoastSat transect loading and LOWESS smoothing for the rate-comparison figures.

From the script's original header:

```text
CoastSat transect loading and LOWESS smoothing for the rate-comparison figures.

Transect-level LRR (linear regression rate) values are loaded from
transect_lrr_full.csv, LOWESS-smoothed at transect resolution (physical
along-coast distance as x), then aggregated to GIS-domain resolution for
comparison against the CASCADE model. This module only computes -- nothing
here touches matplotlib; see cascade_pipeline.plotting.rate_comparison for the
figures that consume its output.
```

Notes that were in the code:

```text
From the geometry's first domain, not from 1: rate_comparison maps it
back with along_m / spacing + first_gis_id, and a reach starting at
GIS 0 (an extended geometry) would otherwise plot one domain off.
```

<details><summary>Function notes (the original docstrings)</summary>

**`CoastSatDataset()`**

```text
One CoastSat transect CSV to load and smooth.

Attributes:
    label: Legend label, e.g. "CoastSat LRR (1984-2004)".
    period_start: Calendar year this dataset's period begins; compared
        against the active run's start year to decide solid vs faded
        styling downstream (build_coastsat_series' active flag).
    csv_path: Path to transect_lrr_full.csv.
    domain_col: Column holding the GIS domain ID.
    rate_col: Column holding the LRR value (m/yr).
    transect_id_col: Column holding a per-transect ID, used to keep
        each domain's transects in a stable order.
```

**`LowessConfig()`**

```text
LOWESS smoothing settings shared by every CoastSat dataset.

Attributes:
    window_domains: One or two window widths, in domain units
        (1 domain ~= 500 m). The largest is treated as the primary
        reference window wherever a single window is needed.
        Hannah, 2026-09-28: 7 domains, the research group's range,
        for every smoothed observation in the project. From 2026-09-10
        to 09-28 the default and every skill target were 10 alone.
    skip_southern_domains: Domains 1..N shown as raw per-domain means
        instead of LOWESS-smoothed -- boundary effects near Oregon Inlet
        dominate this zone and smoothing can obscure the sharp
        gradient there. 0 disables the splice (LOWESS used everywhere).
```

**`estimate_transect_spacing()`**

```text
Median spacing between consecutive transects, in meters.

Args:
    along_coast_m: Along-coast distance for each transect.

Returns:
    Median of the positive consecutive differences, or 50.0 if there
    are fewer than two distinct positions.
```

**`load_transect_data()`**

```text
Load individual transect LRR values for one CoastSatDataset.

Along-coast distance is derived by spreading each domain's transects
evenly across its domain_spacing_m band.

Args:
    dataset: A CoastSatDataset.
    domains: DomainGeometry; first_gis_id/last_gis_id/domain_spacing_m
        are used to filter and space transects.

Returns:
    (domain_ids, lrr_values, along_coast_m): int/float arrays, or
    (None, None, None) if the CSV doesn't exist.
```

**`lowess_transect_values()`**

```text
The LOWESS itself: one smoothed value per TRANSECT, before any domain
averaging.

Factored out of lowess_smooth_transect_to_domains (2026-09-22) so a figure
can draw the smoother at the resolution it is actually fitted at without a
second copy of the frac rule. That function now calls this and aggregates
the result, so there is one implementation and the drawn curve cannot
drift from the graded one.

Args:
    along_coast_m, lrr: Per-transect arrays from load_transect_data.
    window_domains: LOWESS window width, in domain units.
    domains: DomainGeometry; only domain_spacing_m is used.

Returns:
    (smoothed, frac): smoothed is a per-transect array, NaN where lrr is,
    or None when fewer than 5 transects have a rate at all.
```

**`lowess_smooth_transect_to_domains()`**

```text
Apply LOWESS at transect resolution, then aggregate to domain resolution.

Args:
    along_coast_m, lrr, domain_ids: Output of load_transect_data.
    window_domains: LOWESS window width, in domain units.
    domains: DomainGeometry; only domain_spacing_m is used.

Returns:
    (gis_x, smoothed, frac): GIS domain IDs with at least one transect,
    the domain-averaged smoothed LRR (m/yr), and the LOWESS frac used
    (for logging). (None, None, frac) if fewer than 5 valid transects.
```

**`compute_domain_means()`**

```text
Mean LRR per GIS domain within [gis_min, gis_max].

Used to substitute raw per-domain averages for LOWESS smoothing in the
southernmost domains, where boundary effects dominate.

Args:
    domain_ids, lrr_values: Per-transect arrays (e.g. from
        load_transect_data).
    gis_min, gis_max: Inclusive GIS domain ID range.

Returns:
    (gis_x, means): GIS domain IDs with at least one transect, mean
    LRR per domain.
```

**`splice_lowess_with_raw_south()`**

```text
Return (plot_x, plot_y) for one LOWESS window, LOWESS line starting at skip_n+1.

Domains 1..skip_n are omitted from the returned line -- they're shown as
raw scatter only, no smoothed line, since boundary effects near Oregon
Inlet dominate that zone.

Args:
    win_gis_x, win_smoothed: Output of lowess_smooth_transect_to_domains.
    transect_domain_ids, transect_lrr_values: Accepted for interface
        stability with callers that pass per-window transect data; not
        read by the current splice.
    skip_n: Domains 1..skip_n excluded from the LOWESS line.
    is_widest_window: Accepted for interface stability; not read by the
        current splice.

Returns:
    (plot_x, plot_y): arrays to plot, LOWESS domains > skip_n only.
```

**`build_coastsat_series()`**

```text
Load every CoastSat dataset and LOWESS-smooth it at each configured window.

Replaces the per-dataset loop previously inlined in main(): loads
transects, applies lowess_smooth_transect_to_domains at each window in
lowess_config.window_domains, and tags each dataset active/reference by
comparing its period_start to active_period_start (the run's START_YEAR).

Args:
    datasets: Sequence of CoastSatDataset.
    active_period_start: The run's START_YEAR; datasets with a matching
        period_start are drawn solid/full-opacity downstream.
    lowess_config: LowessConfig.
    domains: DomainGeometry.

Returns:
    List of dicts, one per dataset that loaded successfully:
        label, period_start, active (bool), transect_domains,
        transect_rates, transect_along_coast,
        windows: list of dicts (window, gis_x, smoothed, frac).
```

**`scale_coastsat_series()`**

```text
A copy of build_coastsat_series output with every rate times `factor`.

Turns an LRR series (m/yr) into a distance (m): LRR x span years. EXACT
for the LOWESS curves as well as the transects, because statsmodels'
lowess -- robustness iterations included -- is equivariant under a
positive scale (the robustness weights see residual / median |residual|,
which a constant cancels), so smoothing the scaled transects would return
`factor` x the stored curve.

Args:
    cs_series: build_coastsat_series output.
    factor: Multiplier, e.g. 14 for a 14-year window.
    label: Replace every series' label, or None to keep them.
    active: Force every series' `active` flag, or None to keep them. Set
        True when the series is a reference fitted on a different window
        (the 1996-2024 LRR drawn against a 1996-2010 run).
```

**`spliced_lowess_series()`**

```text
Per-domain series: a LOWESS of `values` at transect resolution north of
GIS `skip`, the raw domain means at or below it.

The two steps hindcast.build_target_table applies to the scoring target,
factored out (2026-09-21) so an analysis can build the same series at a
window other than TARGET_WINDOW without going through a CoastSatDataset.
build_target_table itself is untouched -- it is the production scoring
path and still reads its window off a built CoastSatDataset.

Args:
    domain_ids, along_coast_m, values: per-transect arrays.
    window: LOWESS window width in domain units. 0 means no smoothing at
        all: raw domain means everywhere, the unsmoothed comparison.
    skip: domains <= skip keep their raw means. Default the shared
        LowessConfig's skip_southern_domains.
    domains: DomainGeometry.

Returns:
    (Series indexed first_gis_id..last_gis_id, the lowess frac used, or
    nan when window is 0).
```

</details>

### domains.py

Generic CASCADE grid geometry: real domains padded by buffers on each end.

From the script's original header:

```text
Generic CASCADE grid geometry: real domains padded by buffers on each end.

A CASCADE/BRIE run pads its GIS-numbered real domains with buffer domains
on each side so alongshore diffusion doesn't see an artificial edge at the
real-domain boundary. DomainGeometry's field defaults happen to match the
Hatteras Island hindcast (90 real domains, 15-domain buffers, 500 m
spacing, starting at GIS ID 1) -- override them for a different site or
grid resolution; nothing else in this module assumes Hatteras' numbers.
```

<details><summary>Function notes (the original docstrings)</summary>

**`DomainGeometry()`**

```text
Describes the padded CASCADE domain array.

Attributes:
    num_real_domains: Count of real, GIS-numbered domains (defaults to
        90, matching the Hatteras Island hindcast).
    num_buffer_domains: Count of buffer domains padding each end
        (defaults to 15).
    first_gis_id: GIS ID of the first real domain (defaults to 1).
    domain_spacing_m: Alongshore width of one domain, in meters
        (defaults to 500).
```

**`gis_to_pad()`**

```text
Convert a GIS domain ID (e.g. 1-90) to its padded array index.

Args:
    gis_id: GIS domain ID, or a numpy array/sequence of them.

Returns:
    The padded array index (or array of indices).
```

**`pad_to_gis()`**

```text
Convert a padded array index back to its GIS domain ID.

The inverse of `gis_to_pad`. Buffer domains have no GIS ID, so an
index outside the real span returns a number below first_gis_id or
above last_gis_id -- check the range yourself if it matters.

Args:
    pad_index: Padded array index, or a numpy array/sequence of them.

Returns:
    The GIS domain ID (or array of IDs).
```

</details>

### hindcast.py

Shared machinery for the Hatteras hindcast: loaders, diagnostics, build and run.

From the script's original header:

```text
Shared machinery for the Hatteras hindcast: loaders, diagnostics, runner.

WHY THIS MODULE EXISTS
    `HAT_hindcast_1984_2024.ipynb`, its headless mirror
    `HAT_hindcast_1984_2024.py`, and `HAT_groin_sweep_worker.py` each carried
    their own copy of these functions -- 667 lines duplicated character for
    character between the notebook and the .py, and a third copy of
    `build_cascade` in the sweep worker that the worker's own comment warned
    "can drift". It had already drifted in the docstrings.

    They live here now, defined once. The notebook and the .py keep everything
    that describes THIS study -- the paths, the switches, the reports, the
    figures -- and import the machinery that is the same either way.

WHAT IS NOT HERE
    Nothing that decides what is simulated. Every run-selecting value still
    comes from `HAT_hindcast_config` (which reads `hat_run.yaml`, overridden by
    the environment) and every path still comes from section 2 of the calling
    file, so a run remains fully described by the settings file plus the file
    you can read top to bottom. The functions below take those
    values as arguments rather than reaching for module globals, which is the
    only change made to any body while moving it.

THE SYNC RULE
    Editing a function here changes the notebook, the .py, and the sweep at
    once. That is the point. Verify with a full hindcast run, not by reading.
```

Notes that were in the code:

```text
True: sandbox copy with the pre-AST groin hook. False: real package, hook
folded in. Resolved here so the notebook, the .py and the sweep worker cannot
disagree about which Cascade they built.

Read from the environment rather than from `HAT_hindcast_config`, which is
the deliberate direction of the dependency: this package is shared, and the
Hatteras settings file is not something it should know about. The runner
reads `use_sandbox_cascade` from hat_run.yaml and exports it as
CASCADE_USE_SANDBOX in its section 1, BEFORE this module is imported -- the
choice selects an import, so it cannot be made after the fact.
```

```text
scripts/ is the parent of this package, so the resolver is importable
without a path hack. It owns the dune-topo directory AND the array names.
```

```text
--- unit and file-format constants ------------------------------------------
Fixed by Barrier3D's input contract and by the extractor's RUN_MANIFEST, not
by the scenario, so they are the same for every run and every period.
```

```text
BRIE's (1.2 sin^2 - cos^2) changes sign here; above it the
shoreline is anti-diffusive.
```

```text
A zero peak means the two runs are identical, so the threshold is
also zero and every "abs(value) < 0" test is False -- the walk would
run to the end of the grid and report the whole island as fillet.
```

```text
An extended geometry (2026-09-16) reaches beyond the surveyed raw; its
domains are in raw_offsets/ext/<vintage>_duneline_offset_raw_ext.csv
(duneline_to_raw_offsets.py --extension), same columns. Read only when
the geometry needs it, so a base run never touches the file.
```

```text
A sensitivity cell differs from the matrix run beside it in ONE value, and
nothing in SCENARIO_SWITCHES sees that value -- the switches describe which
management modules were built, not what they were forced with. So without a
token of its own, an Hs = 1.2 cell derives the matrix run's exact name, lands
in its directory, and replaces its row in run_index.csv. That is the failure
output/calibration/groin/README.md records for the rig sweep, and it is silent.

The token is emitted ONLY when the value is off its calibration default, so
every run that predates this file is named exactly as it was. A run at the
defaults has no token, which is correct: it IS the matrix run.
```

```text
Field name -> the abbreviation that appears in the run name. Abbreviated
rather than spelled out because all four can move at once and
"..._wave_angle_high_fraction0p3" is longer than the rest of the name.
```

```text
The split is what lets the caller hold a built-but-unstepped Cascade. BRIE's
diffusivity and the groin's fillet prediction are only meaningful as initial
conditions, and a prediction printed after the run is not one.

`run_years` is TRANSITIONS, not states. Barrier3D seeds _x_s_TS = [x_s] at
init and appends one entry per update, so N updates produce N+1 annual states
spanning N years. The original signature took `nt` and looped `range(nt - 1)`
with `nt = END_YEAR - START_YEAR`, which ran 19 updates for a 20-year period
while dividing by 20 -- every rate came out low by 19/20, and storm years 20
and 21 were never applied. Here time_step_count = run_years + 1 and the loop
runs exactly run_years updates.

Nourishment goes to cascade.nourishment_volume, which is where CASCADE reads
it. Writing it onto the manager instead hits the attribute CASCADE overwrites
one line before the manager reads it, so the fill spends the init default.
```

```text
run_years transitions need run_years + 1 states. TMAX is set from
this (brie_coupler.py:117), and the loop runs exactly run_years
updates, so the last write lands on the final valid index.
```

```text
A STANDARD RELOCATION TARGET IS NOW AN ARGUMENT, not a fix-up.

`road_setback` used to do two jobs in CASCADE: it placed the road at
t = 0, AND cascade_groin.py re-read it every year as the distance a
relocated road is rebuilt at. There was no separate parameter, so this
function used to overwrite `cascade._road_setback` after construction --
safe only because the constructor had already consumed it, and only for
as long as that stayed true.

Since 2026-09-14 the model takes `road_relocation_setback` directly and
the yearly update reads that instead, so the measured array keeps its
measured meaning and nothing here reaches into the object afterwards.
Passing None gives the old behaviour: every domain relocates to its own
measured offset.

WHY A STANDARD AT ALL. The measured setback is observed 1984/2004
geometry, not a design standard: 30 distinct values from 0 to 430 m
across the 55 road domains, exactly one of which is 30 m. Using it as the
relocation target makes a domain's rebuild rule an accident of where the
road happened to sit. At GIS 85 and 86 the measured value is 0 m, so a
relocation returns the road to the dune line with no clearance and the
next 10 m of retreat re-fires it -- 7 relocations against 7.3 cells of
retreat at GIS 85, 6 against 6.0 at GIS 86, one per cell.

WHAT IT DOES NOT CHANGE. Every domain's FIRST relocation still fires in
exactly the year it fires today, because that depends only on the initial
setback. What changes is where the road lands, and so every trigger after
the first.

ONE SIDE EFFECT WORTH KNOWING. `road_relocation_checks` refuses to
relocate when `setback + 2 * road_width > average_barrier_width`, and
sets relocation_break, which abandons the road. With measured targets up
to 430 m that guard can fire; with a 20 m standard it effectively cannot.
So a standard does not only move roads, it makes relocation POSSIBLE in
narrow domains where the measured target would have been refused.

WHAT IT STILL DOES NOT DECOUPLE. A prescribed historical relocation adds
its measured DISPLACEMENT to the road's live setback, so the clearance a
road was last rebuilt at still affects where a later historical event
lands. That is correct -- a historical relocation moved the road from
wherever it then was -- and the alternative, an absolute setback in the
1984 frame, would count the dune migration between 1984 and the event
twice. The sensitivity is in the physics, not in the plumbing.
```

```text
One emitter for everything the loop says, so the display mechanism is
the caller's choice and this function does not care which it is.
```

```text
--- historical beach nourishment -----------------------------------
apply_to_cascade rewrites BOTH nourish_now and nourishment_volume in
full every year, so nothing carries over. It writes the volume to
cascade.nourishment_volume, which stock Cascade.update() copies into
each BeachDuneManager before calling it -- see section 6.
```

```text
The .npz model pickle is ~160 MB and is what lets a figure be re-derived
without re-running, so it is written by default. Everything downstream
(rate CSV, shoreline matrix, metadata) is written regardless, so a run
skipped here is still a complete row in run_index.csv -- just one that
cannot be re-plotted from state.
```

<details><summary>Function notes (the original docstrings)</summary>

**`build_domain_file_paths()`**

```text
Builds one elevation and one dune file path per padded domain.

Delegates to hat_topo_version.domain_arrays(), which owns BOTH the
directory and the filename. This function used to spell the names itself -
f"domain_{gis_id}_topography_{init_year}.npy" - and that body was copied
verbatim into the groin sweep worker and the notebook, each with its own
`init_year`. On 2026-08-26 the year was dropped from the names entirely and
the three copies collapsed into one call; see the note at the top of
hat_topo_version.py for why a per-period tag was tried and reverted.

Args:
    geometry: DomainGeometry describing the padded domain array.
    product: topography product, e.g. "1984-start". None resolves the
        resolver's DEFAULT_PRODUCT.
    override: pin a dune-topo version instead of resolving it.

Returns:
    An (elevation_paths, dune_paths) tuple of string lists. Each list is
    geometry.total_domains long and index-aligned with the padded array.
```

**`load_barrier3d_contract()`**

```text
Reads the unit-relevant Barrier3D parameters from the CASCADE YAML.

Mirrors the conversions in barrier3d/load_input.py so the values can be
compared against the raw arrays on disk.

Args:
    parameter_path: Path to the CASCADE parameter YAML.

Returns:
    A dict with barrier_length_cells, mhw_dam and berm_el_dam.
```

**`check_domain_units()`**

```text
Checks one domain's raw arrays against the Barrier3D input contract.

Args:
    elevation_dam: Raw elevation array as stored on disk.
    dune_dam: Raw dune array as stored on disk.
    contract: Mapping from load_barrier3d_contract.

Returns:
    A dict of check name -> bool, True when the array satisfies the check.
```

**`load_island_offset_dam()`**

```text
Loads a padded BRIE island-offset file and converts it to decameters.

Args:
    offset_path: Path to a single-column padded offset CSV, in meters.
    geometry: DomainGeometry the file must match in length.

Returns:
    A 1-D array of offsets in decameters, one per padded domain, ready to
    pass to Cascade() as shoreline_offset.

Raises:
    ValueError: If the file is not one value per padded domain.
```

**`close_offset_ring()`**

```text
Buffer values carrying the shoreline from real_m[-1] back to real_m[0].

BRIE's alongshore domain is PERIODIC -- `x_s[np.r_[1:ny, 0]]` makes the
last padded domain a neighbour of the first -- so the padded offset has to
come back to where it started. Hatteras does not: it runs 5.4 km across
the shore from Cape Point to Oregon Inlet and keeps going. The buffers are
the invented coast that closes the loop.

A cubic Hermite is used, with end tangents matched to the island's own
local slopes, so there is no kink where buffer meets real domain. The path
runs real[-1] -> right buffer -> (wrap) -> left buffer -> real[0], which is
2 * buffers + 1 steps.

Why it matters that this be smooth: BRIE's diffusivity term
`(1.2 sin^2 - cos^2)` changes SIGN near 42 degrees. Below that a shoreline
smooths itself; above it, bumps GROW (high-angle instability). The original
linear-bridge padding leaves a 1154 m step at the seam, which is 66.6
degrees once the offset is in its correct unit -- unstable. It went
unnoticed because the decameter bug divided it to a harmless 13 degrees.

Args:
    real_m: The real-domain offsets, in meters.
    buffers: Buffer domains per side.

Returns:
    A (right_buffer, left_buffer) tuple of arrays, each `buffers` long.
```

**`pad_offset_ring()`**

```text
The padded island offset: buffers, real domains, buffers, closed smoothly.

The one definition of the padding. island_offset_hybrid.py writes the
offset files with it (every v1 since 2026-09-24), and build_island_offset
re-closes with it, so a file built by it is exactly what metres mode hands
Cascade.

Args:
    real_m: The real-domain offsets, in metres, south to north.
    buffers: Buffer domains per side.

Returns:
    A 1-D array of len(real_m) + 2 * buffers values, in metres.
```

**`build_island_offset()`**

```text
Builds the `shoreline_offset` array Cascade is handed.

UNITS. `Cascade(shoreline_offset=...)` must be in METRES.
`brie_coupler.offset_shoreline` adds the values straight onto
`brie.x_t` / `brie.x_s` with no conversion, and those are metres -- see
`brie_coupler.py:390` ("convert from dam to meters") and line 344
(`barrier3d.x_s = brie.x_s / 10`). The measurement file is already metres,
so it needs no conversion at all.

Args:
    offset_path: Padded offset CSV, one value per padded domain, meters.
    geometry: DomainGeometry the file must match in length.
    mode: Which variant to build.
        "asrun"     - meters / 10, reproducing the historical unit error
                      exactly, buffers included -- so only with the
                      build a run was made from (v1 for every run before
                      2026-09-24, now superseded_20260924_pre-metres/v1; the
                      current v1 closes the buffer differently).
        "metres"    - the measurement as-is, with the ring re-closed.
                      Carries the island's full planform including its
                      ~7 degree lean. The default since 2026-09-24. A
                      file built since then (the current v1) holds this
                      closure, so it comes back unchanged; an older
                      file's slope-and-bridge buffers are replaced.
        "detrended" - the measurement with its linear trend removed, ring
                      re-closed. Keeps all 2052 m of real curvature and
                      drops the lean, which is what lets the periodic
                      domain close without a steep buffer.

Returns:
    A 1-D array, one value per padded domain, ready to pass to Cascade.

Raises:
    ValueError: If the file length does not match the geometry, or the
        mode is unknown.
```

**`island_offset_tilts()`**

```text
Local shoreline tilts BRIE will read from an offset array, in degrees.

Reproduces `brie.py:820`, `theta = atan2(diff(x_s), dy)`, including the
periodic wrap, so a run can report the angles it is about to impose rather
than leaving them to be discovered in the output.

Args:
    offset: The padded array passed to Cascade, meters.
    geometry: DomainGeometry.

Returns:
    A dict with island_max_deg, buffer_max_deg, seam_deg and unstable.
```

**`load_storm_series()`**

```text
Loads a Barrier3D storm series into a DataFrame.

Args:
    storm_path: Path to the storm .npy file.

Returns:
    A DataFrame with the raw columns plus Rhigh_m and Rlow_m in meters.

Raises:
    ValueError: If the file does not have the expected column count.
```

**`build_background_erosion()`**

```text
Expands sparse per-GIS background-erosion rates onto the padded array.

Args:
    be_rates: Mapping of GIS domain ID to rate in m/yr. Domains absent
        from the mapping get 0.0.
    geometry: DomainGeometry describing the padded array.

Returns:
    A list of geometry.total_domains rates, ready to pass to Cascade().

Raises:
    ValueError: If a GIS ID falls outside the padded array.
```

**`groin_trapping_schedule()`**

```text
Effective trapping rate for every model year of a period.

Reads the callback's deterioration curve without running it.
`_effective_trapping_rate` is a pure function of the year -- it does not
touch `_call_count` -- so querying it here leaves the callback's year
counter at zero for the actual run.

Args:
    callback: A GroinCallback.
    start_year: First model year.
    end_year: Last model year, exclusive.

Returns:
    A (years, M_eff) tuple of 1-D arrays.
```

**`implied_interception_m3_yr()`**

```text
Sediment volume the dipole transfers across the groin each year.

A shoreline displacement of `m_per_yr` over one domain implies a volume of
`m_per_yr * dy * (d_sf + h_b)`. The dipole moves that from the downdrift
side to the updrift side; the transfer is volume-neutral overall.

Args:
    m_per_yr: Trapping rate M, in meters per year.
    profile_height_m: Active profile height (d_sf + h_b), in meters.
    geometry: DomainGeometry supplying the alongshore domain width.

Returns:
    The implied transfer in cubic meters per year.
```

**`measure_groin_extent()`**

```text
Alongshore extent of the groin's effect, from a paired baseline run.

The effect is this run's final shoreline minus the no-groin baseline's,
sign-flipped so positive means seaward of the baseline. Extent is the
contiguous run of domains, outward from the groin, where the effect holds
at or above `threshold_frac` of its own peak magnitude.

Args:
    shoreline_m: [time x domain] matrix from this run.
    baseline_m: [time x domain] matrix from the no-groin run.
    geometry: DomainGeometry describing the padded array.
    updrift_gis, downdrift_gis: The groin's flanking domains.
    threshold_frac: Fraction of the peak effect defining the edge.

Returns:
    A dict with effect (padded array, m), peak_m, threshold_m, and the
    updrift/downdrift extents in domains and meters.
```

**`build_target_table()`**

```text
Per-domain target rate, labelled with where each value came from.

The target is not one curve: GIS 1..skip_southern_domains are raw
per-domain means (LOWESS is suppressed there), and the rest is the
LOWESS-smoothed reference window. Every row records which.

Args:
    cs: One entry from build_coastsat_series.
    lowess_config: LowessConfig used to build it.
    geometry: DomainGeometry describing the real-domain span.
    window: Reference window width, in domains.

Returns:
    A DataFrame with gis_domain, target_lrr_m_yr, source, n_transects.

Raises:
    ValueError: If `window` was not computed for this dataset.
```

**`groin_differential()`**

```text
Observed updrift-minus-downdrift mean LRR for one CoastSat period.

The observational check on the dipole section 7 imposes: a positive
differential means the updrift domain is retreating more slowly (or
advancing faster) than the downdrift one, which is the signature a
functioning groin should leave.

Args:
    cs: One entry from build_coastsat_series.
    updrift_gis: GIS domain on the updrift side.
    downdrift_gis: GIS domain on the downdrift side.

Returns:
    A dict of label, updrift, downdrift, differential, and per-domain
    transect counts. Rates are NaN where a domain has no transects.
```

**`load_absolute_dune_distance()`**

```text
Mean distance from the offshore datum to the dune line, per GIS domain.

The absolute quantity the padded offset files are built from, BEFORE
island_offset_hybrid.py subtracts that year's own minimum. Absolute
distances share a fixed datum across years, so differencing two of them is
a real shoreline change; differencing the padded files is not.

One row per transect (`LineID`) is kept before averaging, mirroring
island_offset_hybrid.py:108-113, so a transect sampled at many points does
not outweigh one sampled at few.

Args:
    year: Period year (a start or an end). The file read is
        <raw_dir>/<vintage>_duneline_offset_raw.csv, where the vintage
        is hat_topo_version.DUNE_LINE_FOR_YEAR[year] -- the imagery
        year of the digitised line that stands for this period year
        (1996 reads the 1997 line). Since 2026-09-15; before that a
        period year with no line of its own read a copy filed under
        its name.
    geometry: DomainGeometry supplying the real-domain GIS range.
    raw_dir: Directory holding the raw per-vintage CSVs.
    columns: Mapping with "domain", "distance" and "transect" keys.

Returns:
    1-D array of length geometry.num_real_domains, in meters, increasing
    LANDWARD, indexed by GIS domain. NaN where a domain has no transects.

Raises:
    FileNotFoundError: If that year has no raw file.
    KeyError: If an expected column is missing.
```

**`build_shoreline_target()`**

```text
Surveyed end-year shoreline position, in the model's own x_s frame.

Takes the model's year-0 position and adds the OBSERVED start->end change,
so the two share a base and the gap between them is the model's misfit.
Returns positions rather than a change because that is what
make_shoreline_gif's `target_m` expects.

Args:
    model_year0_m: shoreline_m[0], the run's year-0 position (raw x_s_TS
        convention, meters, increasing landward), one value per PADDED
        domain.
    start_year: Run start year; must have a raw offset file.
    end_year: Run end year; returns None if it has no raw offset file.
    geometry: DomainGeometry.

Returns:
    (target_m, observed_change_m):
        target_m: padded array of target positions, NaN in the buffers.
        observed_change_m: real-domain change, + = landward.
    Both None if the end year was never surveyed.
```

**`scenario_run_name()`**

```text
Run name for a variant of the current scenario.

Rebuilds the name from the same switch tokens section 7.5 used, with named
switches overridden. Built from the token list rather than by editing the
name string, so flipping "groin" cannot accidentally match inside "nogroin".

Args:
    switches: SCENARIO_SWITCHES, as (label, value, token) triples.
    stem: RUN_NAME_STEM, the period prefix.
    **overrides: label -> replacement token. A token of None or "" drops
        that switch from the name, matching how 7.5 omits default states.

Returns:
    The run name for that variant.

Raises:
    KeyError: If an override names a switch that does not exist -- a typo
        would otherwise silently produce the unmodified name.
```

**`_number_token()`**

```text
A number, spelled so it can live in a directory name.

"." becomes "p" and a leading "-" becomes "m", because a name is split on
"_" and read by eye: 1.2 -> "1p2", 8.0 -> "8". %g drops the trailing zero
so 8.0 and 8 produce one token, not two names for one model.
```

**`wave_climate_token()`**

```text
The run-name token for a wave climate, or None if it is the default one.

Args:
    values: field name -> value this run is using, for the four fields in
        _WAVE_TOKEN_FIELDS.
    defaults: the same keys -> the calibration default, from
        HAT_hindcast_config.field_default.

Returns:
    e.g. "waveHs1p2", or "waveHs1p2Tp10" if two moved, or None when every
    value is at its default.
```

**`relocation_setback_token()`**

```text
The run-name token for the relocation target, or None if default.

Spelled `rset`, not `reloc`: `reloc` is already the token for whether the
historical 1989/1999 events fire at all, and the two are independent -- a
run can have the events off and still relocate emergently.

Args:
    value: relocation_setback_m for this run; None means "each domain
        relocates to its own measured offset".
    default: the calibration value, from field_default.

Returns:
    e.g. "rset40", or "rsetmeasured" for the None case, or None when the
    value is the default.
```

**`build_cascade()`**

```text
Constructs a Cascade and attaches the groin, without stepping it.

Separate from run_cascade_simulation so section 11 can inspect the model
before any time step runs: BRIE's diffusivity and the shoreface depth are
only meaningful as initial conditions, and the fillet prediction they feed
is only a prediction if it is printed before the run.

Args:
    run_years: Annual transitions the run will simulate. The model is built
        with time_step_count = run_years + 1, giving run_years + 1 annual
        states. See the section 10 comment on why this is not run_years.
    name: Run name, passed to Cascade.
    storm_file, elevation_file, dune_file: Barrier3D input paths.
    alongshore_section_count: Padded domain count.
    num_cores: Cores for the parallel Barrier3D step.
    rmin, rmax: Dune growth rate bounds, per domain.
    dune_design_elevation: Rebuild target, m MHW. roadway_manager raises
        this to berm + 1.0 m on the first step if it is lower.
    dune_minimum_elevation: Rebuild trigger, m MHW. roadway_manager raises
        this to berm + 0.3 m on the first step if it is lower.
    road_ele, road_width, road_setback: Padded roadway forcing.
    overwash_filter, overwash_to_dune: beach_dune_manager forcing.
    nourishment_volume: Per-domain init volume. Every scheduled year
        overwrites it via the schedule; see section 6.
    background_erosion: Padded source/sink rates.
    roadway_management_on, beach_dune_manager_on: Per-domain module flags.
    sea_level_rise_rate, sea_level_constant: RSLR forcing.
    sandbag_management_on, sandbag_elevation: Sandbag forcing.
    enable_shoreline_offset, shoreline_offset: Island orientation, in
        METRES: brie_coupler.offset_shoreline adds it to brie.x_s with
        no conversion (build_island_offset).
    wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction:
        Wave climate. wave_height also sets BRIE's shoreface depth,
        d_sf = 8.9 * Hs (brie.py:270).
    berm_elevation, MHW: Barrier3D datums.
    data_base: Datadir handed to Cascade -- the directory the
        parameter file and the padded input files are resolved
        against.
    parameter_file: Barrier3D parameter yaml name, resolved by
        CASCADE inside data_base.
    groin_callback: A GroinCallback to attach, or None.
    relocation_setback_m: Standard distance behind the dune line, in
        metres, that a relocated roadway is rebuilt at. None leaves every
        domain relocating to its own measured setback, which is CASCADE's
        built-in behaviour. See the block below for why the two are
        separable and why they were not.

Returns:
    The constructed Cascade, before any update().
```

**`run_cascade_simulation()`**

```text
Steps a built Cascade through its period and writes the run artifacts.

Args:
    cascade: A Cascade from build_cascade, not yet stepped.
    run_years: Annual transitions to simulate; the loop runs exactly this
        many updates.
    name: Run name, used for output filenames.
    run_dir: Directory for the saved model and logs.
    start_year: Calendar year of the run's first state.
    geometry: DomainGeometry, for GIS <-> pad translation.
    alongshore_section_count: Padded domain count.
    historical_road_events: RelocationEvent / BridgeEvent sequence.
    relocations_enabled: Global toggle for relocation events.
    setback_check: {gis: measured_setback_m} reported beside relocations.
    nourishment_schedule: A NourishmentSchedule, or None for no fills.
    groin_callback: The attached GroinCallback, for diagnostics output.
    progress: Optional tqdm-like object with update() and write().
        When given, the year counter advances the bar and every event
        message is written through it, so the messages scroll above a
        live bar instead of shredding it. None falls back to a plain
        carriage-return counter, so callers without tqdm still work.
    save_model_state: Whether to write the ~160 MB model pickle.
        The caller passes RUN_CONFIG.save_model_state; everything
        else is written either way.

Returns:
    The same Cascade, after the run.
```

**`brie_r_ipl()`**

```text
BRIE's diffusion number at the freshly built model's initial state.

Reproduces `brie.py:1294`, which computes r_ipl as a local variable inside
update() and never stores it:

    r_ipl = coast_diff[clip(round(90 - theta))] * dt / 2 / dy**2

`_coast_diff` is the wave-climate-averaged shoreline diffusivity, built once
in BRIE's __init__ and deleted by its finalize(). Nothing in CASCADE calls
finalize, but the angle-dependent index changes as the shoreline evolves, so
this is only the initial-condition value.

Args:
    cascade: A built Cascade.
    theta_deg: Shoreline angle to evaluate at. 0 is shore-normal, the
        reference used for the fillet scaling.

Returns:
    The dimensionless diffusion number.
```

</details>

### nourishment.py

Beach and dune forcing for a CASCADE run: nourishment schedule and overwash filter.

From the script's original header:

```text
Beach and dune forcing for a CASCADE run: nourishment schedule, overwash filter.

Everything CASCADE's `beach_dune_manager` needs, prepared before the run and
checked against what the module actually did afterwards. Nothing here is
site-specific -- domain geometry arrives as a `DomainGeometry`, and the
Hatteras instances (project extents, volumes, community zones) live in
`hatteras_site_config`.

Four things about `BeachDuneManager` this module exists to get right:

- **`overwash_filter` is a PERCENT, not a fraction.** `filter_overwash` divides
  it by 100, and the docstring cites 40-90 % from Rogers et al. (2015):
  residential to commercial. A value of 0.4 filters 0.4 % of overwash, which
  is indistinguishable from no filtering. `BeachDuneConfig` rejects the
  fraction scale rather than letting it pass silently.

- **The per-year volume must be written to the Cascade object, not the
  manager.** `cascade.update()` copies its own `_nourishment_volume[iB3D]`
  into each manager on every time step, immediately before calling
  `BeachDuneManager.update()` (cascade.py:772). A volume assigned straight
  onto `cascade.nourishments[i]` is therefore overwritten before it is ever
  used, and every nourishment silently applies the Cascade init default
  instead. `NourishmentSchedule.apply_to_cascade` writes
  `cascade.nourishment_volume`, which survives.

- **Enabling the module is not the same as scheduling a nourishment.** With
  `beach_nourishment_module=True`, a domain runs overwash filtering, the 50 m
  community-width drowning check, and fixed-dune-line dynamics EVERY year of
  the run -- `nourish_now` only controls whether sand is added. There is no
  way to get the fill without the rest.

- **The manager's time series are offset by one.** `update_dune_domain`
  increments `barrier3d.time_index` before the managers run, so a manager
  writing at `time_index - 1` lands on index `year - start_year + 1`.
  `NourishmentSchedule.time_index` is the single place that conversion lives.
```

Notes that were in the code:

```text
Two projects overlapping in one domain-year would silently
clobber each other; sum instead, which is what placing both
volumes on the same beach means.
```

```text
Real projects run roughly 20-2000 m^3/m; outside that is a unit slip
(total m^3 not divided by domain length, or cy never converted).
```

```text
Trailing zeros are unwritten years after the road stopped being
managed, not a setback of zero.
```

<details><summary>Function notes (the original docstrings)</summary>

**`BeachDuneConfig()`**

```text
Percent-scale settings for CASCADE's beach_dune_manager.

Attributes:
    community_overwash_filter_pct: Percent of overwash deposition removed
        from the interior in a developed community and returned to the
        shoreface. Rogers et al. (2015) give 40-90 %, residential to
        commercial.
    default_overwash_filter_pct: Percent for every other domain --
        undeveloped ground filters nothing, so this is normally 0.
    overwash_to_dune_pct: Percent of the remaining overwash bulldozed onto
        the dunes. filter + to_dune must stay under 100.
    minimum_community_width_m: BeachDuneManager's own hard-coded threshold,
        repeated here so pre-flight checks can use it. Below this average
        interior width the community is abandoned and management stops.
```

**`NourishmentProject()`**

```text
One historical beach nourishment project.

Volume is the project total as reported by the permitting agency, in cubic
yards, spread evenly across the project's domains. Even spreading is an
assumption: real fill templates taper at the ends, but the reported totals
are not resolved by domain.

Attributes:
    name: Project name, used in printed summaries and run metadata.
    year: Calendar year the fill was placed.
    gis_domains: GIS domains the project covered.
    volume_cubic_yards: Reported project total.
    note: Human-readable description for run metadata.
    enabled: Whether to include this project in a schedule.
```

**`NourishmentSchedule()`**

```text
Per-year nourishment flags and volumes on the padded array.

Attributes:
    nourish_now: Mapping of calendar year to a 0/1 float array of length
        geometry.total_domains, matching CASCADE's nourish_now.
    volume_m3_per_m: Mapping of calendar year to a float array of length
        geometry.total_domains, in m^3/m.
    projects: Projects that fall inside the period.
    skipped: (project, reason) pairs for projects that do not.
    geometry: DomainGeometry the arrays are indexed by.
    start_year: First model year.
    end_year: Last model year.
```

**`build_schedule()`**

```text
Expands nourishment projects onto per-year padded arrays.

Projects dated outside [start_year, end_year] are skipped rather than
raising, so one project list serves every period -- a 1984-2004 run simply
returns an empty schedule.

Args:
    projects: Iterable of NourishmentProject.
    geometry: DomainGeometry describing the padded array.
    start_year: First model year.
    end_year: Last model year.

Returns:
    A NourishmentSchedule.

Raises:
    ValueError: If a project covers a GIS domain outside the real span, or
        names no domains, or carries a negative volume.
```

**`build_overwash_filter()`**

```text
Per-domain overwash filter percent, for Cascade's `overwash_filter`.

Only developed ground filters overwash, so the community percent applies
inside the community zones and the default (normally 0) everywhere else --
including nourished domains outside a village, which are road corridor and
refuge rather than development.

Args:
    geometry: DomainGeometry describing the padded array.
    community_zones: Iterable of inclusive (first_gis, last_gis) spans.
    config: BeachDuneConfig supplying the two percentages.

Returns:
    A list of geometry.total_domains percentages.
```

**`build_beach_dune_management_on()`**

```text
Which domains run CASCADE's beach_dune_manager.

The union of two footprints drawn from different sources: the permanent
community zones, which need overwash filtering every year, and the
nourishment project extents, which need the module on at all to receive
fill. The project extents are wider than the villages, so the union is
larger than either.

Everything in this list runs the module's always-on behaviour for the whole
run -- fixed dune line, community-width drowning check -- not just in
nourishment years. See the module docstring.

Args:
    geometry: DomainGeometry describing the padded array.
    community_zones: Iterable of inclusive (first_gis, last_gis) spans.
    nourished_gis: Iterable of GIS domains that receive fill.
    enabled: Global switch; False turns the module off everywhere, which
        is what a natural-dynamics run wants. Mirrors the same argument on
        roadway.build_roadway_management_on.

Returns:
    A list of geometry.total_domains booleans.
```

**`find_double_managed()`**

```text
GIS domains where the roadway and beach-dune managers both run.

CASCADE does not guard against this. Both loops run in `cascade.update()`,
roadway first, and the combination has two consequences worth naming:
overwash is removed twice over (bulldozed, then the survivors filtered),
and BeachDuneManager holds `dune_migration_on` False, which leaves
`ShorelineChangeTS` at 0 -- the value RoadwayManager reads to decrement the
road setback. The road therefore stops retreating in these domains.

Args:
    beach_dune_on: Per-domain booleans for beach_nourishment_module.
    roadway_on: Per-domain booleans for roadway_management_module.
    geometry: DomainGeometry describing the padded array.

Returns:
    A sorted tuple of GIS domain IDs, empty if the footprints are disjoint.
```

**`audit_schedule()`**

```text
Pre-flight check that the schedule can actually reach the model.

Catches the failures that produce a plausible-looking run rather than an
error: a domain scheduled for fill whose module is off (the nourishment is
silently dropped), a volume outside any real project's range, or a schedule
year outside the run.

Args:
    schedule: The NourishmentSchedule about to be run.
    beach_dune_on: Per-domain booleans for beach_nourishment_module.
    config: BeachDuneConfig, for the reported drowning threshold.

Returns:
    A dict with 'ok', 'module_off' (scheduled but unreachable, as GIS IDs),
    'implausible_volume' (event rows), 'out_of_period' (event rows), and
    'notes' (human-readable strings).
```

**`verify_nourishment()`**

```text
Checks what BeachDuneManager actually did against the schedule.

Reads each manager's own `_nourishment_TS` and `_nourishment_volume_TS`,
which are written inside the nourishment branch of
`BeachDuneManager.update` from the volume passed to `shoreface_nourishment`
(beach_dune_manager.py:794-795). So this reports the volume the model
used -- not the volume the run script intended, which is what a
schedule-side log records.

That distinction is the point of the function. Writing the volume onto the
manager instead of onto the Cascade object produces a run where every
nourishment applies the Cascade init default while the log looks correct;
only the manager's own record disagrees.

Args:
    cascade: The Cascade instance after the run.
    schedule: The NourishmentSchedule the run was driven by.
    beach_dune_on: Per-domain booleans for beach_nourishment_module.
    tolerance_m3_per_m: Allowed volume difference before it counts as
        wrong.

Returns:
    A dict with 'ok', 'confirmed', 'missing', 'wrong_volume',
    'unexpected', 'abandoned', and 'notes'. Every list holds dicts with
    year, gis, and the expected/actual volumes where relevant.
```

**`verify_setbacks_frozen()`**

```text
Confirms the road setback stopped moving in double-managed domains.

The predicted consequence of running both managers in one domain, checked
against the run rather than assumed: BeachDuneManager holds
`dune_migration_on` False, Barrier3D then leaves `ShorelineChangeTS` at 0,
and RoadwayManager only decrements the setback when that value is non-zero.

Args:
    cascade: The Cascade instance after the run.
    double_managed_gis: GIS domains reported by `find_double_managed`.
    geometry: DomainGeometry describing the padded array.

Returns:
    A list of dicts (gis, setback_start_m, setback_end_m, frozen), one per
    double-managed domain, sorted by GIS ID.
```

**`volume_m3_per_m()`**

```text
Volume per unit alongshore length, the unit CASCADE wants.

Args:
    domain_spacing_m: Alongshore width of one domain, in meters.

Returns:
    Nourishment volume in m^3/m.
```

**`time_index()`**

```text
Index into a BeachDuneManager time series for a calendar year.

Barrier3D's `update_dune_domain` increments `time_index` before the
management modules run, so a manager writing at `time_index - 1` for
the first model year lands on index 1, not 0.

Args:
    year: Calendar year.

Returns:
    The index into `_nourishment_TS`, `_beach_width`, and the other
    per-time-step arrays on BeachDuneManager.
```

**`events()`**

```text
Every scheduled nourishment, one row per domain-year.

Returns:
    A list of dicts with year, time_index, pad, gis, and
    volume_m3_per_m, sorted by year then GIS domain.
```

**`apply_to_cascade()`**

```text
Sets this year's nourishment flags and volumes on a live Cascade.

Call this immediately BEFORE `cascade.update()` for the matching model
year. Both arrays are rewritten in full every year, so nothing carries
over from the previous one.

The volume goes to `cascade.nourishment_volume` rather than to the
BeachDuneManager instances: `cascade.update()` copies its own list into
each manager before calling it, so a value written onto the manager is
overwritten before use. See the module docstring.

Args:
    cascade: A live Cascade instance.
    year: Calendar year about to be stepped.
    default_volume_m3_per_m: Volume carried by domains with no event.
        Unused by the model -- `nourish_now` is 0 there and no
        nourishment interval is set -- so 0 keeps it unambiguous.

Returns:
    A list of dicts (year, pad, gis, volume_m3_per_m) for the domains
    nourished this year, in GIS order. Empty when nothing is scheduled.

Raises:
    ValueError: If the schedule is wider than the Cascade grid.
```

</details>

### plotting/__init__.py

cascade_pipeline.plotting: the rendering modules, animations and static figures.

From the script's original header:

```text
Rendering modules: shoreline_gif (animation) and rate_comparison (static figures).

Both import the shared foundation (domains, run_info, annotations,
coastsat_lowess) but never import each other -- keep it that way. If a
helper is needed by both, it belongs one level up, not duplicated here.
```

### plotting/init_planview.py

Plan-view rendering of a CASCADE initialization surface (t=0).

From the script's original header:

```text
Plan-view rendering of a CASCADE initialization surface (t=0).

Composites per-domain Barrier3D topography onto one alongshore canvas, each
domain shifted cross-shore by its BRIE island offset, so the initial island
can be read as a map before a run starts. Nothing here is site-specific --
domain geometry arrives as a DomainGeometry and the file paths and offsets are
supplied by the caller.

Two conventions this module depends on, both Barrier3D's:

- Elevation arrays are stored in decameters on a fixed-size grid cell, and are
  indexed (cross_shore, alongshore) with row 0 seaward.
- The extractor writes topography trimmed to the island interior, so the
  landward water rows are absent from the .npy files. `PlanViewConfig.topo_rows`
  is the untrimmed frame height; the missing rows are refilled with
  `sentinel_water_m` (see RUN_MANIFEST.txt in the extractor's version folder).
```

Notes that were in the code:

```text
No np.fliplr here. It used to reverse the alongshore cells WITHIN each
domain, inside this very loop, which put every 500 m block backwards
against the ascending domain order. That compensated for the extractor
writing the within-domain alongshore order reversed; the extractor now
fixes it at load (ALONGSHORE_FLIP), so flipping again double-flips.
_warn_if_alongshore_reversed below is the guard.
```

```text
Adaptive so a small multiple (a 7-domain zoom) still gets labels;
a fixed step of 5 gives such a panel only two ticks.
```

<details><summary>Function notes (the original docstrings)</summary>

**`PlanViewConfig()`**

```text
Rendering and unit settings for a plan-view initialization figure.

Attributes:
    topo_rows: Cross-shore rows in the extractor's untrimmed frame.
    sentinel_water_m: Elevation used to refill trimmed water rows.
    cell_size_m: Grid cell size in meters.
    dam_to_m: Decameters-to-meters factor for the stored arrays.
    elev_min_m: Colorbar floor.
    elev_max_m: Colorbar ceiling.
    sea_level_m: Elevation the colormap pivots at.
    sea_level_pos: Colormap position sea level maps to. Below 0.5 pushes
        green onto low above-water terrain instead of burying it in the
        sub-sea-level range.
    ocean_color: Fill for canvas cells no domain covers.
    text_color: Foreground color for titles, labels and ticks.
```

**`pad_cross_shore()`**

```text
Extends a trimmed interior to the full cross-shore frame.

Args:
    elevation_m: 2-D (cross_shore, alongshore) elevations in meters.
    config: PlanViewConfig supplying the frame height and water value.

Returns:
    The array padded (or truncated) to config.topo_rows rows. Rows run
    seaward to landward, so padding is appended on the landward end.
```

**`load_domain_grids()`**

```text
Loads every domain's elevation array, in meters, on the full frame.

Args:
    elevation_paths: Padded-order elevation file paths, one per domain.
    config: PlanViewConfig supplying unit and frame settings.

Returns:
    A list of 2-D arrays in meters, each config.topo_rows tall.
```

**`_warn_if_alongshore_reversed()`**

```text
Warns when the alongshore axis is reversed within each domain.

A per-domain alongshore reversal is nearly invisible in a 45 km plan view --
it reads as roughness -- but it lays every 500 m block down backwards. In
numbers it is unmistakable: the mean jump across domain seams over the mean
jump within a domain sits near 1 for a continuous island, and reached 21 on
the corrected arrays while this function's caller still applied np.fliplr.

Args:
    grids: The per-domain arrays about to be composited.
    config: PlanViewConfig supplying sentinel_water_m.
    warn_ratio: Ratio above which a warning is printed.

Returns:
    The seam/inner ratio, or nan when it cannot be computed.
```

**`build_canvas()`**

```text
Composites domains onto one canvas, offset by the island dune line.

Args:
    domain_grids: Padded-order elevation arrays from load_domain_grids.
    offset_cells: Per-domain cross-shore offset, in whole grid cells.
    geometry: DomainGeometry describing the padded array.
    include_buffers: Whether to draw the buffer domains as well as the
        real ones.
    config: PlanViewConfig supplying the frame height.

Returns:
    A (canvas, domain_col_starts, cells_per_domain, first_real_idx) tuple.
    Cells no domain covers are NaN. first_real_idx is the position of the
    first real domain within the plotted set.
```

**`plot_canvas()`**

```text
Draws a composited canvas as a plan-view map.

Args:
    canvas: Canvas from build_canvas.
    domain_col_starts: Column start of each plotted domain.
    cells_per_domain: Alongshore cell count of each plotted domain.
    first_real_idx: Index of the first real domain within the plotted set.
    geometry: DomainGeometry used to label the real domains.
    title: Axes title.
    ax: Axes to draw on. A new figure is created when omitted.
    include_buffers: Whether buffer domains are present, and so whether to
        bracket and label them.
    xlabel: X axis label.
    colorbar: Whether to draw the elevation colorbar. Off for small
        multiples that share one bar with a parent panel.
    config: PlanViewConfig supplying colors and limits.

Returns:
    The matplotlib Figure the canvas was drawn on.
```

**`plot_initialization()`**

```text
Loads, composites and draws an initialization surface in one call.

Args:
    elevation_paths: Padded-order elevation file paths, one per domain.
    offset_cells: Per-domain cross-shore offset, in whole grid cells.
    geometry: DomainGeometry describing the padded array.
    title: Axes title.
    include_buffers: Whether to draw the buffer domains.
    ax: Axes to draw on. A new figure is created when omitted.
    config: PlanViewConfig supplying unit and rendering settings.
    **kwargs: Passed through to plot_canvas.

Returns:
    The matplotlib Figure.
```

</details>

### plotting/rate_comparison.py

Modelled against CoastSat shoreline-change-rate figures.

From the script's original header:

```text
Modeled-vs-CoastSat shoreline-change-rate figures.

Two entry points:
  plot_rate_comparison           the working REAL-domains/ALL-domains QC
                                  figure (toggle via real_domains_only)
  plot_annotated_rate_comparison the publication/poster figure with the
                                  full geographic annotation layer

Both consume cs_series from cascade_pipeline.coastsat_lowess.build_coastsat_series
-- this module only renders, it doesn't load or smooth CoastSat data itself.

STYLE. Drawn under the house standard (scripts/site_layer/hat_figure_style.py): printed
width, 8-9 pt type, one alongshore axis label, village bands. The one place
these figures depart from it is the TITLE. Every other figure in the project
moves its title sentence into a CAPTIONS.md beside the image; these two are
written automatically into a run folder and opened months later with nothing
around them, so the run's identity -- period, scope, background erosion, run
name -- stays on the canvas. It is set small and muted beside the subject
rather than as a 14 pt banner, but it is not moved off the image and there is
no captions file in a run folder. The wave height and the SLR rate came OFF
that line on 2026-09-10 (Hannah): they crowded it, and both are in the run's
metadata JSON and TXT beside the PNG.
```

Notes that were in the code:

```text
`scripts/` is on sys.path already -- cascade_pipeline lives inside it, so
importing this package at all means the style module is importable too.
```

```text
COLOUR HERE IS NOT THE HOUSE SEMANTIC PAIR, ON PURPOSE (Hannah, 2026-09-10).
This is the figure read most often in the project and it has one reading:
the ORANGE model curve against the BLUE observation. Orange is the site
config's -- HATTERAS_ANNOTATIONS.model_color -- and this module reads it off
the `annotations` object every call site already passes, so the hex is
defined in exactly one place and stays there. The blues below are the
observed layer: the darker for the widest LOWESS window, the lighter for a
narrower one, and a muted blue for the per-transect cloud underneath.
The red/blue VINTAGE pair is a different case -- two periods drawn together
-- and deliberately does not appear here.
```

```text
The sensitivity style names the wave settings here, because its model
legend is just "CASCADE" (Hannah, 2026-09-29).
```

```text
After the limits, so a span outside the view is skipped and a label
is clamped to the visible part of its span.
```

```text
Extra pad: the title row sits above the secondary GIS axis, which
the parent axes' title placement knows nothing about.
```

```text
Scatter/LOWESS transition marker -- only when southern domains show
raw scatter only; marks where dots end and LOWESS lines begin.
```

```text
The unsmoothed means swing wider than the LOWESS at every peak
(GIS 68, 1996-2010: +5.1 against +0.1), so a bound locked to the
smoothed curves alone would cut the new line off (2026-09-22).
```

```text
The southern dots are drawn too, and at Cape Point they run past
every curve (2010-2024 D1: +9.1 against a +7.2 bound), so a bound
without them clipped the top transects off (2026-09-29).
```

```text
A white backing: the erosion label sits at the right edge, which is where
the observed curve runs on this period.
```

```text
The run's identity, kept on the canvas -- see the module docstring. The
separate italic "Model | Obs | SLR | Run" footnote this figure used to
carry beneath the legend said the same things twice; it is folded in
here, where the eye already is.
```

<details><summary>Function notes (the original docstrings)</summary>

**`RateComparisonConfig()`**

```text
Styling for the modeled-vs-CoastSat shoreline-change-rate figures.

Attributes:
    window_colors: {window_domains: color} for each LOWESS window.
    window_color_default: Fallback color for an unlisted window size.
    window_styles: [(linewidth, linestyle, alpha_factor), ...] matched
        to lowess_config.window_domains by position; extra windows fall
        back to (1.5, "-", 0.80).
    raw_color: Color for the individual-transect scatter.
        The modelled curve's colour is NOT here: it is the site config's
        `AnnotationConfig.model_color` (HATTERAS_ANNOTATIONS defines the
        orange), which is the single authority for it.
    plot_domain_means: Draw the UNSMOOTHED per-domain mean of the
        transect rates under the LOWESS curve, the whole island (Hannah,
        2026-09-22). Without it the only observed curve on a model figure
        is the smoothed one, and a reader cannot tell a real alongshore
        feature from one the 5 km window flattened -- the observation
        figures beside these (rates_figures.lrr_figures) draw the same
        two curves, so the two products no longer disagree.
    domain_mean_lw, domain_mean_alpha: Its weight; thinner and fainter
        than the LOWESS it underlies, in the same colour as that window.
    plot_raw_lrr: Show the transect scatter at all.
    raw_lrr_southern_only: True -> scatter only for the domains where
        LOWESS is suppressed (D1-lowess_config.skip_southern_domains).
        False -> scatter for every real domain. Read by
        plot_annotated_rate_comparison always, and by
        plot_coastsat_overlay only when overlay_raw_southern_only is set.
    overlay_raw_southern_only: Make plot_coastsat_overlay (the real-
        domains and with-buffers run figures) honour
        raw_lrr_southern_only too. Off by default so the matrix figures
        keep the whole-island scatter; the sensitivity cells turn it on
        with plot_domain_means=False (Hannah, 2026-09-29: "only keeping
        the 7 domain LOWESS, just use the dots for domains 1-10").
    south_mean_line: Draw the unsmoothed per-domain mean over
        D1..skip_southern_domains as a DASHED line, so the observed curve
        runs the whole island but the stretch that is NOT a LOWESS looks
        different from the stretch that is (Hannah, 2026-09-29: "add the
        line through D1-10 ... but ensure it isnt LOWESS").
    show_features: Draw the shoal zones, piers and groin on the
        real-domains figure (the annotated figure always has them).
    model_label_plain: Legend the model as just run.model_name
        ("CASCADE"); the wave settings go in the subtitle instead.
    show_wave_climate: Add run.wave_climate as a subtitle line.
    quantity: "rate" (m/yr, the default) or "position" -- the caller
        passes the model's end-minus-start change (m) as `change_rate`
        and a CoastSat series already scaled to metres
        (coastsat_lowess.scale_coastsat_series). Only the publication
        labels and caption read it.
    observed_label: Legend text for the observed LOWESS curve in the
        publication style; None keeps "CoastSat LRR (7-domain LOWESS)".
    observed_description: Caption phrase for the observed quantity in
        position mode, e.g. "total change: the per-transect LRR fitted on
        1996-2010 x 14 yr".
    publication_text: Replace the title and provenance lines with ONE
        short tag above the plot (period and wave settings), write the
        rest to CAPTIONS.md beside the PNG, and legend in the house
        wording (Hannah, 2026-09-29: "more academic and professional").
        Both the real-domains and the annotated (with-buffers) figure;
        the annotated one also drops its S/N compass text, which the
        x label already says and which sat where the tag goes.
    These four are off by default; rerender_run_figures.py --lowess-only
    turns them on together (the sensitivity style).
    plot_reference_period: Also show the CoastSat period that doesn't
        match the run's start year (faded).
    raw_scatter_size: Marker area, points^2.
    raw_scatter_alpha: Opacity for the active period; the reference
        period (if shown) uses 0.35x this.
    domain_tick_step: X-axis tick spacing, in GIS domains.
    ylim: (low, high) in m/yr to fix the y axis, or None (the default) to
        fit it to the data. Set when a SET of runs has to be read on one
        scale (Hannah, 2026-09-27: every figure in raw_runs/matrix on one
        y axis); rerender_run_figures.py --ylim passes it. Applies to
        the with-buffers figures.
    ylim_real: the same for the real-domains-only figure, which has no
        buffer swings to hold and can be tighter (Hannah, 2026-09-27:
        "so we are efficient with space"). None falls back to ylim.
        rerender_run_figures.py --ylim-real passes it.
```

**`coastsat_domain_mean()`**

```text
(gis ids, unsmoothed per-domain mean rate) for one CoastSat series.

The plain mean of the transect LRRs in each 500 m domain -- what
domain_lrr_summary.csv holds -- recomputed from the same transect arrays
the LOWESS is fitted to so the two curves on a figure are the same data
two ways, not two files.
```

**`south_domain_mean_line()`**

```text
The raw domain means over D1..skip, dashed, for one CoastSat series.

Not a LOWESS: south of skip the target IS the unsmoothed domain mean, and
the dashes say so against the solid LOWESS line north of it. Returns the
label drawn, or None when there is nothing to draw.
```

**`_features_only()`**

```text
The annotation layer minus the town spans and village lines, for a
figure whose towns are already drawn by town_bands. The groin label goes
beside the line at the bottom, south of it: D1-5 run positive in both
windows, so the lower-left corner is the one empty place near the groin.
```

**`plot_coastsat_overlay()`**

```text
Draw CoastSat raw-transect scatter + LOWESS lines for one axis.

Public because the sensitivity figures draw many model curves on one axis
and need the SAME observed layer underneath them; a second implementation
there would let the two drift apart on styling and on the southern splice.

Shared by both branches of plot_rate_comparison (REAL vs ALL domains);
the annotated figure has extra features (fill_between, legend handle
tracking) and builds its overlay separately.

Args:
    ax: Axes to draw on.
    cs_series: Output of build_coastsat_series.
    lowess_config: LowessConfig (only window_domains is read here, to
        find the widest window).
    config: RateComparisonConfig.
    x_transform: along_coast_m array -> x-axis coordinate array, for
        the raw scatter.
    gis_x_transform: gis-domain-ID array -> x-axis coordinate array,
        for the LOWESS lines. Defaults to identity (GIS-ID x-axis).
```

**`_rate_axis_label()`**

```text
Axis label naming the estimator the modelled curve was built with.

The observed curve is always an LRR -- CoastSat supplies a per-transect
OLS slope -- so a figure that does not say which estimator the MODEL
used is inviting the reader to assume they match. They only match when
estimator is "lrr".

Args:
    estimator: "lrr", "endpoint", or None to leave the estimator
        unnamed (the pre-2026-08-22 label, kept so an old figure can
        be redrawn unchanged).
    title_case: Match the annotated figure's title-case axis labels.

Returns:
    The y-axis label string.
```

**`_run_parameters()`**

```text
The run's identity, as two short lines for the provenance title.

What a reader needs in order to know WHICH run they are looking at, months
later, out of a run folder: the scope of the figure, whether background
erosion was on, and the run name. Plus the alongshore endpoints, which the
axis label no longer carries.

NOT the wave height and NOT the SLR rate (Hannah, 2026-09-10): both were
crowding the line, and both are in the run's metadata JSON and TXT beside
the figure, so nothing is lost by leaving them off the canvas. The period
belongs to the subject line, not here.

`endpoints` off for a figure that already names them at the axis corners.
```

**`_provenance()`**

```text
Subject on the title line, the run's identity stacked under it.

The house rule is that nothing on the canvas belongs in a caption; the
documented exception is a per-run artefact, which has no caption file to
be read beside (see the module docstring). Small and muted, not a banner.

STACKED, not set to the right of the subject: the run name alone is 36
characters, so a right-aligned block ran straight through the title at the
printed width (seen on the first render, 2026-09-10). The title pad is
computed from the number of lines so constrained_layout reserves the room
and the two never touch.
```

**`_tick_step()`**

```text
`step`, widened until `span` carries no more than `max_ticks` labels.

At the printed width a 120-domain axis ticked every 5 is a grey smear;
the config's step is the floor, not the answer, on the widest layout.
```

**`_publication_axes()`**

```text
Axis labels and the one-line tag (period and wave settings) that
replace the title and provenance lines.
```

**`_publication_legend()`**

```text
The legend in house wording (STYLE.md, 2026-09-19): short noun
phrases, no period or dataset repeats -- those are in the caption.
```

**`_publication_caption()`**

```text
What the title and the three provenance lines used to say, as a
caption for CAPTIONS.md (house rule: nothing on the canvas that belongs
in a caption).
```

**`plot_rate_comparison()`**

```text
Modeled vs. observed shoreline-change-rate figure (REAL or ALL domains).

Two layouts, chosen by real_domains_only:
  True  -> x-axis is GIS domain IDs (domains.first_gis_id to
           domains.last_gis_id) only. Community spans are drawn from
           annotations.town_spans, same source of truth as the
           annotated figure -- label text is whatever's stored there.
  False -> x-axis is all domains.total_domains padded indices, buffers
           shaded red, GIS IDs on a secondary top axis.

Args:
    change_rate: 1-D array, length domains.total_domains, m/yr (from
        cascade_pipeline.shoreline.compute_change_rate).
    cs_series: Output of build_coastsat_series.
    run: RunInfo for this run.
    real_domains_only: Selects the layout (see above).
    estimator: Which estimator built change_rate -- "lrr", "endpoint",
        or None to leave it unnamed. Names it on the y axis, since
        the observed curve is always an LRR and a figure that does
        not say invites the reader to assume the two match.
    sea_level_rise_rate_m_yr: Accepted and ignored by the figure
        since 2026-09-10 -- the SLR rate crowded the provenance
        line and is in the run's metadata beside the PNG. Kept in
        the signature so no caller has to change.
    save_path: If given, fig.savefig(save_path, dpi=300, bbox_inches="tight").
    show: Call plt.show() before returning.

Returns:
    (fig, ax, fig_suffix): fig_suffix is "REAL_DOMAINS_ONLY" or
    "ALL_DOMAINS_WITH_BUFFERS", handy for building the output filename.
```

**`plot_annotated_rate_comparison()`**

```text
Publication/poster figure: modeled rate + full geographic annotation layer.

Always uses the real-domains-only (GIS first_gis_id-last_gis_id) x-axis,
regardless of the REAL/ALL toggle used for plot_rate_comparison -- this
figure is for sharing, not QC.

Args:
    change_rate: 1-D array, length domains.total_domains, m/yr.
    cs_series: Output of build_coastsat_series.
    run: RunInfo for this run.
    estimator: Which estimator built change_rate -- "lrr",
        "endpoint", or None to leave it unnamed. See
        plot_rate_comparison.
    sea_level_rise_rate_m_yr: Accepted and ignored by the figure
        since 2026-09-10 -- the SLR rate crowded the provenance
        line and is in the run's metadata beside the PNG. Kept in
        the signature so no caller has to change.
    save_path: If given, fig.savefig(save_path, dpi=300,
        bbox_inches="tight", facecolor="white").
    show: Call plt.show() before returning.

Returns:
    (fig, ax)
```

</details>

### plotting/road_planview.py

Plan view of the roadway on a CASCADE initialization surface.

From the script's original header:

```text
Plan-view of the roadway on a CASCADE initialization surface.

Draws NC-12 where CASCADE will actually put it: on the same plan-view canvas
`init_planview` builds, with each domain shifted cross-shore by its island
offset, so the road's setback is read against the island the model runs on
rather than against a map.

The road is drawn per domain as a horizontal bar at `offset + setback`, because
that is literally what CASCADE does -- `road_start = int(setback / cell)` is one
row index applied to every alongshore profile in the domain. A relocated
position, when supplied, is drawn as a second bar so the displacement is
visible against the barrier that has to absorb it.

Nothing here is site-specific: geometry arrives as a DomainGeometry, and the
setbacks, offsets and events are supplied by the caller.
```

Notes that were in the code:

```text
Aliased: `figsize` is already a keyword argument of plot_roadway_island,
and no caller should have to change name for a restyle.
```

```text
The subject line stays on the canvas: this is a per-run artefact that is
opened out of a run folder with nothing around it.
```

<details><summary>Function notes (the original docstrings)</summary>

**`IslandSection()`**

```text
A named alongshore stretch, drawn as a labelled band.

Attributes:
    name: Label drawn above the band.
    first_gis: First GIS domain in the stretch.
    last_gis: Last GIS domain in the stretch.
    color: Band fill, or None to label the stretch without shading it.
```

**`RoadPlanViewStyle()`**

```text
Colors and sizing for the roadway overlay.

The three overlay colors carry identity and were checked for colour-vision
separation (worst adjacent OKLab dE 25.5 against a target of 8; normal
vision 26.8 against a floor of 15). `relocated` warns on contrast against a
light surface, so it is always accompanied by a legend entry and never used
as the only cue. They are the house colours since 2026-09-10: NC-12 is
C["ROAD"] wherever it is drawn in this project, the relocated position is
the one house orange that shoreline_gif.RELOCATION_COLOR also uses, and the
drowning marker is the cool pole of the vintage pair -- three hues that
still separate on luminance as well as on hue.

Attributes:
    road: Color for the initial road footprint.
    relocated: Color for the post-relocation position.
    drowning: Color for the marker on domains that drown at t=0.
    dune_line: Color for the dune line the setback is measured from.
    bar_linewidth: Line width of a road bar, in points.
    text_color: Foreground for titles and labels.
```

**`crop_to_data()`**

```text
Finds the rows a canvas actually carries data in.

The island offsets span kilometres, so a canvas built from them is mostly
empty. Cropping to the occupied rows is what turns the plan view from a
thin diagonal ribbon into a readable map; a fixed guess slices the south
end off, because that is where the offsets are largest.

Args:
    canvas: The plan-view canvas; uncovered cells are NaN.
    pad_rows: Rows of breathing room to keep either side.

Returns:
    A (first_row, last_row) tuple, suitable for set_ylim.
```

**`draw_sections()`**

```text
Shades and labels the named alongshore stretches.

Args:
    ax: Axes carrying the plan-view canvas.
    sections: IslandSection instances.
    geometry: DomainGeometry describing the padded array.
    col_starts: Canvas column where each plotted domain starts.
    cells_per_domain: Canvas column count for each plotted domain.
    first_real_idx: Position of the first real domain in the plotted set.
    include_buffers: Whether buffers are in the plotted set.
    style: RoadPlanViewStyle.
```

**`road_rows()`**

```text
Converts setbacks into canvas rows, the way bulldoze indexes them.

Args:
    setbacks_m: Padded per-domain setbacks, in meters.
    offset_cells: Padded per-domain cross-shore offset, in whole cells.
    geometry: DomainGeometry describing the padded array.
    config: PlanViewConfig supplying the cell size.

Returns:
    A float array of canvas row positions, one per padded domain. Domains
    with a zero setback (no road) are NaN.
```

**`overlay_roadway()`**

```text
Draws the roadway onto an existing plan-view axes.

Args:
    ax: Axes already carrying an init_planview canvas.
    setbacks_m: Padded per-domain setbacks, in meters.
    offset_cells: Padded per-domain cross-shore offset, in whole cells.
    geometry: DomainGeometry describing the padded array.
    col_starts: Canvas column where each plotted domain starts.
    cells_per_domain: Canvas column count for each plotted domain.
    first_real_idx: Position of the first real domain in the plotted set.
    relocated_m: Optional padded setbacks after relocation, in meters.
    drowning_gis: GIS domains whose road drowns at t=0.
    include_buffers: Whether buffers are in the plotted set.
    config: PlanViewConfig supplying the cell size.
    style: RoadPlanViewStyle.

Returns:
    The list of legend handles added.
```

**`plot_roadway_planview()`**

```text
Renders the initialization surface with the roadway drawn on it.

Args:
    elevation_paths: Padded-order elevation file paths, one per domain.
    offset_cells: Padded per-domain cross-shore offset, in whole cells.
    setbacks_m: Padded per-domain setbacks, in meters.
    geometry: DomainGeometry describing the padded array.
    title: Axes title.
    relocated_m: Optional padded setbacks after relocation, in meters.
    drowning_gis: GIS domains whose road drowns at t=0.
    sections: IslandSection instances for the alongshore bands.
    crop: Whether to trim the view to the rows carrying data. The offsets
        span kilometres, so an uncropped canvas is mostly empty.
    legend: Whether to draw the roadway legend.
    include_buffers: Whether to draw the buffer domains.
    ax: Axes to draw on. A new figure is created when omitted.
    config: PlanViewConfig supplying unit and rendering settings.
    style: RoadPlanViewStyle.
    **kwargs: Passed through to init_planview.plot_canvas.

Returns:
    The matplotlib Figure.
```

**`plot_roadway_island()`**

```text
Renders the island-wide roadway figure, with optional zoom panels.

The full-island panel is cropped to the occupied rows and drawn with
`aspect="auto"`, which exaggerates the cross-shore axis. Hatteras is ~45 km
long and ~1 km wide, so at true aspect the island is a hairline; the
exaggeration is what makes the road's position readable, and the cross-shore
axis should be read as indicative rather than to scale.

Args:
    elevation_paths: Padded-order elevation file paths, one per domain.
    offset_cells: Padded per-domain cross-shore offset, in whole cells.
    setbacks_m: Padded per-domain setbacks, in meters.
    geometry: DomainGeometry describing the padded array.
    title: Figure suptitle.
    relocated_m: Optional padded setbacks after relocation, in meters.
    drowning_gis: GIS domains whose road drowns at t=0.
    sections: IslandSection instances for the alongshore bands.
    zoom_windows: (first_gis, last_gis, label) tuples, each drawn as its
        own panel below the island.
    config: PlanViewConfig supplying unit and rendering settings.
    style: RoadPlanViewStyle.
    figsize: Figure size in inches, or None for the printed double
        column (the house width; it was a 20 x 11 in canvas until
        2026-09-10, which set its 9 pt type at 3 pt on a page).
    xlabel: Alongshore axis label.
    **kwargs: Passed through to init_planview.plot_canvas.

Returns:
    The matplotlib Figure.
```

</details>

### plotting/road_relocation_gif.py

Side-by-side animation of NC-12 against a migrating dune line.

From the script's original header:

```text
Side-by-side animation of NC-12 against a migrating dune line.

WHAT IS DRAWN
    One frame per model year, two panels sharing a year clock and a y-axis:

        left   relocations OFF -- roadway_manager decides on its own
        right  relocations ON  -- the measured historical displacements

    In each panel, per alongshore domain:

        dune line   the modelled shoreline, as displacement from year 0
        road        that line PLUS the domain's current setback
        marker      a domain that relocated in this frame's year

WHY THE ROAD IS DRAWN AS "DUNE LINE PLUS SETBACK"
    `road_setback` is the road's distance LANDWARD of the interior domain's
    seaward edge, and `roadway_manager` decrements it by dune migration every
    year precisely so the road stays geographically put while the dune line
    advances on it. On a landward-positive axis, `dune + setback` reproduces
    that: a road that is not moving traces a FLAT line while the dune line
    climbs toward it, and a relocation shows up as the road stepping up in a
    single frame. A road drawn at an absolute position would hide the very
    mechanism the animation exists to show.

    The y-axis is therefore cross-shore displacement relative to each domain's
    own year-0 dune line, not an absolute cross-shore coordinate. Two domains
    at the same height on the plot are NOT at the same place on the island.

SIGN CONVENTION
    x_s_TS increases landward, and `run.flip_sign_model` turns that into a
    seaward-positive quantity -- which shoreline_gif then draws on an inverted
    axis. This module negates once more instead, so the axis is plainly
    landward-positive and ascending: 0 is the year-0 dune line and larger is
    further landward. Do not pre-flip the matrix before passing it in.

WHERE A LINE STOPS
    A road that drowns stops being drawn. The managed span is dated from
    `_road_ele_TS`, never from the setback, because a setback of exactly
    0.0 m is legitimate -- see `_last_managed`.
```

Notes that were in the code:

```text
THE PALETTE IS THE HOUSE PALETTE (2026-09-10). The glyphs, the rings and the
panel wording are unchanged -- those were settled on 2026-09-09 and are
deliberate -- but every colour now comes from hat_figure_style, so NC-12 is
the same ink here as on the road plan view, the relocation marker is the
same orange as shoreline_gif's, and the ocean shoreline is the cool pole of
the vintage pair rather than a fourth near-identical blue.
```

```text
Deliberately not the star colour: a prescribed move and a module-triggered
one are different claims and must not share a glyph. Purple rather than the
old cyan, which sat a few degrees of hue from the dune line and read as a
marker ON that line at GIF resolution. Now the house purple, which is the
same hue at a hair more saturation.
```

```text
The barrier interior: the style module's 0.5-1.0 m elevation class, so the
island body is the same sand tone as on every elevation figure.
```

```text
One place for the typographic and axis conventions every frame in this module
shares. The VALUES are the house standard's now, applied per-artist as before
-- these functions are called from long analysis scripts that draw their own
figures, and per-artist styling means a frame looks the same whichever entry
point drew it. (apply_style() above sets the typeface and the rc defaults;
it is idempotent, and shoreline_gif calls it too, so importing this module
has been touching rcParams since 2026-09-10 either way.)
```

```text
Frames are drawn at the PRINTED double-column width and the dpi carries the
pixels: an 18 in canvas set this 8-9 pt type at 3-4 pt on a page. The target
pixel width is what it always was, so the GIFs are the same size on screen.
```

```text
the line figures have no room under the rule (the panel titles sit
there); the note goes in the header row, right-aligned against the
clock, so it neither jitters nor runs into a long title
```

```text
The relocation tracker (2026-09-09, Hannah): how many of the historically
relocated domains has each panel relocated BY THIS FRAME, and how many
relocations it made elsewhere - against how many the record says should
have happened by now. Island-wide, from the full road_series, whatever the
window; the window's own share is given beside it.
```

```text
The road stops being drawn where the record stops -- a drowned road
must leave a gap, not a line frozen at its last position.
```

```text
A STAR IS ALWAYS THE MODULE, IN BOTH PANELS. It is driven by
`_road_relocated_TS`, the RoadwayManager's own counter, and a PRESCRIBED
historical move never increments it -- the pipeline applies those as a
displacement before the manager updates. So arm B stars in the years its
module fired, not in 1989/1999, and the prescribed move shows only as a
step in the road. That was silently misleading, hence the separate
prescribed marker -- a RING, drawn around the star rather than over it,
so a year in which BOTH happened in the same domain still reads as both.
The distinction lives in the legend now: it did not fit in these panel
titles, where the two long labels collided across the gutter.
```

```text
One y-axis for both panels, fixed across all frames: a road that jumps
in one panel and not the other must be readable as a difference between
the arms, not as a difference between two autoscaled axes.
```

```text
Include the back-barrier or the island body is drawn clipped, which
reads as the sound being part of the barrier.
```

```text
Landward-positive (see _panel_series), so ocean-at-bottom is the plain
ascending axis rather than an inverted one.
```

```text
Fixed margins (NOT bbox_inches="tight") so every frame is the same
size -- mismatched frame dimensions break GIF assembly. Left margin
is wide enough for the y-label at every window width; the bottom
clears a two-row legend strip.
```

```text
Domains are integers; the default locator offers halves on a
narrow window, which reads as a domain that does not exist.
```

```text
The island itself: everything between the ocean dune line
and the back-barrier shoreline.
```

```text
The gap between the two lines IS the setback; shading it makes
"the dune line is closing on the road" the thing you see.
```

```text
The PRESCRIBED move, marked only in the arm that carries it and
only in its event year. Without this the bottom panel's headline
event was invisible: the road simply steps, with no marker at
all, while the stars nearby are the module doing something else.
```

```text
Left-aligned and lettered, so the panels can be cited as (a) and
(b) in a caption rather than by position.
```

```text
No arrow glyph: the label is rotated 90 degrees and a triangle
rotates with it, so it ends up pointing at the axis, not landward.
```

```text
Every mark on the frame, including the two shaded bands, which
were previously unexplained. Two rows of four, sized to fit the
narrowest window this function draws (9 in): a legend entry that
runs off the canvas is worse than no legend at all.
```

```text
Figure-level and below the axes: an in-axes legend sits on top of
the road wherever the setback is large, which is most of the island.
```

```text
dpi EXPLICIT: the house rcParams put savefig.dpi at 300 for print,
which would render every frame at nearly twice the pixels the GIF
is sized for (and blow up the file).
```

```text
The line version above abstracts the island to two curves. This one paints
Barrier3D's actual interior grids, so the animation shows overwash fans,
the barrier narrowing, and the road sitting on real ground.

GEOMETRY. Every domain carries `DomainTS[t]`, a (cross-shore rows, 50)
elevation grid in dam MHW whose row 0 is the seaward edge of the interior.
The domain's absolute cross-shore position is `x_s_TS[t]` (dam), so row r
lands at `x_s_TS[t] + r`. Alongshore, each domain is 50 cells of 10 m, i.e.
500 m, and the domains abut -- so 90 real domains tile 45 km of island.

WHAT IS WATER. The interior grid extends well past the subaerial barrier
(176 rows where the barrier is only ~32), so most of it is sound. Cells at
or below 0 m MHW are masked and drawn as water rather than as low land.
```

```text
Cross-shore window drawn around the window's reference shoreline, in dam.
The seaward end is fixed (a few cells of ocean for context); the landward
end is chosen per window from how far the land actually reaches, because a
fixed value either clips the wide Tri-Village section or drowns narrow Pea
Island in an empty frame.
```

```text
A PERCENTILE, not the maximum. One wide domain in the window -- a
Tri-Village section beside narrow Pea Island, say -- otherwise sets a
frame tall enough to make the domains under test unreadable. Clipping
the widest domain costs nothing here: the subject is the road corridor,
which sits within a few hundred metres of the ocean shoreline.
```

```text
One reference shoreline for the whole window, so the panel does not
shift under the island as the shoreline retreats.
```

```text
setback is metres landward of the interior's seaward edge; the
raster is in 10 m rows, so /10 puts it in the same units.
```

```text
Fixed for the whole animation and shared by both panels, so a colour or
a height means the same thing in every frame and in both arms.
```

```text
Last flag marks the arm carrying the prescribed moves; see
COLOR_PRESCRIBED and the label note on make_road_relocation_gif.
```

```text
NC-12 AS CASCADE HOLDS IT (Hannah, 2026-09-09): two straight rows
per domain, 20 m wide, spanning the domain's 50 cells - not a
line through the domain centres, which drew slopes between
domains that no cell of the model has.
```

```text
Cells are 10 m, so the axis is labelled in metres: "row 23"
is not a distance anyone can check against a map.
```

```text
The raster keeps all four spines: here the frame IS the edge
of the data, not decoration.
```

```text
This panel had NO legend at all: the star and the dotted lines were
drawn unexplained, and a reader had no way to tell whether a star
meant "the module decided to move the road" or "history did". The
panel titles now carry the per-panel meaning; these entries name the
glyphs themselves.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_style_axes()`**

```text
Applies the shared axis conventions in place.

Args:
    ax: The axes to style.
    grid_axis: "y", "x", "both", or None. Defaults to a HORIZONTAL-only
        grid: the vertical reference lines these panels draw mark real
        domains, and a vertical grid of the same weight competes with them.
    box: True keeps all four spines, for a raster panel where the frame is
        the edge of the data. False drops the top and right, which is the
        convention for the line panels.
```

**`_figure_header()`**

```text
Draws the title block: name on the left, year clock on the right.

Replaces a centred suptitle carrying both. A centred title moves as the
year text changes width and, in a GIF, that jitter is the first thing the
eye tracks; anchoring the two ends to the axes margins holds them still.

Args:
    fig: The figure.
    left: Left axes margin in figure coordinates.
    right: Right axes margin in figure coordinates.
    title: Figure title.
    year: Calendar year of the frame.
```

**`_panel_series()`**

```text
Year-0-referenced dune-line displacement for one arm, LANDWARD-POSITIVE.

shoreline_gif keeps a seaward-positive axis and then inverts the y-limits
to put the ocean at the bottom. That works for a single shoreline but reads
badly here, where the road sits tens of metres landward of the dune line
and would be drawn at large NEGATIVE numbers. This negates once more, so
landward is positive and up, the dune line starts at 0 and climbs as it
retreats, and the road is simply `dune + setback`.

Args:
    shoreline_m: 2-D [n_years, total_domains] raw x_s_TS matrix, metres.
    run: RunInfo carrying flip_sign_model.
    pad_lo: First padded index to keep.
    pad_hi: One past the last padded index to keep.

Returns:
    A 2-D array [n_years, pad_hi - pad_lo]; row 0 is identically zero and
    positive means landward of the year-0 dune line.
```

**`_referenced()`**

```text
Puts a second cross-shore series on the dune line's own axis.

The back-barrier shoreline has to be referenced to the SAME year-0 dune
line as `_panel_series`, not to its own year-0 position; otherwise both
curves start at zero and the island is drawn with no width at all.

Args:
    other_m: 2-D raw x_b_TS matrix in metres, same shape as shoreline_m.
    shoreline_m: The x_s_TS matrix supplying the reference.
    run: RunInfo carrying flip_sign_model.
    pad_lo: First padded index to keep.
    pad_hi: One past the last padded index to keep.

Returns:
    A 2-D landward-positive array on the dune line's axis.
```

**`_last_managed()`**

```text
Index of the last year the manager ran, dated from road ELEVATION.

Not from the setback series: a setback of exactly 0.0 m is legitimate --
it means the road sits on the dune line, which is where six of the ten
historical domains start -- so a nonzero test on the setback would draw no
road at all in precisely the domains this animation exists to show. Road
elevation has no such ambiguity; the module stops the moment it drops
below 0 m MHW.

Args:
    entry: One road_series() value, carrying "elevation".

Returns:
    The index, or -1 if the manager never ran.
```

**`_tracker()`**

```text
Counts for one panel at year index `t`.

A historical domain counts once it has relocated AT ALL by `t`, early or
late (cumulative; the timing error stays visible in the stars and rings).
In a panel that carried the prescribed moves, the applied move is read
OFF THE SETBACK SERIES - a jump of a cell or more at the event year - not
assumed from the event table, so a run that failed to apply an event
shows as a miss rather than a hit by construction. 'Elsewhere' is the
module's own relocations in domains the record does not list, i.e. the
false positives, cumulative.

Returns:
    dict(hist_hit, hist_total, hist_hit_win, hist_total_win, elsewhere,
         observed, observed_by_year)
```

**`_tracker_label()`**

```text
The per-panel line: 'historical 3 of 10 (window 2 of 4) - elsewhere 2'.

Short, as the line above says. The long form ("historical domains ... this
window ...") was 3.6 in of 7.5 pt type, so at the printed width it overran
its half-width panel and ran across the gutter into the other one.
```

**`_road_matrix()`**

```text
Builds road position and relocation-flag grids aligned to `series`.

Args:
    series: Dune-line displacement [n_years, n_domains] for this arm.
    road_series: {gis: {"setback", "relocated", ...}} for this arm.
    gis_lo: First GIS domain in the window.
    gis_hi: Last GIS domain in the window.
    domains: DomainGeometry.
    n_years: Number of animation years.

Returns:
    A (road, relocated) tuple of [n_years, n_domains] arrays. Domains the
    run did not manage are NaN in `road` and False in `relocated`, so an
    unmanaged domain draws no road rather than a road at the dune line.
```

**`make_road_relocation_gif()`**

```text
Writes the two-panel road/dune animation for one alongshore window.

`prescribed_panels` says, per panel, whether that run CARRIED the
prescribed 1989/1999 moves and so earns a ring at the event year. The
default is the relocation comparison's pairing (a free, b prescribed).
A pairing of two versions under ONE scenario passes (False, False) for
an emergent scenario and (True, True) for a prescribed one; before this
flag (2026-09-09) panel b was ringed unconditionally, which drew
"measured move applied" on a v3 run that applied nothing.

Args:
    arm_a: (shoreline_m, RunInfo) for the relocations-OFF run.
    arm_b: (shoreline_m, RunInfo) for the relocations-ON run.
    road_series_a: {gis: series dict} for arm A, from the comparison
        script's road_series().
    road_series_b: The same for arm B.
    gis_lo: First GIS domain to draw.
    gis_hi: Last GIS domain to draw.
    out_path: Destination .gif path.
    back_a: Optional raw x_b_TS matrix for arm A, metres. Supplying it
        draws the barrier body between the two shorelines instead of a
        single dune line, which is what makes the plot read as an
        island rather than as a chart.
    back_b: The same for arm B.
    label_a: Panel title for arm A.
    label_b: Panel title for arm B.
    event_years: {gis: year} marking the historical relocation year, drawn
        as a reference tick on the right panel.
    domains: DomainGeometry.
    annotations: Geographic annotation layer (town spans, groins).
    gif_config: Shared GifConfig.
    fps: Frames per second; defaults to gif_config.
    stride: Year stride; defaults to gif_config.
    title: Figure suptitle prefix.

Returns:
    The written path, or None if the animation was skipped.
```

**`make_all_road_gifs()`**

```text
Writes one animation per alongshore window.

Args:
    arm_a: (shoreline_m, RunInfo) for the relocations-OFF run.
    arm_b: (shoreline_m, RunInfo) for the relocations-ON run.
    road_series_a: {gis: series dict} for arm A.
    road_series_b: {gis: series dict} for arm B.
    windows: Sequence of (name, gis_lo, gis_hi) tuples.
    out_dir: Directory for the .gif files.
    event_years: {gis: year} for the historical relocation reference ticks.
    gif_config: Shared GifConfig.
    **kwargs: Forwarded to make_road_relocation_gif.

Returns:
    A list of written paths, skipped windows omitted.
```

**`_window_reference()`**

```text
Cross-shore origin shared by every panel, domain and year.

Taken as the most SEAWARD shoreline anywhere in the window, over both arms
and every drawn year, minus a few cells of open ocean. An earlier version
used the window's mean year-0 shoreline, which put domains seaward of that
mean at a negative row: `_island_raster` clipped their land to row 0 while
`_road_track` drew the road at its true negative position, so NC-12
appeared to float in the sea off Pea Island. A minimum cannot do that.

Args:
    cascades: Iterable of finished Cascade instances sharing the window.
    gis_lo: First GIS domain in the window.
    gis_hi: Last GIS domain in the window.
    domains: DomainGeometry.
    years: Iterable of year indices that will be drawn.

Returns:
    The cross-shore dam coordinate of raster row 0.
```

**`_cross_shore_rows()`**

```text
Chooses a landward extent that fits the land in this window.

Scans every drawn year and takes the furthest landward row that is still
above 0 m MHW anywhere in the window, then adds headroom. Fixed once for
the whole animation so the frame does not breathe year to year.

Args:
    cascade: A finished Cascade instance.
    gis_lo: First GIS domain in the window.
    gis_hi: Last GIS domain in the window.
    domains: DomainGeometry.
    years: Iterable of year indices that will be drawn.
    x_ref: Cross-shore origin from _window_reference.

Returns:
    The number of cross-shore rows to draw.
```

**`_island_raster()`**

```text
Assembles one year's elevation raster for an alongshore window.

Args:
    cascade: A finished Cascade instance.
    year: Index into the per-domain time series (0 = first model year).
    gis_lo: First GIS domain in the window.
    gis_hi: Last GIS domain in the window.
    domains: DomainGeometry.
    n_cross: Number of cross-shore rows to draw, from
        _cross_shore_rows.
    x_ref: Cross-shore origin from _window_reference.

Returns:
    A [n_cross, n_domains * 50] array in metres MHW, NaN where the model
    has no cell.
```

**`_road_track()`**

```text
Road cross-shore position for one year, in raster row units.

Args:
    cascade: A finished Cascade instance.
    road_series: {gis: series dict} for this arm.
    year: Index into the time series.
    gis_lo: First GIS domain in the window.
    gis_hi: Last GIS domain in the window.
    domains: DomainGeometry.
    x_ref: Cross-shore dam coordinate of raster row 0.

Returns:
    An (x_cells, rows, relocated_mask) tuple. `rows` is NaN wherever the
    domain is unmanaged or its record has stopped, so the road line breaks
    rather than being drawn across a drowned stretch.
```

**`_land_colormap()`**

```text
matplotlib's `terrain`, restricted to its LAND range.

Full `terrain` spends its first quarter on blues and cyans for bathymetry.
Everything drawn here is at or above 0 m MHW -- water is masked out and
painted separately -- so those blues would render low-lying land in the
same colour as the sound beside it, which is exactly the distinction the
animation exists to show. Starting at 0.25 gives the familiar green-to-
brown-to-white land ramp with nothing wasted below sea level.

Returns:
    A Colormap whose "bad" colour is the water fill.
```

**`make_topography_gif()`**

```text
Animated elevation map of the island with NC-12 drawn on it.

`prescribed_panels`: per panel, whether that run carried the prescribed
moves (ring at the event year); see make_road_relocation_gif.

Two panels sharing a colour scale and a year clock. Cells at or below
0 m MHW are drawn as water, so the barrier's real outline, its overwash
fans and its narrowing are all visible rather than implied.

VERTICAL EXAGGERATION. The window is ~640 m across-shore and 500 m per
domain alongshore, so the aspect is set from the data and stated on the
figure. Nothing here is to scale in the way a map is; it is a raster of
model cells.

Args:
    cascade_a: Finished Cascade for the relocations-OFF run.
    cascade_b: Finished Cascade for the relocations-ON run.
    road_series_a: {gis: series dict} for arm A.
    road_series_b: {gis: series dict} for arm B.
    gis_lo: First GIS domain to draw.
    gis_hi: Last GIS domain to draw.
    out_path: Destination .gif path.
    start_year: Calendar year of model year 0.
    label_a: Panel title for arm A.
    label_b: Panel title for arm B.
    event_years: {gis: year} for the historical relocation reference ticks.
    vmax_m: Top of the elevation colour scale, metres MHW. 3 m rather
        than the barrier's true maximum: almost every cell is under 2 m,
        so a scale topped by the few dune crests washes the island out.
    domains: DomainGeometry.
    gif_config: Shared GifConfig.
    fps: Frames per second; defaults to gif_config.
    stride: Year stride; defaults to gif_config.
    title: Figure suptitle prefix.
    planform_note: Optional caption, e.g. recording that the run's
        shoreline offset does not carry the island's true curvature.

Returns:
    The written path, or None if the animation was skipped.
```

</details>

### plotting/setup_qc.py

Pre-run QC figures for the Hatteras hindcast setup.

From the script's original header:

```text
Pre-run QC figures for the Hatteras hindcast setup.

WHY THIS MODULE EXISTS
    These three figures answer "does the initial condition look right" before
    a run is started: which way the island is oriented, what the
    initialization surface looks like in plan view, and how much sea level
    rises over the period. They are worth having and cost ~125 lines of
    matplotlib to define, so the definitions live here and the notebook keeps
    one call each.

    Nothing downstream reads their output. Skipping them changes no result --
    which is exactly why they belong out of the file that describes the run.
```

Notes that were in the code:

```text
The one alongshore axis label the whole project uses. The endpoints this
constant used to carry ("Cape Point to Pea Island") were the only CORRECT
pair anywhere in the repo, so they are not lost: ENDPOINT_NOTE states them
beside each figure's title instead.
```

```text
NOT shoreline change: each year is zeroed on its own most-seaward
domain, so this carries a constant offset. Pattern only; section
9.4 rebuilds the real change from the raw transect files. The mean
stays on the canvas: these figures are handed back to a notebook and
never written to disk, so there is no CAPTIONS.md to hold it.
```

<details><summary>Function notes (the original docstrings)</summary>

**`plot_island_orientation()`**

```text
Plots island offset by GIS domain for every period, active one bold.

Args:
    offsets_by_year: Mapping of start year to a padded offset array (dam).
    active_year: The start year currently selected.
    geometry: DomainGeometry used to slice out the real domains.

Returns:
    The matplotlib Figure.
```

**`plot_initialization_planview()`**

```text
Draws the initialization surface in plan view, with and without buffers.

Args:
    elevation_file_paths: Padded list of domain elevation .npy paths.
    island_offset_dam: Padded offsets in decameters.
    geometry: DomainGeometry.
    start_year: Start year, for the titles.
    config: PlanViewConfig, or None for the extractor defaults
        (200 rows, -3.0 m water).
    verbose: Whether to print each canvas shape.

Returns:
    The matplotlib Figure.
```

**`plot_sea_level_rise()`**

```text
Plots cumulative RSLR for each period from its own start year.

Args:
    periods: Mapping of start year to a period config dict.
    active_year: The start year currently selected.

Returns:
    The matplotlib Figure.
```

</details>

### plotting/shoreline_gif.py

Plan-view shoreline animation (GIF) for a completed CASCADE run.

From the script's original header:

```text
Plan-view shoreline animation (GIF) for a completed CASCADE run.

One frame per model year; current shoreline drawn over a year-0 reference
(dashed grey), shaded blue where seaward of the reference, red where
landward. See make_shoreline_gif's docstring for the three y-axis modes.

STYLE. House standard (scripts/site_layer/hat_figure_style.py) for type, colour and
frame. TWO deliberate departures, both because this is a per-run artefact
that lands in a run folder and is opened later with no context around it:
the run's identity stays on the canvas as a small muted line beside the
subject, and the year clock stays in the corner of every frame. There is no
CAPTIONS.md in a run folder.

FRAME SIZE. A GIF needs every frame to be the same pixel size, so the frames
are drawn at a fixed figure size with fixed margins (never bbox_inches=
"tight") and the pixel width comes from the dpi, not from a 16-inch canvas.
The dpi is passed to savefig EXPLICITLY: the house rcParams set savefig.dpi
to 300 for print, which would otherwise triple every frame.
```

Notes that were in the code:

```text
Shared with road_planview.RoadPlanViewStyle.relocated and the plan-view
animation, so a relocation reads as the same event in all three. The house
orange, so the three of them and road_relocation_gif's relocation star are
now ONE orange across the project instead of three near-misses.
```

```text
Frames are sized to a target pixel WIDTH: the canvas is the printed double
column, and the dpi carries the pixels. Matches road_relocation_gif.
```

```text
-- Optional observed/target line ----------------------------------------
Referenced exactly as the run is, so it shares the run's axis. "position"
subtracts the same year-0 alongshore mean; "displacement" subtracts the
same per-domain year-0 position, which turns the target into the
displacement the run is trying to reproduce. Both lines therefore sit on
the same year-0 base, so the gap between them is the model's misfit even
though that base is itself an initial condition.
```

```text
The run's identity, kept on the canvas -- see the module docstring. No
wave height and no SLR rate (Hannah, 2026-09-10): both crowded the line
and both are in the run's metadata beside the GIF.
```

```text
Fixed margins (NOT bbox_inches="tight") so every frame is byte-identical
in size -- mismatched frame dimensions break GIF assembly.
```

```text
Drawn after the limits are set, below: town_bands clamps its
labels to the visible part of each span and needs the view.
```

```text
-- Shoreline + erosion/accretion shading -----------------------------
The RdBu band fills: light blue where the shoreline lies seaward of
the reference, light red where it lies landward.
```

```text
-- Roadway relocations, in the year they happen ----------------------
An EVENT, not a state: _road_relocated_TS is raised in the year the
roadway manager moves the road, so a frame shows that year's moves
only. Drawn on the shoreline itself because that is the line whose
arrival at the road caused them.
```

```text
The compass ends, on the axis-label row rather than on top of the
tick labels (they landed on "10" and "90" at the printed width).
ONLY on a window that actually reaches both ends of the island: a
groin zoom over GIS 1-15 was being labelled "Pea Island | N" at its
right-hand edge, which is 37 km from Pea Island.
```

```text
-- Year clock + provenance + legend ----------------------------------
Both stay on the canvas: this frame is a run artefact, not a figure
with a caption beside it.
```

```text
Legend outside (below) the axes so it never covers the town / shoal
labels, which sit at fixed axes fractions inside the panel.
THREE columns, not five: five long labels ran off the right edge of
the printed width and the last one was clipped.
```

```text
-- Render to an in-memory PNG (constant size, no disk churn) ----------
dpi EXPLICIT: the house rcParams put savefig.dpi at 300 for print,
and a frame rendered at 300 dpi is three times the pixels the GIF
wants (and no longer matches a frame drawn anywhere else).
```

```text
A "groin" job fans out into one GIF per structure in
annotations.groins, so adding a second groin needs no change here
or in `jobs`. Use "which" to restrict, or range="groin_span" for
one enclosing window.
```

```text
A job's own "target" wins over the shared one; "target": False opts
out, which is how a job stays a clean model-only figure.
```

<details><summary>Function notes (the original docstrings)</summary>

**`GifConfig()`**

```text
Shared settings for every shoreline-animation job.

Attributes:
    fps: Frames per second.
    year_stride: 1 = every model year, 2 = every other year, etc.
    annotate: Draw the geographic annotation layer (GIS-axis modes).
    auto_open: Pop the finished GIF open in the OS default viewer.
    keep_frames: Also write individual frame PNGs to disk.
    save_matrix: Save shoreline_m.npy alongside the GIFs in
        make_all_shoreline_gifs, so this run can serve as a future
        baseline for a "difference" job.
    ocean_at_bottom: True -> ocean/seaward at the bottom of the plot,
        sound/landward at the top (Hatteras' real cross-shore layout).
        False -> seaward at the top.
    baseline_label: Shown in a "difference" job's title/y-axis label.
    target_label: Legend label for the optional observed/target line
        drawn by a job that passes target_m.
```

**`_groin_zoom_window()`**

```text
Inclusive GIS window centred on one groin, clipped to real domains.

ANN_GROINS-style values sit on a domain boundary (e.g. 5.5 = the 5/6
interface), so the window is built around that interface, not a cell
centre.
```

**`_select_groins()`**

```text
Resolve a job's "which" key into [(name, domain), ...] from annotations.groins.

Args:
    which: "all" (default, every groin -- new structures are picked up
        with no config change), a single name, or a list of names.
    annotations: AnnotationConfig.
```

**`_resolve_gif_domain_range()`**

```text
Turn a GIF job's "range" value into concrete plotting coordinates.

Args:
    spec: "real" | "all" | "groin" | "groin_span" | (gis_lo, gis_hi).
    pad: half-width in domains; only read for "groin"/"groin_span".
    groin: (name, domain) tuple, required when spec == "groin". Supplied
        by make_all_shoreline_gifs, which fans one "groin" job out into
        one call per structure.
    domains: DomainGeometry.
    annotations: AnnotationConfig (for "groin"/"groin_span").

Returns:
    (pad_lo, pad_hi, axis_kind, x_lo, x_hi, tag):
        pad_lo, pad_hi: padded-array slice bounds (pad_hi exclusive).
        axis_kind: "gis" (x-axis in GIS domain IDs) or "pad" (padded
            index).
        x_lo, x_hi: inclusive x-axis limits in whichever space
            axis_kind names.
        tag: short string for the output filename.
```

**`make_shoreline_gif()`**

```text
Plan-view animation of the shoreline, one frame per model year.

Current shoreline drawn over a year-0 reference (dashed grey); shaded
blue where seaward of the reference, red where landward. Axis
orientation is set by gif_config.ocean_at_bottom: ocean at bottom /
sound at top by default, matching Hatteras' real cross-shore layout.

SIGN CONVENTION: build_shoreline_matrix() returns raw x_s_TS (x10 for
m), NO sign flip applied -- in BRIE/Barrier3D, x_s_TS INCREASES as the
shoreline retreats landward. run.flip_sign_model is applied internally
here so larger = more seaward, matching compute_change_rate's sign
handling. Do not pre-flip shoreline_m before passing it in.

MODE:
  "displacement" (default) - each domain's change from its OWN year-0
    position; all start at 0, so a few m/yr stays legible across all 90
    domains. Use for the full island.
  "position" - cross-shore position relative to year-0 alongshore mean.
    Keeps planform shape, but x_s_TS is referenced to each domain's own
    local grid, so the alongshore "shape" is partly an artefact of how
    the initial DEMs were extracted, and ~80 m of 40-yr change can be
    nearly invisible against a planform range of several hundred m.
    Fine for a narrow groin-zone window like (1, 30); use
    "displacement" for the full island.
  "difference" - this run minus baseline_m, in displacement terms:
    (run - run_yr0) - (base - base_yr0). Differencing displacements
    (not raw positions) cancels any per-domain initial-condition offset
    between runs, so year 0 is exactly zero. With a no-groin baseline
    this isolates the groin's effect: blue = seaward of baseline
    (updrift fillet), red = landward (downdrift notch).

Args:
    shoreline_m: 2-D array [n_years, domains.total_domains] from
        cascade_pipeline.shoreline.build_shoreline_matrix(cascade,
        to_meters=True); raw x_s_TS convention.
    run: RunInfo for this run.
    domain_range: "real" | "all" | "groin" | "groin_span" | (gis_lo, gis_hi).
    pad: half-width in domains; only read when domain_range is
        "groin"/"groin_span".
    groin: (name, domain) tuple, required for domain_range="groin".
    baseline_m: same-shape matrix from a previous run (same raw x_s_TS
        convention); required for mode="difference".
    target_m: 1-D array of length domains.total_domains giving an
        observed/target shoreline position to draw as a static line in
        every frame -- e.g. the surveyed end-year island position, so
        "how close did we get" is readable off the animation. Same raw
        x_s_TS convention as shoreline_m (meters, increasing landward);
        run.flip_sign_model and the mode's own referencing are applied to
        it exactly as they are to the run, so it lands in the same frame
        as the shoreline. Ignored by mode="difference", where the y-axis
        is a run-minus-run quantity the target has no place in.
    target_label: legend label for that line; defaults to
        gif_config.target_label.
    Everything else defaults to gif_config.

Returns:
    Path to the saved GIF, or None if the job was skipped.
```

**`make_all_shoreline_gifs()`**

```text
Run every job in `jobs` for one completed run.

Also saves the shoreline matrix (when gif_config.save_matrix) so this
run can serve as the baseline for a later difference GIF. Each job is
independent: one failing (bad range, missing baseline) logs and is
skipped rather than taking the others -- or the run -- down with it.

Args:
    shoreline_m: [n_years, domains.total_domains] array, raw x_s_TS
        convention (from build_shoreline_matrix, to_meters=True).
    run: RunInfo for this run.
    jobs: List of dicts, same shape as the original GIF_JOBS list, e.g.
        [dict(range="real", mode="displacement"),
         dict(range="groin", mode="difference", pad=9)]
    baseline_npy: Path to a previous run's *_shoreline_matrix.npy to
        difference against. None skips "difference" jobs.
    target_m: Observed/target shoreline position (raw x_s_TS convention,
        meters, one value per padded domain) drawn as a static line on
        every "position"/"displacement" job. A job may override it with
        its own "target" key, or opt out with "target": False.
    relocations: Optional (n_years, domains.total_domains) boolean array,
        True in the year a domain's roadway is relocated. Passed to every
        job; each marks only the domains inside its own window. None -- the
        default, and what a run without roadway management should pass --
        draws no markers and leaves the GIFs exactly as they were.

Returns:
    List of saved GIF paths.
```

</details>

### reports.py

The printed reports the Hatteras hindcast emits at each section.

From the script's original header:

```text
The printed reports the Hatteras hindcast emits at each section.

WHY THIS MODULE EXISTS
    The reports are the run's evidence: they are what says a switch reached
    the module it names, what the sediment budget implies, and where the
    model disagrees with the survey. They are not debug output, so none of
    them were dropped -- but ~400 lines of f-strings sat duplicated between
    `HAT_hindcast_1984_2024.ipynb` and its headless mirror
    `HAT_hindcast_1984_2024.py`, and every wording change had to be made
    twice or the two runs stopped saying the same thing.

    Every function here prints exactly what the inline block printed. The
    values still come from the calling file -- these take arguments, never
    module globals -- so the notebook remains the place that decides what is
    reported, and this is only where the formatting lives.

ONE RULE
    A report never computes a result the run depends on. `run_units_check`
    is the single exception and it is named for it: it loops the domains and
    returns whether anything failed, because the loop existed only to be
    printed. Everything else takes already-built objects.
```

Notes that were in the code:

```text
PERIOD["enable_nourishment"] describes the period; nothing reads it. It
is printed beside the switch that does reach the model so the two cannot
quietly disagree -- the volume there is a legacy default that section 6
overwrites. Flagged in one direction only: a period that expects fill and
is not getting it is a scenario choice worth restating, whereas fills
left on in a period with no projects withholds nothing.
```

<details><summary>Function notes (the original docstrings)</summary>

**`run_units_check()`**

```text
Runs the units contract over every real domain and reports the tally.

Args:
    elevation_paths, dune_paths: Padded lists of .npy paths.
    contract: Mapping from load_barrier3d_contract.
    geometry: DomainGeometry.

Returns:
    True if every check passed on every domain.
```

**`scenario_report()`**

```text
Prints the scenario, what it expands to, and the predicted run name.

Args:
    scenario: Selected scenario key.
    scenarios: The full scenario preset table.
    departures: Mapping of key -> (wanted, got) for overridden switches.
    roadway_on, relocations_on, beach_dune_on, fills_on, groin_enabled:
        The resolved switches.
    relocations_forced_off, fills_forced_off: Whether the switch was
        forced off by its parent module rather than chosen.
    period_expects_nourishment: PERIOD["enable_nourishment"], printed
        beside the switch that does reach the model.
    run_name_preview: The name predicted from the switches.
    input_files: Sequence of (label, Path) for the period's inputs.
```

**`road_audit_report()`**

```text
Prints which roads would not survive year one, and why.

A domain whose road starts in water is not a warning: roadway_manager
sets _drown_break and returns immediately on every later year, so the
domain becomes an unmanaged barrier wearing a road label.
```

**`groin_report()`**

```text
Prints the groin configuration, its sediment budget, and the overlap.

For a blocking groin (groin_kind "blocking") the schedule is of b, the
intercepted fraction, and the budget is emergent -- it depends on the
transport arriving -- so it is reported after the run, not here.

The double-count note is reported rather than corrected: the calibrated
background-erosion rates were fit to CoastSat spanning the functional-groin
era, so the groin's signature is plausibly already in them.
```

**`scenario_summary_report()`**

```text
Prints every switch, the derived run name, and the two known overlaps.

Args:
    switches: The (label, value, token) list the name was derived from.
    beach_dune_on_updrift: Whether beach_dune_manager runs on the updrift
        domain. Keyed on the module flag, not the schedule: dune_migration
        is held False wherever the manager runs, fill year or not.
```

**`pre_run_report()`**

```text
Prints what was built, the resolved profile, and the fillet prediction.

The prediction is printed BEFORE the run because a prediction printed
after one is not a prediction. Amplitude was tuned; extent was not.
```

**`nourishment_report()`**

```text
Prints the model's own record of what fill it received.

Read from each manager's _nourishment_volume_TS, not from the schedule's
intent -- the failure this catches produced plausible-looking output for
a long time.
```

**`roadway_outcome_report()`**

```text
Prints which managed roads drowned and which relocations were blocked.

Returns:
    The (drowned, relocation_blocked) row lists. Section 12.4 counts both
    into run_index.csv, so they are returned rather than recomputed.
```

</details>

### roadway.py

Roadway forcing for a CASCADE run: setbacks, elevations, events, drowning.

From the script's original header:

```text
Roadway forcing for a CASCADE run: setbacks, elevations, events, drowning.

Everything CASCADE's `roadway_manager` needs, prepared before the run and
checked against the interior it will actually be spent on. Nothing here is
site-specific -- domain geometry arrives as a `DomainGeometry`, and the
Hatteras instances (community zones, historical events, file names) live in
`hatteras_site_config`.

Three Barrier3D conventions this module depends on:

- `road_setback` and `road_width` are METERS; `bulldoze` divides them by the
  cell size to get row indices, so a setback is only ever as precise as one
  cell.
- `road_ele` is METERS, MHW-RELATIVE. `bulldoze` writes it straight into the
  interior grid, which the extractor stores MHW-relative -- see the datum note
  on `load_road_elevations`.
- `bulldoze` tests the rows FLANKING the road, never the road's own cells, and
  drowns the road when either flank is more than `percent_water` water.
  `predict_drowning` reproduces that test rather than approximating it.
```

Notes that were in the code:

```text
The relocation TARGET is deliberately not touched. It used to be set
here "to keep the two in sync", which never survived a year: the
yearly update re-assigns it from the model's own parameter. Since
2026-09-14 that parameter is a real input rather than the initial
setback, so overwriting it here would quietly give this one domain a
different rebuild rule from every other.
```

<details><summary>Function notes (the original docstrings)</summary>

**`RoadwayConfig()`**

```text
Unit and threshold settings matching CASCADE's roadway_manager.

Attributes:
    road_width_m: Roadway width in meters, as passed to Cascade.
    cell_size_m: Grid cell size in meters (Barrier3D's dx/dy).
    dam_to_m: Decameters-to-meters factor for the stored arrays.
    drown_threshold_m: Elevation at or below which a cell counts as water.
        roadway_manager passes 0 and calls it "m MSL", but compares it
        against MHW-relative elevations, so it is effectively 0 m MHW.
    percent_water: Fraction of flanking cells that may be water before the
        road drowns (roadway_manager's percent_water_cells_touching_road).
    sea_level_dam: Barrier3D's SL. Fixed at 0 in the Lagrangian frame --
        sea-level rise lowers the domain rather than raising SL.
    elevation_fallback_m: road_ele for domains with no measured value.
```

**`RelocationEvent()`**

```text
A historical roadway relocation, as a DISPLACEMENT.

The road is moved landward by a measured distance; the resulting setback is
that displacement added to whatever setback the model is carrying at the
time. It is deliberately not an absolute setback: an absolute value
referenced to an older dune line double-counts the shoreline retreat
between that line and the topography's own dune line, which is what put
NC-12 into the sound in earlier versions of this pipeline.

Attributes:
    year: Calendar year the relocation happened.
    displacement_m: Landward displacement per GIS domain, in meters.
    note: Human-readable description for run metadata.
    enabled: Whether to apply this event.
```

**`BridgeEvent()`**

```text
A bridge or alternate route replacing the road surface.

Roadway management stops for the listed domains from this year on, with no
setback change -- the road is gone, not moved.

Attributes:
    year: Calendar year the bridge opened.
    gis_domains: GIS domains whose road is removed.
    note: Human-readable description for run metadata.
    enabled: Whether to apply this event.
```

**`load_padded_series()`**

```text
Reads a 2-row CASCADE forcing file into a padded per-domain array.

The 2-row format is (GIS IDs, values). CASCADE consumes forcings as arrays
indexed by padded position, so the IDs are used to place values rather than
being discarded -- reading positionally would silently shift every domain
after a gap.

Args:
    path: Path to the 2-row CSV.
    geometry: DomainGeometry describing the padded array.
    first_gis: First GIS domain the forcing covers.
    last_gis: Last GIS domain the forcing covers.
    fill: Value for domains the file does not cover.

Returns:
    A (values, missing) tuple. values is a float array of length
    geometry.total_domains; missing lists GIS domains in the expected span
    that the file did not supply.

Raises:
    ValueError: If the file is not 2 rows.
```

**`load_road_setbacks()`**

```text
Loads road setbacks, in meters landward of the dune line.

Args:
    path: Path to the 2-row RoadSetback_<year>.csv.
    geometry: DomainGeometry describing the padded array.
    first_gis: First GIS domain carrying road.
    last_gis: Last GIS domain carrying road.

Returns:
    A (setbacks_m, missing) tuple; setbacks are 0.0 outside the road span.
```

**`load_road_elevations()`**

```text
Loads per-domain road elevation, in meters MHW-relative.

Datum: `bulldoze` writes road_ele into the interior grid after dividing by
dz, and the extractor stores that grid MHW-relative. So road_ele must be
MHW-relative meters -- NOT NAVD88. The file this reads is already in that
frame; do not subtract MHW again.

Not period-dependent. Road elevation is a property of the surveyed surface,
and there is one topography for every period, so one elevation set serves
all of them.

Args:
    path: Path to the 2-row RoadElevation.csv. No year in the name --
        see the note above on why one set serves every period.
    geometry: DomainGeometry describing the padded array.
    first_gis: First GIS domain carrying road.
    last_gis: Last GIS domain carrying road.
    config: RoadwayConfig supplying the fallback elevation.

Returns:
    A (elevations_m, missing) tuple. Domains the file omits, and every
    buffer domain, get config.elevation_fallback_m.
```

**`build_roadway_management_on()`**

```text
Builds the per-domain roadway-management flag CASCADE consumes.

Management runs on the road span minus the permanent community zones.
Buffer domains are always off.

Note CASCADE constructs a RoadwayManager for EVERY domain regardless of
this flag ("always initialize just in case we want to add a road during the
simulation"), so the existence of a manager says nothing about whether the
road is managed -- only this flag does.

Args:
    geometry: DomainGeometry describing the padded array.
    first_gis: First GIS domain carrying road.
    last_gis: Last GIS domain carrying road.
    community_zones: (first_gis, last_gis) pairs to exclude.
    enabled: Global switch; False turns management off everywhere.

Returns:
    A list of bools, length geometry.total_domains.
```

**`interior_widths()`**

```text
Land run from row 0 to the first water cell, per alongshore profile.

Transcribes barrier3d.FindWidths, including its `- 1` and its clamp at
zero. This is the only width Barrier3D uses to decide where the island is;
land behind a water cell is invisible to it.

Args:
    interior_m: Interior elevations in meters MHW, (cross_shore, along).
    config: RoadwayConfig supplying the sea level.

Returns:
    An int array of land-run lengths, one per alongshore profile.
```

**`predict_drowning()`**

```text
Reproduces roadway_manager.bulldoze's drowning test at t=0.

bulldoze checks the rows FLANKING the road -- `road_end + 1` on the bay
side and `road_start - 1` on the sea side -- and never inspects the road's
own cells. The row at `road_end` is skipped entirely: it supplies only the
cell count. All three are reported so a caller can see the difference.

Args:
    interior_m: Interior elevations in meters MHW, (cross_shore, along).
    setback_m: Road setback in meters.
    config: RoadwayConfig supplying widths and thresholds.

Returns:
    A dict with road_start, road_end, border_row, n_rows, sea_water,
    bay_water, road_cells_water, drowns and wall. `wall` is None unless the
    setback would corrupt rather than drown the run: "NEGATIVE" (numpy
    wraps the index to the back of the barrier) or "OVERRUN" (the border
    row is past the array, which raises IndexError inside bulldoze).
```

**`audit_setbacks()`**

```text
Predicts, before the run, which road_offset will not survive year one.

A drowned road is not a warning: roadway_manager sets _drown_break and
returns immediately on every later year, so the domain gets no overwash
removal, no dune rebuilding and no relocation for the rest of the run. It
becomes an unmanaged barrier wearing a road label, which is why these
domains have to be named before any managed-vs-unmanaged comparison.

Args:
    elevation_paths: Padded-order interior .npy paths, one per domain.
    setbacks_m: Padded per-domain setbacks, in meters.
    geometry: DomainGeometry describing the padded array.
    first_gis: First GIS domain carrying road.
    last_gis: Last GIS domain carrying road.
    management_on: Optional padded management flags; when given, the
        `managed` column records whether CASCADE will actually run the
        manager for that domain.
    config: RoadwayConfig supplying widths and thresholds.

Returns:
    A list of per-domain dicts sorted by GIS id, each carrying gis, pad,
    setback_m, the predict_drowning fields, interior-width statistics and
    `managed`.
```

**`summarise_audit()`**

```text
Reduces audit_setbacks output to the counts worth reporting.

Args:
    rows: Output of audit_setbacks.
    config: RoadwayConfig supplying the water threshold, for the message.

Returns:
    A dict with n_domains, n_managed, drowning (GIS ids that will drown AND
    are managed), drowning_unmanaged, and blocking (GIS ids whose setback
    would corrupt the run rather than drown it).
```

**`relocated_setbacks()`**

```text
Applies a relocation event to a padded setback array.

Adds the measured displacement to the CURRENT setback. At t=0 the current
setback is the initial one, so this is what the road's position would be if
the event fired immediately -- the upper bound on the relocated setback,
since CASCADE decrements the setback by dune migration between t=0 and the
event year.

Args:
    setbacks_m: Padded per-domain setbacks, in meters.
    event: A RelocationEvent.
    geometry: DomainGeometry describing the padded array.

Returns:
    A new padded setback array; domains the event does not touch are
    unchanged.
```

**`summarise_road_management()`**

```text
Reports which domains kept their road, from CASCADE's post-run state.

Reads state CASCADE already exposes rather than instrumenting the model:
roadway_management_module says which domains were managed at all,
drown_break and relocation_break say why management stopped, and
_road_ele_TS is written every managed year so its last non-zero entry dates
the last managed year.

Args:
    cascade: A Cascade instance after its run.
    geometry: DomainGeometry describing the padded array.
    first_gis: First GIS domain carrying road.
    last_gis: Last GIS domain carrying road.

Returns:
    A list of per-domain dicts for the managed domains only, each with gis,
    pad, drowned, relocation_blocked, last_managed_year, reason,
    overwash_removed_m3, dunes_rebuilt and relocations.
```

**`apply_historical_event()`**

```text
Applies one historical roadway event to a running CASCADE model.

Called from inside the time loop when `event.year` comes up. Both event
types mutate live model state, and each does so in a way that needs
explaining:

`RelocationEvent` adds its measured DISPLACEMENT to the setback the model
is currently carrying. `roadway_manager` has been decrementing that setback
by dune migration since t=0, so it already holds the modelled retreat;
adding the displacement counts the retreat exactly once. An absolute
setback referenced to an older dune line would count it twice.

`road_ele` is deliberately left alone. A real relocation rebuilds the road
at grade on new ground, so resetting it looks right -- but it would be an
exact no-op. `road_ele` is initialised from the 2004 alignment, which IS
the post-relocation road for events preceding 2004, and `roadway_manager`
decrements it in the Lagrangian frame, so at year t it already holds
`measured_2004 - sum(RSLR[0:t]) * 10`. Rebuilding at grade on that same
alignment gives the identical number. This only breaks if a
post-relocation alignment has a different measured elevation than the
initial one, which is not the case here by construction.

CASCADE's own relocation guards are evaluated but NOT obeyed. The
model-driven relocation path refuses a move when the road would drown or
the island is too narrow; this prescribed path does not, because these
relocations actually happened and refusing them would be historically
wrong. Each refusal is recorded in the returned row so the disagreement
goes on the record rather than disappearing.

`BridgeEvent` switches `cascade.roadway_management_module` off for its
domains. It does NOT set `cascade.roadways[pad] = None`: `Cascade.update()`
reads `self._roadways[iB3D].drown_break` on every domain with no None
check, so nulling the object raises AttributeError on the next step.

Args:
    cascade: A Cascade instance mid-run.
    event: A RelocationEvent or BridgeEvent whose year has arrived.
    geometry: DomainGeometry describing the padded array.
    relocations_enabled: Global toggle for relocation events; a False
        value skips RelocationEvents and reports them as skipped. Bridge
        events ignore it.
    setback_check: Optional {gis: measured_setback_m} for an independent
        cross-check, printed by the caller alongside the new setback.

Returns:
    A list of per-domain dicts. Relocation rows carry gis, pad, kind,
    old_setback_m, displacement_m, new_setback_m, check_m and warnings;
    bridge rows carry gis, pad and kind. An event that is skipped or
    touches no managed domain returns an empty list.
```

**`_apply_bridge()`**

```text
Switches roadway management off for a BridgeEvent's domains.

Flips `cascade.roadway_management_module` rather than tracking the pads in
a set the model never reads -- which is what an earlier version did, so the
bridged domains kept full road management for the rest of the run.
```

</details>

### run_info.py

Identifying information for one completed CASCADE run.

From the script's original header:

```text
Identifying info for one completed CASCADE run.

Bundles the handful of values every plotting function in cascade_pipeline needs
(run name/dir, period, wave height, sign convention) so call sites pass one
object instead of five loose scalars pulled from module globals.
```

<details><summary>Function notes (the original docstrings)</summary>

**`RunInfo()`**

```text
Identifying info for one completed CASCADE run.

Attributes:
    run_name: Run name used in output filenames (e.g. run_name_hs).
    run_dir: Output directory for this run's files.
    start_year: Calendar year the run starts (period START_YEAR).
    end_year: Calendar year the run ends.
    Hs: Fixed significant wave height used for this run (m), shown in
        figure titles/legends. None omits it.
    flip_sign_model: CASCADE's x_s_TS increases landward (erosion).
        True (the convention used throughout cascade_pipeline's plotting)
        flips model output so a larger value means more seaward.
    background_erosion_on: Whether DOMAIN_BE_RATES has any non-zero
        entries for this run; used only for the "BE=on/off" label in
        figure titles and the GIF caption.
    model_name: Model name shown in figure titles/captions, e.g.
        "CASCADE". cascade_pipeline.shoreline is itself CASCADE-specific
        (it reads cascade.barrier3d / b3d.x_s_TS directly), so this
        defaults to "CASCADE" rather than a fully generic placeholder --
        override it if you're comparing against a different model run
        through the same figures.
    wave_climate: The run's wave settings as one line, e.g. "Hs 2.0 m,
        Tp 7.5 s, asym 0.6, high-angle 0.5". Shown in the subtitle only
        when RateComparisonConfig.show_wave_climate is set. None omits it.
```

</details>

### run_layout.py

Where a run's files live inside its run folder, and what they are called.

From the script's original header:

```text
Where a run's files live inside its run folder, and what they are called.

WHY THIS MODULE EXISTS
    A run folder used to be thirteen files in a heap, every one of them
    prefixed with the full run name:

        HAT_1984_2004_calibBE_road_bdm_groin/
            HAT_1984_2004_calibBE_road_bdm_groin.npz
            HAT_1984_2004_calibBE_road_bdm_groin_annotated.png
            HAT_1984_2004_calibBE_road_bdm_groin_shoreline_change_rate.csv
            ... ten more

    Since 2026-09-10 the output is sorted by kind, and the files inside those
    subfolders drop the run-name prefix, which the folder already carries:

        HAT_1984_2004_calibBE_road_bdm_groin/
            HAT_..._run_metadata.json      identity and state stay at the root:
            HAT_..._run_metadata.txt       these are globbed ACROSS runs and
            HAT_....npz                    travel outside the folder, so they
            HAT_..._shoreline_matrix.npy   keep the prefix
            figures/     shoreline_change_rate.png, ..._with_buffers.png
            animations/  displacement_domains_1-90.gif, ...
            tables/      shoreline_change_rate.csv, groin_diagnostics.csv,
                         road_management.csv

    The rename is not cosmetic. The longest path in the runs tree was 229
    characters and Windows gives up at 260; adding a subfolder took it to 240,
    leaving no headroom for a longer run name. The prefix was 56 of those
    characters, and dropping it inside the subfolders brings the worst case
    back to about 184.

THE FALLBACK IS THE POINT
    `resolve()` looks for the new location first and falls back to the old
    flat one. So a half-migrated tree always reads correctly, the migration
    can be interrupted, and a run folder restored from an old backup still
    works. Nothing in the codebase should join a run-folder filename by hand;
    call `resolve()` to read and `write_path()` to write.

USAGE
    from cascade_pipeline.run_layout import resolve, write_path, migrate_run

    p = resolve(run_dir, "rate_csv", run_name)          # read, either layout
    p = write_path(run_dir, "figure_rate", run_name)    # write, new layout
    migrate_run(run_dir)                                # move an old folder
```

Notes that were in the code:

```text
kind -> (subfolder or "" for the run root,
new filename (no run-name prefix unless it is a root file),
legacy filename as a "{run}" template)
```

```text
--- figures ----------------------------------------------------------
"annotated" never said what the figure was: it is the rate comparison
drawn across the buffer domains as well as the real ones.
```

```text
An animation's name is built from its mode and window rather than listed,
because the windows are open-ended (any GIS range is legal). These turn the
plotting module's range tag into the new form:
D1-90                          -> domains_1-90
groinZoom_BuxtonGroin_D1-15    -> buxton_groin_domains_1-15
groinSpan_D1-15                -> groin_span_domains_1-15
ALL120pad                      -> all_120_padded
```

<details><summary>Function notes (the original docstrings)</summary>

**`resolve()`**

```text
Where to READ this file: new layout if present, else the old flat one.

Returns the new-layout path when neither exists, so a caller that is about
to create the file gets the right place. `must_exist=True` raises instead.
```

**`plan_run()`**

```text
(source, destination) for every file in this run folder that moves.

Only files whose OLD name is present and whose NEW path is free are
planned, so running this twice is a no-op and an interrupted migration
resumes cleanly.
```

</details>

### run_registry.py

Run provenance, output-directory guarding, and the cross-run index.

From the script's original header:

```text
Run provenance, output-directory guarding, and the cross-run index.

The hindcast is run as a matrix -- period x source/sink preset x groin -- and
the three failure modes that costs are all bookkeeping ones:

  * a re-run silently overwriting the outputs of the run it was meant to be
    compared against (guard_run_dir),
  * a finished run whose metadata cannot say which code produced it, because
    the sandbox flag or the extractor version was flipped days earlier
    (git_provenance, and the [identity] section callers pass),
  * twelve runs whose results can only be compared by opening twelve
    hand-formatted text files (rebuild_run_index, which derives run_index.csv
    from every run's metadata; append_run_index is the pre-2026-09-16 writer
    and is kept for the groin sweep).

Metadata is written twice from ONE structure: a .txt to read and a .json to
parse. Rendering both from the same `sections` mapping is what keeps them from
disagreeing -- the previous inline version built only the prose, so anything
downstream had to re-parse it.

Used by both HAT_hindcast_1984_2024.ipynb and its headless mirror
HAT_hindcast_1984_2024.py, so the two cannot drift apart.
```

Notes that were in the code:

```text
Files whose presence means a run directory holds real output. A directory
containing only these is treated as empty, so a stray .gitkeep or an
editor's .DS_Store does not block a run.
```

```text
Paths a run WROTE BACK into the repository, so their being modified says
nothing about whether the run's code was committed. See git_provenance.
Since 2026-09-16 the runner copies the parameters yaml into the run
directory and CASCADE rewrites only the copy, so this file should stay
clean; it is still excluded so a tree dirtied by an older run reports right.
```

```text
The route_overwash index fix (2026-09-24). Barrier3D's subaerial test read
Elevation[TS, i, d+1:d+10] -- row and column swapped -- which read the wrong
cells and, on narrow domains, out of bounds (the silent crashes). Fixed on
the local Barrier3D branch fix/route-overwash-axis-swap; see
output/raw_runs/experiments/code-checks/2026-09-24-metres-3-barrier3d-overwash-fix/NOTE.md.
```

```text
The 2026-09-28 adoption (Barrier3D branch hatteras/adopted): the three
overwash fixes -- all present, all absent, or mixed (None) -- and
per-cell dune ceilings (DuneCeilingFromStart).
```

```text
Sorted and formatted rather than hashed off repr(): dict order and float
repr are not things to make a run's identity depend on.
```

```text
One run's outputs are at

<raw_runs>/[<arm>/]<start>_<end>/<preset>/<run_name>/

and the ARM COMPONENT IS ABSENT for the calibration arm. That asymmetry is
deliberate -- every run made before forcing arms existed is at the short
path, and emitting the component unconditionally would rename all of them --
but it does mean the path cannot be built by joining a fixed number of parts.

THIS IS THE ONLY PLACE THAT SPELLING BELONGS. It was previously rebuilt by
hand in six scripts, five of which predate the arm component and so join
<period>/<preset>/<name> with no slot for it: an arm-scoped run is simply
invisible to them. plot_sensitivity.py skipped its target-window check
silently whenever the path did not resolve, which is the quiet-wrong-path
failure hat_topo_version.py exists to end for the domain arrays, one tree
over. A name that is not on disk is an error here, and the error names the
arms the run IS under.
```

```text
A run is filed by what it is FOR, then by period and preset:

matrix/<period>/<preset>/<run_name>/                  the production runs
sensitivity/<axis>/<period>/<preset>/<run_name>_<token>/   sweep cells
experiments/<tag>/<period>/<preset>/<run_name>/       one question each
versions/<tag>/<period>/<preset>/<run_name>/          input-version pairs
archive/<tag>/<period>/<preset>/<run_name>/           superseded, intact

The run NAME still describes the scenario and is derived from the switches
the runner built. What used to be an "arm" -- one string that meant a
forcing value, an input version or an ad hoc experiment label without saying
which -- is now a KIND (one of KINDS) and a TAG. A matrix run has no tag. A
sensitivity cell's tag is its axis folder, derived from the trailing token on
its name, so the value it was swept to is in the NAME and the axis is in the
PATH. An experiment's or version's tag is the folder it was filed under,
"<set>/<member>" for a set of related runs.

Why (Hannah, 2026-09-16): a wave sweep had fanned out into twelve top-level
arms/waveHs<x>/ folders, each holding one run; nineteen of thirty arms were
finished one-off experiments nothing marked as finished; and the layout
could not tell those from the version comparisons that are kept on purpose.

BOTH EARLIER LAYOUTS STILL READ. find_run_dir tries the purpose path, then
the 2026-09-10 layout (arms/<arm>/..., <period>/<preset>/sweeps/<family>/),
then the flat pre-09-10 one, so a tree that is half migrated resolves. The
`arm=` keyword is accepted everywhere as a legacy spelling and translated
through LEGACY_ARMS.
```

```text
The legacy name of the unscoped tree. Still what a pre-09-16 index row says
in its `arm` column, and what old call sites pass.
```

```text
Trailing name tokens that mark a sensitivity cell, and the axis folder each
files under. The runner derives the token (cascade_pipeline.hindcast's
wave_climate_token / relocation_setback_token); a token combining two wave
fields ("waveHs3Tp10") files under the FIRST family it starts with.
```

```text
Where each pre-09-16 arm was filed on 2026-09-16, decided by Hannah in the
same interview: the three version comparisons keep their names under
versions/; everything else was a one-off experiment and is filed under
experiments/ by the date it was run, with the arm's own name as the member.
The twelve waveHs<x> arms are ABSENT on purpose: those were the 1996 wave
cells filed by the 2026-09-01 rule, re-run as sensitivity cells and deleted.
```

```text
A period directory is exactly <4 digits>_<4 digits>. Anything else directly
under a kind folder is a tag, so the two levels can be told apart without a
registry of tag names.
```

```text
A pre-09-16 caller passing the ARM positionally, where kind now
sits (find_run_dir(raw, name, period, preset, "version-pair/v2")).
```

```text
The purpose layout: kind folders, each holding tags (one or two levels)
or, for the matrix, periods directly.
```

```text
run_index.csv is a RESTATEMENT of every run's metadata JSON in one table, so
a question across runs is one read. Since 2026-09-16 no run appends to it:
each run writes the row it would have appended INTO its metadata, under
"index row", and the runner then calls rebuild_run_index, which regenerates
the whole file from every metadata on disk and replaces it atomically. Two
runs finishing at once both rebuild the same complete table, so the last
writer wins nothing -- which is what lets runs be concurrent. Rows for runs
that predate the "index row" section are carried over from the existing
file by their old key, so nothing is lost in the changeover.
```

```text
The old key only means something in a file that still HAS the arm
column. On a rebuilt file every row's arm reads "", so three runs of
one name collapse onto one old key and the last wins -- which handed
the 1996 matrix run a version row's skill on 2026-09-16.
```

```text
BACKFILL: a legacy run's row, once found, is written into its own
metadata JSON so the next rebuild derives it from the run and not
from whatever the index file happens to hold. Without this, a
rebuild that read a file lacking the arm column matched three
runs of one name to one row (2026-09-16). The .txt is left alone.
```

```text
Column order: identity first, then everything else in first-seen order,
so the file stays readable as columns accumulate.
```

```text
csv, not pandas: pandas round-trips every float through repr, which once
rewrote unrelated rows. Atomic: written beside, then replaced, so a
reader never sees a half-written file and two writers cannot interleave.
```

```text
WINDOWS: os.replace fails with "Access is denied" while another
process has the target open -- which is exactly when a parallel run
is replacing it too (2026-09-24: a finished run died here, its
outputs already written). The index is derived, so waiting a moment
and trying again is safe; the last writer's table is complete.
```

```text
Refused rather than handled: every file a run writes is flat, so a
subdirectory here means run_dir is not pointing where it is thought
to be -- a half-built path, a period directory, the output root. The
cost of guessing wrong is a recursive delete of someone's runs.
A run directory holds a KNOWN set of subfolders (run_layout's
figures/ animations/ tables/, plus the gif_frames scratch) and
nothing else. Any OTHER subdirectory still means run_dir is not
pointing where it is thought to be -- a half-built path, a period
directory, the output root -- and the cost of guessing wrong is a
recursive delete of someone's runs, so that stays refused. Before
2026-09-10 EVERY subdirectory was refused, which stopped OVERWRITE
working at all once the layout gained folders.
```

```text
One column width for the whole file, wide enough for the longest key, so
the "=" stays aligned across sections instead of jogging in and out.
```

```text
A COMPOSITE key is allowed because the run name no longer identifies a
run on its own: forcing that is not part of the scenario -- Hs -- scopes
the output DIRECTORY instead of adding a name token, so two runs can share
a name and differ in what they were forced with. Replacing on name alone
would silently drop one of them from the index.
```

```text
AS TEXT. Parsing the existing rows into pandas and writing them
back rewrites every float through repr, which on 2026-09-10 silently
truncated the last digit of five columns in a row this call was not
touching. Only the row being added should change.
```

```text
Stable ordering: the identity columns first, then whatever else exists,
so the file stays readable as columns accumulate.
```

<details><summary>Function notes (the original docstrings)</summary>

**`git_provenance()`**

```text
Records which commit produced a run, and whether the tree was clean.

A dirty tree is not an error -- most runs happen mid-edit -- but it does
mean the commit hash alone will not reproduce the run, so the flag is
recorded beside it rather than inferred later.

ONE PATH IS EXCLUDED, AND WITHOUT IT THE FLAG IS USELESS.
data/hatteras_init/Hatteras-CASCADE-parameters.yaml is TRACKED and is
REWRITTEN BY EVERY CASCADE CONSTRUCTION -- it is the shared file behind
the "never run two sweep orchestrators at once" rule. So the moment any
run starts the tree is dirty and stays dirty, and `dirty` was True on
every report this pipeline had ever written, including runs whose code
was fully committed. A flag that cannot be False carries no information.

Excluding it makes DIRTY TREE mean what it is read as meaning:
uncommitted CODE or INPUTS, not the run's own scratch output. Found
2026-08-31, when five relocation comparisons were stamped DIRTY TREE and
the only run-relevant dirty paths turned out to be this yaml and one
benign helper addition.

Anything else volatile that a run writes back into the repository belongs
in EXCLUDED_FROM_DIRTY too -- otherwise it re-breaks the flag silently.

Args:
    repo_root: Path to the repository root.

Returns:
    A dict with commit, branch and dirty. Values are the string "unknown"
    (and dirty None) if git is unavailable or this is not a repository, so
    provenance capture can never be the thing that fails a run.
```

**`barrier3d_provenance()`**

```text
Which Barrier3D a run used, and whether it carries the overwash fix.

Barrier3D is installed editable, so the branch checked out in its own
repository IS the model -- and CASCADE's git_provenance cannot see it. A
`git checkout master` there would silently put every later run back on the
indexing bug. This records the branch, commit and dirty flag, and checks
the SOURCE OF THE MODULE ACTUALLY IMPORTED, not the file on disk, for the
fixed line.

Returns:
    A dict with branch, commit, dirty and route_overwash_fix (True, False,
    or None if the source could not be read). Never raises: provenance
    capture must not be what fails a run.
```

**`values_digest()`**

```text
Fingerprints a {domain: rate} mapping.

A preset's NAME and its VALUES are separate facts. The source/sink table
is edited between runs -- re-solving an end domain, testing a different
edge value -- while the preset keeps the name "edgeBE", so a run that
records only the name cannot be told apart from the trial before it. Two
scalar columns cover the end domains; this covers everything else,
including the interior domains of a 90-domain calibrated fit that no
column carries.

Args:
    mapping: Mapping of domain id to rate.
    length: Characters of hex digest to keep. 12 is ~10^-14 collision
        odds over any realistic number of runs.

Returns:
    A hex string, or "empty" when the mapping is empty -- a preset that
    imposes nothing has nothing to fingerprint, and saying so directly
    beats a hash of the empty string that looks like a real value.
```

**`period_component()`**

```text
The <start>_<end> path component, from either spelling of a period.

Callers hold a period as a string in some places and as two integers in
others; accepting both is what lets every call site pass what it already
has rather than reformatting at each one.

Args:
    period: Either "1984_2004" or a (start_year, end_year) pair.

Returns:
    The directory component as a string.

Raises:
    ValueError: If the string is not <4 digits>_<4 digits>, or the pair is
        not two values. A malformed period would otherwise build a path
        that cannot exist, and be reported as a missing run.
```

**`sweep_family()`**

```text
The sensitivity axis a run belongs to, or "" if it is a scenario run.

Args:
    run_name: The run's derived name, which is also its directory name.

Returns:
    One of SWEEP_FAMILIES, or "" when the trailing token names none of
    them -- which is every scenario run.
```

**`check_tag()`**

```text
Validates a tag: one to four path components joined by '/'.

Two levels is a set and its members -- the four probes of one experiment,
the two versions of one pair. A third level exists for one thing only:
the steps of an iterative solve inside a member (`<set>/<member>/step2`,
2026-09-16), so a Newton solve's history sits under the member it
belongs to rather than beside it as `<member>-step2`, which put 26
folders at one level for a five-member experiment. A fourth since
2026-09-25: experiments are grouped by theme (Hannah), so a tag is
`<theme>/<study>/<member>[/<step>]`. The value reaches
here from an environment variable and is joined onto the output root, so
it must not escape it.

Args:
    tag: The tag string. None and "" mean no tag.

Returns:
    The tag, stripped; "" for none.

Raises:
    ValueError: If it is not one to four clean path components.
```

**`legacy_arm_to_kind_tag()`**

```text
(kind, tag) for a pre-2026-09-16 arm name.

Args:
    arm: The `arm` column of an old index row, or the value an old call
        site passes as arm=. None, "" and "calibration" are the matrix.

Returns:
    (kind, tag). A wave arm (`waveHs1p2`) maps to a sensitivity cell
    under the family its name starts with; an arm in LEGACY_ARMS maps as
    that table says; anything else is an experiment tagged with the arm
    name itself, so an unknown arm still resolves somewhere sensible.
```

**`_resolve_identity()`**

```text
Normalises the three spellings a caller may use into (kind, tag).

`arm=` wins when given, because a call site still passing it is one that
has not been updated and means the OLD thing. For a sensitivity cell the
tag is the axis folder, derived from the name when not given.
```

**`preset_dir_for()`**

```text
The directory holding every run of one period, preset, kind and tag.

This is the runner's OUTPUT_BASE_DIR. It is the level anything that
ENUMERATES runs works at -- the scenario grid, the relocation comparison,
the source/sink calibration -- as against `run_dir_for`, which is for a run
already named.

Args:
    raw_runs: The output/raw_runs root.
    period: "1984_2004" or (1984, 2004).
    preset: Source/sink preset, e.g. "calibBE".
    kind: One of KINDS; the default is the matrix.
    tag: The experiment/version tag, or the axis for a sensitivity run.
    arm: LEGACY. A pre-09-16 arm name, translated through LEGACY_ARMS.

Returns:
    The preset directory as a Path. Does not check it exists.
```

**`run_dir_for()`**

```text
The directory one run's output belongs in. Does not check it exists.

The inverse of the runner's RUN_DIR, and the only place the layout is
spelled. Use `find_run_dir` to READ a finished run; this builds the path a
run would be WRITTEN to, which is what a writer and a collision guard need
and what a reader should not be doing by hand.

Args:
    raw_runs: The output/raw_runs root.
    run_name: The run's derived name, which is also its directory name.
    period: "1984_2004" or (1984, 2004).
    preset: Source/sink preset, e.g. "calibBE".
    kind: One of KINDS; the default is the matrix.
    tag: The experiment/version tag. For a sensitivity run it is derived
        from the name's token when not given.
    arm: LEGACY. A pre-09-16 arm name.

Returns:
    The run directory as a Path.
```

**`legacy_run_dirs_for()`**

```text
Every place this run could sit under the two earlier layouts.

2026-09-10 layout:  arms/<arm>/<period>/<preset>/<run>  for an arm,
                    <period>/<preset>/sweeps/<family>/<run>  for a cell,
                    <period>/<preset>/<run>  otherwise.
Before that:        <arm>/<period>/<preset>/<run>  and  <period>/<preset>/<run>.

Args:
    raw_runs: The output/raw_runs root.
    run_name: The run's directory name.
    period: "1984_2004" or (1984, 2004).
    preset: Source/sink preset.
    arm: The legacy arm the run was filed under; None/"" for calibration.

Returns:
    Candidate Paths, most recent layout first. Not checked for existence.
```

**`kinds_holding()`**

```text
Every (kind, tag) under which this run exists on disk, either layout.

A run name describes the SCENARIO, so one name can legitimately exist
under several tags -- the matrix run and the version pair made from it.
This is what makes that discoverable rather than a surprise: it is what
`find_run_dir` reports when the place it was asked for holds nothing.

Args:
    raw_runs: The output/raw_runs root.
    run_name: The run's directory name.
    period: "1984_2004" or (1984, 2004).
    preset: Source/sink preset.

Returns:
    Sorted list of (kind, tag) pairs. Empty if the name is nowhere under
    this period and preset.
```

**`find_run_dir()`**

```text
Locates a finished run, raising with what IS on disk if it is absent.

KIND DEFAULTS TO THE MATRIX RATHER THAN SEARCHING. A search would let a
figure silently draw a run forced or versioned differently from the one it
names -- exactly what the kind and tag exist to prevent. Naming nothing
means the matrix, which is also what every call site did before arms
existed, so routing an existing script through this cannot change which
run it reads.

Args:
    raw_runs: The output/raw_runs root.
    run_name: The run's directory name.
    period: "1984_2004" or (1984, 2004).
    preset: Source/sink preset.
    kind: One of KINDS. Defaults to the matrix.
    tag: The experiment/version tag, or the axis of a sensitivity cell
        (derived from the name when not given).
    arm: LEGACY. A pre-09-16 arm name; translated, and also tried at its
        old path.

Returns:
    The run directory as a Path, which exists.

Raises:
    FileNotFoundError: If that directory is absent. The message names the
        places the run DOES exist, so a run filed elsewhere reads as "it
        is over there" rather than as "it was never made".
```

**`run_dir_for_index_row()`**

```text
The run directory named by one run_index.csv row.

The index carries `kind`, `tag`, `start_year`, `end_year` and
`source_sink_preset` for exactly this: a row and a directory can be
matched without either side reconstructing the other's spelling. A row
from before 2026-09-16 carries `arm` instead, which is translated.

Args:
    raw_runs: The output/raw_runs root.
    row: A mapping or pandas Series with run_name, start_year, end_year,
        source_sink_preset, and kind/tag (or the legacy arm).

Returns:
    The run directory as a Path, resolved through find_run_dir so an
    unmigrated tree still answers.
```

**`arm_component()`**

```text
LEGACY. The path component a pre-09-16 arm contributed.

Kept only so an old caller importing it still imports. New code files by
kind and tag; see preset_dir_for.
```

**`_kind_tag_from_path()`**

```text
(kind, tag) read off a run directory's place in the tree.

Either layout: a purpose folder says so directly; an old arms/<arm>/ path
or a token-named cell is translated as legacy_arm_to_kind_tag would.
```

**`_legacy_arm_from_path()`**

```text
The `arm` a pre-09-16 index row spelled for a run at this path.

"calibration" for the unscoped tree (matrix runs AND token-named sweep
cells, either old layout); the arm folder(s) for arms/<arm>/ or a loose
pre-09-10 arm; and for the purpose layout, the arm LEGACY_ARMS maps the
(kind, tag) back to. This is what finds an old row for a moved run.
```

**`run_status()`**

```text
current / superseded / archived, for the index's `status` column.

A matrix or sensitivity run on a topography that is no longer the
product's CURRENT is superseded: comparing it with a run made today would
attribute the pick difference to whatever the figure is about. A version
or experiment run is judged against nothing -- it names its own inputs
deliberately -- and an archived run says so by where it sits.
```

**`rebuild_run_index()`**

```text
Regenerates run_index.csv from every run's metadata; returns the rows.

Args:
    raw_runs: The output/raw_runs root.
    index_path: Where to write; default raw_runs/run_index.csv.
    current_versions: {topo_product: CURRENT dune-topo version}, for the
        status column. None leaves status "current" for everything but
        the archive; pass hat_topo_version's answer to get supersession.

Returns:
    The list of row dicts written, in file order.

Notes:
    A run whose metadata has no "index row" section (made before
    2026-09-16) keeps the row the existing file holds for it, found by the
    old (run_name, Hs_m, arm) key and given kind/tag from its path. A run
    with neither is indexed by its identity columns alone, so it is not
    lost. Rows whose run is gone from disk are dropped here; recording
    them is HAT_index_runs.py's job (retired_runs.csv), which calls this.
```

**`load_run_index()`**

```text
run_index.csv as a DataFrame with kind/tag/status guaranteed present.

A file from before 2026-09-16 has `arm` instead; it is translated so a
reader written against the new columns works on either.
```

**`run_dir_contents()`**

```text
Lists what a run directory holds that counts as output.

The single definition of "this directory holds a result", so the guard
below and anything reporting on it cannot disagree about whether a stray
.gitkeep means the directory is occupied.

Args:
    run_dir: Directory to inspect. Need not exist.

Returns:
    Sorted list of Paths, empty if the directory is absent or holds
    nothing but ignorable files.
```

**`guard_run_dir()`**

```text
Refuses to write into a run directory that already holds output.

Called before the model is stepped, not after, so a name collision costs
nothing. Overwriting is what makes a paired comparison quietly wrong: the
baseline the groin run measures against is resolved by directory name, so
replacing one run's outputs silently redefines the other run's answer.
That is why the default refuses rather than asks.

With overwrite=True the directory is EMPTIED, not written over in place.
Writing over in place leaves behind any file the previous run produced and
the new one does not -- a difference GIF from a run that had a paired
baseline, a road summary from a run whose roadway manager was on -- and a
leftover file in a run directory is indistinguishable from a current one.
Emptying first means the directory holds exactly one run's output.

Args:
    run_dir: Directory this run will write to.
    overwrite: True to empty the directory and reuse it. Intended for
        iterating on one scenario -- tweak a value, re-run, read the
        figures, tweak again -- where only the current state matters. The
        previous run's output is deleted and is NOT recoverable.

Returns:
    The run_dir as a Path, created if it did not exist.

Raises:
    RuntimeError: If the directory holds output and overwrite is False,
        or if overwrite is True and the directory holds a subdirectory
        that is not one of the run's own (run_layout.SUBFOLDERS plus
        gif_frames). An unexpected subdirectory means this is not the
        directory it is taken to be, and deleting its contents is
        refused rather than guessed at.
```

**`render_metadata_text()`**

```text
Renders the metadata sections as the human-readable .txt.

Args:
    sections: Ordered mapping of section name to a mapping of key to
        either a value, or a (value, comment) tuple. The comment is shown
        in the .txt and dropped from the .json.
    header: Comment lines placed at the top of the file, without the
        leading "# ".

Returns:
    The file contents as a string.
```

**`write_run_metadata()`**

```text
Writes the run metadata as both .txt and .json.

Args:
    run_dir: Directory to write into.
    run_name: Run name, used for the filenames.
    sections: As accepted by render_metadata_text.
    header: Comment lines for the top of the .txt.

Returns:
    A (txt_path, json_path) tuple of Paths.
```

**`skill_vs_target()`**

```text
Model-minus-observed skill over the real domains.

Reported over two spans, deliberately. The end domains carry the locked
source/sink values (tens of m/yr), so an island-wide RMSE for a calibBE or
edgeBE run is dominated by two domains that were pinned rather than
predicted -- and a zeroBE run has no such term. Comparing presets on the
island-wide number alone would mostly compare the boundary treatment.
The interior number excludes them, which is the same exclusion the
source/sink QC plot makes for its zoom panel.

Args:
    change_rate: Padded per-domain model rate array, m/yr, already
        sign-flipped so (+) is seaward.
    target_table: DataFrame with gis_domain and target_lrr_m_yr columns
        (COASTSAT_TARGET).
    geometry: DomainGeometry describing the padded array.
    interior_margin: Domains excluded from each end for the interior
        metrics. 1 drops the two locked end domains.
    interior_gis: (first, last) GIS bounds for the interior instead of
        a margin, so an extended geometry (2026-09-16) is scored on the
        same GIS 2-89 as its 90-domain baseline. Overrides the margin.

Returns:
    A dict of mean bias and RMSE in m/yr, island-wide and interior, plus
    the domain count each was computed over. NaN where no domains remain.
```

**`append_run_index()`**

```text
Adds one run to the cross-run index CSV, replacing any earlier row.

Replaces rather than appends on a repeat run name so the index tracks what
is currently on disk. A run directory holds exactly one result, so two
index rows for one name could only ever mean one of them is stale.

New columns are unioned in, so adding a field later does not invalidate an
index written before it existed -- older rows get NaN for it.

Args:
    index_path: Path to the index CSV. Created if absent.
    row: Mapping of column name to value for this run.
    key: Column, or sequence of columns, identifying a run uniquely.

Returns:
    The full index as a DataFrame, as written.

Note:
    The index is a DERIVED view of the runs on disk: every value in it is
    restated from a run's own metadata. This keeps it in step as runs are
    made; `HAT_index_runs.py` rebuilds it from the runs, checks it against
    them, and records anything whose run has been deleted.
```

</details>

### shoreline.py

Shoreline time-series extraction from a completed CASCADE run.

From the script's original header:

```text
Shoreline time-series extraction from a completed CASCADE run.

Pure data extraction -- nothing here touches matplotlib. Feeds both the
shoreline GIF and the rate-comparison figures in cascade_pipeline.plotting.
```

Notes that were in the code:

```text
Evenly spaced in calendar years, so the slope is per year rather than
per state. The two agree for annual states; being explicit means a
sub-annual save spacing would not silently rescale the rate.
```

<details><summary>Function notes (the original docstrings)</summary>

**`get_x_s_TS()`**

```text
Extract the shoreline position time series from a Barrier3D object.

Args:
    b3d: A single-domain Barrier3D instance (element of cascade.barrier3d).

Returns:
    1-D float array, one value per model year, in the raw x_s_TS
    convention (dam, increasing landward with erosion).

Raises:
    AttributeError: Neither x_s_TS nor _x_s_TS is present.
```

**`build_shoreline_matrix()`**

```text
Build a [time x domain] shoreline matrix from a completed CASCADE run.

Args:
    cascade: A Cascade instance after cascade.update() has finished.
    to_meters: Convert dam -> m (CASCADE's native unit is dam).

Returns:
    2-D float array, shape (n_years, n_domains), raw x_s_TS convention
    (increases landward / with erosion). No sign flip is applied here;
    see compute_change_rate for the plotting convention used elsewhere
    in cascade_pipeline.
```

**`compute_change_rate()`**

```text
Modeled shoreline change rate (m/yr) from a [time x domain] matrix.

change_rate = (shoreline_m[-1, :] - shoreline_m[0, :]) / span_years

Args:
    shoreline_m: Array from build_shoreline_matrix(cascade, to_meters=True).
    span_years: Denominator for the rate, e.g. END_YEAR - START_YEAR.
        Defaults to (n_years - 1) if not given.
    flip_sign: CASCADE's x_s_TS increases landward (erosion). True
        (default) flips the sign so a positive rate means
        seaward/accreting -- the convention used throughout
        cascade_pipeline's plotting functions. Pass False to keep the raw
        x_s_TS sign convention.

Returns:
    1-D float array, length n_domains, m/yr.
```

**`compute_lrr()`**

```text
Modeled shoreline linear regression rate (m/yr), with fit quality.

The ordinary-least-squares slope through EVERY annual state, which is the
estimator the observational target is defined by: CoastSat's
transect_lrr_full.csv holds a per-transect OLS slope fit through the full
set of satellite shoreline positions in the period (`lrr_m_yr`, with
`r_squared` and `unc_m_yr` beside it). compute_change_rate's endpoint
difference is a different quantity -- a net displacement divided by a
span -- so plotting one against the other compares two estimators rather
than a model against an observation.

The two also differ in what they are sensitive to. An endpoint difference
reads only years 0 and N, so any single-year excursion still present in
the final state lands in the result at full amplitude. That is not
hypothetical here: a nourishment fill enters as an instantaneous step in
x_s, and BRIE's Crank-Nicolson alongshore solve answers a step with a
grid-scale (2*dy) mode whose amplification factor is (1-4r)/(1+4r) for
r = D*dt/(2*dy^2). At Hatteras r ~= 1.05, so that factor is about -0.61:
the mode alternates sign along the coast and decays only ~39% per year.
A fill late in a run therefore leaves a sawtooth in the final state, and
an endpoint rate reports all of it. The OLS slope spreads the same
excursion over every year and reports about a quarter of it.

r_squared is returned rather than left implicit because a slope is only a
summary of a trend where a trend is what the domain has. A nourished
domain sits flat for most of the period and then steps, which fits a line
badly (r_squared near 0) however the slope is computed -- that is a fact
about the trajectory worth carrying next to the number, not a defect in
the estimator. Note also that r_squared is a variance ratio, so it goes
small wherever the total signal is small, quite apart from linearity.

Args:
    shoreline_m: Array from build_shoreline_matrix(cascade, to_meters=True),
        shape (n_states, n_domains).
    span_years: Elapsed years the states span, e.g. END_YEAR - START_YEAR.
        Used only to space the regressor, so the slope comes out per
        calendar year. Defaults to (n_states - 1), which is the same
        thing whenever the states are annual.
    flip_sign: CASCADE's x_s_TS increases landward (erosion). True
        (default) flips the sign so a positive rate means
        seaward/accreting -- matching compute_change_rate and the
        convention used throughout cascade_pipeline's plotting.

Returns:
    An (lrr, r_squared) tuple of 1-D float arrays, length n_domains. lrr
    is in m/yr; r_squared is NaN for a domain that never moved, where the
    ratio is 0/0 rather than a perfect fit.
```

</details>
