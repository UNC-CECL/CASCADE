# Hatteras hindcast -- methods

The reasoning behind each section of `HAT_hindcast_1984_2024.ipynb`.
It was written in the notebook itself, which pushed the code that
actually runs a long way apart; it is the same text, moved, not
rewritten. The notebook keeps a short statement of what each section
does and links here for why.

### 2.1 Project paths

Everything derives from `SCRIPTS_DIR` (Section 1), so no absolute paths are
hardcoded. Layout under `data/hatteras_init/1-barrier3d-domains/`:

```
2009-dune-topo/<version>/topography/domain_<gis>_topography_<year>.npy
2009-dune-topo/<version>/dunes/domain_<gis>_dune_<year>.npy
2009-buffer/sample_1_topography.npy      # one profile, reused by every buffer
2009-buffer/sample_1_dune.npy
```

`TOPO_DUNE_VERSION` is **resolved, not pinned**: `topo_dirs()` reads `VERSION`
out of `HAT_dune_topo_extractor.py`, so the runner, the dune-start road
setbacks and the audits always describe the same extraction. Bump `VERSION` in
the extractor and everything follows; pass `topo_dirs(override="2009_v3")` to
reproduce an older run on purpose.

It has to work that way. The setback is measured from interior row 0, so a
runner pinned to one extraction while the setbacks are built on another places
the road against a row that does not exist on the grid being run. That is
exactly what happened from 2026-08-19, when the setbacks moved to `2009_v4` and
then `2009_v5` while this pin stayed on `2009_v3`.

`2009_v1` and `2009_v2` are in `2009-dune-topo/incorrect/`; `topo_dirs()` raises
if the resolved version is not on disk, before the check in 2.3 runs.

### 2.4 Units check against Barrier3D's input contract

Barrier3D's `load_input.py` converts the *scalar* YAML parameters but loads the
elevation and dune `.npy` files **verbatim** — `InteriorDomain` is whatever the
file contains, with no unit conversion applied:

```python
params["MHW"] /= 10.0                                  # m   -> dam
params["BarrierLength"] = int(params["BarrierLength"] / 10.0)   # m -> cells
params["BermEl"] = params["BermEl"] / 10.0 - params["MHW"]      # m -> dam above MHW
params["InteriorDomain"] = load_elevation(...)         # <- no conversion
```

So the arrays on disk must *already* be in decameters relative to MHW. A file
written in metres would run without error and silently model an island 10x too
tall. This step verifies that, plus the shape contract `load_input.py` relies
on, across all 90 real domains — deliberately on the **raw** arrays, before any
`DAM_TO_M` scaling.

What `load_input.py` assumes:

- `InteriorDomain.shape[1]` is alongshore. If it exceeds `BarrierLength` the
  array is silently truncated; if it is smaller, `BarrierLength` is silently
  reduced to match. Either way the mismatch is never reported.
- `DuneStart` is sliced `[0:BarrierLength]`, so the dune array must be at
  least that long.
- Dune values are heights *above the berm*, not elevations
  (`newDuneHeight = newDuneElev - BermEl`), so they must be positive.

## 3. Island orientation — set `START_YEAR`

`START_YEAR` is the one flip in this notebook. It selects a period from
`HATTERAS_PERIODS` (in `hatteras_site_config.py`), and everything in Section 4
follows from it: run length, RSLR rate, storm series, background-erosion
preset, plus the road-setback and nourishment settings used later.

The island offset is what "orientation" means here — one value per padded
domain giving that domain's cross-shore starting position. It becomes
`shoreline_offset` on the `Cascade()` call, so it sets where the shoreline
starts, and a wrong `START_YEAR` shows up as the wrong island shape.

### `SCENARIO` — the one place a run is defined

The cell below is where you choose what to simulate. `SCENARIO` names a
management combination and expands to the four `ENABLE_*` switches that
sections 5 and 6 read:

| `SCENARIO` | roadway | beach/dune | fills | what it isolates |
|---|---|---|---|---|
| `natural` | off | off | off | the counterfactual |
| `roadway_only` | **on** | off | off | the NC-12 program alone |
| `beachdune_only` | off | **on** | **on** | the beach program alone |
| `full_management` | **on** | **on** | **on** | the status-quo hindcast |
| `full_no_fill` | **on** | **on** | off | the fill, against `full_management` |

`GROIN_ENABLED` sits in the same cell but outside the table, because §12.3
measures the groin against a paired no-groin baseline identical in every other
token. Every scenario is run twice: `False` first to write the baseline, then
`True`.

Named scenarios are the intended route, but the commented override lines below
them let you depart for a single run. A departure is detected and printed, and
the run name follows the **switches**, never the label — so an overridden run
cannot be filed under the scenario it departed from.

Two combinations cannot both hold, and are resolved in the cell rather than
left to fail later: fills need a `beach_dune_manager` to spend them, and a
relocation event needs a `roadway_manager` to read the setback it moves.

The cell prints the run name it will produce. That preview is advisory — §7.5
derives the real name from the modules sections 5 and 6 actually built and
raises if the two disagree, so a switch that fails to reach its module is an
error rather than a mislabelled directory.

### 4.3 Source/sink (background erosion)

Background erosion is CASCADE's stand-in for shoreline change driven by
alongshore transport gradients that the model does not resolve. It is a per
domain rate in m/yr, passed to Barrier3D as `Rat`. Sign convention, from
`cascade/brie_coupler.py`: **(-) = erosion, (+) = accretion**.

The presets live in `hatteras_site_config.py` and are sparse -- a domain absent
from a preset gets 0.0 m/yr. Each is one hypothesis about where the sediment
budget is unresolved, and they are the source/sink axis of the run matrix:

- **`zeroBE`** imposes nothing anywhere. Whatever the shoreline does is what
  Barrier3D and BRIE produce unaided.
- **`edgeBE`** pins only the two end domains (GIS 1 and 90), which absorb the
  open-boundary artifact at the ends of the modelled reach. Nothing is imposed
  on the interior, so an interior misfit is the model's own. Its values are
  *sliced* from `calibBE`, not retyped, so the two cannot disagree about the
  ends.
- **`calibBE`** is the per-domain fit against the CoastSat LRR target rates in
  Section 8.

All three are plotted so the selected preset is visible against the
alternatives. The old names `base` and `calibrated` still work as deprecated
aliases for `zeroBE` and `calibBE`; Section 3 normalises them to the canonical
key, so a run directory is always named for the canonical preset. Note that
`base` was **all zeros**, with the edge values sitting beside them as
comments -- it was `zeroBE` all along, despite the name.

## 5. `roadway_manager` — setbacks, per-domain elevation, historical events

NC-12's forcing is three things: **where** the road sits (setback, by period),
**how high** it is (elevation, period-independent), and **which domains are
managed** at all.

The setback is measured landward of the digitised dune line of its own year,
then spent against a grid whose row 0 is the 2009 dune crest. Where the island
narrowed between those dates the road can land behind the barrier, so §5.1
audits every setback against the interior *before* the run and names the
domains whose road drowns in year one.

Two conventions worth stating, both Barrier3D's:

- `road_ele` is **metres MHW-relative**, not NAVD88 — `bulldoze` writes it
  straight into the interior grid, which the extractor stores MHW-relative.
- relocation events carry a **displacement**, not an absolute setback. CASCADE
  already decrements the setback by dune migration each year, so adding the
  measured displacement counts the retreat once; an absolute setback referenced
  to an older dune line counts it twice.

Whether any of this is *acted on* is `ENABLE_ROADWAY_MANAGEMENT`, in section 3.
Off, the road stays fully loaded, audited, and plotted, but no `RoadwayManager`
ever bulldozes, rebuilds a dune, or decrements a setback. The road becomes a
record of where NC-12 was, not a thing the model defends.

## 6. `beach_dune_manager` — nourishment schedule + overwash filter

Two different things arrive bundled in one CASCADE module, and the distinction
matters more than the name suggests.

**Always-on**, in every domain where the module is enabled, for every year of
the run: a percentage of overwash deposition is removed from the interior and
returned to the shoreface, the dune line is held fixed in the cross-shore
(`dune_migration_on = False`), and the community is abandoned if the average
interior width falls below 50 m.

**Event-driven**, only where and when `nourish_now` says so: sand is added to
the shoreface, moving the shoreline seaward and widening the beach.

There is no way to get the second without the first. Turning the module on to
deliver a 2022 fill also turns on overwash filtering and a fixed dune line in
that domain from year one.

Three conventions, all of them places this pipeline has been wrong before:

- `overwash_filter` is a **percent**, not a fraction. `filter_overwash` divides
  it by 100, and Rogers et al. (2015) give 40–90 %, residential to commercial.
  A value of `0.4` filters 0.4 % of overwash, which is indistinguishable from
  no filtering. `BeachDuneConfig` now refuses the fraction scale.
- The per-year volume goes to **`cascade.nourishment_volume`**, which is where
  CASCADE reads it from. This is not a way around the module: `beach_dune_manager`
  still performs every fill, unmodified. Stock `Cascade.update()` copies its own
  `nourishment_volume` list into each manager on every step and then calls the
  manager, which spends it (`cascade_groin.py:781-789` ->
  `beach_dune_manager.py:788`). So a volume written onto the `BeachDuneManager`
  instance instead lands on the same attribute CASCADE overwrites one line before
  the manager reads it, and the fill quietly spends the `Cascade()` init default.
  `apply_to_cascade()` writes the two Cascade-level inputs -- `nourish_now` and
  `nourishment_volume` -- and nothing else.
- The manager's time series are **offset by one**. `update_dune_domain()`
  increments `time_index` before the management modules run, so a manager
  writing at `time_index - 1` lands on index `year - start_year + 1`.

### The footprint is a union of two different maps

`beach_dune_manager` runs on the permanent community zones **and** every domain
that receives fill. Those come from unrelated sources — settlement extents
versus engineering project extents — and the projects are wider than the
villages. Buxton 2022 runs north out of the village into the NC-12 corridor,
and the Rodanthe 2014 emergency fill sits north of Tri-Village entirely.

Where that union overlaps the roadway footprint, **both managers run on the
same domain in the same year**. CASCADE permits this — the two loops in
`cascade.update()` are independent — but it is outside what Anarde et al.
published, and it has two consequences the cell below names explicitly:
overwash is removed twice over (bulldozed first, then the survivors filtered),
and the fixed dune line leaves `ShorelineChangeTS` at 0, which is the value
`RoadwayManager` reads to decrement the road setback. NC-12 stops retreating
in those domains.

That is reported, not corrected — §12 asserts it against the finished run
rather than trusting the prediction.

### Two switches, because the module is two things

Section 3 sets both. `ENABLE_BEACH_DUNE_MANAGEMENT` decides whether the module
runs at all: off, the footprint is empty everywhere, `overwash_filter` is handed
to CASCADE as zeros, and the village domains get natural shoreline behaviour --
the dune line migrates, overwash lands where it lands, no width check.

`ENABLE_NOURISHMENT_FILLS` sits inside that: with the module on, it decides only
whether the historical fill is spent. It works by driving an **empty schedule**
for the same period, not by shrinking the footprint, so `BEACH_DUNE_MANAGEMENT_ON`
is identical either way. That is what makes the fills-on / fills-off pair a
single-variable comparison -- the always-on behaviour described above is held
fixed, and the difference between the two runs is the sand.

The distinction matters because of the paragraph above: turning fills off does
**not** give you natural dynamics in the villages. Only turning the module off
does.

## 7. `hard_structures` / groin -- Buxton groin field

`cascade/groin.py`'s `CorrectedGroinCallback` attaches to a run through
`cascade._groin_callback` and is called once per model year from inside
`Cascade.update()`, immediately before the alongshore-transport solve
(`cascade/cascade_groin.py`). It reads BRIE's shoreline and wave-climate tables,
adds to `x_s_dt`, and changes nothing inside BRIE or Barrier3D. Each active year:

1. **It holds back a share of the drift.** BRIE's net alongshore drift for the
   wave climate, read from its own tables at the shoreline angle where the drift
   arrives at the GIS 5|6 face, is converted to shoreline change
   (`q dt / (dy (h_b_crit + d_sf))`); a share `s` of it is moved from the
   downdrift domain (GIS 5) to the updrift one (GIS 6).
2. **It stops BRIE leaking that sand away.** BRIE's diffusive step would smooth
   the step the groin builds and spread the held sand onto the rest of the
   island. Within 10 cells of the face the callback adds back
   `2 (r_i - r_{i-1}) (x_i - x_{i-1})` for this shoreline minus the same term for
   the groin-free shoreline, so the correction acts on the groin's own bends and
   never on the coast's natural ones. The groin-free shoreline is this shoreline
   less the groin's tracked footprint: each year BRIE's alongshore step for the
   previous year is repeated with and without the groin's own contribution, and
   the difference is the footprint. A groin of strength 0 is therefore exactly
   the no-groin run, which a full run confirms to 0 m.

Units and sign are BRIE's: `x_s_dt` reaches the hook in metres, and `x_s`
increases landward, so the updrift advance is a negative `x_s_dt`.

### Two numbers, fitted over the groin's whole life

The groin is set by its **strength**, the share of the arriving drift held
back (`groin.strength`, 0.25), and the **first year it no longer works**
(`groin.end_year`, 2017). Both were fitted on one continuous run from the 1967
wet/dry shoreline, two years before the groins, to 2025, scored on the gap
across the groins against the wet/dry photo dates and annual CoastSat
(`hard-structures/groin/3-hindcast/2-blocking-1996-2025/2026-10-13-pregroin-planform/`).
Fitted era by era, the strength is 0.25 for 1969-1995, 0.25 for 1996-2009 and
0.20-0.25 for 2010-2016 (0.25 is within 0.5 m), so one strength is kept. It is
the minimum under three different stand-ins for the pre-1984 storms. The end
year is the 2017 Buxton fill, placed across the whole groin field: left on, the
groin misses 2017-2025 by 47 m.

**Why not fit the calibration and test windows separately.** Each window's run
starts from the observed shoreline, which already holds the fillet the groins
built. A strength fitted that way measures how fast the model relaxes that
fillet (0.16 on 1996-2009, 0.10 on 2010-2016), not how strong the groin is: run
as one history, those values leave the gap 30-40 m below the record from 1996
to 2017.

**Why the no-groin runs move at GIS 5|6.** The 1996 and 2009 start offsets step
240-264 m seaward from GIS 5 to 6, against 104 m on the 1967 dune line: about
100 m is the coast's bend toward the cape and about 150 m the fillet. Without a
groin BRIE smooths it, GIS 6 eroding and GIS 5 building out with almost no net
change; it is not a source or a sink.

### What is not represented

- **The build-up speed.** At full strength from installation the model builds
  its fillet in about 3 years; the photos show about 9. Ramping the strength in
  over 1969-1978 reproduces it without moving the plateau, which matters only
  for a run that starts before the fillet exists (a forecast of a rebuilt groin).
- **The fillet's width.** The model's change sits in the two domains beside the
  face; the observed one spans three or four, and is mostly downdrift retreat.
- **One groin against three.** 0.25 is the three-groin field's strength; only
  the southernmost groin is rebuilt in 2026.
- **The 2017 fill's footprint.** The fill starts at GIS 6, so it makes a gap
  across the groins of its own (35 m in the no-groin run).

### History

The dipole (a fixed `-M`/`+M` a year, with a deterioration floor `f`) and the
blocking groin (a share `b` of BRIE's diffusive exchange cancelled at the face)
were the runner's groins until 2026-10-14; their fits and why they were replaced
are in `hard-structures/groin/README.md`, and git keeps their code (last at
commit `4d937b79`).

## 8. CoastSat target rates -- LOWESS windows

This section builds the observational target. The `calibBE` source/sink
preset in section 4.3 was fit against the curve produced here, and section 12
draws the model against it -- so what this section decides is what the model is
judged by. It emits `cs_series` for the section 12 figures and
`COASTSAT_TARGET`, a per-domain table of the target rate itself.

### Smooth at transect resolution, then aggregate

`transect_lrr_full.csv` holds one linear-regression rate per CoastSat transect:
906 transects over GIS 1-90 in each period, roughly ten per domain. LOWESS runs
at **transect** resolution using along-coast distance as x, and only then
averages to domain resolution. The order matters -- averaging to domains first
would discard the within-domain structure the smoother exists to use.

Two windows are computed: 7 domains (~3.5 km) and 10 domains (~5.0 km).

### "Primary reference" is a convention, not a setting

The widest window is the target. That is *implicit*: `LowessConfig` has no field
naming a reference window, and `rate_comparison` picks `max(window_domains)`.
The old run script had an explicit `RESIDUALS_LOWESS_WINDOW = 10`, which agreed
only because 10 happened to be the larger of the two. Set `window_domains` to
`(7, 12)` and the reference silently becomes 12, with nothing raising. The cell
below names `TARGET_WINDOW` explicitly so the choice is at least visible at the
point of use.

### Domains 1-10 carry no smoothed line

`skip_southern_domains = 10` suppresses LOWESS across GIS 1-10: boundary effects
near the Cape dominate that zone and smoothing flattens the sharp gradient
there. Those domains are shown as raw transect scatter only.

Worth being explicit, because an earlier config comment in the run script says
otherwise: it describes raw per-domain means being *stitched* onto the widest
LOWESS line so the curve stays continuous. That never ran. Both the old script
and `splice_lowess_with_raw_south` drop domains `1..skip_n` for every window;
the function's `transect_domain_ids`, `transect_lrr_values` and
`is_widest_window` arguments are accepted and ignored.

So the target has two different provenances along its length, and
`COASTSAT_TARGET` labels each row accordingly: **raw per-domain mean** for GIS
1-10, **LOWESS 7-domain** for GIS 11-90 (10-domain until 2026-09-28).

That boundary matters for section 7. The Buxton groin sits at GIS 5.5, so both
of its flanking domains fall in the unsmoothed zone -- the groin's
observational target is a raw mean over a handful of transects, not a point on
the reference curve.

### Suppression is display-only, and it leaks north

`skip_southern_domains` drops GIS 1-10 from the plotted line. It does **not**
remove them from the fit: LOWESS runs over all 906 transects, and only the
result is truncated. So the southern transects the skip exists to discount
still pull the smoothed values immediately north of the cut.

Refitting with GIS 1-10 excluded outright, against the target as built:

| domain | target as built | refit without D1-10 | difference |
|---|---|---|---|
| D11 | -2.10 | -2.78 | **+0.68** |
| D12 | -1.84 | -2.18 | +0.34 |
| D13 | -1.43 | -1.60 | +0.17 |
| D15 | -0.44 | -0.53 | +0.08 |
| D16 and north | | | < 0.01 |

The leak decays within about half a window and is negligible by GIS 16, but
0.68 m/yr at GIS 11 is not nothing on a target spanning -4.8 to +1.6 -- and
GIS 9 is where NC-12 management begins. Left as built, because the calibrated
source/sink preset was fit against this curve and changing the fit domain would
invalidate it. Recorded so the choice is deliberate rather than inherited.

### The along-coast axis is reconstructed, not measured

`transect_lrr_full.csv` carries no along-coast coordinate, so
`load_transect_data` builds one: each domain's transects are spread evenly
across its 500 m band, `(domain - 1) * 500 + (rank + 0.5) * (500 / n)`. Real
positions do exist -- every one of these 906 transects appears in
`groin_analysis_chainage_all.csv` with a measured `alongshore_m` -- and they
are not quite the same. GIS 7's four transects really span 327 m of their 500 m
band, and there is a 124 m gap between the last GIS 6 transect and the first
GIS 7 one.

The reconstruction is monotonic and correct to within a few tens of meters per
domain, so the smoothed curve is close. It is kept because the alternative
makes this pipeline depend on a groin-analysis output; the fix belongs upstream
in `input_prep`, not here. Recorded so the approximation is not mistaken for
measurement.

### The groin differential

The cell prints observed updrift-minus-downdrift LRR for both periods, since
that is the observational check on the dipole direction section 7 imposes. A
positive differential means the updrift domain is doing better than the
downdrift one, which is what a functioning groin should produce.

## 9. Figure configuration

The plotting *functions* are already imported in section 1 -- `make_shoreline_gif`,
`make_all_shoreline_gifs`, `plot_rate_comparison`, `plot_annotated_rate_comparison`,
`build_shoreline_matrix`, `compute_change_rate`. This section is not those. It is
the configuration they need, gathered before section 12 calls them.

### Why it is a section and not four keyword arguments in section 12

Every plotting entry point defaults to the **package** defaults, not Hatteras':

```python
plot_rate_comparison(..., domains=DEFAULT_DOMAINS, annotations=DEFAULT_ANNOTATIONS,
                     lowess_config=DEFAULT_LOWESS, config=DEFAULT_RATE_COMPARISON)
make_all_shoreline_gifs(..., domains=DEFAULT_DOMAINS, annotations=DEFAULT_ANNOTATIONS,
                        gif_config=DEFAULT_GIF_CONFIG)
```

`DEFAULT_ANNOTATIONS` is an `AnnotationConfig()` with every field empty --
`town_spans {}`, `village_lines {}`, `piers {}`, `groins {}`, `shoal_zones {}`.
Forget to pass `annotations=HATTERAS_ANNOTATIONS` at any one of those call sites
and that figure comes out with no villages, no piers, no groin line and no shoal
zones. Nothing raises; the figure just quietly loses its geography.

`lowess_config` is the same shape of hazard with a subtler failure. `DEFAULT_LOWESS`
happens to equal what section 8 configures today, so omitting it currently
changes nothing -- but section 8 owns that setting, and the moment it changes,
section 12 would keep drawing the package default and the figure would no longer
show the target section 8 reported.

`domains=DEFAULT_DOMAINS` is safe, and only by coincidence: `domains.py` notes
its field defaults happen to match this hindcast (90 real, 15 buffer, first GIS
1, 500 m).

So the two bundles below, `RATE_FIG_KWARGS` and `GIF_KWARGS`, exist to be
splatted into every call rather than retyped at each one. Two rather than one
because the GIF and rate-comparison functions take different keywords.

### Difference GIFs resolve their own baseline

A `mode="difference"` job needs a previous run's shoreline matrix to difference
against, and `make_all_shoreline_gifs` silently skips the job when
`baseline_npy` is `None`. Since section 7.5 derives run names from the active
switches, the matching no-groin run's name is computable: same scenario, `groin`
token swapped for `nogroin`. The cell looks for that run's saved matrix and uses
it if it exists, reporting either way -- so a difference GIF is either produced
or explained, never skipped in silence.

This only works because `gif_config.save_matrix` is `True`: each run writes
`<run_name>_shoreline_matrix.npy` into its own output directory, which is what
makes it available as a baseline later.

### Provenance

The plotting code was extracted from the run script into
`cascade_pipeline.plotting` -- `shoreline_gif.py` and `rate_comparison.py`, with
the shared geographic layer in `cascade_pipeline.annotations`. The roughly twenty
module-level constants became explicit dataclass configs (`GifConfig`,
`AnnotationConfig`, `RateComparisonConfig`) plus a `RunInfo` bundling
run_name/run_dir/period/Hs/sign convention, passed in rather than read off
globals -- which is what makes them importable here independent of the run
script's state. Public function names are unchanged, so this is a drop-in for any
other script that already calls them.

### 9.4 Validation target -- where the island actually ended up

A `"position"` or `"displacement"` GIF already carries a dashed grey year-0
reference. This adds the other end of the comparison: a static line showing the
**surveyed** island position in the run's end year, so how close the run gets is
readable off the animation instead of only off the rate figure.

#### Why this cannot be built from the padded offset files

The obvious construction -- `Island_Dune_Offsets_2004` minus
`Island_Dune_Offsets_1984` -- is wrong, and quietly so. `island_offset_hybrid.py`
computes each year independently as

```python
baseline         = min(mean_distances)      # :125
relative_offsets = np.subtract(mean_distances, baseline)
```

so **each year is zeroed on its own most-seaward domain**. That domain (GIS 78,
in both years) itself retreated 56.4 m between 1984 and 2004, so differencing the
two padded files subtracts a 56.4 m constant that has nothing to do with any
domain's motion. It turns a real mean retreat of +15.1 m into an apparent 41.3 m
of accretion -- the wrong sign for Hatteras.

The fix is to go back to the raw transect CSVs, whose `ORIG_LEN` is measured from
a **single fixed offshore datum** (`offshore_datum_line.geojson`, a straight N-S
line at easting 460198.45 shared by every year). Absolute distances differenced
across years are a true shoreline change; min-subtracted ones are not.

Sign follows `ORIG_LEN` being a distance *from an offshore line*: larger = farther
landward, which is the same direction as `x_s_TS`, so the change adds straight
onto the model's year-0 position with no flip.

#### What the line is, and what it is not

It is the **dune line**, the same feature that set `shoreline_offset` -- not
CoastSat's wet/dry shoreline, which section 8 already targets separately. Using
the dune line keeps the target in the same frame as the initial condition, so the
year-0 base cancels and the gap between the two lines is the model's misfit.

Buffer domains get `NaN`: their offsets are extrapolated ramps, not surveys, so
there is nothing to validate against and the line simply stops at the real span.

Only produced when a raw file exists for the run's end year. The 2004-2024 period
has no 2024 survey, so it reports "no target" and every figure falls back to the
model-only form.

## 10. Define `run_cascade_simulation()`

Two functions. `build_cascade` constructs the model and attaches the groin;
`run_cascade_simulation` steps it through the period, applying historical
management, and writes the run's artifacts. Both take every input explicitly
rather than closing over notebook state, so the sweep scripts can call them
too.

The split is what lets section 11 hold a built-but-unstepped `Cascade`.
BRIE's diffusivity and the groin's fillet prediction are only meaningful as
initial conditions, and a prediction printed after the run is not one.

### `run_years` is transitions, not states

The original signature took `nt` and looped `range(nt - 1)`, and `main()`
passed `nt = END_YEAR - START_YEAR`. That combination is off by one, and it
cost a model year on every run to date:

- Barrier3D seeds `_x_s_TS = [x_s]` at init (`barrier3d.py:1289`) and appends
  one entry per update, so N updates produce N+1 annual states spanning N
  years.
- `nt = 20` therefore ran 19 updates -> 20 states -> **1984 to 2003**, while
  the rate was divided by `END_YEAR - START_YEAR` = 20. Rates came out low by a
  factor of 19/20, about 5 %.
- The storm file carries time indices 1-21 and Barrier3D selects with
  `StormSeries[:,0] == time_index`, which only ever reached 19. Storm years 20
  and 21 were never applied.

This function takes **`run_years`**, the number of annual transitions, and
derives the rest: `time_step_count = run_years + 1`, and the loop runs exactly
`run_years` updates. Passing `RUN_YEARS` now means what it says, and the
denominator section 12 uses matches the span actually simulated.

That fits the preallocation exactly. `brie_coupler.py:117` writes
`TMAX = time_step_count`, so arrays hold `run_years + 1` slots; after
`run_years` updates `time_index` is `run_years + 1` and the last write lands at
`time_index - 1 = run_years`, the final valid index. One more update would
overflow.

**Both periods therefore gain a model year relative to every previous run.**
1984-2004 and 2004-2024 each now simulate 20 transitions instead of 19. The
2004 state is shared between the two periods as a boundary condition -- it is
the end state of one run and the observed initial condition of the next -- but
no transition is simulated twice.

### Nourishment is written where CASCADE reads it

`beach_dune_manager` does the nourishing, unmodified. Nothing in CASCADE is
patched. The only question is which attribute the per-year volume is written
to, and stock `Cascade.update()` answers it:

```python
# cascade_groin.py:781-789, unmodified
self._nourishments[iB3D].nourishment_volume = self._nourishment_volume[iB3D]
[self._nourish_now[iB3D], self._rebuild_dune_now[iB3D]] =     self._nourishments[iB3D].update(...)
```

CASCADE copies its own `nourishment_volume` list into each manager on every
step, then calls the manager, which spends it at `beach_dune_manager.py:788`
(`self._nourishment_volume / 100`, m^3/m -> dam^3/dam). So
`cascade.nourishment_volume` is the **designed input**, not a workaround, and
`BN_SCHEDULE.apply_to_cascade(cascade, year)` writes exactly there.

The original loop wrote the volume onto the manager instead:

```python
nourishment_obj._nourishment_volume = float(nourish_vol[iB3D])   # clobbered
```

That is the same private attribute CASCADE overwrites one line before the
manager reads it, so the value never survived to be spent.

The failure was narrower than "no nourishment happened". `nourish_now` was set
at the Cascade level, which is correct, so fills fired on the right domains in
the right years -- they just spent the `Cascade()` init default. For 2004-2024
that default was 100 m^3/m, against actual project volumes of 619.3 (Rodanthe
2014), 183.5 (Buxton 2022) and 841.0 (Avon 2022) m^3/m.

The init volume is now `[0.0] * total_domains` from section 6, deliberately:
if a scheduled fill ever fails to reach the model again, that year nourishes
nothing and section 12's verification catches it, rather than a plausible-
looking default hiding the failure. `nourishment_interval` stays `None`, so
`beach_dune_manager.py:779` (`if nourish_now or self._nourishment_counter == 0`)
can only be triggered by the historical schedule -- no automatic
renourishment enters a hindcast.

Section 12 verifies delivery by reading each manager's own
`_nourishment_volume_TS`, the volume the model actually spent, and asserts on
a mismatch.

### Historical roadway events live in the package

The event logic moved to `cascade_pipeline.roadway.apply_historical_event`,
beside the setback loading and auditing it belongs with. It returns rows; this
cell prints them. Three behaviours documented there and worth repeating:

- A relocation is a **displacement** added to the setback the model currently
  holds, which already carries the modelled retreat. An absolute setback
  double-counts it -- that is what put NC-12 into the sound in earlier
  versions.
- `road_ele` is deliberately not reset. It is initialised from the 2004
  alignment, which is already the post-relocation road, and the manager
  decrements it in the Lagrangian frame, so rebuilding at grade would be an
  exact no-op.
- CASCADE's own relocation guards are probed and **reported, then overridden**.
  These relocations happened; refusing them would be historically wrong. Every
  refusal is printed.

A `BridgeEvent` now actually switches management off, by flipping
`cascade.roadway_management_module`. The original recorded the pads in a
`terminated_road_pads` set that nothing ever read, so GIS 82-88 kept full road
management for the rest of the run after the 2022 Jug Handle Bridge. It still
does not null `cascade.roadways[pad]` -- `Cascade.update()` reads
`.drown_break` on every domain with no None check.

### Artifacts

The function is handed `run_dir` rather than deriving it, and writes the saved
model, the nourishment log and the groin diagnostics there. Section 12 does
figures only.

## 11. Initialize Cascade -- single config, no sweep

Section 10 defined the procedure; this section supplies the values sections
2-9 do not already produce and builds the model. Section 12 steps it.

The split exists so a built-but-unstepped `Cascade` can be inspected. Two of
the numbers below are only meaningful as initial conditions:

- **BRIE's diffusion number.** `r_ipl = coast_diff * dt / (2 * dy^2)` is a
  local variable inside BRIE's `update()` (`brie.py:1294`), recomputed every
  step from the local shoreline angle and never stored. `_coast_diff` itself is
  set once at init, so the cell below evaluates `r_ipl` at shore-normal
  incidence from the freshly built model. Reading it after the run would give
  the end-of-run geometry.
- **The fillet prediction.** `predict_fillet` turns the groin's `M` into a
  predicted amplitude and alongshore extent. Printing it here, before any time
  step, is what makes the extent an actual prediction rather than something
  fitted after the fact. Section 12 checks the emergent extent against it.

### Wave height sets the shoreface

`Hs = 2.5 m` is the calibration value, and it does more than force storms:
BRIE derives the shoreface depth from it directly, `d_sf = 8.9 * Hs`
(`brie.py:270`), giving 22.25 m. That is exactly the `DShoreface: 22.25` in
`Hatteras-CASCADE-parameters.yaml`, which confirms that value is **metres** --
one of the two units section 7 had to leave open. The cell resolves section 7's
implied-interception volume with it.

### Dune thresholds are floored by CASCADE, so both numbers are printed

`dune_design_elevation` and `dune_minimum_elevation` are documented as m MHW
(`cascade_groin.py:236-238`). `roadway_manager.py:723-729` raises them on the
first step:

```python
dune_design_elevation  = max(passed, BermEl * 10 + 1.0)
dune_minimum_elevation = max(passed, BermEl * 10 + 0.3)
```

The original run script passed `0.01` for the minimum, commented `# dam` --
the wrong unit for the argument, and far below the floor, so CASCADE replaced
it every time. The value below is stated in the documented unit instead, and
the cell prints both what was passed and the floor CASCADE will apply, so the
threshold that actually governs is visible rather than implicit.

That printed floor is also a live check on an unresolved units question.
`brie_coupler.py:183` writes `BermEl` as `[m NAVD88]`, Barrier3D pops it
unconverted, and then uses it in both decametre contexts
(`Dmax = Dmaxel - BermEl  # (dam)`) and metre contexts
(`MaxAvgSlope = BermEl / 10`). If `BermEl * 10 + 0.3` prints as roughly 2 m the
metre reading holds; if it prints near 17 m it does not. Flagged, not chased.

## 12. Run the loop, verify, then figures

Steps the model section 11 built, checks what it actually did, and only then
draws anything.

### Verification comes before figures, deliberately

Every check below runs against the finished model, and the nourishment check
**asserts**. A fill that never reached the model invalidates the run, and that
exact failure went unnoticed for a long time precisely because it produced
plausible-looking output -- so it stops the notebook rather than quietly
feeding a figure.

Nothing is lost when it fires. `run_cascade_simulation` saves the model and its
logs before returning, so a failed assert costs the figures, not the run.

Four checks, each testing something a previous section predicted rather than
assumed:

- **Run length.** `run_years + 1` annual states, or the run ended early.
- **Nourishment delivery** (`verify_nourishment`), read from each manager's own
  `_nourishment_volume_TS` -- the volume the model spent, not the volume the
  schedule intended. Section 6 explains why only the manager's record can tell
  those apart.
- **Frozen setbacks** (`verify_setbacks_frozen`) in the double-managed domains
  section 6 identified. `beach_dune_manager` holds `dune_migration_on` False,
  Barrier3D then leaves `ShorelineChangeTS` at 0, and `roadway_manager` only
  decrements the setback when that is non-zero -- so NC-12 should stop
  retreating there. Predicted in section 6, confirmed here.
- **Roadway survival** (`summarise_road_management`), from CASCADE's own
  `drown_break` / `relocation_break` state.

### Two estimators, and why the figures draw the LRR

The run writes two columns for the same thing, and they are not
interchangeable:

| column | definition | reads |
|---|---|---|
| `change_rate_m_yr` | `(x[-1] - x[0]) / RUN_YEARS` | 2 of the 21 states |
| `lrr_m_yr` | OLS slope through every state | all 21 |

The observational target is an LRR. `transect_lrr_full.csv` holds a
per-transect ordinary-least-squares slope fit through the full set of CoastSat
shoreline positions in the period -- `lrr_m_yr`, with `r_squared`, `p_value`
and `unc_m_yr` beside it -- and section 8 LOWESS-smooths those slopes
alongshore. Plotting an endpoint difference against it compares a net
displacement to a trend, so `RATE_ESTIMATOR = "lrr"` is what the figures, the
skill metrics, and the `calibBE` fit all use. `SKILL_ENDPOINT` is reported
next to `SKILL` only because `calibBE` and the groin `M` were fit before this
distinction was drawn.

**Expect the LRR skill to be slightly worse.** Modelled erosion decelerates
across a period, so an all-years slope is more erosive than an endpoint
difference: interior RMSE rises in most runs (zeroBE 2004-2024
`road_bdm_nourish`: 2.048 to 2.136 m/yr). That is the estimator changing, not
the model getting worse, and the old number is kept in the metadata so the two
can be told apart.

#### What the endpoint difference was actually reporting

CASCADE applies a nourishment fill as an **instantaneous step** in `x_s`
(`shoreface_nourishment`: `x_s - 2V / (2 h_b + D_sf)`; the 841 m^3/m Avon
project moves the shoreline 68 m seaward in one year). BRIE answers a step
with ringing. Its alongshore solve is Crank-Nicolson --
`brie.py` builds `mat = I + 2r`, `rhs = x + r d^2x/dy^2 + dx_dt` with
`r = D dt / (2 dy^2)` -- which is unconditionally stable but **not monotone**.
At Hatteras `D ~= 5.2e5 m^2/yr`, `dy = 500 m`, `dt = 1 yr`, so `r ~= 1.05` and
the grid-scale (2*dy) mode has amplification `(1 - 4r) / (1 + 4r) ~= -0.61`:
it alternates sign along the coast and decays only ~39% per year.

Measured on the runs, from the 2014 Rodanthe fill, the year-on-year decay of
the grid-scale roughness is 0.49, 0.66, 0.55, 0.65, 0.58, 0.66, 0.60, 0.69 --
the predicted 0.61. That fill had ten years to damp out and is invisible. The
2022 Buxton and Avon fills have **two**, so 35 m of ringing is still in the
2024 state, and an endpoint rate divides it by 20 and reports a +-1.7 m/yr
sawtooth through GIS 6-15 and 22-28. The LRR spreads the same excursion over
every year and reports about a quarter of it; what survives is the real fill
edge. Grid-scale roughness of the plotted profile, `max |d2 rate / dy2|`:

| runs | endpoint | LRR |
|---|---|---|
| 1984-2004 (no fills) | 0.66-0.73 | 0.19-0.21 |
| 2004-2024, fills off | 0.87-1.07 | 0.17-0.32 |
| 2004-2024, fills on | 1.76-2.18 | 0.40-0.51 |

This is a reporting artefact, not a broken module. `beach_dune_manager` is
applied exactly as section 6 describes, `nourishment_ok` is True, and a
`bdm`-on/fills-off run is as smooth as a `nobdm` run.

#### `lrr_r2` is written, and read

A slope only summarises a trend where the domain has one. A nourished domain
sits flat for eighteen years and then steps, which fits a line badly however
the slope is computed -- GIS 23-26 come out at `r2 ~= 0.00`. That is a fact
about the trajectory, not a defect in the estimator, so it ships as a column
and the run prints how many domains fall below `LRR_R2_FLOOR`. Note `r2` is a
variance ratio, so a domain that barely moved lands there too; the two cases
are told apart by whether the domain took a fill.

For scale on the other side: the observed transect LRRs have a median
`r_squared` of 0.20 and a median `unc_m_yr` of 0.24. The target curve is a
weak trend through noisy positions. That is normal for CoastSat, and it is
the reason both sides need the same estimator rather than one side using a
cleaner-looking two-point difference.

### The rate denominator is asserted, not assumed

`compute_change_rate` divides by `span_years`. Section 10's `run_years` is the
transition count, so `shoreline_m` has `run_years + 1` rows and the correct
denominator is `nt - 1 == RUN_YEARS`. That equality is checked rather than
trusted, because the version of this pipeline that divided a 19-year change by
20 looked exactly as healthy as this one.

### The groin extent check

Section 11 printed a predicted fillet extent **before** the run. Here it is
measured. Extent needs the groin's effect isolated, which needs the paired
no-groin run: this run's final shoreline minus the baseline's, with the span
where that difference exceeds `GROIN_EXTENT_THRESHOLD_FRAC` of its own peak
taken as the extent, updrift and downdrift separately.

The threshold is a choice, and it is stated in the output rather than buried.
Amplitude was tuned; extent was not, so extent is the part that can actually
fail. When no baseline exists on disk the check is skipped with a message
saying how to produce one -- run once with `GROIN_ENABLED = False`.

### Re-running this cell

The cell refuses to step a model that has already been stepped. Running it
twice would walk past `TMAX` and raise an `IndexError` from inside Barrier3D,
which is an unhelpful way to discover you needed to re-run section 11.

### 12.4 Write

The change-rate CSV, the road-management summary, the run metadata, and one
row appended to the cross-run index.

Metadata is written twice from a single structure: a `.txt` to read and a
`.json` to parse. Building both from one mapping is what keeps them from
disagreeing -- the earlier version emitted only prose, so anything comparing
runs had to re-parse it.

`[identity]` records the three things that drift silently across a batch run
over several days: `USE_SANDBOX_CASCADE` (which `Cascade` class ran),
`TOPO_DUNE_VERSION` (which extractor built the init surface), and the git
commit with a dirty flag. A dirty tree is not an error -- most runs happen
mid-edit -- but it does mean the commit alone will not reproduce the run, so
the flag is recorded rather than inferred later.

`run_index.csv`, at the root of `output/raw_runs/`, carries one row
per run: the switches, the skill metrics, and whether the run completed. It is
keyed on `run_name`, so re-running a scenario replaces its row instead of
adding a second -- a run directory holds one result, so two rows for one name
could only mean one is stale.

Skill is reported over two spans. The end domains carry the locked source/sink
values (tens of m/yr), so an island-wide RMSE for a `calibBE` or `edgeBE` run
is dominated by two domains that were pinned rather than predicted, while a
`zeroBE` run has no such term. Ranking the presets on the island-wide number
alone would mostly rank their boundary treatment, so the interior number
(GIS 2-89) is carried beside it.
