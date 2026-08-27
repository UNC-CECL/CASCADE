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

`cascade/groin.py` attaches to a run through `cascade._groin_callback` and is
called once per model year from inside `Cascade.update()`, immediately before
the alongshore-transport solve (`cascade/cascade_groin.py:600`). Each active
year it adds `-M` to the updrift domain and `+M` to the downdrift domain of
`x_s_dt`, then hands the array on. BRIE's implicit diffusion solve spreads that
dipole in the same step, so the fillet's taper and alongshore extent are
**emergent** -- only its amplitude is imposed.

Units and sign are the ones BRIE already uses, verified rather than assumed:
`brie_coupler.py:57` converts `x_s_dt` to meters before the hook sees it and
line 344 converts back, so `M` is meters per year with no conversion in the
module. `x_s` increases landward, so `-M` updrift is a seaward advance.

### It is a forcing, not a barrier

Nothing in this module blocks alongshore transport. Littoral drift crosses GIS
5.5 in the model exactly as it would with no structure there. What the dipole
reproduces is the *shoreline signature* of blockage -- accretion one side,
starvation the other -- by moving shoreline displacement between two adjacent
cells. The pair is volume-neutral by construction: `-M` and `+M` cancel.

Three consequences, stated here rather than discovered later:

- **Trapping never saturates on state.** A real groin fills to its tip and
  bypasses; trapping shuts off because the fillet grew. Here `M` is applied at
  full strength every year, modulated only by a calendar schedule. The
  amplitude does level off -- diffusion balances injection at
  `A ~= M / (4 * r_ipl)` -- but for a reason unrelated to the physical one.
- **The year counter is the call counter.** `groin.py:230` computes
  `year = start_year + _call_count`. The callback never reads `cascade` and
  never checks `time_index`, so it assumes it is called exactly once per model
  year, in order. It cannot detect the off-by-one that section 6 documents for
  `BeachDuneManager`.
- **`install_year` is inert here.** 1969 precedes both 1984 and 2004, so
  `active` is `True` at every step of either period and `active_TS` carries no
  information. It mattered for the 1967 test runs.

One dipole also stands in for the whole Buxton groin field. At 500 m domains
the field is sub-domain scale, so that is a resolution limit rather than a
modelling choice.

### Two independent estimates of `M` disagree, and that is reported

`M` is the single tunable knob, and two lines of evidence give answers an order
of magnitude apart.

The sensitivity sweep fit `M` against a **shoreline-position** target and found
`M = 50 m/yr` (RMSE 24.0, `HAT_groin_sweep_results.csv`) -- from the 1967-2017
test setup, not this period.

The **sediment budget** says that is unaffordable. Moving `M` meters of
shoreline over a 500 m domain transfers `M * dy * (d_sf + h_b)` per year. At
`M = 50` that is 3-6 x 10^5 m3/yr depending on the active profile height,
against a reach-integrated transport-gradient loss of 5.9 x 10^5 m3/yr for the
whole of Oregon Inlet to Cape Hatteras (Inman & Dolan 1989, via Moore et al.
2010, doi:10.1029/2009JF001299). The groin would be moving something close to
the island's entire annual budget across one domain boundary, every year, while
deteriorating from 1996 and storm-damaged in 2003.

That comparison is order-of-magnitude only: 5.9 x 10^5 m3/yr is a *divergence*
over a ~60 km reach, not a gross flux past Buxton. The cell prints that caveat
with the number on every run.

The run is configured at the sweep's `M = 50` and the breach is **reported, not
corrected** -- the same treatment section 6 gives the double-managed domains.
If both estimates are right, the extra shoreline signal at D5/D6 is not groin
trapping; it is the Cape Point shoal dynamics the calibrated source/sink rates
already claim there (`+1.2` at D5, `+2.0` at D6, both accretion). Those rates
were fit to observed CoastSat LRR spanning the functional-groin era, so the
groin's signature is plausibly counted twice. Left standing pending the
source/sink re-analysis; section 12 checks it against the finished run.

### What the fitted differential actually measures

`differential = LRR[D6] - LRR[D5]` is algebraically the **OLS slope of the
fillet**, `x_s[D5] - x_s[D6]`. Its trend across the window, not its size. Three
consequences, all counter-intuitive, all measured on the 2026-08-22 sweep:

- The fillet **saturates** -- ~18 m at `M = 40`, reached by 1990, where BRIE's
  diffusion balances the trapping. A healthy, strong, *constant* groin therefore
  scores near **zero**.
- A **degrading** groin scores **negative**, because its fillet is relaxing
  during the window.
- Sign reports building-vs-failing, not strong-vs-weak. Two very different
  groins can score identically.

| year | fillet, `f = 0.00` | fillet, `f = 1.00` | applied `M` |
|---|---|---|---|
| 1984 | 0.0 | 0.0 | 40.0 |
| 1990 | 17.4 | 17.4 | 40.0 |
| 1996 | 17.6 | 17.6 | 40.0 |
| 2000 | 10.6 | 18.3 | 17.1 |
| 2004 | 0.7 | 18.3 | 0.0 |

#### The "period 2 is unreachable" claim was wrong

Earlier versions of `HAT_groin_sweep_config.py` argued that a negative
differential could not be produced at any `M >= 0`, since the callback adds `-M`
updrift and `+M` downdrift. That reasoning holds for the fillet's *size*, not
its *slope*. Falsified twice:

1. Period-1 cell `M40_beNA_f0.00` scores **-0.713**, more negative than the
   no-groin baseline's **-0.229** -- a groin at `M = 40` driving the
   differential *below* no groin at all.
2. The period-2 seed run inherits a **25.3 m fillet at t = 0** (the real Buxton
   fillet is in the 2004 initial shoreline) and relaxes to 21.9 m, scoring
   **-0.168 m/yr** at `M = 50, f = 0.9`. Negative, with a groin attached.

So the observed 2004-2024 differential of **-2.47 m/yr** is reachable, and is
what a fillet relaxing after the 2003 storm damage looks like.

#### Why the fit ranks on period 1 alone

Cumulative trapping is `M*(16 + 4f)` in 1984-2004 but `20*M*f` in 2004-2024.
Period 2 sits entirely past the ramp, so `M` and `f` enter **only as their
product** -- that window cannot separate them. Every bit of information
distinguishing the two lives in period 1, the one window straddling the
1996-2003 ramp.

Summing both legs let a window that cannot identify the parameters pull them
anyway, and in the wrong direction: holding period-1 trapping fixed, driving `f`
from 0.9 to 0 forces `M` from 50 to **61.3 m/yr** -- 52% to 63% of the reach
sediment budget. Period 2 is now reported out-of-sample, where its failure to
match is informative rather than guaranteed.

### Deterioration is scheduled, not simulated

Last repair 1996, storm damage 2003, encoded as a linear ramp between them
(`deterioration_delay_years = 27`, `ramp_years = 7`). For a 1984 start the ramp
falls at run-years 12-19; for a 2004 start it has already completed, so the
groin sits at its floor from year 0. The floor itself, `0.9`, is the sweep's
paired value for `M = 50` -- it leaves the structure at 90 % strength after
failure, which the documented 1996/2003 history does not support. Printed as an
open tension rather than silently reconciled.

### What the 1984 sweep found

The 1984 `zeroBE` sweep, 43 cells, is the only leg run to completion
(2026-08-23). `zeroBE` imposes **no** source/sink anywhere, and `edgeBE`
imposes it only at GIS 1 and 90, so neither preset puts anything at D5/D6 --
the groin is fitted against the whole unexplained signal there. `calibBE` is
excluded from `PRESETS` for the opposite reason: fitting a groin against a
field that already absorbed the groin would count it twice.

**On the fillet's slope, the fit has no optimum.** Error falls monotonically in
`f` at every `M`, and monotonically in `M` along `f = 1.0`. The best cell is
the top-right corner, `M = 110, f = 1.0`, railed on both axes:

| M \ f | 0.0 | 0.4 | 0.8 | 1.0 |
|---|---|---|---|---|
| 40 | -0.713 | -0.418 | -0.121 | +0.026 |
| 70 | -1.086 | -0.565 | -0.029 | +0.249 |
| 110 | -1.616 | -0.773 | +0.102 | **+0.555** |

Observed is **+0.986**, the no-groin baseline **-0.229**. The `f = 1.0`
contribution is linear in `M` to R^2 = 0.998, so the grid ceiling delivers 65 %
of what is needed and the line reaches the target only at `M ~= 168 m/yr` --
**175 % of the reach sediment budget**, moving 1.13x the entire Buxton 2022
nourishment across one domain boundary every year. Matching the observed
differential would need a groin that never degrades *and* is unaffordable, on
two independent counts.

**On the fillet's size, the same 43 runs give an interior optimum.** Rescored
with `measure_fillet` against the detrended observation (13.80 m):

| ranked on | best cell | error | railed | mean applied `M` | budget |
|---|---|---|---|---|---|
| slope | M=110, f=1.0 | 0.431 m/yr | both axes | 110.0 | 115 % |
| **size** | **M=50, f=0.6** | **0.068 m** | neither | 46.0 | 48 % |

The ridge is well formed -- M=50/f=0.6 at 13.87 m, M=70/f=0.4 at 13.46 m,
M=40/f=0.8 at 14.77 m -- and `f = 0.6` is consistent with the 1996 repair and
2003 storm the schedule already encodes, where `f = 1.0` flatly is not. This is
the number to quote, with the caveat below.

### The caveat that limits the fit: the observed shape is not a dipole

Detrending the observed LRR across D1-D12 with the groin pair excluded gives
the pair's local anomaly. It is **one-signed**:

```
observed anomaly:   D5 +0.786    D6 +1.476     sum +2.26
model anomaly:      D5 -1.543    D6 -0.759     sum -2.30
```

Both observed domains sit *above* the regional trend. The groin dipole is
volume-neutral by construction, so it must place a **notch** at D5 -- and the
model does, at every `(M, f)` on the grid. The model therefore cannot reproduce
the observed *shape* at any parameter, regardless of amplitude.

D5's `+0.786` is also within the scatter of the regional fit (D1 `+0.869`, D3
`-0.812`, D8 `+0.774`); only D6's `+1.476` is a genuine outlier. What the data
shows is an accretional bulge centred near Cape Point, not a groin dipole. A
well-conditioned fit to a feature whose shape the model cannot produce is a
better-behaved number, not a truer one -- so report the fillet-size `(M, f)`
*and* this mismatch together.

### Bypassing is missing, and it matters at high `M`

`GroinCallback` is a pure forcing: "trapping never saturates on state". A real
groin stops trapping once its fillet reaches the seaward tip, which also stops
the downdrift starvation, because sand resumes rounding the end. Without that,
the source saturates on diffusion while the sink saturates on nothing:

| M | downdrift D5 | updrift D6 | fillet |
|---|---|---|---|
| 40 | +10.5 m erosion | -7.8 m accretion | 18 m |
| 70 | +25.2 m | -7.6 m | 33 m |
| 110 | **+52.2 m** | **-0.8 m** | 53 m |

52 m of notch against 0.8 m of advance is not a groin; it is a runaway sink.
**The differential fit rails straight into that regime**, which is a second
reason not to quote its `M = 110`.

`bypass_fillet_m` was added to `GroinCallback` (default `None`, so every prior
run reproduces) as a **capacity limit** -- the year's trapping is capped at the
fillet the structure can still accept, `cap - fillet`. A hard switch and a
linear taper were both tried first and both degenerate at an annual timestep,
because one year at `M = 70` moves the fillet ~140 m and overshoots any
plausible cap; their shared signature is mean applied `M` pinned to exactly
`M/2`. Verified at M=70, f=1.0:

| cap | mean applied `M` |
|---|---|
| 15 m | 18.60 |
| 25 m | 25.44 |
| 40 m | 35.82 |
| 100 m | 70.00 (never binds; identical to uncapped) |

Capping at a plausible fillet cuts effective trapping ~4x at high `M`. **The
full capped sweep was not completed** -- three attempts were lost to the machine
suspending mid-sweep, which kills in-flight workers on `WORKER_TIMEOUT_S`
wall-clock -- so the cap is a demonstrated mechanism and a documented model gap,
not a fitted parameter. The cap value itself is unsourced: the 1969 Navy groins
have no published length in anything this project cites.

### What is deliberately not here

No `predict_fillet` call. Its `r_ipl` is a local variable inside BRIE's
`update()` (`brie.py:1294`), computed per domain per step from `_coast_diff`
and never stored -- and `finalize()` deletes `_coast_diff` outright. Section 11
reads it off the constructed model and prints the predicted amplitude and
extent **before** the time loop, so the emergent-extent check is pre-registered
rather than retrofitted.

## 8. CoastSat target rates -- LOESS windows

This section builds the observational target. The `calibBE` source/sink
preset in section 4.3 was fit against the curve produced here, and section 12
draws the model against it -- so what this section decides is what the model is
judged by. It emits `cs_series` for the section 12 figures and
`COASTSAT_TARGET`, a per-domain table of the target rate itself.

### Smooth at transect resolution, then aggregate

`transect_lrr_full.csv` holds one linear-regression rate per CoastSat transect:
906 transects over GIS 1-90 in each period, roughly ten per domain. LOESS runs
at **transect** resolution using along-coast distance as x, and only then
averages to domain resolution. The order matters -- averaging to domains first
would discard the within-domain structure the smoother exists to use.

Two windows are computed: 7 domains (~3.5 km) and 10 domains (~5.0 km).

### "Primary reference" is a convention, not a setting

The widest window is the target. That is *implicit*: `LoessConfig` has no field
naming a reference window, and `rate_comparison` picks `max(window_domains)`.
The old run script had an explicit `RESIDUALS_LOESS_WINDOW = 10`, which agreed
only because 10 happened to be the larger of the two. Set `window_domains` to
`(7, 12)` and the reference silently becomes 12, with nothing raising. The cell
below names `TARGET_WINDOW` explicitly so the choice is at least visible at the
point of use.

### Domains 1-10 carry no smoothed line

`skip_southern_domains = 10` suppresses LOESS across GIS 1-10: boundary effects
near the Cape dominate that zone and smoothing flattens the sharp gradient
there. Those domains are shown as raw transect scatter only.

Worth being explicit, because an earlier config comment in the run script says
otherwise: it describes raw per-domain means being *stitched* onto the widest
LOESS line so the curve stays continuous. That never ran. Both the old script
and `splice_loess_with_raw_south` drop domains `1..skip_n` for every window;
the function's `transect_domain_ids`, `transect_lrr_values` and
`is_widest_window` arguments are accepted and ignored.

So the target has two different provenances along its length, and
`COASTSAT_TARGET` labels each row accordingly: **raw per-domain mean** for GIS
1-10, **LOESS 10-domain** for GIS 11-90.

That boundary matters for section 7. The Buxton groin sits at GIS 5.5, so both
of its flanking domains fall in the unsmoothed zone -- the groin's
observational target is a raw mean over a handful of transects, not a point on
the reference curve.

### Suppression is display-only, and it leaks north

`skip_southern_domains` drops GIS 1-10 from the plotted line. It does **not**
remove them from the fit: LOESS runs over all 906 transects, and only the
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
                     loess_config=DEFAULT_LOESS, config=DEFAULT_RATE_COMPARISON)
make_all_shoreline_gifs(..., domains=DEFAULT_DOMAINS, annotations=DEFAULT_ANNOTATIONS,
                        gif_config=DEFAULT_GIF_CONFIG)
```

`DEFAULT_ANNOTATIONS` is an `AnnotationConfig()` with every field empty --
`town_spans {}`, `village_lines {}`, `piers {}`, `groins {}`, `shoal_zones {}`.
Forget to pass `annotations=HATTERAS_ANNOTATIONS` at any one of those call sites
and that figure comes out with no villages, no piers, no groin line and no shoal
zones. Nothing raises; the figure just quietly loses its geography.

`loess_config` is the same shape of hazard with a subtler failure. `DEFAULT_LOESS`
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
and `unc_m_yr` beside it -- and section 8 LOESS-smooths those slopes
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
