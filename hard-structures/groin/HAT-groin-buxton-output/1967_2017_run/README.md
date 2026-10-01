# 1967_2017_run - the 1967-2017 groin rig

The reduced-reach groin rig: GIS D2-D12 (11 real domains inside 41 padded),
1967 to 2017, the Buxton groin as a dipole at D5/D6 with the 1995-2003
deterioration ramp, and the 1971 and 1973 nourishments. The (M, f) fit, its
sweep and the edge correction are run from here; the fitted pair the
production hindcast uses is recorded in `output/calibration/groin/`, and
`GROIN_PLAN.md` (one folder up from `HAT-groin-buxton-output/`) is the
authority on it.

| script | what it does | writes |
|---|---|---|
| `HAT_groin_hindcast_1967_2017.py` | the rig runner: every key in `RUN_MATRIX`; imported by the sweep, its worker and the edge solve | `output/calibration/groin_rig/<run>/` |
| `HAT_groin_sensitivity_sweep.py` | coarse then fine grid over M and the deterioration fraction, one process per cell, resumable | `../../HAT-buxton-hindcast-groin-test/sensitivity_sweep/` |
| `HAT_groin_sweep_single_combo.py` | one sweep cell, launched by the sweep | a profile `.npy` in `sensitivity_sweep/profiles/`, one stdout line |
| `HAT_solve_edge_be.py` | solves the D2 / D12 edge background-erosion rates, groin off | `edge_be_solved.json` beside it |
| `HAT_groin_trajectory_target.py` | the observed fillet trajectory 1967-2023, and the scorer against it | prints only; a library |
| `HAT_groin_threeway_hindcast_1967_2017.py` | baseline / nourishment only / nourishment + groin in one go | `output/raw_runs/<run>/` |
| `HAT_plot_groin_runs.py` | figures of saved rig runs; the runners import its `fig_` functions | `PLOT_*.png` in the first run's folder |
| `HAT_check_nourishment_applied.py` | did the 1971/1973 fills reach x_s inside CASCADE? | prints only |
| `comparison/HAT_groin_effect_comparison.py` | no-groin vs groin against the observed shoreline at checkpoint years | see `comparison/README.md` |

All observed targets come from
`../shoreline_position_output/Change_from_wetdry_1967_D2_D12.csv`, made by
`HAT-groin-buxton-input/input_prep/shoreline_position/HAT_geometric_distance_sanity_check.py`.

The two runners descend from `HAT_groin_hindcast_1967_1997.py`, the original
30-year (1967-1997) groin-only test. That script was deleted on 2026-10-01
with the rest of `1967_1997_run/`; recover it from git (see
`../README.md`).

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning, the
choices behind it and its history are here, one section per script. Moved out
of the scripts on 2026-10-01, when they were brought in line with
`scripts/STYLE.md`. Passages from the scripts' own headers are kept word for
word under "From the script's original header".

### HAT_groin_hindcast_1967_2017.py

Extended groin-deterioration-test hindcast (1967-2017, GIS D2-D12). Built
from `HAT_groin_hindcast_1967_1997.py` (deleted 2026-10-01; in git) -- same
inputs and conventions -- and extended from a 30-year (1967-1997) groin-only
test to a 51-year (1967-2017) test that also exercises the 1995->2003
deterioration ramp, which needs the longer window to actually play out (both
dates fall after 1997, so the original 30-year test never reached them).

    RUN_MATRIX = ["no_groin", "groin"]
      -> HAT_1967_2018_M60_deterioration_no_groin   (groin OFF -- erosive baseline)
      -> HAT_1967_2018_M60_deterioration_groin      (groin ON, WITH deterioration)

Run keys: `"no_groin"` (erosive baseline), `"groin"` (dipole at D5/D6, with
deterioration; needs the hook), `"groin_be"` (groin plus regional background
erosion, `REGIONAL_BE_RATE_M_YR`). Put one key in the list to run a single
experiment; list several to run them in sequence. Names match the comparison
folders and the plotter's `RUNS` list. (The run names actually come out as
`HAT_1967_2018_edge_calibrated_<key>`, from `RUN_NAME_SUFFIX`.)

**Which Cascade.** While testing the groin, the runner uses a sandbox copy of
`cascade.py`, `cascade/cascade_groin.py`, inside the package beside the real
one, with the 3-line groin hook added; the real `cascade.py` stays
untouched. Once the groin is proven, fold the hook into the real `cascade.py`
and set `USE_SANDBOX_CASCADE = False`. The `"no_groin"` run needs nothing
beyond the working base setup; the `"groin"` run needs `cascade.groin`
importable and the inert pre-AST hook (it sets `cascade._groin_callback`).
If the hook is missing, the groin run warns loudly (the diagnostics stay
empty), so a no-op is never mistaken for a real groin run.

**The groin is the shared module.** This used to import
`scripts.groin_module.hindcast_groin_test.version_control.HAT_groin_module`,
which no longer exists: the groin became part of the package as
`cascade/groin.py`. Importing the shared definition also means this 1967 fit
and the production hindcast fit the SAME model rather than two copies that
can drift apart.

**END_YEAR is exclusive.** `RUN_YEARS = END_YEAR - START_YEAR`, so the
original 1967-1997 run, with `END_YEAR = 1997`, only ever simulated through
1996. To reach 2017 as the final modelled year, `END_YEAR` is 2018 (51 model
years, 1967-2017 inclusive), matching the storm file built by
`HAT_build_1967_2017_storms.py`.

**Deterioration.** Last repair 1995 (after Hurricane Gordon damage in 1994)
-> linear decline -> Hurricane Isabel 2003 locks in the new deteriorated
state. The delay is relative to `GROIN_INSTALL_YEAR` (this script's own 1970,
not the true historical 1969), so the module resolves the right calendar
years whichever install-year convention the test uses.

The floor fraction was originally M/3 (~0.333, Katherine's "one groin still
functional out of three" framing; see the discussion with Laura and
Katherine on the Coastal Sediments 2027 abstract), but that removed too much
of the groin signal; 0.5 (50% of trapping retained) was then tested as a
less severe assumption, still folded into the ramp rather than run as a
separate experiment. It is **0.60**, the decided pair, 2026-08-30. It was
0.50, the 2026-08-24 sweep answer on the PRE-FIX topography. Re-run on
1984-start/v1, the rig's own sweep returns f = 0.6 (RMSE 23.78 against 24.21
at 0.5), bracketed on both sides -- and 0.6 is what production uses. The
sweep overrides this per cell; it matters only for a standalone run, which is
exactly what the full-life figure plots.

**Edge source/sink correction (Section 3b).** Mirrors the main 1984-2024
hindcast's edge correction at GIS 1 / GIS 90: the outermost REAL domains
(here D2 and D12) sit directly against buffer padding, and the buffer's
flat/repeated orientation can introduce an artificial alongshore-transport
signal right at that boundary. A background_erosion value (m/yr, same sign
convention as `REGIONAL_BE_RATE_M_YR`: negative = erosive) at the edge
domain(s) corrects for it. This is a STRUCTURAL fix (same role as GIS 1 / 90
in the main script), not a scientific choice like `REGIONAL_BE_RATE_M_YR` --
so it is applied to EVERY run in `RUN_MATRIX` (no_groin, groin, groin_be
alike), and it STACKS additively with `REGIONAL_BE_RATE_M_YR` on groin_be
runs rather than being overwritten by it. Set both edge values to 0.0 to
disable. The rates 5.0 and 10.0 are placeholders marked "(solve for this)";
`HAT_solve_edge_be.py` solves them.

**Historical beach nourishment (Section 3c: 1971, 1973).** Two documented NPS
nourishment projects fall within this 1967-2017 window, both targeting the
erosion embayment adjacent to the (Navy, 1969) groin field this script
models. Applied to EVERY run in `RUN_MATRIX` -- same treatment as the edge BE
correction -- since these are real historical events, not an experimental
variable.

Sources:

- Machemehl (1973): reports the 1971 volume as 300,000 cy, sourced from a
  man-made lake at Cape Point, pumped via 14-in cutterhead dredge (JA LaPort
  Dredging Co.) ~3.5 mi to a discharge point near the Hatteras Court Motel,
  left to migrate south under normal littoral drift.
- NPS (1980), pg 48: reports the 1971 volume as 200,000 cy (borrow material
  "proved insufficient to have any significant impact"); reports the 1973
  volume as 1,300,000 cy from an interior Cape Point borrow area (basin still
  visible today as altered vegetation), discharged via a 16-in dredge + 3
  boosters ~4 mi north near the Hatteras Court Motel, widening the beach
  ~500 ft over a cited 5,000-ft reach.
- Dolan, Hayden, Riddel & Ponton (1974), "1973 Buxton Beach Nourishment
  Project: An Annotated Photographic Atlas," NPS Contract No. CX5000031059
  (the primary, station-surveyed source behind NPS 1980's 1973 figures):
  independently confirms 200,000 cy for 1971 (agreeing with NPS 1980, not
  Machemehl); gives exact engineering stations for the 1973 project (south
  limit STA 2235+00, explicitly excluding the groin cells -- "no material was
  pumped into the groin cells, southern sediment drift caused a build-up" --
  north limit near STA 2164+80, MP41.0) and confirms the Navy's 1969 groin
  field is this script's modelled groin.

Volume choice: NPS (1980)'s 200,000 cy for 1971 rather than Machemehl's
300,000 cy -- two independent sources (NPS 1980, Dolan et al. 1974) agree on
200,000 cy against one source for 300,000 cy.

Domain range (derived, NOT yet confirmed against the project's own GIS domain
shapefile -- verify before treating these as final):

- 1973: templated engineered fill, ~D6-D10. Anchored on the atlas's own
  station data (south limit STA 2235 ~= 0.24 domains north of the groin; north
  limit STA 2164+80 ~= 4.5 domains north of the groin), converted at
  500 m/domain and 100 ft/station, using the script's GROIN_UPDRIFT /
  DOWNDRIFT_GIS (D5/D6) as the shared anchor point with the atlas's "Navy
  Groins" landmark. 1,300,000 / 5 x 0.764555 = 198,784.3 m^3 per domain.
- 1971: NOT a templated fill -- the atlas gives no station range for it, only
  that it targeted the same "north embayment" and was discharged as a point
  source near the motels, then left to migrate south under natural littoral
  drift. Modelled as a SINGLE-domain point injection at D8 (the motels, per
  the derived station crosswalk), deliberately NOT spread across a template --
  letting BRIE's own alongshore transport carry it south is the same "extent
  is emergent, not prescribed" logic used for the groin dipole, and matches
  what the historical record says happened. 200,000 cy x 0.764555 =
  152,911.0 m^3, all in D8.

D8 needs entries for BOTH years (the 1971 point injection AND its share of
the 1973 template), so its list is built by combining the two rather than
appearing twice as a dict key.

**Units reminder.** Every volume is stored and printed in m^3, NOT cubic
yards. 1 cy = 0.764555 m^3, so the m^3 figures look ~24% SMALLER than the cy
numbers quoted in the source reports (e.g. 1973's "1,300,000 cy" becomes
993,921.5 m^3 total). If a number looks low compared to a report, check units
before assuming an error -- it is very likely cy vs m^3, not a mistake.

**Where the fill volume goes.** In the time loop the runner sets Cascade's OWN
`nourishment_volume` list, NOT `cascade.nourishments[iB3D]._nourishment_volume`
directly. `cascade.update()` overwrites `nourishments[iB3D].nourishment_volume`
from `cascade._nourishment_volume[iB3D]` every call (see `cascade_groin.py`,
right before `nourishments[iB3D].update()` is called) -- setting the object's
attribute directly gets silently discarded the moment `update()` runs, which
is exactly why the requested volume never reached x_s.
`HAT_check_nourishment_applied.py` checks that it does.

`build_nourishment_arrays_from_manual_inputs()` builds per-year
nourishment-on and volume arrays for the CASCADE time loop from
`HAT_BN_YEARS` + `HAT_BN_VOLUME_BY_DOMAIN`, mirroring the main 1984-2024
hindcast's function of the same name exactly. Years outside
[START_YEAR, END_YEAR] are silently skipped, so
`ENABLE_HISTORICAL_NOURISHMENT = False` returns all-zero arrays with no other
change needed. It returns `nourishment_on_by_year` ({year: array of
TOTAL_DOMAINS}, 1/0) and `nourishment_volume_by_year` ({year: list of
TOTAL_DOMAINS}, m^3/m).

**Paths.** `PROJECT_BASE_DIR` is derived, not hardcoded. It was `r"/"`, which
silently resolved every data path to the filesystem root -- the scripts
imported fine and then reported every input as missing. Walking up to
`pyproject.toml` survives the reorganisation that moved this tree from
`scripts/groin/` to `hard-structures/`.

Runs go to `output/calibration/groin_rig/`, not `raw_runs/`, since
2026-08-31: the rig is a 41-domain grid and production is 120, M is
grid-specific, and the rig files no run_index row -- so its runs sat in
raw_runs unindexed, next to production runs they must not be compared with.
One of them was found holding an unstable M = 70 cell while named as though
it were the calibrated run.

**Topography.** The rig is a 1967-2017 window; `1984-start` is the nearer
product in time and the one the production period-1 groin fit reads, so the
two routes agree on a surface. Chosen 2026-08-30; previously this resolved to
2004-start by default. `TOPO_DUNE_INIT_YEAR` and `TOPO_DUNE_SUBFOLDER` are
legacy labels, no longer used to build filenames.

`build_file_lists()` is RESOLVED, NOT HARDCODED. It used
`HATTERAS_DATA_BASE/topography/2009/` and `/dunes/2009/`, a flat layout that
no longer exists: the domains moved under
`1-barrier3d-domains/2009-dune-topo/<version>/` and are now VERSIONED. Every
combination of the 2026-08-24 sweep died on the missing files.
`site_layer/hat_topo_version.py` is the project's single resolver for which
version is current -- pinning a version string here is what created this
breakage in the first place, so it is deliberately not pinned. It was
REPOINTED on 2026-08-30 for the period-first restructure of 2026-08-25.
Three things had gone stale and none of them errored:

1. `topo_dirs()` with no product resolves `DEFAULT_PRODUCT` ("2004-start").
   The rig is a 1967 window, so it should read the 1984 product -- and that is
   also what the production period-1 fit uses, so the two routes share a
   surface. This is the same omission that put the production groin sweep on
   the wrong island until 2026-08-30.
2. Array names lost their year suffix: `domain_N_topography_2009.npy` is now
   `domain_N_topography.npy`. `array_name()` owns that spelling.
3. `"2009-buffer"` became `"buffer"`; `BUFFER_DIR` owns that path.

**Figures.** `_save_run_figures()` generates and saves each run's figures into
its own folder, reusing the plotting functions from `HAT_plot_groin_runs.py`
(single source of truth -- no duplicated plot code): position change, rate,
trajectories, model vs observed, planform. It never blocks the run: any
plotting error is caught and reported. `_save_run_gif()` animates the
modelled shoreline over the run (real domains D2-D12), year by year, in
POSITION mode matching the main hindcast: real planform relative to the
year-0 alongshore mean, ocean at bottom (seaward downward). It shows the
fillet growing at the groin against the real island orientation.

Plot with `comparison/HAT_groin_effect_comparison.py`:

    RUN_NO_GROIN = "HAT_1967_2018_M60_deterioration_no_groin"
    RUN_GROIN    = "HAT_1967_2018_M60_deterioration_groin"

(no_groin first = baseline for the difference figure.)

**Cited by line number in archived configs.** The archived 7-source-sink
site configs cite `HAT_groin_hindcast_1967_2017.py:280`. The restyle of
2026-10-01 moved that line; it means `RIG_TOPO_PRODUCT`. The live READMEs
that cited `:76` and `:280` now name `_gis_to_pad()` and `RIG_TOPO_PRODUCT`
instead.

### HAT_groin_sensitivity_sweep.py

Grid search over `GROIN_TRAPPING_RATE_M_YR` (M) and
`GROIN_DETERIORATION_FRACTION` for the 1967-2017 groin hindcast, scored
against the observed wet/dry 2018 target (D2-D12, full range).

**Redesigned for crash-safety.** An earlier in-process version of this sweep
(running all 30 simulations back-to-back inside ONE Python process) crashed
partway through with a Windows access violation (0xC0000005) -- the signature
of accumulated state (Cascade/Barrier3D objects, joblib worker pools, or
similar) building up over many repeated simulations in a single long-lived
process, not a bug in any individual run. Fix: EVERY (M, fraction)
combination now runs in its OWN fresh subprocess
(`HAT_groin_sweep_single_combo.py`), guaranteeing the OS fully reclaims
memory and handles between every simulation, and results are written to CSV
IMMEDIATELY after each combo, not batched at the end -- so a crash costs at
most one combination's worth of time, not the whole sweep.

**Resumable.** If the script is interrupted (crash, closed terminal, etc.),
run it again. Already-successful combinations (found in the results CSV) are
skipped; failed ones (recorded with RMSE = NaN) are retried automatically. A
missing results CSV starts a typed empty frame: dtypes are declared
explicitly (not left to infer from an empty frame) -- concatenating a real
numeric row onto an empty, dtype-unspecified DataFrame is exactly what
triggers pandas' "concatenation with empty/all-NA entries" FutureWarning, and
silently produces object-dtype columns even on pandas versions that no longer
warn about it. Each result is appended and the CSV rewritten at once, so
progress is never lost even if the very next combination crashes.

A cell counts as failed on ANY failure -- non-zero exit code (including an
access violation), or a clean exit with unparseable output -- so the
orchestrator logs it and moves on rather than losing the whole sweep.

**Fit metric.** RMSE over the FULL D2-D12 range (both updrift and downdrift)
-- matching the observed shoreline position as closely as possible overall,
not isolating the groin's own signal. M and the fraction mainly move domains
near the groin (roughly D5-D12); the downdrift domains furthest from it
(D2-D4) are dominated by Cape Point dynamics the groin barely touches, so the
"best" combo found here reflects a balance between the groin-sensitive
domains and a residual error elsewhere that no groin parameter can close.

**Strategy** -- efficient by construction, not brute force. Stage 1: a COARSE
grid (few points per axis, wide range) to find the promising region cheaply.
Stage 2: a FINE grid zoomed into the neighbourhood of the coarse best. The
fine grid was narrowed to a clean 3 x 3 = 9 (half-width == step, so it lands
exactly on best-step / best / best+step) -- down from the original
7 x 7 = 49, to cut total runtime roughly in half (30 coarse + 9 fine = 39,
vs 30 + 49 = 79). The coarse ranges are illustrative, not validated.

Outputs: `HAT_groin_sweep_results.csv` (every combination's RMSE, written
incrementally -- safe to inspect mid-run), `HAT_groin_sweep_heatmap.png`
(fine grid, M x fraction, coloured by RMSE), and the single best combination
printed. Two profile figures follow: the best-fit profile against observed
(the direct visual answer to "how close does the best combo actually get",
not just its aggregate RMSE), and the top-N profiles overlaid (whether there
is one clear winner or a family of comparably good fits, i.e. M and fraction
trading off against each other, which the heatmap's single best-marker
cannot show). `load_observed_target()` here is a reference copy for printing
only -- each worker loads its own copy, so this script never needs the real
CASCADE inputs itself.

From the script's original header: "this script has not been run end-to-end
against a real CASCADE install (not available in the environment that built
it) -- the resume/subprocess/incremental-save logic was tested with a mocked
worker, but treat the first real run as a shakedown." Run it in this folder,
beside `HAT_groin_hindcast_1967_2017.py` and `HAT_groin_sweep_single_combo.py`.

### HAT_groin_sweep_single_combo.py

Runs ONE (M, fraction) combination of the sweep in its own fresh Python
process, then exits. Launched as a subprocess by the sweep (its old header
called it `HAT_groin_sensitivity_sweep_v2.py`; the file is
`HAT_groin_sensitivity_sweep.py`), NOT run directly in a loop -- the point
is that each combination gets a clean process, so accumulated state
(Cascade/Barrier3D objects, joblib worker pools, whatever else builds up over
30 sequential simulations in one process) can never carry over from one
combination to the next. That is the fix for the 0xC0000005 (Windows access
violation) crash that showed up partway through an earlier in-process sweep
-- process isolation guarantees a full OS-level cleanup between every run,
regardless of the underlying cause.

It prints exactly one line to stdout on success, in an easily parsed format,
`RESULT_RMSE=<value>`. Everything else printed is normal `run_one()` logging,
which the orchestrator ignores. A non-zero exit code (including an access
violation) is treated as "this combination failed" -- the sweep moves on.

The profiles go to the same tree the sweep READS from. This said
`scripts/groin/` -- the pre-reorg location -- so every profile landed in a
directory nothing looked at, and the best-fit / top-N figures silently
plotted whatever stale arrays were left in the real one. The path is split
across two lines, which is why the earlier bulk path fix missed it.

### HAT_solve_edge_be.py

Solves the D2 / D12 edge background-erosion correction for the 1967 rig.

**Why this exists.** `HAT_groin_hindcast_1967_2017.py` carries
`EDGE_BE_RATES_GIS = {2: 5.0, 12: 10.0}   # "(solve for this)"` --
placeholders that were never solved. They are not a minor detail here. The
reduced reach is 11 real domains inside 41, so D2 sits THREE domains from the
groin pair at D5/D6: an imposed rate at the edge diffuses into the pair
within a few years and lands directly on the signal the sweep is fitting. The
observed edge rates over 1967-2023 are +1.18 and +1.35 m/yr, so the
placeholders are 4-7x too strong.

**What is solved, and against what.** The same method the production hindcast
uses for GIS 1 / GIS 90 -- the site config records those as "solved on the
edgeBE road_bdm base run". A trial edge rate is imposed, the base run is
stepped, and the modelled change at that domain is compared against the
surveyed change from the fixed 1967 datum. The rate is then adjusted and the
run repeated.

GROIN OFF, ALWAYS. If the edge were solved with the groin attached, the
correction could absorb groin signal and the sweep would afterwards be
fitting M against a background that had already eaten part of its effect.
The whole point of an edge correction is that it is STRUCTURAL -- a fix for
the buffer's artificial orientation -- so it must be solved with the
structure of interest switched off.

Both edges are solved together rather than one at a time. They are nearly
independent (opposite ends of the reach) but not exactly, since each one's
signal diffuses inward, so a joint secant step converges without pretending
the coupling is zero.

**Method.** Secant iteration on each edge. The shoreline response to an
imposed background rate is close to linear over this range, so two
evaluations bracket the root and each further step refines it. Converges in
2-3 iterations; the cap exists so a non-convergent case stops rather than
spinning. Each trial runs 1967-2024 on the 1967-2024 storm file
(`END_YEAR = 2025`, exclusive; 58 states) so the 2023 survey is inside the
run. The 1971/73 fills are history, applied to every run in the matrix -- the
runner requires them, and omitting them would push their signal into the
solved edge rate. `run_one` returns a run NAME and writes the matrix to disk;
it does not hand back a Cascade, so the solve reads the matrix the same way
the sweep worker does.

In the wet/dry column names the SECOND year is the survey year; the first is
the 1967 datum. Matching the first one silently returns the datum for every
column.

Writes the solved rates, the observed and modelled change, the fit year and
"groin: off" to `edge_be_solved.json` for the sweep to read. Nothing else may
run while this does -- every CASCADE construction writes the shared
`Hatteras-CASCADE-parameters.yaml`, and concurrent writers corrupt it.

The header's first line is kept on the docstring's opening line on purpose:
`--help` prints `__doc__.split("\n")[0]` as its description.

### HAT_groin_trajectory_target.py

The observed groin-fillet TRAJECTORY, 1967-2023, and how to score against it.

**Why a trajectory and not an end state.** M and f are not separable from a
single fillet measurement. A groin with a high trapping rate and a low
deterioration floor, and one with a low rate and a high floor, reach a
similar fillet by 2024 along completely different paths -- the first builds
fast and then relaxes, the second builds slowly and holds. Every scalar
target tried on the 1984-2004 and 2004-2024 windows produced a RIDGE of
equally good (M, f) pairs for exactly this reason. The path is what separates
them, and the path is measured: 24 dated wet/dry surveys between 1967 and
2023.

**Why 1967 and not 1984.** The fillet's whole life, measured from
`Change_from_wetdry_1967_D2_D12.csv`:

    build    1967-1978     0 -> 117 m
    plateau  1978-2004   117 -> 150 m
    decline  2004-2023   150 ->  74 m

Both hindcast windows begin AFTER the build. A run starting in 1984 inherits
the fillet in its initial condition rather than predicting it, so its fillet
CHANGE carries almost no information about M. The 1967 window is the only
one containing the growth that M controls, while the post-2004 decline is
what constrains f. One run, both knobs, separated by different stretches of
the same curve.

**No paired baseline is needed.** Elsewhere the modelled fillet is
differenced against an M = 0 run, because a run starting in 1984 begins with
the real fillet already in its initial shoreline and the baseline is what
removes it. Here the observations are themselves changes from a fixed 1967
datum, and the model starts before the groin existed, so both sides are
already referenced to the same zero. The fillet is simply the change in
(D5 - D6) since year 0.

**Sign.** The wet/dry table is landward-positive (+ = erosion), matching
Barrier3D's x_s. The fillet is D5 minus D6: positive when the downdrift
domain has retreated further than the updrift one, which is what a groin
builds.

The functions, from their original docstrings:

- `observed_trajectory(table)`: observed fillet against the fixed 1967 datum,
  per dated survey. Returns a DataFrame indexed by year with a `fillet_m`
  column, sorted by year. Raises `FileNotFoundError` if the table is absent,
  `ValueError` if neither groin domain has any dated column. The datum year
  is 0 by construction; it is carried explicitly so a model trajectory that
  starts at 0 is compared against a target that does too.
- `model_trajectory(shoreline_m, start_year, updrift_pad, downdrift_pad)`:
  modelled fillet per simulated year, referenced to year 0. `shoreline_m` is
  a [state x padded domain] array, metres, landward-positive; `start_year` is
  the calendar year of state 0. Returns a DataFrame indexed by year with a
  `fillet_m` column.
- `score(model, observed)`: RMSE between a modelled and observed trajectory,
  on shared years. Only years the survey actually sampled are compared -- the
  model has a value every year, the observations do not, and interpolating
  the survey onto the model's grid would invent data and weight the
  well-sampled 2010s the same as the sparse 1970s. Returns a dict with
  `rmse_m`, `n_years`, `bias_m` and the paired frame.

### HAT_groin_threeway_hindcast_1967_2017.py

A sibling of the rig runner, from the same 1967-1997 original and with the
same END_YEAR convention, deterioration ramp, edge correction, nourishment
sources and units as `HAT_groin_hindcast_1967_2017.py` above. Those passages
are written once, in that section, and apply here word for word: the
sandbox Cascade and the hook, END_YEAR, the deterioration ramp and floor
history, the edge correction (Section 3b), the nourishment sources, volume
choice and domain range (Section 3c), the units reminder, where the fill
volume goes, and the run figures and GIF.

`RUN_MODE = "three_way"` (the default) runs THREE hindcasts in one execution,
driven by `RUN_RECIPES`: baseline (no groin, no nourishment),
nourishment-only (no groin), and the full model (nourishment + groin) -- all
other parameters (storms, edge BE correction, groin M and fraction) held
identical, only nourishment on/off and groin on/off varying. This is what
feeds `HAT-buxton-hindcast-groin-test/comparison/HAT_three_run_comparison.py`
directly. Each recipe is (label, run_key, enable_nourishment,
run_name_suffix). The baseline gets a DISTINCT suffix (`_no_BN`) from the
other two -- it shares run_key `"no_groin"` with the nourishment-only run,
and without a distinct suffix the two would produce the identical run name
and overwrite each other on disk.

`RUN_MODE = "single"` falls back to the old behaviour: run whatever is in
`RUN_MATRIX`, nourishment controlled by the single module-level
`ENABLE_HISTORICAL_NOURISHMENT`. `RUN_MODE` only affects what `main()` does.

`GROIN_TRAPPING_RATE_M_YR` (50) and `GROIN_DETERIORATION_FRACTION` (0.90) are
marked placeholders -- fill in the sweep's winning combination before running
the full-model recipe for real; the values were carried over from before the
sweep existed. The floor-fraction history is the rig runner's (M/3, then 0.5).

**Nourishment is togglable per run** through
`build_nourishment_arrays_from_manual_inputs(enable=...)`. Unlike earlier
versions, `ENABLE_HISTORICAL_NOURISHMENT` no longer empties `HAT_BN_YEARS` /
`HAT_BN_VOLUME_BY_DOMAIN` at import time -- the full schedule stays intact,
and the function takes an explicit `enable=` override so one script run can
produce nourishment-on and nourishment-off runs. `enable=None` (the default)
falls back to `ENABLE_HISTORICAL_NOURISHMENT`, which is what every existing
caller that passes no argument gets, unchanged. Pass True/False to override
per call, as `main()`'s `RUN_RECIPES` loop does.

The module-level `NOURISHMENT_MANAGEMENT_ON` is superseded per run inside
`run_one()`, which derives it fresh from whatever schedule was passed to that
call; the module-level one is kept only for backward compatibility with
anything that might reference it. Deriving it per call means a run built with
an empty (nourishment-disabled) schedule gets `beach_nourishment_module=False`
everywhere, not just zero volume. This matters: even with `nourish_now` never
firing, `beach_nourishment_module=True` still runs the whole
`NourishmentManager.update()` every year (narrow_break checks, beach-width
tracking, etc.) -- a clean "no human intervention machinery at all" baseline
needs this False, not just quietly zeroed.

Differences from the rig runner that matter: it writes to `output/raw_runs/`,
not `calibration/groin_rig/`; it still builds file lists from the flat
`data/hatteras_init/dunes/2009/` and `topography/2009/` layout, which no
longer exists (see the rig runner's `build_file_lists()` history), so it
stops at its missing-file check as written; and its groin import still names
`scripts.groin_module.hindcast_groin_test.version_control.HAT_groin_module`,
which no longer exists (the rig runner imports `cascade.groin`).

Its old header also said `HAT_groin_sweep_single_combo.py` relied on this
file's `run_one()` / `build_nourishment_arrays_from_manual_inputs()`; the
worker imports `HAT_groin_hindcast_1967_2017.py`, not this one.

Plot with `HAT_three_run_comparison.py` (RUN_MODE "three_way") or
`comparison/HAT_groin_effect_comparison.py` (RUN_MODE "single", comparing
exactly two named runs).

### HAT_plot_groin_runs.py

Plots for the 1967-2017 groin-deterioration-test runs. Loads saved shoreline
matrices (the `*_shoreline_matrix.npy` files the runners write) and produces:

    FIG 1  Shoreline CHANGE (m, end - start) per domain          [position change]
    FIG 2  Shoreline change RATE (m/yr) per domain               [LRR-style]
    FIG 3  Shoreline POSITION over time, selected domains        [trajectories]
    FIG 4  model vs observed: positions above, change since 1967 below
    FIG 5  planform against the 1967 alongshore mean
    FIG 6  (multi-run only) run-vs-run DIFFERENCE, e.g.
           groin_be minus no_groin -> the isolated groin signal

It matches the main hindcast's conventions: model orange (#FF8C00), groin red
(#B71C1C) at the D5/D6 boundary, GIS domain x-axis, real-domain focus. No
CoastSat dependency -- it is for inspecting the model runs themselves. Edit
`RUNS` to list the run folder name(s) under `output/raw_runs/`; one run gives
the single-run figures, two or more also give the difference figure against
the first as baseline. `RUNS` was set on 2026-08-30 to the run the sweep's
best cell corresponds to, produced by `HAT_groin_hindcast_1967_2017.py` on
1984-start/v1 at M = 60, f = 0.6. (The rig runner now writes to
`output/calibration/groin_rig/`, so a rig run is not where this script looks
unless it is moved.)

**END_YEAR is EXCLUSIVE** (`RUN_YEARS = END_YEAR - START_YEAR`), so
`END_YEAR = 2018` means the run's actual final modelled year is 2017, not
2018. Every figure label derives the true final year from the loaded data's
own shape (`_final_model_year()`: START_YEAR + nt - 1) rather than trusting
END_YEAR -- END_YEAR is bookkeeping only, never a figure title. If runs
disagree in length it warns, and uses the first run's.

**Sign and rate conventions** (match the main hindcast exactly). Main script:
`change_rate = (x_s[-1] - x_s[0]) / (nt-1)`; then `*= -1` if FLIP. Raw x_s
increases LANDWARD (retreat). With `FLIP_SIGN_MODEL = True`, the plotted
rate is + = seaward/accretion, - = landward/erosion.

**Observed comparison.** `fig_model_vs_observed()` reads the validated wet/dry
change table (`Change_from_wetdry_1967_D2_D12.csv`, built by
`HAT_geometric_distance_sanity_check.py`) instead of the raw dune-line CSVs --
see `comparison/README.md` for the reasoning (x_s is conceptually closer to a
water-line proxy than to the dune line, and the wet/dry line was confirmed
seaward of the dune line at every domain). `OBSERVED_YEARS` are the checkpoint
years overlaid, roughly every ~10 years, chosen from years with full (or
near-full) D2-D12 coverage -- 2017 itself only has 6/11 domains (D8-D12
missing), so 2018 (11/11, one year later than the model's actual endpoint) is
used instead. The two-panel figure: TOP the absolute shoreline position,
start (1967) and modelled end; BOTTOM the change since 1967, modelled vs
observed at every year in `OBSERVED_YEARS` (the validity check), coloured
chronologically (blue = older, red = newer) so several decades overlay
readably. The change panel is the apples-to-apples overlay: both are metres
of change from 1967, + = seaward/accretion (flipped from raw).

The table path was REPOINTED on 2026-08-30. It named
`HAT-buxton-hindcast-groin-test/input_prep/shoreline_position/output/`, which
does not exist -- the table lives under `HAT-groin-buxton-output`. The script
did not error: it fell through to "observed wet/dry change table not found --
model only" and drew the model against nothing, which looks like a finished
figure. The sweep reads the same table from its own correct path
(`HAT_groin_sweep_config.py:601`, as cited then), so the FIT was never
affected; only these figures were. The script now refuses to start
without the table rather than draw model-only figures that look like
comparisons.

The planform figure is the position-mode view (main-hindcast convention): the
REAL island planform, with the year-0 (1967) shoreline as the 0 reference, and
how it changes -- cross-shore position relative to the year-0 alongshore mean,
ocean at bottom. "Real island orientation as position 0."

`PROJECT_BASE_DIR` is derived, not hardcoded, for the same reason as in the
rig runner: it was `r"/"`, which silently resolved every data path to the
filesystem root.

### HAT_check_nourishment_applied.py

Checks whether historical beach nourishment (1971, 1973) actually reached x_s
inside CASCADE, or was silently skipped by NourishmentManager's narrow_break
early-return (`beach_dune_manager.py`, around line 650) -- which exits BEFORE
the nourish block, regardless of `nourish_now` / volume having been set
correctly by the runner.

Console prints from the hindcast runner only confirm what the RUNNER requested
(set before `cascade.update()` is called that year) -- they cannot see whether
Cascade's internal `update()` actually executed the nourish branch. This
script checks the one signal that can: `nourishment_volume_TS`, which is only
written on the exact line inside the nourish branch, so a zero there when a
nonzero volume was requested means it was blocked. The expected requested
volumes (m^3/m) are from the actual console log -- used only to print a
side-by-side comparison, not to change the check itself.

Point `RUN_DIR` at the saved run folder (the one containing the
`{run_name}.npz` file), then run. The home-directory path it used to name was
anchored on the repo root on 2026-09-14 (ORGANIZATION.md rule 5). The run it
names, `HAT_1967_2018_edge_calibrated_no_groin`, is looked for under
`output/raw_runs/`; rig runs are now written to `output/calibration/groin_rig/`.
