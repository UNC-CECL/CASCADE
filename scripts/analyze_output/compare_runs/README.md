# compare_runs - comparisons of finished runs, grouped by question

Nothing here runs the model. Each subfolder answers one question and keeps its
scripts together because they import one another; each has a README with the
reasoning behind its scripts.

```
compare_runs.py              general tool: overlay up to four runs against the
                             smoothed CoastSat LRR. RUNS_TO_COMPARE ships
                             EMPTY; fill it first -> comparisons/<COMPARISON_NAME>/
hindcast_vs_observed/        the hindcast against the shoreline and dune line
    rate_windows.py              every window          -> comparisons/model_vs_observed/
    target_comparison.py         which target to grade -> comparisons/target_comparison/
    smoothing_scale.py           does the LOWESS width matter
    smoothed_lowess7_with_cascade.py   model over the smoothed sheets
matrix_vs_observed/          every run-matrix scenario against the observations
    matrix_vs_observed.py        per-run figures + scores -> comparisons/matrix_vs_observed/
    matrix_management_ladder.py  management added rung by rung -> raw_runs/matrix/figures/
adoption_2026-09-28/         the matrix before vs after the 09-28 adoption (one-off)
    adoption_before_after.py     score, then report
    adoption_before_after_figures.py
offset_source/               dune-line vs shoreline island offset
    offset_source_comparison.py  -> comparisons/offset_source/
```

Where to start: `hindcast_vs_observed/` is the live model-vs-observation work;
`matrix_vs_observed/` is the scenario-by-scenario view; the other two are
studies that answered one question.

## The script in detail

Moved here from `scripts/analyze_output/README.md` on 2026-09-30.

## compare_runs.py

A general run-vs-run and run-vs-CoastSat tool: loads up to four finished
CASCADE runs from their saved shoreline change rate CSVs, overlays them on one
figure, and adds the LOWESS-smoothed CoastSat LRR for reference.
**`RUNS_TO_COMPARE` ships empty**: fill it (and `COMPARISON_NAME`) before
running.

Each run must already have been run by the hindcast runner, which writes
`<run_dir>/tables/shoreline_change_rate.csv` (before the 2026-09-10 layout
change, `<run_dir>/{run_name}_shoreline_change_rate.csv`; `run_layout`
resolves either). Columns: `gis_domain | change_rate_m_yr | lrr_m_yr | lrr_r2`.
The run folder is resolved from (run_name, period, preset, arm) by
`cascade_pipeline.run_registry.find_run_dir`, never joined by hand; it raises
naming the arms a run *is* under, which reads as "it is over there" rather than
"it never existed". The domain ids are read from the file's own `gis_domain`
column: this used to slice padded indices out of a 120-row CSV and pair them
positionally with a hand-built 1-90, so a run written in the other alongshore
order would have shifted every rate against its label with nothing raised.

**Which rate column is read is a method choice, not a spelling.** The file
carries `lrr_m_yr` (an OLS slope through the annual shoreline positions) and
`change_rate_m_yr` (the endpoint rate). This script reads `lrr_m_yr` because
the CoastSat target it is plotted against is an OLS rate too; putting an
endpoint rate up against an OLS one would move the difference between two
estimators into the residual panel, where it reads as model error.

**`RUNS_TO_COMPARE` fields:**

| field | meaning |
|---|---|
| `run_name` | folder name and filename prefix for this run; must match the name the run was made under |
| `period` | the run's period directory, e.g. `"1996_2010"` |
| `preset` | source/sink preset the run was made under (`calibBE`, `edgeBE`, `zeroBE`), the directory below the period |
| `arm` | optional forcing arm. Omit for the calibration arm. Naming no arm never silently searches: a run forced off the calibration wave climate must be asked for by name |
| `label` | legend label |
| `start_year` | which CoastSat period is drawn solid for this run, and which panel it appears in in the two-period figure |
| `sort_key` | optional; orders the run along the colour gradient (lower = lighter), typically the swept parameter. Without it, list order is used |
| `color` | optional explicit colour that overrides the gradient for this run only |
| `run_dir` | optional escape hatch: full path to a run outside the raw_runs tree (another drive, a collaborator's folder, an archive). Overrides the resolver for this run only. Prefer period/preset: a hand-typed path is how a figure ends up drawing a run other than the one it names |

Colours are sampled from `RUN_COLORMAP` (`YlOrRd`) at each run's rank along
`sort_key`, within `RUN_COLORMAP_RANGE` (0.35-0.95, avoiding near-white and
near-black), so reordering or adding runs never means re-picking hex codes.
`plasma` or `inferno` keep more contrast for more than five runs.

**The four runs the list used to hold are gone.** `HAT_1984_2004_SStest`,
`HAT_1984_2004_FinalSS`, `HAT_2004_2024_SStest` and `HAT_2004_2024_FinalSS`
were addressed by absolute paths into `output/raw_runs/source&sink_tests/`
and the flat `raw_runs/<name>` level; neither exists any more, and none of the
four names appears in `run_index.csv` or in `superseded_20260828/` (now
`output/archive/2026-08-28_full-tree/`). Their four figures sat at
`output/comparisons/source_sink_zones/` (2026-06-19) until 2026-09-17, when
they were deleted as unregenerable; nothing cited them (checked 2026-09-02 and
2026-09-17).

**CoastSat:** each `COASTSAT_DATASETS` entry points to a
`transect_lrr_full.csv` (one row per ~50 m transect). LOWESS is applied at
transect resolution, then averaged to domains. The period whose start matches
most runs is drawn solid (active) and the other faded (reference);
`ACTIVE_PERIOD_START` overrides that. Domains 1-10
(`LOWESS_SKIP_SOUTHERN_DOMAINS`) show raw transect scatter instead of LOWESS,
because Oregon Inlet boundary effects dominate there and smoothing hides the
sharp gradient; this matches the hindcast script. `LOWESS_WINDOW_DOMAINS`
lists one or two widths (7 since 2026-09-28; [7, 10] before), each with its
own style.

**Colour palette:** CoastSat is the cool blue family, light to dark = less to
more processed (transect scatter and domain means very light blue #9ECAE1,
7-domain LOWESS #6BAED6, 10-domain LOWESS #08519C). CASCADE runs are the warm
orange-red gradient. Blue against orange-red is the most colour-blind-safe
pairing (deuteranopia and protanopia); within each family lighter to darker
encodes less to more smoothing or a lower to higher parameter; the period is
carried by line style, not colour.

**Output:** `output/comparisons/<COMPARISON_NAME>/`:
`<name>_diagnostic.png` (quick multi-run check), `<name>_annotated.png` (the
publication figure with geographic annotations), the two-period figure
(1984-start runs left, 2004-start right, one y axis and one legend; skipped
for a period with no runs), and `<name>_residuals.png` (each run minus the
active CoastSat LOWESS). Per-run fit statistics are shown in the legend.

**History of fixes recorded in the code:**

- The 2004 CoastSat dataset pointed to `2004_2024_specific_dates`, a folder
  that did not exist, so the load failed silently and every 2004 run was
  compared against the 1984 CoastSat data in every figure (residuals reaching
  +9 m/yr instead of ±1-2, no transect scatter in the 2004 panel). The script
  now stops if a run's CoastSat period failed to load.
- The pier label heights were 85 and 70 (meant as percentages) in
  axes-fraction coordinates, which placed labels 85x the axes height above the
  plot; the tight bounding box grew to 300+ inches on save and crashed with an
  out-of-memory error. They are 0.85 and 0.70 now.
- Legends and captions placed below the canvas (negative y) made
  `bbox_inches="tight"` expand to include them, which on the wide two-period
  figure produced a ~490-inch canvas and another out-of-memory crash. Every
  figure now reserves its margins first and places legend and caption in
  figure-fraction coordinates on the canvas, saving with the figure's own
  bounds. The two-period figure's bottom margin is 0.20 (xlabel ~0.04, legend
  ~0.10, caption ~0.02); it was 0.30 before the legend was deduplicated.
  `loc="best"` used to put the diagnostic legend on top of the data.
- The two-period legend is deduplicated by visual style: both panels draw
  their CoastSat in the same colours, so one entry per style is enough and the
  period is dropped from the label (each panel's title names its period).
  Deduplicating by label text used to list every CoastSat entry twice.

**Dead code** (left in place; it runs but nothing uses it): `TRANSECT_DATASETS`
is not referenced anywhere (`COASTSAT_DATASETS` already carries the transect
data), `_gis_to_pad` is never called, and `sys`/`Path` are imported twice.

<details><summary>Function notes (the original docstrings)</summary>

**`assign_run_colors()`**

```
Auto-generate a light->dark color gradient for a list of run configs.

Runs are ranked by `sort_key` if every run provides one; otherwise by
their position in the list (so simply listing runs low-to-high parameter
value works without setting sort_key explicitly). Colors are sampled
evenly across RUN_COLORMAP_RANGE of RUN_COLORMAP, lightest first.

Any run with an explicit, non-None `color` field keeps that color
untouched and is excluded from the rank-based assignment - useful for
pinning one run (e.g. a baseline) to a fixed color while the rest of the
sweep auto-generates.

Parameters
----------
runs : list of dict - entries from RUNS_TO_COMPARE (or equivalent)

Returns
-------
list of dict - same entries, each with a resolved 'color' key (hex string)
```

**`load_run_rates()`**

```
Load the shoreline change rate CSV produced by HAT_hindcast_1984_2024_old version.py.

Parameters
----------
run_name : str       — used to resolve the default folder AND, inside it,
                        the rate CSV, whose name run_layout knows in both
                        the current and the pre-2026-09-10 layout. Which
                        is used is unchanged regardless of how the folder
                        was resolved.
period   : str, optional — "1984_2004" or "2004_2024".
preset   : str, optional — source/sink preset the run was made under.
arm      : str, optional — forcing arm; None means the calibration arm.
run_dir  : str, optional — ESCAPE HATCH. Full path to the folder holding
                        that CSV, for a run outside the raw_runs tree.
                        If given it OVERRIDES the resolver entirely.
                        Otherwise period and preset are both required and
                        find_run_dir locates the run, raising with what
                        IS on disk when it is absent.

Returns
-------
gis_ids   : int array — GIS domain IDs, READ FROM the CSV's gis_domain
                        column rather than assumed. Normally 1-90; a
                        short array means the run wrote fewer rows, which
                        is warned about rather than padded over.
rates_myr : float array, same shape — the run's LRR in m/yr, from
                        RUN_RATE_COL, aligned to gis_ids by construction.
run_dir   : str                                   — full path to run folder (resolved)
```

**`load_transect_data()`**

```
Load individual transect LRR values from transect_lrr_full.csv and derive
along-coast distance by spreading each domain's transects evenly across its
500 m band (mirrors 6-scr-smooth/lowess_method_comparison.py: load_transect_csv).

Returns
-------
domain_ids    : int array   — CASCADE domain ID for each transect
lrr_values    : float array — LRR (m/yr) for each transect
along_coast_m : float array — cumulative along-coast distance (m)
All three arrays share the same length (one entry per transect).
Returns (None, None, None) on load failure.
```

**`lowess_smooth_transect_to_domains()`**

```
Apply LOWESS at transect resolution using physical along-coast distance (m) as x,
then aggregate smoothed values to CASCADE domain resolution by averaging within
each domain.  Mirrors smooth_transect_df() + aggregate_to_domains() from
6-scr-smooth/lowess_method_comparison.py.

window_domains is converted to km (× DOMAIN_SPACING_M) so the physical window
is consistent regardless of transect density.

Returns
-------
gis_x    : int array   — domain IDs that have at least one transect
smoothed : float array — domain-averaged smoothed LRR (m/yr), same length as gis_x
frac     : float       — LOWESS frac used (for logging)
```

**`splice_lowess_with_raw_south()`**

```
Trim a LOWESS curve so it starts north of the southernmost `skip_n`
domains, leaving that southern zone to show raw transect scatter only.

Ported from HAT_hindcast_1984_2024_old version.py's function of the same name for
visual consistency between the two scripts - domains 1-skip_n are
boundary-affected (Oregon Inlet dynamics) and LOWESS smoothing there can
obscure the sharp gradient rather than clarify it, so the LOWESS line is
simply not drawn there; the raw scatter (already restricted to this same
zone via RAW_LRR_SOUTHERN_ONLY) carries the signal instead.

Parameters
----------
win_gis_x    : int array   - GIS domain IDs from the LOWESS result
win_smoothed : float array - LOWESS-smoothed LRR (m/yr)
skip_n       : int         - domains 1..skip_n are excluded from the
                              returned line. Defaults to
                              LOWESS_SKIP_SOUTHERN_DOMAINS.

Returns
-------
plot_x, plot_y : arrays - the LOWESS curve restricted to domains > skip_n
```

**`add_geographic_annotations()`**

```
Draw the standard Hatteras geographic annotation layer onto ax.
X-axis must be in GIS domain IDs (1–90).
```

**`plot_diagnostic()`**

```
Quick diagnostic plot — all runs + CoastSat + geographic annotations on
one panel. Uses the same drawing helper as plot_annotated/plot_two_period
so all comparison figures stay visually consistent.

Previously this plot had NO geographic annotations and used
loc="best" for the legend, which with 4+ runs landed on top of the
data (see uploaded screenshot) - both fixed below.
```

**`_draw_comparison_panel()`**

```
Draw geographic annotations + CoastSat curves + model run lines onto a
single axis. Shared by plot_annotated (one panel) and plot_two_period
(two panels, called once per period) so both stay visually identical.

Parameters
----------
show_reference_period : bool
    True  -> the inactive CoastSat period is still drawn, faded, with a
             "(ref)" label (useful in plot_annotated's single-panel view,
             where showing the other period for context is informative).
    False -> only the active period (cs["period_start"] == active_period)
             is drawn at all; the inactive period is skipped entirely.
             This is what plot_two_period uses, since each panel is
             already dedicated to one period - showing the other period
             there just duplicates what the OTHER panel is for.

Returns
-------
model_handles, cs_handles : lists of Line2D proxies for the legend
```

**`plot_two_period()`**

```
Two-panel figure: left = all runs with start_year=1984 (vs. 1984-2004
CoastSat), right = all runs with start_year=2004 (vs. 2004-2024 CoastSat).
Shares one y-axis range and one combined legend below both panels.

If every run shares the same start_year, the empty panel is skipped and
a single-panel figure is produced instead (so this is always safe to call
regardless of what's in RUNS_TO_COMPARE).
```

**`plot_residuals()`**

```
Residual panel — each model run minus the active CoastSat LOWESS curve.
Helps identify where each run over- or under-predicts relative to observations.
Only produced if PLOT_RESIDUALS = True.
```

</details>
