# analyze_output - questions that span more than one run

Nothing here runs the model. Each script reads finished runs and writes to
`output/comparisons/`, resolved as `hat_figure_style.COMPARISONS_ROOT`
(2026-09-18); a figure finished for the manuscript also goes to
`output/figures/` (numbered layout, map in its README.md), from `scripts/figure_making/`.

```
compare_runs/            comparisons of finished runs, one subfolder per
                         question (map + details in compare_runs/README.md)
    compare_runs.py          general tool: up to four runs vs CoastSat.
                             RUNS_TO_COMPARE is EMPTY; fill it first
    hindcast_vs_observed/    LIVE. The hindcast vs shoreline and dune line
                             (rate_windows.py, target_comparison.py, ...)
                             -> comparisons/model_vs_observed/, target_comparison/
    matrix_vs_observed/      every run-matrix scenario vs the observations
    adoption_2026-09-28/     the matrix before/after the 09-28 adoption
    offset_source/           dune-line vs shoreline island offset
overwash/
    compare_overwash_figures.py   a run's overwash (Qow) against the observed
    compare_overwash_observed.py  record; both resolve their run through
                                  run_registry (HAT_1984_2004_calibBE_road_bdm_groin,
                                  read from raw_runs/archive/2026-09-24-pre-metres/
                                  since 2026-09-30) -> comparisons/overwash/
    superseded_20260918/          plot_overwash.py (dead absolute path) and an
                                  early copy of Roya's Pea Island script
smoothing_vs_cascade/
    smoothing_vs_cascade.py   what the LOWESS smoothing does to the
                             comparison. Names HAT_1984_2004_SQ_BE_Hs2p0, which
                             no longer exists anywhere under raw_runs/: point it
                             at a current run before running
```

The observed overwash record and its own figures are
`scripts/input_prep/8-overwash-analysis/`; this folder only compares a RUN
against it.

The smoothing figures and the source/sink-zone figures were deleted from
`output/comparisons/` on 2026-09-17 (unregenerable, uncited); running either
script recreates its folder, so give it a current run first.

The run tree is addressed through `cascade_pipeline.run_registry`, never by
building a path: a run's NAME describes its scenario and its PATH describes its
forcing, and the two are easy to get wrong by hand.

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning, the
choices behind it and its history are here, one section per script. Moved out
of the scripts on 2026-09-30, when they were brought in line with
`scripts/STYLE.md`.

Shared history: every script here drew in matplotlib's defaults until
2026-09-17, when each started calling `apply_style()` from
`site_layer/hat_figure_style.py` at import. Their paths used to be absolute
literals typed on one machine; they are all found by searching upward for the
repo root now (ORGANIZATION.md rule 5), so they follow the checkout and survive
a file changing depth.

### compare_runs/

The compare_runs scripts' details moved on 2026-09-30 into the README of the
subfolder each now lives in: `compare_runs/README.md` is the map.

### overwash/compare_overwash_observed.py and overwash/compare_overwash_figures.py

A run's overwash (Qow) against the observed record from the imagery. Both
read one run's `.npz` and the observation workbook (`Overwash_Matrix` sheet,
resolved by `site_layer/hat_overwash.py`) and write to
`output/comparisons/overwash/`. `compare_overwash_figures.py` is the fuller
of the two: the same two figures in a revised palette, plus the spatial
figures. The two are largely duplicates, and a future cleanup could keep one.

| figure | switch | shows |
|---|---|---|
| 1, stacked | `PLOT_STACKED` | model Qow heatmap above (every model year, imagery years highlighted), observed overwash below (binary, imagery years; hatched where there is no image); shared domain axis and island-section bar |
| 2, contingency | `PLOT_CONTINGENCY` | each image x domain classified hit / miss / false alarm / correct rejection; summary statistics printed |
| 3 and 4, spatial | `PLOT_SPATIAL` (figures script only) | overwash frequency alongshore, model (mean with IQR shading) against observed (bars): normalised on one panel, and as two panels |

**The run:** `HAT_1984_2004_calibBE_road_bdm_groin`, 1984-2004, calibBE,
resolved through `cascade_pipeline.run_registry`. *Fixed 2026-09-30:* the run
moved to `raw_runs/archive/2026-09-24-pre-metres/` when the matrix was
archived on 2026-09-24, and both scripts crashed on import until they were
pointed at it there (`RUN_KIND`, `RUN_TAG`). Whether to compare a current
matrix run instead is an open choice.

**Threshold:** a model cell counts as overwash if Qow > `QOW_THRESHOLD`
(dam³/yr). Start at 0 (any non-zero flux) and raise it if there are too many
false alarms.

**Imagery-year remap** (figures script): some imagery post-dates the storm it
captures, so `YEAR_REMAP` maps an image year to the model year that best
represents the event visible in it. The May 2004 imagery shows Hurricane
Isabel's (September 2003) deposits, so 2004 is remapped to 2003. 1996 imagery
is flagged as poor quality (`POOR_QUALITY_YEARS`).

**Palette** (figures script): muted coastal earth and ocean tones, with
consistent warmth across all elements; the section bar uses desaturated
versions of the warm and cool families in the contingency legend. Section
bar: warm linen #CABB9E for villages, dusty maritime blue #8DAFC2 between
them. Contingency: hit deep maritime green #2E7D5A, miss deep maritime blue
#2A5F8F, false alarm warm brick #A85840, correct rejection warm cream #F2EFE9,
no imagery warm greige #BEB9B4. Observed overwash uses the same brick as a
false alarm. The look was modelled on AGU / Nature Geoscience figures.

**To use on another run:** set `RUN_NAME`, `RUN_PERIOD`, `RUN_PRESET` (and
`RUN_KIND`/`RUN_TAG` if it is not in the matrix), then adjust `START_YEAR`,
`END_YEAR`, the domain constants and `SECTIONS` to match it, and choose
`QOW_THRESHOLD`. If the run is not where you say, `find_run_dir` raises naming
where it *is* on disk.

**History:** the paths used to be absolute literals under a folder spelling
that no longer exists (`input_preperation`), so the scripts could not read
their observations and crashed on the first save; they are anchored on the
repo root now. The workbook left `scripts/input_prep/8-overwash-analysis/`
on 2026-09-10 and is resolved by `site_layer/hat_overwash.py` since
2026-09-18. The `.npz` loader substitutes stand-in classes for anything the
pickled model references but cannot import, and is duplicated from
`plot_overwash.py` so each script stands alone.

<details><summary>Function notes (the original docstrings)</summary>

**`plot_contingency()`**

```
For each imagery year × domain where obs data exists, classify as:
  Hit (H)             : Qow > T  AND  obs = 1   (green)
  Miss (M)            : Qow <= T AND  obs = 1   (blue)
  False Alarm (FA)    : Qow > T  AND  obs = 0   (orange)
  Correct Rejection(C): Qow <= T AND  obs = 0   (light grey)
  Not Assessed (NA)   : obs = NaN               (medium grey)
Rows: imagery years only. Columns: all 90 domains.
```

</details>

<details><summary>Function notes (the original docstrings)</summary>

**`compute_spatial_data()`**

```
Returns per-domain summary arrays used by both spatial profile figures.

qow_mean  : mean annual Qow per domain across all model years (dam³/yr)
qow_p25   : 25th percentile of annual Qow per domain
qow_p75   : 75th percentile of annual Qow per domain
obs_freq  : fraction of assessed imagery years with overwash per domain
            (NaN where no assessed imagery exists for that domain)
bar_colors: colour per domain based on section type (village vs inter-village)
```

</details>

### smoothing_vs_cascade/smoothing_vs_cascade.py

Applies LOWESS smoothing to the CoastSat shoreline change rates of both
periods (1984-2004, 2004-2024) and draws a CASCADE run against them. LOWESS is
applied to each period's series independently; it keeps the large-scale
alongshore pattern and removes per-domain noise.

| figure | shows |
|---|---|
| `overview_smoothed.png` | both periods, raw CoastSat (faded) with the LOWESS overlay |
| `smoothed_only_comparison.png` | both periods, LOWESS lines only (for presentations) |
| `combined_periods.png` | both periods' smoothed lines on one panel |
| `smoothing_sensitivity_<period>.png` | one panel per LOWESS window |
| `window_comparison.png` | raw CoastSat with every window in `COMPARE_WINDOWS_DOMAINS` overlaid, for choosing a window |
| `cascade_vs_lowess.png` | the CASCADE run against raw and LOWESS CoastSat, per period |
| `cascade_by_window.png` | the CASCADE run against one LOWESS window per panel |

Also `coastsat_smoothed_table.csv`. Output:
`output/comparisons/smoothing_vs_cascade/1984_2004/`.

**Bandwidth:** `LOWESS_FRAC` is the fraction of the data in each local fit:
0.10 is ~9 domains (more local), 0.167 ~15 domains (the default), 0.20 ~18
domains (smoother, loses finer structure). Window sizes in domains convert to
frac as n / 90, over the whole island regardless of how many domains a CSV
has data for: 5 domains ≈ 0.056 (2.5 km), 10 ≈ 0.111 (5 km), 15 ≈ 0.167
(7.5 km).

**`CASCADE_RUNS`:** one `dict(label, period, csv)` per run, where `period`
matches a CoastSat period label ("1984–2004") and `csv` is
`run_rate_csv("<run_name>")`, which finds the rate CSV in either run-folder
layout (`base=` for a run outside `output/raw_runs`). Add one colour to
`C_CASCADE` per run. An empty list skips the model figures.

**Cannot draw the model today** (open, needs your choice):

- `CASCADE_RUNS` names `HAT_1984_2004_SQ_BE_Hs2p0`, which no longer exists
  anywhere under `raw_runs/`, so the model is skipped and only the CoastSat
  figures are drawn.
- Even with a current run, `load_cascade_rate` expects columns
  `gis_domain_id` and `model_rate_m_per_yr`, but run tables have carried
  `gis_domain`, `lrr_m_yr` and `change_rate_m_yr` since the 2026-09-10 layout
  change. Which rate to read (OLS or endpoint) is a method choice, so this was
  not changed blind.

**Fixed 2026-09-30:** the script drew every figure and then crashed printing
its summary on a Windows console (arrows that cp1252 cannot encode); it now
switches its output to UTF-8 first.

**Styling note:** after `apply_style()` it sets its own rcParams (Arial,
dotted grid, no top/right spines) and its own period colours, so its figures
are not fully in the house style. The smoothing windows run cool to warm, teal
(5 domains) to amber (10) to crimson (15), so more smoothing reads as a warmer
colour, and `cascade_by_window.png` uses the same colour per window as
`window_comparison.png` so the figures can be read against each other. Legends
use `loc="best"` because the free corner differs between periods.

**CSV quirk:** some exports prepend the file name to the first column
("domain_lrr_summary.csvdomain_number"); the loader strips it.

**History:** every path used to be an absolute literal; the output one had
lost its drive and wrote figures to `C:\scripts\`, and the input ones still
spelled the folder `input_preperation`. They are anchored on the repo root
now, and the CoastSat CSVs are resolved through `hat_observed_rates.py`
(2026-09-18).

<details><summary>Function notes (the original docstrings)</summary>

**`run_rate_csv()`**

```
The shoreline change rate CSV inside one run folder.

RESOLVED, NOT JOINED. That file is tables/shoreline_change_rate.csv in the
new run layout and {run_name}_shoreline_change_rate.csv in the old one;
run_layout.resolve returns whichever is on disk, so a half-migrated tree
reads either way. `base` is the folder the run folder sits in -- pass it
for a run outside output/raw_runs.
```

**`load_cascade_rate()`**

```
Load a CASCADE shoreline change rate CSV produced by HAT_hindcast_1984_2024_old version.py.

Expects columns:
    gis_domain_id       – integer 1–90 for real domains, NaN for buffer rows
    model_rate_m_per_yr – shoreline change rate in m/yr

Returns a DataFrame with columns [domain, model_rate] filtered to
DOMAIN_MIN–DOMAIN_MAX, or None if the file is missing.
```

**`domains_to_frac()`**

```
Convert a window size in CASCADE domains to a LOWESS frac value.

Always divides by the total island domain range (DOMAIN_MAX - DOMAIN_MIN + 1)
so that fracs are consistent regardless of how many domains have valid data
in a given CSV.  This ensures:
  5  domains → frac ≈ 0.056
  10 domains → frac ≈ 0.111
  15 domains → frac ≈ 0.167  (matches the LOWESS_FRAC default)
```

**`add_annotations()`**

```
Add all geographic reference annotations to an axis.

Layer order (bottom → top):
  1. Wimble Shoals influence zone  (hatched amber fill, bottom label)
  2. Community shaded spans        (steel-blue fill, top labels)
  3. Village center lines          (dashed gray,  y=0.88)
  4. Pier lines                    (dash-dot blue, y=0.76, rotated)
  5. Groin lines                   (dotted red,    y=0.76, rotated)

All label y-positions use blended axes-fraction coordinates so they
stay fixed relative to the panel height regardless of data range.
```

**`annotation_legend_handles()`**

```
Return proxy artists explaining the annotation layer types.
Append these to a plot's legend handle list so readers can decode
all reference marks without scanning every label individually.
```

**`plot_smoothed_only()`**

```
2-panel: LOWESS smoothed lines only, no raw data.
Cleanest version for presentations or dissertation figures.
```

**`plot_combined_periods()`**

```
Single panel: both periods of smoothed CoastSat on one axis.
Good for directly comparing the two periods.
```

**`plot_smoothing_sensitivity()`**

```
3-panel showing the effect of each LOWESS window on CoastSat data.
One smoothed line per panel so individual window behavior is clear.
Fracs are derived from COMPARE_WINDOWS_DOMAINS / 90 domains, matching
exactly what plot_window_comparison uses.
```

**`plot_window_comparison()`**

```
2-panel figure (one per period) showing the raw CoastSat LRR (faded)
overlaid by LOWESS-smoothed lines for each window size in
COMPARE_WINDOWS_DOMAINS.  All smoothed lines share one panel so spatial
patterns and differences between window choices are directly visible.

Window sizes are specified in CASCADE domains and converted to LOWESS fracs
using the actual number of valid domains in each dataset.

Parameters
----------
cs_1984, cs_2004 : pd.DataFrame or None
    Loaded CoastSat data for each period.
out_path : str
    Full path for the comparison PNG.
window_domains : list of int
    Number of CASCADE domains for each smoothing window to compare.
    Defaults to COMPARE_WINDOWS_DOMAINS from CONFIG.
```

**`plot_cascade_vs_lowess()`**

```
One panel per CoastSat period (1984–2004, 2004–2024).

Each panel shows:
  • Raw CoastSat LRR (faded, period color)
  • Three LOWESS-smoothed CoastSat curves (green / orange / purple)
  • CASCADE modeled change rate(s) for that period (thick black line)

CASCADE runs are matched to panels by their 'period' key in CASCADE_RUNS.
If no CASCADE run exists for a period, that panel still shows the CoastSat
smoothed curves alone.

Parameters
----------
cs_1984, cs_2004 : pd.DataFrame or None
    CoastSat data for each period.
cascade_runs : list of dict
    Loaded CASCADE rate DataFrames, each with an extra 'label' and
    'period' key (same structure as CASCADE_RUNS but with 'df' added).
out_path : str
    Full path for the comparison PNG.
window_domains : list of int
    LOWESS window sizes in CASCADE domains (from COMPARE_WINDOWS_DOMAINS).
```

**`plot_cascade_by_window()`**

```
3-panel figure (one per LOWESS window) for a single CoastSat period.

Each panel shows:
  • Raw CoastSat LRR (faded, period color)
  • ONE LOWESS-smoothed CoastSat curve (bold, period color)
  • CASCADE modeled rate (thick black)

This isolates the model-vs-observation comparison for each smoothing
choice so you can judge fit quality without the visual clutter of
seeing all three windows simultaneously.

Panels share the same y-axis limits so differences in smoothing
level — not axis scaling — drive the visual comparison.

Parameters
----------
cs_df : pd.DataFrame or None
    CoastSat data for the target period.
cascade_runs : list of dict
    Loaded CASCADE run dicts (with 'df', 'label', 'period' keys).
    Only runs whose 'period' matches period_label are plotted.
period_label : str
    e.g. "1984–2004".  Used for title and CASCADE run matching.
out_path : str
    Full path for the comparison PNG.
window_domains : list of int
    LOWESS window sizes in CASCADE domains.
```

</details>

