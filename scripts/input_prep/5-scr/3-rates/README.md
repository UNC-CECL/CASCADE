# 3-rates — the model targets

The rate fits the model is graded against. These are **model inputs**: the
hindcast runner reads the LRR table on every run, and the calibrated
source/sink preset is fitted against it. Everything else here is a reading of
the same record taken a different way.

Grouped by source, then by what kind of number it is.

```
coastsat/
    lrr/           the OLS rate fit — the thing runs are actually graded on
    endpoint/      net change between +/-6-month means at the dune-line dates
    5yr_bins/      the OLS in successive 5-year bins: WHEN did it change?
    window_convergence/  the OLS on NESTED window families pinned at each
                   end: WHICH WINDOWS recover the long-term rate?
                   1-rate_profiles/ (whole profile) and 2-settling_window/
                   (each location)
    total_change/  the rate as a distance, and the projections
    extension/     the same fit beyond the 90 surveyed domains
duneline/
    endpoint       end line minus start line
rates_figures.py   one house-style figure per window, for all of the above
```

## A window is an interval, and the ends are not equivalent

Folders are named `<start>_<end>` because a rate fit spans an interval (rule 2
of `ORGANIZATION.md`). The five windows are **not** peers: `1996 -> 2010 ->
2024` is the main chain, `1984_2004` and `2004_2024` are the older one, and
`1996_2024` is context that nothing is graded against.
`data/hatteras_init/5-scr/WINDOWS.md` says which is which, generated from
`hat_observed_rates.WINDOW_ROLE` by `../tools/windows_index.py` so the
table cannot drift from the code.

## Why the figure script is not inside a product folder

`rates_figures.py` draws the figure for *every* product here, beside its
tables. Putting it in any one product folder would misfile it four times over.
Run it after whichever product you rebuilt. It replaced the autoscaled
quick-looks the LRR fit used to draw itself, which are archived under
`5-scr/archive/coastsat_lrr_quicklooks_20260918/`.

Two figures stay at a window folder's top level with `lrr_<w>.png` —
`smoothing_windows_<w>.png`, from
`coastsat/lrr/coastsat_lrr_smoothing_windows.py`. By contrast
`coastsat_lrr_transect_zoom.py` writes into a `transects/` subfolder, because
that family grows a file per reach and per variant.

## The scripts

```
coastsat/lrr/coastsat_domain_lrr.py
    The fit. One window per run; reads the lookup from ../2-transect-frame/
    and the chainage from 1-observations/coastsat_timeseries/. Tables only
    since 2026-09-18 — the figure is drawn by rates_figures.py.
    (Was coastsat_domain_lrr_fixed.py until 2026-09-22; the "_fixed" had no
    unfixed sibling left to distinguish it from.)

coastsat/lrr/coastsat_lrr_windows.py
    Every window side by side on ONE y axis. Each window's own figure is
    autoscaled, so four of them together draw a 2 m/yr swing as tall as a
    7 m/yr one and the eye reads the wrong story. Also the module the other
    figure scripts import for the shoals, fills and structure drawing.

coastsat/lrr/coastsat_lrr_smoothing_windows.py
    The unsmoothed field plus the LOWESS at 3, 5 and 10 domains (1.5, 2.5 and
    5.0 km) on one axis, darkest being the window runs are graded at.

coastsat/lrr/coastsat_lrr_transect_zoom.py
    One window at transect resolution over a short reach, transects on the
    x axis instead of domains.

coastsat/endpoint/coastsat_endpoint.py
    The mean CoastSat position within +/-6 months of each dune-line survey
    date, end minus start, in metres and m/yr. Centred on the dune dates so
    that a gap between the shoreline and the dune line is a real difference
    and not an artefact of comparing different moments.

coastsat/5yr_bins/coastsat_5yr_bins.py         the table
coastsat/5yr_bins/coastsat_5yr_bins_figure.py  the figure
    Bin edges snap to whole years and bins under 3.75 yr are dropped, so a
    short tail bin cannot manufacture a large rate.

coastsat/window_convergence/coastsat_window_convergence.py
    Which windows recover the long-term rate? The same OLS on two families of
    NESTED windows, one pinned at each end -- forward walks the END out from
    1996, backward walks the START back from 2024 -- so the pair brackets the
    answer rather than giving one side of it. Both converge on the 1996-2024
    rate BECAUSE that rate is the longest window of each family, so what the
    sweep gives is the window at which the curve stops leaving a tolerance,
    not a match test. Three tolerances are scored side by side rather than one
    being chosen -- but the HEADLINE, what the figures draw and the READMEs
    lead with, is CI OVERLAP: a window passes when its own 95% interval reaches
    the reference's. It is the only one that is not a number somebody picked,
    and it is the only one shaped right -- a funnel, wide where the record is
    short, narrowing onto the reference -- because it asks a five-year fit to
    be indistinguishable from the long-term rate rather than as precise as it.
    It barely moves the answer (27 -> 26 years forward, 25 -> 22 backward), and
    that is the finding: the windows disagree because the shoreline changed,
    not because the fits were noisy. Three scales: `a-eight_sites/` draws eight transects in full,
    `b-every_transect/` runs all ~906 and profiles them alongshore, and
    `c-domain_means/` aggregates that sweep to the unit the model is actually
    graded on -- a groupby, not a refit, because the target is the MEAN of
    transect slopes. Averaging turned out to buy nothing: the domain mean
    settles at the same median window as a single transect, which says the
    window-to-window disagreement is shoreline behaviour and not per-transect
    noise. It draws its own figures -- per-unit panel families, not the
    alongshore profile rates_figures.py covers.
    Output: `window_convergence/2-settling_window/<direction>_from_<year>/`;
    `--ref-end 2020` files under `experiments/record_cut_2020/` (that run was
    deleted 2026-09-29; the option still works).

coastsat/window_convergence/coastsat_window_profiles.py
    The same nested families drawn as whole alongshore profiles, from two-year
    windows, over the 1996-2024 rate, each scored by the alongshore Pearson r
    against it (shape only). Reuses the settling sweep's per-transect fit.
    Output: `window_convergence/1-rate_profiles/<direction>_from_<year>/`.

    WHY the windows disagree is answered next door, in
    1-observations/detrended_position/: one island-wide excursion, which on its own
    biases the fitted rate by -0.37 m/yr over 1996-2010 and +0.88 m/yr over
    2010-2024.

coastsat/total_change/coastsat_total_change.py
    A rate turned into a distance, named for the window it was FITTED on:
      --product total      LRR(W) x the years of W      -> total_change/
      --product projected  the 1996-2024 LRR onto a window it was not fitted
                           on                           -> projected/
    The distinction is the whole point of the split; see the docstring.

coastsat/extension/coastsat_extension_lrr.py
    The same fit for the coast beyond GIS 90, for the Pea Island extension
    experiment. Numbers those transects through hat_extension_domains rather
    than the 90-polygon lookup, which stops at 90.

duneline/duneline_endpoint.py
    End dune line minus start line, per 100 m transect and per domain.
    Endpoint and not LRR on purpose: the quantity is net change, and an OLS
    through the intermediate lines would answer a different question. It
    replaced a dune-line LRR product on 2026-09-18.
```

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### coastsat/5yr_bins/coastsat_5yr_bins.py

The CoastSat LRR in successive 5-year bins, per GIS domain, for each window of the canonical chain.

From the script's original header:

```text
CoastSat Shoreline Change Rate — Discrete Interval Line Plots
Divides each analysis period into non-overlapping sequential bins of
fixed width, computes LRR (m/yr) per domain per bin, and plots all bins
as overlaid spatial profiles.

Key quality controls (fixes for inflated rate values):
  1. Bin edges are snapped to whole years — no tiny partial tail bins.
  2. Bins shorter than MIN_BIN_FRACTION (default 0.75) of the interval
     width are dropped entirely.
  3. LRR is computed via compute_lrr() from coastsat_lrr.py,
     the same function used in all other CoastSat scripts.
  4. Per-transect results are filtered by p-value and R² before
     domain aggregation — unreliable fits are excluded.
  5. Domain aggregation uses the MEDIAN by default — one bad transect
     cannot dominate the domain value.
  6. A physical hard cap (MAX_RATE_M_YR) rejects per-transect LRR
     values that are physically implausible before aggregation.

Expected line counts per figure:
  Period 1 (1984-2004, 20 yr)  |  5-yr intervals ->  4 lines
  Period 1 (1984-2004, 20 yr)  |  3-yr intervals ->  6 lines
  Period 1 (1984-2004, 20 yr)  |  1-yr intervals -> 20 lines
  Full     (1984-2024, 40 yr)  |  5-yr intervals ->  8 lines

Inputs
    transect_domain_lookup.csv    (from coastsat_domain_mapping.py)
    CoastSat time-series CSVs     (one per transect, standard format)
    coastsat_lrr.py      (same directory or PYTHONPATH)

Outputs  (OUTPUT_DIR / <window>/, since 2026-09-18; tables only)
    1996_2010/lrr_bins_5yr.csv   rows = bins, cols = GIS domains, LRR m/yr
    2010_2024/lrr_bins_5yr.csv
    1996_2024/lrr_bins_5yr.csv
    Figures: coastsat_5yr_bins_figure.py -> 4-comparisons/coastsat_5yr_bins/

Usage
    Edit CONFIG section, then:
        python coastsat_lrr_interval_lines.py
```

Notes that were in the code:

```text
Resolved through hat_observed_rates.py (2026-09-18); the typed path named
coastsat_lrr/ and used _PATH_REPO before it was defined.
```

```text
Anchored on this file 2026-09-12. The literals here were
drive-rooted and had never resolved; the data they name also
moved out of the scripts tree on that date.
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
--- Time periods: (file_tag, figure_title, file_stem, start_date, end_date) ---
The canonical chain since 2026-09-17 (1996 -> 2010 -> 2024), rebuilt here
2026-09-18; the 1984-2004 / 2004-2024 / 1984-2024 run of 2026-06-02 is in
5-scr/archive/coastsat_5yr_bins/. One folder per WINDOW, the coastsat lrr
naming, calendar years inclusive like every other window.
```

```text
--- Interval sizes (years): one subfolder per entry ---
On a dynamic barrier island like Hatteras, 5-yr is the minimum window
that starts to average through storm-recovery cycles. 1-yr captures
mostly interannual noise; 3-yr is marginal. Recommended: [5] or [3, 5].
```

```text
--- Minimum observations per transect per bin ---
Only gate that must be passed to compute LRR for a transect.
Lower values include more transects; raise to demand denser coverage.
```

```text
--- Statistical quality filters (applied before domain aggregation) ---
Set MAX_PVALUE = 1.0 and MIN_R2 = 0.0 to disable (recommended default).
Short bins (1-3 yr) rarely achieve p <= 0.10 with only 4-8 observations
even when the underlying signal is real — these filters will silently
empty most domains and should be left off unless you have a specific
reason to restrict to high-confidence fits only.
```

```text
--- Physical plausibility cap (m/yr) ---
Per-transect LRR values with |LRR| > MAX_RATE are excluded before
domain aggregation.
```

```text
--- Minimum bin width relative to the nominal interval ---
A bin covering less than this fraction of the full interval is dropped.
Example: with 5-yr intervals, 0.75 drops bins shorter than 3.75 years.
This prevents partial tail bins (e.g., "2024-2024*") from polluting results.
```

```text
--- Domain aggregation method: "mean" or "median" ---
"mean" matches coastsat_domain_lrr_fixed.py and is the correct default
for consistency across the pipeline. The per-transect LRR is already
computed via a full linear regression on all observations in the bin,
so the domain value is the mean of those individual transect slopes.
```

```text
--- Color palette ---
"custom" with the list below gives maximally distinct colors for 4-8 bins,
running cool (oldest) -> warm (most recent).
```

```text
--- Y-axis range for LRR panel ---
None = auto (98th percentile clip). Or fix e.g. (-6, 6) for cross-period comparison.
```

```text
--- Community span annotations (kept for backward compatibility) ---
The full annotation system below supersedes these — leave as-is.
```

```text
SECTION 4b: GEOGRAPHIC ANNOTATION STYLING
Matches the style used across all CoastSat figures in this project.
```

```text
--- LOWESS spatial smoothing overlay ---
LOWESS_OVERLAY = True  : smoothed line drawn on top of the raw line.
LOWESS_ONLY    = True  : raw line is drawn faintly (alpha * RAW_ALPHA_SCALE)
so the smoothed signal is the dominant visual.
Set both True to suppress nearly all raw noise.
LOWESS_FRAC    : fraction of domains used for each local fit.
7-domain window over 90 domains -> frac = 7/90 ≈ 0.078 (10/90 until 2026-09-28).
Matches the window used in the cross-period LOWESS comparison
and preserves community-scale signals (Avon, Wimble Shoals)
while filtering sub-kilometer noise.
Increase toward 0.20 to smooth more aggressively.
```

```text
--- Physical plausibility cap (m/yr) ---
20 m/yr passes genuine extreme signals at Oregon Inlet margin (domain 3)
and the Rodanthe / post-Isabel zone (domains 76-77) while still rejecting
clear CoastSat detection errors.
```

```text
Use the identical compute_lrr as every other CoastSat script in the project.
The companion file sits beside this one (script_dir is on the path above);
the package path it was imported by, scripts.input_preperation..., has not
existed since the tree was renamed, so this script could not run (2026-09-18).
```

```text
Label shows the INCLUSIVE year range: end year minus 1 makes
clear that e.g. "1984-1988" ends at Dec 31 1988 and
"1989-1993" starts at Jan 1 1989 — no overlap, no gap.
```

```text
Auto Y-limit: use the TRUE data range (not percentile clipping)
so all values are visible. A small padding keeps lines off the edges.
```

```text
--- Raw line ---
When LOWESS_ONLY is active, draw the raw line at a much reduced alpha
so the smoothed signal reads as the primary trace.  Setting
RAW_ALPHA_SCALE = 0.0 hides the raw line entirely.
```

```text
One folder per window, overwritten in place (deterministic); no
timestamped run folders and NO FIGURES since 2026-09-18 -- 3-rates holds
data only. The figure is drawn in the house style by
coastsat_5yr_bins_figure.py into 4-comparisons/coastsat_5yr_bins/.
```

<details><summary>Function notes (the original docstrings)</summary>

**`build_color_list()`**

```text
Return exactly n colors interpolated through the chosen palette.
Handles built-in named palettes, matplotlib colormaps, and "custom".
```

**`make_bins()`**

```text
Divide a period into non-overlapping sequential bins of interval_yr years.

Bin edges are snapped to whole years (from the integer start year of
the period), which prevents tiny partial tail bins. Any bin shorter
than min_bin_fraction * interval_yr is silently dropped.

Example — period_start="1984-01-01", period_end="2004-12-31",
           interval_yr=5, min_bin_fraction=0.75:
    -> [(1984, 1989, "1984-1989"),
        (1989, 1994, "1989-1994"),
        (1994, 1999, "1994-1999"),
        (1999, 2004, "1999-2004")]   ← 4 clean bins, no partial tail

Returns
-------
list of (bin_start_yr: int, bin_end_yr: int, label: str)
```

**`compute_domain_bin_lrr()`**

```text
Build a (bins x domains) LRR matrix.

For each bin and each CASCADE domain:
  1. Collect all CoastSat observations within the bin date range.
  2. Compute LRR using compute_lrr() from coastsat_lrr.py
     (same function as coastsat_custom_range_dates_plot.py and all
      other scripts in this project).
  3. Filter transect results by min_obs, max_pvalue, min_r2, and max_rate.
  4. Aggregate remaining valid transects by median (or mean).

Parameters
----------
all_data      : {transect_id: DataFrame}
lookup        : DataFrame with 'transect_id' and 'domain_number'
bins          : list of (start_yr, end_yr, label) from make_bins()
period_start  : ISO date string — data restriction for the whole period
period_end    : ISO date string
min_obs       : minimum obs per transect per bin
max_pvalue    : maximum p-value to accept a transect fit as valid
min_r2        : minimum R² to accept a transect fit as valid
max_rate      : physical cap — |LRR| > this is excluded (noise/error)
agg           : "median" or "mean" for domain aggregation
buffer_domains: domain numbers to exclude

Returns
-------
lrr_df : DataFrame, index=bin_label, columns=domain numbers
```

**`compute_overall_lrr()`**

```text
Single full-period LRR per domain for the bottom summary bar panel.
Uses the same quality filters as the per-bin computation.
```

**`lowess_smooth()`**

```text
Apply LOWESS smoothing along the domain (spatial) axis.
Only fits on non-NaN points; NaN gaps remain NaN in comparison.
```

**`draw_annotations()`**

```text
Apply the full geographic annotation suite to an axes object,
matching the Section 4b style used across all CoastSat figures.

Parameters
----------
ax           : matplotlib Axes to annotate
domains      : np.ndarray of domain numbers present in the figure
show_labels  : if False, draw bands/lines but suppress all text.
               Use False for the bottom panel to avoid collision
               with the panel's own set_title() label.
```

**`plot_interval_lines()`**

```text
Two-panel figure for one period and interval size.

Top    : LRR (m/yr) by domain, one line per time bin. Legend shows
         the date range of each bin (e.g. "1989-1994").
Bottom : Overall LRR across the full period as a reference bar chart.
```

</details>

### coastsat/5yr_bins/coastsat_5yr_bins_figure.py

When inside a window did the shoreline change? The CoastSat LRR in successive 5-year bins, one panel per bin.

From the script's original header:

```text
When, inside a window, did the shoreline change? The CoastSat LRR in
successive 5-year bins, one panel per bin, per GIS domain, in the house style.
Written 2026-09-18 when the bins were rebuilt on the 1996 -> 2010 -> 2024
chain; it writes beside the table, as every 3-rates product does.

READS    3-rates/coastsat/5yr_bins/<window>/lrr_bins_5yr.csv
         (coastsat_5yr_bins.py: per transect an OLS over the bin's
         positions, domain MEAN, |rate| > 50 m/yr dropped, bins under 3.75 yr
         dropped)
DRAWS    one stacked panel per bin, filled blue where the shoreline moved
         seaward and red where it moved landward (the coastsat_lrr_windows
         panel drawing, imported); ONE y axis for all three figures, the largest
         |rate| over GIS 2-90 plus 1 m rounded up, with GIS 1 clipped and
         labelled where it runs past (Hannah, 2026-09-18). Village bands, the groin and
         piers, the offshore shoals as faint hatched boxes; a model-input
         beach fill is marked above the panel of the bin it falls in.
WRITES   3-rates/coastsat/5yr_bins/<window>/lrr_5yr_bins_<window>.png
         (+ PDF and CAPTIONS.md under supporting/)

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/5yr_bins/coastsat_5yr_bins_figure.py
```

Notes that were in the code:

```text
Beside its table since 2026-09-18 (every 3-rates product carries its own
figure; Hannah). It sat in 4-comparisons/coastsat_5yr_bins/ for an hour.
```

```text
THE AXIS IS CAPPED (Hannah, 2026-09-18). Cape Point (GIS 1) reaches about
+35 m/yr in the 2020s bins -- the shoal attaching -- and at full scale that
one domain set a +/-36 axis that flattened every other bin, all of which stay
within +/-13. The bound is taken over GIS 2-90 of EVERY window, so the three
figures share it; a value beyond it is clipped at the edge, marked with a
triangle and labelled with its value, never dropped.
```

```text
draw_panel ticks every 2 m/yr, which is unreadable once one bin's
spike (Cape Point, GIS 1, 2020-2024) sets a +/-36 axis
```

<details><summary>Function notes (the original docstrings)</summary>

**`mark_clipped()`**

```text
A triangle at the axis edge for every domain beyond it, labelled with
the domain and its value. Returns [(gis, value), ...].
```

</details>

### coastsat/endpoint/coastsat_endpoint.py

CoastSat net shoreline change per window: the mean position around the end survey minus that around the start.

From the script's original header:

```text
The stored CoastSat net shoreline change: for each CoastSat transect, the mean
shoreline position in a +/-6-month window about the window's END dune-line
survey date minus the mean about its START survey date, in METRES and m/yr,
per transect and per GIS domain. Built 2026-09-18 (Hannah, by interview) as
the counterpart of 3-rates/duneline/endpoint/, so the shoreline and the dune
line difference like for like: same windows, same survey dates, same sign.

WHY THE DUNE DATES
    A dune line is a survey at a moment. Centring the CoastSat windows on the
    same moments (1997-10-12, 2009-05-30, 2023-07-01 assumed; the 1984 and
    2004 lines for the older windows) means a gap between the two changes is
    beach-width change, not a date mismatch. Each end averages one full
    seasonal cycle of satellite positions. The window means are the ones
    coastsat_vs_duneline.py uses (window_mean / endpoint_by_transect,
    imported, not copied), so its CoastSat endpoint and this agree.

    change_m = end mean - start mean. CoastSat chainage grows SEAWARD, so
    SEAWARD IS POSITIVE, as in the dune product. rate_m_yr divides by the
    survey interval and inherits any assumed date; change_m does not depend on
    the interval (the windows are still centred on the assumed date).

OUTPUT   data/hatteras_init/5-scr/3-rates/coastsat/endpoint/<start>_<end>/
    transect_endpoint.csv         per CoastSat transect: both window means,
                                  how many positions fell in each, the first
                                  and last date inside each, change_m,
                                  rate_m_yr, the survey dates
    domain_endpoint_summary.csv   per domain: n, mean/std/min/max change_m,
                                  mean rate, pct_landward, the median
                                  positions per end window, and how many of
                                  its transects had an empty end window
    PROVENANCE.md
    Read through hat_observed_rates.coastsat_endpoint_csv(start, end, level).

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/endpoint/coastsat_endpoint.py
    python scripts/input_prep/5-scr/3-rates/coastsat/endpoint/coastsat_endpoint.py --windows 1996_2024
```

### coastsat/extension/coastsat_extension_lrr.py

CoastSat rates for the coast beyond the 90 surveyed domains, for the Pea Island extension experiment.

From the script's original header:

```text
CoastSat rates for the coast BEYOND the 90 surveyed domains
The Pea Island extension experiment (2026-09-16) models GIS 1-115 (and
0-115) instead of 1-90, and solves the end-domain source/sink at the new
ends against CoastSat, as the matrix does at GIS 1 and 90. The committed
rate tables stop at GIS 90 only because the transect-to-domain lookup was a
polygon join onto 90 polygons; the CoastSat record itself runs to Oregon
Inlet (263 transects on disk between GIS 90 and 115, ~10 per domain, a
median 260 observations each over 1996-2010).

This script numbers those transects the way the surveyed ones were
numbered -- the transect's origin point within a 500 m domain polygon, now
Hannah's whole-island polygons (hat_extension_domains.join_origins) -- and
fits them with the same LRR as coastsat_domain_lrr_fixed.py, for one window.
A transect no polygon covers is left out, as the surveyed mapping leaves
them out. It writes, beside the surveyed products and never into them:

    5-scr/2-transect-frame/transect_domains/transect_domain_lookup_ext.csv
    5-scr/3-rates/coastsat/lrr/<start>_<end>/ext/transect_lrr_full.csv
    5-scr/3-rates/coastsat/lrr/<start>_<end>/ext/domain_lrr_summary.csv
    5-scr/3-rates/coastsat/lrr/<start>_<end>/ext/transect_lrr_with_base.csv

The last is the surveyed table with the extension rows appended: what an
extended-geometry run loads as its active dataset. The window's own
transect_lrr_full.csv must already exist (coastsat_domain_lrr_fixed.py).

Usage
    python coastsat_extension_lrr.py --start-year 1996 --end-year 2010
```

Notes that were in the code:

```text
The origin point within a polygon, the rule of
coastsat_domain_mapping.py, onto Hannah's whole-island polygons
(2026-09-16). A transect no polygon covers is left out.
```

<details><summary>Function notes (the original docstrings)</summary>

**`fit_window()`**

```text
coastsat_domain_lrr_fixed.compute_all_lrr, for the extension rows.
That script parses its arguments at import, so the ten lines are
repeated here rather than imported.
```

</details>

### coastsat/lrr/coastsat_domain_lrr.py

Step 2 of 2: per-domain CoastSat LRR for one window, from the transect-to-domain lookup.

From the script's original header:

```text
CoastSat Domain-Level LRR Summary
Step 2 of 2 in the domain-level LRR workflow.

Requires:
  1. transect_domain_lookup.csv   – from coastsat_domain_mapping.py
  2. CoastSat time-series CSVs    – one per transect
  3. coastsat_lrr.py     – in the same directory (or on PYTHONPATH)

Outputs (tables only since 2026-09-18; the window's figure is drawn beside
them by scripts/input_prep/5-scr/3-rates/rates_figures.py):
  domain_lrr_summary.csv  –  one row per domain with aggregated LRR stats
  transect_lrr_full.csv   –  full transect-level results with domain assignments

Usage
Edit the CONFIG section below, then run:
    python coastsat_domain_lrr.py
```

Notes that were in the code:

```text
Path to the lookup table produced by coastsat_domain_mapping.py
THE WINDOW IS GIVEN AS A PERIOD, and every path is anchored on this file
(2026-09-11). The three literals here were a machine-specific absolute
path into "input_preperation", a tree that no longer exists, plus two
drive-rooted paths from before the 5-scr rename -- none of them resolved.

python coastsat_domain_lrr_fixed.py --start-year 1996 --end-year 2010
```

```text
Root folder containing all site subfolders (e.g. usa_NC_0032_timeseries, usa_NC_0033_timeseries, ...)
The script will automatically find every CSV in every subfolder one level down.
Example: r"C:/Users/hahenry/Downloads"
```

```text
Optional: only include subfolders whose names contain this string.
Set to "" to include ALL subfolders under ROOT_DATA_DIR.
```

```text
CASCADE buffer domains to EXCLUDE from summaries
(e.g., the 15 buffer domains on each end of your 90+30 setup)
Set to an empty list [] to include all domains
```

```text
---- No figures here (2026-09-18) ----
The figure is drawn by scripts/input_prep/5-scr/3-rates/rates_figures.py. The
domain_lrr_bar.png / transect_lrr_scatter.png quick-looks this used to
draw were autoscaled, titled and coloured by magnitude, a second picture
of the numbers that clashed with the house-style figures; they are in
5-scr/archive/coastsat_lrr_quicklooks_20260918/. The window's figure is
python scripts/input_prep/5-scr/3-rates/rates_figures.py
-> 3-rates/coastsat/lrr/<window>/lrr_<window>.png (since 2026-09-19). plot_domain_lrr and
plot_transect_scatter are kept, unused, for a one-off look.
```

<details><summary>Function notes (the original docstrings)</summary>

**`collect_csv_map()`**

```text
Auto-discover all time-series CSVs under root_dir.

Walks one level of subfolders (e.g. usa_NC_0032_timeseries/) and
collects every CSV inside them. Optionally filters to subfolders
whose names contain site_filter.

Returns:
    { transect_id_stem : full_filepath }
    e.g. { 'usa_NC_0032_0011' : 'C:/Downloads/usa_NC_0032_timeseries/usa_NC_0032_0011.csv' }
```

**`compute_all_lrr()`**

```text
For every transect in the lookup table, find its CSV, compute LRR,
and return a merged DataFrame with domain assignments.
```

**`domain_summary()`**

```text
Aggregate transect-level LRR results to domain level.

Returns one row per domain with:
  n_transects    – total transects in domain
  n_valid        – transects with valid LRR
  mean_lrr       – mean LRR across transects (m/yr)
  median_lrr     – median LRR
  std_lrr        – standard deviation
  min_lrr        – most erosional transect
  max_lrr        – most accretionary transect
  pct_eroding    – % of transects with negative LRR
```

**`plot_domain_lrr()`**

```text
Bar chart of domain-level LRR coloured by magnitude.
metric: 'mean_lrr' or 'median_lrr'
```

**`plot_transect_scatter()`**

```text
Scatter of individual transect LRRs coloured by domain,
sorted by domain number.  Good for seeing within-domain spread.
```

</details>

### coastsat/lrr/coastsat_lrr_smoothing_windows.py

The three LOWESS windows overlaid on one LRR rate field, in the house style.

From the script's original header:

```text
The three LOWESS windows overlaid on ONE LRR rate field, in the house style
(Hannah, 2026-09-21). The rate is the thing smoothed here -- this is the field
itself, not a projection into metres and not a model comparison.

    coastsat/lrr/<w>/smoothing_windows_<w>.png

WHAT IT SHOWS
    The unsmoothed per-domain mean rate as the palest line, then the LOWESS
    curve at 3, 5 and 10 domains (1.5, 2.5, 5.0 km) on a light-to-dark ramp,
    darkest being the window every run is actually graded at. Village bands,
    groin and piers, the hatched shoal boxes and the model-input beach fills
    come from the coastsat_lrr_windows panel, so this figure reads directly
    against lrr_<w>.png beside it -- same y bound, same marks.

WHY THE SIGN FILL IS NOT HERE
    lrr_<w>.png colours the rate blue seaward / red landward. Four curves
    share this panel, so that pair is not available: the smoothing width is an
    ORDERED variable and takes a sequential ramp instead, the one
    smoothing_scale.py uses, anchored on the house shoreline blue. Sign is
    read off the zero line, which is drawn.

THE SPLICE
    Every curve keeps the raw domain means over GIS 1-10 (the Oregon Inlet
    boundary treatment, cascade_pipeline.coastsat_lowess.DEFAULT_LOWESS), so all
    four are identical there by construction -- the same splice the scoring
    target is built through.

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_smoothing_windows.py
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_smoothing_windows.py --window 1996_2024
```

Notes that were in the code:

```text
Domain units. 0 is the unsmoothed domain means; 10 is the grading window
(cascade_pipeline.hindcast TARGET_WINDOW). 3 is roughly the alongshore
decorrelation scale of the domain-mean rate, 1.5 km -- the set
smoothing_scale/ sweeps (Hannah, 2026-09-21).
```

```text
The two shoal-fronted peaks of the full-period field, the ones the widest
window flattens. Read off the figure in smoothing_scale/PROVENANCE.md; the
caption states the share of removed signal that actually falls in them
rather than asserting the attribution.
```

```text
The transect cloud the curves are fitted to (Hannah, 2026-09-22): the same
size and alpha as the dots on lrr_<w>.png, but ONE neutral grey rather than
the sign pair -- sign is the other figure's variable, smoothing width is
this one's, and a second colour scale on the same panel reads as a third.
```

```text
Does the widest window take its structure out of the shoal-fronted
peaks, or evenly along the island? Squared removed signal per domain,
peaks against the rest, over the unspliced reach only.
```

```text
mean_lrr all-NaN draws the frame, grid, village bands and structures with
no sign fill: the ramp below carries the ordered variable instead.
```

<details><summary>Function notes (the original docstrings)</summary>

**`win_label()`**

```text
One smoothing width, as the legend says it: the physical width alone.

Which width the runs are graded at is a caption matter, not a legend one
-- spelling it here ran the four entries past the figure edge. The domain
count went the same way when the transect cloud took a fifth entry
(2026-09-22); the caption gives both units for every width.
```

**`transects()`**

```text
(frame, domain ids, along-coast metres, rate) for one LRR window.

along-coast metres follow the convention every target build uses: each
domain's transects spread evenly across its 500 m band, ordered within the
domain by transect_id.
```

**`dots()`**

```text
The individual transect rates under the curves. Returns how many fell
outside the shared y bound, drawn as open markers at the edge as they are
on lrr_<w>.png.
```

**`structure()`**

```text
What each width takes out of THIS field, for the caption.

Everything is in m/yr of alongshore structure removed, never a share of
variance: the domain-mean variance is dominated by the long-wavelength
swings, so a wiggle that is plainly visible on the figure reads as a few
per cent of it and the percentage badly undersells the effect (Hannah,
2026-09-21). These reproduce the target row of
output/comparisons/target_comparison/smoothing_scale/tables/
target_structure.csv exactly, and are recomputed here so the caption
cannot drift from the table.
```

**`lrr_half()`**

```text
The y bound of every LRR window figure, so this one matches the figure
already beside it (rates_figures.lrr_figures).
```

</details>

### coastsat/lrr/coastsat_lrr_standalone.py

Shoreline change rate (LRR) for every CoastSat transect in a folder.

From the script's original header:

```text
Shoreline change rate (LRR) for every CoastSat transect in a folder.

Input:  CoastSat time-series CSVs, one per transect (columns "dates UTC",
        "chainage (m)"), anywhere under FOLDER.
Output: CSV of transect_id, lrr_m_yr, unc_m_yr (95% CI half-width), n_obs.

Rate = ordinary least-squares slope of position against time, using every
position from 1 Jan START to 31 Dec END. No filtering, no weighting; at
least 3 positions. Positive = seaward (accretion), negative = erosion.

    python coastsat_lrr_standalone.py FOLDER START END OUT.csv
    python coastsat_lrr_standalone.py coastsat_timeseries 1984 2025 rates.csv

START and END are whole calendar years, both included -- the window
convention of every 5-scr script and of 5-scr/template/, whose
shoreline_rates_template.py is the fuller version of this file.
```

### coastsat/lrr/coastsat_lrr_transect_zoom.py

One window's LRR at transect resolution over a short reach, transects on the x axis.

From the script's original header:

```text
One window's LRR at TRANSECT resolution over a short reach, with the transect
itself on the x axis (Hannah, 2026-09-22: "x being the transects instead of
the domains").

    coastsat/lrr/<w>/transects/lrr_transects_<w>_gis<lo>-<hi>.png
    coastsat/lrr/<w>/transects/slides/  the --slide versions

    Its own subfolder because the family grows a file per reach and per
    variant, and the window folder's convention is tables plus the one figure
    for the window (2026-09-22). smoothing_windows_<w>.png stays at the top
    level with lrr_<w>.png: coastsat_total_change.py and five PROVENANCE.md
    files cross-reference it by that path.

WHY IT EXISTS
    Every other figure in 3-rates puts the GIS domain on x and shows the
    transects as a cloud behind the domain mean. That is the right axis for
    comparing against the model, whose cell IS the domain, but it hides the
    step the averaging actually takes: which transects fall in which domain,
    how far apart they sit, and how much of the domain mean is one transect.
    This figure is that step, drawn -- the example for explaining the chain,
    not a product any run reads.

WHAT IT SHOWS
    One marker per CoastSat transect, in order south to north, coloured by its
    own sign (the house blue / red pair) with its 95 % regression confidence
    interval as a whisker. The domain each transect belongs to is the shaded
    band behind it, labelled along the top; the flat segment across each band
    is that domain's mean -- the plain average of the markers above it, the
    value domain_lrr_summary.csv holds. That is the whole figure by default.

    --target adds the scoring curve as a second set of segments, written to
    <stem>_with_target.png so it never overwrites the plain version. It is not
    an average of anything inside the band: the LOWESS reads 5 km either way,
    so it can sit off the mean, which is the point of drawing it -- but it is
    a third quantity on the panel and needs explaining before it can be read,
    so the teaching version leaves it off.

THE REACH
    Default GIS 64-72, the Tri-Village gradient: the domain mean crosses zero
    at 64 and reaches +2.8 m/yr at 67, so a real alongshore signal is resolved
    transect by transect rather than a flat stretch where every marker is the
    same number. --gis takes any span.

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_transect_zoom.py
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_transect_zoom.py --window 1996_2024 --gis 28 36
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_transect_zoom.py --target
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_transect_zoom.py --slide
```

Notes that were in the code:

```text
--slide: the same figure on a canvas that fits a slide beside text. The
house font sizes are absolute, so a smaller canvas makes the type
relatively LARGER, which is what a projected figure needs; only the marker
size, the tick spacing and the legend columns have to come down with it
(Hannah, 2026-09-22 -- same layout, just smaller).
```

```text
A slide panel is a third the width of the page figure, so the same nine
domains would put 90 markers in 3.4 in and the gradient would read as a
smear. --slide therefore narrows the REACH as well as the canvas, to four
domains across the peak, unless --gis says otherwise (Hannah, 2026-09-22).
```

```text
The along-coast axis covers anything from a 4 km reach to the whole 45 km
island, so neither its unit nor its tick step can be a constant: 1 km ticks
over the island put 45 labels on top of each other.
```

```text
Off by default (Hannah, 2026-09-22): this figure's job is the transect
-> domain step, and the scoring curve is a third thing on the panel
that has to be explained before it can be read. --target puts it back,
under its own name, so the two versions never overwrite each other.
```

```text
Every other domain shaded, so a band is a domain without a legend entry.
Past ~20 domains the bands are a comb and the numbers collide, so a long
reach gets neither and reads as the plain scatter it is.
```

```text
The per-domain values, each flat across the domain it belongs to. The
target takes a mid-ramp blue when the LOWESS curve is drawn too, so the
curve and the domain-resolution version of it are not one colour.
```

```text
No bands to align to, so the panel holds the markers and a little
air, not the 500 m boundaries they happen to sit between.
```

```text
Bounds over EVERY series drawn, not just the markers: the LOWESS reads
5 km beyond the reach and routinely sits below everything in it, so
bounds taken from the markers alone would clip the curve.
```

```text
The page figure always holds zero -- the line sign is read against is
worth the white space in a document. A slide panel is a third the width
and cannot spare it, so it crops to the data and the caption says when
zero fell off (Hannah, 2026-09-22).
```

```text
The page label does not fit a 3.4 in panel -- it ran off both edges --
so the slide keeps only what the axis cannot be read without, and the
shading and the domain numbers are left to the caption (2026-09-22).
```

```text
No domain vocabulary on the panel at all, not even the reach: the
figure is the transect fits and nothing else, and where on the
island it sits is the caption's job (Hannah, 2026-09-22).
```

<details><summary>Function notes (the original docstrings)</summary>

**`load()`**

```text
Transect rows inside GIS lo..hi, south to north, with an x position.

x is the transect's ORDER in the reach, not a distance: one unit per
transect is what puts the transect on the axis. The transects are ~50 m
apart and evenly spread inside a domain, so the two differ only where a
domain holds an unusual number of them, and the band edges below are drawn
from the counts rather than assumed.

With `lowess_window`, a `lowess` column carries the smoother's value at each
transect. It is fitted over the WHOLE island first and sliced to the reach
afterwards, never fitted to the reach alone: a LOWESS reads 5 km either
way, so a fit stopping at the reach edge would be a different curve from
the one the target is built through.
```

**`along_axis()`**

```text
(scale, unit, tick step) for an along-coast axis spanning `span_m`.

scale divides the metre values for plotting, so the numbers on the page
stay short; the coordinate itself is unchanged.
```

**`bands()`**

```text
[(domain, left, right)] per domain, in the coordinate on the x axis.

In metres the edges are the domain's true 500 m boundaries. On the
transect-order axis there is no such thing -- a domain is however many
transects fell in it -- so the edges are half a transect either side of
its first and last.
```

</details>

### coastsat/lrr/coastsat_lrr_windows.py

Observed shoreline change rate along the island, one panel per window, every panel on the same y axis.

From the script's original header:

```text
Observed shoreline change rate along the island, one panel per rate window,
every panel on the SAME y axis.

WHY
    Each window under coastsat_lrr/<start>_<end>/ ships its own
    domain_lrr_bar.png, autoscaled to that window. Put four of them side by
    side and a 2 m/yr swing in 1984-2004 is drawn as tall as a 7 m/yr swing in
    2010-2024, so the eye reads the wrong story. These figures pin one y range
    across the windows (Hannah, 2026-09-15: "the y axis among these plots must
    be held to the same bounds so they are easily comparable").

WHAT IS DRAWN
    The per-domain mean LRR (the model's grading target, `mean_lrr` in
    domain_lrr_summary.csv) as a LINE between domain centres over the 90 GIS
    domains (a step at the bin edges was tried 2026-09-15 and Hannah's
    advisor asked for a line). The line is coloured by sign and filled to
    zero: RdBu blue where the shore accreted, red where it eroded - the same
    reading as the per-window domain_lrr_bar.png, and Hannah's call over a
    BrBG pair (2026-09-15). Two dotted lines are +/- one standard deviation
    across the CoastSat transects inside each domain.
    Village spans are the house bands; the Buxton groin and the two piers are
    hairlines from the site config. Nothing else is on the canvas; the
    captions carry the method.

Y BOUNDS
    Symmetric: the largest |mean| over all windows plus a 1 m pad, rounded up
    to the next metre, mirrored about zero, ticks every 2 m/yr. The std lines
    are NOT in the bound - one wide domain at Cape Point was pushing every
    panel to +/-9 with nothing above 7 (the 2026-09-15 first cut). The value
    used is written to supporting/y_bounds.txt.

OUTPUT   data/hatteras_init/5-scr/3-rates/coastsat/lrr/  (since 2026-09-19)
    ONLY the 2 x 2 is drawn now. Until 2026-09-19 this wrote
    4-comparisons/coastsat_windows/ with one folder per window and the
    --overlay halves figure too; those duplicated rates_figures.py's window
    and chain figures and were archived (--overlay is retired). The listing
    below is as it was.
    lrr_four_windows.png           2 x 2: the 1984-start period in the left
                                   column, the 1996-start period in the right
    supporting/
        lrr_windows_wide.csv       the four means and stds side by side
        y_bounds.txt               the bounds every panel uses
        lrr_four_windows.pdf, CAPTIONS.md
    <start>_<end>/lrr_<start>_<end>.png    one figure per window, its PDF and
                                           caption under its own supporting/

    --overlay 1996_2024 (2026-09-18) draws ONLY, into 1996_2024/:
    lrr_1996_2024_halves.png       two stacked panels: (a) the long window
                                   filled by sign, the model-input fill
                                   footprints as bars above it; (b) the two
                                   default windows that chain across it
                                   (1996-2010 grey, 2010-2024 black). No std
                                   lines; the y axis is the tightest whole
                                   metre holding every line, NOT the shared
                                   bound. Also published to
                                   output/figures/2-observations/shoreline/
    supporting/lrr_1996_2024_halves.csv   the three means side by side
    The long window is context, not a grading target, so it is NOT added to
    the default four or to the 2 x 2.

USAGE
    python coastsat_lrr_windows.py                      # the four default windows
    python coastsat_lrr_windows.py --windows 1984_2004 2004_2024
    python coastsat_lrr_windows.py --overlay 1996_2024
```

Notes that were in the code:

```text
The 2 x 2 layout is by model period: each column is one CHAIN of windows,
the second starting where the first ends (1984-2004 then 2004-2024 is the
1984-start period; 1996-2010 then 2010-2024 the 1996-start). Lettered across
then down. Any window set that is not two chains of two falls back to one
column.
```

```text
The RdBu poles and their light fills: red erosion, blue accretion, as the
per-window bar charts already read. In this figure the pair means SIGN, not
vintage; no vintage is drawn here, so the two readings never meet.

The fills are the house light pair lightened by a third toward white
(Hannah, 2026-09-15): at full strength the fill carried more weight than the
lines over it, and on the comparison figure the black model line is what
should read first. Both figures take these constants, so they move together.
```

```text
structures() moved to hat_figure_style 2026-09-15 (shared with
coastsat_vs_duneline); imported above.
```

```text
TWO STACKED PANELS (Hannah, 2026-09-18). The first cut laid the halves over
the long window as dotted and dashed ink lines in one panel; they were more
variable than the long window, crossed it everywhere, and pulled the eye off
the quantity the figure is about. Now (a) is the long window alone, filled
by sign, and (b) the two model periods as plain lines.

The halves are a LUMINANCE pair, not the house vintage pair: panel (a)
already spends red/blue on SIGN, and red "earlier" in (b) directly under red
"landward" in (a) would read as one thing. Lighter grey is the earlier
period, ink the later -- it survives greyscale and colour deficiency.
```

```text
SHOALS (Hannah, 2026-09-18: "lighter so it doesn't take away from the
shoreline change, I just want to see where it is"). Full-height BOXES, a
thin amber outline over a sparse, faint amber hatch, with NO fill: the
other shoal figures' solid wash would tint the sign fill in (a) and merge
with the village greys, while a hatch stays readable over both. A thin
bottom strip was tried first; Hannah asked for hatched boxes instead.
The house amber (C["ADDED"]) the other alongshore figures use for shoals.
```

```text
The single-window figures' bound, named in the caption so a reader
knows this figure's axis differs from theirs.
```

```text
RETIRED 2026-09-19: the halves overlay duplicated the 3-rates chain
figure (lrr/chains/lrr_chain_1996_2010_2024) and was archived.
```

```text
Since 2026-09-19 only the 2 x 2 is drawn, into 3-rates/coastsat/lrr/
(OUT_DIR). The per-window figures are rates_figures.py's; drawing them
here too would overwrite 3-rates/coastsat/lrr/<w>/lrr_<w>.png.
```

<details><summary>Function notes (the original docstrings)</summary>

**`shared_bounds()`**

```text
Half-range: the largest |mean| over every window plus the pad, rounded
up to the next whole metre per year.
```

**`signed_segments()`**

```text
The line as (segment, colour) pairs. A segment that crosses zero is
split where it crosses, so each half carries its own sign and the colour
changes exactly where the fill does.
```

**`draw_panel()`**

```text
The observed panel. `std`, `line_lw` and `fill_y` exist for the
comparison figure that lays a scoring target over this
(scripts/analyze_output/compare_runs/rate_windows.py):
with `fill_y` given, THAT series takes the fill and a light outline, and
the per-domain means are only the thin line over it, so the reference is
the shape and the data the line (Hannah, 2026-09-15, option A).
```

**`_save()`**

```text
PNG in the folder, PDF under supporting/ -- the house save() does that
itself since 2026-09-15 (it put the PDF beside the PNG, and this script
kept a pdf/ folder of its own, before then).
```

**`_chains()`**

```text
Windows linked end-to-start, each chain sorted by start, chains by
their first start: [(1984,2004),(2004,2024)], [(1996,2010),(2010,2024)].
```

**`grid_figure()`**

```text
2 x 2 when the windows form two chains of two (a column per chain,
the earlier window above); one column otherwise.
```

**`fills_in()`**

```text
(year, first_gis, last_gis) of every ENABLED model-input fill placed
inside the window, read from the site config the hindcast reads -- the
model-input footprint, not the wider or narrower record span (Hannah,
2026-09-18).
```

**`draw_fills()`**

```text
A bar just ABOVE the frame over each fill footprint, the year on it.
A mark, not a shade: the village bands already shade. Above rather than
along the bottom because the pier labels stand at the bottom, and Avon's
ran through the 2022 label there (first cut, 2026-09-18).

Placed in AXES fractions, not data units, so the same call works on the
m/yr and the metre figures (duneline_windows.py reuses it); `half` is
kept for the old callers and not used.
```

**`draw_shoals()`**

```text
Each shoal zone of the site config as a hatched, outlined box the
full height of the panel, behind the data, named at the bottom when
`label`. Hatch and outline are two patches because matplotlib 3.9 takes
the hatch colour from the edge colour.
```

**`tight_bound()`**

```text
The smallest whole metre per year that holds every line drawn, no pad.
Hannah asked for a tighter axis than the shared +/-8 (2026-09-18) and
named +/-6; 2010-2024 reaches +6.9 at GIS 1, so a fixed 6 would clip a
measured value. The bound is computed so it can never clip.
```

**`draw_halves()`**

```text
Panel (b): the two periods as solid lines over the same furniture as
(a), unlabelled -- (a) carries the village and structure names.
```

**`_save_both()`**

```text
Into the comparison folder, and published to output/figures/2-observations/shoreline/
(finished figures publish by subject, 2026-09-18). Each copy gets its
CAPTIONS.md entry through the caption() wrapper.
```

</details>

### coastsat/lrr/coastsat_obx_lrr.py

Full-record CoastSat LRR for every transect from Cape Point to the Virginia line: the Murray lab hand-off.

From the script's original header:

```text
Full-record CoastSat LRR for every transect from Cape Point to the Virginia
line: the table handed to the Murray lab (2026-09-27) for their diffusivity
work, plus a README and one figure.

WHY A SIBLING OF coastsat_domain_lrr.py
    That script is driven by transect_domain_lookup.csv, so it only sees the
    transects inside the 90 GIS domains. This hand-off covers 153 km of coast,
    ~90 km of it north of the model domain where no domain exists, so it is driven by the CoastSat
    transect layer instead. The FIT is the same one -- coastsat_lrr.compute_lrr
    on the calendar window, every point, no filter, no weights, >= 3 points --
    so a transect inside the domain gets the number the model target uses.

THE SPEC (decided 2026-09-27)
    window    1984-01-01 -> 2025-12-31 (whole record, last full calendar year)
    extent    usa_NC_0032_0021 (south end of the model domain, GIS 1) through
              usa_NC_0049_0230 (ends at the NC/VA line, 36.550 N)
    columns   the hand-off CSV is what was asked for: transect ID, alongshore
              km, origin lon/lat, LRR, its 95% CI, flag. Every other fit stat
              and CoastSat's beach slope go in supporting/..._full.csv
    flags     flagged, never withheld: fewer than 50 positions, under 20 yr
              between first and last position, within 2 km of Oregon Inlet or
              of the south end at Cape Point (only the two location flags
              fire on this record)

ALONGSHORE DISTANCE
    km from the origin of usa_NC_0032_0021, accumulated transect to transect
    along the coast. Each step between neighbouring origins is projected onto
    the local shore-parallel direction (perpendicular to the two transects'
    mean bearing), because the origins wander cross-shore by up to ~300 m on
    the Currituck Banks and a straight point-to-point sum would add that
    wander as length. Oregon Inlet is the one real gap (~1.1 km).

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_obx_lrr.py

OUTPUT  data/hatteras_init/5-scr/3-rates/coastsat/lrr/1984_2025_obx/
    coastsat_lrr_obx_1984_2025.csv        the hand-off table (7 columns)
    README.md                             written for the recipients
    lrr_obx_1984_2025.png                 rate vs alongshore km, place names,
                                          1 km median line (a guide only)
    supporting/coastsat_lrr_obx_1984_2025_full.csv   all 19 columns
    supporting/ PDF and CAPTIONS.md
    The maps are drawn from the full table by coastsat_obx_lrr_maps.py; a
    one-file version of the fit for others is coastsat_lrr_standalone.py.
```

Notes that were in the code:

```text
The hand-off table is the short one (what was asked for: rate by transect
ID); every fit and CoastSat field is kept in supporting/ (Hannah, 09-27).
```

```text
Places named along the coast, by latitude. The one list for both this
script's profile and coastsat_obx_lrr_maps.py, which imports it.
```

```text
CoastSat's own per-transect beach slope, the one its tidal correction
used, carried through as given (added 2026-09-27 at the recipients' use).
```

```text
compute_lrr rounds p to 6 dp, which printed 2,436 of them as 0.0.
Same x and y as its fit (years since the first position), unrounded.
```

```text
Flagged transects are drawn like the rest (Hannah, 2026-09-27: the hollow
markers came off the figures); the flag lives in the table.
```

```text
Running median over +/-0.5 km of coast, each side of the inlet on its
own. A count window (21 transects) was tried first: it spans well over
1 km where transects are sparse and reached across the inlet.
```

```text
Place names along the top edge, at the transect nearest each latitude;
Oregon Inlet at the gap itself, and the two ends of the reach.
```

```text
every named place gets the same faint dashed line (Oregon Inlet was the
only one until 09-27); the ends of the reach are the axis limits
```

```text
vertical: at 40 degrees the pairs 4 km apart (Cape Point/Buxton,
Carova/NC-VA line) ran into each other
```

```text
Coloured by rate on the maps' scale (Hannah, 2026-09-27: erosion red,
accretion blue, deeper with severity); saturates at +/-V_HALF.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_transects()`**

```text
The CoastSat layer's transects for SITES, south to north, from FIRST_ID.

Order is site then transect number; both run south to north here (checked
2026-09-27: four steps go <= 31 m south, all local wiggles).
```

</details>

### coastsat/lrr/coastsat_obx_lrr_maps.py

Maps of the full-record CoastSat LRR, Cape Point to the Virginia line: an overview and four regional zooms.

From the script's original header:

```text
Maps of the full-record CoastSat LRR, Cape Point to the Virginia line: one
overview and four regional zooms. Reads the table coastsat_obx_lrr.py wrote;
fits nothing.

EACH FIGURE IS TWO PANELS ON ONE NORTHING AXIS
    (a) the map: every transect drawn at its true position and length,
        coloured by rate on a fixed diverging scale (RdBu, +/-3 m/yr), over
        Esri World Shaded Relief recoloured to light grey (contextily, so it
        needs the internet). That basemap carries no labels: every name on
        a map is placed here from PLACES (coastsat_obx_lrr.py) and
        MAP_AREAS. Latitude ticks on the left edge, exact at the coast; the
        NC/VA state line drawn where it is in view; north arrow top right.
    (b) the same rates against northing, sharing (a)'s y axis, so a colour on
        the map reads straight across to its value; the black line is the
        1 km median of the transects, a guide drawn over the points.
    The coast here runs within ~25 degrees of north, which is what makes a
    shared northing axis honest. Cape Point, where it turns, is flagged in
    the table anyway.

THE REGIONS  split at towns, ~35-40 km each, 1 km overlap
    A  Cape Point to Salvo
    B  Rodanthe to South Nags Head (Oregon Inlet in the middle)
    C  Nags Head to Duck
    D  Corolla to the Virginia line

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_obx_lrr_maps.py

OUTPUT  beside the table in 5-scr/3-rates/coastsat/lrr/1984_2025_obx/
    lrr_obx_1984_2025_map_overview.png
    lrr_obx_1984_2025_map_<A-D>_<region>.png
    (+ supporting/ PDFs and CAPTIONS.md entries)
    Reads supporting/coastsat_lrr_obx_1984_2025_full.csv (it needs the
    seaward ends), so run coastsat_obx_lrr.py first.
```

Notes that were in the code:

```text
Label-free (09-27): Esri's Gray Canvas printed "BODIE ISLAND" beside Kitty
Hawk, its historical name for the whole Nags Head-Duck barrier, which read
as a clash with the Bodie Island spit at Oregon Inlet. Every name on the
maps is now one this script places. Esri's shaded relief carries no labels;
it is turned to light grey in draw_map (its blue water fought the blue
accretion colours). CARTO's no-label tiles were tried and now need an API
key -- they come back watermarked.
```

```text
Latitudes of the places named on the maps and used to split the regions.
Towns are placed at the transect nearest their latitude; labels go landward.
```

```text
Stretches rather than towns, named on the maps only (italic): the reviewer
asked for them in map B, where the flagged inlet transects are.
```

```text
B runs ~15 km past Oregon Inlet, so not "to Oregon Inlet"; and not "to
Bodie Island", which also names the Nags Head-Duck barrier (09-27)
```

```text
Equal aspect by widening x to fill the panel ("datalim"), so the map
keeps the profile's full height and northing lines up across the two.
"box" shrank a diagonal region (C) vertically and broke that. The
limits only settle on a draw, so draw, then fetch tiles for them.
```

```text
above the line: the overview's box D has its top edge on the line
itself, and a label below it was cut by that edge (09-27)
```

```text
every transect drawn alike; the flags are in the table only (09-27)
thin grey edge so rates near zero (near-white) stay visible (09-27)
s=14: at s=3-4 the points hid under the median line (Hannah, 09-27)
```

```text
small and in the scale's end colour: at Oregon Inlet ~40 of them
stack, and full-size black triangles merged into one heavy bar
```

```text
Size the map panel from the region's own shape. With a fixed panel the
equal-aspect map of a diagonal stretch (C) could only fit by trimming
x, which clipped transects at the edge.
```

<details><summary>Function notes (the original docstrings)</summary>

**`greyscale_basemap()`**

```text
Recolour the basemap tiles to light grey: luminance stretched to
[lo, hi], so water sits a shade below land and nothing competes with the
rate colours.
```

**`latitude_ticks()`**

```text
Latitude on the map's left edge (the reviewer: "neither panel shows
coordinates"). The axis is UTM northing, so each tick is the northing of
that latitude at the panel's median coastline longitude: exact at the
coast, which is where the transects are.
```

**`north_arrow()`**

```text
A cartographic north arrow: a solid black notched arrowhead with N
above the tip, on a white box with a hairline border, in the map's
top-right corner. Sized and inset in inches from that corner, so it is
identical in every map whatever the panel's shape. Local to these maps
(Hannah, 2026-09-27: "more academic", "all black and in the upper
right", then "a little larger" with "a white box behind it"); the
house-style _north_arrow used elsewhere is unchanged.
```

</details>

### coastsat/total_change/coastsat_total_change.py

The CoastSat LRR turned into a distance, beside the distance the shoreline actually moved.

From the script's original header:

```text
The CoastSat linear regression rate turned into a DISTANCE, beside the
distance the shoreline actually moved. Built 2026-09-19 (Hannah, by
interview, for her advisor's "total change in shoreline position from the
long-term rate"); split into two named products 2026-09-21 (Hannah, by
interview) after the folder called `lrr_projected/` turned out to hold no
projections at all.

THE VOCABULARY.  A rate turned into a distance is named by the window it was
FITTED on, never by the arithmetic:

  TOTAL SHORELINE CHANGE   --product total  ->  3-rates/coastsat/total_change/
      The rate is evaluated over the SAME window it was fitted on.
      LRR(1996-2010) x 14 yr, LRR(2010-2024) x 14 yr, LRR(1996-2024) x 28 yr.
      Nothing is extrapolated, so nothing is projected. This is what the
      whole of the old `lrr_projected/` tree actually was.

  PROJECTED SHORELINE      --product projected  ->  3-rates/coastsat/projected/
  CHANGE
      The 1996-2024 rate carried onto a window it was NOT fitted on:
      LRR(1996-2024) x 14 yr over 1996-2010, and the same over 2010-2024.
      1996_2024 is deliberately absent -- there it would BE the total change.
      This is the pairing the model's CoastSat target uses in both halves
      (output/comparisons/target_comparison/projected/), here on the
      observations alone.

  OBSERVED CHANGE          both products, unchanged
      No rate anywhere: per transect, the mean position over the whole END
      calendar year minus the mean over the whole START calendar year (all of
      2010 minus all of 1996). Both means are centred mid-year, so the span
      is the same as the rate's multiply, and both use only data inside the
      window. Not the dune-date endpoint in 3-rates/coastsat/endpoint
      (1997-10 to 2023-07, 25.7 yr), which is a shorter span.

    observed - (total or projected) is how far the actual change departs from
    the trend: positive where the shoreline ended up more seaward than the
    trend predicts, negative where more landward. Under --product projected
    it is the more interesting residual of the two, because the rate there
    was never fitted to the window it is being judged over. SEAWARD IS
    POSITIVE throughout, as in every 3-rates product.

SMOOTHED    the same comparison after an alongshore LOWESS (Hannah, by
            interview, 2026-09-21). The rate the MODEL is graded against is not
            the raw rate: it is raw over GIS 1-10 and a 10-domain LOWESS of the
            transect values beyond (cascade_pipeline.coastsat_lowess). The raw
            comparison above therefore tests the fairness of a quantity nobody
            uses; this one tests the target as it is actually applied.

            Note LOWESS commutes with the x years multiply -- the weights depend
            only on the transect positions and the robust reweighting is scale
            equivariant -- so smoothing the RATE and smoothing the DISTANCE
            give the same number to machine precision. Nothing here turns on
            the order; what matters is that BOTH sides are smoothed, at the
            same window, so the residual is not a smoothed quantity minus an
            unsmoothed one.

            Windows 3, 5 and 10 domains (1.5, 2.5, 5.0 km) are all built. The
            sweep is the point: if the residual collapses as the window
            widens, the departures from trend are transect-scale estimation
            noise; if a departure survives 10 domains, the rate genuinely
            fails there, at the scale the model resolves. Read the bias and
            the RMS residual, NOT r -- smoothing strips high-frequency
            variance that is uncorrelated between the two sides, so r rises
            whether or not the smoothing is telling the truth.

FIGURE TITLES carry quantity, window and method, so a figure pulled out of
its folder still says which of the two products it is (Hannah, 2026-09-21):
"Total shoreline change, 1996-2010 (CoastSat LRR 1996-2010 x 14 yr)" against
"Projected shoreline change, 1996-2010 (CoastSat LRR 1996-2024 x 14 yr)" --
the method names the window the rate was FITTED on, so the reader never has
to trust the folder.

OUTPUT   data/hatteras_init/5-scr/3-rates/coastsat/<product>/<start>_<end>/
    transect_<product>.csv             per transect: lrr_m_yr, its uncertainty,
                                       <product>_change_m (and its uncertainty),
                                       the two calendar-year means with their
                                       counts, observed_change_m,
                                       observed_minus_<total|projected>_m
    domain_<product>_summary.csv       per domain: the means of those, std,
                                       pct_landward of each
    <product>_<start>_<end>.png        the rate's distance as the house-style
                                       fill and dots, observed as a black
                                       line; PDF and caption under supporting/
    PROVENANCE.md
    smoothed/<product>_smoothed_<start>_<end>_w<NN>.png
                                       one per LOWESS window, same axis as each
                                       other so the windows can be read side by
                                       side; PDFs and captions under
                                       smoothed/supporting/
    smoothed/tables/domain_smoothed.csv        long: one row per domain per
                                       window (0 = the raw product above)
    smoothed/tables/residual_by_scale.csv      the sweep: bias, RMS, range,
                                       sign agreement and r per window
    smoothed/PROVENANCE.md

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/total_change/coastsat_total_change.py
    python ... --product projected
    python ... --product both
    python ... --windows 1996_2024 2010_2024
    python ... --smooth-windows 3 5 10 | --no-smoothed
```

Notes that were in the code:

```text
The model target's own smoother, imported rather than re-implemented so the
7-domain figure here IS the treatment the runs are graded under (10 until 2026-09-28).
```

```text
The canonical 1996 -> 2010 -> 2024 chain: the full period and its two halves.
Which of these a product is defined for is on the Product below, because
`projected` has no 1996_2024 (see PROJECTED).
```

```text
The observed CoastSat line, PURPLE since 2026-09-22 (Hannah). It was
black, and a black line in 4-comparisons/shoreline_vs_duneline is the
DUNE LINE -- same glyph, two meanings across trees, which is exactly
how this one got read as the dune line. The house ACCENT purple, so it
is a colour the project already uses rather than a new one.
```

```text
ONE fixed metre axis on every figure here (Hannah, 2026-09-22), shared
with 4-comparisons/shoreline_vs_duneline/total_change and
output/comparisons/target_comparison, so a figure from any of the three
can be laid beside another without rescaling by eye. It is FIXED, not a
floor: a window whose data exceeds it is marked at the edge and named in
the caption (hat_figure_style.mark_offaxis) rather than given its own
axis, which would defeat the point. Only Cape Point does, in practice.
```

```text
The 1996-2024 rate carried onto a window it was not fitted on. 1996_2024 is
absent on purpose: there the rate window IS the change window, so the answer
is TOTAL, and building it here would put the same numbers under two names.
```

```text
LOWESS window widths in domain units (1 domain = 500 m). 7 is the model
target's window since 2026-09-28 (coastsat_lowess.LowessConfig.window_domains;
the group's range); 10 was until then and is kept for what still reads it;
3 and 5 are there to show how fast the residual collapses with scale.
```

```text
GIS 1..SPLICE_DOMAINS keep their raw domain means instead of the LOWESS --
coastsat_lowess.LowessConfig.skip_southern_domains, the boundary treatment at
Oregon Inlet. Applied to the OBSERVED side too, so the two never differ in
treatment at any domain.
```

```text
Internal names are neutral (`rate`); Product.cols() renames them to the
product's own on the way out, so no CSV is ambiguous about which it is.
```

```text
Quantity, window, method (Hannah, 2026-09-21). The method is what tells
total from projected at a glance -- "LRR × 14 yr" against "1996–2024 LRR
× 14 yr" -- so it is on the canvas, not left to the folder. draw_fills
puts its bars at 1.025 in axes fractions with the year above them, so
the title has to clear those when the window contains a fill.
```

```text
Both series here are CoastSat. 4-comparisons is where a CoastSat
series meets a dune-line one; this tree never mixes sources, and
after the observed line was read as the dune line it says so; it is
purple now for the same reason.
```

```text
Named for what it is: a symmetric smoother strips variance that is
uncorrelated between the two sides, so this climbs with the window
whether or not the smoothing is right. It is not a score.
```

```text
mean all-NaN draws the frame, grid, village bands and structures with no
sign fill: four curves share the panel, so the blue/red pair is not
available and the ordered ramp below carries the width instead.
```

```text
Quantity, window, method, like every other figure in the tree since
2026-09-21. draw_fills puts its bars at 1.025 in axes fractions and the
year above them, so the title has to clear that when the window contains
a fill -- at the default pad it lands on the 2022 labels. The window's
ROLE in the 1996-2010-2024 chain used to be the title; it is in the
caption now, because the product is the thing a reader cannot recover.
```

```text
At most three columns: five entries on one row ran off both edges of the
canvas once the 7-domain curve joined the overlay.
```

```text
The rate figure exists only where coastsat_lrr_smoothing_windows.py has been run;
a cross-reference to a file that is not there is worse than none.
```

```text
Quantity, window, method -- plus the smoothing width, which is the
only thing that separates these panels from each other.
```

```text
Both series here are CoastSat. 4-comparisons is where a CoastSat
series meets a dune-line one; this tree never mixes sources, and
after the observed line was read as the dune line it says so; it is
purple now for the same reason.
```

```text
A window a product is not defined for is skipped loudly rather than
built: `projected/1996_2024` would be `total_change/1996_2024` under
another name, and two folders of identical numbers is the exact
confusion this rename was done to end.
```

```text
The overlays share one y bound across every window built for this
product, so they are drawn in a second pass once all series exist.
```

<details><summary>Function notes (the original docstrings)</summary>

**`Product()`**

```text
One of the two named products. The ONLY thing that differs between them
is which window the rate is read from; everything downstream -- the
observed side, the figures, the LOWESS sweep -- is identical, which is the
point of building both from one script.

Attributes:
    key: folder name under 3-rates/coastsat/ and the file stem.
    noun: how the quantity is named in a title, caption or legend.
    tok: the token that replaces `rate` in the written column names, so
        every CSV says which product it is without its path.
    windows: the change windows this product is defined for.
```

**`build()`**

```text
The product's rate turned into a distance over `start`-`end`, beside the
observed change over the same years.

The rate is read from `prod.rate_window(start, end)`, which is the window
it was FITTED on -- the same window for TOTAL, always 1996-2024 for
PROJECTED. That single line is the whole difference between the two
products; everything below is shared.
```

**`_smooth_series()`**

```text
One alongshore LOWESS pass at transect resolution, averaged to domains,
with GIS 1..SPLICE_DOMAINS put back to their raw domain means -- the
scoring target's own two steps, shared with the smoothing-scale sweep in
analyze_output/compare_runs/smoothing_scale.py.

Args:
    dom_ids, along_m, values: per-transect domain id, along-coast distance
        in metres, and the quantity to smooth (projected or observed).
    window: LOWESS window width in domain units.

Returns:
    (Series indexed 1..N_DOMAINS, the lowess frac used).
```

**`smooth()`**

```text
The rate-vs-observed comparison repeated under the model target's
alongshore LOWESS, at each window in `windows`. BOTH sides get the same
pass and the same splice, so no window compares a smoothed quantity with
an unsmoothed one. Window 0 in the output is the raw product.
```

**`_smooth_bounds()`**

```text
The fixed metre axis, as everywhere else here (Hannah, 2026-09-22).
Kept as a function so the call sites read the same as before.
```

**`overlay_bounds()`**

```text
One y half-range and tick for EVERY window's overlay, from the
rate-derived series alone.

The per-window panels bound on projected AND observed together; the
overlay draws no observed side, so bounding it that way would set the axis
from a series that is not on the figure. These three are meant to be read
against each other -- 28 yr against two 14 yr halves -- so they share one
bound, and it is not the panels' bound. Said in each caption.
```

**`smooth_overlay_figure()`**

```text
Every LOWESS width's distance on ONE panel, no observed side (Hannah,
2026-09-21).

The per-window panels above each answer "does the trend hold HERE"; this
one answers "what does the window do to the target", which needs the
curves on top of each other and nothing else competing for the eye.

Note this is the RATE figure of input_prep/5-scr/coastsat_lrr_smoothing_windows.py
in metres: LOWESS commutes with the x years multiply, so the curves have
the same shape and only the units differ. It is drawn because metres is
the unit the model and the dune line are read in, not because it shows a
different field.
```

**`method()`**

```text
The parenthetical in a figure title: WHERE THE RATE CAME FROM and
how the distance was made -- "CoastSat LRR 1996-2010 x 14 yr".

The fit window is in the string on purpose (Hannah, 2026-09-21): the
reader compares it against the window in the title, and the two being
equal or not IS the difference between total change and a projection.
So the two products read in parallel and differ in one number:
    Total shoreline change, 1996-2010 (CoastSat LRR 1996-2010 x 14 yr)
    Projected shoreline change, 1996-2010 (CoastSat LRR 1996-2024 x 14 yr)
```

</details>

### coastsat/window_convergence/coastsat_window_convergence.py

Which windows recover the long-term rate, and which are too short?

From the script's original header:

```text
WHICH WINDOWS RECOVER THE LONG-TERM RATE, AND WHICH ARE TOO SHORT?

The model is graded on the LRR over 1996-2010. Is that enough record for the
rate to represent the shoreline, or is it still an artefact of where the record
happens to stop? And if it is not, what window would be?

A NOTE ON COUNTING YEARS. A window here is counted in CALENDAR YEARS OF
RECORD, both ends included, because that is what the fit consumes: 1996-2010
is fifteen years of satellite positions. Rule 2 of ORGANIZATION.md counts the
same window as fourteen SIMULATED years, because the model spends 1996..2009
and the end year is a boundary. Both are right about different things; the
tables and figures here are always the fit's count.

TWO SWEEPS, ONE FROM EACH END (Hannah, 2026-09-23). Both fit the same OLS the
target uses (`coastsat_lrr.compute_lrr`) on a family of NESTED windows, and
both converge on the same reference, the 1996-2024 rate -- but they approach
it from opposite sides, and each answers a different question:

    forward_from_1996/   the START is pinned at 1996 and the END walks out:
                         1996-2000, 1996-2001 ... 1996-2024.
                         HOW MUCH RECORD DO YOU NEED from the start of the
                         chain before the rate settles?

    backward_from_2024/  the END is pinned at 2024 and the START walks back:
                         2020-2024, 2019-2024 ... 1996-2024.
                         HOW LATE CAN A WINDOW BEGIN and still recover the
                         long-term rate?

The pair brackets the answer. Forward gives a window 1996-YYYY; backward gives
a window YYYY-2024. In both the reference is the LONGEST window of that sweep,
which is 1996-2024 either way, fitted in the same loop as every other window
so it cannot drift from a stored product.

The moving year 2010 is marked in both, and in both it is a real model window:
forward that is 1996-2010, the graded window, and backward it is 2010-2024,
the second leg of the canonical chain.

TWO SCALES, in two folders under each direction:

    a-eight_sites/  one transect at the middle of each of eight evenly spaced
                    domains. The readable case: every position, every fit, one
                    panel per site. Eight transects cannot speak for an
                    island, but they show WHAT is happening.
    b-every_transect/  every CoastSat transect on the island, ~906 of them, the
                    same sweep. Answers whether the eight were representative:
                    the convergence year as an alongshore profile, and the
                    spread within each domain.

WHAT THIS CAN AND CANNOT SAY
    The windows are NESTED, so a curve converges on the reference BECAUSE the
    reference is its endpoint. "Does 1996-2010 match 1996-2024" therefore has
    no yes/no answer, and it does not need one -- what the sweep gives is the
    shape of the approach and the window at which the rate stops leaving a
    tolerance. Read it as "this window recovers the long-term rate", never as
    a match test.

    It also cannot separate a rate that was WRONG from one that was merely
    EARLY. A site whose shoreline genuinely changed behaviour in 2010 and a
    site that simply needed more observations draw the same curve. Disjoint
    windows are what separate those -- 3-rates/coastsat/5yr_bins/ is the
    product that asks when the change happened. Running BOTH directions is
    the cheapest partial answer: a real change of behaviour shows as a
    forward sweep and a backward sweep that disagree about where the good
    window is.

THREE TOLERANCES, reported side by side rather than chosen, so the sensitivity
to the choice is visible in the table:

    ci    within the 1996-2024 fit's own 95% confidence half-width. Scales
          with how well constrained each site is; nothing arbitrary to defend.
    abs   within +/-0.25 m/yr. The same band everywhere, so sites compare;
          ~7 m of shoreline over 28 yr, small against rates of 2-3 m/yr.
    rel   within +/-20% of the reference. Sensible where the rate is fast,
          unreachable where it is near zero -- which the table will show.

And two readings of each, because they differ exactly where the answer is
interesting:

    first entry    the shortest window inside the band. May be a crossing the
                   curve then leaves again.
    stable entry   the shortest window after which every LONGER window is
                   also inside. This is the convergence window; the gap
                   between the two is a site that wandered back out.

Inputs
    transect_domain_lookup.csv      2-transect-frame/, via hat_observed_rates
    CoastSat time-series CSVs       1-observations/coastsat_timeseries/
    coastsat_lrr.compute_lrr        5-scr/lib/, via scr_paths

Outputs  (hat_observed_rates.window_convergence_dir(direction, anchor))
    2-settling_window/<direction>_from_<anchor>/
        a-eight_sites/
            shoreline_position_window_fits_*.png   the record: positions, the
                                                   annual median, six fits
            window_convergence_*.png               the sweep as an ERROR
            window_convergence_transects.csv       a row per site per window
            convergence_summary.csv                a row per site
        b-every_transect/
            convergence_alongshore_*.png           the island profile
            domain_convergence_summary.csv         a row per GIS domain
            convergence_summary_all_transects.csv  a row per transect
            window_convergence_transects_all.csv   the full sweep
        c-domain_means/
            the transect sweep grouped to the 90 domains (no refit)
    README.md beside each, supporting/ for PDFs and CAPTIONS.md
    A record other than 1996-2024 files under experiments/record_cut_<end>/.

Usage
    python .../coastsat_window_convergence.py                 both, both scales
    python .../coastsat_window_convergence.py --direction forward --scale sites
    (--abs-tol, --rel-tol, --domains to vary it)
```

Notes that were in the code:

```text
THE RECORD THE SWEEP MAY SEE, and the two pins. The longest window of either
sweep is this one, so both directions converge on the same number, and every
error in the product is measured against it.

TRUNCATING IT IS AN EXPERIMENT, NOT A SETTING. The island-wide detrended
position anomaly steps +17 m between 2020 and 2021 and holds there, in the
same year the satellite record goes from 12.8 to 30.4 observations per
transect per year. Those four years sit at the highest-leverage end of the
1996-2024 fit, so if the step is an artefact of the sensor mix rather than
the shoreline, the reference itself is biased -- and the reference is the
model's grading target. `--ref-end 2020` refits everything on the record
before the step, and the two record spans file side by side under
`experiments/record_cut_<end>/`, apart from the main result in
`2-settling_window/`, so they can be differenced rather than confused.
```

```text
The shortest window either sweep fits. Below five years an OLS through ~90
positions is describing a storm cycle, not a trend.
```

```text
EIGHT EVENLY SPACED DOMAINS (Hannah, 2026-09-23). The 90 domains run south
(1, Cape Point) to north (90, Pea Island); these are every 12th or 13th,
both ends included. Not chosen for behaviour -- an even spread has no
selection argument to defend, and the alongshore gradient is covered.
```

```text
ONE LITERAL TRANSECT PER SITE, not the domain mean (Hannah, 2026-09-23):
the noisiest case, and so an upper bound on how long convergence takes. The
transect is the MIDDLE one of its domain by alongshore order -- a centroid
proxy that needs no geometry, and each domain holds 9-13 transects.
```

```text
SEVEN TOLERANCES, AND NONE OF THEM IS THE ANSWER (Hannah, 2026-09-23, on
being shown how strict the first three were: "would it be possible to make it
more forgiving so it is not an exact match but is close?").

The first three ask a window to land inside a band set by the REFERENCE, and
the reference is the most precise fit in the record -- a median 95% half-width
of 0.151 m/yr over 906 transects. A 15-year window is itself uncertain to
about 0.44 m/yr, so those three judge an imprecise estimate as though it were
as precise as the thing it is being compared to. That asymmetry, not the
shoreline, is why almost nothing passed them.

`overlap` is the principled loosening and the one to reach for first: a window
passes when its OWN 95% interval reaches the reference's, i.e. when the two
fits are not statistically distinguishable. It is a FUNNEL, not a band -- wide
where the record is short and narrowing onto the reference -- which is exactly
the shape the question deserves and the reason it cannot be expressed as a
number. The rest are arbitrary widths, stated openly as such so the answer's
sensitivity to the choice is on the page rather than in a decision no one
wrote down.

ONE THAT WAS TRIED AND REJECTED: the nested-sample standard error,
sqrt(u_window^2 - u_ref^2), the Hausman variance of the difference between an
efficient estimator and a less efficient one. It is correct for nested windows
and more forgiving where the record is short -- but its band COLLAPSES TO ZERO
as the window approaches the full record, and "stable entry" is decided at the
long end, so it comes out STRICTER than everything here (1996-2024 in both
directions). Kept in this comment because a plausible-sounding statistic that
fails for a structural reason is worth one paragraph.

Each entry: (tag, label for figures and prose, test(diff, unc_window,
unc_ref, ref_lrr) -> bool array).
```

```text
THE ONE THE PRODUCT ANSWERS WITH (Hannah, 2026-09-23). Every other tolerance
is still scored and tabled; this is the one the figures draw, the panel
titles name and the READMEs lead with.

CI overlap, because it is the only criterion here that is not a number
somebody chose. It asks whether a window's rate is STATISTICALLY
DISTINGUISHABLE from the long-term rate -- the two 95% intervals meet -- so a
five-year fit is judged against what a five-year fit can actually resolve
(+/-1.8 m/yr) and a twenty-five-year fit against what that can (+/-0.2). The
strict bands ask both to land inside +/-0.15, which is the reference's own
precision and no one else's.

It is worth knowing it barely changes the answer: the island median moves
from 27 to 26 years forward and 25 to 22 backward. That is the finding, not a
disappointment -- the windows disagree because the shoreline changed, not
because the fits were noisy, so no defensible loosening rescues a 15-year
window. The strict bands stay in the tables so that claim can be checked.
```

```text
THE MARKED YEAR, and it is the same number in both directions because in
both it names a real model window: forward it is the end of 1996-2010, the
window the model is graded on; backward it is the start of 2010-2024, the
second leg of the canonical chain. See [[cascade-canonical-periods]].
```

```text
THE FIGURES DO NOT SINGLE OUT THE MODEL WINDOW for now (Hannah, 2026-09-29):
she is choosing which window to use, so nothing is drawn in amber -- no
15-year line on the years-needed figure, no amber fit or dashed continuation
on the eight-site figure; the marked window is drawn like any other. The
tables and README text still score it. Set True to bring the highlight back.
```

```text
SIX WINDOWS ARE DRAWN, NOT TWENTY-FIVE (Hannah, 2026-09-23: "too many lines
to distinguish anything"). The sweep still FITS every window -- the tables
carry all of them -- but a panel that draws them all is a smear in which the
two that matter, the marked window and the reference, are lost among
twenty-three that differ from their neighbour by one year of record.
```

```text
Dropped to naive UTC purely so searchsorted has a datetime64 array to
bisect: a tz-aware column comes back as objects and compares against
nothing. The column itself, which compute_lrr reads, stays tz-aware.
```

```text
The intercept the slope implies, so the straight line can be
DRAWN without a second fit disagreeing with the stored one. x is
years since the window's FIRST observation, as compute_lrr sets it.
```

```text
The UNIT the sweep is about. A transect names itself; a domain
mean names its domain. Everything downstream -- scoring, the
convergence walk, the figures -- keys on this, so the same code
serves a transect and a 10-transect mean without a branch.
```

```text
The reference is always drawn. Its moving year is REF_END going forward
but REF_START going back, which DRAWN_YEARS does not hold, so until
2026-09-29 the backward figure listed the reference and never drew it.
```

```text
THREE PLAIN THRESHOLDS, NOT SEVEN TOLERANCES (Hannah, 2026-09-29, on the old
figures: too many of them, an abstract "rate minus reference" axis, and
tolerance jargon). The figure asks one thing -- how many years of record does
each transect need before its rate stays within X m/yr of the 1996-2024 rate
-- at three values of X, in m/yr and nothing else. All seven tolerances are
still scored in the tables. Strictest last, so it stacks on top.
```

```text
x IS THE TRANSECT COUNT, 1-906 south to north (Hannah, 2026-09-29): every
transect one unit wide, nothing stretched to fit a domain. The domains
are only labels, on a second axis under panel (b).
```

```text
GIS domain labels under the transect axis: piecewise-linear between the
domains' centre transects, so each label sits over its own transects.
```

```text
The transect id only adds something when the unit IS a transect: at the
domain scale it is "GIS 36 mean" beside "GIS 36", which reads as a typo.
```

<details><summary>Function notes (the original docstrings)</summary>

**`windows_for()`**

```text
The nested family, SHORTEST FIRST, so the reference is always last.

Sorting by length rather than by year is what lets one convergence walk
serve both directions: "stable entry" always means "the shortest window
after which every longer one is also inside the band".

Returns [(start, end, moving_year), ...]. The moving year is the end that
is not pinned -- the end year going forward, the start year going back.
```

**`window_label()`**

```text
`1996-2010` going forward, `2010-2024` going back.

The year is rounded because callers pass medians as well as single
windows, and a median over 906 transects arrives as 2022.0 -- which then
printed "1996-2022.0" into a README (found 2026-09-23).
```

**`sweep_one()`**

```text
Every window of one family for one transect, as a list of row dicts.

The series is loaded and sorted once and each window is a SLICE of it, so
a 906-transect run is a couple of minutes rather than an afternoon.
`compute_lrr` still does every fit, so the estimator is the target's.

The reference row (the longest window) is fitted in the same loop as the
rest, so it cannot disagree with them.
```

**`score()`**

```text
Add the reference, the difference and the three in-band flags in place.

The reference is the LAST row -- the longest window of the family, which
is 1996-2024 in either direction.
```

**`convergence_entry()`**

```text
(first entry, stable entry) for one transect under one flag.

`sub` is ordered SHORTEST WINDOW FIRST. First entry is the shortest window
inside the band. Stable entry is the shortest window after which every
longer one is also inside -- the convergence window. Both are reported as
the MOVING year, so forward gives an end year and backward a start year.
The longest window is always inside (it IS the reference), so both exist.
```

**`summarise_domains()`**

```text
One row per GIS domain: the median and spread across its transects.

The median, not the mean: a single transect that never settles would drag
a mean to the end of the record and say the whole domain did.
```

**`domain_mean_sweep()`**

```text
The transect sweep aggregated to the model's own unit: the domain MEAN.

WHY A MEAN OF SLOPES, NOT A FIT THROUGH POOLED POSITIONS. This has to be
the same construction as the grading target, and `coastsat_domain_lrr.py`
builds that by fitting each transect and averaging the slopes. Pooling the
positions first would be a different number, because each transect's
chainage sits on its own arbitrary origin.

No refitting happens here. Every window of every transect is already in
`sweep`; this is a groupby.

TWO UNCERTAINTIES, and the honest one is the wider. `unc_m_yr` is the MEAN
of the transects' 95% half-widths -- the typical fit uncertainty in the
domain. `unc_if_independent_m_yr` is what the half-width of the mean would
be if the ~10 transects were independent samples, which they are not:
they are 10-metre-spaced views of the same shoreline and move together.
Propagating as if they were would divide the band by about sqrt(10) and
manufacture a much later convergence window out of an assumption. The
scoring uses the mean, which is conservative -- it makes convergence look
EARLIER, not later -- and the independent figure is carried as a column so
the size of that choice is visible rather than buried.
```

**`_annual_median()`**

```text
The year-by-year median position: the trajectory under the scatter.

A CoastSat transect gets ~18 positions a year and their spread is tens of
metres, so the raw cloud hides the very shape the fits are arguing about.
The median is per CALENDAR year, matching the window convention.
```

**`draw_fits()`**

```text
The record itself: position against time, with six windows' fits.

The companion to the error figure, and the thing it is derived FROM. y is
shoreline position in metres, x is time across the whole record; each
straight line is one nested window's OLS drawn over the span it was fitted
on, so the spread between them IS the disagreement the sweep measures.

Going forward the marked window is also continued as a DASHED line to the
end of the record: not a claim about the future, but the plainest way to
see what the marked window's record would have had you believe about 2024.
```

**`alongshore_x()`**

```text
Each transect's x on the domain axis: its domain, with the domain's
transects spread evenly across [d - 0.5, d + 0.5] in alongshore order, so
906 transects and the village bands share one axis.
```

**`draw_years_needed()`**

```text
THE figure of the settling sweep: years needed, per transect, alongshore.

One panel per direction, read from the stored every-transect summaries, so
it can be redrawn without refitting. The three thresholds are NESTED -- a
rate that stays within 0.25 m/yr also stays within 0.5 -- so they stack as
shaded bands rather than crossing as lines. The model's 15-year window is
the amber line: wherever a band reaches above it, 15 years was not enough
at that threshold.
```

**`one_direction()`**

```text
Run and write one direction at whatever scales were asked for.

The transect sweep is run ONCE and serves both the all_transects product
and the domain means, which are a groupby of it rather than a refit.
```

</details>

### coastsat/window_convergence/coastsat_window_profiles.py

At what window does the alongshore rate profile start to look like 1996-2024?

From the script's original header:

```text
AT WHAT WINDOW DOES THE ALONGSHORE RATE PROFILE START TO LOOK LIKE 1996-2024?

A companion to `coastsat_window_convergence.py`, drawn the other way round.
That script asks, point by point, when each transect's rate settles into a
tolerance of the reference. This one draws the WHOLE PROFILE -- rate against
alongshore position -- once per window, over the 1996-2024 reference, and
scores each window by one number: the alongshore Pearson r against the
reference, i.e. whether the pattern of erosion and accretion hotspots is right
regardless of offset (Hannah, 2026-09-29, by interview).

THE SAME NESTED FAMILIES, FROM TWO YEARS (Hannah, 2026-09-29):

    1-rate_profiles/forward_from_1996/    1996-1997, 1996-1998 ... 1996-2024
    1-rate_profiles/backward_from_2024/   2023-2024, 2022-2024 ... 1996-2024

The convergence sweep stops at five years because below that an OLS is
describing a storm cycle, not a trend. Hannah asked to see those windows
anyway, on a full y-axis: seeing how far off the short windows are is part of
the question. The minimum is therefore a SEPARATE constant here, and this
script never writes into the convergence sweep's folders or tables, so its
five-year scoring is untouched.

THE UNIT is the transect, all 906 on the island, unsmoothed. The x-axis is the
GIS domain with each domain's transects spread evenly across it in alongshore
(sorted id) order, so the village bands and the domain numbers read as in
every other alongshore figure.

WHAT r CAN AND CANNOT SAY
    The windows are NESTED: each contains the one before and the reference
    contains them all, so r goes to 1 at the reference BY CONSTRUCTION. Read
    where the curve gets there and how steadily, not whether it does. r is
    blind to offset and scale -- a window that has every hotspot in the right
    place at twice the rate scores 1.0 -- and the convergence sweep's bias and
    tolerance tables are the magnitude side of the same question.

Each window is fitted by `sweep_one` from the convergence script, which calls
the target's own `coastsat_lrr.compute_lrr`, so the estimator is the target's.
A window with fewer than MIN_OBS positions on a transect gets no fit there,
and r for that window is over the transects that have one (n in the table).

Outputs  (obs.window_profiles_dir(direction, anchor))
    window_profiles_<direction>_from_<year>.png         (a) every window over
                                                        the reference, (b) r
    window_profiles_panels_<direction>_from_<year>.png  one panel per window
    window_profiles_transects.csv                       a row per transect per window
    window_profiles_correlation.csv                     a row per window: r, n
    README.md, supporting/ (PDFs, CAPTIONS.md)

Usage
    python .../coastsat_window_profiles.py                     both directions
    python .../coastsat_window_profiles.py --direction forward
```

Notes that were in the code:

```text
The sibling sweep, for its loader and its per-transect fit, so the two
products cannot fit a window differently.
```

```text
TWO YEARS, not the sweep's five (Hannah, 2026-09-29): the short windows are
drawn so their distance from the reference can be seen. Two is the least an
OLS through a year boundary can mean anything at all.
```

```text
Windows light (short) to dark (long) on the shoreline-blue ramp, coloured by
LENGTH so both directions read the same way. The reference is the purple
ACCENT, off the ramp, so it cannot be mistaken for a long window.
```

```text
THE PANEL FIGURE IS ZOOMED (Hannah, 2026-09-29: "much more zoomed in"). The
reference spans about -3..+3 m/yr, but the full axis ran to +80 for the
short windows at Cape Point and flattened every panel. The panels share a
fixed +/-PANEL_Y_HALF; a short window that runs past it is cut at the edge,
and the caption says so. The overlay figure keeps the full axis.
```

```text
No window is singled out in the panels (Hannah, 2026-09-29): she is choosing
the window, so the 2010 model window is not titled in amber any more.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_profile()`**

```text
x and rate for one window, in alongshore order, NaN where no fit, so
the line breaks instead of bridging a gap.
```

</details>

### duneline/duneline_endpoint.py

Net dune-line change between the two dune lines that bound a window, per transect and per GIS domain.

From the script's original header:

```text
The stored dune-line observation: net change between the two dune lines that
bound a window, per 100 m transect and per GIS domain, in METRES and in
m/yr. Replaced the dune-line LRR product (3-rates/duneline_lrr/, an OLS
through every line inside a window) on 2026-09-18 (Hannah: "these should not
be lrr, they would just be endpoint, we are tracking net change").

WHAT IS MEASURED
    A window <start>_<end> reads one line per period year through
    hat_topo_version.DUNE_LINE_FOR_YEAR (1996 -> the 1997 line, 2010 -> 2009,
    2024 -> 2023). Each line's per-transect station is ORIG_LEN from
    2-brie-offset/raw_offsets/<vintage>_duneline_offset_raw.csv, the first row
    per transect, exactly as the hindcast's end-year target loader reads it:
    distance from a fixed offshore datum, growing LANDWARD. So per transect
        change_m  = start station - end station       (SEAWARD POSITIVE)
        rate_m_yr = change_m / interval_yr
    and per domain the mean over its ~5 transects (every domain has the same
    transects in both lines, so the mean of the transect changes IS the
    change of the domain means).

    change_m needs no dates. rate_m_yr divides by the interval between the
    two SURVEY DATES (coastsat_vs_duneline.KNOWN_SURVEY_DATES); a vintage
    with no known date is centred on 1 July of its year and flagged in the
    date_assumed columns and the PROVENANCE -- today that is the 2023 line.

OUTPUT   data/hatteras_init/5-scr/3-rates/duneline/endpoint/<start>_<end>/
    transect_endpoint.csv         per transect: positions, change_m, rate_m_yr
    domain_endpoint_summary.csv   per domain: n, mean/std/min/max of both
    PROVENANCE.md                 lines, raw files, dates, island summary
    Read through hat_observed_rates.dune_endpoint_csv(start, end, level).
    rate_windows.py draws the dune line from here and nowhere else.

USAGE
    python scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py          # every window
    python scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py --windows 1996_2010
```

### rates_figures.py

One house-style figure per window for every product under 5-scr/3-rates/, written beside its tables.

From the script's original header:

```text
One house-style figure per window for every product under 5-scr/3-rates/,
written beside its tables (Hannah, 2026-09-18: "I wanted all of these to have
figures"). This replaces the autoscaled quick-looks the LRR fit used to draw
(archived in 5-scr/archive/coastsat_lrr_quicklooks_20260918/).

WHAT EACH FIGURE SHOWS
    The per-domain value as the sign-coloured line and fill of the
    coastsat_lrr_windows panel (blue seaward, red landward; the drawing is
    imported), the individual transects behind it as small dots coloured by
    their OWN sign (the same blue / red; Hannah, 2026-09-18), and the
    village bands, groin and piers, the offshore shoals as faint hatched boxes,
    and the model-input beach fills inside the window as bars above the panel.

        coastsat/lrr/<w>/lrr_<w>.png                    m/yr, +/-1 std dotted
        coastsat/endpoint/<w>/coastsat_endpoint_<w>.png m
        duneline/endpoint/<w>/duneline_endpoint_<w>.png m
        coastsat/5yr_bins/<w>/lrr_5yr_bins_<w>.png      m/yr, one panel per bin
                                                        (coastsat_5yr_bins_figure.py)

    and per MODEL CHAIN (1984-2004-2024, 1996-2010-2024) the chain's two
    windows stacked, earlier above, on the same axis (2026-09-18):
        <product root>/chains/<stem>_chain_<y0>_<y1>_<y2>.png
    for lrr, coastsat endpoint and duneline endpoint. 5yr_bins has no chain
    figure: it covers only the 1996 chain, and its 1996_2024 figure is it.

Y AXES
    lrr        the bound the window figures use: the largest |domain mean|
               over every window plus 1 m, rounded up (+/-8 m/yr today)
    endpoint   ONE bound for both endpoint products, the smallest multiple of
               10 m holding every domain mean of every window, so a shoreline
               figure reads against its dune-line figure directly
    The transect dots are NOT in the bound; a dot outside it is drawn at the
    edge as an open marker and counted in the caption.

USAGE
    python scripts/input_prep/5-scr/3-rates/rates_figures.py
```

Notes that were in the code:

```text
Each transect dot takes the colour of its own sign, the line's blue / red
(the model_vs_observed per-domain dots use the same pair). Grey until 2026-09-18.
```

```text
The curve the runs are actually graded against, laid over the domain means
it is built from (Hannah, 2026-09-22: the observation figure and the model
figure should not show different observed curves). Same darkest blue as the
graded curve on smoothing_windows_<w>.png, so the two read as one series.
```

```text
<start>_<end> only: lrr/1984_2025_obx/ (the all-OBX hand-off, 09-27)
starts with a year but is not a model window and has no domain table.
```

```text
Sentence case, the quantity named by its method; "domain mean" and the
seaward / landward colours are stated in each caption (Hannah, 2026-09-19).
```

```text
Short on purpose: which smoothing width, and that it is the scoring
curve, are caption matters -- a fourth long entry runs past the edge
(lrr_smoothing_windows.win_label learned this the same way).
```

<details><summary>Function notes (the original docstrings)</summary>

**`_target()`**

```text
The scoring curve for one window's transect table: a LOWESS at
TARGET_WINDOW north of the splice, the raw domain means at and below it.

Built through cascade_pipeline.coastsat_lowess, the module the run's own
target is built through, so this figure cannot drift from what the model
is graded against.
```

**`chain_figures()`**

```text
One figure per chain per product: the chain's two windows stacked,
earlier above, on the product's shared axis (the same bound as its
single-window figures). Written to <product root>/chains/.
```

</details>
