# coastsat_shoreline — the shoreline analysis around the groin field

One script measures every shoreline source around the Buxton groin field on
CoastSat's own transects and fits change rates per transect: before the groin,
through its functional life, and after it deteriorated. Its tables are the
observed fillet the hindcast fits read (`../structure_history.md`).

```
HAT_groin_shoreline_analysis.py   the analysis; settings in its CONFIG block
shoreline_output_coastsat/        the 2026-07-13 run. THIS is what the fit reads:
                                  scripts/hatteras_ms/groin-sweep/HAT_fullperiod_target.py
                                  and HAT_groin_position_figure.py load it by path
shoreline_output_grid100m/        the 2026-07-14 run, and where the script writes now
                                  (OUTPUT_DIR); re-running does not touch the folder above
```

Most outputs are git-ignored (the chainage CSV alone is 121 MB, the GIFs tens of
MB); the metadata text and a few tables are tracked. Its GIS inputs are in `../gis_data/`;
`transects_100m.geojson` there
is also read by the 2-brie-offset stage (`data/hatteras_init/2-brie-offset/transects/README.md`).

The script still calls itself `HAT_groin_shoreline_analysis_v2.py` inside the
metadata it writes; that text is part of its output and was left unchanged.

## Restyled 2026-10-01

Brought in line with `scripts/STYLE.md`: a short header with the author block,
one-line comments, a CONFIG block, section titles in place of `# ===` banners,
and every function's docstring moved here. `scripts/repo_tools/style_equivalence_check.py`
proved the restyle changed no code against the previous commit. One change was
made after that proof, on purpose: the unused setting
`BASE = str(_PATH_REPO / "scripts" / "groin_module")`, which named a folder that
no longer exists, was removed. Nothing read it.

Every other path is found from the repo root by searching upward, and each
input file exists. The paths in the two metadata files (`C:\Users\...\scripts\groin_module\...`)
are a record of where the runs read from in July, not settings.

Not checked by the restyle, and worth knowing: the script prints `→`, `Δ` and
similar characters, which the Windows cp1252 console cannot encode
(`scripts/STYLE.md`, "Print UTF-8 on Windows").

## The script in detail

The reasoning that used to sit in the script's comments and docstrings, kept
word for word. The section names follow the script's own section titles.

### HAT_groin_shoreline_analysis.py

#### The original header

```text

HAT_groin_shoreline_analysis_v2.py
================================
Buxton Groin Field effectiveness analysis using multi-source shoreline
observations, all measured against CoastSat's own transect network.

ARCHITECTURE NOTE: wet-dry lines, NC Coastal Management shorelines, and
CoastSat shoreline positions are all extracted onto CoastSat's own
transect network (COASTSAT_TRANSECT_GEOM) -- CoastSat chainage is read
directly from its per-transect CSVs (load_coastsat_chainage()), and
wet-dry/NC state chainage comes from geometric intersection with the
same transects (extract_chainage_by_intersection()). Alongshore
position is built from a trusted along-coast ID order (see
_compute_alongshore_positions()), and CASCADE domain assignment comes
from the authoritative HAT_domains.json reference (see
load_domain_reference() / assign_domain_from_northing()), not a
formula.

Analyses (each per-transect, aggregated alongshore):

    1. PRE-INSTALLATION LRR (all shorelines dated before
       PRE_INSTALLATION_YEAR_CUTOFF, currently 1970)
       Establishes background shoreline change rate before groin
       installation, using every historical/aerial shoreline available
       pre-1970 rather than a fixed hand-picked list of years, pooled
       across every source that has pre-1970 data (a per-source
       breakdown is printed). The 1966 Buxton nourishment does bias the
       1966/1967 shorelines at Buxton-area transects -- those
       observations are KEPT in the regression (not dropped) and only
       FLAGGED (nourishment_zone=True), so a reader can see exactly
       which points are nourishment-influenced rather than having them
       silently removed. See NOURISHMENT_EXCLUSIONS.

    2. FULL-PERIOD REGRESSION with SIMPLE PIECEWISE CHANGE-POINT DETECTION
       Fits a single-slope regression over EVERY observation a
       transect has (no year cutoff), and separately searches
       BREAKPOINT_SEARCH_START..END for the best two-segment
       (piecewise linear) breakpoint year. This is a data-driven check
       independent of the assumed era boundaries below: does the
       record itself suggest a regime change, and if so, roughly when.
       (The RSS-reduction fraction reported is descriptive, not a
       significance test -- an exhaustive breakpoint search will
       always find *some* "best" breakpoint even under a true
       single-slope null, so a large reduction is suggestive, not
       confirmatory, without a permutation/Davies-type test.)

    3. PER-ERA LRR
       One rate per transect per structural era. Three eras (see
       ERAS below and HAT_groin_zone_investigation.py for the full
       cited maintenance history): Pre-install (no structure),
       Functional groin (1970-1995: built, damaged/repaired during
       construction, 1975 steel-pile repair, 1980/82 anti-flanking
       extensions, 1994 Gordon damage, 1995 south-groin repair -- the
       LAST documented repair), and Deteriorated (1996-2024: no
       further repairs after 1995, including through Hurricane Isabel
       in 2003, which is annotated on the decadal plots as a notable
       event but is not an era boundary since no maintenance bracket
       it). Shown as per-transect scatter only -- no binning, no
       smoothing, no connecting line, for any era.

    4. DECADAL LRR + ANOMALY (LRR minus pre-install baseline)
       At each transect, fits LRR within fixed, NON-OVERLAPPING
       DECADE_LENGTH_YEARS-year bins starting at DECADE_START_YEAR
       (1960s, 1970s, ...) -- replaces an earlier 10-year window that
       slid one year at a time (adjacent "center years" shared ~90% of
       their underlying observations, which read as smoother/more
       continuous than the survey frequency actually supports). Each
       decade here is a distinct, non-overlapping slice of time.
       The regional baseline subtracted to get the anomaly connects the
       sorted per-transect pre-install points directly (linear
       interpolation, no binning/smoothing) -- see
       interp_preinstall_baseline(). The decadal/anomaly/distance-band/
       signal-extent CSVs are all still written even though their
       companion plots were dropped as unhelpful (see Outputs below);
       the numbers may still be useful even without a chart.

    5. SIGNAL EXTENT OVER TIME (computed, not plotted -- see Outputs)
       The CONTIGUOUS alongshore distance, measured outward from the
       groin, over which the anomaly exceeds SIGNAL_ANOMALY_THRESHOLD_M_YR
       -- i.e. "how far" the groin's influence reaches on each side,
       tracked by decade. The threshold is checked against an empirical
       noise floor (far-field pre-install rate standard deviation,
       printed by compute_preinstall_lrr) rather than assumed outright
       -- see the config comment above SIGNAL_ANOMALY_THRESHOLD_M_YR.

    6. SHORELINE EVOLUTION GIF (groin area)
       Animated view of RAW shoreline position through time, zoomed to
       domain 1 to PLOT_UPDRIFT_MAX_KM updrift of the groin, using the
       same multi-source
       chainage table as everything else (so it covers the full
       historical record, not just the CoastSat era). No geometric
       correction of any kind -- CoastSat's transects are shore-normal
       and run parallel to the coast, so chainage is used exactly as
       provided. GIF_FRAME_MODE picks "date" (default -- one frame per
       unique observation date, using every observation) or "year"
       (pools all sources/dates within a calendar year -- fuller
       per-frame coverage, coarser resolution). Adapted from the
       whole-island version in shoreline_chainage_alldata_evolution.py.

Data sources (identical to shoreline inventory):
    Aerial wet-dry lines   - user-digitized (GeoJSON)
    NC Coastal Management  - historical shoreline shapefile (GeoJSON)
    CoastSat               - satellite-derived, one CSV per transect
    (CoastSat imagery doesn't reach back before ~1984, so the
    pre-installation baseline is necessarily wet-dry + NC state only --
    that's expected, not a bug; see the per-source breakdown printed by
    compute_preinstall_lrr.)

Chainage extraction:
    CoastSat  - read directly from its per-transect CSVs
    Wet-dry   - geometric intersection with transect line (shapely)
    NC state  - same intersection approach

All shoreline sources are treated as equal-weight observations in
regression. Wet-dry (~1-2 m precision) is more accurate than CoastSat
(~10 m), so weighted regression is a reasonable extension for later.

Distance from groin is measured ALONG THE SHORE (curvilinear, the same
cumulative alongshore distance used everywhere else), not as a
straight-line northing difference -- see assign_distance_from_groin().
The origin (0 km) is the NORTHERNMOST individual groin feature, since
updrift = north; the whole groin field (first to last groin) is shaded
on the profile plots. CASCADE domain number is the PRIMARY x-axis on
the profile plots and the GIF, with distance-from-groin (km) on a
secondary top axis -- see set_domain_primary_axis().

Curves vs. smoothing: no binning or averaging anywhere -- every series
(pre-install baseline, full-period, post-install, each era, era-to-era
differences) is real per-transect data, connected in distance order so
the trend reads as a line rather than a loose cloud of dots. Binning
was deliberately removed: even a non-smoothing fixed-width bin pools
multiple transects together and can shift or mute a real, spatially
localized groin signal.

The per-era profile plot (plot_era_lrr_profile) additionally draws a
LOWESS-smoothed overlay (SMOOTHED_OVERLAY_LOWESS_FRAC bandwidth, a
lighter shade of each line's own color) as a supplementary visual aid
-- it is drawn ALONGSIDE the raw connected points, never in place of
them, so the actual per-transect data is always visible underneath.

The regional pre-install baseline used for the DECADAL ANOMALY
calculation (a different thing from the plots above) connects the
sorted per-transect pre-install points directly via linear
interpolation -- see interp_preinstall_baseline().

Outputs
-------
  groin_analysis_chainage_all.csv                 - unified (transect, date, source, chainage) records
  groin_analysis_preinstall_lrr.csv               - per-transect pre-install LRR (nourishment-flagged, not dropped)
  groin_analysis_full_period_lrr.csv              - per-transect full-period LRR + piecewise breakpoint
  groin_analysis_post_install_lrr.csv             - per-transect post-install-only (1970-present) LRR
  groin_analysis_decadal_lrr.csv                  - long-format decadal (fixed, non-overlapping) LRR [no plot]
  groin_analysis_era_lrrs.csv                     - per-transect per-era LRR
  groin_analysis_decadal_anomaly.csv              - decadal LRR minus interpolated pre-install baseline [no plot]
  groin_analysis_signal_extent.csv                - contiguous signal extent by decade [no plot]
  groin_analysis_lrr_distance_band_series.csv     - decadal LRR aggregated into distance bands x decade [no plot]
  groin_analysis_anomaly_distance_band_series.csv - decadal anomaly aggregated into distance bands x decade [no plot]
  groin_analysis_alongshore_profile.png           - LRR vs. CASCADE domain / distance from groin (pre-install vs. full-period)
  groin_analysis_era_lrr_profile.png              - LRR vs. CASCADE domain / distance from groin, one curve per era
                                                     (+ pre-2022-nourishment subset, + LOWESS-smoothed overlay per line)
  groin_analysis_era_difference_profile.png       - era-to-era Δ LRR, both comparisons together (Functional-Pre-install, Deteriorated-Functional)
  groin_analysis_diff_pre_to_functional.png       - Δ LRR, Functional groin minus Pre-install baseline (focused single-comparison version)
  groin_analysis_diff_functional_to_deteriorated.png - Δ LRR, Deteriorated minus Functional groin (focused single-comparison version)
  groin_analysis_diff_deteriorated_to_preinstall.png - Δ LRR, Deteriorated minus Pre-install baseline (whole historical arc, skipping the middle era)
  groin_analysis_decade_increment_lrrs_{N}yr.csv  - per-transect LRR for each N-year window since installation (one file per DECADE_PLOT_INCREMENTS_YEARS entry)
  groin_analysis_decade_increment_profile_{N}yr.png - LRR profile per post-install window (see above), color = chronological order
  groin_analysis_zone_panels.png                  - one stacked panel per era, downdrift/updrift analysis zones shaded with mean-rate labels
  groin_analysis_shoreline_evolution.gif          - animated shoreline position (chainage) through time, groin area --
                                                     water/island shaded around the curve, y-range fixed to the full record
  groin_analysis_shoreline_evolution_zoomed.gif    - same, but zoomed to domain 1-GIF_ZOOM_WINDOW_DOMAIN_MAX for a closer look around the groin
  groin_analysis_metadata.txt                     - run metadata

(The distance-band, anomaly-profile, signal-extent, and transect-
diagnostic PLOTS were removed entirely -- not just unused -- since they
weren't helpful. Their underlying compute_* functions and CSV exports
are still here; only the plotting functions are gone.)

Usage
-----
Edit CONFIG below, then run:
    python HAT_groin_shoreline_analysis_v2.py
```

#### Settings

Each note sat above the named setting in the CONFIG block.

**`_PATH_REPO`**

```text
Anchored 2026-09-14 (rule 5 of ORGANIZATION.md). The anchoring pass pasted
this INSIDE the docstring above, so it never ran and every _PATH_REPO below
was a NameError; moved out 2026-09-18.
```

**`NC_STATE_PATH`**

```text
Repointed 2026-09-14 at the data tree. gis_data/ held a BYTE-IDENTICAL second
copy of this 20 MB file (md5 330d5f32...), so every clone carried it twice and
an edit to one would have silently diverged from the other. The data tree is
the owner -- shoreline_inventory.py reads the same file from there -- and
this script already crosses into data/hatteras_init for the CoastSat transect
geometry two lines below, so the direction is not new.
```

**`COASTSAT_TRANSECT_GEOM`**

```text
These two named coastsat_lrr/ and had not resolved since the 2026-09-12 move.
```

**`DOMAINS_JSON_PATH`**

```text
CASCADE domain reference (authoritative)
A GeoJSON of the 90 real CASCADE domain boxes (D1-D90), each a 500 m
polygon with a "domain_id" property. This REPLACES the original formula-
based domain assignment (CAPE_POINT_NORTHING + fixed 500 m spacing),
which was measurably wrong -- domain 1's real midpoint northing is
3,899,029 per this file, not the ~3,897,750 the original formula assumed
(over 2 domains off). See load_domain_reference() /
assign_domain_from_northing().
```

**`PLOT_UPDRIFT_MAX_KM`**

```text
Plot x-axis limit
Distance-from-groin plots (alongshore profile, era profile) are
clipped to start exactly at CASCADE domain 1 (the southernmost real
domain -- see domain_dist_km()) on the downdrift side, and to
PLOT_UPDRIFT_MAX_KM on the updrift side.
```

**`GROIN_GEOJSON_PATH`**

```text
Groin location
The Buxton Groin Field is provided as a shapefile/GeoJSON with
multiple LineString features (the individual groins). The script
loads the geometry and uses:
- the NORTHERNMOST groin's northing as the origin for "distance
from groin" (0 km), since updrift = north here
- the full north-south extent (southernmost to northernmost groin)
as a shaded band on plots, showing the whole groin field's footprint
```

**`GROIN_NORTHING_FALLBACK`**

```text
Fallback northing if the geojson can't be loaded. Used only if the
geojson is missing. Set to None to error out instead.
```

**`UPDRIFT_DIRECTION`**

```text
Direction of "updrift" — the side where sediment approaches from.
On Hatteras Island, net alongshore transport is southward, so updrift
is NORTH of the groin. Signed distance from the groin is then
positive = north = updrift, negative = south = downdrift.
```

**`GROIN_INSTALLATION_YEAR`**

```text
Groin installation year (used as fixed breakpoint reference,
and as the cutoff for the pre-installation baseline below)
```

**`PRE_INSTALLATION_YEAR_CUTOFF`**

```text
Pre-installation baseline
ALL shoreline observations dated strictly before GROIN_INSTALLATION_YEAR
are used as the pre-groin baseline -- not a hand-picked list of years,
and pooled across EVERY source that has pre-1970 data (aerial wet-dry
lines and/or NC Coastal Management shorelines; CoastSat doesn't reach
back this far). Nothing is dropped for being pre-1970 -- the code
prints a per-source breakdown so it's visible exactly what went in.
```

**`NOURISHMENT_EXCLUSIONS`**

```text
Nourishment zones by (name, year, y_min, y_max): the 1966 Buxton
nourishment placed sand directly in front of the groin field, so the
1966 and 1967 shorelines there reflect a fill event, not the natural
pre-install trend. Per-committee direction: keep this data in the
regression (don't silently drop real observations) -- these
transects are only FLAGGED (nourishment_zone=True in the comparison
tables, shown as a distinct marker on the profile plot) so a reader
can see exactly which points are nourishment-influenced and judge for
themselves, rather than having them removed upstream.
```

**`ERAS`**

```text
Analysis periods (for era-averaged LRR profiles)
Three eras derived from groin structural/maintenance condition
(CSE 2013 documented history -- see HAT_groin_zone_investigation.py
for full citations):
Pre-install     - 1849 to 1969 (baseline, no structure)
Functional groin - 1970 to 1995 (structure actively maintained:
built 1970, storm/construction damage same
year, 1975 steel-pile repair on groin #1,
1980/82 anti-flanking extensions, 1994 Gordon
damage, 1995 south-groin repair)
Deteriorated    - 1996 to 2024 (LAST documented repair was the
1995 south-groin steel sheet piling -- nothing
after that. Hurricane Isabel (2003) is a real,
notable event WITHIN this era, not the era
boundary itself, since no maintenance happened
either right before or after it; it's annotated
on the decadal plots instead of used as a split)

NOTE: this merges what were previously two separate eras
("Post-Gordon" 1995-2003 and "Deteriorated" 2004-2024) into one,
since the maintenance history doesn't actually support treating them
as functionally distinct -- both are "no more repairs happened."
If the piecewise-breakpoint search (compute_fullperiod_lrr) turns up
a real, data-supported rate change clustering near 2003 across many
transects, that would be evidence FOR re-splitting at Isabel --
worth checking rather than assuming either way.
```

**`NOTABLE_EVENTS`**

```text
Notable storm/repair events -- shown as annotated vertical lines on
the decadal plots even though they no longer define era boundaries.
```

**`DECADE_START_YEAR`**

```text
Decadal window parameters (replaces the original sliding rolling window) ──
Fixed, NON-OVERLAPPING 10-year bins are much easier to read and
report than a continuously-slid 1-year-step window: in the original
scheme, adjacent "center years" shared ~90% of their underlying
observations (a 10-year window stepped by 1 year overlaps 9 of its
10 years with its neighbor), which made the resulting time series
look smoother and more continuous than the underlying survey
frequency actually supports, and made it hard to say what was
genuinely new information between e.g. center-year 1994 and 1995.
Decades start at DECADE_START_YEAR and run every DECADE_LENGTH_YEARS
until the data runs out -- each number on a plot is now a real,
distinct slice of time with no double-counted observations.
```

**`SIGNAL_ANOMALY_THRESHOLD_M_YR`**

```text
Signal-extent detection
The alongshore extent of the groin signal is defined as the largest
CONTIGUOUS run of domain-bins, starting immediately at the groin and
working outward, whose |anomaly| (decadal LRR - pre-install baseline)
exceeds SIGNAL_ANOMALY_THRESHOLD. Requiring contiguity from the groin
(rather than "furthest point anywhere that crosses the threshold")
stops a single noisy far-field transect from being reported as part
of the groin signal. Signal extent is computed separately updrift and
downdrift. See compute_signal_extent_over_time().

Is 1.0 m/yr justified? It's a round, literature-typical number, and
it IS checked against this site's own noise floor: compute_preinstall_lrr()
computes the standard deviation of pre-install LRR at transects
beyond SIGNAL_MAX_SEARCH_DISTANCE_M (far enough that the groin can't
plausibly reach) and prints how many multiples of that std dev this
threshold represents. A defensible threshold is roughly 1-2x that
natural background standard deviation, so the detector isn't flagging
ordinary shoreline variability as "groin signal." READ THAT PRINTED
NUMBER on your first real run -- if the ratio isn't in a reasonable
1-2x range, adjust this constant to match your site's actual noise
floor rather than leaving it at the round default.
```

**`SIGNAL_EXTENT_BIN_WIDTH_M`**

```text
End-of-line comment was: bin width for contiguity check
```

**`SIGNAL_EXTENT_BIN_WIDTH_M`**

```text
End-of-line comment, continued: (matches CASCADE domain length)
```

**`SIGNAL_EXTENT_MAX_GAP_BINS`**

```text
End-of-line comment was: allow this many below-threshold
```

**`SIGNAL_EXTENT_MAX_GAP_BINS`**

```text
End-of-line comment, continued: bins in a row before stopping
(absorbs single-bin noise)
If a transect lacks a pre-install LRR (no historical shorelines on
it), fall back to a regional baseline built from the transects that
do have pre-install LRR -- see interp_preinstall_baseline(). This
connects the sorted per-transect pre-install values directly (linear
interpolation), not a binned step function or a smoothed curve. No
LOWESS/kernel smoother is used anywhere in this script.
```

**`DISTANCE_BAND_EDGES_M`**

```text
Distance bands used for the "signal through time" curve plots that
replace the transect x year heatmaps. Each band gets its own curve
(decadal LRR, or anomaly, vs. decade) applied
symmetrically updrift and downdrift.
```

**`BREAKPOINT_SEARCH_START`**

```text
Piecewise breakpoint search
The piecewise-linear regression per transect uses ALL observations
(1849–present). This range is only the SEARCH space for the single
breakpoint year. Extending it later lets the analysis catch regime
changes near the recent record.
```

**`TRANSECT_INTERSECTION_LENGTH_M`**

```text
Transect chainage extraction
Length of transect line (from origin, seaward) used for shoreline
intersection. CoastSat transects are typically ~500 m; a longer value
handles cases where the actual shoreline sits further seaward than
the CoastSat mean line.
```

**`CHAINAGE_MAX_ABS_M`**

```text
Maximum absolute chainage retained. Anything beyond this is treated
as a spurious intersection (e.g., two shorelines happened to overlap
far from the intended cross-shore location).
```

**`GIF_FRAME_MODE`**

```text
Shoreline evolution GIF (groin area)
Animated view of RAW shoreline position (chainage) through time,
zoomed to the groin area -- no geometric correction of any kind:
CoastSat's transects are shore-normal and run parallel to the coast,
so chainage is used exactly as provided. Adapted from the whole-
island version in shoreline_chainage_alldata_evolution.py, but
built on THIS script's multi-source unified chainage table (wet-dry +
NC state + CoastSat) so it can show the full historical record, not
just the CoastSat era (1984+). See create_groin_evolution_gif().

GIF_FRAME_MODE controls what one "frame" is:
"date" -- one frame per unique observation DATE, using every
observation the data supports (not collapsed into
years). Maximum temporal resolution, but a single
wet-dry or NC-state date, or an early/partial CoastSat
pass, often only covers PART of the window, so
individual frames can look spatially sparse -- that's
real data sparsity on that date, not a bug.
"year"  -- one frame per calendar year, pooling every source/date
within that year (median). Fuller per-frame spatial
coverage, coarser temporal resolution.
```

**`GIF_MIN_TRANSECT_FRACTION`**

```text
A date/year needs to clear BOTH of these bars to get a frame --
lowered the fraction from the original 0.15, which was excluding most
of the historical record: a single wet-dry/NC-state survey, or an
early CoastSat pass, usually covers a small fraction of the window on
its own, so a 15%-of-window bar was throwing out the majority of real
observations. The absolute floor just guarantees a frame has at
least a few real points to draw a meaningful curve from, even when
the fraction bar alone would let through 1-2 lone transects.
```

**`GIF_PRE_COASTSAT_CUTOFF_YEAR`**

```text
Before CoastSat existed, every historical shoreline is precious and
inherently sparse -- so ANY frame with at least GIF_MIN_TRANSECT_ABS
transects gets kept regardless of coverage fraction, for any date
before this year. GIF_MIN_TRANSECT_FRACTION only applies from this
year onward, where CoastSat's dense modern coverage makes a
meaningful fraction bar reasonable to enforce.
```

**`GIF_ZOOM_WINDOW_DOMAIN_MAX`**

```text
Second GIF, zoomed to domain 1 through this domain (instead of the
full domain-1-to-PLOT_UPDRIFT_MAX_KM window) for a closer look at
change right around the groin.
```

**`GIF_REFERENCE_YEAR`**

```text
Last pre-install shoreline in the record -- drawn as a fixed,
always-visible reference line on every frame dated AFTER this year
(not on frames before/during it, since it's meant as a baseline for
what came after, not a comparison for the earlier record). Built the
same way as any other frame -- median chainage per transect for this
specific year, gaps left where no observation exists that year -- so
it's exactly as sparse/honest as the real 1967 coverage actually is.
```

**`GIF_WATER_COLOR`**

```text
Background shading so each frame reads as a real coastal cross-
section rather than an abstract line chart: everything seaward of
the shoreline (larger chainage) shaded like water, everything
landward (smaller chainage, back toward the transect origin) shaded
like the island itself.
```

**`GIF_TRANSECT_SOURCE`**

```text
GIF transect source
Which transects define alongshore position and/or measure shoreline
position, for the GIF SPECIFICALLY -- every other analysis in this
script always uses CoastSat's own shore-normal transects, since those
measure true cross-shore distance regardless of which way they
happen to face, which is what makes them correct for computing rates.

"coastsat" (default) -- CoastSat's own transects for BOTH the
chainage measurement (correct, shore-normal) and the alongshore
x-axis position (via its trusted ID order).
"hybrid" -- keeps CoastSat's own (correctly-measured) chainage
values exactly as they are, but repositions each one along the
x-axis using the nearest transects_100m.geojson transect's own
alongshore position instead of CoastSat's ID-order-based one --
perfectly even, precisely-known 100 m spacing on the x-axis,
without changing how the shoreline itself was measured.
"grid100m" -- fully re-measures the shoreline against the 100 m
grid's OWN transects (parallel, not shore-normal): each
chainage observation gets reconstructed into a real (x, y)
point using its original CoastSat transect, then re-projected
onto the nearest 100 m-grid transect. WARNING: the 100 m grid's
transects all point the same fixed direction rather than
following the coast's curve, which reintroduces a known
geometric artifact (a spurious drift in raw position wherever
the true coast diverges from that fixed direction -- see the
"IMPORTANT" note this caused several turns of debugging
earlier in this project). This option exists so you can see
that for yourself, not because it's recommended.
```

**`NUM_REAL_DOMAINS`**

```text
CASCADE domain assignment
Domains are ~500 m long, numbered from south (D1) to north (D90).
Real boundaries come from DOMAINS_JSON_PATH (see load_domain_reference()
/ assign_domain_from_northing()) -- NOT a formula anymore. The original
formula (D1 midpoint at a fixed "CAPE_POINT_NORTHING") was measurably
wrong: real D1 midpoint is at northing 3,899,029 per the domains file,
over 2 domains away from what the formula assumed.
```

**`PRE_NOURISHMENT_PERIOD`**

```text
Extra display-only sub-period: Deteriorated era before the 2022
nourishments (Buxton 2022, Avon 2022) -- shown alongside the full
Deteriorated line on the era profile plot, NOT part of the official
ERAS list (so it doesn't touch the anomaly/signal-extent pipeline,
which is built on the 3 official eras only).
```

**`DECADE_PLOT_START_YEAR`**

```text
Decade-increment LRR evolution plot
A finer-grained complement to the 3-era profile: shows the LRR
profile across successive fixed-length windows since installation
(e.g. 1970-1980, 1980-1990, ...), to see whether/how the groin's
effect has changed continuously over its lifespan, not just as a
single before/after-deterioration split.
```

**`ZONE_PANEL_DOWNDRIFT_DOMAINS`**

```text
Zone-panel profile plot
Adapted from HAT_groin_zone_investigation.py's alongshore_profile_panels
figure -- one stacked panel per era (using THIS script's own era
definitions/data, not that script's separate rate pipeline), with
shaded analysis zones and a mean-rate label per zone per panel.
```

**`ZONE_PANEL_FUNCTIONAL_EARLY_WINDOWS_YEARS`**

```text
Extra zone-panel rows: the Functional groin era restricted to just
its first N years (starting at GROIN_INSTALLATION_YEAR), to see
whether the groin looked most effective early on rather than across
its whole 1970-1995 span. Empty list to skip.
```

**`DELTA_COLOR`**

```text
Used for any DIFFERENCE/delta line (era-to-era Δ LRR) -- deliberately
NOT any color already used for a raw era curve (grey/green/red/
orange), so a delta is never visually confusable with one of the two
raw curves it was computed from.
```

**`SMOOTHED_OVERLAY_LOWESS_FRAC`**

```text
LOWESS bandwidth (fraction of points) for the smoothed-overlay lines
on most profile plots -- a VISUAL AID only, drawn as a lighter
shade of each line's own color, in addition to (not instead of) the
raw per-transect points + connecting line. Lowered from 0.08 so the
smoothed line tracks the raw per-transect data more closely (less
aggressive smoothing) while still reading as a clean line.
```

**`ERA_PROFILE_LOWESS_FRAC`**

```text
Dedicated (lower) bandwidth for the structural-era profile
specifically -- the sharp rate jumps right around the groin get
smeared out at the shared 0.04 bandwidth, and this plot in
particular is the one where those sharp local jumps are the point.
```

**`ANN_TOWN_SPANS`**

```text
SECTION 6: GEOGRAPHIC ANNOTATION STYLING (annotated publication figure)
Community/feature positions given as CASCADE domain numbers (or
domain ranges), converted to distance-from-groin (km) via
build_domain_to_dist_km() so they land in the right place regardless
of which plot or x-window they're drawn on. See
add_geographic_annotations().
```

**`ANN_TOWN_SPANS`**

```text
Tri-Village removed -- it and Rodanthe sit right at/beyond the
edge of the 30 km window and were rendering intersected with the
plot boundary rather than cleanly inside it.
```

**`ANN_MAX_DOMAIN`**

```text
General cutoff: any annotation (village line, town span, shoal zone)
past this domain is skipped or clipped -- e.g. "Waves" (domain 74)
sits past this and is dropped entirely; a span like Wimble Shoals
that starts before the cutoff but extends past it gets clipped to
end AT the cutoff rather than being dropped outright.
```

**`ANN_WIMBLE_SHOALS`**

```text
NOTE: the Buxton Groin itself is NOT re-added here via ANN_GROINS --
it's already the central, uniquely-marked feature of every plot in
this script (the "Groin" line + "Groin field" shaded band), so a
second label for the same structure would just be a duplicate.
```

**`LABEL_ACCRETION_Y`**

```text
Accretion / Erosion side labels - set to None for auto-computed midpoint,
or a 0-1 axes fraction to pin to a fixed position.
```

**`ANN_C_TOWN_SPAN_EDGE`**

```text
Darker outline colors for the town-span / shoal shaded regions --
their fill colors are close to the GIF's ocean/island background, so
a plain fill alone can be hard to see; a darker border makes the
extent of each region clear regardless of what's underneath it.
```

#### Functions

In file order. The docstring first, then the notes from inside the body, each
with the line it stood above in the previous version of the script.

##### `_require_file()`

```text
Check a configured file path exists before handing it to
geopandas/pyogrio, which otherwise raises a deep, unhelpful GDAL
traceback for a plain "you haven't saved this file here yet" typo.
```

##### `fix_crs()`

```text
Reproject to PROJECTED_CRS, handling missing or invalid CRS.
Logs native CRS + reprojected bounds to _CRS_LOG for the
consolidated spatial-reference confirmation (see print_crs_summary()).
```

##### `print_crs_summary()`

```text
Consolidated confirmation of the spatial reference system for
every shoreline/geometry file loaded this run: native CRS as read
from the file, confirmation it was reprojected to PROJECTED_CRS,
and the reprojected bounding box so a wrong/garbled reprojection
(wildly out-of-range coordinates) is obvious at a glance rather
than silently corrupting downstream distances.
```

##### `parse_date_column()`

```text
Auto-detect date encoding and parse to pandas datetime.
Handles: string dates, Unix milliseconds since epoch, year-only ints.
```

##### `decimal_year()`

```text
Convert datetime to decimal year (e.g., 1985-07-02 → ~1985.50).
```

##### `load_domain_reference()`

```text
Load the authoritative CASCADE domain boundaries (DOMAINS_JSON_PATH)
-- 90 real domain polygons (D1-D90), each ~500 m alongshore, with a
domain_id property. Replaces the original CAPE_POINT_NORTHING-based
formula, which was measurably wrong (real D1 midpoint is at northing
3,899,029 per this file, over 2 domains away from what the formula
assumed).

Returns DataFrame: domain_id, y_min, y_max, y_mid -- sorted by
domain_id.
```

##### `assign_domain_from_northing()`

```text
Assign northing(s) y (scalar or array-like) to a CASCADE domain
number using the authoritative domain midpoints from
load_domain_reference() -- nearest-midpoint assignment, so the
small (~few meter) gaps/overlaps between the source polygons don't
leave any northing unassigned. Values beyond the outermost real
domain (buffer-zone transects) clip to D1 or D90 rather than
extrapolating.
```

##### `domain_dist_km()`

```text
Mean dist_from_groin_m (km) of transects assigned to the given
CASCADE domain -- used to anchor a plot's x-limit to a specific
domain (e.g. "start the x-axis at domain 1") instead of wherever
the data happens to start or stop. Returns None if that domain has
no transects in this dataset.

Prints the transect count and computed position so it's directly
checkable in the run log -- if only one or two transects happen to
be assigned to that domain (e.g. near the Cape Point bend, where
domain assignment can be less certain), the resulting position may
look further out than a naive "domain_number x 500 m" estimate
would suggest; that's a real reflection of sparse/awkward coverage
there, not necessarily a bug.
```

##### `build_domain_to_dist_km()`

```text
Build a callable domain_number -> dist_from_groin_m (km) mapping
(linear interpolation across whichever domains have transects), and
its inverse dist_km -> domain_number. Shared by
set_domain_primary_axis() and add_geographic_annotations() so every
domain-based reference on every plot uses the exact same mapping.

Returns (domain_to_dist, dist_to_domain), both vectorized.
```

##### `set_domain_primary_axis()`

```text
Make CASCADE domain number the PRIMARY (bottom) x-axis label,
moving distance-from-groin (km) to a secondary (top) axis instead.
The underlying plotted x-COORDINATE is unchanged -- still distance-
from-groin in km, which is what all the binning/curve positions use
-- this only changes which unit is the main axis label.

Uses the domain assignment already computed per transect
(transects['domain'], from the authoritative HAT_domains.json --
see assign_domain_from_northing()), so it stays perfectly
consistent with every other domain reference in this analysis.

Call this AFTER the caller has already set the final x-limits
(e.g. via apply_distance_xlim_and_autoscale_y), so only domains
actually visible in that range get a tick.
```

Above `edge_margin_km = 0.6`:

```text
Exclude ticks within a small margin of the plot's edges -- a
label sitting right at the boundary (e.g. D70 at the very edge
of a 30 km window) tends to collide with geographic annotation
text placed near the edge, rather than reading as a normal tick.
```

Above `domain_min, domain_max = agg.index.min(), agg.index.max()`:

```text
D65 is explicitly always shown if it's anywhere in the true
visible range (not shrunk by edge_margin_km) -- the margin
exclusion above is meant for ticks landing right at the very
edge (like D70), not for D65, which sits well within the window
for the standard 30 km-updrift view. IMPORTANT: domain_to_dist()
clips (doesn't extrapolate) queries outside the domains actually
present in `transects` -- on a narrower window (e.g. the zoomed
GIF, which only has domains up to ~25), domain_to_dist(65) would
silently return the position of whatever the real max domain is,
mislabeling it "D65". Guard against that by checking 65 is
actually within this transects table's own domain range first.
```

Above `all_positions_km = dist_km  # one position per domain (order`:

```text
Thin vertical line at EVERY individual domain's position (not
just the labeled ticks above, which only show every 5th/10th),
so each ~500 m domain boundary is visible on the plot itself.
```

##### `lighten_color()`

```text
Return a lighter shade of hex_color, blended toward white by
`amount` (0 = unchanged, 1 = white). Used to give a smoothed-curve
overlay a distinguishable shade of the same base color as its raw/
connected-dots series, without hand-picking a second hex code for
every line.
```

##### `add_geographic_annotations()`

```text
Add community/feature context to a distance-from-groin plot:
shaded town spans (Buxton, Avon), pier markers, and shoal-zone
shading (Wimble Shoals, Avon Shoals) -- so a reader can see what's
actually around a given stretch of coast, not just an abstract
domain number. Positions come from ANN_* config (domain numbers/
ranges), mapped to km via build_domain_to_dist_km() so they land
correctly on whatever x-window this particular plot is showing.

Styling matches the established convention from
HAT_groin_hindcast_1984_2024.py: labels use a blended transform (x in
data coordinates, y in AXES-FRACTION coordinates), which keeps
label height fixed relative to the panel regardless of the y-data
range or whether the axis is inverted (e.g. the GIF), and a white
bbox behind every label so it stays legible even where a groin/pier
line passes directly behind it.

gif_mode=True switches to the bolder styling requested specifically
for the GIF (higher alpha, and a higher zorder than the water/
island fill_between so the shading renders ON TOP of it instead of
being hidden underneath) -- the static plots use the subtler,
hindcast-matching style by default.

The groin itself is NOT re-labeled here -- it's already the
central, uniquely-marked feature on every plot in this script.

Call this AFTER the x-limits are set (same requirement as
set_domain_primary_axis), and safe to call even if none of the
annotated features fall within the current visible window (those
are just skipped).
```

Above `for name, (d_lo, d_hi) in ANN_TOWN_SPANS.items():`:

```text
Town spans (shaded band + label) -- centered by default, but
nudged off-center for any span where a pier sits near enough to
the midpoint to collide with a centered label (e.g. Avon Pier at
domain 26, the middle of Avon's 21-31 span).
```

Above `for name, d in ANN_VILLAGE_LINES.items():`:

```text
Village center lines
```

Above `for name, (d, label_y_frac) in ANN_PIERS.items():`:

```text
Piers
```

##### `annotation_legend_handles()`

```text
Consolidated proxy legend artists for add_geographic_annotations()
-- one entry per LAYER TYPE (not one per named feature), matching
the established convention from HAT_groin_hindcast_1984_2024.py. Append
these to a plot's own legend handles rather than relying on labels
attached to the individual axvspan/axvline calls, which would
otherwise add one legend entry per town/shoal/pier and make an
already-busy legend worse.
```

##### `_iter_axes_real_xy()`

```text
Yield (xdata, ydata) arrays for every REAL data artist on ax --
line plots and scatter points via their offsets, fill_between bands
via their polygon vertices -- while skipping axhline/axvline
reference lines (2-point constant lines using an axes-fraction
blended transform on one axis, not real data).
```

Above `verts = np.concatenate([p.vertices for p in coll.get_paths()`:

```text
fill_between / PolyCollection: real data is the polygon
vertices themselves
```

##### `apply_distance_xlim_and_autoscale_y()`

```text
Clip a distance-from-groin plot's x-axis and set the y-axis from
the data actually VISIBLE within that window -- rather than an
arbitrary fixed clip like the original +/-8 m/yr, or a symmetric +/-
x_max_km that wastes space on the downdrift (short-coverage) side.

x_max_km defaults to PLOT_X_MAX_KM (updrift coverage runs close to
the full island). x_min_km defaults to None, which means: auto-
detect the most negative (downdrift) x-value actually present in
the plotted data, minus a small margin -- i.e. start the axis right
where the real data starts instead of extending all the way to
-x_max_km when downdrift coverage is much shorter than updrift.
Pass an explicit x_min_km (e.g. -10) to override the auto-detect.

y_percentile: the y-range is set from these percentiles of the
visible raw data (not the strict min/max) -- a single extreme
outlier point (which raw, unsmoothed per-transect data can have)
would otherwise stretch the axis far past where nearly everything
else sits, wasting most of the figure on empty space. Pass
(0, 100) to restore strict min/max if you want every point
guaranteed in view.

Call this AFTER all plot/fill_between calls on the axes.
```

##### `load_groin_geometry()`

```text
Load the groin GeoJSON/shapefile and derive:
    reference_x, reference_y - centroid of the NORTHERNMOST
                   individual groin feature -- used as the origin
                   for distance-from-groin (0 km), since updrift =
                   north here
    y_min, y_max - alongshore extent across ALL groin features
                   (southernmost to northernmost), used to shade
                   the whole groin field's footprint on plots
    geometry    - the unioned geometry (for map overlays)
    n_features  - number of individual groin features
Falls back to GROIN_NORTHING_FALLBACK if the file is missing.

If domain_table is provided (see load_domain_reference()), prints
which CASCADE domain the groin actually landed in -- a quick sanity
check against wherever you expect the groin to be (e.g. "the groin
is around domain 6/7"). A mismatch here usually means
GROIN_GEOJSON_PATH points at the wrong file.
```

Above `centroids = gdf.geometry.centroid`:

```text
Per-feature centroid, to find the NORTHERNMOST individual groin
feature specifically -- this becomes the distance-from-groin
origin (0 km), since updrift = north here. The overall span
(miny to maxy, across every feature) is kept separately for
shading the whole groin field's footprint.
```

Above `order = np.argsort(centroids.y.values)`:

```text
Southernmost and "middle" individual groin features too -- used
for lighter reference lines showing where each groin sits within
the field, not just the northernmost origin and the overall span.
```

##### `assign_distance_from_groin()`

```text
Add a signed distance-from-groin column (in meters), measured
ALONG THE SHORE (curvilinear alongshore_m), not as a straight-line
northing difference. Hatteras Island bends sharply through Cape
Point/Buxton; a straight north-south (northing) distance is not
monotonic with true alongshore position there, which corrupts the
x-axis of every distance-from-groin plot near the cape (visible as
a fold-back/box artifact in profile plots). Using alongshore_m --
the same cumulative along-path distance used everywhere else in
this script -- keeps "distance from groin" consistent and
monotonic with alongshore position.

Positive = updrift, negative = downdrift, per UPDRIFT_DIRECTION.
```

Above `if "shore_x" in tx.columns:`:

```text
Locate the groin's own alongshore position: the transect whose
NEAR-SHORE reference point (shore_x/shore_y) is closest to the
NORTHERNMOST groin feature's centroid (groin["reference_x"/"_y"]),
not the centroid of the whole groin field.
```

##### `load_study_area_filter()`

```text
Load and buffer the study area filter polygon.
```

##### `load_shoreline_lines()`

```text
Load a shapefile/GeoJSON of shoreline lines, project, spatially
filter, and parse dates. Returns one row per feature with columns:
    geometry, date, year, source
Multiple features on the same date are kept (needed for intersection).
```

##### `load_coastsat_transects()`

```text
Load CoastSat's own transect layer and use it as the analysis
backbone for EVERY source (wet-dry, NC state, and CoastSat itself)
-- per your direction to use CoastSat's transects for all shoreline
change analyses, replacing the earlier shared-100m-grid approach.

Returns GeoDataFrame with columns:
    transect_id, geometry, origin_x, origin_y, dir_x, dir_y,
    shore_x, shore_y, alongshore_m, domain
```

Above `gdf["shore_x"] = gdf["origin_x"]`:

```text
Near-shore reference point for groin-proximity matching (see
assign_distance_from_groin) -- for this short, beach-hugging
network, origin itself is already near shore.
```

##### `load_100m_grid_transects()`

```text
Load the alternate 100 m transect grid (transects_100m.geojson)
-- evenly-spaced, PARALLEL transects along a straight offshore
reference line -- used ONLY by the GIF's "hybrid"/"grid100m"
transect-source options (see GIF_TRANSECT_SOURCE). Not used
anywhere else in this script; every other analysis always uses
CoastSat's own shore-normal transects.

Returns GeoDataFrame: transect_id, geometry, shore_x, shore_y
(seaward end, on the straight reference line), origin_x, origin_y
(landward end), dir_x, dir_y, alongshore_m (this grid's own
northing-based alongshore position -- since it's evenly spaced
along a STRAIGHT line, alongshore position is just northing raw_offset
from the southernmost transect, no curvilinear ordering needed),
domain (via the same authoritative domain_table used everywhere
else in this script, so domain labels/annotations stay consistent
regardless of which transect source the GIF is using).
```

Above `dx, dy = sx - lx, sy - ly`:

```text
Direction points from the landward origin TOWARD the sea --
matching CoastSat's own convention (origin=landward,
direction=seaward), so "+chainage = accretion" means the same
thing for both transect sources. (seaward - landward), NOT
the reverse -- that sign was backwards before, which is why
grid100m chainage came out inverted.
```

##### `build_hybrid_position_map()`

```text
"hybrid" GIF transect-source: for each CoastSat transect, find
the alongshore position it would have under the 100 m grid's
coordinate system -- via nearest-neighbor matching (by real
distance between each CoastSat transect's shore-side point and the
100 m grid's shore-side points), NOT by reprojecting or re-
measuring anything. This lets CoastSat's own (correctly-measured)
chainage values be repositioned onto an evenly-spaced 100 m
alongshore axis for display, without touching how the shoreline
itself was measured.

Returns a Series indexed by transect_id (str): dist_from_groin_m
in the 100 m grid's coordinate system.
```

##### `reconstruct_and_remeasure_on_grid100m()`

```text
"grid100m" GIF transect-source: reconstruct each chainage
observation into a real-world (x, y) point using its ORIGINAL
(CoastSat) transect's origin/direction, then re-measure chainage
against the NEAREST 100 m-grid transect's own origin/direction --
restoring the shared-100m-grid approach from earlier in this
project. WARNING: reintroduces the known parallel-transect
geometric artifact (see GIF_TRANSECT_SOURCE config comment).

Returns a new DataFrame with "transect_id" replaced by the matched
100 m-grid transect ID and "chainage_m" replaced by the re-measured
value.
```

Above `best_idx = int(tree.nearest(pt))`:

```text
True nearest-neighbor query (no distance cap) -- a buffer-
limited search (the previous approach) can silently drop
points that are legitimately far from any 100m-grid line,
which is exactly what happens near Cape Point/Buxton where
the parallel grid transects diverge most from CoastSat's
true shore-normal ones. This always finds SOME match.
```

##### `_compute_alongshore_positions()`

```text
Cumulative alongshore distance for each transect.

PRIMARY METHOD: if id_col parses cleanly as a numeric sequence, sort
transects by that ID and accumulate REAL geometric distance between
ID-CONSECUTIVE origins. CoastSat generates transect IDs in order
along the digitized reference shoreline, so ID order IS true
alongshore order -- trusting it sidesteps the failure mode below
entirely.

FALLBACK: greedy nearest-neighbor path from the southernmost
origin. This is a heuristic over a raw point cloud and CAN go
wrong: if the path runs out of nearby unvisited points it jumps
across the island instead of along it, silently inflating
alongshore_m for whatever transects get swept up afterward (this
produced a real artifact previously — a handful of transects
reporting distances-from-groin of 100+ km on an island that's
actually ~45 km long). Any large jump is now flagged below
regardless of which method was used, so it's visible in the run
log rather than silently corrupting downstream plots.
```

Above `id_str = gdf[id_col].astype(str).reset_index(drop=True)`:

```text
CoastSat's global transect IDs are often compound,
"{site}-{index}" (e.g. "usa_NC_0031-0160"), since a long
coastline gets split into multiple site chunks -- the
whole string doesn't parse as a plain number, but sorting
by (site prefix, numeric suffix) still recovers true
alongshore order, since CoastSat's own site numbering
runs along the coast too. This is trusted the same way
the plain-numeric case above is.
```

Above `if len(segment_lengths) > 0:`:

```text
Flag suspiciously large jumps so they're visible in the run log
rather than silently inflating alongshore_m for later transects.
```

##### `collect_csv_map()`

```text
Walk one level of subfolders under root_dir and return
{csv_stem: full_filepath} for every CSV found.
```

##### `_detect_coastsat_columns()`

```text
Return (date_col, chainage_col), auto-detecting CoastSat's
(slightly version-dependent) column names, falling back to
"first column is date, second is chainage".
```

##### `load_coastsat_chainage()`

```text
Read every CoastSat per-transect CSV directly, keyed by
CoastSat's own transect ID -- straightforward again now that
CoastSat's transects ARE the analysis backbone, so no point-
reconstruction/re-snapping step is needed.

Returns long-format DataFrame: transect_id, date, chainage_m, source
```

##### `extract_chainage_by_intersection()`

```text
For each (shoreline feature, transect) pair whose bounding boxes
overlap, compute the intersection point and the chainage (signed
distance from transect origin along transect direction).

transect_length_m / valid_range_m default to
TRANSECT_INTERSECTION_LENGTH_M / ±CHAINAGE_MAX_ABS_M (CoastSat's
short, beach-hugging transects).

Returns long-format DataFrame:
    transect_id, date, chainage_m, source
```

Above `result = (result.groupby(["transect_id", "date", "source"])`:

```text
If multiple features from same date intersect same transect
(e.g., NC state segmented shorelines), collapse to per-date mean
```

##### `_choose_intersection_point()`

```text
A single shapely intersection can be a Point, MultiPoint, LineString,
or MultiLineString. Pick the single 2D coordinate best interpreted as
the transect–shoreline crossing:
  - Point: use directly
  - MultiPoint: pick the one nearest to transect origin
  - LineString: pick midpoint (shoreline coincident with transect)
```

##### `build_unified_chainage_table()`

```text
Concatenate all three sources and enrich each row with alongshore
position, domain, and decimal year. Also drops transects with too few
observations to be useful.
```

##### `regress_lrr()`

```text
Simple unweighted OLS linear regression.
Returns dict with slope (m/yr), intercept, r_squared, n, se_slope.
```

##### `compute_preinstall_lrr()`

```text
For each transect, fit LRR using ALL observations dated before
PRE_INSTALLATION_YEAR_CUTOFF (not a fixed list of years), pooled
across every source that has pre-1970 data. Nourishment-affected
observations (1966 Buxton fill, and the 1967 shoreline it biases)
are NOT dropped -- they're included in the regression and simply
FLAGGED (nourishment_zone=True) so a reader can see exactly which
points are nourishment-influenced rather than having them silently
removed.

Returns DataFrame:
    transect_id, alongshore_m, domain, slope, intercept, r_squared,
    n, se_slope, nourishment_zone
```

Above `sub = chainage[chainage["year"] < PRE_INSTALLATION_YEAR_CUTO`:

```text
Filter to everything before the cutoff -- this picks up whatever
pre-1970 years actually exist in the data, rather than a
hand-picked subset.
```

Above `def in_any_nourishment_zone(y):`:

```text
Flag (do NOT drop) observations in a documented nourishment zone,
so downstream plots can mark them distinctly instead of hiding them.
```

Above `control = out[(out["n_obs"] >= 3) & (~out["nourishment_zone"`:

```text
Noise-floor diagnostic, to help justify SIGNAL_ANOMALY_THRESHOLD_M_YR
Far-field (beyond SIGNAL_MAX_SEARCH_DISTANCE_M, i.e. outside where
the groin could plausibly reach) and not nourishment-flagged: this
is "ordinary" pre-install shoreline variability, unrelated to the
groin. Compare SIGNAL_ANOMALY_THRESHOLD_M_YR against this.
```

##### `compute_post_install_lrr()`

```text
For each transect, fit LRR using ALL observations dated on/after
GROIN_INSTALLATION_YEAR (1970-present) -- the direct complement to
compute_preinstall_lrr(). This isolates the post-groin rate
cleanly, which compute_fullperiod_lrr() does NOT do: that function
blends pre- and post-install years into a single slope per
transect, which dilutes the apparent post-groin signal whenever a
transect has pre-install observations mixed in with its post-
install ones. Comparing THIS to the pre-install baseline is a
cleaner before/after contrast than comparing the full-period line
to the pre-install baseline.

Returns DataFrame:
    transect_id, alongshore_m, domain, dist_from_groin_m, n_obs,
    slope_m_yr, intercept, r_squared, se_slope
```

##### `compute_fullperiod_lrr()`

```text
Per transect: single-slope OLS over ALL observations, and a
simple piecewise-linear (two-segment) breakpoint search across
BREAKPOINT_SEARCH_START..END. Returns both fits.
```

##### `_piecewise_breakpoint()`

```text
Search for the two-segment breakpoint year (integer) that
minimizes total residual sum of squares (RSS). Both segments must
have BREAKPOINT_MIN_POINTS_PER_SEGMENT observations to be valid.
```

##### `compute_decadal_lrr()`

```text
For each transect, fit LRR within each fixed, NON-OVERLAPPING
DECADE_LENGTH_YEARS-year bin starting at DECADE_START_YEAR (default:
1960s, 1970s, 1980s, ...) -- replaces a continuously-slid 1-year-
step rolling window. The original version reported one number per year
from 1975-2020, but adjacent years shared ~90% of their underlying
observations (a 10-yr window stepped by 1 yr overlaps 9 of its 10
years with its neighbor), so the resulting series looked smoother
and more continuous than the actual survey frequency supports, and
it was hard to say what was genuinely new information between e.g.
center-year 1994 and 1995. Fixed decades are simpler to read and
report: each number is a distinct, non-overlapping slice of time.

Returns long-format DataFrame:
    transect_id, decade_start, decade_label, n_obs, slope_m_yr, r_squared
```

##### `compute_era_lrrs()`

```text
For each (transect × era in ERAS), fit an OLS linear regression
of chainage vs. decimal year using ONLY observations within that
era's year range. This is the structured complement to piecewise
detection: instead of asking "when did the rate change", it asks
"what was the rate during each pre-defined era?"

extra_periods: optional list of additional (name, y_lo, y_hi)
tuples computed the SAME way but NOT part of the official ERAS
list -- e.g. a narrower display-only sub-period (see
PRE_NOURISHMENT_PERIOD) shown for visual comparison in one plot,
without disturbing the anomaly/signal-extent pipeline built on the
3 official eras.

Returns long-format DataFrame:
    transect_id, alongshore_m, dist_from_groin_m, domain, era,
    y_lo, y_hi, n_obs, slope_m_yr, intercept, r_squared, se_slope
Rows are omitted where n_obs < 3 (regression not meaningful).
```

##### `interp_preinstall_baseline()`

```text
Build the regional pre-install baseline by linearly interpolating
("connecting the dots") between the sorted per-transect pre-install
LRR values themselves -- no binning, no smoothing. Query distances
beyond the data's range are clamped to the nearest end value.

Used for BOTH the anomaly calculation and the reference curve drawn
on the profile plots, so there's exactly one baseline definition,
not two. No LOWESS/kernel smoother is used anywhere in this script.
(An earlier version used a 500 m binned-median step function
instead of interpolating; the sharp box-like jumps at each bin
edge read poorly on the profile plots, so this connects the actual
per-transect points directly instead.)

Returns a function baseline(dist_from_groin_m) -> float (vectorized).
```

##### `compute_decadal_anomaly()`

```text
Enrich the decadal LRR table with distance-from-groin and the
anomaly relative to the pre-install regional baseline at each
transect. Anomaly = LRR - baseline. Positive anomaly means
the transect is accreting (or eroding less) faster than baseline;
negative means the opposite.
```

##### `_contiguous_extent_one_side()`

```text
Walk outward from the groin (distance 0) in bin_w-wide bins,
taking the median anomaly per bin. The signal is considered to
extend through the bin as long as direction*median > thresh,
tolerating up to max_gap_bins consecutive non-exceeding bins before
stopping (absorbs single-bin noise). direction=+1 looks for
anomaly > +thresh (updrift/accretion signal); direction=-1 looks
for anomaly < -thresh (downdrift/erosion signal).
Returns (extent_m, peak_anomaly_within_extent).
```

##### `compute_signal_extent_over_time()`

```text
For each decade, find how far the groin
signal extends updrift (positive anomaly = accretion beyond
baseline) and downdrift (negative anomaly = erosion beyond
baseline).

Method: CONTIGUOUS run from the groin outward. Distances are
binned (SIGNAL_EXTENT_BIN_WIDTH_M) and the run of bins whose
median |anomaly| exceeds SIGNAL_ANOMALY_THRESHOLD_M_YR, starting
at the groin, defines the extent (small gaps up to
SIGNAL_EXTENT_MAX_GAP_BINS are tolerated). This is a deliberate
change from "the single furthest transect anywhere that happens to
cross the threshold" -- that definition let one noisy far-field
transect get reported as the edge of the groin's influence, which
produced a signal-extent time series with implausible sudden km-
scale jumps unrelated to any physical process near the groin.

Returns DataFrame:
    decade_start,
    updrift_extent_m,       (contiguous extent where anomaly > +threshold)
    updrift_peak_anomaly,   (peak anomaly WITHIN that extent)
    downdrift_extent_m,     (contiguous extent where anomaly < -threshold)
    downdrift_peak_anomaly,
    n_updrift_transects, n_downdrift_transects
```

##### `plot_alongshore_lrr_profile()`

```text
LRR (m/yr) vs distance from groin. Every series is per-transect
data connected in distance order -- no binning, no smoothing, no
averaging -- just the real points joined so the trend reads as a
line instead of a loose cloud of dots.
    - Grey: pre-install baseline, per-transect LRR
      (nourishment-flagged transects are included here too, not
      singled out with a different marker, though still flagged
      in the data)
    - Dark grey/black: full-period (1849-2024, every observation
      blended into one slope per transect)
    - Blue: post-install only (1970-2024) -- a cleaner
      before/after contrast against the pre-install baseline than
      the full-period points, since it isn't diluted by any
      pre-install years mixed into the same transect's regression
    - Vertical shaded band: groin footprint
    - Top axis: CASCADE domain number
```

Above `fs = fullperiod_lrr[fullperiod_lrr["slope_full"].notna()].co`:

```text
Full-period per-transect points + connecting line
(still every real point, in order -- nothing averaged/binned)
```

Above `ok_pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &`:

```text
Pre-install per-transect points + connecting line (nourishment-
flagged transects included, same marker as everything else -- still
kept in the regression and in preinstall_lrr's nourishment_zone
column for anyone who wants to filter on it, just not called out
with a different marker here)
```

Above `ax.axhline(0, color="#999", lw=0.5, linestyle=":")`:

```text
Reference lines
```

Above `if groin.get("geometry") is not None:`:

```text
Groin footprint band -- spans the real southernmost-to-northernmost
extent of the groin field, relative to the northernmost-groin
origin (0 km) -- generally NOT symmetric around 0, since 0 sits at
the north end of the field, not its middle.
```

Above `d1_km = domain_dist_km(transects, 1)`:

```text
Clip x-axis to start exactly at domain 1 (downdrift) and run to
PLOT_UPDRIFT_MAX_KM updrift, rather than wherever the data
happens to start/stop
```

##### `_lowess_overlay()`

```text
LOWESS-smoothed (x, y), sorted by x -- a VISUAL AID overlay
only, drawn as a lighter shade of a line's own color, alongside
(never instead of) the raw per-transect points + connecting line.
Returns (x_smooth, y_smooth), or (None, None) if statsmodels isn't
available or there's too little data to smooth meaningfully.
```

##### `plot_era_lrr_profile()`

```text
Alongshore LRR profile with one color per era. Every series is
per-transect data connected in distance order -- no binning, no
averaging -- just the real points joined so the trend reads as a
line. A LOWESS-smoothed overlay (lighter shade of the same color,
ERA_PROFILE_LOWESS_FRAC bandwidth, deliberately tighter than the
shared SMOOTHED_OVERLAY_LOWESS_FRAC so it captures the sharp jumps
right around the groin) is drawn alongside each line as a visual
aid, never replacing the raw connected points.

Includes PRE_NOURISHMENT_PERIOD (1996-2021, pre-2022 nourishment)
as an extra display-only line alongside the full Deteriorated era,
to check whether the 2022 nourishments are skewing that era's
overall rate.
```

Above `periods_to_plot = [(n, lo, hi) for n, lo, hi in ERAS if n !=`:

```text
Post-install eras (official) + the pre-nourishment sub-period:
per-transect points + connecting line, both in one faint call
(light background texture) + smoothed overlay (bold, the primary
signal to read).
```

Above `ax.plot(xs, ys, color=color, lw=2.8, zorder=6,`:

```text
Pre-nourishment overlaps Deteriorated's range and color
family -- dash it so the two smoothed lines don't blur
together where they nearly coincide.
```

##### `_build_decade_periods()`

```text
Build a list of (name, y_lo, y_hi) non-overlapping period tuples
from start_year to end_year in increment_years steps -- e.g.
1970-1980, 1980-1990, ... The LAST period may be shorter than
increment_years if end_year isn't an exact multiple of steps away
from start_year -- that's expected (the data just runs out there),
not an error.
```

##### `plot_decade_lrr_profile()`

```text
How has the groin's effect changed continuously over its
lifespan? Shows the LRR profile for each fixed-length post-install
window in `periods` (see DECADE_PLOT_START_YEAR /
DECADE_PLOT_INCREMENTS_YEARS) as its own line, colored along a
sequential colormap (DECADE_PLOT_COLORMAP) so chronological order
reads directly from color -- a finer-grained complement to the
3-era profile (which only splits pre/functional/deteriorated),
for the specific question of whether effectiveness has been
trending rather than just stepping between two states.
```

Above `colors = [cmap(0.45 + 0.55 * (i / max(n - 1, 1))) for i in r`:

```text
Sample from 0.45-1.0: both 0.15 and 0.28 were still reported too
light for the earliest period's line to read clearly against a
white background. 0.45 gives a solid, clearly-visible medium
green as the lightest color, still distinguishable from the
darkest (1.0) end.
```

##### `_zone_mean_and_label()`

```text
Draw a heavy horizontal segment at a zone's mean rate, plus a
numeric label above/below it -- adapted from
HAT_groin_zone_investigation.py's draw_zone_mean_rates(), using
domain-number zone boundaries (converted to km via domain_to_dist)
instead of that script's alongshore-meters convention.
```

##### `plot_zone_panels()`

```text
Multi-panel alongshore profile, one stacked panel per era --
adapted from HAT_groin_zone_investigation.py's
alongshore_profile_panels.png, but using THIS script's own era
definitions and data (Pre-install, Functional groin, Deteriorated
-- already computed elsewhere, not recomputed here) and a
simplified two-zone structure (ZONE_PANEL_DOWNDRIFT_DOMAINS,
ZONE_PANEL_UPDRIFT_DOMAINS) instead of that script's four-zone
system, since these are the specific zones asked for.

extra_functional_panels: optional list of (name, df) tuples --
e.g. "Functional groin (first 5yr)" -- inserted right after the
full "Functional groin" panel, for a closer look at whether the
groin's early years looked different from its full functional era.

Each panel shows the per-transect points connected in order, a
LOWESS-smoothed overlay (matching this script's own convention
elsewhere, not the reference script's rolling mean), a faint
pre-install reference line (on every panel except its own), and a
mean-rate horizontal segment + label for each zone.
```

##### `plot_era_difference_profile()`

```text
Difference profile: the CHANGE in shoreline change rate between
consecutive eras, computed PER TRANSECT (matched by transect_id) --
no binning of any kind. A transect only appears in a given
difference if it has a valid rate in BOTH eras being compared.
Both comparisons use the same "later minus earlier" convention:
    - Functional groin − Pre-install baseline
      (did the groin change the rate relative to background?)
    - Deteriorated − Functional groin
      (positive = the rate was higher after deterioration than
      while the groin was still functional; negative = the rate
      was higher during the functional era, i.e. effectiveness
      declined)

Points are connected in distance order, same as the other profile
plots. Worth keeping in mind here specifically: a difference of two
independently-noisy per-transect estimates can swing more sharply
transect-to-transect than either estimate alone, so this line may
look jumpier than the era profile's -- that's the data, not a
rendering choice.
```

Above `m1 = pre.merge(func.drop(columns=["dist_from_groin_m"]),`:

```text
Functional − Pre-install: only transects present in BOTH.
func's own dist_from_groin_m is dropped before merging -- it's
the same transect's same position as pre's copy, so keeping both
would just get them suffixed apart (dist_from_groin_m_pre/_func)
instead of leaving one plain dist_from_groin_m column to plot.
```

##### `_diagnose_era_overlap()`

```text
Print how many transects overlap between two era-filtered sets
before merging them for a difference calculation. ids_a/ids_b
should largely be the SAME underlying transect set (just filtered
to whichever era each has enough observations in) -- if one era's
valid-transect count is close to the full transect total, the
overlap with ANY other era's valid set is mathematically guaranteed
to be at least that other set's own size. If the printed overlap is
well below that, something's off in how the two sets are being
matched, not just "this is how the data happens to look" -- this
print exists to make that visible rather than silently producing a
plot that only covers part of the coast.
```

##### `_diagnose_alongshore_gaps()`

```text
After merging two era-filtered tables on transect_id, check for
ALONGSHORE GAPS in the result -- stretches with no matched data at
all. On a connected-line plot, a gap like this shows up as a
misleadingly smooth, dead-straight "bridge" jumping from the last
real point before it to the first real point after it, which can
look like a plotting bug even though it's actually just a gap in
the merged data being connected across.

For each gap found, checks whether name_a's or name_b's FULL
(pre-merge, not-yet-intersected) table has any data in that same
stretch -- so instead of leaving "why is there a gap" a mystery,
this identifies whether ONE era genuinely lacks coverage there
(e.g. pre-install has no historical shorelines pre-dating CoastSat
in that stretch) or whether both sides have data that simply isn't
matching by transect_id (which WOULD point to a real bug worth
investigating further).
```

Above `a_has = int(((full_a[dist_col] / 1000 > lo) &`:

```text
STRICTLY inside the gap (exclusive bounds) -- lo and hi are
themselves merged/intersected points, so they're guaranteed
to exist in BOTH full_a and full_b; inclusive bounds would
always count at least those 2 boundary points and falsely
suggest "both sides have data" regardless of the interior.
```

##### `_plot_single_difference()`

```text
Shared rendering for a single (one-line) era-difference plot --
used by plot_pre_to_functional_difference() and
plot_functional_to_deteriorated_difference(), which each isolate
ONE of the two comparisons already shown together on
plot_era_difference_profile(), as their own focused figure.

`color` is the dedicated DELTA color (see DELTA_COLOR) -- used for
the fill between the two curves, the Δ line, and its smoothed
overlay. It's deliberately NOT the same as either curve_a_color or
curve_b_color, so the difference is never visually confusable with
one of the two raw curves it was computed from.

When curve_a/curve_b are provided, the two ORIGINAL curves being
subtracted are drawn too (thin reference lines in their own era
colors), with the gap between them shaded -- so the difference is
visible both as an explicit Δ line AND as "how far apart these two
curves actually sit", which is what makes a delta concrete.
```

Above `if curve_a_values is not None and curve_b_values is not None`:

```text
The two original curves + shaded gap between them (drawn first,
low zorder, so the Δ line and its smoothed overlay stay on top)
```

##### `plot_pre_to_functional_difference()`

```text
Single-comparison plot: Functional groin era LRR minus pre-
install baseline LRR, matched per transect (both required). Isolates
just this ONE comparison from plot_era_difference_profile() as its
own focused figure -- did the groin change the rate relative to
background?
```

##### `plot_deteriorated_to_preinstall_difference()`

```text
Single-comparison plot: Deteriorated era LRR minus pre-install
baseline LRR (later minus earlier, same convention as the other two
focused difference plots), matched per transect (both required).
Spans the WHOLE historical arc -- before the groin existed at all,
to its current deteriorated state -- skipping over the middle
"Functional groin" era entirely. Positive means the rate is higher
now than it was before the groin was ever built; negative means
the opposite.
```

##### `plot_functional_to_deteriorated_difference()`

```text
Single-comparison plot: Deteriorated era LRR minus Functional
groin era LRR (later minus earlier, matching the same convention as
Functional − Pre-install), matched per transect (both required).
Isolates just this ONE comparison from
plot_era_difference_profile() as its own focused figure -- positive
means the rate was higher after deterioration than while the groin
was still functional; negative means the rate was higher during the
functional era (effectiveness declined).
```

##### `compute_distance_band_series()`

```text
Groups transects into named alongshore distance bands
(DISTANCE_BAND_EDGES_M), separately updrift and downdrift, and
computes the median value_col per band per decade. CSV-only now
(the companion plot was removed as unhelpful) -- kept because the
numbers may still be useful even without a chart.

Returns long-format DataFrame:
    decade_start, side ('updrift'/'downdrift'), band_label,
    band_lo_m, band_hi_m, n, median
```

##### `create_groin_evolution_gif()`

```text
Animated GIF of RAW shoreline position (chainage, m -- no
correction/centering of any kind) through time, zoomed to
domain 1 to PLOT_UPDRIFT_MAX_KM updrift by default (same window as
the static profile plots). Pass window_domain_max (e.g. 25) for a
tighter, domain-bounded window instead -- domain 1 to whatever
km position that domain maps to -- e.g. for a closer look at
change right around the groin rather than the whole island stretch.
output_suffix (e.g. "_zoomed") distinguishes the comparison filename/
frame folder when generating more than one GIF in the same run.
Adapted from the whole-island
version in shoreline_chainage_alldata_evolution.py, but built
on this script's multi-source unified chainage table (wet-dry + NC
state + CoastSat, already merged in `chainage`) instead of CoastSat
alone, so the animation can show the full historical record instead
of starting at 1984.

GIF_TRANSECT_SOURCE controls which transects define alongshore
position and/or measure shoreline position (see its config comment
for the full explanation of "coastsat" / "hybrid" / "grid100m").
domain_table and grid100m are only needed for "hybrid"/"grid100m"
-- pass None (default) when using "coastsat".

CoastSat's transects are shore-normal and run parallel to the
coast, so the chainage values are used exactly as provided -- no
geometric correction of any kind.

GIF_FRAME_MODE picks "date" (default -- one frame per unique
observation date, using every observation the data supports) or
"year" (pools all sources/dates within a calendar year into one
frame -- fuller per-frame spatial coverage, but coarser temporal
resolution). In "date" mode, a single wet-dry or NC-state date, or
an early/partial CoastSat pass, often only covers PART of the
window, so individual frames can look spatially sparse -- that's
real data sparsity on that date, not a bug.

x-axis is CASCADE domain number (primary), with distance from groin
(km) on a secondary top axis. y-axis is fixed to the TRUE min/max
chainage across the entire animation (not a per-frame or percentile
range), so the real shoreline position is always fully in view and
the axis never rescales frame to frame -- and INVERTED, so the
island (smaller chainage, landward) renders at the top and the
ocean (larger chainage, seaward) at the bottom. Everything seaward
of the current shoreline curve is shaded like water
(GIF_WATER_COLOR), everything landward like the island itself
(GIF_LAND_COLOR) -- the curve is the actual land/water boundary at
that alongshore position and moment in time, not just an abstract
line, and the shading updates every frame as that boundary moves.

The current frame's line is colored by structural era (matches
ERA_COLORS elsewhere); the preceding GIF_TRAIL_FRAMES frames are
drawn as a fading grey trail (unshaded, so they don't visually
compete with the current frame's water/island fill). A fixed,
always-visible dashed reference line shows the GIF_REFERENCE_YEAR
(1967) shoreline -- the last shoreline before the groin was built --
on every frame dated after it.

Requires Pillow (pip install pillow). Writes nothing and returns
None if Pillow is unavailable or there's too little data in the
window to animate.
```

Above `if window_domain_max is not None:`:

```text
x_max_km: either the standard PLOT_UPDRIFT_MAX_KM (full-island
window, same as the static plots) or, if window_domain_max is
given, wherever that specific domain maps to -- e.g. a tighter
domain 1-25 window zoomed in on the groin area.
```

Above `d1_km = domain_dist_km(transects, 1)`:

```text
Restrict to the window computed above: domain 1 (downdrift) ──
to x_max_km (updrift) of the groin, in CoastSat's own coordinate
system (transects already has this applied from main()'s Cape
Point exclusion).
```

Above `print("  Repositioning CoastSat transects onto the 100 m gri`:

```text
CoastSat's own (correctly-measured) chainage values, `sub`,
are untouched -- only WHERE each transect sits on the x-axis
changes, via nearest-neighbor matching to the 100 m grid.
```

Above `gx, gy = groin.get("reference_x"), groin["reference_y"]`:

```text
Window/order the 100 m-grid transects the same way: domain 1
to PLOT_UPDRIFT_MAX_KM, converted directly via the domain
table's authoritative boundary -- exact here (no curvilinear
ambiguity), since this grid's alongshore position IS northing.
```

Above `if GIF_FRAME_MODE == "year":`:

```text
One row of raw chainage per frame (date or year)
```

Above `ref_data = sub[sub["date"].dt.year == GIF_REFERENCE_YEAR]`:

```text
Persistent reference baseline: the GIF_REFERENCE_YEAR shoreline ──
Built exactly like any other frame -- median chainage per transect
for that specific year, NaN (gap) where no observation exists that
year -- so it's exactly as sparse/honest as the real coverage.
Drawn on every frame dated after GIF_REFERENCE_YEAR as a fixed
comparison line (not part of the fading trail).
```

Above `coverage_counts = sub.groupby("transect_id").size()`:

```text
Diagnostic: per-transect total observation count across the
WHOLE historical record in this window. A gap that recurs across
many frames at the same alongshore position is usually because
those specific transects simply have much sparser coverage than
their neighbors (a real data characteristic worth knowing about),
not a processing bug -- printed here so that's checkable
directly rather than left a mystery.
```

Above `all_vals = np.concatenate(`:

```text
Shared y-axis: TRUE min/max across every frame in the whole time
period (not a percentile-trimmed range), so the shoreline's real
movement is always fully in view, in every single frame -- never
rescaling between frames, and never clipping a real extreme.
Includes the reference row too.
```

Above `groin_span = None`:

```text
Groin field footprint (southernmost to northernmost groin),
relative to the northernmost-groin origin (0 km) -- generally NOT
symmetric around 0, since 0 sits at the north end of the field.
```

Above `if n_ref > 0 and frame_year > GIF_REFERENCE_YEAR:`:

```text
Persistent reference baseline (drawn every frame after the
reference year, not part of the fading trail)
```

Above `row_series = pd.Series(row)`:

```text
Water (seaward of the shoreline, larger chainage) / island
(landward of it, smaller chainage) shading -- the shoreline
curve itself is the boundary between them. Uses an
INTERPOLATED version of row for the shading only, so missing-
data gaps don't flash white between frames -- the actual line
below still uses the raw row and keeps its honest gaps.
```

##### `write_metadata()`

```text
Write a plain-text metadata file summarizing the run.
```

##### `main()`

Above `domain_table = load_domain_reference()`:

```text
2. Load the authoritative CASCADE domain reference (replaces the
original CAPE_POINT_NORTHING formula) -- loaded before the groin so its
position can be sanity-checked against a known domain.
```

Above `transects = load_coastsat_transects(filter_gdf, domain_table`:

```text
5. Load CoastSat's own transect network -- the analysis backbone
for EVERY source below (wet-dry, NC state, and CoastSat itself),
per direction to use CoastSat's transects for all shoreline
change analyses.
```

Above `domain1_row = domain_table[domain_table["domain_id"] == 1]`:

```text
5b. Exclude transects south of domain 1's TRUE boundary (per the
authoritative domains.json, not just "nearest domain by midpoint"
-- assign_domain_from_northing() has no distance cutoff, so a
transect far south of domain 1 still gets labeled domain=1 simply
because domain 1 is the closest available bucket). This is a
REAL data exclusion, not a display clip: it happens before any
chainage extraction or regression below, so it affects every
curve, every smoothing pass, and the GIF -- not just where the
plots' x-axis starts. Per direction: this stretch is close enough
to Cape Point's sharp curvature that alongshore ordering and
chainage baselines are much less reliable there, and it isn't
useful for this analysis anyway.

IMPORTANT: the cutoff is applied to dist_from_groin_m (alongshore
position, built from the trusted CoastSat ID-order) -- the SAME
coordinate every plot's x-axis actually uses -- not to northing
directly. Northing and alongshore position are exactly the two
things that disagree near Cape Point's sharp curve (that's the
reason this exclusion exists at all), so filtering on northing
while plotting against alongshore position let some transects
through the filter that still showed up south of domain 1 on the
actual plots. The domain-1 northing boundary is converted to an
alongshore-position threshold via a linear fit across ALL
transects (robust to the local Cape Point anomalies, since the
fit is dominated by the well-behaved majority of the island),
rather than by averaging just the handful of transects nearest-
midpoint-labeled "domain 1" (which is itself an unstable quantity
right in the area being excluded).
```

Above `coastsat_chainage = load_coastsat_chainage(transects, COASTS`:

```text
6. CoastSat chainage -- read directly, keyed by CoastSat's own
transect ID (no reconstruction/snapping needed now that CoastSat's
transects ARE the backbone).
```

Above `wet_dry_chainage = extract_chainage_by_intersection(`:

```text
6. Wet-dry & NC state chainage via intersection with the SAME
CoastSat transects (short, ~800 m; default TRANSECT_INTERSECTION_
LENGTH_M / CHAINAGE_MAX_ABS_M).
```

Above `wd_check = wet_dry_chainage.merge(`:

```text
Sanity check: how much wet-dry coverage exists on each side of
the groin?
```

Above `print_crs_summary()`:

```text
Consolidated spatial-reference confirmation for every file loaded
so far (study area, groin, domains, CoastSat transects, wet-dry,
NC state)
```

Above `print(f"\n{'='*72}\nRegional baseline (connect-the-dots, not`:

```text
9. Regional baseline: connects the sorted per-transect pre-install
points directly (linear interpolation), used for the anomaly
calculation specifically -- see interp_preinstall_baseline().
```

Above `postinstall_lrr = compute_post_install_lrr(chainage, transec`:

```text
10b. Post-install-only LRR (1970-present) -- cleaner post-groin
rate than the full-period line, since it isn't diluted by any
pre-install years mixed into the same transect's regression.
```

Above `decadal_anomaly = compute_decadal_anomaly(decadal_lrr, trans`:

```text
12. Decadal anomaly (LRR minus the interpolated baseline)
with distance-from-groin
```

Above `lrr_band_series = compute_distance_band_series(decadal_lrr.m`:

```text
13b. Distance-band time series (curve replacement for the original
transect x year heatmaps)
```

Above `print("\n" + "=" * 72)`:

```text
14. Plots
Per your feedback: only the two profile plots below were actually
useful. The distance-band, anomaly-profile, signal-extent, and
transect-diagnostic PLOTS have been dropped (their underlying CSVs
from steps 12-13b are still written, in case the numbers are
useful later even without a chart).
```

Above `for increment_years in DECADE_PLOT_INCREMENTS_YEARS:`:

```text
Decade-increment LRR evolution: how has the groin's effect
changed continuously since installation, not just stepped
between two states? One plot per increment in
DECADE_PLOT_INCREMENTS_YEARS (e.g. both 5-yr and 10-yr).
Reuses compute_era_lrrs() with the decade windows as
extra_periods (the 3 official eras get computed alongside but
aren't used by this plot -- harmless, avoids a second near-
duplicate per-transect loop function).
```

Above `early_functional_periods = [`:

```text
Zone-panel profile: one stacked panel per era (Pre-install,
Functional groin, Deteriorated), with the downdrift/updrift
analysis zones shaded and a mean-rate label per zone per era --
adapted from HAT_groin_zone_investigation.py's
alongshore_profile_panels.png. Also inserts a panel for each
early-functional-era window in
ZONE_PANEL_FUNCTIONAL_EARLY_WINDOWS_YEARS (e.g. "first 5yr",
"first 10yr" of the groin's functional life), to see whether it
looked most effective early on.
```

Above `grid100m = (load_100m_grid_transects(domain_table)`:

```text
14b. Shoreline evolution GIF -- re-enabled: now shaded as water/
island around the actual shoreline curve, with a y-range fixed to
the full historical extent so nothing goes off-frame at any point
in the animation.
```
