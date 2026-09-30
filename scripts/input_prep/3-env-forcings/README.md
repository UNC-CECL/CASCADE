# 3-env-forcings - sea level, storms and waves

The environmental forcing the hindcast reads: the Duck water-level record, the
relative sea-level trend, and the storm series built from water levels and
WIS waves.

```
1-records/
    HAT_download_water_levels.py    the Duck gauge, hourly, gap-checked
    HAT_storm_record_1984_2021.py   the named-storm timeline figure
2-rslr/
    duck_rslr_analysis.py           the RSLR trend per hindcast window
3-storms/
    historical_storm_creation_v3_HAT.py   the storm series the model reads
    storm_figures.py                      the storm record figures, 1996-2024
    storm_validation/                     the named-storm catalogue and the validator
    from_Hannah/                          an earlier generator (HAT_create_storms.py)
    from_lexi/, from_roya/                colleagues' generators, kept as they were (not restyled)
    superseded_20260922/                  the retired storm_check scripts (see WHY.md)
```

`historical_storm_creation_v3_HAT.py` is also read by
`hatteras_ms/experiments/HAT_storm_max_duration.py`, which compiles four of its
top-level functions by name: keep them top-level and keep their names.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 1-records/HAT_download_water_levels.py

Download the Duck, NC (8651370) hourly water levels from NOAA CO-OPS, gap-checked and completeness-verified.

Notes that were in the code:

```text
HAT_download_water_levels.py
NOAA CO-OPS Water Level Downloader — Duck, NC (Station 8651370)
```

```text
Purpose:
Produce a gap-checked, completeness-verified hourly water level CSV for use
as the `water_levels_file` input to historical_storm_creation_v3_HAT.py.

Output columns match that script's config exactly:
t = datetime (GMT)          -> t_name_water = "t"
v = water level [m NAVD88]  -> water_name   = "v"

```

```text
WHY THIS EXISTS:
The previous Duck CSV was missing the record's three largest events
(Gloria 1985, Halloween 1991, Isabel 2003). Root cause was not the API —
it was that download_noaa_water_levels_resumable.py accepted ANY non-empty
response as a complete month, cached it permanently, wrote the comparison CSV
even when months were missing, and checked completeness at the YEAR level
(< 6000 records) — too coarse to ever see a missing month (8760 - 720 =
8040, well above the threshold). v3's load_data() then dropna()'d the
absent hours in silence, so the storm file lost Isabel with no error.

WHAT CHANGED, AND WHY:
Monthly chunks         Kept from the resumable script — one month (~720
records) is far lighter on CO-OPS than a large
window and times out much less. This was the right
call and is preserved.

Cache reuse            Same folder and same {station}_{datum}_{YYYYMM}.csv
naming, so months already downloaded are NOT re-
fetched. Existing cache carries over.

Cache VALIDATION       NEW. Every month — cached or freshly fetched — is
checked against its expected hour count. A truncated
month is re-fetched rather than trusted forever.
This is the fix for the actual bug.

Short-month handling   A month that stays short across retries is treated
as a genuine gauge outage (not API truncation) and
logged, not retried forever.

Fail loudly            Output CSV is NOT written if any month is empty or
if any checked storm window is uncovered.

Conservative fill      Only gaps <= INTERP_LIMIT_HR are interpolated.
Linear fill across a multi-hour outage at a storm
peak flattens the peak and biases Rhigh low — the
quiet version of the bug we just found.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
```

```text
--- Paths ---
CACHE_DIR intentionally matches download_noaa_water_levels_resumable.py so the
existing cache is reused. Point it at your real cache folder.
The record and its download cache moved into the data tree 2026-09-12;
only this downloader lives under scripts/. Anchored on this file, since
the two literals here were drive-rooted and had never resolved.
```

```text
Derived from BEGIN/END so the filename always states its own span. This is a
guard, not cosmetics: the previous broken CSV was named for the full span, and
writing a different span under that name is how a stale file gets read as fresh.
```

```text
Leave True — this is the bug fix. Set False only
to skip the (cheap, local) recount.
```

```text
--- Storm windows that MUST be covered ---
If the gauge is missing here, Rhigh is wrong in a way nothing downstream sees.
```

```text
Storms known to be uncoverable from this gauge. The gate still fires for anything
NOT listed here — so a new gap can never pass silently — but a documented outage
does not block the run. The reason string is your methods-section note.
```

```text
--- known-short registry -----------------------------------------------------
A month that stays short across all retries is a real gauge outage, not API
truncation. Without a record of that, REVALIDATE_CACHE re-fetches it (5 attempts
with backoff, ~100 s) on EVERY run, forever. This registry remembers the
confirmed count so the month is accepted from cache next time — but re-fetches
if the count ever changes, in case CO-OPS backfills the record.
```

```text
NOTE: raw.index is a DatetimeIndex NAMED 't'. Passing the Index
and a Series into DataFrame({...}) makes pandas adopt the Series'
index as the frame index — giving a frame with BOTH a column 't'
and an index named 't', which makes sort_values('t') ambiguous.
.values strips both, leaving a clean RangeIndex.
```

```text
Re-fetch if the cached month is short — UNLESS we already confirmed this
exact count is all CO-OPS has. Then it's a real outage, so accept it.
```

```text
ignore_index=True discards every part's index at the stitch, so no part can
reintroduce an index level named 't' and make sort_values('t') ambiguous.
```

```text
Fill short gaps only.

CAREFUL: pandas interpolate(limit=N) fills the first N NaNs of EVERY run,
including long ones — a 37 h outage would get 3 fabricated hours at its
leading edge. So interpolate, then restore NaN across every run longer
than the limit. Only wholly-short runs survive as filled.
```

```text
--- stale-comparison guard ---------------------------------------------------
If we are not going to write out_name this run, any pre-existing file with
that name is stale and will be silently read by v3_HAT.py. Quarantine it.
```

<details><summary>Function notes (the original docstrings)</summary>

**`fetch_month()`**

```text
Download one month with retry/backoff, keeping the LARGEST response seen.

Retrying matters because CO-OPS under load returns a truncated month rather
than an error. Keeping the largest response distinguishes a transient
truncation (a later attempt returns more) from a genuine gauge outage
(every attempt returns the same short count).

Returns (DataFrame or None, n_valid).
```

**`download_all()`**

```text
Walk every month: reuse valid cache, re-fetch invalid or missing.
Returns (DataFrame[t, v], status DataFrame).
```

**`check_storm_coverage()`**

```text
Returns (unexpected_missing, acknowledged_missing) — lists of storm names.
Only unexpected_missing should ever block a run.
```

</details>

### 1-records/HAT_storm_record_1984_2021.py

The named storms that affected Hatteras, 1984-2021, on a timeline with a detail sidebar.

From the script's original header:

```text
Named-storm record for Hatteras Island, NC, 1984-2021.

Data sources:
  - NOAA Historical Hurricane Tracks search (60 nm buffer of Hatteras, Dare
    County, NC) -> "start"/"end"/"cat" per HISTORICAL_STORMS below, plus
    max sustained wind speed / min pressure pulled from the same search.
  - "Hatteras, NC Hurricane History Since 1985" -> local landfall/impact
    detail (surge, evacuations, damage) for the more notable storms.
  - Bertha, Fran (1996) and Isaias (2020): NOT in the 60 nm buffer search
    (their tracks passed outside it), but documented as local impacts in
    the hurricane-history doc. Wind/pressure for these three were cross-
    checked against NHC/NWS tropical cyclone reports (source="local"
    below). This distinction is kept in the data for provenance but is no
    longer flagged on the figure itself (previously a dagger + footnote).

Known remaining gaps (not yet added -- flagged to Hannah, not included
without confirmation): Earl 2010, Sandy 2012, Maria 2017, Florence 2018,
and Michael 2018 all appear in the hurricane-history doc but, like
Bertha/Fran/Isaias, are absent from the 60 nm buffer search. Unlike
Bertha/Fran/Isaias their NC landfalls (or, for Sandy/Michael, their
tracks) were far enough from Hatteras that local impact was minor/indirect
per the doc's own text, so they were left out pending a decision on
whether to include them.

Nor'easters / non-tropical winter storms are still not included (no
dataset available for these).

Chart labels: storms with a documented local-impact note get a plain
asterisk (*); the writeup is in the "Storm Details" sidebar, matched by
name/year rather than a numbered reference (numbers next to the year
abbreviation, e.g. Gloria '85, read as confusingly similar to the year
itself).
```

Notes that were in the code:

```text
cat        = max Saffir-Simpson category reached over the storm's full
lifetime (NHC best track), not necessarily at landfall
wind       = max sustained wind speed (mph) corresponding to `cat`
pressure   = min central pressure (mb) corresponding to `cat`
landfall   = True  -> doc confirms direct Hatteras/Dare Co. impact or evac
False -> doc/NHC indicates storm passed by with minor/no
local impact (or made landfall well south)
None  -> no local-impact narrative available (NOAA-only)
source     = "noaa"  -> from the 60 nm buffer search
"local" -> from the hurricane-history doc only; wind/
pressure cross-checked against NHC reports
```

```text
De-clutter: several years have 2-5 storms only weeks apart (e.g. five in
2016). Rather than a running left-to-right push (which can cascade a
whole cluster's positions into a neighboring year), group consecutive
storms whose gap is below MIN_SEP and evenly redistribute *within* each
cluster, centered on the cluster's true mean date. Labels still show each
storm's own year, so nothing displayed becomes inaccurate.
```

```text
Cluster-level spreading can leave the tail of one cluster close to the
head of the next (e.g. Colin '16 vs Hermine '16). Clean up any residual
close pairs with a small, local nudge.
```

```text
Axes formatting (set BEFORE label placement so pixel<->data mapping used
for collision checks below matches the final rendered chart)
```

```text
Collision-aware label placement

Rather than assuming above/below alternation is enough (it isn't once 3+
storms land within the same year), place each label, measure its actual
rendered bounding box, and check it against every box placed so far
(marker diamonds + other labels). If it collides, try the next vertical
tier out. Names are placed most-intense-storm-first, so e.g. Isabel/
Dorian keep their close, prominent labels.

Notes no longer live on the chart at all -- they're too long to fit
without crowding, even with tiering. Instead, storms with a note get a
plain asterisk, and readers find the matching writeup by name in the
chronological sidebar list (no number-matching needed).
```

```text
Sidebar: chronological detail list for every storm with a note.

Positions are anchored to axes-fraction (0, 1) -- the top-left of the
panel -- with a running raw_offset in points, so the physical panel-pixel
height doesn't depend on any data/ylim choice. Each entry's vertical
footprint is measured directly from its actual rendered bounding box
(same renderer-based approach used for the main chart), so header/note
spacing is exact rather than a guessed line-height.
```

```text
Products land in the data tree beside the record they describe
(2026-09-12). These were absolute paths into "input_preperation",
a folder renamed long ago, so neither had resolved since.
```

### 2-rslr/duck_rslr_analysis.py

Relative sea-level rise at Duck, NC: an OLS trend per hindcast window from the monthly gauge record.

From the script's original header:

```text
Relative sea level rise (RSLR) rate for the Duck, NC gauge
(NOAA CO-OPS station 8651370), fitted over each hindcast window.

Reads the NOAA monthly mean sea level file (seasonal cycle removed), fits an
OLS trend inside each window, and writes everything to
data/hatteras_init/3-env-forcings/2-rslr/ (2026-09-15 layout; the folder was rslr/ until 2026-09-18):

    record/   duck_8651370_meantrend.csv       the NOAA download, untouched
    fits/     duck_rslr_rates.csv              ONE ROW PER WINDOW: slope, CI,
                                               n, and the value the site
                                               config carries
              duck_rslr_timeseries_<w>.csv     the monthly record inside the
                                               window with its fitted trend
                                               and residual
    figures/  duck_rslr_full_record.png/.pdf   the record with the windows
              duck_rslr_windows.png/.pdf       one panel per window
              duck_rslr_residuals.png/.pdf     residuals per window
              CAPTIONS.md                      written by caption()

THE RATES FILE IS NEW (2026-09-15). Until then the fitted slopes existed only
as annotations burned onto the figures and as hand-typed literals in
scripts/site_layer/hatteras_site_config.py (HATTERAS_PERIODS[...]["sea_level_rise_rate"],
rounded to 0.001 m/yr). The config still carries those literals -- this
script does NOT feed the model -- but the file is the record they were read
from, and its `config_m_yr` column is the rounded value so the two can be
diffed.

STYLE. Drawn under the house standard (scripts/site_layer/hat_figure_style.py) since
2026-09-15: printed width, Arial, panel letters, nothing on the canvas that
belongs in a caption. Windows are drawn in the vintage pair -- the earlier of
two windows in red, the later in blue -- and the two pairs (1984-2004 with
2004-2024; 1996-2010 with 2010-2024) never share a panel, so the pair rule
holds everywhere. The NOAA full-record trend is the reference green.

Units: all computed rates are in metres per year [m/yr]; the CSVs carry
mm/yr beside them.

Date:   5/4/2026; restyled and split into record/fits/figures 2026-09-15
```

Notes that were in the code:

```text
The gauge record and every product of this script live in the DATA tree;
only the script lives here (2026-09-12). Anchored on this file so it follows
the checkout.
```

```text
--- Windows ---
(start_year, end_year) inclusive of both calendar years, so 1984-2004 fits
on 21 years of monthly values. The end year is the survey that closes the
window, which is why the windows overlap at their joints.
Windows 3 and 4 added 2026-09-11. They OVERLAP windows 1 and 2 on purpose --
these are four hindcast windows over one record, not a partition of it.

Each window has a PAIR INDEX (0 = earlier of its pair, 1 = later) that picks
its colour, and a PAIR that decides which panel of the full-record figure
it is drawn on.
```

```text
One legend for the figure: the record and reference once, then the
window trends of each panel in order.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_noaa_meantrend()`**

```text
Load a NOAA CO-OPS monthly mean trend CSV file.

These files have 4 metadata header lines, a blank line, then a column
header line, then data. Values are metres relative to the station's MSL
datum (the file's own header says so; an earlier version of this script
labelled the axis MLLW, which was wrong).
```

**`fit_linear_trend()`**

```text
OLS trend on Monthly_MSL within [start_year, end_year].

Returns slope, intercept, R², p-value, the 95% CI half-width on the slope,
dense predicted arrays with the CI band on the MEAN, and the subset.
```

**`export_rates()`**

```text
One row per window. `config_m_yr` is the slope at the precision the site
config stores, so a diff against HATTERAS_PERIODS is one column.
```

**`_draw_window()`**

```text
Trend line and its 95% band for one window, on top of whatever record
the caller has already drawn.
```

**`plot_full_record()`**

```text
The whole record, twice: (a) with the 1984-2004 / 2004-2024 pair,
(b) with the 1996-2010 / 2010-2024 pair. The two pairs overlap in time,
so on one panel the four bands would sit on top of each other.
```

**`plot_windows()`**

```text
One panel per window: the monthly record inside it, the trend and its
band. Shared y so the slopes compare by eye.
```

</details>

### 3-storms/from_Hannah/storm_creation/HAT_create_storms.py

An earlier CASCADE storm-file generator for Hatteras: NOAA Duck water levels plus WIS waves.

Notes that were in the code:

```text
HAT_create_storm_file.py
CASCADE-Formatted Storm File Generator — Hatteras Island
```

```text
Description:
Generates a CASCADE-compatible storm file (.npy) from NOAA tide gauge data
and WIS hindcast wave data. Storm events are identified as periods where
the computed Total Water Level (TWL) exceeds a threshold. For each event,
five parameters are extracted and stored in CASCADE's native unit system:

[Year_Index, Rhigh, Rlow, Wave_Period, Duration]

Physics (Stockdon et al., 2006):
Rhigh = eta_obs + R2%          (2% exceedance total runup, in dam)
Rlow  = still water level at storm peak, floored at MHW (in dam)
Wave Period = mean peak Tp during storm (s)
Duration = storm length in hours (h), capped at MAX_STORM_DURATION_HR
Year_Index = calendar year of storm peak minus START_YEAR

Units note: CASCADE uses decameters (dam) internally. All elevations are
divided by 10 before saving. Wave period and duration remain in their
native units (seconds and hours respectively).

```

```text
KEY PARAMETER DECISIONS (informed by two working Outer Banks reference files):

MAX_STORM_DURATION_HR = 36
CASCADE ran successfully with max durations of exactly
36 hours. Hannah's file crashed with durations up to 189 hours.
triggering a C-level access violation (0xC0000005).

MIN_INTER_STORM_GAP_HR = 48
Increasing to 48 hours keeps separate
storm events distinct, producing realistic individual durations.
MAX_STORM_DURATION_HR acts as a backstop for any events that still
exceed the ceiling after the gap change.

STORM_THRESHOLD = 1.7 m
This produces a realistic distribution of forcing intensities.

Rlow = max(water_level_at_peak_TWL, MHW_M) / 10
Rlow is the still water level (tide + surge, no runup or setup) at the
storm's peak intensity moment, floored at MHW. This prevents near-zero
or negative Rlow values. Rlow_m and Rlow_dam in the readable CSV are
always consistent: Rlow_dam = Rlow_m / 10.

MAX_STORMS_PER_YEAR = 5
Highest-Rhigh events are retained when trimming is needed.

Adapted from: Storm_Creation.ipynb
Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
--- Time period ---
Run once per hindcast period, changing the four variables below each time.
Period 1: BEGIN_DATE="19840101", END_DATE="20041231", START_YEAR=1984, OUTPUT_NAME="storms_1984_2004.npy"
Period 2: BEGIN_DATE="20040101", END_DATE="20241231", START_YEAR=2004, OUTPUT_NAME="storms_2004_2024_base.npy"
```

```text
--- NOAA Tide Gauge ---
Duck, NC (Station 8651370) — closest long-record gauge on the Outer Banks.
Downloaded automatically in annual chunks.
```

```text
--- WIS Wave Data ---
CSV from WIS Generic Export (https://wisportal.erdc.dren.mil/), station ST63228.
One file covering the full 1984–2024 period works for both hindcast periods.
```

```text
--- Hatteras site parameters ---
MHW_M must match MHW_ELEVATION in HAT_hindcast_1984_2024_old version.py.
Set at berm crest so only events that reach or overtop the berm are included.
```

```text
FIX: increased from 24 → 48 to prevent separate
nor'easters from merging into 100–189 h mega-events
```

```text
--- Duration cap ---
Hard ceiling on storm duration in the saved .npy file.
Both working reference files (colleague and Benton/Ocracoke) have max
duration of exactly 36 hours. Set here as the cap; any event longer is
truncated. Acts as a backstop for the gap merging fix above.
```

```text
--- Per-year storm cap ---
Retains highest-Rhigh events when a year exceeds MAX_STORMS_PER_YEAR.
Barrier3D's internal limit is ~5 storms/year.
```

```text
--- Surge multiplier ---
Scales non-tidal residual only — tidal signal unchanged:
water_level = eta_A  +  (SURGE_MULTIPLIER * eta_NTR)
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_noaa_water_levels()`**

```text
Download hourly NOAA water levels and decompose into tidal prediction
(eta_A) and non-tidal residual (eta_NTR = storm surge proxy).

Returns DataFrame with columns:
    observed_wl : raw observed water level (m NAVD88)
    eta_A       : tidal prediction (m)
    eta_NTR     : non-tidal residual / storm surge (m)
```

**`load_wis_data()`**

```text
Load WIS hindcast CSV from wisportal.erdc.dren.mil Generic Export.
Expected columns: time, waveHs, waveTp, waveMeanDirection.

Returns DataFrame with columns: Hs (m), Tp (s), WAVD (deg).
```

**`compute_twl()`**

```text
Merge tide and wave data; compute Stockdon runup components and TWL.

Key columns produced:
    water_level : tide + scaled surge, no runup (m NAVD88) — used for Rlow
    Rhigh_m     : full 2% exceedance TWL (m NAVD88)
    TWL         : same as Rhigh_m (used for storm identification threshold)
```

**`identify_storms()`**

```text
Identify discrete storm events as contiguous periods where TWL > threshold.

1. Find all hours where TWL exceeds threshold.
2. Merge events separated by fewer than min_gap_hr hours.
   min_gap_hr = 48 prevents separate nor'easters from merging into
   unrealistically long events that crash Barrier3D.
3. Discard events shorter than min_duration_hr hours.

Returns list of (start_timestamp, end_timestamp) pairs.
```

**`extract_storm_params()`**

```text
For each storm event, extract CASCADE storm parameters.

Rlow (KEY FIX):
    Rlow = max(water_level_at_peak_TWL, mhw_m) / 10
    where water_level = tide + surge only (no runup, no setup).
    Flooring at MHW ensures Rlow is never below the tidal datum.
    Rlow_m and Rlow_dam are always consistent: Rlow_dam = Rlow_m / 10.

Duration (KEY FIX):
    Duration = min(actual exceedance hours, max_duration_hr)
    Both working reference files cap at 36 h. Raw duration before capping
    is preserved in the readable CSV as Raw_Duration_h for QA.

Returns:
    arr      : (N, 5) float64 array: [Year_Index, Rhigh_dam, Rlow_dam, Tp_s, Duration_h]
    readable : DataFrame with human-readable QA columns
```

**`cap_storms_per_year()`**

```text
Limit storms to max_storms per Year_Index, retaining highest-Rhigh events.

Barrier3D's internal per-year array limit is ~5 storms. When trimming,
the most morphologically significant storms (highest Rhigh) are kept.
```

</details>

### 3-storms/historical_storm_creation_v3_HAT.py

Build the historical storm series CASCADE reads for one hindcast window, from Duck water levels and WIS waves.

From the script's original header:

```text
Creating a historical storm series for CASCADE
CASCADE requires an input storm series consisting of the model year the storm
occurs, maximum run-up, minimum run-up, wave period, and storm duration.
Synthetic storms can be created using the multi-variate sea storm model of
Wahl et al. (2016), but this script converts recorded water level and wave
data into a CASCADE storm input.

Input files should be csv files with headers. Any missing data will be filled
with nans and dropped. Water levels can be pulled from a NOAA tide gauge
(https://tidesandcurrents.noaa.gov/) and wave information can be obtained
from a WIS Station (https://wisportal.erdc.dren.mil/).

Input files include:
- water levels file: should have datetimes and water levels in m NAVD88
- wave information file: should have datetime, significant wave height (m),
  and wave period (s) columns.

Input variables include:
- beach slope (to calculate runup)
- berm elevation in m NAVD88
- conversion from m NAVD88 to m MHW for your tidal gauge (check the NOAA
  datums for your gauge)
```

Notes that were in the code:

```text
THE WINDOW IS GIVEN ONCE, AS A PERIOD (2026-09-11). Everything below that
names a year -- the dates, the output directory, the file name -- is derived
from it, so they cannot disagree. In the committed state they did: save_dir
said 1984_2004 while save_name said 2004_2024, a pair that cannot both be
right about what was produced.

The two source paths were drive-rooted literals from before the tree was
reorganised (str(_PATH_REPO / "scripts" / "..."), str(_PATH_REPO / "data" / "hatteras_init" / "3-env-forcings" / "storms" / "...")), so they
resolved only if the interpreter happened to start at the drive root. They
are anchored on this file now, like every other path in the input tree.

python historical_storm_creation_v3_HAT.py --start-year 1996 --end-year 2010
```

```text
Anchored 2026-09-14: this named a home directory, or a tree renamed since.
Rule 5 of ORGANIZATION.md.
```

```text
THE LONG-EVENT RULE (2026-09-28, Hannah adopted "trim24"). Grouping runs
less than weather_grouping hours apart can make one event of 100-190 h, and
the original rule DROPPED any event longer than max_storm_dur -- which at
72 h removed Isabel 2003, March 2018, Florence, Dennis and Nor'Ida. The 72
existed only because the pre-49fd069 Barrier3D crashed on long storms
(experiments/storms-and-overwash/2026-09-28-storm-max-duration). "trim"
keeps every event and cuts one longer than max_storm_dur to the
max_storm_dur hours around its peak TWL; Rhigh, Rlow and the period are taken
from what is kept. The defaults (drop, 72) reproduce the v3_72 files.
```

```text
THE SPLIT RULE (2026-09-29, Hannah adopted "split12"). The 24 h grouping
chains storms a week apart, because the berm is overtopped at most high tides
in an active spell: Edouard + Fran 1996 became one 119 h event and Jose +
Maria 2017 one of 186 h, and the trim then kept only the larger storm, so
Fran and Jose were not in the model. --split-gap H splits each grouped event
wherever the water stays below the berm for >= H hours; a piece shorter than
the minimum duration is folded into the piece before it (after, for the
first), so no hour the grouped event counted is lost. Each piece is then an
event of its own (trimmed, dated by its own start). Omitted, events are not
split and the older series rebuild exactly.
(experiments/storms-and-overwash/2026-09-29-event-splitting)
```

```text
DO NOT EDIT CODE BLOCKS BELOW - all inputs are updated automatically. To
change results, please change your inputs in the block above. However, you
should review the comparison below to ensure it makes sense.
```

```text
loop through the NaN times and check if the difference between times is more than 1 hour
if so, that means this is the transition from one consecutive time period of nans to another
add the end value plus the next start value to our list of datetimes
```

```text
separate out by start and end datetime
start dates will be odd indeces, end dates will be even indeces
```

```text
check for missing data
water levels
```

```text
split_gap: cut each grouped event where the water stays below the berm
for >= split_gap hours; a piece shorter than min_storm_dur joins the
piece before it (after, for the first). See --split-gap above.
```

```text
"trim": keep a long event, cut to the max_storm_dur hours above the
berm centred on its peak (clamped to the event's ends)
```

```text
MODEL STEP 1 IS THE WINDOW'S FIRST YEAR, not the first year that happens
to hold a storm (2026-09-11). Anchoring on the data means a window whose
opening year is quiet silently shifts the WHOLE series one year earlier,
and nothing downstream can detect it: CASCADE indexes storms by step, so
the run would simply spend the wrong year's storms with no error.

window_start_year=None keeps the old data-anchored behaviour, so this
reproduces the committed series exactly where the first year is stormy --
which is the case for both of the periods already on disk.
```

<details><summary>Function notes (the original docstrings)</summary>

**`find_time_gaps()`**

```text
this function finds all values that are NaN and groups them by consecutive time periods
NOTE: this function does not fill or modify missing data. it just informs you where the data is missing
df: dataframe
    has datetime indeces and a column with data to evaluate
col_name: string
    name of col to evaluate in df
```

**`load_data()`**

```text
This function loads the data, combines everything into a single dataframe with the specified start and stop values, 
and drops rows with nans. It returns the combined dataframe.
If data is missing from water level, wave height, or wave period, it will call the function above.

start_time: string
    date to start the storms
end_time: string
    date to end the storms
water_levels_file: string
    file that contains the water levels from NOAA gauge
wis_file: string
    file that contains the wave height and period from WIS gauge
t_name_water: string, optional
    name of the column that contains the datetimes in the water levels file
water_name: string, optional
    name of the column that contains the water levels in the water levels file [m NAVD88]
t_name_wis: string, optional
    name of the column that contains the datetimes in the WIS file    
waveHs_name: string, optional
    name of the column that contains the significant wave heights in the WIS file [m]
waveTp_name: string, optional
    name of the column that contains the wave periods in the WIS file [s]
```

**`create_storms()`**

```text
df_merged: dataframe
    comparison dataframe from the previous function
berm_elevation: float
    berm elevation of your location of interest [m NAVD88]
weather_grouping: int, optional
    if storms occur within the specified limit, they are assumed part of same weather system and are grouped into one event [hrs]
MHW: float, optional
    conversion from m NAVD88 to m MHW for your tidal gauge [m]
min_storm_dur: int, optional
    minimum storm duration that is considered a storm event [hrs]
max_storm_dur: int, optional
    maximum duration to include in storm events [hrs]
save_dfs: bool, optional
    determine whether to save the dataframes as csv and npy files
savedir: string, optional
    directory to save the storm files. if blank, use the current working directory
save_name: string, optional
    name of this storm run
window_start_year: int, optional
    calendar year that becomes model step 1. None anchors on the first
    year that holds a storm, which is what this script did before 2026-09-11.
```

</details>

### 3-storms/storm_figures.py

Storm record figures for the canonical 1996 -> 2010 -> 2024 chain, coloured by HURDAT2 type.

From the script's original header:

```text
Storm record figures for the canonical 1996 -> 2010 -> 2024 chain.
Reads the storm series the hindcast runs on (hat_env_forcings'
DEFAULT_STORM_VARIANT, v3_split12_trim24 since 2026-09-29) for the two windows and
draws them as one 1996-2024 record:

    storm_record_1996_2024.png           (a) events per year by storm type,
                                         (b) every event's Rhigh through time,
                                         coloured by storm type (HURDAT2)
    storm_characteristics_1996_2024.png  the two periods compared: Rhigh
                                         exceedance, seasonality, event length
                                         before trimming, Rhigh vs wave period

Output: data/hatteras_init/3-env-forcings/3-storms/figures/1996_2024/, PNGs at the top,
PDFs + CAPTIONS.md + the labelled-event table under supporting/.

THE WINDOWS OVERLAP BY ONE YEAR. Both summary files hold calendar 2010 (the
generator writes start..end inclusive), but the model loop runs start..end-1,
so a 1996-2010 run never spends its 2010 storms. The record here takes
1996-2009 from the 1996_2010 file and 2010-2024 from the 2010_2024 file, which
is exactly what the two runs consume. The two files' 2010 rows are identical.

UNITS. The summary stores Rhigh/Rlow in decametres above MHW (the generator's
(TWL - MHW) / 10 with MHW = 0.36 m NAVD88). Figures show metres above MHW (add 0.36
for NAVD88); the berm threshold is 1.7 m NAVD88 = 1.34 m MHW.

STORM TYPE (tropical cyclone within 500 km / other high-water event) comes from the NHC
HURDAT2 best tracks in 1-records/hurdat2/; see TC_NEAR_KM below.

    python storm_figures.py                # 1996-2024, the canonical chain
    python storm_figures.py --variant v3_72
```

Notes that were in the code:

```text
STORM TYPE, FROM THE NHC BEST TRACKS (HURDAT2). Two classes:
tropical  a tropical or subtropical fix (TD, TS, HU, SD, SS) within
TC_NEAR_KM of Cape Hatteras inside +/-TC_WINDOW_H of the event's
peak. Checked 2026-09-29 for the labelled storms: same storm at a
+/-72 h window, and each sits next to Hatteras at the peak hour.
other     every other event, "other high-water event". Most are nor'easters
(198 of the 240 at the time peak October-April), but the class is
what is LEFT when no cyclone matches, not a positive
identification, so it is not called nor'easter (Hannah,
2026-09-29).
A distant-cyclone class (500-1,000 km, or a swell match by wave direction up
to 2,000 km) was tried the same day and removed: it rested on thresholds
that could not be checked, so Noel 2007 (589 km) and the swell events are
now "other". Tropical labels are the HURDAT2 names; other events carry no
official name and are unlabelled. Every event's class and nearest cyclone
are written to supporting/storm_types_1996_2024.csv.
```

```text
House colours (Hannah, 2026-09-29): the style's amber C["ADDED"] for the
cyclones and the theme blue for the rest, taken from the style's blue ramp
(SMOOTH_RAMP[1], a step lighter than C_1997) so the amber still stands out
on top of it. Tried and replaced the same day: orange/teal, amber/grey,
plum/sage.
```

```text
The generator's total water level, rebuilt to find each event's peak HOUR
(the summary keeps only start and end): Duck water level + Stockdon (2006)
R2% from WIS Hs/Tp at the generator's beach slope. It reproduces every
event's Rhigh to <1 cm, and load_record() checks that it still does.
```

```text
EVERY tropical event is labelled (Hannah, 2026-09-29), once per storm: a
storm with two events (Ophelia 2005) is named at its higher one.
```

```text
The highest non-tropical event (March 2018) is the second-highest of the
record and a documented nor'easter; it has no official name, so it is
marked by type, not by date (Hannah, 2026-09-29).
```

```text
Label placement: each label tries these offsets (points from its event) and
takes the first that overlaps no placed label, no event dot and nothing
already written on the panel; failing that, the least crowded one. Labels
are placed highest event first, so the big storms get the best positions.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_hurdat2()`**

```text
Tropical/subtropical fixes since `first_year`, with the distance (km)
of each fix from Cape Hatteras.
```

**`label_text()`**

```text
A named storm's label. Nor'easters have no official name and are left
unlabelled (Hannah, 2026-09-29: a month/year label reads as a name).
```

</details>

### 3-storms/storm_validation/HAT_storm_catalog.py

The named-storm record for Hatteras Island, NC: one source for every script that needs it.

Notes that were in the code:

```text
HAT_storm_catalog.py
Named-storm record for Hatteras Island, NC — SINGLE SOURCE OF TRUTH
```

```text
WHY THIS FILE EXISTS:
This catalog previously lived in three places: HAT_storm_record_1984_2021.py
(42 storms, with wind/pressure/landfall/note) and both validator scripts
(39 storms, name/start/end/cat only). They had drifted — BERTHA 1996,
FRAN 1996 and ISAIAS 2020 appeared on the record chart but were never
tested against any storm file. Those three are exactly the source="local"
entries added by hand from the hurricane-history doc, i.e. the ones most
in need of testing.

Import from here everywhere. Do not copy the list.

from HAT_storm_catalog import HISTORICAL_STORMS, storms_in_period

HAT_storm_record_1984_2021.py should also be edited to import from here
and delete its own copy.

```

```text
FIELDS:
name      storm name ("Unnamed" for the 2017 ET event)
start/end HURDAT2 full active period — NOT the Hatteras influence window.
A storm active two weeks may affect Hatteras for 24–48 h.
cat       max Saffir-Simpson category over full lifetime (not at landfall)
wind      max sustained wind (mph) corresponding to `cat`
pressure  min central pressure (mb) corresponding to `cat`
landfall  True  -> documented direct Hatteras/Dare Co. impact or evacuation
False -> passed by, minor/no local impact
None  -> no local-impact narrative available (NOAA-only)
note      local-impact detail, where documented
source    "noaa"  -> NOAA Historical Hurricane Tracks, 60 nm buffer of
Hatteras, Dare County NC
"local" -> hurricane-history doc only; wind/pressure cross-checked
against NHC/NWS tropical cyclone reports. Track fell
outside the 60 nm buffer.

KNOWN GAPS (deliberate, pending a decision):
Earl 2010, Sandy 2012, Maria 2017, Florence 2018, Michael 2018 appear in the
hurricane-history doc but are absent from the 60 nm buffer search, and their
NC landfalls/tracks were far enough from Hatteras that local impact was
minor/indirect per the doc's own text.

Nor'easters and non-tropical winter storms are NOT included — no catalog
available. At Hatteras these are the dominant morphological forcing, so
unmatched events in a TWL-derived storm file are EXPECTED, not errors.
```

<details><summary>Function notes (the original docstrings)</summary>

**`storms_in_period()`**

```text
Return catalog entries whose START year falls in [begin_year, end_year],
with 'start_ts'/'end_ts' Timestamps added.

NOTE: filtering on start year (not overlap) matches the original validator
behaviour. A storm starting 2003-12-28 and ending 2004-01-03 belongs to the
year it started in, so it cannot be double-counted across the two periods.
```

</details>

### 3-storms/storm_validation/HAT_validate_storms.py

Validate a CASCADE storm file against the documented named-storm record, per hindcast period.

Notes that were in the code:

```text
HAT_validate_storms.py
Unified storm-record validation — Hatteras Island
```

```text
Compares a CASCADE storm file against the documented record of named tropical
and extratropical storms near Hatteras Island, NC.

REPLACES: HAT_validate_storms.py + HAT_validate_storms_windowchange.py
(which differed only in period config and raw_offset convention)

```

```text
WHAT CHANGED, AND WHY:

Catalog moved out          The list lived in three files and had drifted:
BERTHA/FRAN 1996 and ISAIAS 2020 were on the
record chart but tested by nothing. Now imported
from HAT_storm_catalog.py. Do not paste it back.

Format auto-detection      Reads BOTH schemas:
- v3 summary (Lexi's): StartTime, EndTime,
calendar_year, Rhigh [dam rel. MHW]
- HAT_create_storms readable: Storm_Start,
Storm_End, Peak_TWL_Time, Rhigh_m [m NAVD88]
The original validator required Peak_TWL_Time, so it
could not read v3 comparison at all.

Datum normalisation        v3 stores (TWL - MHW)/10 dam; HAT_create_storms
stores TWL/10 dam. Everything is converted to
m NAVD88 internally so numbers are comparable
across formats. Get MHW right or the axis lies.

Overlap matching           Default MATCH_MODE="overlap": the model storm
WINDOW must intersect the named storm window
(+/- MATCH_WINDOW_DAYS). Works for every format,
needs no Peak_TWL_Time, and is the better test:
a storm is captured if the model has an event
running at the same time, not if a single
instant lands inside a fuzzy box.
MATCH_MODE="peak" reproduces the original behaviour
where Peak_TWL_Time exists.

Shared-match reporting     A merged mega-event can span several named
storms (the 2003 event runs 09-09 -> 09-19 and
could claim more than Isabel). The original code
silently let the first claimant win while still
reporting every storm as matched. Shared matches
are now counted and printed.

Usage:
Set STORM_FILE / BEGIN_YEAR / END_YEAR below and run.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
```

```text
--- Periods ------------------------------------------------------------------
Each period validates independently and writes to its own subfolder under
validation/ subfolder, so the two never overwrite each other and each can
be checked on its own.

BOUNDARY — DECIDED BUT NOT YET APPLIED.
2004 currently appears in BOTH storm files, with the same 10 events. If
Period 2 initialises from Period 1's final morphological state, those 10 are
forced twice. Agreed fix: Period 1 ends 2003-12-31 23:00, Period 2 begins
2004-01-01 — matching historical_storm_creation_HAT.py (v1), which ended
Period 1 at 2003-12-31 23:00 for exactly this reason.

Config below matches the files AS THEY EXIST TODAY (Period 1 -> 2004), so the
script runs against them unchanged. The CROSS-PERIOD SUMMARY will keep
reporting the 2004 overlap every run — that is deliberate, and it stops when
the input is actually regenerated.

WHEN YOU REGENERATE, change in TWO places:
1. historical_storm_creation_v3_HAT.py -> end_time = '2003-12-31 23:00:00'
save_name = "1984_2003_storms_v3"
2. here -> "file": "1984_2003_storms_v3_summary.csv",  "end": 2003
Note this also drops Alex/Bonnie 2004 from Period 1's catalog slice, so its
denominator goes 22 -> 20 and the capture rate shifts for that reason alone.

Point these at the *_summary.csv files, NOT the bare CASCADE files. The bare
files have only a model-year index and cannot be date-matched.
```

```text
Paths are DERIVED from each period's name to match the folder layout:

<STORM_ROOT>/
1984_2004/
1984_2004_storms_v3_summary.csv      <- input  (SUMMARY_TEMPLATE)
validation/                          <- comparison (VALIDATION_SUBFOLDER)
2004_2024/
2004_2024_storms_v3_summary.csv
validation/

So renaming a period (e.g. 1984_2004 -> 1984_2003) updates the input path AND
the comparison folder together. Override per period with explicit "file" /
"outdir" keys if a run ever sits outside this convention.
```

```text
--- Datum --------------------------------------------------------------------
MHW conversion for the Duck gauge [m]: 0 m NAVD88 = MHW m MHW.
Used to convert v3's (TWL - MHW)/10 dam back to m NAVD88. Must match the MHW
in the storm-creation config or every Rhigh here is shifted.
```

```text
"auto"  -> infer from columns (recommended)
"v3"    -> StartTime/EndTime/calendar_year, Rhigh in dam relative to MHW
"readable" -> HAT_create_storms.py readable CSV, Rhigh_m in m NAVD88
```

```text
--- Matching -----------------------------------------------------------------
"overlap" -> model storm window intersects [named_start - W, named_end + W]
"peak"    -> Peak_TWL_Time falls in that range (needs the readable format)
```

```text
Days of slack either side of the named storm's HURDAT2 active period.
(1) Duck gauge (8651370) is ~40 km north of Hatteras centre — surge timing
differs by hours depending on approach direction.
(2) HURDAT2 reports full lifetime at sea, not the Hatteras influence window.
Recommended 3-5. Larger windows inflate the capture rate; see the sensitivity
sweep printed at the end.
```

```text
--- Output -------------------------------------------------------------------
Each period writes into <STORM_ROOT>/<name>/validation/ :
storm_check_<name>.png     figure
match_table_<name>.csv     per-named-storm match detail
report_<name>.txt          the full console report
```

```text
Beside the storm_check output for the same window, not inside the
model-input folder (2026-09-18; was <window>/validation/).
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_storms()`**

```text
Normalise any supported schema to:
    start_ts, end_ts, peak_ts (NaT if absent), year,
    Rhigh_m, Rlow_m [m NAVD88], duration_h, period_s
```

**`match_storms()`**

```text
For each named storm, find the best-matching model event.
Best = highest Rhigh among candidates. Returns (model, results).
```

</details>
