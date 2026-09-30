# lib — what the 5-scr scripts share

Two modules. Neither produces a product, and neither has a data twin, which is
why this folder is not numbered.

```
coastsat_lrr.py   The OLS. load_timeseries, filter_dates, compute_lrr,
                  _empty_lrr — imported by the LRR fit, the 5-year bins, the
                  extension fit and the raw DSAS comparison. One copy so the
                  four cannot drift apart; there were two byte-identical
                  copies of it until 2026-09-22.
                  (Was coastsat_lrr_analysis.py, under CoastSat/.)

scr_paths.py      Where the shared modules live. Import it and every
                  module-bearing folder in 5-scr goes onto sys.path.
```

## Using scr_paths

A script that imports a sibling module carries these two lines, and no others:

```python
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
```

After that, `import rates_figures` or `from coastsat_lrr import compute_lrr`
resolves from anywhere in the tree.

**To move a module, edit its row in `MODULE_DIRS` and nothing else.**
`python lib/scr_paths.py` checks every row against disk and names the ones that
are wrong — which is the thing the twelve hand-built `sys.path` inserts this
replaced could never do. A stale insert does not fail where it is written: it
succeeds against a directory that no longer exists, and the error arrives forty
lines later naming the module instead of the path.

Six modules are registered. `coastsat_lrr` and `scr_paths` live here; the other
four live beside the product they belong to, because they are working scripts
that other scripts happen to import:

| module | lives in | because |
|---|---|---|
| `rates_figures` | `3-rates/` | it draws every 3-rates product |
| `coastsat_lrr_windows` | `3-rates/coastsat/lrr/` | it is the multi-window LRR figure |
| `duneline_endpoint` | `3-rates/duneline/` | it is the dune-line product |
| `coastsat_vs_duneline` | `4-comparisons/shoreline_vs_duneline/` | it is the base comparison |
| `total_change_vs_duneline` | `4-comparisons/shoreline_vs_duneline/` | same folder, same question |

Pulling those four in here would file them away from the product they make,
which is the larger harm. Only code that exists *solely* to be imported
belongs in `lib/`.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### coastsat_lrr.py

The linear regression rate (LRR) for CoastSat transects: the fit every rate script in 5-scr imports.

From the script's original header:

```text
CoastSat Transect LRR (Linear Regression Rate) Analysis
Loads CoastSat time-series CSVs from one or more folders,
filters to a user-defined date range OR a set of specific years,
and computes the Linear Regression Rate (LRR) for every transect.

Expected file naming convention:
    <site>_<transect_id>.csv   e.g.  usa_NC_0033_0002.csv

Expected CSV format (CoastSat standard):
    Column 1 – "dates UTC"     : ISO-8601 datetime string
    Column 2 – "chainage (m)"  : cross-shore distance in metres

Usage
Edit the CONFIG section below, then run:
    python coastsat_lrr.py
```

Notes that were in the code:

```text
Date range filter (inclusive). Used when MATCH_YEARS is empty.
Set either to None to include all dates.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_timeseries()`**

```text
Load a CoastSat time-series CSV.
Returns a DataFrame with columns ['date', 'chainage_m'].
```

**`filter_dates()`**

```text
Clip DataFrame to a continuous [start, end] date range (inclusive).

Use this when you want all CoastSat observations between two dates,
e.g. all imagery from 1997-01-01 to 2019-12-31.

For matching specific USGS shoreline years instead, use filter_to_years().

A date-only `end` ("2023-12-31") includes that WHOLE day. Until
2026-09-30 it was compared as a timestamp, i.e. midnight, so a pass later
that day fell outside a window the docstring called inclusive. No window
in use ends on a year with a 31 December pass (those are 2003 and 2023),
so no stored product changed.
```

**`filter_to_dates()`**

```text
Retain only CoastSat observations that fall within ±window_days of
one or more specific USGS shoreline survey dates.

This is the preferred filter for direct DSAS comparisons because it
anchors the CoastSat window to the actual survey date of each
digitized shoreline, not just the calendar year.  For example, if
your DSAS Period 1997–2019 uses shorelines digitized on
1997-09-15, 2008-04-22, and 2019-10-03, this function will collect
CoastSat imagery within ±30 days of each of those dates, ensuring
seasonal and tidal conditions are as comparable as possible.

Parameters
----------
df : pd.DataFrame
    Output of load_timeseries() — must have a tz-aware 'date' column.
survey_dates : list[str]
    ISO-8601 date strings for each USGS shoreline used in DSAS,
    e.g. ["1997-09-15", "2008-04-22", "2019-10-03"].
window_days : int, optional
    Half-width of the search window in days around each survey date.
    Default 30 (i.e. ±1 month).  Increase if CoastSat has sparse
    coverage in your area (cloud cover, satellite gaps).

Returns
-------
pd.DataFrame
    Filtered copy sorted by date, index reset.
    Column 'survey_date' is added to show which anchor date each
    observation was matched to (nearest anchor wins if windows overlap).

Notes
-----
If the windows for two survey dates overlap, an observation is
attributed to the nearest anchor date.  This avoids double-counting
a single image in the regression.

Examples
--------
# Period 1997–2019 with 3 USGS shorelines, ±30-day window
df_matched = filter_to_dates(
    df_raw,
    survey_dates = ["1997-09-15", "2008-04-22", "2019-10-03"],
    window_days  = 30,
)

# Period 1978–1997 with 3 USGS shorelines, ±30-day window
df_matched = filter_to_dates(
    df_raw,
    survey_dates = ["1978-06-10", "1986-08-01", "1997-09-15"],
    window_days  = 30,
)
```

**`filter_to_years()`**

```text
Backwards-compatible wrapper.  Prefer filter_to_dates() for new runs.
If window_days > 0, anchors windows on Jan 1 of each year.
If window_days == 0, keeps any observation whose calendar year is in the list.
```

**`compute_lrr()`**

```text
Compute Linear Regression Rate (LRR) from a time-series DataFrame.

Returns a dict with:
    lrr_m_yr   – slope in m/yr
    r_squared  – R² of the regression
    p_value    – p-value of the slope
    unc_m_yr   – 95 % confidence interval half-width (m/yr)
    n_obs      – number of observations used
    start_date – earliest date in filtered series
    end_date   – latest date in filtered series
```

**`inspect_transect()`**

```text
Load, filter, compute LRR, and plot a single transect CSV.

Pass match_years to use year-matching instead of a continuous range.

Examples:
    # Full date range
    inspect_transect("data/usa_NC_0033_0002.csv", "1997-01-01", "2019-12-31")

    # Year-matched (DSAS-comparable)
    inspect_transect("data/usa_NC_0033_0002.csv", match_years=[1997, 2019])
```

</details>

### scr_paths.py

Where the 5-scr modules live: import this and any 5-scr module imports by name.

From the script's original header:

```text
Where the 5-scr modules live -- the one place that knows.

WHY THIS EXISTS
    Twelve scripts under 5-scr import a module from a SIBLING folder:
    rates_figures for the house drawing, coastsat_vs_duneline for the chainage
    loader, coastsat_lrr for the OLS fit. Until 2026-09-22 each one built that
    folder's path by hand --

        sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr"
                               / "coastsat_vs_duneline"))

    -- so the folder layout was written down in twelve places, in five
    different spellings, and the reorganisation on 2026-09-22 would have had to
    edit all twelve. Worse, a stale one does not raise where it is written: the
    insert succeeds against a directory that no longer exists, and the failure
    surfaces forty lines later as `ModuleNotFoundError: coastsat_vs_duneline`,
    naming the module rather than the path that is wrong.

    Rule 6 of ORGANIZATION.md: a location is decided once, in a resolver, and
    everything else asks. This is that resolver for 5-scr's own modules, as
    site_layer/hat_observed_rates.py is for its data.

USAGE  -- the two lines every 5-scr script that imports a sibling carries:

        sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
        import scr_paths  # noqa: E402,F401

    The import is the whole point: bringing the module in runs the loop at the
    bottom, which puts every module-bearing folder on sys.path. After it, a
    plain `import rates_figures` resolves no matter which folder the importing
    script sits in.

MOVING A MODULE
    Edit its row in MODULE_DIRS below. Nothing else changes. A name whose
    folder has gone missing is reported by check() with the path that is
    wrong, which is the thing the old hand-built inserts could never say.
```

Notes that were in the code:

```text
Every module another 5-scr script imports BY NAME, and the folder holding it.
The comment is who imports it, so a move can be checked against real callers.
```

```text
Importing this module is what does the work. dict.fromkeys keeps the folders
unique and in declaration order; two modules share a folder.
```

<details><summary>Function notes (the original docstrings)</summary>

**`check()`**

```text
Report any module in MODULE_DIRS whose file is not where it is claimed.

Returns a list of complaint strings, empty when the table is accurate.
Run as `python lib/scr_paths.py` after moving anything.
```

</details>
