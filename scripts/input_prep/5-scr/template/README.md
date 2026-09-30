# template — a bare version of this process, to share

Standalone scripts that do what 5-scr does, with nothing specific to Hatteras
Island in them. Meant to be **copied out of this repository** and handed to
someone starting the same kind of work somewhere else.

Numbered in run order. Step 1 first; the two step-2 scripts read its lookup
and can run in either order.

```
1-zone-join/
    transect_zone_join_template.py    transect lines + zone polygons -> lookup
2-shoreline-change/
    shoreline_rates_template.py       time series -> OLS rate (LRR) per zone
    shoreline_endpoint_template.py    time series -> net change (end - start)
```

The scripts keep their comments to a minimum; the reasoning behind each
choice is in this file.

None of them imports anything from this repository. Send the files you need:
the two time-series templates need only `pandas`, `numpy`, `scipy` and
`matplotlib`; the join also needs `geopandas`.

They are meant to be used together, and they share one contract:

| | the rule | why |
|---|---|---|
| **window** | `--start-year 1996 --end-year 2010` means **1 Jan 1996 through 31 Dec 2010**, both years whole | a period is its calendar years. Compared on the year, never against a `YYYY-12-31` timestamp, which is midnight and would drop a pass later that day |
| **transect id** | the time-series **filename without `.csv`** (`usa_NC_0032_0021`) | it is exactly the `id` in CoastSat's transect layer, so the join's lookup matches with no renaming. Not the trailing number: CoastSat restarts numbering at every site, so `0032_0001` and `0033_0001` would collide |
| **input** | a folder of CSVs, searched recursively; col 1 date, col 2 cross-shore position (m) | CoastSat downloads arrive one subfolder per site |
| **lookup** | `transect_id, zone_id`, written by the join | a transect **not in it** is fitted and kept in the per-transect table but belongs to **no zone** — the join refused it, and averaging it in somewhere would undo that |
| **sign** | seaward positive (`SEAWARD_POSITIVE`) | CoastSat chainage grows seaward |

The same window convention is what `../3-rates/coastsat/lrr/coastsat_lrr_standalone.py`
uses: the minimal forty-line rate script, for someone who wants one number per
transect and nothing else.

## The steps, which are the transferable part

**1 · Join** (`transect_zone_join_template.py`): reduce each transect line to one
point (its origin, by default), keep it only if it lies inside **exactly one**
zone polygon. No nearest-neighbour snapping: a transect left out is visible in
`transect_zones_problems.csv` and on the map; a transect snapped into the
wrong zone is visible nowhere. It also reports zones that got no transect.

**2 · Rate** (`shoreline_rates_template.py`):

| | step | why it is there |
|---|---|---|
| 1 | **load** | one series per transect: a date and a cross-shore position |
| 2 | **window** | a rate means nothing without the interval it was fitted over |
| 3 | **fit** | OLS of position against time; the slope is the rate |
| 4 | **screen** | drop bad *measurements* before averaging — but not small *signals* |
| 5 | **group** | average surviving transects into zones you reason about |
| 6 | **write** | the per-zone table, the per-transect table, and a figure |

**2 · Endpoint** (`shoreline_endpoint_template.py`): the same load, window, screen,
group and write, but step 3 is *end mean minus start mean*. By default each end
is the whole first and last calendar year (a year of passes averages one
seasonal cycle), and the rate divides by the interval between the two window
centres: 14.0 yr for 1996–2010. The gap between the passes' own mean dates is
written beside it as `interval_obs_yr`, not used to adjust it.
`--end-window-years 3` averages more years at each end;
`--start-date/--end-date` centres ±6-month windows on survey dates instead, for
differencing against a dune line or photo taken on a day.

Rate and endpoint are **different questions**, not two estimates of one
number. The LRR is the trend through every pass; the endpoint is how far the
shoreline actually went. They agree on a steadily moving coast and disagree
when something happened early or late in the window — a nourishment, a storm,
a jump. That disagreement is a finding. Compute both.

Every output writes both tables on purpose. The per-transect file is the only
way to tell a real alongshore signal from one loud transect dragging a zone
mean, and it is the first thing to ask for when a result looks surprising.

## Try it before trusting it

Build a synthetic coast with rates you already know, and check the scripts
recover them:

```python
import numpy as np, pandas as pd, pathlib
rng = np.random.default_rng(0)
ts = pathlib.Path("data/timeseries"); ts.mkdir(parents=True, exist_ok=True)
dates = pd.date_range("1995-01-01", "2025-01-01", freq="30D")
rows = []
for i in range(1, 25):
    rate = -2.0 + 0.1 * i                       # a known alongshore gradient
    yrs  = (dates - dates[0]).days / 365.25
    pos  = 100 + rate * yrs + rng.normal(0, 3, len(dates))
    tid  = f"site_{i:04d}"                      # the filename IS the id
    pd.DataFrame({"date": dates, "position_m": pos}).to_csv(
        ts / f"{tid}.csv", index=False)
    rows.append({"transect_id": tid, "zone_id": f"Z{(i-1)//6 + 1}"})
pd.DataFrame(rows).to_csv("data/transect_zones.csv", index=False)
```

```
python 2-shoreline-change/shoreline_rates_template.py    --start-year 1996 --end-year 2024
python 2-shoreline-change/shoreline_endpoint_template.py --start-year 1996 --end-year 2024
```

Zone rates should come back near **-1.65, -1.05, -0.45, +0.15 m/yr**, which
are the true means of each group of six (verified 2026-09-30: -1.65, -1.06,
-0.45, +0.15). The endpoint should give the same rates within noise
(-1.62, -1.06, -0.43, +0.16) and changes of about 28 times them. If they do,
the machinery is right and anything odd in your real output is in the data or
the CONFIG.

The join has no synthetic set here because it needs geometry; test it on your
own layers and **look at `transect_zones_map.png`** before using the lookup.
Most of your transects in `no zone` usually means the polygons stop short of
the transect origins (`--join-point midpoint`) or the two files are in
different places.

## The one trap worth reading before you edit

`MAX_P_VALUE` is off (`1.0`) in the rates CONFIG block, and it should stay off.

Turning on significance screening looks like rigour and quietly biases the
answer. A genuinely **stable** transect has a true slope near zero, so its
trend cannot be distinguished from flat, so it fails the test and is dropped —
while its eroding neighbours pass. The zone mean is then an average over the
transects that moved, biased away from zero, with nothing in the output saying
so.

You can watch it happen: set `MAX_P_VALUE = 0.05`, run the synthetic set above,
and the single transect it discards is `site_0020` — the one whose true rate is
exactly 0.0 m/yr.

Screen on evidence of a **bad measurement** (too few observations, a physically
impossible rate). Do not screen on the size of the signal you are measuring.
The endpoint template follows the same rule: it screens on passes per end
window, never on the size of the change. This is also the position the main
pipeline takes: `coastsat_5yr_bins.py` ships with `MAX_PVALUE = 1.0` and
`MIN_R2 = 0.0` and calls that the recommended default.

## What to change, and what to leave

Change **CONFIG**. Change **`load_timeseries`** if your CSVs are shaped
differently — it reads columns by position, not by name, because every provider
names them something else.

Everything else is the method. Changing it changes what the number means, which
is fine as long as the change is written down next to the output.

## What these assume and do not do

- **Tidally corrected positions.** CoastSat's published series already are
  (FES tide model, per-transect beach slope). If you ran CoastSat yourself from
  raw imagery, correct for tide and remove outliers **before** these scripts;
  they do neither.
- **No alongshore smoothing.** Zone means are plain means. A LOWESS across
  zones is a separate, later decision.
- **No window-sensitivity check.** A rate depends on the window it was fitted
  over; on Hatteras no window under about 25 years recovers the long-term
  rate. Try more than one window before believing a short one.

## How this differs from the real thing

The production scripts in `../2-transect-frame/` and `../3-rates/` add what a
specific study needs, and each addition is a decision these templates do not
make for you:

- paths resolved through `site_layer/hat_observed_rates.py` rather than typed
- domain numbering and the 90-domain frame, a fixed projected CRS, and a
  CRS fallback for ArcGIS exports (the template refuses a missing CRS instead)
- the endpoint anchored on the dune-line survey dates by default, so it
  differences like for like with the dune product
- a house figure style, village bands, structures and shoals
- other readings of the same record — 5-year bins, rate as a distance, window
  convergence — each answering a different question about the same series

If you want those, read `../README.md`. If you want to get a defensible number
out of a folder of CSVs this week, these files are enough.

---

Hannah A. Henry, Coastal Environmental Change Lab, University of North Carolina at Chapel Hill  
hahenry@unc.edu · version 2026-09-30
