# Shoreline change templates

These are stripped-down versions of the scripts I use in 5-scr to get shoreline
change rates from CoastSat for Hatteras Island. I took out everything specific
to Hatteras so they can be copied out of this repository and used on another
site. Nothing in them imports from the rest of the repo.

The folders are numbered in the order you run them. Run step 1 first. The two
step-2 scripts both read the lookup that step 1 writes, and they can be run in
either order.

```
1-zone-join/
    transect_zone_join_template.py    transect lines + zone polygons -> lookup
    example/                          35 real transects, 3 zones, expected output
2-shoreline-change/
    shoreline_rates_template.py       time series -> OLS rate (LRR) per zone
    shoreline_endpoint_template.py    time series -> net change (end - start)
    example/                          the same 35 transects' time series, expected output
```

Each `example/` folder is a small piece of Hatteras Island (three of my model
domains, 1.5 km of coast) with the input files, the command to run, and the
output it should produce. I'd run the examples before using your own data. If
your output matches `expected_output/`, your setup is working. The step-2
example uses the lookup produced by the step-1 example, so you can run them
in order.

The step-2 scripts need `pandas`, `numpy`, `scipy` and `matplotlib`. The join
also needs `geopandas`.

The scripts themselves only have short comments. I explain the reasoning
behind the choices here instead.

## How the scripts fit together

All three scripts follow the same conventions:

| | convention | reason |
|---|---|---|
| **window** | `--start-year 1996 --end-year 2010` means 1 Jan 1996 through 31 Dec 2010, including both full years | I treat a period as its calendar years. The scripts never compare against a `YYYY-12-31` timestamp, because that timestamp is midnight and would drop any image taken later on 31 Dec. |
| **transect id** | the time-series filename without `.csv` (e.g. `usa_NC_0032_0021`) | This is the same as the `id` field in CoastSat's transect layer, so the lookup matches without any renaming. I don't use the trailing number on its own because CoastSat restarts the numbering at each site, so `0032_0001` and `0033_0001` would clash. |
| **input** | a folder of CSVs, searched recursively; column 1 is the date, column 2 is cross-shore position in metres | CoastSat downloads come as one subfolder per site. |
| **lookup** | `transect_id, zone_id`, written by the join | If a transect isn't in the lookup, it still gets a rate in the per-transect table but isn't counted in any zone. The join left it out for a reason, and averaging it into a zone anyway would undo that. |
| **sign** | seaward is positive (`SEAWARD_POSITIVE`) | CoastSat chainage increases seaward. |

`../3-rates/coastsat/lrr/coastsat_lrr_standalone.py` uses the same window
convention. It's a forty-line script for when you just want one rate per
transect.

## The steps

### 1. Join transects to zones (`transect_zone_join_template.py`)

Each transect line is reduced to one point (its origin by default). A transect
is kept only if that point falls inside exactly one zone polygon. I don't snap
transects to the nearest zone. If a transect is left out, you can see it in
`transect_zones_problems.csv` and on the map. If it were snapped into the
wrong zone, nothing would show it. The script also reports any zones that
didn't get a transect.

### 2a. Rate (`shoreline_rates_template.py`)

| | step | why |
|---|---|---|
| 1 | load | one time series per transect: a date and a cross-shore position |
| 2 | window | a rate only means something alongside the period it was fitted over; whole calendar years by default, or `--from/--to` (below) |
| 3 | fit | ordinary least squares of position against time; the slope is the rate (LRR) |
| 4 | screen | remove bad measurements before averaging, but not small signals (see below) |
| 5 | group | average the transects that pass into zones |
| 6 | write | a per-zone table, a per-transect table, and a figure |

If you need a window that doesn't start or end on a year boundary, for
example starting the day after a storm, use `--from` and `--to` instead of the
years:

```
python shoreline_rates_template.py --from 2003-09-19 --to 2024-12-31
```

Both days are included in full, so an image taken in the afternoon of the
`--to` day still counts. The output files are named with the two dates
(`rates_per_zone_2003-09-19_2024-12-31.csv`). `--from 1996-01-01 --to
2024-12-31` gives exactly the same numbers as `--start-year 1996 --end-year
2024`. You can't mix the two forms in one run.

These are not the same as the endpoint script's `--start-date/--end-date`.
Those are the centres of two ±6-month averaging windows, not the edges of the
period, which is why the rates flags have different names.

### 2b. Endpoint (`shoreline_endpoint_template.py`)

This does the same load, window, screen, group and write steps, but step 3 is
the mean position at the end minus the mean position at the start. By default
each end is a whole calendar year, so each end averages over one seasonal
cycle. The rate is the change divided by the time between the centres of the
two windows, which is 14.0 yr for 1996–2010. The time between the actual mean
image dates is written next to it as `interval_obs_yr`, but it isn't used to
adjust anything.

`--end-window-years 3` averages more years at each end.
`--start-date/--end-date` uses ±6-month windows centred on specific survey
dates instead. I use that when comparing against a dune line or an aerial
photo from a known date.

### Rate vs endpoint

The LRR and the endpoint answer different questions, so I don't treat them as
two estimates of the same number. The LRR is the trend through every image.
The endpoint is how far the shoreline actually moved. On a coast that moves
steadily they agree. They disagree when something happened early or late in
the window, like a nourishment, a storm, or a sudden jump, and that difference
is worth looking into. I compute both.

Both scripts always write a per-transect table as well as the per-zone one.
That's the only way to tell a real alongshore pattern from one noisy transect
pulling a zone mean, so it's the first thing I check when a result looks odd.

## Testing with synthetic data

Before trusting the scripts, you can build a fake coast with rates you already
know and check that they come back out:

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

The zone rates should come back close to -1.65, -1.05, -0.45 and +0.15 m/yr,
which are the true means of each group of six transects. When I ran it on
2026-09-30 I got -1.65, -1.06, -0.45 and +0.15. The endpoint should give
about the same rates (I got -1.62, -1.06, -0.43, +0.16), and changes of about
28 times those. If yours match, the scripts are working, and anything strange
in your real results is coming from the data or the CONFIG settings.

The join can't be tested this way because it needs real geometry, so
`1-zone-join/example/` is its test. On your own data, look at
`transect_zones_map.png` before using the lookup. If most of your transects
end up as `no zone`, it usually means the zone polygons stop short of the
transect origins (try `--join-point midpoint`), or the two files aren't
actually in the same place.

## Don't screen on p-value

`MAX_P_VALUE` in the rates CONFIG is set to `1.0`, which turns it off, and I'd
leave it that way.

Screening on significance seems like the careful thing to do, but it biases the
result. A transect that is actually stable has a slope near zero, so its trend
isn't significantly different from flat and it gets dropped, while the eroding
transects around it pass. The zone mean then only includes the transects that
moved, so it's pushed away from zero, and nothing in the output tells you this
happened.

You can see this with the synthetic data above. Set `MAX_P_VALUE = 0.05` and
the one transect it drops is `site_0020`, which is the one with a true rate of
exactly 0.0 m/yr.

So I screen for bad measurements (too few observations, or a rate that isn't
physically possible), not for the size of the signal. The endpoint script
works the same way: it screens on the number of images at each end, never on
how big the change is. The main pipeline does this too:
`coastsat_5yr_bins.py` uses `MAX_PVALUE = 1.0` and `MIN_R2 = 0.0` as its
defaults.

## What to change

Change the CONFIG block at the top of each script. If your CSVs are laid out
differently, change `load_timeseries`. It reads the columns by position rather
than by name, because every data source names them differently.

The rest of each script is the method. You can change it, but that changes
what the numbers mean, so write down what you changed alongside the output.

## Assumptions and limitations

- The positions are assumed to already be tidally corrected. CoastSat's
  published time series are (FES tide model and a per-transect beach slope).
  If you ran CoastSat yourself from raw imagery, do the tide correction and
  outlier removal before using these scripts, because they don't do either.
- There is no alongshore smoothing. Zone means are plain means. Smoothing
  across zones (e.g. LOWESS) is a separate step if you want it.
- There is no check on how sensitive the rate is to the window. On Hatteras,
  no window shorter than about 25 years gets back to the long-term rate, so I
  would try more than one window before trusting a short one.

## Differences from my full scripts

The scripts in `../2-transect-frame/` and `../3-rates/` add things that are
specific to my study. These templates leave those decisions to you:

- paths come from `site_layer/hat_observed_rates.py` instead of being typed in
- the 90 numbered model domains, a fixed projected CRS, and a fallback CRS for
  ArcGIS exports (the template refuses a file with no CRS instead)
- the endpoint is anchored on the dune-line survey dates by default, so it
  compares like for like with the dune-line data
- my figure style, plus village bands, structures and shoals on the plots
- other ways of looking at the same data (5-year bins, rates converted to a
  distance, window convergence), each answering a different question

For those, see `../README.md`. If you just need a defensible number from a
folder of CoastSat CSVs, these templates are enough.

---

Hannah A. Henry, Coastal Environmental Change Lab, University of North Carolina at Chapel Hill  
hahenry@unc.edu · version 2026-10-01
