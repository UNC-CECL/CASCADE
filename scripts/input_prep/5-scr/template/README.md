# template — a bare version of this process, to share

One standalone script that does what 5-scr does, with nothing specific to
Hatteras Island in it. Meant to be **copied out of this repository** and handed
to someone starting the same kind of work somewhere else.

```
shoreline_rates_template.py   satellite time series -> rates -> table + figure
```

It imports nothing from this repository. Send the single file; it needs only
`pandas`, `numpy`, `scipy` and `matplotlib`.

## What it shows

The six steps, which are the transferable part:

| | step | why it is there |
|---|---|---|
| 1 | **load** | one series per transect: a date and a cross-shore position |
| 2 | **window** | a rate means nothing without the interval it was fitted over |
| 3 | **fit** | OLS of position against time; the slope is the rate |
| 4 | **screen** | drop bad *measurements* before averaging — but not small *signals* |
| 5 | **group** | average surviving transects into zones you reason about |
| 6 | **write** | the per-zone table, the per-transect table, and a figure |

Step 6 writes both tables on purpose. The per-transect file is the only way to
tell a real alongshore signal from one loud transect dragging a zone mean, and
it is the first thing to ask for when a result looks surprising.

## Try it before trusting it

Build a synthetic coast with rates you already know, and check the script
recovers them:

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
    pd.DataFrame({"date": dates, "position_m": pos}).to_csv(
        ts / f"site_{i:04d}.csv", index=False)
    rows.append({"transect_id": f"{i:04d}", "zone_id": f"Z{(i-1)//6 + 1}"})
pd.DataFrame(rows).to_csv("data/transect_zones.csv", index=False)
```

```
python shoreline_rates_template.py --start 1996-01-01 --end 2024-01-01
```

Zone means should come back near **-1.65, -1.05, -0.45, +0.15 m/yr**, which are
the true means of each group of six. If they do, the machinery is right and
anything odd in your real output is in the data or the CONFIG.

## The one trap worth reading before you edit

`MAX_P_VALUE` is off (`1.0`) in the CONFIG block, and it should stay off.

Turning on significance screening looks like rigour and quietly biases the
answer. A genuinely **stable** transect has a true slope near zero, so its
trend cannot be distinguished from flat, so it fails the test and is dropped —
while its eroding neighbours pass. The zone mean is then an average over the
transects that moved, biased away from zero, with nothing in the output saying
so.

You can watch it happen: set `MAX_P_VALUE = 0.05`, run the synthetic set above,
and the single transect it discards is number 20 — the one whose true rate is
exactly 0.0 m/yr.

Screen on evidence of a **bad measurement** (too few observations, a physically
impossible rate). Do not screen on the size of the signal you are measuring.
This is also the position the main pipeline takes: `coastsat_5yr_bins.py` ships
with `MAX_PVALUE = 1.0` and `MIN_R2 = 0.0` and calls that the recommended
default.

## What to change, and what to leave

Change **CONFIG**. Change **`load_timeseries`** if your CSVs are shaped
differently — it reads columns by position, not by name, because every provider
names them something else.

Everything else is the method. Changing it changes what the number means, which
is fine as long as the change is written down next to the output.

## How this differs from the real thing

The production scripts in `../3-rates/` add what a specific study needs, and
each addition is a decision this template does not make for you:

- paths resolved through `site_layer/hat_observed_rates.py` rather than typed
- transects joined to domains by point-in-polygon, with nearest-neighbour
  snapping deliberately disabled (`../2-transect-frame/`)
- a house figure style, village bands, structures and shoals
- other readings of the same record — endpoint change, 5-year bins, rate as a
  distance — each answering a different question about the same series

If you want those, read `../README.md`. If you want to get a defensible number
out of a folder of CSVs this week, this file is enough.
