# How model years work

Written 2026-10-04, when the 1996 run was found to stop a year short of its data.

## The two numbers a period carries

Each period in `HATTERAS_PERIODS` (`scripts/site_layer/hatteras_site_config.py`) has two year fields:

| Field | What it is | What reads it |
|---|---|---|
| `end_year` | The window **label**. | File and folder names: the storm file (`<start>_<end>_storms_*.npy`), the CoastSat target folder (`3-rates/coastsat/lrr/<start>_<end>/`), run folders and run names (`HAT_<start>_<end>_...`), comparison scripts. |
| `last_model_year` | The last calendar year the model **steps through**. | The run length, the nourishment schedule, the road-event filter. |

Read them through the helpers, not the dict:

```python
from site_layer.hatteras_site_config import last_model_year, run_years
run_years(1996)        # 20 = last_model_year - start_year + 1
```

`last_model_year()` raises if `last_model_year` is not `end_year` or `end_year - 1`, so changing a window's label without its run length fails loudly.

## What one model year is

The run makes one transition per calendar year. In the time loop (`cascade_pipeline/hindcast.py`, `run_cascade_simulation`), step `k = 0 .. run_years - 1` is calendar year `start_year + k`.

- **States.** There are `run_years + 1` saved states. State 0 is the initial topography, taken as 1 January of `start_year`. State `k + 1` is the end of calendar year `start_year + k`. The final state is **1 January of `last_model_year + 1`**.
- **Storms.** Barrier3D picks the storms whose year column equals its `time_index`, which runs 1 .. `run_years`. Storm year `k` is calendar year `start_year + k - 1`. The storm file's `_summary.csv` has a `calendar_year` column to check against.
- **Fills and road events.** A project or event dated year `Y` is applied in the step for calendar year `Y` (`time_index = Y - start_year + 1`). Anything dated after `last_model_year` is skipped by `build_schedule`. Before 2026-10-04 the schedule used `end_year`, so a fill dated in the label year was listed as in the period but never fired.
- **Rates.** `compute_lrr` fits a line through all `run_years + 1` states, spaced one year apart.

## The current periods

**Since 2026-10-05 the periods run DEM to DEM** (advisor plan): calibrate 1996 → 2009 and test 2009 → 2025. The 2010 key became 2009, since the period starts in its DEM's year. Both labels are exclusive. Each run starts on 1 Jan of the DEM year and steps whole years: 13 and 16. That puts the final state about 7.5 months before the DEM date (2009-08-17) or its anniversary; this offset is reported, not corrected. The targets are net change between DEM-centred mean shorelines, not the LRR. Rows for the earlier 1996–2015 / 2010–2026 windows are kept below the table as history.

| Period | `end_year` (label) | `last_model_year` | `run_years` | Final state | Storm file holds | CoastSat target spans |
|---|---|---|---|---|---|---|
| 1996 | 2009 | 2008 | 13 | 1 Jan 2009 | 1996–2009 (one spare year) | start ±1 yr of 1996-10-12; end ±1 yr of 2009-08-17 |
| 2009 | 2025 | 2024 | 16 | 1 Jan 2025 | 2009–2025 (one spare year) | start ±1 yr of 2009-08-17; end 2025-08-17 ±6 months (data stops 2026-01-13) |
| 1984 (legacy) | 2004 | 2003 | 20 | 1 Jan 2004 | 21 years | — |
| 2004 (legacy) | 2024 | 2023 | 20 | 1 Jan 2024 | 21 years | — |

Before 2026-10-05 (1996–2015 and 2010–2026) the labels followed different rules: **1996–2015 is inclusive** (the data runs through December 2015), **2010–2026 is exclusive** (the data stops in early January 2026). That is why the label alone can't set the run length. The legacy 1984 and 2004 periods keep the length they always had; their storm files hold one more year than they use.

Until 2026-10-04 the run length was `end_year - start_year`. That was right for 2010–2026 but made 1996–2015 a 19-year run (1996–2014): the 2015 storms were never used, and the model's rate covered a year less than the CoastSat target it was graded against. Every 1996–2015 run made before that date has the short length.

## Changing or adding a window

1. Decide which calendar years the observations cover. CoastSat window folders are named for the calendar years they include.
2. Set `last_model_year` to the last full calendar year of data, and `end_year` to whatever label the window's files use.
3. Check the storm file holds exactly `run_years` years (`calendar_year` in its summary CSV).
4. Re-solve the edgeBE ends if the run length changed.
