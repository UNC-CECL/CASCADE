# storm_check/1995-10-12_1997-10-12 -- were there big storms around this window mean?

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_storm_check.py`
on 2026-09-29. A check, not an input: nothing here changes
the mean line in the folder above.

The mean line is 1995-10-12 to 1997-10-12, ±1 yr of the 1996 fall East Coast NOAA/NASA ALACE lidar (flown 1996-10-09 to 1996-10-16). The question is whether a
big storm inside that window, or just before it, pulled the mean landward.

## Answer

2 major storm(s) peaked inside the window or within 90 days before it. Dropping the passes in the 90 days after them moves the island-median transect mean by **+1.5 m** (positive: the mean would sit further seaward without the post-storm passes), against a median standard error of a transect mean of ~2.4 m and the 10 m Barrier3D cell. The window held **415 storm-hours** above the berm; of the 2-yr spans of 1984–2024 it ranks 362 of 475 (stormier than 24% of them).

## 1. The major storms, 3 yr either side of the window

**Major** = at or above the median annual maximum of 1984–2024 in
**height** (Rhigh ≥ 2.88 m above MHW, 3.24 m NAVD88) or in
**length** (≥ 86 h above the berm): a level the record reaches in half its
years. Ranks are among the 413 events of 1984–2024 (1 = highest / longest).
Hours are above the berm before the 24 h trim.

| peak | name | type | Rhigh (m MHW) | rank, Rhigh | hours | rank, hours | major by | vs window |
|---|---|---|---|---|---|---|---|---|
| 1992-12-14 | — | other | 3.06 | 17 | 108 | 10 | both | before |
| 1993-08-31 | Emily 1993 | tropical | 2.88 | 25 | 35 | 105 | height | before |
| 1994-11-18 | Gordon 1994 | tropical | 3.18 | 13 | 41 | 89 | height | before |
| 1995-08-16 | Felix 1995 | tropical | 2.98 | 23 | 117 | 5 | both | before |
| 1996-09-01 | Edouard 1996 | tropical | 3.03 | 19 | 77 | 31 | height | inside |
| 1998-02-05 | — | other | 2.27 | 109 | 105 | 12 | length | after |
| 1999-08-30 | Dennis 1999 | tropical | 3.29 | 11 | 116 | 6 | both | after |

The five longest events in the same span (a long nor'easter can move more
sand than a higher, shorter storm):

| peak | name | type | Rhigh (m MHW) | rank, Rhigh | hours | rank, hours | major by | vs window |
|---|---|---|---|---|---|---|---|---|
| 1992-12-14 | — | other | 3.06 | 17 | 108 | 10 | both | before |
| 1995-08-16 | Felix 1995 | tropical | 2.98 | 23 | 117 | 5 | both | before |
| 1995-09-09 | — | other | 2.39 | 76 | 83 | 27 | — | before |
| 1998-02-05 | — | other | 2.27 | 109 | 105 | 12 | length | after |
| 1999-08-30 | Dennis 1999 | tropical | 3.29 | 11 | 116 | 6 | both | after |

Inside the window itself: 14 events, 1 of them major;
highest 3.03 m (1996-09-01).

## 2. Was the window stormy?

415 storm-hours above the berm fell inside the window. Over every
2-yr span of 1984–2024 (stepped monthly) the median is
573 h and the range 162–1174 h,
so this window ranks 362 of 475.

## 3. Did a storm move the mean?

Each transect's window mean recomputed without the CoastSat passes that
fall within **90 days after** a major storm's peak (storms
peaking inside the window, or within 90 days before it). The
shift is (mean without) − (mean as built); positive means the post-storm
passes had pulled the mean **landward**. A transect left with fewer than
10 passes gets no shift.

| storm (peak) | passes removed per transect (median) | shift, median (m) | shift, 5th to 95th pct (m) | largest |shift| (m) | transects left < 10 passes |
|---|---|---|---|---|---|
| Felix_1995_1995-08-16 | 1 | +0.0 | -0.7 to +0.9 | +2.0 | 0 |
| Edouard_1996_1996-09-01 | 3 | +1.4 | +0.0 to +3.2 | +4.3 | 0 |
| all major storms together | 4 | +1.5 | -0.1 to +3.4 | +4.9 | 0 |

Per transect in `storm_check_1995-10-12_1997-10-12_mean_shift.csv`.

## 4. The shoreline through the window

The island-median CoastSat position, each transect relative to its
own window mean, from 3 yr before to 3 yr after. 97 image
dates have at least 50% of the 905 transects
(20 of them inside the window). A storm that moved the shoreline
shows as a step down after its peak that does not recover.
Series in `supporting/storm_check_1995-10-12_1997-10-12_island_series.csv`.

## Caveats

- The storm water levels are the model's estimate (Duck gauge + Stockdon
  R2% from WIS waves, beach slope 0.06), not observations at Hatteras.
- "Major" and the 90-day recovery are choices; re-run with
  `--major-rhigh`, `--major-hours` or `--recovery-days` to test them.
- The shift is the storm's effect *through the passes that followed it*. A
  storm that moved the shoreline for good moves every later pass too, and
  that part is not removable by dropping passes: look at the island series for it.
- Landsat 5 alone before 1999, so the 1996 window has fewer passes to drop.

## Files

| file | what it is |
|---|---|
| `storm_check_1995-10-12_1997-10-12.png` | every storm around the window, major ones circled |
| `storm_check_1995-10-12_1997-10-12_events.csv` | every event in the context span: peak, Rhigh, hours, type, HURDAT2 name and distance, ranks in 1984–2024, major, inside window |
| `storm_check_1995-10-12_1997-10-12_mean_shift.csv` | per transect: window mean, passes removed, mean without, shift |
| `supporting/storm_check_1995-10-12_1997-10-12_island_series.csv` | per image date: transects, median and quartiles of the anomaly, coverage |
| `supporting/storm_check_1995-10-12_1997-10-12.pdf`, `supporting/CAPTIONS.md` | vector figure, caption |

Storm series: `v3_split12_trim24` (`hat_env_forcings.DEFAULT_STORM_VARIANT`).
