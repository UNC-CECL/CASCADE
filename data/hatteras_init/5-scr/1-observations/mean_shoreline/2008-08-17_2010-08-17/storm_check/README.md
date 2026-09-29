# storm_check/2008-08-17_2010-08-17 -- were there big storms around this window mean?

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_storm_check.py`
on 2026-09-29. A check, not an input: nothing here changes
the mean line in the folder above.

The mean line is 2008-08-17 to 2010-08-17, ±1 yr of the 2009 USACE NCMP topobathy lidar (CHARTS) (flown 2009-08-10 to 2009-08-24). The question is whether a
big storm inside that window, or just before it, pulled the mean landward.

## Answer

1 major storm(s) peaked inside the window or within 90 days before it. Dropping the passes in the 90 days after it moves the island-median transect mean by **+0.9 m** (positive: the mean would sit further seaward without the post-storm passes), against a median standard error of a transect mean of ~2.8 m and the 10 m Barrier3D cell. The window held **858 storm-hours** above the berm; of the 2-yr spans of 1984–2024 it ranks 78 of 475 (stormier than 83% of them).

## 1. The major storms, 3 yr either side of the window

**Major** = at or above the median annual maximum of 1984–2024 in
**height** (Rhigh ≥ 2.88 m above MHW, 3.24 m NAVD88) or in
**length** (≥ 86 h above the berm): a level the record reaches in half its
years. Ranks are among the 413 events of 1984–2024 (1 = highest / longest).
Hours are above the berm before the 24 h trim.

| peak | name | type | Rhigh (m MHW) | rank, Rhigh | hours | rank, hours | major by | vs window |
|---|---|---|---|---|---|---|---|---|
| 2006-11-22 | — | other | 3.63 | 8 | 43 | 77 | height | before |
| 2007-11-03 | — | other | 3.04 | 18 | 40 | 91 | height | before |
| 2009-11-13 | — | other | 2.85 | 28 | 103 | 13 | length | inside |
| 2010-09-03 | Earl 2010 | tropical | 3.83 | 4 | 46 | 70 | height | after |
| 2010-11-12 | — | other | 3.32 | 10 | 87 | 20 | both | after |
| 2011-08-27 | Irene 2011 | tropical | 3.80 | 5 | 35 | 105 | height | after |
| 2012-10-28 | Sandy 2012 | tropical | 3.64 | 7 | 56 | 53 | height | after |
| 2013-03-09 | — | other | 3.11 | 15 | 128 | 4 | both | after |

The five longest events in the same span (a long nor'easter can move more
sand than a higher, shorter storm):

| peak | name | type | Rhigh (m MHW) | rank, Rhigh | hours | rank, hours | major by | vs window |
|---|---|---|---|---|---|---|---|---|
| 2008-09-24 | — | other | 2.73 | 37 | 84 | 23 | — | inside |
| 2009-11-13 | — | other | 2.85 | 28 | 103 | 13 | length | inside |
| 2010-11-12 | — | other | 3.32 | 10 | 87 | 20 | both | after |
| 2011-11-05 | — | other | 2.59 | 49 | 70 | 39 | — | after |
| 2013-03-09 | — | other | 3.11 | 15 | 128 | 4 | both | after |

Inside the window itself: 31 events, 1 of them major;
highest 2.85 m (2009-11-13).

## 2. Was the window stormy?

858 storm-hours above the berm fell inside the window. Over every
2-yr span of 1984–2024 (stepped monthly) the median is
573 h and the range 162–1174 h,
so this window ranks 78 of 475.

## 3. Did a storm move the mean?

Each transect's window mean recomputed without the CoastSat passes that
fall within **90 days after** a major storm's peak (storms
peaking inside the window, or within 90 days before it). The
shift is (mean without) − (mean as built); positive means the post-storm
passes had pulled the mean **landward**. A transect left with fewer than
10 passes gets no shift.

| storm (peak) | passes removed per transect (median) | shift, median (m) | shift, 5th to 95th pct (m) | largest |shift| (m) | transects left < 10 passes |
|---|---|---|---|---|---|
| event_2009-11-13 | 2 | +0.9 | -0.7 to +2.4 | +4.1 | 1 |

Per transect in `storm_check_2008-08-17_2010-08-17_mean_shift.csv`.

## 4. The shoreline through the window

The island-median CoastSat position, each transect relative to its
own window mean, from 3 yr before to 3 yr after. 122 image
dates have at least 50% of the 906 transects
(27 of them inside the window). A storm that moved the shoreline
shows as a step down after its peak that does not recover.
Series in `supporting/storm_check_2008-08-17_2010-08-17_island_series.csv`.

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
| `storm_check_2008-08-17_2010-08-17.png` | every storm around the window, major ones circled |
| `storm_check_2008-08-17_2010-08-17_events.csv` | every event in the context span: peak, Rhigh, hours, type, HURDAT2 name and distance, ranks in 1984–2024, major, inside window |
| `storm_check_2008-08-17_2010-08-17_mean_shift.csv` | per transect: window mean, passes removed, mean without, shift |
| `supporting/storm_check_2008-08-17_2010-08-17_island_series.csv` | per image date: transects, median and quartiles of the anomaly, coverage |
| `supporting/storm_check_2008-08-17_2010-08-17.pdf`, `supporting/CAPTIONS.md` | vector figure, caption |

Storm series: `v3_split12_trim24` (`hat_env_forcings.DEFAULT_STORM_VARIANT`).
