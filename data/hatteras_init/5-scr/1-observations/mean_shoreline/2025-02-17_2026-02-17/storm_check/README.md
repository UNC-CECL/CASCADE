# storm_check/2025-02-17_2026-02-17 -- were there big storms around this window mean?

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_storm_check.py`
on 2026-10-06. A check, not an input: nothing here changes
the mean line in the folder above.

The mean line is 2025-02-17 to 2026-02-17. The question is whether a
big storm inside that window, or just before it, pulled the mean landward.

## Answer

3 major storm(s) peaked inside the window or within 90 days before it. Dropping the passes in the 90 days after them moves the island-median transect mean by **+5.7 m** (positive: the mean would sit further seaward without the post-storm passes), against a median standard error of a transect mean of ~1.9 m and the 10 m Barrier3D cell. The window held **453 storm-hours** above the berm; of the 1-yr spans of 1984–2024 it ranks 83 of 488 (stormier than 83% of them).

## 1. The major storms, 3 yr either side of the window

**Major** = at or above the median annual maximum of 1984–2024 in
**height** (Rhigh ≥ 2.88 m above MHW, 3.24 m NAVD88) or in
**length** (≥ 86 h above the berm): a level the record reaches in half its
years. Ranks are among the 413 events of 1984–2024 (1 = highest / longest).
Hours are above the berm before the 24 h trim.

| peak | name | type | Rhigh (m MHW) | rank, Rhigh | hours | rank, hours | major by | vs window |
|---|---|---|---|---|---|---|---|---|
| 2022-05-10 | — | other | 2.68 | 39 | 95 | 18 | length | before |
| 2023-03-12 | — | other | 2.21 | 125 | 94 | 19 | length | before |
| 2023-09-15 | — | other | 3.02 | 21 | 73 | 36 | height | before |
| 2023-12-18 | — | other | 3.09 | 16 | 67 | 41 | height | before |
| 2024-09-23 | — | other | 1.95 | 208 | 102 | 14 | length | before |
| 2025-08-21 | Erin 2025 | tropical | 3.63 | 9 | 105 | 12 | both | inside |
| 2025-09-30 | — | other | 2.77 | 33 | 107 | 12 | length | inside |
| 2025-10-12 | — | other | 2.58 | 52 | 111 | 9 | length | inside |

The five longest events in the same span (a long nor'easter can move more
sand than a higher, shorter storm):

| peak | name | type | Rhigh (m MHW) | rank, Rhigh | hours | rank, hours | major by | vs window |
|---|---|---|---|---|---|---|---|---|
| 2022-05-10 | — | other | 2.68 | 39 | 95 | 18 | length | before |
| 2024-09-23 | — | other | 1.95 | 208 | 102 | 14 | length | before |
| 2025-08-21 | Erin 2025 | tropical | 3.63 | 9 | 105 | 12 | both | inside |
| 2025-09-30 | — | other | 2.77 | 33 | 107 | 12 | length | inside |
| 2025-10-12 | — | other | 2.58 | 52 | 111 | 9 | length | inside |

Inside the window itself: 8 events, 3 of them major;
highest 3.63 m (2025-08-21).

## 2. Was the window stormy?

453 storm-hours above the berm fell inside the window. Over every
1-yr span of 1984–2024 (stepped monthly) the median is
285 h and the range 64–647 h,
so this window ranks 83 of 488.

## 3. Did a storm move the mean?

Each transect's window mean recomputed without the CoastSat passes that
fall within **90 days after** a major storm's peak (storms
peaking inside the window, or within 90 days before it). The
shift is (mean without) − (mean as built); positive means the post-storm
passes had pulled the mean **landward**. A transect left with fewer than
10 passes gets no shift.

| storm (peak) | passes removed per transect (median) | shift, median (m) | shift, 5th to 95th pct (m) | largest |shift| (m) | transects left < 10 passes |
|---|---|---|---|---|---|
| Erin_2025_2025-08-21 | 9 | +1.6 | -0.5 to +4.8 | +8.3 | 0 |
| event_2025-09-30 | 13 | +3.5 | -0.5 to +9.0 | +14.6 | 0 |
| event_2025-10-12 | 14 | +3.8 | -0.8 to +10.0 | +14.9 | 0 |
| all major storms together | 17 | +5.7 | +0.1 to +13.0 | +19.4 | 0 |

Per transect in `storm_check_2025-02-17_2026-02-17_mean_shift.csv`.

## 4. The shoreline through the window

The island-median CoastSat position, each transect relative to its
own window mean, from 3 yr before to 3 yr after. 135 image
dates have at least 50% of the 906 transects
(34 of them inside the window). A storm that moved the shoreline
shows as a step down after its peak that does not recover.
Series in `supporting/storm_check_2025-02-17_2026-02-17_island_series.csv`.

## Caveats

- The storm water levels are the model's estimate (Duck gauge + Stockdon
  R2% from WIS waves, beach slope 0.06), not observations at Hatteras.
- "Major" and the 90-day recovery are choices; re-run with
  `--major-rhigh`, `--major-hours` or `--recovery-days` to test them.
- The shift is the storm's effect *through the passes that followed it*. A
  storm that moved the shoreline for good moves every later pass too, and
  that part is not removable by dropping passes: look at the island series for it.
- Landsat 5 alone before 1999, so the 1996 window has fewer passes to drop.
- The storm record (Duck gauge and WIS) ends 2025-12-31, before the window does (2026-02-17): storms of the last 48 days are not counted, and the after-window context is empty.

## Files

| file | what it is |
|---|---|
| `storm_check_2025-02-17_2026-02-17.png` | every storm around the window, major ones circled |
| `storm_check_2025-02-17_2026-02-17_events.csv` | every event in the context span: peak, Rhigh, hours, type, HURDAT2 name and distance, ranks in 1984–2024, major, inside window |
| `storm_check_2025-02-17_2026-02-17_mean_shift.csv` | per transect: window mean, passes removed, mean without, shift |
| `supporting/storm_check_2025-02-17_2026-02-17_island_series.csv` | per image date: transects, median and quartiles of the anomaly, coverage |
| `supporting/storm_check_2025-02-17_2026-02-17.pdf`, `supporting/CAPTIONS.md` | vector figure, caption |

Storm series: `v3_split12_trim24` (`hat_env_forcings.DEFAULT_STORM_VARIANT`).
