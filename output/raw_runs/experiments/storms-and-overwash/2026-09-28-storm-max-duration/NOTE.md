# Storm series maximum duration: why longer series failed, and what they do (2026-09-28)

**Question (Hannah).** The storm series keeps events of 8–72 h. Longer limits were abandoned because "the barriers kept drowning and the simulation would end". Why? And what do the longer series do now?

**Background.** The builder (`historical_storm_creation_v3_HAT.py`) groups hours above the berm into one event when they fall less than 24 h apart. It then DROPS any event longer than the limit; it does not shorten it. At 72 h that removes 29 events from 1996–2024, among them:

- Isabel 2003: 133 h above the berm, Rhigh 5.22 m MHW, the largest event in the record
- the March 2018 nor'easter: 4.30 m
- Florence 2018: 3.61 m
- November 2010: 3.32 m
- Dennis 1999: 3.29 m
- Nor'Ida 2009

The pipeline figure shows it: `output/figures/pipeline/3-storms/storm_construction_steps.png`.

**Nothing in the main code changed.**

- **Storm variants.** They come from the builder's own functions, read out of its source with `ast`, and are written to `storms/` here, never to `hindcast_storms/`. The 72 h rebuild equals the committed series exactly in both windows.
- **Runs.** Each is the unchanged hindcast runner, with that period's storm file pointed at the variant inside the run's own process. Barrier3D is the current version (49fd069).

Driver: `scripts/hatteras_ms/experiments/HAT_storm_max_duration.py` (`build`, `run`, `cause`, `diagnose`, `compare`).

Variants (storm count, 1996–2010 / 2010–2024):

| variant | rule | storms |
|---|---|---|
| 72 | the committed series | 139 / 159 |
| 96 | drop events over 96 h | 142 / 171 |
| 120 | drop events over 120 h | 147 / 174 |
| 240 | drop events over 240 h (identical to no limit; the longest event is 193 h) | 150 / 178 |
| 72trim | keep every event, cutting the longer ones to the 72 h around their peak | 150 / 178 |

## Why the longer series failed: the route_overwash bug, not the storms

- **Today's model:** all 16 variant runs finish (both windows, natural and full_management, 96/120/240/72trim). No domain drowns (`tables/drowning.csv`).
- **The cause test (`cause`):** the same series on Barrier3D **before** the 2026-09-24 route_overwash axis-swap fix (commit ce36866, a detached worktree at `../Barrier3D-prefix-ce36866`):

| window | series | pre-fix Barrier3D |
|---|---|---|
| 1996–2010 | 72 | finishes |
| 1996–2010 | 240 | **crashes silently in model year 8 (2003, Isabel)**, exit code 0xC0000005 (access violation) |
| 2010–2024 | 72 | finishes |
| 2010–2024 | 240 | finishes |

The unfixed router read `Elevation[TS, i, d+1:d+10]`, which goes out of bounds whenever `i >= rows`. A long storm routes through far more steps, so Isabel reaches the bad read and the process dies with no message. That matches "the simulation would end". The 72 h limit worked only because it happened to drop the storms long enough to hit the bug.

The older runs also used the ÷10 offset and earlier topography, so a separate drowning then cannot be ruled out. But the crash is reproduced, it lands exactly on Isabel, and it is gone on the fixed code.

## What the longer series do (`tables/comparison.csv`)

Each run against its matrix control (72 h), edgeBE, option A waves:

| window, scenario | series | observed-overwash hit rate | mean shoreline change | interior RMSE (m/yr) |
|---|---|---|---|---|
| 1996–2010 natural | 72 | 16% | −10.2 m | 1.02 |
| | 120 | 75% | −20.4 m | 1.16 |
| | 72trim | 80% | −44.5 m | 3.18 |
| | 240 | 80% | −51.1 m | 3.83 |
| 1996–2010 managed | 72 | 13% | −4.0 m | 1.05 |
| | 120 | 70% | −7.1 m | 1.07 |
| | 72trim | 73% | −16.3 m | 1.58 |
| | 240 | 73% | −20.6 m | 1.96 |
| 2010–2024 natural | 72 | 87% | −37.5 m | 4.14 |
| | 72trim | 89% | −74.5 m | 7.16 |
| 2010–2024 managed | 72 | 78% | −5.7 m | 2.25 |
| | 72trim | 82% | −19.9 m | 3.41 |

The hit rate dates each run's overwash by its own storm file, with the 7-day grace of the observed record.

- **Overwash.** The missing storms are what the imagery shows: the 1996–2010 hit rate goes from about 15% to about 70–80%.
- **Shoreline change.** The same storms make the shoreline retreat much more, so skill against CoastSat drops. The edge rates and wave settings were calibrated on the 72 h series, so those numbers are not a fair test of the longer series; a longer series would need the ends re-solved.
- **The 120 h series** (keeps Dennis 1999 and November 2010; Isabel is still dropped at 133 h) gets most of the overwash gain in the managed 1996 run at almost no RMSE cost (1.05 → 1.07).
- **The 72trim series** keeps every storm at the 72 h the model was built around, and adds less retreat than 240 h.

**Status: record, decision pending (Hannah).** Which series to adopt, and whether to re-solve the ends on it, is Hannah's call. Nothing in `hindcast_storms/` or the builder has changed.

Files: `storms/`, `runs/<variant>_<scenario>/`, `runs/prefix_<variant>_natural/` (pre-fix cause test), `logs/` (with `launches.jsonl`), `tables/drowning.csv`, `tables/comparison.csv`.
