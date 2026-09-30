# Splitting back-to-back storms (2026-09-29)

**Question (Hannah: "test splitting the events").** The storm builder groups above-berm spells less than 24 h apart into one event. At Hatteras the berm is overtopped at most high tides during an active spell, so storms a week apart chain together:

- Edouard and Fran 1996 became one event of 119 h above the berm.
- Jose and Maria 2017 became one of 186 h.

The adopted `v3_trim24` series keeps only the 24 h around each merged event's highest peak, so **Fran 1996 (2.62 m MHW, 33 h above the berm) and Jose 2017 (2.81 m, 35 h) are not in the model**. Across 1996–2024, 56 spells of ≥8 h above the berm are cut out of merged events (41 of them peak above 2 m). Does putting them back change the model?

**Setup** (6 runs, all complete, on `hatteras/adopted@d6546c7`; driver `scripts/hatteras_ms/experiments/HAT_storm_event_splitting.py`):

- **Storms, every one trimmed to 24 h as adopted:**
  - `trim24`: the adopted series (the control). It was rebuilt from the builder's own functions and is identical in both windows.
  - `g12`: the builder with a 12 h grouping gap.
  - `split12`: the 24 h grouping kept, each system split where the water stays below the berm ≥12 h, and pieces shorter than 8 h folded into a neighbour.
- **Runs:** managed (`full_management`), edgeBE, both windows, the site config's end rates (not re-solved), the unchanged runner with the storm file swapped in its own process.
- **Scoring:** as in `../2026-09-28-trim-length-adopted/`.

**The series** (`tables/series.csv`; the calendar years each run spends):

| window | series | events | storm-hours | Fran 1996 | Jose 2017 |
|---|---|---|---|---|---|
| 1996–2009 | trim24 | 132 | 2,524 | absent | |
| | g12 | 125 | 2,460 | 2.62 m | |
| | split12 | 138 | 2,626 | 2.62 m | |
| 2010–2024 | trim24 | 178 | 3,345 | | absent |
| | g12 | 174 | 3,260 | | 2.81 m |
| | split12 | 191 | 3,542 | | 2.81 m |

g12 recovers Fran and Jose but has FEWER events: a merged storm's tidal fragments shorter than 8 h become separate events and fall under the minimum-duration rule. split12 keeps every hour the adopted series counts and adds the second storms.

**The runs** (`tables/scores.csv`; "events" there counts the whole .npy, including the 2010 year the 1996 run does not spend):

| window | storms | PSS | POD | POFD | timing r | space r | RMSE (m/yr) | bias (m/yr) | overwash total (m³/m) | mean shoreline change (m) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1996–2010 | **trim24** | 0.60 | 0.78 | 0.19 | 0.84 | 0.63 | 1.17 | +0.06 | 1,614 | 5.71 |
| | g12 | 0.59 | 0.78 | 0.19 | 0.83 | 0.64 | 1.17 | +0.06 | 1,615 | 5.73 |
| | split12 | 0.60 | 0.78 | 0.19 | 0.84 | 0.63 | 1.17 | +0.06 | 1,614 | 5.78 |
| 2010–2024 | **trim24** | 0.17 | 0.57 | 0.40 | 0.38 | 0.21 | 2.07 | −1.36 | 3,575 | 5.08 |
| | g12 | 0.17 | 0.57 | 0.40 | 0.38 | 0.21 | 2.07 | −1.37 | 3,670 | 5.25 |
| | split12 | 0.17 | 0.57 | 0.40 | 0.38 | 0.21 | 2.07 | −1.37 | 3,648 | 5.27 |

**What it shows.**

1. **Putting the hidden storms back changes neither where nor when the model overwashes.** PSS, POD, POFD and the timing and space correlations are the same to two decimals in both windows. Fran and Jose (2.6–2.8 m MHW) do not overtop the per-cell dune crests (about 5 m) anywhere the imagery was checked.
2. **The shoreline barely moves.** Interior RMSE is unchanged, the bias shifts by ≤0.01 m/yr, and mean shoreline change differs by at most 0.2 m over 14–15 yr.
3. **Overwash volume rises slightly, and only in 2010–2024:** +2–3% (Jose and the extra 2010s pieces). In 1996–2010 the difference is under 0.1%.

**Reading.** The grouping hides real storms from the storm *record*, but on the adopted dunes it does not matter to the *model*. That is the same pattern as the trim-length check: storm peaks against dune heights decide overwash, and the extra hours only add volume. split12 is the more faithful series (every storm present, nothing the adopted series counts is lost). Switching to it would be a correctness fix with no measurable effect on the scores; keeping trim24 costs nothing measurable either. The decision is Hannah's.

**Incident.** The first launch of the four non-control runs completed the simulation but failed at the runner's metadata step. The runner stores the storm file relative to `data/hatteras_init` (`STORM_FILE.relative_to(HATTERAS_DATA_BASE)`, HAT_hindcast_1984_2024.py), and these files live here. Those runs lost their shoreline matrix and metadata, and were deleted and re-run with a `..`-relative path. `logs/launches_attempt1_failed-metadata.jsonl` is the record. A future experiment that swaps in a storm file from outside `data/` needs the same.

**Status: current.** Hannah adopted split12 on 2026-09-29 ("use split12 as the storm series going forward"): the builder has `--split-gap`, the files are `hindcast_storms/<window>/<window>_storms_v3_split12_trim24`, and `hat_env_forcings.DEFAULT_STORM_VARIANT` points at them. The matrix runs and the solved ends are still from `v3_trim24`.
