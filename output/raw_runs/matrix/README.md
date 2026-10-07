# matrix - the production runs

Only current windows live here. Superseded windows go to `D:\CASCADE_offload\output\raw_runs\matrix\` (see the bottom).

## Current windows (the 2026-10-05 calibration/test plan)

| window | role | presets | target |
|---|---|---|---|
| `1996_2009/` | **calibration**: 1996 ALACE DEM to the 2009 USACE DEM | zeroBE, edgeBE (nogroin), domainBE (blocking groin) | CoastSat net change 1996 start line to 2009 line, LOWESS-7 |
| `2009_2025/` | **test**: 2009 DEM to 2025-08-17, fills on | zeroBE, edgeBE (nogroin), domainBE (blocking groin, BE set 1) | CoastSat net change to the 2025-02-17..2026-02-17 mean |
| `1996_2025/` | **full window**: one 30-yr run, through Dec 2025 | zeroBE, edgeBE (nogroin) | CoastSat LRR 1996-2025 and net change to the calendar-2025 mean; figures in `output/comparisons/full_window/` |

The plan, scores and open decisions are in `scripts/hatteras_ms/DEM_TO_DEM_CALIBRATION.md`. Each run draws its own net-change figure under `<run>/figures/`.

## Moved off C: on 2026-10-07

`1996_2010/`, `2010_2024/` (the 14-yr windows, superseded 2026-10-03), `1996_2015/`, `2010_2026/` (the candidate windows, superseded by the calibration/test plan 2026-10-05) and `figures/` (management-ladder figures on 1996_2015 / 2010_2026). All of them are on `D:\CASCADE_offload\output\raw_runs\matrix\`, same paths, and their tracked text files are in git history before the commit that removed them.
