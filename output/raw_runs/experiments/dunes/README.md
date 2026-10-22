# dunes

The dune ceiling (Dmaxel), the growth rate (rmin, rmax) and the design height, and how the modelled dunes compare with the DEMs and the lidar surveys. Split out of `storms-overwash-and-dunes/` on 2026-10-20 (`../RENAMES.md`); the storm and overwash scoring studies are in [`../storms-overwash/`](../storms-overwash/README.md). The 10-13 studies are written up in `scripts/hatteras_ms/experiments/README.md` under their scripts' names; the 10-19 and 10-20 ones in their own `NOTE.md`.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| `2026-10-13-dune-design-height` | Does raising the design height from 3.0 to 3.94 m MHW change the calibration run? | Skill 0.548 → 0.575, fewer false alarms (121 → 113), more rebuilds (32 → 44). | **current**; 3.94 m MHW in the runner and notebook |
| `2026-10-13-dune-ceiling-uniform` | Does one island-wide ceiling do as well as the per-cell ceilings, and at what height? | Raises calibration skill only by cutting false alarms; a cliff between 5.5 and 6.0 m; halves the test hit rate. | record |
| `2026-10-13-dune-ceiling-domain-value` | Does one ceiling per domain, from its own starting crests, beat the per-cell ceilings? | Keeps the alongshore pattern but loses hits where the observed overwash went through. Not adopted. | record |
| `2026-10-13-dune-ceiling-callback` | Can a CASCADE callback reproduce the per-cell ceilings on base Barrier3D? | Close on scores (calibration skill 0.553 vs 0.575), but 0.4-0.6 m crest differences (P95) and more test rebuilds. | record |
| `2026-10-13-dune-ceiling-range` | Across island-wide ceilings 4.5-6.0 m, how much overwash, and where? | One Dmaxel of 5.5 m NAVD88 adopted 10-13 (skill 0.594 / 0.120). | superseded 2026-10-20 (Dmaxel 12.2 m) |
| `2026-10-13-dune-growth-rate` | How much do the growth-rate bounds change the dunes and the overwash? | Little on skill or shoreline; 0.55-0.95 kept. | superseded 2026-10-20 (r 0.005-0.01) |
| `2026-10-19-dune-crest-dem-vs-model` | Is the ceiling an elevation or a height, and does the model build the dunes up early? | Elevation (Dmax = Dmaxel − 1.7 m berm). Low dunes reach the ceiling by year 3; cells under the run-up fall from 22% to 7% (2012). | record |
| `2026-10-19-dune-ceiling-dem-growth` | Ceiling at the measured P99 crest with slower growth? | Skill falls (0.59 → 0.24-0.44) and the 2009 crest sits 1.6-3.0 m above the lidar. | record |
| `2026-10-19-dune-ceiling-growth-grid` | Ceiling × growth grid, picked on 1996-2009 skill then the 2009 lidar. | Winner 5.0 m, r 0.25-0.45 (skill 0.612). | record; not adopted |
| `2026-10-19-dune-growth-low-rates` | Low growth rates, ranked on the surveyed dune trajectory. | No growth fits the trajectory best (P10 RMSE 0.76 vs 2.20 m for the control), at a cost of skill 0.51 vs 0.59. | record |
| `2026-10-19-dune-crest-vs-surveys` | Modelled dunes against the 2014-2019 surveys. | The real low dunes stayed low; the control fills them in 3-4 years (alongshore r 0.18 by 2019); the grid winner is closer (0.59-0.76). | record |
| `2026-10-20-dune-findings-figures` | Three summary figures of the 10-13 and 10-19 findings. | Figures in `figures/`, numbers in `tables/`. | record |
| `2026-10-20-stable-dune-growth-rate` | With Dmax 10.5 m above the berm, what growth rate keeps stable dunes stable 1996-2009 (advisor's method)? | Mean r 0.005 (0-0.01) closest on the seven stable domains (−0.11 m). | superseded 2026-10-20 by the `-topo-v3` re-run on the corrected 1996 start; its NOTE.md still sets out the selection |
| `2026-10-20-stable-dune-growth-rate-topo-v3` | The same on the 1996 start with the ALACE offset corrected. | r 0.02 fits (+0.01 m); 0.01-0.03 adopted. | **current**, but runs only: the write-up, tables and figures moved to `data/hatteras_init/1-barrier3d-domains/1996-start/dune-growth-rate/` on 2026-10-21 |
| `2026-10-20-test-period-growth-rate` | Does 2009-2025 want a faster growth rate? | 0.02-0.04 tracks the 2014-2019 median crest within 0.16 m; adopted for 2009. | **current**, but runs only: the write-up and table moved to `data/hatteras_init/1-barrier3d-domains/2009-start/dune-growth-rate/` on 2026-10-21 |

**Status**: **current** means its answer is in use now. **record** means a finished check, kept so the number can be traced. **superseded** means a later study replaced the answer.

Every 10-13 and 10-19 study ran on the old storms (`v3_split12_trim24`), beach slope 0.04, and dune settings that changed on 10-20. Their `.npz` model states were deleted on 2026-10-20 (`../../README.md`); the tables, logs and metadata are kept. Studies from before 2026-10-05 are on `D:\hatteras_ch1\model_runs\output\raw_runs\experiments\storms-overwash-and-dunes\`.

Back to [the map](../README.md).
