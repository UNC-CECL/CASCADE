# groin — the Buxton groin field

Everything about the groin at Buxton (GIS 5|6): what was observed, what the
module does on its own, how it was fitted to Buxton, and later, what it implies
for the future. The folders are grouped by **what kind of claim** they hold and
numbered in the order the work builds on itself. Inside 2- and 3-, steps run
from simplest to most complex.

```
1-observations/            true whatever the model does; nothing here is fitted
    structure_history.md       install 1969, last repair 1995, damage 2003; the fillet in numbers
    gis_data/                  the groins, wet/dry lines, 100 m transects, datum line, domains
    wetdry_photo_positions/    the 24-survey wet/dry and dune-line change tables (GIS 2-12, from 1967)
    coastsat_shoreline/        CoastSat around the field: era rates, profiles, GIFs
    coastsat_groin_condition/  when the gap stopped widening (1995, not 2004)
    figures/                   the two shorelines and the gap between them
2-module-tests/            what the module does on its own; no Buxton fit
    TEST_PLAN.md               the 2026-09-11 design for the idealized rig
    1-straight-coast/          four module versions on straight coasts, -40 to +40 deg
    2-solver-audit/            the dipole in BRIE's alongshore solve alone
    3-real-planform/           the dipole and blocking groin on the real planform, option A waves
    figures/                   the dipole's arithmetic as a schematic
3-hindcast/                the module fitted to Buxton in real model runs
    1-dipole-1967-2017/        SUPERSEDED: the dipole (M, f) fit on the 1967-2017 rig
    2-blocking-1996-2025/      the blocking groin on the calibration and test periods
4-forecast/                not started
```

**Where things stand (2026-10-08).**
- **Pinned in the runner:** the blocking groin, b 0.6, f 0.6, instant failure at the 2004 step (`output/calibration/groin/joint_fit.json`, read by `HAT_run_all.py`; evidence in `3-hindcast/2-blocking-1996-2025/2026-10-05-blocking-fit-calibration/`).
- **Open decision:** the failure schedule. CoastSat puts the end of trapping at the 1995 repair, and refitting on it gives a weaker groin failing from 1996 with the same post-failure strength (b × f 0.36), so test and forward runs are unchanged (`3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/`).
- **Open problem:** neither module is a physically correct groin yet. On a straight coast the pinned groin conserves sand but blocks only the tilt-driven part of the transport, not the net drift; a drift-blocking groin traps correctly but BRIE's alongshore solve does not conserve sand around it (`2-module-tests/1-straight-coast/`).

**Model output is not here.** Runs are under `output/raw_runs/experiments/groin/<study>/`; the pinned `joint_fit.json` is in `output/calibration/groin/`. The dipole (M, f) sweeps, `SELECTED_M60_f0.60/` and the 1967-2018 rig runs moved to `D:\CASCADE_offload\output\calibration\` on 2026-10-07.

## To do: move to the offload drive (D: attached 2026-10-09)

D: was not connected on 2026-10-08. When it is, move these to `D:\CASCADE_offload\hard-structures\groin\<same path>`, check that the file count and size match, delete them here, and give each parent README a line saying where they went (the 2026-10-07 `output/` offload is the model to follow):

| path | size | what it is |
|---|---|---|
| `1-observations/coastsat_shoreline/shoreline_output_coastsat/gif_frames_groin_area/` | 182 MB, 1,248 PNGs | frames behind `groin_analysis_shoreline_evolution.gif`; the analysis script rewrites them |
| `1-observations/coastsat_shoreline/shoreline_output_grid100m/gif_frames_groin_area/` | 66 MB, 390 PNGs | the same for the 100 m grid run |
| `1-observations/coastsat_shoreline/shoreline_output_grid100m/gif_frames_groin_area_zoomed/` | 63 MB, 501 PNGs | the zoomed GIF's frames |
| `3-hindcast/1-dipole-1967-2017/results/sensitivity_sweep/archive_july_20260824_081742/` | 31 KB, 24 files | an earlier dipole rig sweep |
| `3-hindcast/1-dipole-1967-2017/results/sensitivity_sweep/archive_pre1984start_20260830/` | 582 KB, 42 files | the dipole rig sweep before the 1984 start |
| `3-hindcast/1-dipole-1967-2017/runs/1967_1997_run/` | 1.7 MB, 6 PNGs | figures left from the precursor runs deleted on 2026-10-01 |

All untracked except one results CSV in each sweep archive; `git rm` those two. The GIFs themselves stay.

## Conventions

- **Code sits beside its data.** Each experiment folder holds its script, tables, figures and README, unlike the project's split into `scripts/`, `data/` and `output/`. This is deliberate, so an experiment can be read and rerun in one place.
- **Scripts follow `scripts/STYLE.md`** (header and author block, one-line comments, CONFIG block, reasoning in the folder README), since 2026-10-01.
- **Figures** follow the house style; each PNG has its PDF in `supporting/` or beside it, and its caption in `CAPTIONS.md`.
- **Experiments are dated** (`YYYY-MM-DD-<question>/`) inside a numbered step, and never edit the main code: they drive the unchanged runner with an in-process patch.

## Layout history

Reorganized on 2026-10-08 from the study's original folders. Old name → new place:

| was | now |
|---|---|
| `GROIN_PLAN.md` | split: `1-observations/structure_history.md` and `3-hindcast/1-dipole-1967-2017/dipole_fit_notes.md` |
| `HAT-groin-gis-analysis/` | `1-observations/coastsat_shoreline/`, its `gis_data/` to `1-observations/gis_data/` |
| `HAT-groin-condition-analysis/` | `1-observations/coastsat_groin_condition/` |
| `HAT-groin-buxton-output/shoreline_position_output/` and `HAT-groin-buxton-input/input_prep/shoreline_position/` | `1-observations/wetdry_photo_positions/` |
| `HAT-groin-figures/` | one `figures/` per category, each beside what it shows |
| `groin-module-test/0-solver-audit/` | `2-module-tests/2-solver-audit/`, `3-real-planform/`, `1-straight-coast/` |
| `groin-module-test/1-dem-to-dem/` | `3-hindcast/2-blocking-1996-2025/` |
| `HAT-groin-buxton-input/`, `HAT-groin-buxton-output/`, `HAT-buxton-hindcast-groin-test/` | `3-hindcast/1-dipole-1967-2017/inputs/`, `runs/`, `results/` |

Code paths were rewritten in the same change. Run metadata written before it still names the old folders.
