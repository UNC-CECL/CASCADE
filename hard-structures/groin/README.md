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
    wetdry_photo_positions/    the 24-survey wet/dry and dune-line change tables (GIS 2-12, from 1967), and
                               the groin in one picture: the two shorelines and the gap between them
    shoreline_rates_by_era/    wide view, ±60 km: rates by era and decade, how far the effect reaches
    gap_across_groins/         narrow view: the gap across the groins each year, when it stopped widening
2-module-tests/            what the module does on its own; no Buxton fit
    TEST_PLAN.md               the 2026-09-11 design for the idealized rig
    1-straight-coast/          six module versions on straight coasts, -40 to +40 deg
    2-solver-audit/            the dipole in BRIE's alongshore solve alone
    3-real-planform/           the dipole and blocking groin on the real planform, option A waves
    4-groin-field/             four groins on one straight coast: which module behaves like a groin field
    5-fine-grid-reference/     a real-sized field on a 25 m grid as the yardstick for the modules
    6-long-term/               the modules over 200 years: stable, conserving, levelling off?
    figures/                   the dipole's arithmetic as a schematic
3-hindcast/                the module fitted to Buxton in real model runs
    1-dipole-1967-2017/        SUPERSEDED: the dipole (M, f) fit on the 1967-2017 rig
    2-blocking-1996-2025/      the blocking groin on the calibration and test periods
4-forecast/                not started
```

**Where we left off (2026-10-09).** Commit `6c8cde10`.

- **Pinned in the runner, unchanged:** the blocking groin, b 0.6, f 0.6, instant failure at the 2004 step (`output/calibration/groin/joint_fit.json`, read by `HAT_run_all.py`).
- **Module design, settled in principle:** build on BRIE and Barrier3D, never edit them. Of the six versions tested on straight coasts (`2-module-tests/1-straight-coast/`), the drift-blocking groin ("option 1") is the one to carry forward. It moves a fraction b of the wave-driven drift arriving at the groin across GIS 5|6 each year, read from BRIE's own wave climate. Its known limitation is BRIE's sand leak at the bends the groin makes, about +100 m of summed shoreline over 1996-2009; it is reported, not corrected. The corrected version (module 6) conserves sand but changes BRIE's behaviour at Buxton's natural bend (`2-module-tests/3-real-planform/corrected_groin/`).
- **Option 1 in the full model** (`3-hindcast/2-blocking-1996-2025/2026-10-08-drift-groin/`): best combined fit on 1996-2009 at b 0.1, f 0.4 (joint RMSE 15.8 m against the pin's 20.9 m). The annual CoastSat score and the photo-date score want different b (0.06-0.07 against 0.12). On 2009-2025 it beats the pin but not no groin: it holds the gap too wide (−18 m against CoastSat's −50 m).
- **Next, in order:**
  1. **Decide the failure schedule** (options A-C in `3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/`). It now limits both groin versions on 2009-2025.
  2. **Decide which observation drives the calibration:** annual CoastSat, the photo dates, or the joint score.
  3. **Refit option 1 under that schedule and score,** then rerun the 2009-2025 test.
  4. **Then decide whether option 1 replaces the pin.** If it does, move the drift groin from `drift_groin.py`'s in-process patch into `cascade/groin.py` as its own callback, alongside `BlockingGroinCallback`, and refit for the runner.
- **Check before quoting:** the pinned 2009-2025 test run reads −416 m of island change against no groin. It is a 10-05 run and may not share the baseline's setup.
- **Still open:** the D: offload below.

**Model output is not here.** Runs are under `output/raw_runs/experiments/groin/<study>/`; the pinned `joint_fit.json` is in `output/calibration/groin/`. The dipole (M, f) sweeps, `SELECTED_M60_f0.60/` and the 1967-2018 rig runs moved to `D:\archive\CASCADE_offload\output\calibration\` on 2026-10-07.

## To do: move to the offload drive (D: attached 2026-10-09)

D: was not connected on 2026-10-08. When it is, move these to `D:\archive\CASCADE_offload\hard-structures\groin\<same path>`, check that the file count and size match, delete them here, and give each parent README a line saying where they went (the 2026-10-07 `output/` offload is the model to follow):

| path | size | what it is |
|---|---|---|
| `3-hindcast/1-dipole-1967-2017/results/sensitivity_sweep/archive_july_20260824_081742/` | 31 KB, 24 files | an earlier dipole rig sweep |
| `3-hindcast/1-dipole-1967-2017/results/sensitivity_sweep/archive_pre1984start_20260830/` | 582 KB, 42 files | the dipole rig sweep before the 1984 start |
| `3-hindcast/1-dipole-1967-2017/runs/1967_1997_run/` | 1.7 MB, 6 PNGs | figures left from the precursor runs deleted on 2026-10-01 |

All untracked except one results CSV in each sweep archive; `git rm` those two. (The shoreline GIF frames that were on this list were removed on 2026-10-08, when the plan-view GIFs replaced them.)

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
| `HAT-groin-gis-analysis/` | `1-observations/shoreline_rates_by_era/` (first `coastsat_shoreline/`), its `gis_data/` to `1-observations/gis_data/` |
| `HAT-groin-condition-analysis/` | `1-observations/gap_across_groins/` (first `coastsat_groin_condition/`) |
| `HAT-groin-buxton-output/shoreline_position_output/` and `HAT-groin-buxton-input/input_prep/shoreline_position/` | `1-observations/wetdry_photo_positions/` |
| `HAT-groin-figures/` | `figures/` beside what each shows: the two-shorelines figure in `1-observations/wetdry_photo_positions/figures/` (it had its own `1-observations/figures/` until 2026-10-09), the module schematic in `2-module-tests/figures/`, the timeline in `3-hindcast/1-dipole-1967-2017/figures/` |
| `groin-module-test/0-solver-audit/` | `2-module-tests/2-solver-audit/`, `3-real-planform/`, `1-straight-coast/` |
| `groin-module-test/1-dem-to-dem/` | `3-hindcast/2-blocking-1996-2025/` |
| `HAT-groin-buxton-input/`, `HAT-groin-buxton-output/`, `HAT-buxton-hindcast-groin-test/` | `3-hindcast/1-dipole-1967-2017/inputs/`, `runs/`, `results/` |

Code paths were rewritten in the same change. Run metadata written before it still names the old folders.
