# groin — the Buxton groin field

Everything about the groin at Buxton (GIS 5|6): what was observed, what the
module does on its own, how it was fitted to Buxton, and later, what it implies
for the future. The folders are grouped by **what kind of claim** they hold and
numbered in the order the work builds on itself. Inside 2- and 3-, steps run
from simplest to most complex.

```
1-observations/            true whatever the model does; nothing here is fitted
    data/                      what was recorded or digitized
        structure_history.md       install 1969, last repair 1995, damage 2003, decay 2015-21; sources
        gis/                       the groins, wet/dry lines, 100 m transects, datum line, domains
    analysis/                  what was measured from it, each with its scripts, tables and figures
        1-wetdry_photo_positions/  the 24-survey wet/dry and dune-line change tables (GIS 2-12, from 1967)
        2-gap_across_groins/       narrow view: the gap each year, when it stopped widening, the state in 1996
        3-shoreline_rates_by_era/  wide view, ±60 km: rates by era and decade, how far the effect reaches
2-module-tests/            what the module does on its own; no Buxton fit
    TEST_PLAN.md               the 2026-09-11 design for the idealized rig
    1-straight-coast/          version 1 (dipole) vs version 6 (corrected) on straight coasts, -40 to +40 deg
    2-solver-audit/            the dipole in BRIE's alongshore solve alone
    3-real-planform/           the dipole and blocking groin on the real planform, option A waves
    4-groin-field/             four groins on one straight coast: which module behaves like a groin field
    5-fine-grid-reference/     a real-sized field on a 25 m grid as the yardstick for the modules
    6-long-term/               the modules over 200 years: stable, conserving, levelling off?
    figures/                   the dipole's arithmetic as a schematic
    archive/                   the retired versions 2-5 and 7 (version 2 is still the runner's groin)
3-hindcast/                the module fitted to Buxton in real model runs
    1-dipole-1967-2017/        SUPERSEDED: the dipole (M, f) fit on the 1967-2017 rig
    2-blocking-1996-2025/      the blocking groin on the calibration and test periods
4-forecast/                not started
```

**Where we left off (2026-10-12).** All numbers on the 10-12 setup (commit `6219a55b`: test period 2009-2025, CoastSat's gap closing by 62.0 m).

- **Pinned in the runner, unchanged:** the blocking groin, b 0.6, f 0.6, instant failure at the 2004 step (`output/calibration/groin/joint_fit.json`, read by `HAT_run_all.py`).
- **Module design, settled:** build on BRIE and Barrier3D, never edit them. Of the seven versions tested on straight coasts, the corrected groin (version 6: the drift groin plus a correction for BRIE's sand leak at the groin's own bends) is the one to carry forward; versions 2-5 and 7 are in `2-module-tests/archive/`. It is about 26 lines of physics plus about 15 for the reference below.
- **No no-groin run needed** (`3-hindcast/2-blocking-1996-2025/2026-10-09-corrected-groin/`, findings 15-17): the correction's reference is now the groin's own tracked footprint, not a matching no-groin run. It reproduces no groin exactly at b 0 and stays within 5 m of the old method; it is the default.
- **Best groin so far: the two-number groin**, one constant strength s 0.12 (fitted to annual CoastSat on 1996-2009) and an end year of 2017 (the Buxton fill). Test joint RMSE 14.6 m, against no groin 29.2, the pin 31.6, the drift groin 28.1 and the three-number corrected groin 16.1; best on both halves, including the independent 2009-2016 (10.3 m). In calibration the corrected groin edges out the pin (20.3-20.4 m against 22.9 m).
- **What the decisions now look like:**
  1. **Failure schedule:** the calibration cannot tell the schedules apart; what the test period needs is the groin off at the 2017 fill. The decline the record describes (2015-2021) fits worse (finding 10). The 2017 end was suggested by the test misfit, so only 2009-2016 is fully independent.
  2. **Calibration target:** annual CoastSat generalises (s 0.12); the photo dates do not (s 0.25, pulled by the 2004 survey where photos and CoastSat disagree by 42 m). The joint minimum is 0.18.
- **Next, in order:**
  1. **Decide** the strength's target (annual CoastSat, recommended by findings 12-13) and the 2017 end.
  2. **Decide whether the corrected groin replaces the pin.** If it does, move it from `corrected_groin.py`'s in-process patch into `cascade/groin.py` as its own callback, alongside `BlockingGroinCallback`, and point the runner at it.

**Model output is not here.** Runs are under `output/raw_runs/experiments/groin/<study>/`; the pinned `joint_fit.json` is in `output/calibration/groin/`. The dipole (M, f) sweeps, `SELECTED_M60_f0.60/` and the 1967-2018 rig runs moved to `D:\archive\CASCADE_offload\output\calibration\` on 2026-10-07.

## Conventions

- **Code sits beside its data.** Each experiment folder holds its script, tables, figures and README, unlike the project's split into `scripts/`, `data/` and `output/`. This is deliberate, so an experiment can be read and rerun in one place.
- **Scripts follow `scripts/STYLE.md`** (header and author block, one-line comments, CONFIG block, reasoning in the folder README), since 2026-10-01.
- **Figures** follow the house style; each PNG has its PDF in `supporting/` or beside it, and its caption in `CAPTIONS.md`.
- **Experiments are dated** (`YYYY-MM-DD-<question>/`) inside a numbered step, and never edit the main code: they drive the unchanged runner with an in-process patch.

## Layout history

Reorganized on 2026-10-08 from the study's original folders. Old name → new place:

| was | now |
|---|---|
| `GROIN_PLAN.md` | split: `1-observations/data/structure_history.md` and `3-hindcast/1-dipole-1967-2017/dipole_fit_notes.md` |
| `HAT-groin-gis-analysis/` | `1-observations/analysis/3-shoreline_rates_by_era/` (first `coastsat_shoreline/`), its `gis_data/` to `1-observations/data/gis/` |
| `HAT-groin-condition-analysis/` | `1-observations/analysis/2-gap_across_groins/` (first `coastsat_groin_condition/`) |
| `HAT-groin-buxton-output/shoreline_position_output/` and `HAT-groin-buxton-input/input_prep/shoreline_position/` | `1-observations/analysis/1-wetdry_photo_positions/` |
| `HAT-groin-figures/` | `figures/` beside what each shows: the two-shorelines figure in `1-observations/analysis/1-wetdry_photo_positions/figures/` (it had its own `1-observations/figures/` until 2026-10-09), the module schematic in `2-module-tests/figures/`, the timeline in `3-hindcast/1-dipole-1967-2017/figures/` |
| `groin-module-test/0-solver-audit/` | `2-module-tests/2-solver-audit/`, `3-real-planform/`, `1-straight-coast/` |
| `groin-module-test/1-dem-to-dem/` | `3-hindcast/2-blocking-1996-2025/` |
| `HAT-groin-buxton-input/`, `HAT-groin-buxton-output/`, `HAT-buxton-hindcast-groin-test/` | `3-hindcast/1-dipole-1967-2017/inputs/`, `runs/`, `results/` |

Code paths were rewritten in the same change. Run metadata written before it still names the old folders.
