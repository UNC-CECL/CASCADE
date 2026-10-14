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

**Where we left off (2026-10-14).** Numbers on the calendar-year setup (commit `e9c5e663`).

- **In the runner since 2026-10-14: the corrected groin**, `cascade.groin.CorrectedGroinCallback`, the only groin left. Set in `scripts/hatteras_ms/hat_run.yaml` by the share of the arriving drift held back, per window since the afternoon of 2026-10-14 (`groin.strength_1996` 0.16, `groin.strength_2009` 0.10; `groin.strength` overrides both), and `groin.end_year` 2017. At its defaults the runner reproduces the refit runs to 0 m in both windows. The dipole (M, f) and blocking (b, f) groins, the M/f sweeps (`scripts/hatteras_ms/groin-sweep/`) and the driver's sweep and joint-fit stages were removed the same day; the last commit with them is `4d937b79`.
- **Module design:** build on BRIE and Barrier3D, never edit them. Of the seven versions tested on straight coasts, the corrected groin (version 6: the drift groin plus a correction for BRIE's sand leak at the groin's own bends) is the one carried forward; versions 2-5 and 7 are in `2-module-tests/archive/`. Its reference is the groin's own tracked footprint, so it needs no no-groin run (`3-hindcast/2-blocking-1996-2025/2026-10-09-corrected-groin/`, findings 15-17).
- **The groin's strength, from one continuous run 1967-2025** (`3-hindcast/2-blocking-1996-2025/2026-10-13-pregroin-planform/`): 0.25 from installation to 2009 (perhaps 0.30 after the 1995 repair), 0.20-0.25 to 2016 (0.25 within 0.5 m, so one strength is kept), off from the 2017 fill. It follows the observed gap to 17 m over the whole record, and 0.25 is the minimum under three different pre-1984 storm sequences.
- **The matrix runs use the window strengths, not 0.25.** Each hindcast window starts from its observed shoreline, which already holds the fillet. From there 0.25 overshoots (1996-2009 30.4 m against 13.4 m at 0.16; 2009-2025 28.1 m, worse than no groin's 21.9 m, against 9.7 m at 0.10). Refit on today's runner in `3-hindcast/2-blocking-1996-2025/2026-10-14-window-refit/`, matching the 10-12 fit (`2026-10-09-corrected-groin/`, findings 18-22) to 0.1 m. The window strengths measure how fast the model relaxes the observed fillet, not the groin itself: run as one history from 1967 they leave the gap 30-40 m below the record. 0.25 stays the groin's whole-life strength, for runs that start before the fillet exists. The 2009 value is fitted on the test window, so for the groin that window is not an independent test.
- **The jump at GIS 5|6 in no-groin runs is real planform**, not a source or sink: about 100 m of the start offsets' 240-264 m step there is the natural bend toward the cape and about 150 m the fillet the groins built (the 1967 dune line steps 104 m). Without a groin, BRIE smooths it.
- **Next, in order:**
  1. **Rerun the groin matrix** (`HAT_run_all.py --stages 3`, started 2026-10-14 at 0.16 / 0.10) and point BE set 1 (`HAT_be_domain_solve_net_change.py`, whose start run is still a blocking-groin run) at a corrected-groin start.
  2. **Forecast** (`4-forecast/`, `groin.strength` set): the 2026 rebuild at about 0.25; at full strength from day one the model builds the fillet in about 3 years against about 9 observed, so ramp it in over 5-9 years if the timing matters. Only the southernmost of the three groins is rebuilt, so 0.25 is an upper bound.
  3. Open: the model's fillet sits in the two domains beside the groin where the observed one spans three or four; the 2017 fill starts at GIS 6 and makes its own 35 m gap in the no-groin run.

**Model output is not here.** Runs are under `output/raw_runs/experiments/groin/<study>/`; the old blocking pin's `joint_fit.json` is in `output/calibration/groin/` (no longer read by anything). The dipole (M, f) sweeps, `SELECTED_M60_f0.60/` and the 1967-2018 rig runs moved to `D:\archive\CASCADE_offload\output\calibration\` on 2026-10-07.

## Conventions

- **Code sits beside its data.** Each experiment folder holds its script, tables, figures and README, unlike the project's split into `scripts/`, `data/` and `output/`. This is deliberate, so an experiment can be read and rerun in one place.
- **Scripts follow `scripts/STYLE.md`** (header and author block, one-line comments, CONFIG block, reasoning in the folder README), since 2026-10-01.
- **Figures** follow the house style; each PNG has its PDF in `supporting/` or beside it, and its caption in `CAPTIONS.md`.
- **Experiments are dated** (`YYYY-MM-DD-<question>/`) inside a numbered step, and never edit the main code: they drive the unchanged runner with an in-process patch.
- **Code these studies ran against.** The dipole (`GroinCallback`) and blocking (`BlockingGroinCallback`) groins that the scripts here use or patch were removed from `cascade/groin.py` on 2026-10-14, when the corrected groin became the runner's only groin. The scripts ran against the code at commit `4d937b79`; check that commit out to rerun one.

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
