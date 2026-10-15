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
        1-aerial_photo_gap/  the 24-survey wet/dry and dune-line change tables (GIS 2-12, from 1967)
        2-coastsat_gap/       narrow view: the gap each year, when it stopped widening, the state in 1996
        4-shoreline_rates_by_era/  wide view, ±60 km: rates by era and decade, how far the effect reaches
        5-shoreline_maps/          the shoreline each year as true-scale maps: groin window, whole coast
2-module-tests/            what the module does on its own; no Buxton fit
    TEST_PLAN.md               the 2026-09-11 design for the idealized rig
    1-straight-coast/          version 1 (dipole) vs version 6 (corrected) on straight coasts, -40 to +40 deg
    2-solver-audit/            the dipole in BRIE's alongshore solve alone
    3-real-planform/           the dipole and blocking groin on the real planform, option A waves
    4-groin-field/             four groins on one straight coast: which module behaves like a groin field
    5-fine-grid-reference/     a real-sized field on a 25 m grid as the yardstick for the modules
    6-long-term/               the modules over 200 years: stable, conserving, levelling off?
    7-stress/                  the runner's groin pushed until it breaks; where its bounds come from
    8-mirror/                  drift from either direction: the groin mirrors exactly, BRIE nearly
    figures/                   the dipole's arithmetic as a schematic
    archive/                   the retired versions 2-5 and 7 (version 2 is still the runner's groin)
3-hindcast/                the module fitted to Buxton in real model runs
    1-dipole-1967-2017/        SUPERSEDED: the dipole (M, f) fit on the 1967-2017 rig
    2-blocking-1996-2025/      the blocking groin on the calibration and test periods
4-forecast/                not started
```

**Where we left off (2026-10-14; bounds added 2026-10-15).** Numbers on the calendar-year setup (commit `e9c5e663`).

- **Decision (2026-10-15): two separate window fits.** The groin keeps its own strength in each hindcast window, 0.16 for 1996-2009 and 0.10 for 2009-2025, off from the 2017 fill, because it kept deteriorating (2003 damage, 2015-21 decay). The 2009 value is fitted on that window, so for the groin 2009-2025 is a second fit, not an independent test; the split-sample result below is the groin's validation. A strength carried from the calibration window into the test (before/after the 2003 damage) was considered and left for later.
- **Validated by split sample (2026-10-15).** Calibrated on 1996-2009 (0.16) and run unchanged on 2009-2025, the groin cuts the gap error from 21.9 m (no groin) to 14.1 m, skill 0.59 against no groin; the reverse transfer has skill 0.84. Only the groin's strength is held out (`3-hindcast/2-blocking-1996-2025/2026-10-15-split-sample/`). A mirror check shows the groin handles drift from either direction identically (`2-module-tests/8-mirror/`).
- **Since 2026-10-15 the module has stated bounds and stops a run outside them** (`GroinOutsideTestedRange`; see "The module's bounds" below). Both Buxton windows run unaffected and bit-identical. The groin-on wave sweep's asymmetry 0.7 and 0.8 cells (`output/raw_runs/sensitivity/waveasym/`, run 2026-10-14/15) are outside the bounds: they show the groin breaking, not the coast's sensitivity to asymmetry.
- **In the runner since 2026-10-14: the corrected groin**, `cascade.groin.CorrectedGroinCallback`, the only groin left. Set in `scripts/hatteras_ms/hat_run.yaml` by the share of the arriving drift held back, per window since the afternoon of 2026-10-14 (`groin.strength_1996` 0.16, `groin.strength_2009` 0.10; `groin.strength` overrides both), and `groin.end_year` 2017. At its defaults the runner reproduces the refit runs to 0 m in both windows. The dipole (M, f) and blocking (b, f) groins, the M/f sweeps (`scripts/hatteras_ms/groin-sweep/`) and the driver's sweep and joint-fit stages were removed the same day; the last commit with them is `4d937b79`.
- **Module design:** build on BRIE and Barrier3D, never edit them. Of the seven versions tested on straight coasts, the corrected groin (version 6: the drift groin plus a correction for BRIE's sand leak at the groin's own bends) is the one carried forward; versions 2-5 and 7 are in `2-module-tests/archive/`. Its reference is the groin's own tracked footprint, so it needs no no-groin run (`3-hindcast/2-blocking-1996-2025/2026-10-09-corrected-groin/`, findings 15-17).
- **The groin's strength, from one continuous run 1967-2025** (`3-hindcast/2-blocking-1996-2025/2026-10-13-pregroin-planform/`): 0.25 from installation to 2009 (perhaps 0.30 after the 1995 repair), 0.20-0.25 to 2016 (0.25 within 0.5 m, so one strength is kept), off from the 2017 fill. It follows the observed gap to 17 m over the whole record, and 0.25 is the minimum under three different pre-1984 storm sequences.
- **The matrix runs use the window strengths, not 0.25.** Each hindcast window starts from its observed shoreline, which already holds the fillet. From there 0.25 overshoots (1996-2009 30.4 m against 13.4 m at 0.16; 2009-2025 28.1 m, worse than no groin's 21.9 m, against 9.7 m at 0.10). Refit on today's runner in `3-hindcast/2-blocking-1996-2025/2026-10-14-window-refit/`, matching the 10-12 fit (`2026-10-09-corrected-groin/`, findings 18-22) to 0.1 m. The window strengths measure how fast the model relaxes the observed fillet, not the groin itself: run as one history from 1967 they leave the gap 30-40 m below the record. 0.25 stays the groin's whole-life strength, for runs that start before the fillet exists. The 2009 value is fitted on the test window, so for the groin that window is not an independent test.
- **The jump at GIS 5|6 in no-groin runs is real planform**, not a source or sink: about 100 m of the start offsets' 240-264 m step there is the natural bend toward the cape and about 150 m the fillet the groins built (the 1967 dune line steps 104 m). Without a groin, BRIE smooths it.
- **Next, in order:**
  1. **Point BE set 1** (`HAT_be_domain_solve_net_change.py`, whose start run is still a blocking-groin run) at a corrected-groin start. The groin matrix was rerun at 0.16 / 0.10 on 2026-10-14 (22 of 22 runs). Against the CoastSat net change (full management, edgeBE) the groin helps in 1996-2009 near the groin (GIS 2-12: 21.7 against 25.5 m without it; island-wide 19.0 against 19.6 m) and makes no difference in 2009-2025, where it is weak and off from 2017 (`2026-10-14-window-refit/matrix_groin_effect_figure.py`). The record's 1996-2009 erosion at GIS 7-14 (15-38 m) is in no run, and is not a groin-reach problem (`2026-10-14-fillet-shape/`, finding 5).
  2. **Forecast** (`4-forecast/`, `groin.strength` set): the 2026 rebuild at about 0.25; at full strength from day one the model builds the fillet in about 3 years against about 9 observed, so ramp it in over 5-9 years if the timing matters. Only the southernmost of the three groins is rebuilt, so 0.25 is an upper bound.
  3. Open: the model's fillet sits in the two domains beside the groin where the observed one spans three or four (the groins all sit in GIS 6 and the record's build-out peaks at GIS 7; a field layout across GIS 6 fits the 1967-1995 build-up better but not the hindcast windows, so the runner keeps its 5|6 face: `3-hindcast/2-blocking-1996-2025/2026-10-14-face-6-7/`); the 2017 fill starts at GIS 6 and makes its own 35 m gap in the no-groin run.

## The module's bounds

The groin is a simplified module: it holds back a share of BRIE's own net drift at one face and cancels the
sand BRIE's row-scaled solve would create or destroy at the bends it makes. It was tested to destruction on
2026-10-15 (`2-module-tests/7-stress/`, about 600 paired runs of the runner's own callback on a straight
coast, checked against the full model at Buxton). It was kept simple rather than extended: outside the range
below it stops the run with `cascade.groin.GroinOutsideTestedRange` instead of returning a wrong answer.

**What the testing found.**

- **At option A it works.** Strength 0 is exactly the no-groin run. Up to strength 0.5 it keeps the sand
  balance on every coast angle tested (−30° to +30°); at Buxton's 0.25 the largest leak is 34 m of summed
  shoreline in 100 years, and nothing runs away over 500 years. It is unaffected by a sharp starting step at
  Buxton's angles (±300 m; Buxton's is 256 m), a single switch-off, the cross-shore change Barrier3D passes in,
  and run length.
- **It breaks in two ways, both set by BRIE's alongshore diffusion at the groin** (its diffusion number
  r = D dt / 2 dy²):
  - *Diffusion too strong* (high-angle fraction 0.3 or less; Hs 2.5 and above on coasts turned with the
    drift): the leak correction, an explicit term in front of BRIE's implicit step, overshoots, and the groin
    loses from about 0.4 km to tens of km of summed shoreline in 100 years.
  - *Diffusion switched off* (asymmetry 0.7 and above, or any groin face steep enough that BRIE clips its
    diffusivity to zero): nothing carries sand back across the groin, the groin keeps holding drift, and the
    gap grows without limit, 2-3.4 km in 100 years. The module has no groin length, so nothing stops it.

  These are limits of BRIE's alongshore diffusion, which any groin built on BRIE inherits.
- **Year-to-year changes in the wave climate break it at any setting** (35 of 36 runs leak), because the
  correction assumes one climate.
- **Repeated switching on and off leaks** (up to 39 m at strength 0.25 and 458 m at 1.0, cycling every 5
  years).
- **The correction zone does not wrap the ends of the domains**, so a groin near either end leaks.

**The bounds.**

| | tested range | Buxton hindcast |
|---|---|---|
| wave climate | one climate for the whole run; one parameter at a time from option A, these ran without breaking in the full model: Hs 0.75-3.0, Tp 6-12 s, high-angle fraction 0.45-0.55, asymmetry 0.5-0.65 (high-angle 0.4 not run with the guard) | option A, fixed |
| groin face | at most 45° | 16-29° |
| strength | up to 0.5 | 0.16 (1996), 0.10 (2009), 0.25 whole life |
| switching | one switch-off | off from 2017 |
| placement | face at least `half_width - 1` (9) cells from either end of the domains | pads 19\|20 of 120 |

**What the module checks** (constants on `CorrectedGroinCallback`). Each year the groin is installed, the run
stops if the wave climate has changed, if the groin has made or lost more than 10 m of summed shoreline (or 5%
of its footprint's gross, whichever is larger: `SAND_TOLERANCE_M`, `SAND_TOLERANCE_FRACTION`), or if the groin
face is steeper than 45° (`MAX_FACE_ANGLE_DEG`). When it is built, it refuses a correction zone that runs past
the ends. On the stress test these stop 225 of the 227 runs that break, mostly by year 2. Both Buxton windows
run with the checks silent and identical to the matrix runs (0.0 m); asymmetry 0.8 is stopped in 2000 and
high-angle 0.3 in year 3 (`output/raw_runs/experiments/2026-10-15-groin-guard-check/`).

**Caveats for reporting.**

- **The thresholds are empirical.** 10 m / 5% and 45° were chosen from the same tests, and the full-model
  check is four runs. Describe them as an applicability range established by testing, not derived limits.
  45° sits between every well-behaved full-model run (31° at most; asymmetry 0.65 reached 41°) and the broken
  ones, which pass 45° within four years.
- **Option A sits near BRIE's own edge.** At high-angle fraction 0.5 BRIE's diffusivity is already zero on a
  coast at 0°; the groin works at Buxton because its faces sit at −16° to −29°, where diffusion is still on.
  The tested window in high-angle fraction is narrow (0.45-0.55).
- **The wave sensitivity loses its upper asymmetry range with the groin on.** Report asymmetry 0.7-0.8 as
  outside the module's range, or present that axis with the groin off.
- **Forecasts.** A forecast with a changing wave climate or a much stronger groin is outside the bounds. It
  would need a groin length (so the fillet stops at the tip and bypasses) and a correction that holds when
  the waves change.

**Validation and symmetry.** Split-sample: the strength fitted on 1996-2009 keeps skill 0.59 against no groin
on 2009-2025, and the reverse 0.84 (`3-hindcast/2-blocking-1996-2025/2026-10-15-split-sample/`). Mirror: the
groin's response to drift from the other direction is its exact mirror image; whole runs mirror only to within
BRIE's own directional bias, which the groin does not add to at Buxton-like settings by more than 10 m
(`2-module-tests/8-mirror/`).

**A sentence for the methods**, as a starting point: *"The groin module is applicable while BRIE's alongshore
diffusivity at the groin is non-zero and moderate: a fixed wave climate near the calibrated one, a groin face
under 45° and a strength up to 0.5. Within this range it conserves sand and reaches a steady fillet; outside it
the model stops the run."*

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
| `HAT-groin-gis-analysis/` | `1-observations/analysis/4-shoreline_rates_by_era/` (first `coastsat_shoreline/`), its `gis_data/` to `1-observations/data/gis/` |
| `HAT-groin-condition-analysis/` | `1-observations/analysis/2-coastsat_gap/` (first `coastsat_groin_condition/`) |
| `HAT-groin-buxton-output/shoreline_position_output/` and `HAT-groin-buxton-input/input_prep/shoreline_position/` | `1-observations/analysis/1-aerial_photo_gap/` |
| `HAT-groin-figures/` | `figures/` beside what each shows: the two-shorelines figure in `1-observations/analysis/1-aerial_photo_gap/figures/` (it had its own `1-observations/figures/` until 2026-10-09), the module schematic in `2-module-tests/figures/`, the timeline in `3-hindcast/1-dipole-1967-2017/figures/` |
| `groin-module-test/0-solver-audit/` | `2-module-tests/2-solver-audit/`, `3-real-planform/`, `1-straight-coast/` |
| `groin-module-test/1-dem-to-dem/` | `3-hindcast/2-blocking-1996-2025/` |
| `HAT-groin-buxton-input/`, `HAT-groin-buxton-output/`, `HAT-buxton-hindcast-groin-test/` | `3-hindcast/1-dipole-1967-2017/inputs/`, `runs/`, `results/` |

Code paths were rewritten in the same change. Run metadata written before it still names the old folders.
