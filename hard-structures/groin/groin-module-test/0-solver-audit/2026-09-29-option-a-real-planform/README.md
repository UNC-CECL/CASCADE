# Groin stability under option A, on the real planform (2026-09-29)

## The scripts

| script | question | writes, beside the script |
|---|---|---|
| `groin_stability_option_a.py` | is the dipole stable under option A, on the real planform? | `stability_summary.csv`, `stability_traces.csv` |
| `launch_groin_grid.py` | the full-model grid, through the unchanged runner | runs under `output/raw_runs/experiments/<study>/`, logs in `grid_logs/` |
| `score_groin_grid.py` | the option A grid against the observed fillet change | `grid_scores.csv` |
| `diagnose_no_groin_relaxation.py` | why the no-groin model erases the fillet | prints; saved as `diagnose_no_groin_relaxation.txt` |
| `blocking_groin_emulator.py` | approach 1: a groin that blocks a fraction b of transport | `blocking_groin_scores.csv` |
| `failure_schedule_test.py` | does an instant failure at the 2003 storm fit both windows? | `failure_schedule_scores.csv` |
| `trajectory_check.py` | the instant-failure candidates, at the observed dates | `trajectory_check.csv` |
| `score_instant_grid.py` | the full-model instant-2004 grids against the observed gap | `instant_grid_scores.csv` |

The emulators import each other: `blocking_groin_emulator.py` takes the
planform, climates and schedule constants from `groin_stability_option_a.py`,
which takes the solve from `../HAT_groin_solver_audit.py`;
`failure_schedule_test.py` and `trajectory_check.py` build on
`blocking_groin_emulator.py`. `WHERE_WE_LEFT_OFF.md` is the hand-over note of
2026-09-29.

**Question.** M = 60, f = 0.6 was fitted with the /10 planform at Hs 2.5 / Tp 8.0 / asym 0.7 /
high-angle 0.1. Is it still runnable under option A (Hs 2.0 / Tp 7.5 / asym 0.6 / high-angle 0.5,
metres planform), and what M does the old fillet now cost?

**Method.** `groin_stability_option_a.py` uses the solver audit's emulator (BRIE's own
`coast_diff` and sparse solve, no Barrier3D, no source/sink, no storms). It starts from the
BRIE shoreline after year 1 of the adopted `road_nobdm_nogroin` edgeBE runs (1996 and 2010), and
applies GroinCallback's dipole (GIS 6 updrift, GIS 5 downdrift) with the 1996→2003 linear ramp.
Fillet = groin run minus a no-groin run from the same planform. The old climate is judged on the
same planform ÷10, which is what it was calibrated on. The audit's straight-coast rig can't be
used: at 0° option A's diffusivity is −119 m²/yr, so it shuts down in year 1 at any M.
Runtime is about 30 s.

Outputs: `stability_summary.csv` (one row per case), `stability_traces.csv` (per year),
`run_log.txt`.

## Result

| | old /10 setup | option A |
|---|---|---|
| diffusivity at the groin cells (m²/yr) | 513,000 / 490,000 | 10,500 / 28,600 |
| groin-cell angles at start | −1° / −3° | −12° / −25° |
| 14-yr fillet at M = 60, f = 0.6 (1996) | 17 m | **216 m** |
| downdrift angle at end, M = 60 | −2.7° | **−37°** (peaks −42° in yr 4) |
| M giving the 22.11 m period-1 fillet | ~77 (emulator) | **~3** |

- **M = 60 does not run away, but it is ~10× too strong.** No cell shuts down for more than a
  year. The fillet overshoots to 244 m by year 3, then relaxes. The downdrift angle reaches the
  ~42° band where BRIE turns anti-diffusive. A single-cell shutdown appears at M ≥ 40 (f = 1) and
  at M = 60 (f = 0.6).
- **The same fillet now costs about 25× less M.** This is the `M / r_ipl` scaling the audit
  found, driven by the weaker diffusivity. For scale, the emulator puts the old setup's fillet at
  M ≈ 77 against the fitted 60. So the absolute M values here are approximate; the ratio between
  the two setups is what carries over.
- **Low M is well behaved.** M = 1–5 stays at −18° to −23° at the downdrift cell and is
  saturated by year 14 (fillet growing ≤ 0.6 m/yr). Both windows behave the same way.

## Consequences for the groin sweep

- The old grid (M 0–160) is entirely out of range. A re-sweep needs roughly M ∈ {0, 1, 2, 3, 4,
  5, 7, 10}.
- The sediment-budget objection ([[cascade-groin-M-not-affordable]]) largely goes away. M ≈ 3 is
  about 18,000 m³/yr, about 3% of the reach budget, against 52% at M = 50.
- The 22.11 m target is the 1984–2004 fillet. The 1996 target has to come from
  `observed_fillet_m(1996)`.
- Not tested here: cross-shore feedback, edgeBE, and storms. Before the sweep, a handful of full
  CASCADE runs at M ≈ 2–5 should confirm these numbers.

## Full-model grid (launched 2026-09-29 ~01:05)

`launch_groin_grid.py` runs the unchanged runner over M {1, 2, 3, 4, 5, 7, 10} × f {0.2, 0.4,
0.6, 0.8, 1.0} × 1996/2010. That is 70 runs: edgeBE, full_management, no historical relocations,
everything else at the code defaults. Output goes to
`output/raw_runs/experiments/groin/2026-09-29-option-a-grid/M<M>_f<f>/`, logs to `grid_logs/`.
M = 0 is the adopted matrix `road_bdm[_nourish]_nogroin` run. A relaunch skips finished cells.

**Why not the sweep worker.** `HAT_groin_sweep_worker.py` hardcodes the old waves
(2.5/8/0.7/0.1). Its parameter-file repair restores `.parameters_pristine.yaml`, which has no
per-cell dune ceilings. Its cells would not be the adopted model.

**Runner fix needed first.** Every groin run crashed in `predict_fillet`, because
`brie_r_ipl` read the diffusion number at 0°, which is negative under option A. Since 2026-09-29
the runner (.py and notebook) reads it at the groin cell's starting angle when a groin is
attached, and records the angle as `r_ipl_theta_deg`. No-groin runs are unchanged.

**Confirmation (smoke cell, 1996, M = 3, f = 0.6).** Full CASCADE gives an 18.9 m fillet in
14 yr against the emulator's 20.5 m. Nothing outside GIS 3-8 moves by more than 0.08 m. Most
of it is updrift advance (−17.9 m); the downdrift notch heals to 1 m.

**Scoring is still to do.** The runner's own extent check reports "no paired baseline", because
it looks for the baseline under the experiment tree. Score against the matrix runs. The
observed fillet CHANGE targets from `HAT_groin_sweep_config.observed_fillet_m` are
**−4.3 m (1996–2010) and −60.4 m (2010–2024)**. The fillet shrank in both windows. A module that
only traps (M ≥ 0) can match that only through the no-groin relaxation, so expect the best fit
near the low-M edge of the grid.

## Grid result (scored 2026-09-29 ~03:00, `score_groin_grid.py` -> `grid_scores.csv`)

106 cells: the 70-cell grid, plus an extension of M 15–40 × f 0/0.1/0.2. All ran clean. Scored as
total fillet change (the groin run's own D5−D6 OLS trend × 14 yr, built the same way as
`observed_fillet_m`). The old `measure_fillet` groin-minus-baseline reading is in the CSV too.

| | observed | no groin | groin needed |
|---|---|---|---|
| 1996–2010 | −4.3 m | **−73.9 m** | strong, persistent: M 10 / f 0.8 → −4.5 m (and along a ridge in M·f) |
| 2010–2024 | −60.4 m | −58.8 m | none: this window sees only M·f, and any M·f > 0 is worse |

**No (M, f) fits both windows.** The best joint cell (M 4, f 1.0) still misses by −34 m and
+39 m, with opposite signs. The "strong before 2003, gone after" corner (f ≈ 0) does not rescue it.
At f = 0 the fillet built in 1996–2003 collapses after the ramp, and the 1996 OLS trend gets
WORSE than no groin (−77 to −89 m).

**The problem is the no-groin model in 1996, not the groin.** With no groin the model erases the
Buxton fillet at ~5 m/yr in 1996–2010, while the observed fillet held. The 1996 downdrift cell
starts at −25°, the steepest angle in the reach, so this is plausibly option A's alongshore
diffusion flattening a real notch the model has no structure to hold. The same relaxation is
right in 2010–2024. Not yet diagnosed.

## Why the no-groin model erases the fillet (diagnosed 2026-09-29)

`diagnose_no_groin_relaxation.py` → `diagnose_no_groin_relaxation.txt`.

1. **It is alongshore diffusion.** The emulator (BRIE's alongshore solve alone, no Barrier3D)
   relaxes D5−D6 by −81 m (1996) and −72 m (2010). The full model gives −74 and −76 (natural).
2. **Nothing else touches it in 1996.** All 12 matrix runs (edge forcing zero/edge ×
   natural / road only / road + beach/dune × relocations) give −73.9 to −74.0 m.
3. **2010 matches only because of the Buxton nourishment.** Natural 2010 relaxes −76 m, like
   1996. The fills cut D6's retreat from +63 to +38 m, which gives −58.8 against the observed −60.4.
4. **The step is real geography, not an input artefact.** D5 sits 256 m landward of D6 in the
   dune-line build and 269 m in the CoastSat 1995–97 shoreline build. The neighbouring steps are
   60–140 m. GIS 5/6 is where the coast turns at Cape Hatteras.
5. **The observations say the cape flank changes much more slowly.** Model LRR 1996–2010: D5
   +1.67, D6 −3.61 m/yr, a gap closing at 5.3 m/yr. CoastSat: D5 +0.10, D6 −1.35, closing at
   1.45 m/yr. The wet/dry table: 0.3 m/yr. Same direction, but 4–17× too fast.
6. **Why the old M = 60 fit didn't see this.** Under the ÷10 offset the step was 25 m, so the
   relaxation was a tenth as large and a groin could hold it. The units bug hid the cape.

**Reading.** BRIE under option A has no process that keeps a cape. The real Cape Hatteras
curvature is sustained by things the model lacks: high-angle wave shaping, Diamond Shoals, and
the groins themselves. A groin fitted in 1996 (M 10, f 0.8) would be standing in for that missing
cape process. This is the double-counting [[cascade-groin-M-not-affordable]] warned about, now
running the other way. No (M, f) choice is made here; the options go to Hannah.

## In the code (2026-09-29)

- `cascade/groin.py`: adds `BlockingGroinCallback` (kind "blocking"; the callback form above, reading
  `cascade._brie_coupler._brie`). The schedule and validation moved into shared helpers with the same
  arithmetic, and the dipole is bit-identical to HEAD over both modes and 80 years. Tests are in
  `tests/test_groin_blocking.py` (11 pass).
- Runner `.py` and notebook, edited identically: the schedule is now **instant, from the 2004
  step, onset delay 35 yr** (was the 1996→2003 linear ramp). There's a `groin.kind` switch
  (`HAT_GROIN_KIND`) and `groin.blocking_b` (`HAT_GROIN_BLOCKING_FRACTION`, default 0.6). Blocking
  runs are named `groinblock`, so they can't collide with dipole runs. The index gains
  `groin_kind`/`groin_blocking_b`. For blocking, `groin_trapping_m_yr` is the mean applied rate.
  `reports.py` has blocking branches, with no a-priori amplitude or extent.
- The earlier `2026-09-29-option-a-grid` runs used the OLD ramp schedule. The code no longer
  reproduces them.

**Full-model check, 1996** (`experiments/groin/2026-09-29-instant-2004-check/`), gap change at
1997 / 2004 / 2008, observed +1.7 / +16.0 / −9.9:

| run | 1997 | 2004 | 2008 | emulator said |
|---|---|---|---|---|
| blocking b = 0 | −18.8 | −69.3 | −80.6 | **bit-identical to the matrix no-groin run** |
| blocking b 0.6, f 0.3 | +3.8 | **+14.8** | −31.9 | +6 / +25 / −19 |
| dipole M 10, f 0.2 | −1.0 | +4.3 | −30.2 | +4 / +22 / −12 |

The blocking groin holds the gap through 2004 in the full model, as intended. After the failure,
both kinds relax 12–20 m more than the emulator predicted: the additive Barrier3D offset is only
approximate. The fit has to be made in the full model; the emulator only picks where to look.

## Full-model grid, instant 2004 failure (2026-09-29, `score_instant_grid.py` → `instant_grid_scores.csv/.txt`)

100 runs, all clean: blocking b 0.4–0.8 × f 0.2–0.6 and dipole M 6–14 × f 0.1–0.5, in both
windows. Filed under `experiments/groin/2026-09-29-instant-2004-grid/`. Scored by date-RMSE of
the D5−D6 gap change against the wet/dry observations.

| | 1996 RMSE | 2010 RMSE | joint | equivalent trapping |
|---|---|---|---|---|
| no groin | 65.0 | 30.2 | 50.7 | — |
| **blocking b 0.6, f 0.4** | 9.1 | 15.8 | **12.9** | 4.5 m/yr (1996), 1.85 (2010), emergent |
| **dipole M 12, f 0.3** | 2.4 | 17.4 | **12.4** | 12 m/yr, then 3.6 |

Both optima are interior. They are equally good: the 0.5 m gap in joint RMSE is far below the
observations' 11–14 m year-to-year scatter.
- 1996: both follow the observed hold-then-drop (obs +2/+16/−10; blocking +4/+15/−26; dipole +3/+20/−11).
- 2010: both roughly halve the no-groin error, but both miss the late drop. At 2023 the model has
  −28 (blocking) and −21 (dipole) against −48 observed; the OLS change is −35/−25 against −60. The
  no-groin run gets the trend but starts falling a decade early.
- Budget, at the runner's profile height (9,384 m³ per m/yr): blocking ≈ 42,000 m³/yr in 1996
  (7% of the reach budget); dipole ≈ 113,000 m³/yr (19%). The blocking groin gets the same fit
  while moving ~2.5× less sand. (The earlier "M 3 ≈ 18,000 m³/yr" in this README used the
  wrong height; the runner reports 28,000.)

## The scripts in detail

The header of each script says what it does and how to run it; the reasoning
behind it is here, with the original header kept word for word. Moved out of
the scripts on 2026-10-01, when the groin study was brought in line with
`scripts/STYLE.md`.

### groin_stability_option_a.py

Why it exists, what it is and is not, and the result are in the sections above
("Question", "Method", "Result"). In the code:

- `PLANFORM_RUNS` are the adopted natural-dunes-with-road (`road_nobdm`)
  edgeBE runs. Any matrix run gives the same BRIE planform at year 1 to within
  the year's cross-shore change.
- The groin wiring is GroinCallback's, as in `HAT_groin_sweep_config`: GIS 6
  updrift, GIS 5 downdrift, 15 buffer domains, first GIS 1.
- The old climate was calibrated on the /10 planform, so it is judged there:
  every alongshore difference one tenth as large.
- The M that reproduces the period-1 fillet size is found per climate by
  linear interpolation over the stable hindcast cells at f = 0.6.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
Is the groin dipole stable under the option A wave climate, on the real coast?

WHY THIS EXISTS
    M = 60, f = 0.6 was fitted under the /10 planform at Hs 2.5 / Tp 8.0 /
    asym 0.7 / high-angle 0.1. Option A (adopted 2026-09-27) is Hs 2.0 / Tp 7.5 /
    asym 0.6 / high-angle 0.5 on the metres planform. The solver audit
    (../HAT_groin_solver_audit.py) found the fillet is bought by M / r_ipl and
    that the high-angle fraction sets where BRIE's diffusivity goes to zero, so
    both changes bear directly on whether the old M is even runnable.

    The audit's rig starts from a STRAIGHT coast. Under option A BRIE's
    diffusivity at 0 deg is -119 m^2/yr, so a straight rig shuts down in year 1
    at any M and says nothing about Hatteras. This check starts instead from
    the real BRIE shoreline of the adopted matrix runs (x_s after year 1), whose
    local angles (median about -9 deg) are what set the diffusivity.

WHAT IT IS AND IS NOT
    The same emulator as the audit: BRIE's own coast_diff table and sparse
    indices, the same row-scaled implicit solve, the same clip at zero. It adds
    only the groin dipole, with GroinCallback's linear-ramp deterioration
    schedule. No Barrier3D, no source/sink, no storms: every number here is the
    alongshore solve's response to the dipole alone, measured against a no-groin
    run from the same planform. Anything a full CASCADE run does differently is
    cross-shore feedback.

    Nothing in the main code is changed; brie and the audit are imported as-is.

Author: Hannah A. Henry, UNC CECL
```

<details><summary>Function notes (the original docstrings)</summary>

**`effective_M()`**

```text
GroinCallback._effective_trapping_rate, linear_ramp mode.
```

**`load_planform()`**

```text
BRIE x_s (metres, landward +) after the first model year.
```

**`solve()`**

```text
Emulated BRIE alongshore solve with one dipole; returns per-year rows.
```

**`run_case()`**

```text
Groin run minus the no-groin run from the same planform.
```

**`summarise()`**

```text
One row: end fillet, whether anything shut down, when.
```

</details>

### launch_groin_grid.py

Why the runner and not `HAT_groin_sweep_worker.py`, what is set, and how cells
are filed: see "Full-model grid" above and the original header below. In the
code, `STUDY` and `KIND` default to the first grid (dipole, with the
1996->2003 linear-ramp schedule the runner had until 2026-09-29); later studies pass
`--study` and `--kind`. Cells run low M first in both windows, so an early
stop still leaves the plausible end of the grid.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
The option A groin grid, M x f x window, run through the unchanged runner.

WHY THE RUNNER AND NOT HAT_groin_sweep_worker. The worker is a second copy of
the runner, and as of 2026-09-29 it is not the adopted model: it hardcodes the
old waves (2.5 / 8 / 0.7 / 0.1), and its parameter-file repair restores a
snapshot with no per-cell dune ceilings. Driving the runner means every cell
is the adopted model by construction. Runs write their own parameters file,
so they are safe alongside other sessions' runs (Hannah, 2026-09-29).

WHAT IS SET. Only the groin and the filing; every other value is the code
default (HAT_IGNORE_SETTINGS=1, stray HAT_* dropped, as HAT_run_all does):
edgeBE, full_management without historical relocations, option A waves,
v3_trim24 storms. M = 0 is not run: the paired baseline is the adopted matrix
run with the same tokens, which the runner resolves by itself.
    1996  HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin
    2010  HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin

The M = 2-5, f = 0.6 cells are the confirmation runs for the emulator in
groin_stability_option_a.py.

FILING. Run names carry `groin` but not M or f, so each cell is its own
experiment member: raw_runs/experiments/groin/2026-09-29-option-a-grid/M<M>_f<f>/.
A cell whose metadata already exists is skipped, so a relaunch resumes.

    python launch_groin_grid.py [--streams 3] [--only 1996:3:0.6]

Author: Hannah A. Henry, UNC CECL
```

### score_groin_grid.py

The two readings and why both are kept: see "Grid result" above and the
original header below. The sign is landward-positive throughout: + means the
downdrift domain sits landward of the updrift one.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
Score the option A groin grid against the observed Buxton fillet.

Two readings, both reported, because they answer different questions:

  groin_contribution_m   end-year (D5 - D6) of the groin run minus the same for
                         its no-groin baseline. This is HAT_groin_sweep_config.
                         measure_fillet, the metric the M = 60 fit ranked on.
                         It is >= 0 for any trapping groin.
  total_change_m         the groin run's OWN (D5 - D6) change over the window,
                         OLS slope x window length. This is how the observation
                         is built (observed_fillet_m: OLS across the wet/dry
                         dates x window), so it is the like-for-like reading.
                         It includes the relaxation the coast does with no
                         groin, which is what an observed shrinking fillet needs.

Both are scored against observed_fillet_m(period). The sign is landward-positive
throughout: + = the downdrift domain sits landward of the updrift one.
M = 0 is the adopted matrix no-groin run.

Author: Hannah A. Henry, UNC CECL
```

### diagnose_no_groin_relaxation.py

The diagnosis it supports is "Why the no-groin model erases the fillet" above.
It has no `main()`: it runs top to bottom, from the repo root, because its
paths are relative to the working directory. Part 2 imports
`groin_stability_option_a.py` for the emulator.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
Why the 1996 no-groin run erases the Buxton fillet (2026-09-29). Run from the repo root.
Part 1: D5-D6 change across every adopted matrix run (edge forcing x management).
Part 2: the same with alongshore diffusion alone (the emulator, no Barrier3D).
```

### blocking_groin_emulator.py

The dipole groin imposes +/-M metres a year whatever the state, so it cannot
hold the 256 m step BRIE flattens at GIS 5/6. A physical groin intercepts the
transport arriving at it; this emulator removes a fraction b of the D5|D6
face's share of BRIE's alongshore solve, in two implementations (exact, which
needs a BRIE change, and callback, which does not). The full-model no-groin
gap change in `FULL_NOGROIN` is from the adopted matrix road_bdm[_nourish]
runs. In the solve, each row couples to both neighbours with its own r[i]
(BRIE's row scaling): the face D|U is row D's upper link and row U's lower
link. The matrix is assembled exactly as BRIE does, but from the two weight
vectors, so the unblocked case reproduces BRIE's matrix term for term.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
Approach 1: a groin that BLOCKS a fraction b of alongshore transport.

The dipole groin (GroinCallback) imposes +/-M metres a year regardless of
state, so it cannot hold the 256 m step BRIE flattens at GIS 5/6 (see the
diagnosis in README.md). A physical groin intercepts the transport arriving at
it. Here the groin removes a fraction b of the D5|D6 face's share of BRIE's
alongshore solve, with b ramping b0 -> b0*f over 1996-2003 (the same schedule
GroinCallback uses for M).

TWO IMPLEMENTATIONS, because only one of them can go into CASCADE unchanged:

  exact      Scales the face's coupling in BOTH halves of BRIE's
             Crank-Nicolson step (the explicit Laplacian on the right-hand
             side and the implicit matrix) by (1 - b). The upper bound on
             what the idea can do. Needs a BRIE change.
  callback   What GroinCallback's hook can do: a pre-solve x_s_dt correction.
             It cancels b times the explicit estimate of the face's full
             step, 2 * r_i * (x_j - x_i) on each side (the explicit half
             doubled to stand in for the implicit half). Needs no BRIE change.

SCORING. The D5-D6 gap change over 14 yr, OLS through 15 states (t = 0..14),
against observed_fillet_m: -4.3 m (1996-2010), -60.4 m (2010-2024). The
emulator has no Barrier3D and no nourishment, so each window also gets an
ADJUSTED score: the emulator change plus (full-model no-groin minus emulator
no-groin). That offset is +2 m in 1996 and +16 m in 2010, the 2010 one being
the Buxton fills. It assumes those processes add linearly, which the
full-model grid has to confirm.

Author: Hannah A. Henry, UNC CECL
```

### failure_schedule_test.py

Tests whether an intact groin that fails at once after the 2003 storm fits both
windows, where the linear 1996->2003 wear-down put the decline inside the 1996
window, in which the data show none.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
Does an INSTANT failure at the 2003 storm let one groin fit both windows?

The observed D5-D6 gap (wet/dry table, 24 dates) grows 1967-1995, holds at
134-155 m through 2004, then falls to 125 (2008), 104 (2016), 63-74 (2019-23).
That is an intact groin failing after the 2003 storm, NOT the linear
1996 -> 2003 wear-down both groin emulators were given -- which puts the decline
inside the 1996 window, where the data show none.

Tests both representations under GroinCallback's existing "instant" mode (full
strength until the failure year, x f from then on; no new code needed for the
schedule):
  blocking  approach 1, callback implementation (blocking_groin_emulator.py)
  dipole    today's GroinCallback, M m/yr

Scored as there: adjusted OLS gap change vs observed_fillet_m.

Author: Hannah A. Henry, UNC CECL
```

<details><summary>Function notes (the original docstrings)</summary>

**`dipole()`**

```text
GroinCallback's +/-M dipole on the same solve (x_s_dt only).
```

</details>

### trajectory_check.py

It has no `main()`: it runs top to bottom. It imports
`blocking_groin_emulator.py` and `failure_schedule_test.py`, and swaps the
instant schedule into the emulator's `b_of_year`.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
Trajectory check for the instant-failure candidates (see failure_schedule_test.py).

An OLS trend can be matched by the wrong shape. The observed D5-D6 gap is FLAT
1996-2004 and then declines; a dipole builds a fillet and then drops it. So the
candidates are scored against the observed gap at its own dates, as change since
the window start (the 2010 start is interpolated between 2008 and 2014), with
the full-model no-groin offset spread linearly over the window.
```

### score_instant_grid.py

The scores behind "Full-model grid, instant 2004 failure" above.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
Score the full-model groin grids (instant 2004 failure) on the observed gap.

For every run under raw_runs/experiments/groin/2026-09-29-instant-2004-grid/
(blocking b<b>_f<f>, dipole M<M>_f<f>) and the matrix no-groin baseline:

  date RMSE    the D5 - D6 gap change since the window start, sampled at the
               wet/dry table's own dates inside the window, against the
               observed change (2010 start interpolated between 2008 and 2014).
               The primary score: an OLS trend alone is matched by
               build-then-collapse shapes the data do not show.
  OLS change   trend x 14 yr, as observed_fillet_m builds its target.

Both windows, both kinds; joint = RMS of the two windows' date RMSEs.
Sign: landward-positive, + = downdrift sits further landward of updrift.

    python score_instant_grid.py  ->  instant_grid_scores.csv + printed tables

Author: Hannah A. Henry, UNC CECL
```
