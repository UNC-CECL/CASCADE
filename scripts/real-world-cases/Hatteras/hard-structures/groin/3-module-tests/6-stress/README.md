# 6-stress — trying to break the runner's groin

**Question.** Tests 1-6 asked whether the groin module behaves under the conditions it was built for. This one
asks where it stops behaving: push each input past where the Buxton hindcast uses it, and see whether the
groin keeps the sand balance, stays stable, and settles. Asked for on 2026-10-15 ("try to break it"). What
came out of it, the module's tested range and the guard that now enforces it, is at the end.

**What is tested.** The runner's own groin, `cascade.groin.GroinManager`, imported as it stands in
`cascade/` (not a copy, unlike tests 1-6): the drift groin plus the leak correction against the tracked
groin-free shoreline, with the runner's 10-cell correction zone. It is driven on the straight-coast rig of
`../2-straight-coast/1-run/`: 61 domains of 500 m, one groin at the middle face, BRIE's implicit alongshore
step alone on an infinitely long coast at angle θ0 (the tilt is added only where angles are read, so the
periodic wrap never meets the groin). The groin reads a stand-in for BRIE carrying a real Brie
instance's tables and the rig's shoreline, and gets `x_s_dt` as a list, as CASCADE hands it over. No
Barrier3D, storms or sea level, except where the background family imitates what Barrier3D passes in.

**Rig check.** At strength 0.6 and a 5-cell zone, on a straight coast, the groin reproduces the wave-driven groin of
`1-run/` to 0.0 m at all 17 orientations, in both the updrift cell and the coast's total
(`tables/stress_rig_check.csv`). On a straight coast the tracked groin-free shoreline is the straight line, so
the two must agree; they do.

```
stress_test.py      the rig, every family, the scores, and when the module's guard would have stopped each
                    run (about 3 minutes)
stress_figures.py   the three figures, from the tables
tables/             stress_runs.csv (one row per run: scores, r at the groin, guard_year and guard_reason),
                    stress_by_year.csv.gz (every year of every run, groin and no-groin), stress_rig_check.csv
figures/            stress_1_where_it_breaks, stress_2_other_stresses, stress_3_two_failure_modes;
                    PDFs and CAPTIONS.md in supporting/
logs/run_log.txt    the 2026-10-15 run's printed summary
```

    python stress_test.py [--only waves]
    python stress_figures.py

## The families and the rules

Each case runs at strengths 0.25 (Buxton's whole-life value) and 1.0 (all of the arriving drift) at seven
orientations, −30° to +30°, for 100 years unless stated, and is paired with the same run with the groin off.
Buxton's two groin cells start at about −12° and −25°.

| family | what changes |
|---|---|
| waves | one wave parameter at a time: high-angle fraction 0.1-0.7, asymmetry 0.5-0.9, Hs 1.0-3.5 |
| strength | 0 to 1 at option A |
| notch | a step already at the groin when it starts, −300 to +300 m (Buxton's is 256 m) |
| switch-off | on for 50 years, off for 100 |
| cycling | on and off every 5 years |
| variable | a new climate each year (asymmetry 0.40-0.80, Hs 1.5-2.5; below 0.5 the drift reverses); the groin as written (drift table read once) and with the table refreshed each year; 3 seeds, at −20°, 0°, +20° |
| background | `x_s_dt` arriving non-zero, as from Barrier3D: 2 m/yr retreat; ±10 m/yr noise per domain; a 30 m cut every 10 years |
| placement | the groin face near the end of the array (lower domain 0, 1, 2, 5, 9, 58, 59 of 61) |
| long | 500 years |

The rules, set before the runs: a run **passes** if it finishes; its **net sand** against the no-groin run is
within 5 m or 5% of the groin's gross change, whichever is larger; no face is steeper than 60° (BRIE's table
ends at 90°); and, for constant forcing, the gap moves less than 10% (+5 m) over the last fifth of the run.

## Findings

These are the module's results with its guard recording, not stopping: the guard notes when it would have
stopped each run and lets it finish, so the outcome can be checked against it.

1. **At option A it holds.** Strength 0 is the no-groin run exactly. Up to strength 0.5 it conserves at every
   orientation, −30° to +30°. At 0.25 the largest net change is +34 m in 100 years at +10° (5% of its gross).
   It settles, and over 500 years nothing runs away. At strength 1 it leaks 67 m at +20° in 100 years and
   157 m at 0° over 500.
2. **It breaks two ways, both set by BRIE's diffusion number r = D dt / 2 dy² at the groin** (figure 3):
   - **r too large: the leak correction overshoots.** The correction is an explicit term in front of BRIE's
     implicit step. At high-angle fraction 0.3 the groin loses 434 to 3,250 m of summed shoreline in 100
     years at −30° to −10°. At 0.1 it runs until an angle leaves BRIE's diffusivity table and the groin
     crashed (`IndexError: index 180`; fixed, below). Hs 2.5 and above does the same on coasts turned with
     the drift: 22 to 67 km. On the rig, losses start once r at the start angle is above about 0.13.
   - **r = 0: the coast is frozen.** At asymmetry 0.7 and above, BRIE's diffusivity is clipped to zero at
     every negative angle. Nothing carries sand back across the face, but the groin keeps moving its share of
     the drift every year, so the gap grows without limit: 2,070 m (asymmetry 0.7) and 3,360 m (0.8) by
     year 100 at −20°. The same happens on any coast once the groin's own face steepens into the angles where
     BRIE's diffusivity is zero. This is the fixed-shift groin's runaway from `../TEST_PLAN.md` §2.1, now in the drift
     groin; a real groin fills to its tip and bypasses, and this module has no length.
3. **Waves that change from year to year break it at any r.** 35 of 36 variable-wave runs leak: up to
   1,046 m at strength 0.25 and 104 to 3,841 m at 1.0. Refreshing the drift table each year does not help,
   and neither does widening the correction zone to the whole coast or re-stepping last year's footprint with
   last year's tables (checked by hand, not in the tables). Varying Hs alone conserves (net 2-21 m against
   540-640 m gross); varying asymmetry, even within 0.5-0.65, leaks 40-330 m. The likely reason is that years
   with r ≈ 0 at the groin pile the fillet up steeply and the next high-r year meets it. Not proven.
4. **Starting from a sharp step is fine on the coasts that matter.** From −300 to +300 m it conserves from
   −30° to 0°, which includes Buxton's angles. At +10° to +30° it leaks 20-92 m at strength 0.25 and up to
   222 m at 1.0.
5. **Switching off once is clean; switching on and off repeatedly is not.** Off after 50 years it conserves
   at strength 0.25 except at +10° (23 m). Cycled every 5 years it leaks up to 39 m (0.25) and 458 m (1.0).
6. **What Barrier3D passes in mostly does not matter.** Steady retreat and the storm pulse give exactly the
   no-background result (they add linearly). Per-domain noise leaks 14-109 m at 0°, +20° and +30°.
7. **The correction zone does not wrap the array's ends.** A face within 10 cells of either end leaks up to
   1,528 m. A face at 9 passes. Buxton's face sits at pads 19|20 of 120.

**Checked against the full model.** The groin-on wave sweep (`output/raw_runs/sensitivity/`, run 2026-10-14/15
before the guard) shows the same two modes at Buxton. Its asymmetry 0.7 and 0.8 cells freeze: r at the groin
face is 0 in every year, the face steepens from 28° to 56° and 60° over 1996-2009, and the offset across it
grows from −264 to −725 and −857 m. Its high-angle 0.3 cells reach r 0.32 at the face and 0.44 in the zone.
Asymmetry 0.65 comes close (r at the face down to 0.003, the face to 41°). Every other cell keeps r at the
face between 0.012 and 0.185, and the face at or under 31° (the option-A matrix: r 0.032-0.068 at the face,
0.114 at most in the zone, face 16-29°). Per-run values: `tables/r_at_buxton.csv` (`r_at_buxton.py`).

## The bounds, and the guard (2026-10-15)

The module is meant to be a simplified groin, so rather than being extended to work everywhere it now stops
the run, with `GroinOutsideTestedRange`, when the run leaves the range it was tested in. Each year the groin
is installed:

| check | stops the run when | why this limit |
|---|---|---|
| wave climate | BRIE's diffusivity table differs from the first year's | finding 3; the drift table is read once |
| sand balance | the groin's footprint (the shoreline with it less the tracked groin-free shoreline) sums to more than 10 m, or 5% of its gross, whichever is larger | finding 2, first mode. 10 m rather than the scoring rule's 5 m: at 5 m a Buxton-like 300 m starting step trips it in years 2-3 while the groin-free shoreline relaxes its own step |
| groin face | the face between the two groin cells is steeper than 45° | finding 2, second mode. Every well-behaved full-model run stays at or under 31°, asymmetry 0.65 reaches 41°, and the broken cells pass 45° in their fourth year and keep going. 60° (the scoring rule) let asymmetry 0.8 run to 59.7° unstopped |

When it is built, it also refuses a correction zone that runs past either end of the domains (finding 7).

**How well it separates** (`tables/stress_runs.csv`, `guard_year` and `guard_reason`). It would have stopped
225 of the 227 runs that fail, at a median of year 2. The two it lets through are borderline: high-angle 0.1
at −10° (the balance closes within 2%, but the gap is still moving) and a 150 m starting step at +30° at
strength 1 (8 m against a 79 m gross). It also stops 90 of the 369 runs that pass: 23 by the zone rule, 40 by
the face rule (all at strength 0.75-1 on −20° and −30° coasts, where the face settles at 45-51°), and the rest
by the balance rule, in runs that broke it for some years and came back within tolerance by the end. Of the
42 Buxton-like rig runs (strength 0.25 or less, −30° to −10°) it stops one, which starts with a 50° face (a
300 m step on a −30° coast).

**In the full model** (`output/raw_runs/experiments/code-checks/2026-10-15-groin-guard-check/`, NOTE.md there): both
Buxton windows at option A run with the guard silent, and their shorelines are identical to the matrix runs
(0.0 m). Asymmetry 0.8 is stopped in 2000 (face 46.3°); high-angle 0.3 in its third year (10.8 m lost
against 111.2 m gross).

**The angle clip.** The groin's re-run of BRIE's step clipped the angle index to `_wave_climl`, one past
the end of BRIE's diffusivity table, and crashed at faces steeper than about 89.5°. It now clips to the
table's last entry, as the rest of the groin already did; results are otherwise bit-identical (0.0 m).
BRIE's own step (brie.py line 1297) has the same clip and is not changed here.

## What this means for reporting

Within these bounds the module conserves sand and settles; outside them a run stops rather than returning a
wrong answer:

- **One wave climate for the whole run.** No year-to-year waves.
- **BRIE's diffusion at the groin weak but not zero.** In wave terms, one parameter at a time from option A,
  at Buxton: Hs 0.75-3.0, Tp 6-12 s, high-angle fraction 0.45-0.55 and asymmetry 0.5-0.65 ran without
  breaking in the full model. High-angle 0.3 and asymmetry 0.7-0.8 break it, and the guard stops them.
  High-angle 0.4 sits between and was not run with the guard.
- **A groin face no steeper than 45°.** Buxton's is 16-29°.
- **Strength up to 0.5** conserves on every coast angle tested; 0.75-1 leaks on coasts turned with the drift
  (+20°) and drives the face past 45° on steep coasts. Buxton's fitted strengths are 0.10-0.25.
- **One switch-off.** Repeated on-off cycles leak.
- **The groin at least `half_width - 1` (9) cells from either end of the domains.**

The groin-on wave sweep's asymmetry 0.7 and 0.8 cells (`output/raw_runs/sensitivity/waveasym/`) ran before
the guard and measure the groin breaking, not the coast's sensitivity to asymmetry. Rerun now, they stop.
