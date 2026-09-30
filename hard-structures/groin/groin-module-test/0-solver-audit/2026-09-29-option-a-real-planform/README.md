# Groin stability under option A, on the real planform (2026-09-29)

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
