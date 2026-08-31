# Calibration figures — what was tested on 2026-08-30, and what it showed

Five summary figures documenting the source/sink and groin calibrations, and
**how they were tested**.

The figures themselves are written to `output/groin_sweep/figures/`, which
`.gitignore` does not track — this file lives here so the reasoning survives
even when the PNGs do not. Regenerate them with:

```
python scripts/hatteras_ms/groin-sweep/HAT_calibration_summary_figures.py
```

The conclusions live in prose, in `scripts/hatteras_site_config.py`. These
figures exist because two things were recorded nowhere else at all — see
"What was at risk" below.

## The figures

| file | shows |
|---|---|
| `fig_three_targets.png` | the same 61 cells scored three ways: the fillet says M = 95, the D1–D12 profile says M = 0, D4–D8 demeaned says M = 60 |
| `fig_Mf_identifiability.png` | D4–D8 RMSE over the (M, f) grid with iso-M·f contours — **refutes** the "only the product is identified" claim |
| `fig_top_profiles.png` | the top cells and the no-groin baseline against the observed change profile, fit window marked |
| `fig_be_convergence.png` | interior RMSE per calibration pass, both periods, GIS 90 re-solve marked |
| `fig_period2_and_bug.png` | why period 2 is not fitted, and what the topography-product bug was worth |

Pre-existing, not from this script: `period1_top_n_profiles.png` (the
production D4–D8 ranking) and the four `joint_*` surfaces.

## The result these were built to show

**The fitted groin is set by the target, not by the data.** The same model runs
give M = 95, M = 0 or M = 60 depending only on what they are scored against:

* **the fillet (D5−D6 scalar)** — no admissible M matches it on this grid
  (`HAT_groin_timeseries_check.py:29`). Fitting it anyway rails `f` at the grid
  bound and makes the best M swing 95 → 160 between adjacent `be1` values.
* **the D1–D12 change profile** — ranks M = 0 best, monotonically. Not because
  there is no groin: D2–D4 is Cape Point accretion the parameterisation does
  not represent, and D6–D7 is an erosion trough peaking one domain *north* of
  the structure, which a groin actively worsens by pushing D6 seaward. A ~17 m
  groin signal is swamped. Narrowing to D5–D7 does not help — D6–D7 is inside
  that window too.
* **D4–D8 demeaned** — a clean interior minimum at M = 60, rising on both sides
  through M = 160. Demeaning removes the level offset the source/sink term
  owns; D1 is excluded because the cape's 81–104 m change is ~5× the groin's
  signal.

**M = 60, f = 0.6 stands.** Best cell 10.09 m against 10.5 m for the production
pair, versus 12.2 m with no groin. It was briefly marked void on 2026-08-30
because the sweep had been building 1984 cells on the 2004 island; that bug was
real and is fixed, but re-scoring on the corrected island moved the answer by
less than the width of the ridge it sits on.

**A claim that did not survive testing — and a correction that also failed.**
`HAT_period1_top_n_figure.py`'s caption says the metric "identifies a product,
not a pair". `fig_Mf_identifiability.png` was built to test that and refutes it:

```
corr(RMSE, M·f)  −0.066        equal-product cells at M·f ≈ 40:
corr(RMSE, M)    +0.610          M40/f1.0 = 10.4   M50/f0.8 = 10.7
corr(RMSE, f)    −0.493          M70/f0.6 = 11.2   M95/f0.4 = 12.5
```

But the replacement claim — "M and f are each weakly constrained" — is wrong
too. `hard-structures/groin/GROIN_PLAN.md` has the answer: the ridge is in
**period-1 cumulative trapping, M(15.5 + 4.5f)**, not in M·f. Period 1 mostly
precedes the 1996–2003 deterioration ramp, so f moves it by only 29% across its
whole range, while period 2 is 20·M·f where f = 0 gives zero. Hence **M is set
from period 1 and f from the 1967 rig and period 2** — and at f = 1.0 the best
M is 50 while at f = 0.6 it is 70, so a poor score at M = 50, f = 0.6 means M
was too low for that f, not that f = 0.6 is wrong.

**The groin does not reproduce the shape.** Inside the fit window the observed
profile peaks at D6, dips at D7 and peaks again at D8; every cell draws a
smooth monotonic rise. The RMSE gain comes from matching the overall D4–D8
slope. Read the residual as the split between what the groin explains and what
the source/sink calibration absorbs — not as a successful shape fit.

**Period 2 is not fitted — but the groin still runs in it.** `GroinCallback`
carries an absolute calendar timeline (install 1969, deterioration onset 1996,
ramp to 2003, then hold at M·f), so no period-specific configuration exists.
Running it in period 2 is right for consistency of the structure's timeline,
not because it explains that period's shoreline.

Period 2 records a **release** the module cannot produce: −76 m, of which
`GROIN_PLAN.md` attributes ~85% to the **updrift** side eroding once the
structure failed, not to impounded sand draining downdrift. Trapping is bounded
at ≥ 0, so the groin can stop adding sand but cannot drain the fillet. The groin
supplies ≈ +17.2 m of period 1's observed +52 m (33%) and ≈ 0 of period 2's
−76 m; the rest is carried by the source/sink calibration and the Cape Point
dynamics the dipole does not represent.

So the 2026-08-30 joint fit railing at M = 160 was not a failure to repair —
fitting period 2 is the wrong thing to attempt.

## What was at risk, and is now recorded

* **The BE convergence sequence.** `--overwrite` replaces a run's row in
  `run_index.csv`, so only the final pass survives there, and
  `convergence_history.json` still carries 2026-08-24 baselines from the
  pre-restructure topography. The pass-by-pass numbers existed only as text in
  a comment. `fig_be_convergence.png` is now their durable record — and the
  numbers in it are transcribed from that comment, not re-derived, which is
  stated in the script.
* **The three-target comparison.** Nothing on disk showed it.

## Two caveats on the figures themselves

1. **~~Panel (b) of `fig_three_targets.png` is not a controlled comparison.~~
   WITHDRAWN 2026-08-31 — it is.** This said panel (b) came from the 40-year
   window on the PRE-FIX topography while (a) and (c) were period 1 on the
   corrected one. True when written; the note then outlived the problem. The
   fullperiod sweep was re-run 2026-08-30 18:20, five hours **after** the
   worker's topography fix (`562c75c`, 13:01), and the figure was rebuilt at
   23:52 from those cells.

   Verified rather than assumed: re-running cell `M60_f0.50` through the worker
   reproduces its stored result to **8.3e-05** on rates of ~2.9 m/yr — the same
   noise floor a period-2 cell shows against *itself* (1.0e-04 between two
   re-runs), which the 1984–2024 window inherits because it contains period 2.
   A wrong-island run would differ in the first or second decimal, not the
   fifth. All three panels are on the same corrected topography.
2. **Absolute RMSE here differs from `period1_top_n_profiles.png`** (10.1 vs
   11.4 m for the best cell). This script approximates change as
   `rate × 20 yr` from `sweep_results.jsonl`; the production script reads the
   shoreline matrix directly. The *ranking and the shape of the curves* are
   unaffected, and the production figure is the one to quote.
