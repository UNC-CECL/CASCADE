# `output/calibration/groin` — what is in here, and what it concluded

This was `output/groin_sweep/` until 2026-09-18. Every script builds its path from
`HAT_groin_sweep_config.GROIN_SWEEP_ROOT`, so the next move is a one-line change.

Map of this directory, written 2026-08-30. **If you only want the answer, read
`SELECTED_M60_f0.60/README.md`.**

---

## The answer

**M = 60 m/yr, f = 0.6.** Locked 2026-08-30, single value, no ensemble.

| | value | set by |
|---|---|---|
| `trapping_M` | **60** | period-1 D4–D8 demeaned profile fit at be1 = −42.6, reproduced independently by D3–D9 |
| `deterioration_f` | **0.6** | the 1967–2018 rig, the only window spanning the deterioration ramp |

They come from **different windows on purpose** — period 1 fixes M and barely
sees f (its cumulative trapping is `M(15.5 + 4.5f)`, which f moves by 29%
across its whole range); period 2 is `20·M·f` and carries no information about
M at all (RMSE identical to four significant figures for every M from 0 to 160).

Conditional on be1 = −42.6, void without edgeBE, and M ranges 40–95 across
defensible fit windows. Quote it as a **calibrated parameter pair with those
conditions**, not as a measured property of the structure.

Authority: `hard-structures/groin/GROIN_PLAN.md`. The prose reasoning lives in
`scripts/site_layer/hatteras_site_config.py`; the figure-by-figure record is
`scripts/hatteras_ms/groin-sweep/CALIBRATION_FIGURES.md`.

---

## Directory map

| path | what it is | read it for |
|---|---|---|
| **`SELECTED_M60_f0.60/`** | **the decision, curated** | the case for the pair, in order, with figures and gifs |
| `joint_fit.json` | **the file the pipeline reads** | what stage 6 and the LOESS analysis will actually run |
| `1984_2004_edgeBE/` | period-1 sweep, 488 cells over (M, be1, f) | **the fit that sets M** |
| `1984_2004_zeroBE/` | period 1 without edge forcing | the negative control — the groin just absorbs the missing source/sink term and improves to the grid edge |
| `2004_2024_edgeBE/` | period-2 sweep | that period 2 carries no information: every M scores 14.90 |
| `2004_2024_zeroBE/` | period 2 without edge forcing | same, at 22.42 |
| `fullperiod_1984_2024/` | continuous 40-year sweep, 43 cells | **a bound on M from above, not a fit** — see the warning below |
| `_validation_1984/`, `_validation_2004/` | drift guards | that a sweep cell reproduces its reference run |
| `figures/` | every figure any script has written here | the working set; `SELECTED_` holds the curated subset |

Nothing inside was renamed during the 2026-08-30 cleanup or the 2026-09-18 move.
Figure scripts open these subdirectories by name (e.g. `HAT_d4d7_window_figure.py`
opens `1984_2004_edgeBE/sweep_results.jsonl`), so renaming the sweep
directories would break them. The organisation is additive: `SELECTED_` and
this README sit on top of the existing layout.

---

## Three traps in this directory

**1. `fullperiod_1984_2024/results.csv` ranks M = 0 best. Do not read that as
evidence against the groin.** That window NETS period 1's fillet build (+52 m)
against period 2's collapse (−76 m), and a module whose trapping is bounded at
≥ 0 can only widen the gap — so it can never win there however it is scored.
Its RMSE floor is ~90 m, five times the groin's ~17 m signal. The sweep is
still useful as an upper bound on M. It is not a fit.

**2. `joint_fit.json` is PINNED, not fitted, and stage 5 will overwrite it.**
Its own ranking returned edgeBE M = 160 / f = 0.8 and zeroBE M = 140 / f = 1.0,
both railed at a grid bound, because it scores both periods jointly and period 2
records a release the module cannot produce. That ranking is kept at
`output/archive/2026-08-30_groin-railed-ranking/joint_fit_RAILED_ranking_20260830.json`. `HAT_run_all.py` stage 6
reads `joint_fit.json` and passes whatever it holds to every groin run in the
matrix — **re-run stage 5 and you must re-pin afterwards.**

**3. A run directory does not name its own parameters.** The rig sweep writes
every cell into one run name, so whatever survives in `output/raw_runs/` is the
last cell that finished, not the production run. On 2026-08-30
`HAT_1967_2018_edge_calibrated_groin` was found holding **M = 70, f = 0.6** —
an unstable cell whose fillet runs away to 444 m — while being named as though
it were the calibrated run. It is now
`output/archive/2026-08-30_rig-M70-unstable/`. Check
`*_groin_diagnostics.csv` (`trapping_rate_applied_m_yr`) before believing any
rig run's label.

---

## One operational rule

**Never run two sweep orchestrators at once.** Every CASCADE construction
writes the shared `data/hatteras_init/Hatteras-CASCADE-parameters.yaml`;
concurrent writers corrupt it, which cost 25 cells of a running sweep on
2026-08-24.
