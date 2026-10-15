# Groin under option A: where we left off (2026-09-29)

The full record is in `README.md` in this folder, in the order it happened. This page is the
short version: what's settled, what's waiting on a decision, and what's loose.

## Settled

- **M = 60 / f = 0.6 is dead under option A.** It was fitted under the ÷10 offset, whose 25 m cape
  step hid the problem. On the metres planform it builds a ~216 m fillet.
- **Why no groin fitted both windows.** BRIE's alongshore diffusion flattens the real 256 m step at
  GIS 5/6 (the cape turn), and the old 1996→2003 linear wear-down put a decline inside the 1996 window
  that the data don't show. The observed gap holds at 134–155 m through 2004 and falls after it.
- **In the code (uncommitted):**
  - failure schedule is now **instant from the 2004 step** (runner `.py` + notebook);
  - `BlockingGroinCallback` in `cascade/groin.py`, selected by `groin.kind` / `HAT_GROIN_KIND`,
    with `groin.blocking_b` / `HAT_GROIN_BLOCKING_FRACTION`;
  - runner reads r_ipl at the groin cell's own angle (the 0° value crashed every groin run under
    option A);
  - 11 tests in `tests/test_groin_blocking.py`; dipole bit-identical to before; blocking b = 0
    bit-identical to the matrix no-groin run.
- **Full-model grid** (100 runs, `experiments/groin/2026-09-29-instant-2004-grid/`, scores in
  `instant_grid_scores.csv`), date-RMSE vs the observed gap:

  | | 1996 | 2010 | joint |
  |---|---|---|---|
  | no groin | 65.0 | 30.2 | 50.7 |
  | blocking b 0.6, f 0.4 | 9.1 | 15.8 | 12.9 |
  | dipole M 12, f 0.3 | 2.4 | 17.4 | 12.4 |

  A tie within noise. Blocking moves ~42,000 m³/yr (7% of budget) against the dipole's ~113,000 (19%).

## Waiting on Hannah

1. **Pin the groin.** Recommended: blocking b = 0.6, f = 0.4 as the default, with the dipole
   (M 12, f 0.3) reported as the comparison. f is ONE field shared by both kinds, so pinning moves
   the default f 0.6 → 0.4 for the dipole too, and the dipole's default M (60) is stale either way.
   Pinning touches:
   - `hat_run.yaml` (`groin.kind`, `blocking_b`, `deterioration_f`, `trapping_M`);
   - the config defaults in `HAT_hindcast_config.py`;
   - `output/calibration/groin/joint_fit.json`, which still says M 60 / f 0.6, and which
     `HAT_run_all.py` stage 6 passes to every groin matrix run. It has no blocking field, so it
     needs a new shape or the driver needs to read the config.
2. **Commit.** Nothing from this work is committed. A few study files show as staged ("A"); I didn't stage them, so another
   session may have.

## Loose ends, not started

- **The 2010 late drop.** Both groins miss the 2018–23 fall: −28 / −21 m against −48 m observed at 2023.
  Candidates: the 2018 storms, the Buxton 2017/2022 fills, or the CoastSat 2021 step (a different
  dataset, but worth checking).
- **Groin-on matrix runs** once a kind and values are pinned.
- **`HAT_groin_sweep_worker.py` is stale:** it hardcodes the old waves (2.5/8/0.7/0.1) and restores a
  parameter snapshot with no per-cell dune ceilings. Retire it or bring it up to date before anything
  uses `HAT_run_all.py` stages 3–5 again.
- **`output/calibration/groin/README.md`** still presents M 60 as the answer (a pointer to this page
  has been added at its top).
- **The old grids here are on superseded settings.** `2026-09-29-option-a-grid` (106 runs) and the
  emulator scans used the old linear ramp, which the code no longer reproduces. They are kept as the
  record of why the schedule changed.
- **The cape itself.** BRIE has no process that keeps it, so at b = 0.6 the groin partly stands in for
  the cape. That caveat belongs in the writeup whichever kind is pinned.
