# `output/calibration/groin` - the pinned groin

**PINNED 2026-10-05: blocking groin, b = 0.6, f = 0.6**, fitted on the 1996-2009 calibration period by date
RMSE of the D5-D6 gap at the 1997/2004/2008 photos (4.0 m, against 69.2 m with no groin).

- `joint_fit.json` is a pipeline INPUT: `HAT_run_all` reads it, and it carries the pin and its note.
  M is carried only because the runner reads it; the blocking groin ignores it.
  **Re-running stage 5 (`HAT_groin_joint_fit.py`) overwrites this file.**
- Evidence: `hard-structures/groin/3-hindcast/2-blocking-1996-2025/2026-10-05-blocking-fit-calibration/README.md`
  and the runs in `output/raw_runs/experiments/groin/2026-10-05-blocking-fit-dem-to-dem/`.

Every path is built from `HAT_groin_sweep_config.GROIN_SWEEP_ROOT`.

## Moved off C: on 2026-10-07

The dipole (M, f) sweep and everything from it are on `D:\CASCADE_offload\output\calibration\groin\`,
same paths: the period sweep cells, `fullperiod_1984_2024/`, `_validation_*`, `figures/`, `archive/`
(including the dipole pin `joint_fit_dipole_M60_pin_20260830.json`), `joint_fit.csv`,
`joint_fit_ranking.json`, `SELECTED_M60_f0.60/` and the previous version of this README.
