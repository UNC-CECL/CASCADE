# calibration - pipeline inputs written by a calibration

```
groin/joint_fit.json   the pinned groin (blocking, b 0.6, f 0.6, fitted on 1996-2009,
                       2026-10-05). A pipeline INPUT: HAT_run_all reads it
groin/README.md        what the pin is and where its evidence lives
```

`groin/` paths are built from `HAT_groin_sweep_config.GROIN_SWEEP_ROOT`. The
fit's evidence is `hard-structures/groin/3-hindcast/2-blocking-1996-2025/2026-10-05-blocking-fit-calibration/`
and `raw_runs/experiments/groin/2026-10-05-blocking-fit-dem-to-dem/`. The
current per-domain BE (set 1) is `raw_runs/experiments/source-sink/2026-10-05-be-domain-solve-1996_2009/`.

## Moved off C: on 2026-10-07

Superseded calibrations, now on `D:\CASCADE_offload\output\calibration\`, same
paths (tracked READMEs and DECISION.md are in git history):

- `hs/`: the Hs 3.0 vs 2.5 test on the ÷10 offset (replaced by option A waves).
- `sensitivity/`: the one-forcing-at-a-time sweep's manifests and figures.
- `groin_rig/`: the 1967-2018 dipole rig.
- `groin/` sweep cells, figures, validation runs, `joint_fit.csv`,
  `joint_fit_ranking.json` and `SELECTED_M60_f0.60/` (the dipole M 60 / f 0.6
  pin of 2026-08-30).
