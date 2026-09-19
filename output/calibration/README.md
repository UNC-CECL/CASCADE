# calibration - how the forcing and structure parameters were chosen

Collected here on 2026-09-18. Each folder holds a decision and the evidence
behind it. None of them is part of the production matrix.

```
groin/         the (M, f) sweep, joint_fit.json (a pipeline INPUT) and      (was output/groin_sweep/)
               SELECTED_M60_f0.60/README.md, the M/f decision record
hs/            the Hs 3.0 vs 2.5 test. DECISION.md: keep Hs = 2.5        (was output/hs_experiment/)
sensitivity/   the one-forcing-at-a-time parameter sweep: manifests,      (was output/sensitivity_analysis/)
               logs and figures. figures/README.md is the decision record.
               The sweep's RUNS are in output/raw_runs/sensitivity/.
groin_rig/     the 1967-2018 groin rig, the only window spanning the      (was output/rig_runs/)
               deterioration ramp. It brackets f; it does not fit M
```

`groin/` paths are built from `HAT_groin_sweep_config.GROIN_SWEEP_ROOT`. Two callers
outside `scripts/hatteras_ms/groin-sweep/` spell the path out instead:
`HAT_be_zone_residual_fit.py` and `hard-structures/groin/HAT-groin-figures/HAT_groin_module_logic_figure.py`.
The groin study itself (code, inputs, `GROIN_PLAN.md`, the authority on the fit)
is `hard-structures/groin/`. What is here is its model output.
