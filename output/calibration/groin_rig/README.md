# calibration/groin_rig - the 1967 to 2018 groin rig (was output/rig_runs/ until 2026-09-18)

The only window that spans the Buxton groin's whole life: installation, the
1996 repair, the 2003 storm damage, and the decline after it. That is why it
exists - no hindcast period contains the deterioration ramp end to end.

```
HAT_1967_2018_edge_calibrated_groin/      with the groin
HAT_1967_2018_edge_calibrated_no_groin/   the paired baseline
```

The unstable M = 70 cell once left here under the calibrated run's name is in
`output/archive/2026-08-30_rig-M70-unstable/`.

**The rig brackets the deterioration fraction and only rails on the trapping
rate.** It cannot distinguish M = 60 from any larger value, because the solver
goes unstable above it. Quote the pair as fitted on the production windows, and
cite the rig for f alone. `hard-structures/groin/GROIN_PLAN.md` is the
authority.

The rig runs 1967 off a 1984 island: a deliberate seventeen-year anachronism in
the initial condition, accepted because the target is a shoreline differential.
