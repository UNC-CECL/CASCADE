# archive - retired output, one folder per retirement

**Do not use anything in here for analysis.** Everything that has been
superseded goes here, and only here, as `YYYY-MM-DD_<what>/`. Each folder holds a
`README.md` or `WHY.md` saying what replaced it. Started 2026-09-18, when the
retired material scattered around `output/` was collected in one place.

```
2026-08-28_full-tree/            the whole output/ tree before the re-measured-forcing restart   (was output/superseded_20260828/)
2026-08-30_groin-railed-ranking/ the joint-fit ranking before M/f were pinned                    (was output/groin_sweep/archive/)
2026-08-30_rig-M70-unstable/     the unstable M = 70 rig cell left under the calibrated name     (was output/rig_runs/HAT_1967_2018_rig_M70_f0.6_UNSTABLE_sweep_leftover/)
2026-09-17_figures/              4 figures retired when output/figures/ was sorted by subject    (was output/figures/superseded_20260917/)
2026-09-18_hindcast-calibrated/  the last calibBE render on the 1984/2004 chain                   (was output/comparisons/hindcast_calibrated/superseded_20260918/)
2026-09-18_rate-windows/         the rate_windows figures before the rename to model_vs_observed (was output/comparisons/rate_windows/)
```

**One exception: archived model runs stay in `output/raw_runs/archive/`.**
`archive` is a run kind in `cascade_pipeline.run_registry`, so those runs stay
indexed in `run_index.csv` and still resolve by name. Retired *runs* go there.
Everything else goes here.

Nothing reads this folder. The only references to it are comments that record
where a retired file went.
