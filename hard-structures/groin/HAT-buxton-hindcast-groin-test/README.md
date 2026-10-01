# HAT-buxton-hindcast-groin-test - the 1967 rig's sweep results and its three-run figure

The products of the Buxton groin rig (GIS 2-12, 1967-2017), apart from the
runs themselves, which are in `output/calibration/groin_rig/`. The scripts
that make the runs and the sweep are in
`../HAT-groin-buxton-output/1967_2017_run/`.

```
comparison/          HAT_three_run_comparison.py: no groin / nourishment only /
                     nourishment + groin, side by side (its README says how to
                     make the three runs); figures in comparison/three_run_comparison/
sensitivity_sweep/   the (M, f) sweep on the rig, written by
                     ../HAT-groin-buxton-output/1967_2017_run/HAT_groin_sensitivity_sweep.py:
                     HAT_groin_sweep_results.csv (tracked), the heatmap and
                     profile figures, and profiles/<M>_frac<f>.npy per cell
    archive_july_20260824_081742/    an earlier sweep, results and profiles
    archive_pre1984start_20260830/   the sweep as it stood before the 1984 start
```

`GROIN_PLAN.md` cites `sensitivity_sweep/` as the 1967 rig's sweep (under the
older name `HAT-hindcast-groin-test`).
