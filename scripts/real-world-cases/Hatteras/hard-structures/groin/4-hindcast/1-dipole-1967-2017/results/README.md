# results — the 1967 rig's sweep results and its three-run figure

The products of the Buxton groin rig (GIS 2-12, 1967-2017), apart from the
runs themselves, which are on `D:\hatteras_ch1\model_runs\output\calibration\groin_rig\`. The scripts
that make the runs and the sweep are in
`../runs/1967_2017_run/`.

```
comparison/          HAT_three_run_comparison.py: no groin / nourishment only /
                     nourishment + groin, side by side (its README says how to
                     make the three runs); figures in comparison/three_run_comparison/
sensitivity_sweep/   the (M, f) sweep on the rig, written by
                     ../runs/1967_2017_run/HAT_groin_sensitivity_sweep.py:
                     HAT_groin_sweep_results.csv (tracked), the heatmap and
                     profile figures, and profiles/<M>_frac<f>.npy per cell
```

Two earlier sweeps, `sensitivity_sweep/archive_july_20260824_081742/` and
`sensitivity_sweep/archive_pre1984start_20260830/`, were moved on 2026-10-09 to
`D:\hatteras_ch1\model_runs\hard-structures\groin\3-hindcast\1-dipole-1967-2017\results\sensitivity_sweep\`
(same folder names; their `HAT_groin_sweep_results.csv` left git then).

`../dipole_fit_notes.md` cites `sensitivity_sweep/` as the 1967 rig's sweep.
