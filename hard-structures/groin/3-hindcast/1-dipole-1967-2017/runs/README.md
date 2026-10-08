# runs — the 1967-2017 dipole rig

```
1967_2017_run/               the rig runner, the (M, f) sweep and its worker, the
                             edge solve, the fillet-trajectory target, plots
    comparison/              no-groin vs groin against the observed shoreline
1967_1997_run/comparison/    six PNGs left from the deleted 1967-1997 precursors (untracked)
```

Each folder holding scripts has its own README with the reasoning behind
them. The rig's runs were written to `output/calibration/groin_rig/`, now on
`D:\CASCADE_offload\output\calibration\groin_rig\` (2026-10-07). The observed target every
script here reads is `../../../1-observations/wetdry_photo_positions/Change_from_wetdry_1967_D2_D12.csv`.

## Deleted 2026-10-01: the 1967-1997 precursors

`1967_1997_run/` and `1967_1997_no_BE_run/` held five scripts from the
original 30-year (1967-1997) groin-only test, which the 1967-2017 rig replaced:

| deleted | what it was |
|---|---|
| `1967_1997_run/HAT_groin_hindcast_1967_1997.py` | the original 30-year rig runner; both 1967-2017 runners were built from it |
| `1967_1997_no_BE_run/HAT_groin_hindcast_1967_1997_noBE.py` | the 1967-1997 no_groin / groin pair, its no-BE variant |
| `1967_1997_run/HAT_plot_groin_msweep.py` | figures of the 1967-1997 M sweep |
| `1967_1997_run/HAT_plot_groin_runs.py` | the 1967-1997 copy of the run plotter |
| `1967_1997_run/comparison/HAT_groin_effect_comparison.py` | the 1967-1997 copy of the effect comparison |

They were deleted because `GROIN_PLAN.md` and the fit never cite them, and
nothing imports them; other scripts mention them only as where they were
adapted from. Retired code is deleted, not parked (ORGANIZATION.md rule 4).
Their figures were git-ignored and were kept: they are still in
`1967_1997_run/comparison/` (`base_vs_60M/`, `M_sweep_noBE/`).

To recover any of them:

```
git log --diff-filter=D --oneline -- hard-structures/groin/3-hindcast/1-dipole-1967-2017/runs/<path>
git show <commit>^:hard-structures/groin/3-hindcast/1-dipole-1967-2017/runs/<path>
```
