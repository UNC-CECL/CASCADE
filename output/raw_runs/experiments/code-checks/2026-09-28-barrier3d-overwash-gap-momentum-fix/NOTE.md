# Barrier3D overwash fixes: gap cells and inundation momentum (2026-09-28)

**Question.** How much do three overwash fixes change the hindcast?

**The fixes** are on Barrier3D branch `fix/overwash-gaps-momentum` (2.0.2.dev1, commit db0ba30). They are local only and have not been pushed (Hannah, 2026-09-28). The storm replay found all three (`scripts/figure_making/model/storm_replay.py`, `output/figures/model/storm_routing_*.png`):

| commit | defect | since |
|---|---|---|
| 990c3bd | `DuneGaps` dropped the last overtopped cell of the last gap, and any single-cell gap | upstream, 2020 |
| 015f11e | gap discharge was set on `start:stop`, but `stop` is inclusive, so each gap lost a cell of water and one-cell gaps got none | upstream, 2020 |
| e929e65 | the inundation momentum constant `C` was reset to 0 before routing | upstream b11b880, the 2024 Numba refactor |

**How the fixed code is run.** The main Barrier3D checkout stays on `fix/route-overwash-axis-swap` (49fd069), the code every matrix run uses. The fixes live in a git worktree at `../Barrier3D-overwashfix`, which each fixed run puts first on `PYTHONPATH`. The hindcast runner, notebook and config are unchanged. Each run's metadata records the Barrier3D it imported, and the driver refuses any run that did not use the fix branch.

**Controls.** The matrix runs of 2026-09-27 (edgeBE, option A waves). Re-running the natural 1996 run on today's unchanged code reproduced its shoreline matrix to 0.0 m.

Driver: `scripts/hatteras_ms/experiments/HAT_barrier3d_gap_momentum_fix.py` (`run`, `compare`).

## Result

| run | net change, fixed − control (mean / max) | overwash domain-years | observed hit rate | interior RMSE (m/yr) |
|---|---|---|---|---|
| natural 1996–2010 | −0.61 m / 10.9 m | 419 → 480 | 16% → 22% | 1.020 → 1.021 |
| managed 1996–2010 | −0.14 m / 1.9 m | 290 → 315 | 13% → 15% | 1.049 → 1.048 |
| natural 2010–2024 | −2.58 m / 20.1 m | 1048 → 1071 | 87% → 88% | 4.14 → 4.37 |
| managed 2010–2024 | −0.38 m / 2.4 m | 655 → 701 | 78% → 80% | 2.25 → 2.28 |

Negative means more landward retreat. The observed hit rate uses the rule in `8-overwash-analysis/4-vs-model`.

- **Direction.** The fixes add overwash (2–15% more domain-years) and so a little more retreat. They are consistent with the storm-strength response figure: they matter most for marginal storms.
- **Size.** The effect is small in the managed runs and moderate in the natural ones. The largest local change is in the natural 2010–2024 run, where Buxton (GIS 11–14) retreats up to 20 m more.
- **Observed overwash.** Agreement with the imagery improves slightly in every run. The 1996–2010 hit rate is still low for another reason: the storm series is missing Isabel. See `../../storms-and-overwash/2026-09-28-storm-max-duration/`.
- **Skill.** The skill against CoastSat is unchanged in 1996–2010 and slightly worse in 2010–2024. That is expected: the edge rates were solved on the old code, so adopting the fixes would mean re-solving the ends.

**Status: record.** The fixes are not adopted. Adopting them means checking the branch out (or merging it) in the Barrier3D repository, re-solving the edge rates, and re-running the matrix. That is Hannah's call.

Files: `tables/comparison.csv`, `tables/per_domain.csv`, `figures/fix_effect.png` (caption in `figures/supporting/CAPTIONS.md`), `logs/`, `runs/fixed_<member>/`.
