# sensitivity: the matrix wave-climate sensitivity record

One wave parameter is moved at a time around the edgeBE full-management matrix runs. Everything else is as published. The values were set with Hannah on 2026-09-28. The matrix value (bold) is not re-run; it is the baseline.

| axis folder | parameter | values |
|---|---|---|
| `waveHs` | Hs (m) | 0.75, 1.0, 1.25, 1.5, 1.75, **2.0**, 2.25, 2.5, 2.75, 3.0 |
| `waveTp` | Tp (s) | 6, 7, **7.5**, 8, 9, 10, 12 |
| `waveasym` | asymmetry | 0.5, 0.55, **0.6**, 0.65, 0.7, 0.8 |
| `waveahf` | high-angle fraction | 0.3, 0.4, 0.45, **0.5**, 0.51, 0.52, 0.53, 0.54, 0.55 |

Both windows, 28 cells each: the 24 of the main sweep plus the fine high-angle cells 0.51-0.54 (see below).
- Layout: `<axis>/<window>/edgeBE/<matrix run name>_<token>/`.
- Figures: `figures/<window>_edgeBE/`, written by `scripts/sensitivity_analysis/plot_sensitivity.py`. `01_skill_overview.png` summarises all four parameters, and `summary.csv` holds the numbers. `06_best_setting_optionA.png` is the option A matrix baseline on its own (the best setting; the pink line in `02`-`05`). It is a copy of that matrix run's `figures/shoreline_change_rate.png`, redrawn with `rerender_run_figures.py --loess-only`, so re-copy it if the matrix run changes. `figures/position_change/{total,projected}/<window>_edgeBE/` holds the same 02-06 as shoreline POSITION change (m): each model run's end-minus-start change against CoastSat LRR x 14 yr, where `total` is the window's own LRR and `projected` the 1996-2024 LRR. Drawn by `plot_sensitivity.py --quantity position --reference total|projected`; the per-run versions are in each run's `figures/position_change/<reference>/` from `rerender_run_figures.py --loess-only --position-change`. No 01 there: the skill numbers are rate RMSE.

## Status: done (2026-09-28 23:12 to 09-29 00:19)

Run on the current setup:
- Barrier3D `hatteras/adopted`, storms `v3_trim24`, option A waves.
- The beach/dune manager's 4 m cap limited to the sand it adds.
- Ends 1996 +4.3509 / +19.0935, 2010 +8.0 / +21.2582.

48 of 48 cells completed. Hannah asked for the re-run after the dune-cap fix ("yes, re-run the wave sweep after the matrix"). The queue checked the baselines (Barrier3D branch, cap mode, ends equal to the config) before starting; its status file is `output/calibration/sensitivity/logs/queue_dunecap_20260928.status`. Driver logs are `sweep_<year>_edgeBE_dunecap_20260928.log` beside it, and manifests are `output/calibration/sensitivity/sensitivity_<year>.jsonl`.

Earlier runs of this sweep, both superseded:
- `archive/2026-09-28-pre-ceiling/sensitivity/` (16:24-18:57, before the dune ceiling and storm adoption)
- `archive/2026-09-28-pre-dunecap/sensitivity/` (20:03-21:43, adopted setup but the cap still clipped whole dune cells)

## Fine high-angle sweep (2026-09-29, done 09:04)

Four extra high-angle cells per window (0.51, 0.52, 0.53, 0.54), run to see whether the minimum sits just above 0.5. Status: `output/calibration/sensitivity/logs/ahffine_20260929.status`. Logs: `sweep_<year>_edgeBE_ahffine_20260929.log`. The cells went into the same manifests and figures as the main sweep.

Interior RMSE (m/yr), with mean bias in brackets:

| high-angle fraction | 1996-2010 | 2010-2024 |
|---|---|---|
| **0.50 (baseline)** | **1.172 (+0.06)** | **2.068 (−1.36)** |
| 0.51 | 1.231 (+0.13) | 2.091 (−1.33) |
| 0.52 | 1.269 (+0.11) | 2.133 (−1.36) |
| 0.53 | 1.283 (+0.10) | 2.149 (−1.38) |
| 0.54 | 1.298 (+0.08) | 2.159 (−1.38) |
| 0.55 | 1.298 (+0.09) | 2.167 (−1.38) |

- **0.5 is still the minimum in both windows.** No fine cell beats it, so option A stands.
- **The error rises steeply just above 0.5.** In 1996-2010, going from 0.50 to 0.51 costs +0.06 m/yr, which is half the whole 0.50 → 0.55 rise. The minimum is sharp on that side, so the baseline sits at a narrow bottom rather than a flat one.
- **Bias doesn't explain the rise.** In 2010-2024 it stays at −1.33 to −1.38 across the fine cells (the 2021 CoastSat step). In 1996-2010 it moves only +0.13 → +0.08.
- **Roads are unchanged.** Every 2010-2024 cell still drowns one NC-12 domain, like the baseline; 1996-2010 drowns none.

## What it found

Interior RMSE (m/yr) across each parameter's range. The matrix baseline is 1.172 in 1996-2010 (bias +0.06) and 2.068 in 2010-2024 (bias −1.36). The previous run is in brackets.

| parameter | 1996-2010 | 2010-2024 |
|---|---|---|
| high-angle fraction | 1.30 (0.55) → 3.60 (0.3) [1.28 → 3.63] | 2.17 (0.55) → 4.84 (0.3) [2.22 → 4.88] |
| Hs | 1.25 (0.75) → 1.13 (3.0) [1.25 → 1.13] | 2.36 (0.75) → 2.02 (2.75) [2.43 → 2.07] |
| asymmetry | 1.18-1.23, lowest at 0.6 [1.19-1.24] | 2.10-2.25, lowest at 0.6 [2.16-2.30] |
| Tp | 1.17-1.20 [1.17-1.19] | 2.04-2.10 [2.10-2.15] |

- **The cap fix did not change the picture.** Every curve has the same shape as in the previous run. 2010-2024 sits about 0.05 lower throughout, which is the baseline improvement.
- **High-angle fraction is the only parameter the fit strongly depends on**, in both windows. Below 0.5 the error climbs fast and the bias goes negative. The option A value (0.5) is the minimum in both windows; 0.55 is slightly worse, and the fine cells 0.51-0.54 confirm it (see above).
- **Asymmetry has its minimum at 0.6** in both windows.
- **Higher Hs scores slightly better** (≤ 0.05 m/yr at 2.75-3.0). The ends were solved at Hs 2.0 and no longer match away from it; the 2026-09-27 Hs check showed that this gain is the mismatch (`experiments/wave-climate/2026-09-27-wave-recommendation/`).
- **Tp barely matters.** In 2010-2024, Tp 6 is 0.03 better than 7.5, within the same ends-mismatch margin.
- **No wave setting fixes the 2010-2024 bias**, which stays at −1.2 to −1.5 m/yr at every Hs, Tp and asymmetry value. The cause is the 2021 CoastSat step.
- **Every 2010-2024 cell drowns one NC-12 domain**, and so does the matrix baseline (`roads_drowned` = 1). This is a property of the setup in that window, not of the waves. 1996-2010 drowns none.
