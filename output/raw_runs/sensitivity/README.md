# sensitivity: the matrix wave-climate sensitivity record

One wave parameter is moved at a time around the edgeBE full-management matrix runs. Everything else is as published. The values were set with Hannah on 2026-09-28. The matrix value (bold) is not re-run; it is the baseline.

| axis folder | parameter | values |
|---|---|---|
| `waveHs` | Hs (m) | 0.75, 1.0, 1.25, 1.5, 1.75, **2.0**, 2.25, 2.5, 2.75, 3.0 |
| `waveTp` | Tp (s) | 6, 7, **7.5**, 8, 9, 10, 12 |
| `waveasym` | asymmetry | 0.5, 0.55, **0.6**, 0.65, 0.7, 0.8 |
| `waveahf` | high-angle fraction | 0.3, 0.4, 0.45, **0.5**, 0.55 |

Both windows, 24 cells each.
- Layout: `<axis>/<window>/edgeBE/<matrix run name>_<token>/`.
- Figures: `figures/<window>_edgeBE/`, written by `scripts/sensitivity_analysis/plot_sensitivity.py`. `01_skill_overview.png` summarises all four parameters, and `summary.csv` holds the numbers.

## Status: done (2026-09-28 23:12 to 09-29 00:19)

Run on the current setup:
- Barrier3D `hatteras/adopted`, storms `v3_trim24`, option A waves.
- The beach/dune manager's 4 m cap limited to the sand it adds.
- Ends 1996 +4.3509 / +19.0935, 2010 +8.0 / +21.2582.

48 of 48 cells completed. Hannah asked for the re-run after the dune-cap fix ("yes, re-run the wave sweep after the matrix"). The queue checked the baselines (Barrier3D branch, cap mode, ends equal to the config) before starting; its status file is `output/calibration/sensitivity/logs/queue_dunecap_20260928.status`. Driver logs are `sweep_<year>_edgeBE_dunecap_20260928.log` beside it, and manifests are `output/calibration/sensitivity/sensitivity_<year>.jsonl`.

Earlier runs of this sweep, both superseded:
- `archive/2026-09-28-pre-ceiling/sensitivity/` (16:24-18:57, before the dune ceiling and storm adoption)
- `archive/2026-09-28-pre-dunecap/sensitivity/` (20:03-21:43, adopted setup but the cap still clipped whole dune cells)

## What it found

Interior RMSE (m/yr) across each parameter's range. The matrix baseline is 1.172 in 1996-2010 (bias +0.06) and 2.068 in 2010-2024 (bias −1.36). The previous run is in brackets.

| parameter | 1996-2010 | 2010-2024 |
|---|---|---|
| high-angle fraction | 1.30 (0.55) → 3.60 (0.3) [1.28 → 3.63] | 2.17 (0.55) → 4.84 (0.3) [2.22 → 4.88] |
| Hs | 1.25 (0.75) → 1.13 (3.0) [1.25 → 1.13] | 2.36 (0.75) → 2.02 (2.75) [2.43 → 2.07] |
| asymmetry | 1.18-1.23, lowest at 0.6 [1.19-1.24] | 2.10-2.25, lowest at 0.6 [2.16-2.30] |
| Tp | 1.17-1.20 [1.17-1.19] | 2.04-2.10 [2.10-2.15] |

- **The cap fix did not change the picture.** Every curve has the same shape as in the previous run. 2010-2024 sits about 0.05 lower throughout, which is the baseline improvement.
- **High-angle fraction is the only parameter the fit strongly depends on**, in both windows. Below 0.5 the error climbs fast and the bias goes negative. The option A value (0.5) is the minimum in both windows; 0.55 is slightly worse.
- **Asymmetry has its minimum at 0.6** in both windows.
- **Higher Hs scores slightly better** (≤ 0.05 m/yr at 2.75-3.0). The ends were solved at Hs 2.0 and no longer match away from it; the 2026-09-27 Hs check showed that this gain is the mismatch (`experiments/wave-climate/2026-09-27-wave-recommendation/`).
- **Tp barely matters.** In 2010-2024, Tp 6 is 0.03 better than 7.5, within the same ends-mismatch margin.
- **No wave setting fixes the 2010-2024 bias**, which stays at −1.2 to −1.5 m/yr at every Hs, Tp and asymmetry value. The cause is the 2021 CoastSat step.
- **Every 2010-2024 cell drowns one NC-12 domain**, and so does the matrix baseline (`roads_drowned` = 1). This is a property of the setup in that window, not of the waves. 1996-2010 drowns none.
