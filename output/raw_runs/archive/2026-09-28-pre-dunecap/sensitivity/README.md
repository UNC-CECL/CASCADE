# sensitivity — the matrix wave-climate sensitivity record

One wave parameter moved at a time around the edgeBE full-management matrix
runs, everything else as published. Values (set with Hannah 2026-09-28; the
matrix value in bold is not re-run, it is the baseline):

| axis folder | parameter | values |
|---|---|---|
| `waveHs` | Hs (m) | 0.75, 1.0, 1.25, 1.5, 1.75, **2.0**, 2.25, 2.5, 2.75, 3.0 |
| `waveTp` | Tp (s) | 6, 7, **7.5**, 8, 9, 10, 12 |
| `waveasym` | asymmetry | 0.5, 0.55, **0.6**, 0.65, 0.7, 0.8 |
| `waveahf` | high-angle fraction | 0.3, 0.4, 0.45, **0.5**, 0.55 |

Both windows, 24 cells each. Layout `<axis>/<window>/edgeBE/<matrix run name>_<token>/`;
figures in `figures/<window>_edgeBE/` (`01_skill_overview.png` is the summary of
all four parameters), written by `scripts/sensitivity_analysis/plot_sensitivity.py`.

## Status: done (2026-09-28, 20:03–21:43)

Run on the adopted setup (Barrier3D `hatteras/adopted`, storms `v3_trim24`,
ends 1996 +4.3509 / +19.0935, 2010 +8.0 / +22.4937). 48 of 48 cells completed;
none drowned the barrier. Launched by a queue that waited for both edgeBE
baselines (`output/calibration/sensitivity/logs/queue_adopted_20260928.status`);
driver logs `sweep_<year>_edgeBE_adopted_20260928.log` beside it; manifests
`output/calibration/sensitivity/sensitivity_<year>.jsonl`.

The first run of this sweep (16:24–18:57, on the setup before the dune ceiling
and storm adoption) is in `archive/2026-09-28-pre-ceiling/sensitivity/`.

## What it found

Interior RMSE (m/yr) across each parameter's range; matrix baseline 1.168
(1996–2010, bias −0.03) and 2.123 (2010–2024, bias −1.45).

| parameter | 1996–2010 | 2010–2024 |
|---|---|---|
| high-angle fraction | 1.28 (0.55) → 3.63 (0.3) | 2.22 (0.55) → 4.88 (0.3) |
| Hs | 1.25 (0.75) → 1.13 (3.0) | 2.43 (0.75) → 2.07 (3.0) |
| asymmetry | 1.19–1.24 | 2.16–2.30 |
| Tp | 1.17–1.19 | 2.10–2.15 |

- High-angle fraction is the only parameter the fit depends on strongly, in
  both windows: below 0.5 the error climbs fast and the bias goes negative.
- Higher Hs scores slightly better, but the ends were solved at Hs 2.0; away
  from it they no longer match, and the 2026-09-27 Hs check showed that gain is
  the mismatch (`experiments/wave-climate/2026-09-27-wave-recommendation/`).
- Asymmetry has its minimum at 0.6 in both windows; Tp barely matters.
- 2010–2024 bias stays −1.3 to −1.6 m/yr at every Hs, Tp and asymmetry value:
  no wave setting fixes that window (the 2021 CoastSat step).
- Every 2010–2024 cell drowns NC-12 at GIS 14 in year 9, and so does the matrix
  baseline: a property of the adopted setup in that window, not of the waves.
  1996–2010 drowns none.
