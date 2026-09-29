# sensitivity — the wave sensitivity sweep on the pre-ceiling matrix (ARCHIVED)

> **Archived 2026-09-28 19:15 with the matrix it was run around** (`../README.md`).
> It ran on Barrier3D `fix/route-overwash-axis-swap`, storms `v3_72` and the
> LOESS-7 option A ends; the adopted setup changed all three. Do not use for
> analysis. The current sweep is in `raw_runs/sensitivity/`.

One wave parameter moved at a time around the matrix runs, everything else as
published: option A waves (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5),
edgeBE on the LOESS-7 ends, full management, no groin. Run 2026-09-28, both
windows, 48 cells, all completed; none drowned the barrier.

## Figures

In `figures/` beside the runs (moved from `output/calibration/sensitivity/figures/`):

| window | summary plot (all four parameters) | per-parameter alongshore | numbers |
|---|---|---|---|
| 1996–2010 | [`01_skill_overview.png`](figures/1996_2010_edgeBE/01_skill_overview.png) | `02`–`05_alongshore_<parameter>.png` beside it | [`summary.csv`](figures/1996_2010_edgeBE/summary.csv) |
| 2010–2024 | [`01_skill_overview.png`](figures/2010_2024_edgeBE/01_skill_overview.png) | `02`–`05_alongshore_<parameter>.png` beside it | [`summary.csv`](figures/2010_2024_edgeBE/summary.csv) |


## Layout

```
<axis>/<window>/edgeBE/<matrix run name>_<token>/
  waveHs    Hs 0.75, 1.0, 1.25, 1.5, 1.75, 2.25, 2.5, 2.75, 3.0 m
  waveTp    Tp 6, 7, 8, 9, 10, 12 s
  waveasym  asymmetry 0.5, 0.55, 0.65, 0.7, 0.8
  waveahf   high-angle fraction 0.3, 0.4, 0.45, 0.55
```

The baseline is the matrix run of the same name without the token, now in
`../matrix/<window>/edgeBE/`. Driver:
`scripts/sensitivity_analysis/hindcast_sensitivity.py --start-year <1996|2010>
--param <axis> --preset edgeBE --no-groin`; manifests
`../manifests/sensitivity_<year>.jsonl`, logs `logs/`.

## What it found

Interior RMSE (m/yr) across each parameter's range; matrix baseline 1.19
(1996–2010) and 2.31 (2010–2024).

| parameter | 1996–2010 | 2010–2024 |
|---|---|---|
| high-angle fraction | 1.19 → 3.45 at 0.3 | 2.31 → 4.74 at 0.3 |
| Hs | 1.27 (0.75) → 1.13 (3.0) | 2.78 (0.75) → 2.18 (3.0) |
| asymmetry | 1.20–1.25 | 2.32–2.47 |
| Tp | 1.19–1.21 | 2.27–2.32 |

- High-angle fraction dominates both windows; below 0.5 the error climbs fast
  and the bias goes negative.
- Higher Hs scores better, but the end rates were solved at Hs 2.0: away from
  it they no longer match, and the 2026-09-27 Hs check showed that gain is the
  mismatch, not the waves (`experiments/wave-climate/2026-09-27-wave-recommendation/`).
- Tp and asymmetry barely move the fit.
- 2010–2024 bias stays −1.5 to −1.8 m/yr at every Hs, Tp and asymmetry value (high-angle below 0.5 only makes it worse): no wave setting fixes that
  window (the 2021 CoastSat step).
- 2010–2024 Hs 1.75 drowned NC-12 at GIS 14 in year 9; Hs 1.5 and 2.0 did not.

Why these runs are here and the 2026-09-24..27 wave studies are not: those
searched for the waves, on other end rates, the LOESS-10 target and tokens
relative to the old defaults. See `experiments/wave-climate/README.md`.
