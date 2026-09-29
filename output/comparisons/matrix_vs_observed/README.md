# matrix_vs_observed: every option A matrix run (adopted model) against CoastSat and the dune line

Two figures per run, made 2026-09-27 (Hannah).

**Redrawn 2026-09-28 on the adopted model**: the whole matrix, zeroBE included, was rebuilt on Barrier3D `hatteras/adopted` (the three overwash fixes and per-cell dune ceilings), storms `v3_trim24` (every event kept, trimmed to 24 h around its peak), with ends re-solved on it (1996 +4.3509 / +19.0935, 2010 +8.0 / +22.4937 m/yr; `experiments/end-domain-boundaries/2026-09-28-ends-resolved-adopted/`). The pre-adoption matrix is in `output/raw_runs/archive/2026-09-28-pre-ceiling/`. Against it, 2010–2024 is less erosive in every scenario (natural position bias −50.7 → −25.8 m, full management −18.6 → −14.8 m); 1996–2010 moves by about 1 m.

**Earlier the same day** the tree was redrawn on the edgeBE matrix re-run at the ends re-solved against the LOESS-7 target (1996 +4.8394 / +18.2545, 2010 +18.8657 / +24.2358 m/yr), with the target and the observed change smoothed at 7 domains. Against the LOESS-10 version, rate bias moved by ≤ 0.02 m/yr and rate RMSE rose by 0.04–0.14 m/yr in both presets (the sharper target, not the ends). The previous edgeBE runs are in `output/raw_runs/archive/2026-09-28-loess10-ends/`. The runner draws only rate figures, so neither existed for the matrix.

The runs are the option A no-groin matrix on the adopted model:
- Barrier3D `hatteras/adopted` (the three overwash fixes and per-cell dune ceilings), storms `v3_trim24` (every event kept, trimmed to 24 h around its peak);
- metres offset;
- Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5;
- 1996–2010 and 2010–2024;
- edgeBE and zeroBE, every scenario.

Each figure also sits in its run's own `figures/vs_observed/`.

## The two figures

**`rate_and_position_change/`** ("Where are these position plots?")
- (a) The model's OLS rate against the CoastSat LRR scoring target (7-domain LOESS, raw means GIS 1–10; 10 until 2026-09-28).
- (b) The model's position change (endpoint rate × 14 yr) against the observed CoastSat change: the mean position in the last calendar year minus the first, smoothed at 10 domains.
- `scenarios_rate_and_position_<preset>_<window>.png` puts every scenario on one pair of panels. Start there.

**`start_and_end_positions/`** ("the starting island position with the end modeled position and the observed end position", with both CoastSat and the dune line)
- **Positions are drawn relative to the start line** (Hannah's choice). The island's own position swings about 6 km along the reach while the changes are tens of metres, so on an absolute axis the lines would lie on top of each other.
- The model's start position is the zero line. The other lines are:
  - black: the modelled end position;
  - blue: CoastSat, the observed change added to the start (domain means, with the 7-domain LOESS faint);
  - red: the dune line, the change between the start and end vintages (1997 → 2009 for 1996–2010, 2009 → 2023 for 2010–2024). This is the runner's own end-year target. Its survey interval (11.6 and 14.1 yr) is not the calendar window and is not rescaled.

Seaward is positive. Scores cover the interior, GIS 2–89.

## Layout

```
rate_and_position_change/<window>/scenarios_rate_and_position_<preset>_<window>.png
rate_and_position_change/<window>/<preset>/rate_and_position_<preset>_<scenario>[_reloc]_<window>.png
start_and_end_positions/<window>/<preset>/start_and_end_positions_<preset>_<scenario>[_reloc]_<window>.png
scores.csv      per run: rate bias/RMSE (m/yr); position-change bias/RMSE (m);
                end-position bias/RMSE against CoastSat and against the dune line (m)
y_bounds.txt    the fixed y ranges
```

## One y axis per panel type

The same range is used across every run, so any two figures compare directly:

| panel | range |
|---|---|
| rate | −7.5 to +7.5 m/yr, as each run's `shoreline_change_rate.png`. The run's `_with_buffers` figure is ±10. |
| position change | −130 to +130 m |
| start and end positions | −130 to +130 m |

- The runner's own rate figures were redrawn with `rerender_run_figures.py --arm matrix --ylim=-10,10 --ylim-real=-7.5,7.5`.
- Every model line, CoastSat LOESS and domain mean fits inside ±7.5 m/yr.
- 8 of about 1,800 CoastSat transect dots fall outside ±7.5 and are cut off. The `_with_buffers` figure shows them.

## What the end positions show (edgeBE, full management)

| window | model vs CoastSat end | model vs dune-line end |
|---|---|---|
| 1996–2010 | bias +4.8 m, RMSE 25.1 | bias +10.6 m, RMSE 23.1 |
| 2010–2024 | bias −13.9 m, RMSE 28.2 | bias −5.0 m, RMSE 22.8 |

Before adoption: +5.7 / 25.7 and +11.6 / 23.4 m in 1996–2010, −17.7 / 30.2 and −8.8 / 23.2 m in 2010–2024.

In 2010–2024 the model ends closer to the dune line than to CoastSat, the same reading as `model_vs_observed/`.

Drawn by `scripts/analyze_output/compare_runs/matrix_vs_observed.py`. It was named `matrix_rate_and_position.py`, writing to `matrix_rate_and_position/`, for a few hours on 2026-09-27.
