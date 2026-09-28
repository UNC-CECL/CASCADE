# model_vs_observed: the option A matrix against the shoreline and the dune line

> **Units: this tree is in m/yr.** The other comparison trees (`target_comparison/`,
> `4-comparisons/shoreline_vs_duneline/`, `3-rates/coastsat/{total_change,projected}/`)
> are in **metres**. [`FIGURES.md`](../../../FIGURES.md) indexes all of them.

**Redrawn 2026-09-27 on the option A matrix** (Hannah: "make the model vs observed figures for the
new matrix runs"). The model line comes from these runs:

- metres island offset;
- option A wave climate: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5;
- end rates solved for those waves against the LOESS-7 target (re-solved and the edgeBE matrix re-run 2026-09-28): 1996 GIS 1 +4.8394 / GIS 90 +18.2545, 2010 +18.8657 / +24.2358 m/yr (LOESS-10 values were +17.545 and +18.8 / +24.535; those runs are in `raw_runs/archive/2026-09-28-loess10-ends/`);
- edgeBE, full management, groin off;
- Barrier3D with the route_overwash fix.

| window | run |
|---|---|
| 1996–2010 | `matrix/1996_2010/edgeBE/HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin` |
| 2010–2024 | `matrix/2010_2024/edgeBE/HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin` (nourishment on) |
| 1984–2004, 2004–2024 | no metres run yet; each panel draws the observation and says "model not yet run" |

**Each observation is drawn with the run solved on it**, so each comparison is fair at GIS 1 and 90 by construction:
- The shoreline (CoastSat) figures draw the matrix run above, with ends solved on CoastSat.
- The dune-line figures draw the run with ends solved on the dune line under option A: `experiments/end-domain-boundaries/2026-09-27-ends-solved-on-duneline-option-a/`, mean3 reading. Its ends are 1996 −3.0 / +7.6 and 2010 +3.4 / +15.0 m/yr.
- The `sensitivity/` levels are drawn again:
  - `ends-swapped/`: each observation against the other solve;
  - `dune-raw-solve/`: the raw reading;
  - `mixed-estimator/`.
- Between the first option A redraw and the dune-line re-solve (both on 2026-09-27), the dune-line figures briefly drew the CoastSat-solved run. The switch is `DUNE_SOLVE_CURRENT` in the script, now True.

The previous tree (/10 offset, Hs 2.5, all four windows, the dune-solved runs) is in
`output/archive/2026-09-27_model-vs-observed-div10/`.

## Which figure answers which question

Start with `vs_shoreline/smoothed/model_vs_shoreline_smoothed_grid.png`: the scoring target.

| question | folder | observation | model line |
|---|---|---|---|
| Does the model follow the shoreline? | `vs_shoreline/domain_means/` | CoastSat per-domain mean LRR, ±1 std | OLS rate |
| ...against the smoothed target the run index scores? | `vs_shoreline/smoothed/` | the scoring target (10-domain LOESS, raw means D1–10) as the fill, the means as dots | OLS rate |
| Does the model reproduce the dune line's net change? | `vs_duneline/endpoint_net_change/` | two dune-line surveys differenced per domain | endpoint rate |
| ...smoothed like the CoastSat target? | `vs_duneline/net_change_smoothed/` | the same, 10-domain LOESS | endpoint rate |
| Both observations as net change at the same two dates | `vs_shoreline_and_duneline/` | CoastSat blue, dune line red, both LOESS | endpoint rate |
| Endpoint observation against the model's OLS rate | `sensitivity/mixed-estimator/` | dune-line net change | OLS rate |

Each folder holds one figure per window (`<stem>_<start>_<end>.png`) plus a 2 × 2 `_grid`.
PDFs and `CAPTIONS.md` are under `supporting/`.

## Skill, GIS 2–89 (`tables/skill.csv`)

| window | target | model estimator | bias (m/yr) | RMSE (m/yr) |
|---|---|---|---|---|
| 1996–2010 | CoastSat LRR, smoothed (the scoring target) | OLS | +0.10 | 1.19 |
| 1996–2010 | dune line net change, smoothed (dune-solved run) | endpoint | +0.99 | 1.73 |
| 2010–2024 | CoastSat LRR, smoothed (the scoring target) | OLS | −1.66 | 2.31 |
| 2010–2024 | dune line net change, smoothed (dune-solved run) | endpoint | −0.67 | 1.51 |

- **Smoothed at 7 domains since 2026-09-28**, following the runner's scoring target.
- At 10 domains (`output/archive/2026-09-28_option-a-loess10-comparisons/`) the four RMSEs were 1.05, 1.59, 2.25 and 1.43.
- The narrower window leaves more alongshore detail in the target for the model to miss, so RMSE rises. Bias barely moves.

- At 10 domains the two scoring-target rows reproduced the matrix and the wave-recommendation numbers (RMSE 1.05 and 2.25). At 7 they are 1.19 and 2.31. The runs' own scores in `run_index.csv` were made at 10.
- In 2010–2024 the model sits closer to the dune line than to CoastSat. CoastSat's 2021 +17 m step is the part of the target the model does not make.

## Other files

- `tables/domain_rates_<w>.csv`: every reading and the model per domain, with the residuals.
- `runs_used.csv`: run, folder, commit, topography, offsets, and the dune line's vintages and dates.
- `y_bounds.txt`: the shared y range (±10 m/yr) and the rule behind it.

Drawn by `scripts/analyze_output/compare_runs/rate_windows.py`:
`python scripts/analyze_output/compare_runs/rate_windows.py`
