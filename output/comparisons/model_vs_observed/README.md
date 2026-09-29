# model_vs_observed: the option A matrix, adopted model, against the shoreline and the dune line

> **Units: this tree is in m/yr.** The other comparison trees (`target_comparison/`,
> `4-comparisons/shoreline_vs_duneline/`, `3-rates/coastsat/{total_change,projected}/`)
> are in **metres**. [`FIGURES.md`](../../../FIGURES.md) indexes all of them.

**Redrawn 2026-09-29 after the dune-cap fix** (the other session, 2026-09-28: `experiments/end-domain-boundaries/2026-09-28-ends-resolved-dunecap/`). The fix reran every beach/dune-managed run. It moves the managed domains at Buxton (GIS 3–9), Avon (18–34) and Tri-Village (67–86), not the end domains. The matrix's 2010 GIS 90 end went from +22.4937 to +21.2582 m/yr, and the dune-line and 1996–2024 LRR ends were re-solved on the fixed runs (Hannah, 2026-09-29). The model is 1–2 m more seaward in 1996–2010 and about 1 m in 2010–2024. The pre-fix numbers below are recoverable from git (commit 7de36886); the pre-fix runs are in `raw_runs/archive/2026-09-28-pre-dunecap/`.

**Redrawn 2026-09-28 on the adopted model** (Hannah: "include the overwash fixes, keep option A, go ahead"; then "redraw the model vs observed figures"). The matrix was rebuilt on Barrier3D `hatteras/adopted` (the three overwash fixes and per-cell dune ceilings), storms `v3_trim24` (every event kept, trimmed to 24 h around its peak), its ends re-solved on that model, and the dune-line ends re-solved on it too. The option A pre-adoption version, first drawn 2026-09-27, is recoverable from git (commit 70f2efb3); its runs are in `raw_runs/archive/2026-09-28-pre-ceiling/`.

The model line comes from these runs:

- metres island offset;
- option A wave climate: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5;
- end rates solved on the adopted model against the LOESS-7 target (`experiments/end-domain-boundaries/2026-09-28-ends-resolved-adopted/`): 1996 GIS 1 +4.3509 / GIS 90 +19.0935, 2010 +8.0 / +21.2582 m/yr (after the dune-cap fix; +22.4937 before it). Before adoption they were +4.8394 / +18.2545 and +18.8657 / +24.2358;
- edgeBE, full management, groin off;
- Barrier3D `hatteras/adopted` (the three overwash fixes and per-cell dune ceilings), storms `v3_trim24` (every event kept, trimmed to 24 h around its peak).

| window | run |
|---|---|
| 1996–2010 | `matrix/1996_2010/edgeBE/HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin` |
| 2010–2024 | `matrix/2010_2024/edgeBE/HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin` (nourishment on) |
| 1984–2004, 2004–2024 | no metres run yet; each panel draws the observation and says "model not yet run" |

**Each observation is drawn with the run solved on it**, so each comparison is fair at GIS 1 and 90 by construction:
- The shoreline (CoastSat) figures draw the matrix run above, with ends solved on CoastSat.
- The dune-line figures draw the run with ends solved on the dune line on the dune-cap-fixed model: `experiments/end-domain-boundaries/2026-09-29-ends-solved-on-duneline-dunecap/`, mean3 reading. Its ends are 1996 −2.8 / +7.9 and 2010 +2.1 / +13.6 m/yr (before the fix +2.1 / +14.2 in 2010; before adoption −3.0 / +7.6 and +3.4 / +15.0).
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
| Both observations as net change at the same two dates | `vs_shoreline_and_duneline/change_rate/` (m/yr) | CoastSat blue, dune line red, both LOESS | endpoint rate |
| ...as net change in position (m) over the window | `vs_shoreline_and_duneline/net_change/` | the same two lines x the window's calendar years (measured over the survey interval, scaled to 14 / 20 yr, as `target_comparison/`) | last annual shoreline minus first |
| Endpoint observation against the model's OLS rate | `sensitivity/mixed-estimator/` | dune-line net change | OLS rate |

Each folder holds one figure per window (`<stem>_<start>_<end>.png`) plus a 2 × 2 `_grid`.
PDFs and `CAPTIONS.md` are under `supporting/`.

## Skill, GIS 2–89 (`tables/skill.csv`)

| window | target | model estimator | bias (m/yr) | RMSE (m/yr) |
|---|---|---|---|---|
| 1996–2010 | CoastSat LRR, smoothed (the scoring target) | OLS | +0.06 | 1.17 |
| 1996–2010 | dune line net change, smoothed (dune-solved run) | endpoint | +1.00 | 1.72 |
| 2010–2024 | CoastSat LRR, smoothed (the scoring target) | OLS | −1.36 | 2.07 |
| 2010–2024 | dune line net change, smoothed (dune-solved run) | endpoint | −0.31 | 1.42 |

After the dune-cap fix. Before it: −0.03 / 1.17, +0.91 / 1.70, −1.45 / 2.12, −0.39 / 1.45.

On the adopted model, 2010–2024 is less erosive against both observations: bias −1.66 → −1.45 m/yr against CoastSat and −0.67 → −0.39 against the dune line; the dune-cap fix takes them to −1.36 and −0.31. 1996–2010 barely moves.

- **Smoothed at 7 domains since 2026-09-28**, following the runner's scoring target.
- At 10 domains (`output/archive/2026-09-28_option-a-loess10-comparisons/`) the four RMSEs were 1.05, 1.59, 2.25 and 1.43.
- The narrower window leaves more alongshore detail in the target for the model to miss, so RMSE rises. Bias barely moves.

- At 10 domains the two scoring-target rows reproduced the matrix and the wave-recommendation numbers (RMSE 1.05 and 2.25). At 7, before adoption, they were 1.19 and 2.31; on the adopted model 1.17 and 2.12; after the dune-cap fix 1.17 and 2.07. The runs' own scores in `run_index.csv` were made at 10.
- In 2010–2024 the model sits closer to the dune line than to CoastSat. CoastSat's 2021 +17 m step is the part of the target the model does not make.

## Other files

- `tables/domain_rates_<w>.csv`: every reading and the model per domain, with the residuals.
- `runs_used.csv`: run, folder, commit, topography, offsets, and the dune line's vintages and dates.
- `y_bounds.txt`: the shared y range (±10 m/yr) and the rule behind it.

Drawn by `scripts/analyze_output/compare_runs/rate_windows.py`:
`python scripts/analyze_output/compare_runs/rate_windows.py`
