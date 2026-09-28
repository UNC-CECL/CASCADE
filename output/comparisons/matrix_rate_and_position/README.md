# matrix_rate_and_position: every option A matrix run, rate and position change

Made 2026-09-27 (Hannah: "Where are these position plots?"). The runner draws only the rate figures,
so until now the matrix had no position-change figure.

The runs are the option A no-groin matrix:
- metres offset;
- Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5;
- 1996–2010 and 2010–2024;
- edgeBE and zeroBE, every scenario.

Each figure has two panels:
- **(a) rate:** the model's OLS rate against the CoastSat LRR scoring target (10-domain LOESS, raw means GIS 1–10).
- **(b) position change:** the model's endpoint rate × 14 yr against the observed CoastSat change. The observed change is the mean position in the last calendar year minus the first, smoothed at 10 domains (`5-scr/3-rates/coastsat/total_change/<w>/smoothed`).

Seaward is positive. Scores cover the interior, GIS 2–89.

**Every figure of one kind uses the same y axis**, so any two can be compared directly:

| figure | y range |
|---|---|
| rate panel here, and each run's `shoreline_change_rate.png` (real domains) | −7.5 to +7.5 m/yr |
| each run's `shoreline_change_rate_with_buffers.png` | −10 to +10 m/yr (the buffers reach −9.6) |
| position panel here | −130 to +130 m |

- The runner's figures were redrawn with `rerender_run_figures.py --arm matrix --ylim=-10,10 --ylim-real=-7.5,7.5`.
- Every model line, CoastSat LOESS and domain mean fits inside ±7.5.
- 8 of about 1,800 CoastSat transect dots fall outside it and are cut off at the edge: GIS 1 in 2010–2024 (up to +9.1) and near GIS 80 in 1996–2010 (down to −7.2). The with-buffers figure shows them.
- `model_vs_observed/` keeps its own ±10 rule. `y_bounds.txt` gives the ranges.

| file | what |
|---|---|
| `<window>/scenarios_rate_and_position_<preset>_<window>.png` | every scenario of one preset on one pair of panels (relocation arms left out; they match their twins to 0.001 m/yr) |
| `<window>/by_scenario/rate_and_position_<preset>_<scenario>[_reloc]_<window>.png` | one run per figure. The same figure is in that run's own `figures/rate_and_position_change.png` |
| `scores.csv` | per run: rate bias and RMSE (m/yr), position bias and RMSE (m) |

Start with `2010_2024/scenarios_rate_and_position_edgeBE_2010_2024.png` and its 1996 twin.

Drawn by `scripts/analyze_output/compare_runs/matrix_rate_and_position.py`.
