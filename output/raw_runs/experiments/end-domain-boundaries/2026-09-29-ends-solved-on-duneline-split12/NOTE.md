# end-domain-boundaries/2026-09-29-ends-solved-on-duneline-split12

**Question.** After the storms changed to `v3_split12_trim24` (2026-09-29; grouped events split at ≥12 h below the berm, so Fran 1996 and Jose 2017 are back), what end rates at GIS 1 and 90 make the model match the dune line?

**Why again.** The matrix was re-run on the new storms, with its CoastSat ends re-solved in `../2026-09-29-ends-resolved-split12/`. The dune-line solve starts from the matrix's zeroBE and edgeBE full-management runs, so it follows the matrix. Hannah, 2026-09-29: "re-run the matrix and re-solve the ends", then "redraw the model_vs_observed figures when the matrix finishes". The figures draw this solve beside the matrix.

**Method.** Unchanged from `../2026-09-29-ends-solved-on-duneline-dunecap/`:

    be_dune_edgesolve_loop.py --exp end-domain-boundaries/2026-09-29-ends-solved-on-duneline-split12
        --windows 1996 2010 --smooth mean3 raw --tol 0.02
    be_dune_edgesolve_results.py --exp ... --solved 1996:mean3:3 1996:raw:3 2010:mean3:5 2010:raw:4

**Answer** (GIS 1 / GIS 90, m/yr):

| window | reading | step | ends | residuals | trim24 storms (09-29 dunecap solve) |
|---|---|---|---|---|---|
| 1996–2010 | mean3 | 3 | −2.7 / +8.0 | **+0.032** / +0.023, accepted | −2.8 / +7.9 |
| 1996–2010 | raw | 3 | +0.6 / +10.7 | +0.005 / +0.004, converged | +0.6 / +10.6 |
| 2010–2024 | mean3 | 5 | +2.2 / +13.8 | −0.004 / +0.019, converged | +2.1 / +13.6 |
| 2010–2024 | raw | 4 | +2.0 / +13.7 | **+0.037** / −0.001, accepted | +1.9 / +13.5 |

- Every end moved by at most 0.1–0.2 m/yr, the same order as the CoastSat ends (+0.04).
- Two stalls sit at the 0.1 m/yr probe resolution: 1996 mean3 GIS 1 and 2010 raw GIS 1. Their later steps repeated the same probe. The 09-29 solve accepted 2010 mean3 at 0.032 the same way. The best step is the one used.

**Reading.** As before, the dune-line ends barely change the interior. Interior RMSE is within 0.03 m/yr of the CoastSat-solved run against both observations (`RESULTS.md`).

**Runs.** They stay on disk only (`.gitignore`). Every probe is in `loop_log.csv`, and the driver output is `output/logs/driver/duneline_endsolve_split12_20260929.log`.
