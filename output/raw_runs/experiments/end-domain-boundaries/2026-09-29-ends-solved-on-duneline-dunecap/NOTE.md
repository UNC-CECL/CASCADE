# end-domain-boundaries/2026-09-29-ends-solved-on-duneline-dunecap

**Question.** After the dune-cap fix (2026-09-28, the other session; `../2026-09-28-ends-resolved-dunecap/`), what end rates at GIS 1 and 90 make the model match the dune line?

**Why again.** The fix reran every beach/dune-managed run. The 09-28 dune-line solve (`../2026-09-28-ends-solved-on-duneline-adopted/`) is full management, so its runs carried the uncapped interior. The fix does not reach the end domains: GIS 1–2 and 88–90 moved by under 0.0001 m/yr in the 1996 matrix run. It does move Buxton (GIS 3–9), Avon (18–34) and Tri-Village (67–86) by up to 0.85 m/yr, and `target_comparison/` scores that interior. (Hannah, 2026-09-29: "Re-solve both, then redraw".)

**Method.** Unchanged from 09-28:

- `be_dune_edgesolve_loop.py --exp end-domain-boundaries/2026-09-29-ends-solved-on-duneline-dunecap --windows 1996 2010 --smooth mean3 raw --tol 0.02`.
- It starts from the fixed matrix's zeroBE and edgeBE full-management runs.
- `be_dune_edgesolve_results.py --solved 1996:mean3:3 1996:raw:2 2010:mean3:4 2010:raw:4` wrote `solved.csv`, `skill.csv` and `RESULTS.md`.

**Answer** (GIS 1 / GIS 90, m/yr):

| window | reading | step | ends | residuals | 09-28, before the fix |
|---|---|---|---|---|---|
| 1996–2010 | mean3 | 3 | −2.8 / +7.9 | −0.002 / −0.016, converged | −2.8 / +7.9 |
| 1996–2010 | raw | 2 | +0.6 / +10.6 | +0.010 / −0.015, converged | +0.6 / +10.6 |
| 2010–2024 | mean3 | 4 | +2.1 / +13.6 | **−0.032** / −0.009, accepted | +2.1 / +14.2 |
| 2010–2024 | raw | 4 | +1.9 / +13.5 | +0.010 / +0.016, converged | +1.9 / +14.0 |

- 1996 is unchanged, as the ends of that run are.
- In 2010 GIS 90 drops by 0.5–0.6 m/yr. The matrix's own 2010 GIS 90 end dropped too, from +22.49 to +21.26.
- 2010 mean3 GIS 1 stalled at 0.032, the 0.1 m/yr probe resolution, as on 09-27 and 09-28. Steps 5–6 repeated the probe.

**Reading.** The dune-line ends still barely change the interior. Interior RMSE is within 0.03 m/yr of the CoastSat-solved run against both observations (`RESULTS.md`).

**Runs.** They stay on disk only (`.gitignore`). Logs are in `logs/`, and every probe is in `loop_log.csv`.
