# 2026-09-18-dune-edgesolve

**Question.** The same as `../2026-09-16-dune-edgesolve/`, re-asked after the 1997, 2009
and 2023 dune lines were re-digitized on 2026-09-18: what do the two locked
end domains (GIS 1 and 90) carry when they are solved against the DUNE LINE
rather than CoastSat? Same method and the same two readings (`raw`, `mean3`),
with the endpoint estimator on the model side.

**Re-solved:** 1996-2010, 2004-2024 and 2010-2024, the windows whose dune
target changed. For 1996 and 2010 the island offsets changed too (1996/v3,
2010/v2; renumbered v1 and v1 on 2026-09-19).

**Carried over from 09-16, not re-run:** 1984-2004. Neither of its lines
(1984, 2004) changed, and its inputs are the same.

**Brackets (step 0):** the re-run 1996 and 2010 matrix zeroBE and edgeBE
full-management runs (09-18, on the new offsets), and the 09-16 2004 brackets
under `experiments/2026-09-16-dune-edgesolve/brackets/`. The 2004-start inputs
did not change.

**Target:** read from the stored product `5-scr/3-rates/duneline/endpoint/`,
the same numbers `HAT_rate_windows.py` draws.

**Driven by** `scripts/input_prep/7-source-sink/2-calibrate/HAT_be_dune_edgesolve_loop.py`:
lockstep steps, stopping when both ends are within 0.01 m/yr of the target.
Each step's residuals and next probe are in `loop_log.csv`. The results are
written by `HAT_be_dune_edgesolve_results.py --exp 2026-09-18-dune-edgesolve`
into `solved.csv`, `skill.csv` and `RESULTS.md`, with the 1984 rows carried
from 09-16.

**Layout:** `<reading>/step<k>/<window>/edgeBE/<run>/`, `logs/`.
