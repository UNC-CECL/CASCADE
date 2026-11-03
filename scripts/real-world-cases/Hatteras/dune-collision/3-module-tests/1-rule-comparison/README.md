# 1-rule-comparison — which collision erosion rule fits the lidar?

`HAT_collision_rule_comparison.py`. Offline, no model runs. It takes no arguments: running it refits and
rewrites its output. Rules: none, Long et al. (2014) fractions, Larson per storm or hourly, each with
and without the beach factor, fit to the lidar crest change on the calibration cells of
`../../4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py`, which it imports.
Output: the workspace's `dune-collision/3-module-tests/1-rule-comparison/` (offline, no runs; moved there
from `output/raw_runs/dune-storm-response/2026-10-29-collision-rule-comparison/` 2026-11-03). Its `NOTE.md`
also tabulates every in-model collision run.

Result: the beach factor matters more than the equation. Larson per storm x beach, the module, is best
(CV RMSE 0.49 m, Dennis r 0.86); Long's 10% per storm x beach is the simplest runner-up (0.50, 0.78);
10% as published is worst.
