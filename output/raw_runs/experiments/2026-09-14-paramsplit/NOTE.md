# 2026-09-14-paramsplit

**Question.** Did splitting the site configuration on 2026-09-14 (the
`hatteras_site_config_prebe_20260914_*` backups mark the steps) change the
model? Four runs of the same cell, `HAT_1984_2004_calibBE_road_reloc_bdm_nogroin`
on 1984-start v1, under the old code (`oldcode`), the split with the sync on
(`syncon`), a check (`check`) and the final state (`final`).

**Answer.** All four carry interior RMSE 0.523026: bit-identical, the split
changed nothing. Compare with `tools/HAT_compare_rerun.py`.

**Runs deletable?** Yes; the answer is the four identical numbers in
`run_index.csv`, which `retired_runs.csv` keeps.
