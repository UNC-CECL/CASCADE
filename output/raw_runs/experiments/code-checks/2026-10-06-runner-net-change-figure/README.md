# 2026-10-06-runner-net-change-figure

**Question.** Does the runner draw the net-change figure for a DEM-to-DEM run?

**Answer.** Yes. Every run in a window listed in `NET_CHANGE_WINDOWS` now draws `figures/shoreline_position_change_with_buffers.png` (both sides LOWESS-7, ±100 m axis). The code is `scripts/cascade_pipeline/plotting/net_change_comparison.py`.

**Folder names.** `1996_2009/edgeBE/<run>/`: one edgeBE, no-groin calibration run, re-run so the figure could be checked.

**Status.** Record.
