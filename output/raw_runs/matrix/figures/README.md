# matrix figures: the management ladder

> **Redrawn 2026-10-03 on 1996–2015 and 2010–2026.** The matrix was re-run on 2026-10-03 on the windows 1996–2015 and 2010–2026 (the 2017 Buxton fill and the 1.62 M cy Rodanthe volume, storms `v3_split12_trim24`). Every comparison is against **each window's own CoastSat LRR**. In every smoothed reading, **both sides are smoothed alike**: the model's per-domain values get the target's own 7-domain LOWESS, with GIS 1–10 left raw (`HAT_metres_1_offset_units.smooth_like_target`). No dune line exists for 2015 or 2026, so nothing dune-line based is drawn on these windows. The figures they replace are in `output/archive/2026-10-03_comparisons-14yr-windows/`. Position change is the endpoint rate × each window's years. Both windows now have a fills rung (1996–2015 has Rodanthe 2014). Values beyond the axis are marked at the edge. The table below describes the 14-yr ladder.

The matrix runs of each window and preset in order of increasing management (Hannah,
2026-09-29). There are two versions per window and preset:

- `<window>/management_ladder_rate_<preset>_<window>.png`: the OLS rate against the CoastSat LRR scoring target, ±7.5 m/yr.
- `<window>/management_ladder_position_<preset>_<window>.png`: the position change over the window (endpoint rate × 14 yr) against the observed CoastSat change, ±130 m (the axis of `output/comparisons/matrix_vs_observed/`).

| panel | 1996–2010 | 2010–2024 |
|---|---|---|
| (a) | natural | natural |
| (b) | + road management | + road management |
| (c) | + beach and dune management | + beach and dune management (no fills) |
| (d) | | + nourishment fills |

- **Each panel is cumulative.** It draws every run up to that rung, each in its own blue: lighter means less managed, and the newest rung is the heaviest line. A rung keeps its colour down the figure, so a curve can be followed from panel to panel.
- **Scores:** the row title gives the newest rung's interior (GIS 2–89) RMSE and bias.
- **Missing rungs:** 1996–2010 has no fills, so its ladder stops at (c). Beach/dune-only (a side branch) and the historical-relocation runs are not on the ladder: in 1996–2010 relocation matched full management to 0.0002 m/yr, because it moves the road, not the shoreline, and 2010–2024 has no relocation event.

- Captions: `<window>/supporting/CAPTIONS.md`.
- Numbers: `supporting/ladder_scores.csv` (per rung: rate and position RMSE and bias).

Drawn by `scripts/analyze_output/compare_runs/matrix_vs_observed/matrix_management_ladder.py` from whatever
matrix runs are current. Re-run it after a matrix re-run.
