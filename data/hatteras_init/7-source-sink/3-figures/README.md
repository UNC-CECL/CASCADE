# 3-figures — the source/sink zone rework, in pictures

One folder per period pair. `1-field/` shows the residual itself, twice: `1-rate/` in m/yr over the run and `2-net-change/` in metres from the starting shoreline; `2-method/` shows how the zones are chosen. Captions are in each folder's `supporting/CAPTIONS.md`.

| folder | figure | the question it answers |
|---|---|---|
| `1996_2009__2009_2025/1-field/1-rate/` and `2-net-change/` | `zone_residuals_<view>_edgeBE` | Step 1: where, and by how much, does the edgeBE full-management run miss the observed net change in each period? |
| | `zone_residuals_<view>_edgeBE_calibration_1996_2009`, `..._test_2009_2025` | Step 1: the same, one period per figure |
| | `zone_residuals_<view>_overlap_edgeBE` | Step 1: where do the two periods agree on a missing source or sink, and where do they disagree? |
| `1996_2009__2009_2025/2-method/` | `zone_candidates_edgeBE` | Step 2: which runs of domains are candidate zones, shared (scenario A) and per period (scenario B)? |

The earlier calibration's figures (the 1984/2004 and 1996/2010 pairs) were removed on 2026-10-12; their PDFs are in git history.
