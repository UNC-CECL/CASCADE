# 2026-09-28 — REBUILD: dune line vs shoreline as the island offset, on the OLD ÷10 offset

**A reconstruction, not a record.** It re-runs the 09-25 design with the ÷10 offset and the
÷10-era waves swapped back in; nobody ran it at the time. The ÷10 experiment actually run
then is `../2026-09-22-div10-offset-shoreline-trial-original/` (edgeBE with the ÷10 ends, full
management, the unfixed Barrier3D). Renamed `-rebuild` 2026-09-28 so the name says so.

Hannah, 2026-09-28: rebuild the offset-source experiment on the ÷10 island
offset from before the metres fix, clearly labelled, to track the changes made
while experimenting. This is the **before** state. The same question on the
metres offset is in `../2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8/`
(09-25 waves) and `../2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/`
(option A waves).

| | |
|---|---|
| offset | **÷10**: offset_mode `asrun` (offset / 10, the historical unit error), on `1996/duneline/superseded_20260924_pre-metres/v1` and `1996/shoreline/superseded_20260924_pre-metres/v1`, the builds the ÷10 runs read. `asrun` divides the buffers too, so only these builds reproduce a ÷10 run |
| waves | the ÷10-era defaults, before 2026-09-24: Hs 2.5 m, Tp 8 s, asymmetry 0.7, high-angle 0.1 (Hannah's choice) |
| ends | zeroBE (no source/sink correction); relocations and groins off |
| scope | natural and full management, 1996–2010; 4 runs |
| Barrier3D | the current one, with the route_overwash fix. The ÷10-era runs used the unfixed router, which segfaulted on some storms |
| figures | the same form as the 09-25 and option A studies: net change in metres, observations smoothed with LOWESS over 7 domains (the southern 10 raw), model unsmoothed, no scores on the figures |

**Both the offset and the waves differ from the 09-25 study**, so a
difference between the two folders cannot be put on either alone.

## Layout

```
tables/all_runs.csv            the 4 runs: offset, scenario, scores, status, run folder
figures/duneline_offset_vs_duneline_change_full_management.png
figures/shoreline_offset_vs_coastsat_total_change_full_management.png
figures/shoreline_offset_vs_coastsat_projected_change_full_management.png
figures/total_change_difference_shoreline_minus_duneline_full_management.png
logs/<offset>_<scenario>/<settings>.log, logs/driver.log
runs/<offset>_<scenario>/1996_2010/zeroBE/<run_name>/     on disk only
```


**2010–2024 added 2026-09-28** (full-management pair at the same settings, `run-2010`), so the figures show full management × both periods like every island-offset study: (a) 1996–2010, (b) 2010–2024. The natural runs stay in `tables/all_runs.csv`, not in the figures.

## Reproduce

```
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py run --jobs 4
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py run-2010
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_div10.py plot
```

The driver is the 09-25 one (`HAT_offset_source_comparison.py`) with the run
environment switched to the ÷10 offset and the settings above.
