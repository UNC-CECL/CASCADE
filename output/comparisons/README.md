# comparisons - figures and tables that span more than one run

Nothing here comes from a single run; that lives in the run's own directory
under `output/raw_runs/`. Each folder answers one question, names the script
that writes it and the runs it read, and is only as current as those runs:
the run index records the topography version and git commit of every run,
which is what tells you whether a figure predates a change.

```
hindcast_calibrated/    the headline: modelled rate against the CoastSat
                        target, both windows; edgeBE and zeroBE, 1996-2015 and
                        2010-2026 since 10-03, model smoothed like the target
model_vs_observed/      the model against each window's CoastSat LRR; vs_shoreline/
                        on 1996-2015 + 2010-2026 since 10-03; the dune-line folders
                        are still the 1996-2010 / 2010-2024 record
matrix_vs_observed/     every matrix run on 1996-2015 + 2010-2026 (10-03): rate and
                        position change, start/end positions, against CoastSat
target_comparison/      total_change/runs_vs_coastsat on 1996-2015 + 2010-2026
                        (10-03); projected/ and the dune-line sets are the
                        1996-2010 / 2010-2024 record (09-19)
relocation/             does the model relocate NC-12 where and when history
                        did: per window, per event, across topography
                        versions, and the 20 m rebuild clearance
scenario_grid/          every preset and management scenario on one page
offset_source/          dune line vs shoreline as BRIE's island offset: how much
                        the model changes, and whether orientation drives it (09-28);
                        shoreline_v1_vs_v2/ inside it: the shoreline offset's
                        averaging window, calendar vs DEM-centred (10-01)
```

| folder | script | runs |
|---|---|---|
| `hindcast_calibrated/` | `scripts/figure_making/model_output/hindcast_final_figure_lowess.py` | the edgeBE and zeroBE full-management nogroin matrix arms, 1996-2015 and 2010-2026 |
| `model_vs_observed/` | `scripts/analyze_output/compare_runs/hindcast_vs_observed/rate_windows.py` | the edgeBE full-management nogroin matrix, 1996-2015 and 2010-2026 (`--windows 14yr`: the 14-yr matrix and the dune-line solve), `runs_used.csv` |
| `matrix_vs_observed/` | `scripts/analyze_output/compare_runs/matrix_vs_observed/matrix_vs_observed.py` | all 26 nogroin matrix runs on 1996-2015 and 2010-2026, `scores.csv` |
| `target_comparison/` | `scripts/analyze_output/compare_runs/hindcast_vs_observed/target_comparison.py` | the edgeBE and zeroBE full-management matrix, 1996-2015 and 2010-2026 (`--windows 14yr`: the dune-line and 1996-2024 LRR solves), `runs_used.csv` |
| `relocation/` | `scripts/hatteras_ms/experiments/HAT_relocation_comparison.py` and the three readers named in its README | each set's `report.txt` header |
| `scenario_grid/` | `scripts/figure_making/model_output/scenario_grid.py` | every matrix arm, both periods |
| `offset_source/` | `scripts/analyze_output/compare_runs/offset_source/offset_source_comparison.py` | the full-management duneline/shoreline pairs of `experiments/island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a`, 1996 and 2010, `tables/summary.csv` |
| `offset_source/shoreline_v1_vs_v2/` | `scripts/analyze_output/compare_runs/offset_source/shoreline_v1_vs_v2_comparison.py` | the full-management shoreline_v1/shoreline_v2 pairs of `experiments/island-offset/2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup`, 1996 and 2010, `tables/summary.csv` |

## Layout rule

Question first, then the window, then the version or preset. A window is an
inner axis, never a top-level folder (`relocation/1984_2004/`, not
`relocation_1984_2004/`), so one question is one folder however many windows
it is asked over. A comparison ACROSS versions or windows sits beside the
per-version sets it reads (`relocation/versions/`, `relocation/events/`), not
inside one of them.

## Tracked

One rule for the whole tree (`.gitignore`, the comparisons block): every
README, CAPTIONS.md, `.csv` and `.txt` is tracked; images, PDFs, GIFs and
`.npz` are regenerable from those plus the code and stay ignored. A number
quoted in a results document therefore survives a clone with the table it
came from.

## Removed 2026-09-17

Three sets were deleted rather than archived because nothing outside the
folder read them and none could be regenerated:

- `smoothing_vs_cascade/` (2026-04-06): LOWESS smoothing of the CoastSat
  rates against a run named `HAT_1984_2004_SQ_BE_Hs2p0`, a pre-archive
  name with no run behind it. The smoothing itself is now the scoring
  target's treatment (`cascade_pipeline.hindcast.build_target_table`) and is
  drawn in `model_vs_observed/vs_shoreline/smoothed/`.
- `source_sink_zones/` (2026-06-19): four figures from runs deleted before
  the 2026-08-28 archive; the script's own note said so.
- `relocation_1984_2004/v1/` (2026-09-01, superseded topography since
  09-04): the numbers are quoted in
  `scripts/hatteras_ms/experiments/RELOCATION_COMPARISON_RESULTS.md`, and
  five of the six sets regenerate from the calibration tree.

`overwash/` was an empty folder; the overwash comparison lives at
`data/hatteras_init/8-overwash-analysis/`.
