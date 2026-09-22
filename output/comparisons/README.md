# comparisons - figures and tables that span more than one run

Nothing here comes from a single run; that lives in the run's own directory
under `output/raw_runs/`. Each folder answers one question, names the script
that writes it and the runs it read, and is only as current as those runs:
the run index records the topography version and git commit of every run,
which is what tells you whether a figure predates a change.

```
hindcast_calibrated/    the headline: modelled rate against the CoastSat
                        target, both periods; edgeBE and zeroBE on the
                        1996/2010 chain since 09-18 (calibBE is not solved there)
model_vs_observed/      the four rate windows on one y axis, the model against
                        the shoreline (CoastSat) and the dune line, each with
                        the run solved on it (rate_windows/ until 09-18)
target_comparison/      CoastSat vs the dune line as the target, with the model,
                        as net change (m) over 1996-2010 and 2010-2024 (09-19)
relocation/             does the model relocate NC-12 where and when history
                        did: per window, per event, across topography
                        versions, and the 20 m rebuild clearance
scenario_grid/          every preset and management scenario on one page
```

| folder | script | runs |
|---|---|---|
| `hindcast_calibrated/` | `scripts/figure_making/model_output/hindcast_final_figure_loess.py` | the edgeBE and zeroBE full-management nogroin matrix arms, 1996 and 2010 |
| `model_vs_observed/` | `scripts/analyze_output/compare_runs/rate_windows.py` | the edgeBE nogroin matrix plus the dune-solved experiment, `runs_used.csv` |
| `target_comparison/` | `scripts/analyze_output/compare_runs/target_comparison.py` | the edgeBE matrix and the 09-18 dune edge solve, 1996 and 2010, `runs_used.csv` |
| `relocation/` | `scripts/hatteras_ms/experiments/HAT_relocation_comparison.py` and the three readers named in its README | each set's `report.txt` header |
| `scenario_grid/` | `scripts/figure_making/model_output/scenario_grid.py` | every matrix arm, both periods |

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

- `smoothing_vs_cascade/` (2026-04-06): LOESS smoothing of the CoastSat
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
