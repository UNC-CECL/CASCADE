# comparisons - figures and tables that span more than one run

Nothing here comes from a single run; that lives in the run's own directory
under `output/raw_runs/`. Each folder answers one question, names the script
that writes it and the runs it read, and is only as current as those runs:
the run index records the topography version and git commit of every run,
which is what tells you whether a figure predates a change.

Only comparisons on the current windows live here (see `raw_runs/matrix/README.md`).

```
full_window/            the 30-yr 1996-2025 run (shoreline v2 start, full
                        management, no groin): LRR rate and net change to the
                        calendar-2025 mean, plus a yearly change GIF (10-05)
```

| folder | script | runs |
|---|---|---|
| `full_window/` | `scripts/analyze_output/compare_runs/hindcast_vs_observed/full_window.py` | the zeroBE and edgeBE full-management nogroin runs in `raw_runs/matrix/1996_2025/`, `tables/scores.csv` |

The calibration (1996_2009) and test (2009_2025) runs are compared to the
observed net change by the runner itself: each run's
`figures/shoreline_position_change_with_buffers.png`.

## Layout rule

Question first, then the window, then the version or preset. A window is an
inner axis, never a top-level folder, so one question is one folder however
many windows it is asked over.

## Tracked

One rule for the whole tree (`.gitignore`, the comparisons block): every
README, CAPTIONS.md, `.csv` and `.txt` is tracked; images, PDFs, GIFs and
`.npz` are regenerable from those plus the code and stay ignored.

## Moved off C: on 2026-10-07

Every comparison on a superseded window: `hindcast_calibrated/`,
`model_vs_observed/`, `matrix_vs_observed/`, `target_comparison/`,
`alternative_targets/` (1996-2015 / 2010-2026 and the 14-yr record),
`relocation/` (1984-2004, 1996-2010), `scenario_grid/`, `offset_source/`,
`adoption_2026-09-28/` and the empty `overwash/`. They are on
`D:\CASCADE_offload\output\comparisons\`, same paths; their tracked tables and
READMEs are in git history before the commit that removed them. The scripts
that wrote them are unchanged, so any of them can be redrawn on the current
windows.
