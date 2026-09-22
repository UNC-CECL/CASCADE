# analyze_output - questions that span more than one run

Nothing here runs the model. Each script reads finished runs and writes to
`output/comparisons/`, resolved as `hat_figure_style.COMPARISONS_ROOT`
(2026-09-18); a figure finished for the manuscript also goes to
`output/figures/<subject>/`, from `scripts/figure_making/`.

```
compare_runs/
    rate_windows.py      LIVE. Observed vs modelled shoreline-change rate,
                             every window, the three end-solve model sets
                             -> comparisons/model_vs_observed/
    compare_runs.py      a general run-vs-run / run-vs-CoastSat tool.
                             RUNS_TO_COMPARE is EMPTY: every example in it is
                             commented out. Fill it before running
                             -> comparisons/<COMPARISON_NAME>/
overwash/
    compare_overwash_figures.py   a run's overwash (Qow) against the observed
    compare_overwash_observed.py  record; both resolve their run through
                                  run_registry (HAT_1984_2004_calibBE_road_bdm_groin
                                  resolves today) -> comparisons/overwash/
    superseded_20260918/          plot_overwash.py (dead absolute path) and an
                                  early copy of Roya's Pea Island script
smoothing_vs_cascade/
    smoothing_vs_cascade.py   what the LOESS smoothing does to the
                             comparison. Names HAT_1984_2004_SQ_BE_Hs2p0, which
                             no longer exists anywhere under raw_runs/: point it
                             at a current run before running
```

The observed overwash record and its own figures are
`scripts/input_prep/8-overwash-analysis/`; this folder only compares a RUN
against it.

The smoothing figures and the source/sink-zone figures were deleted from
`output/comparisons/` on 2026-09-17 (unregenerable, uncited); running either
script recreates its folder, so give it a current run first.

The run tree is addressed through `cascade_pipeline.run_registry`, never by
building a path: a run's NAME describes its scenario and its PATH describes its
forcing, and the two are easy to get wrong by hand.
