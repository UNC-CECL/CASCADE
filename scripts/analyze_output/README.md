# analyze_output - questions that span more than one run

Nothing here runs the model. Each script reads finished runs and writes to
`output/comparisons/`, resolved as `hat_figure_style.COMPARISONS_ROOT`
(2026-09-18); a figure finished for the manuscript also goes to
`output/figures/` (numbered layout, map in its README.md), from `scripts/figure_making/`.

```
compare_runs/            comparisons of finished runs, one subfolder per
                         question (map + details in compare_runs/README.md)
    compare_runs.py          general tool: up to four runs vs CoastSat.
                             RUNS_TO_COMPARE is EMPTY; fill it first
    hindcast_vs_observed/    LIVE. The hindcast vs shoreline and dune line
                             (rate_windows.py, target_comparison.py, ...)
                             -> comparisons/model_vs_observed/, target_comparison/
    matrix_vs_observed/      every run-matrix scenario vs the observations
    adoption_2026-09-28/     the matrix before/after the 09-28 adoption
    offset_source/           dune-line vs shoreline island offset
overwash/
    compare_overwash_figures.py   a run's overwash (Qow) against the imagery
                                  record: stacked, contingency, spatial figures.
                                  Reads the archived 1984-2004 calibBE run
                                  -> comparisons/overwash/ (details in its README)
```

The observed overwash record and its own figures are
`scripts/input_prep/8-overwash-analysis/`; this folder only compares a RUN
against it.

The run tree is addressed through `cascade_pipeline.run_registry`, never by
building a path: a run's NAME describes its scenario and its PATH describes its
forcing, and the two are easy to get wrong by hand.

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning, the
choices behind it and its history are in the README of the folder it lives
in (`compare_runs/README.md` is the map of that folder and its subfolders;
`overwash/README.md`). They were moved out of the scripts on 2026-09-30, when
the scripts were brought in line with `scripts/STYLE.md`, and out of this
README on 2026-09-30 and 2026-10-01, when each folder was organized.

Shared history: every script here drew in matplotlib's defaults until
2026-09-17, when each started calling `apply_style()` from
`site_layer/hat_figure_style.py` at import. Their paths used to be absolute
literals typed on one machine; they are all found by searching upward for the
repo root now (ORGANIZATION.md rule 5), so they follow the checkout and survive
a file changing depth.

## Deleted 2026-10-01: smoothing_vs_cascade/

`smoothing_vs_cascade/smoothing_vs_cascade.py` (2026-04-06) LOWESS-smoothed
the CoastSat LRR of 1984-2004 and 2004-2024 and drew one CASCADE run against
it. Hannah chose to delete it rather than park it (git keeps it; a departure
from ORGANIZATION.md rule 4, as for 5-scr on 2026-09-22):

- the run it drew, `HAT_1984_2004_SQ_BE_Hs2p0`, no longer exists anywhere
  under `raw_runs/`, so only its observation-only figures could still be made;
- it used the retired 1984/2004/2024 windows and a ~15-domain LOWESS width
  (`LOWESS_FRAC = 0.167`), not the 1996 -> 2010 -> 2024 chain and the group's
  7 domains;
- both its questions are answered on current settings in
  `compare_runs/hindcast_vs_observed/`: whether the smoothing width matters
  (`smoothing_scale.py`), and the model over the smoothed observations
  (`smoothed_lowess7_with_cascade.py`).

Its figures were already deleted from `output/comparisons/` on 2026-09-17
(see that folder's README). To recover the script:

```
git log --diff-filter=D --oneline -- scripts/analyze_output/smoothing_vs_cascade/smoothing_vs_cascade.py
git show <commit>^:scripts/analyze_output/smoothing_vs_cascade/smoothing_vs_cascade.py
```

Its last committed README section (figure list, bandwidth notes) is in the
same commit's parent: `git show <commit>^:scripts/analyze_output/README.md`.
