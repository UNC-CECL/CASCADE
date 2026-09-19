# target_comparison — CoastSat or the dune line as the CASCADE target?

The two candidate targets and the hindcast, as **net change in position (m)
over each 14-yr model window**, 1996–2010 and 2010–2024. Built 2026-09-19
(Hannah, by interview) by `scripts/analyze_output/compare_runs/HAT_target_comparison.py`,
which reuses `HAT_rate_windows.py`'s loaders, so the observations and runs are
the ones `model_vs_observed/` draws as rates.

| line | what | over 14 yr how |
|---|---|---|
| blue | CoastSat target: the window's LRR (the runner's scoring series) | LRR x 14 |
| red | dune-line target: the MEASURED net change between the digitized lines (1997-10→2009-05, 11.6 yr; 2009-05→2023-07, 14.1 yr, end date assumed) | its rate x 14 (projected to the model years) |
| black | CASCADE: the run's own net change over the window, unchanged | endpoint rate x 14 |

The space between the two targets is the beach-width change they imply
(solid grey widened, hatched narrowed). Lines are raw domain means;
`tables/skill.csv` also scores against the LOESS-smoothed targets, the form
the runs are graded in.

**Two versions of the CoastSat target (2026-09-19, Hannah):**

- `coastsat_full_period_lrr/` — **the target in use**: the 1996–2024 LRR x 14 yr
  in BOTH windows, paired with runs whose ends were solved against it
  (`raw_runs/experiments/2026-09-19-edgesolve-lrr1996_2024/`: GIS 1 / 90 =
  +28.5 / +24.5 in 1996–2010, +37.1 / +25.4 in 2010–2024).
- `coastsat_subperiod_lrr/` — kept for the record: each window's own LRR (what
  the runner grades against), paired with the matrix runs.

**Every figure here uses one y axis, ±80 m** (Hannah, 2026-09-19); each caption names the few domain values that run off it (Cape Point in the sub-period 2010–2024 figures, and the dune line at GIS 12 in 1996–2010).

The dune-line target is the same in both (sub-period: 1997→2009 and
2009→2023, each scaled to 14 yr). Each version has the layout below.
`python ... HAT_target_comparison.py --coastsat-target full|subperiod`.

```
ends_solved_on_coastsat/   the matrix edgeBE run (end domains solved on CoastSat)
ends_solved_on_duneline/   the 09-18 dune edge-solve run, mean3 (solved on the dune line)
    target_comparison_1996_2010_2024.png   1996–2010 above 2010–2024; PDF, CAPTIONS under supporting/
paired/                    START HERE: each target with ITS OWN run only (CoastSat with the
                           CoastSat-solved run, the dune line with the dune-solved run), both
                           one figure per window, all on the same ±80 m axis; per window
                           (a) the CoastSat target as the house fill with its run in black,
                           (b) the dune-line target the same way; end values in the figure title.
                           Chosen 2026-09-19 over dots+lines and residuals (style B)
    target_and_own_run_1996_2010.png
    target_and_own_run_2010_2024.png
tables/
    domain_values_<w>.csv  per domain: both targets raw and LOESS (m), the dune line's
                           measured change and interval, both model sets (m)
    skill.csv              model minus target, GIS 2–89: bias, RMSE (m), r; per window,
                           model set, target (coastsat, coastsat_loess, duneline, duneline_loess)
runs_used.csv              the runs, with their index provenance
```

Hannah had not chosen the target when this was built; the two subfolders are
the two answers side by side.

**No other source/sink correction.** Every run here is edgeBE: the only
source/sink terms are the two end domains, GIS 1 and GIS 90 (checked from the
run index, `be_nonzero_domains == 2`; the paired figure refuses to draw
otherwise). GIS 2–89 carry none, so the interior of each model line is the
model's own response.

| window | ends solved on CoastSat (GIS 1 / 90, m/yr) | ends solved on the dune line |
|---|---|---|
| 1996–2010 | +32.2 / +10.0 | −16.0 / −3.3 |
| 2010–2024 | +72.6 / +31.3 | +20.1 / +19.2 |
