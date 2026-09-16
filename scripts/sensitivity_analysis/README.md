# sensitivity_analysis - varying one forcing at a time

Drives the hindcast repeatedly with one parameter moved off its calibrated
value, and reads the result.

| Script | Does |
|---|---|
| `HAT_hindcast_sensitivity.py` | the sweep driver |
| `HAT_plot_sensitivity.py` | the sweep's figures |
| `HAT_plot_hs_experiment.py` | the wave-height experiment specifically |
| `plot_sensitivity_vs_coastsat.py` | each cell against the observed rates |

Products go to `output/sensitivity_analysis/`. Its `figures/README.md` is one
of the three decision records in the output tree.

**A cell that moves a forcing earns a name token**, so it lands in its own
directory. Without one it would derive the same name as the matrix run it is
being compared against, and the last to finish would wear the production name.
Cells file under `output/raw_runs/sensitivity/<axis>/<period>/<preset>/` (the
driver sets `HAT_RUN_KIND=sensitivity`; the axis is read off the token), and
the index keys a cell on (run_name, kind, tag) so it and its matrix baseline
are two rows. No model state is written for a cell.
