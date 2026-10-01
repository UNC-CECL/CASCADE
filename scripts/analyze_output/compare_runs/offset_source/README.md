# offset_source - does the island offset's source change the model?

The dune-line-built and the shoreline-built island offset, run otherwise
identically (the 2026-09-28 option-A offset study), compared. Nothing here runs
the model; the study's runs are made by
`scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py`.

| script | writes to `output/comparisons/` |
|---|---|
| `offset_source_comparison.py` | `offset_source/` |
| `shoreline_v1_vs_v2_comparison.py` | `offset_source/shoreline_v1_vs_v2/` |

## The scripts in detail

Moved here from `scripts/analyze_output/README.md` on 2026-09-30, when the
scripts were grouped by question.

## offset_source_comparison.py

Asked by Hannah on 2026-09-28: a simple comparison of the model output
started from the dune-line offset against the model output started from the
shoreline offset, for 1996-2010 and 2010-2024. The main question is how much
the island's *orientation* in the offset affects the outcome.

**Runs** (no new runs): the full-management pair of each period from
`output/raw_runs/experiments/island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/runs/{duneline,shoreline}_full_management/<period>/zeroBE/<run>/`.
Option A waves, zeroBE ends, relocations and groins off; the only thing that
differs within a pair is the island offset.

**Offsets:** `2-brie-offset/<start>/{duneline,shoreline}/<version>/*_unpadded.csv`,
at the version each run's metadata records. BRIE adds the offset to x_s, so
larger = more *landward*; everything uses the *seaward* position s = -offset,
with the mean removed (a uniform shift does not change what BRIE does, and the
builds are not on a common datum).

**Quantities,** per GIS domain (500 m), shoreline start minus dune-line start:

- orientation, theta = atan(ds/dx), degrees;
- turning, d(theta)/dx, degrees per km; positive where the line bends
  landward (an embayment), negative at a bulge;
- model change, each run's LRR x 14 yr (m, seaward positive).

Statistics are on the interior, GIS 2-89.

**Target** (second version of the change figure, Hannah 2026-09-29):
projected shoreline change, the CoastSat LRR fitted on 1996-2024, LOWESS over
7 domains (southern 10 raw), x 14 yr: one profile, the same in both periods.
The model stays unsmoothed.

**Label placement** (Hannah, 2026-09-28): every label stays out of the data.
The villages and shoals are named in a strip above the highest line, and the
groin and piers are drawn below that strip and named in the legend. Draw
order, bottom to top: village and shoal bands, grid, zero line, groin and
piers, the two runs, labels; the grid has no ticks inside the label strip.

Output, `output/comparisons/offset_source/`:

| file | shows |
|---|---|
| `offset_source_model_change_full_management.png` | the two runs, both periods |
| `offset_source_model_change_vs_projected_full_management.png` | the same with the projected target on top (2026-09-29) |
| `offset_source_difference_full_management.png` | profiles |
| `offset_source_orientation_vs_model_full_management.png` | scatter |
| `tables/summary.csv`, `per_domain.csv`, `vs_projected.csv` | the numbers |

## shoreline_v1_vs_v2_comparison.py

Asked by Hannah on 2026-10-01: the same comparison for shoreline offset v1
(calendar window) against v2 (DEM-centred window), from the full-management
pairs of `output/raw_runs/experiments/island-offset/2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup/`.
It imports its offsets, orientation and projected target from
`offset_source_comparison.py`, so the two use the same definitions. It checks
each run's metadata offset version against the arm it is filed under. It writes two
figures (profiles with the target on top; scatter against offset and turning
difference) and three tables to `output/comparisons/offset_source/shoreline_v1_vs_v2/`.

<details><summary>Function notes (the original docstrings)</summary>

**`projected_target()`**

```
Projected shoreline change (m): the 1996-2024 CoastSat LRR target, built
as the runner builds it at LOWESS_DOMAINS, x YEARS.
```

**`vs_target()`**

```
Each run against the projected target, interior GIS 2-89: bias and RMS
residual are the numbers to read; explained and r beside them.
```

**`fig_change_only()`**

```
Panels (a, b) of the profile figure on their own: the two runs' change.
With `obs`, the second version: the projected target drawn on top of them.
Every label stays out of the data (Hannah, 2026-09-28): the villages and
shoals are named in a strip above the highest line, and the groin and
piers are drawn below that strip and named in the legend, not on the lines.
```

</details>
