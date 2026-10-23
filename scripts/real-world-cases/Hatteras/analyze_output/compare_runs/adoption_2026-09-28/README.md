# adoption_2026-09-28 - the run matrix before and after the 2026-09-28 adoption

What the 2026-09-28 adoption (`hatteras/adopted` Barrier3D: per-cell dune
ceilings and overwash fixes; 24 h storms) changed, scored over the whole run
matrix. A one-off study: it is tied to that date's archive and will not need
re-running unless the archive is.

| script | step | writes to `output/comparisons/adoption_2026-09-28/` |
|---|---|---|
| `adoption_before_after.py` | 1. score each side, then report | tables |
| `adoption_before_after_figures.py` | 2. figures from those tables | `figures/` |

## The scripts in detail

Moved here from `scripts/analyze_output/README.md` on 2026-09-30, when the
scripts were grouped by question.

## adoption_before_after.py

Hannah, 2026-09-28: "show me the comparison when the matrix finishes". It
compares the matrix before and after the 2026-09-28 adoption:

| side | runs | Barrier3D and forcing |
|---|---|---|
| before | `output/raw_runs/archive/2026-09-28-pre-ceiling/matrix/` | `fix/route-overwash-axis-swap` (49fd069), Dmaxel default (3.4 m NAVD88), storms `v3_72`, the pre-adoption LOWESS-7 ends |
| after | `output/raw_runs/matrix/` | `hatteras/adopted` (overwash fixes + per-cell dune ceilings), storms `v3_trim24`, the ends re-solved on it (`end-domain-boundaries/2026-09-28-ends-resolved-adopted`) |

For every matrix run (both windows, both presets, every scenario) it scores:

- **shoreline:** interior (GIS 2-89) RMSE and bias of the model LRR against the
  CoastSat LOWESS-7 target (`run_registry.skill_vs_target`, as the runner
  scores), and the spatial correlation r;
- **overwash:** against the imagery (`8-overwash-analysis`), each run dated by
  its own storm file: POD, POFD, PSS, timing r, space r;
- **dunes:** for 1996-2010 runs, the 2010 dune crest minus the 2009 lidar.

Beside the scores it writes per-domain tables: `cells_<side>.csv` (every image
x domain, observed and model overwash) and `crest_<side>.csv` (end-of-run crest
per GIS domain, m MHW). `crest_lidar_2009.csv` is the crest from the 2010-start
dune file. The figures are `adoption_before_after_figures.py`.

**Each side is scored under the Barrier3D it ran on**, so the storm sharing
uses that version's DuneGaps and DuneGrowth. `score --side before` must run
with `PYTHONPATH=<Barrier3D at 49fd069 + the ceiling feature, off>` (the
worktree `../Barrier3D-dune-ceiling`); `score --side after` uses the editable
install. The script refuses to score a side under the wrong one.

## adoption_before_after_figures.py

Hannah, 2026-09-28: "make figures of the before and after comparison". It
reads the tables `adoption_before_after.py` writes (`scores_`, `cells_`,
`crest_<side>.csv`, `crest_lidar_2009.csv`) and the runs' own
`shoreline_change_rate.csv`. edgeBE only: zeroBE is within 0.1 of it on every
score (see the README beside the tables). The relocation arms are left out;
in both windows they score the same as their non-relocation twins.

| figure | shows |
|---|---|
| `adoption_scorecard.png` | every score, before -> after, per scenario |
| `adoption_shoreline_alongshore` | model LRR vs CoastSat LOWESS-7, managed + natural |
| `adoption_overwash_map_<scenario>` | image x domain: hit / miss / false alarm |
| `adoption_overwash_by_image` | domains overwashed per image, grouped bars |
| `adoption_dune_crest_2010` | the 1996-2010 runs' 2010 crest vs the 2009 lidar |

Output: `output/comparisons/adoption_2026-09-28/figures/`. On each run it also
deletes the retired `adoption_overwash_by_domain` figure if one is left over.
