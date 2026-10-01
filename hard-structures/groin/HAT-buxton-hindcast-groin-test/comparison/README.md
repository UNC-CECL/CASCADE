# comparison - the three-run figure for the groin rig

`HAT_three_run_comparison.py` draws the 1967-2017 Buxton rig three ways,
side by side: no groin and no nourishment, nourishment only, and nourishment
plus the best-fit groin. Each panel holds the 1967 shoreline, the modelled
2017 position and the observed 2018 target. It writes two static figures
(titled, and `_v2` without the figure title) and two GIFs to
`three_run_comparison/` beside the script. Nothing here runs the model.

    python HAT_three_run_comparison.py

**Make the three runs first**, with
`../../HAT-groin-buxton-output/1967_2017_run/HAT_groin_threeway_hindcast_1967_2017.py`
(or the single-run script with these settings), then set `RUN_BASELINE`,
`RUN_NOURISHMENT_ONLY` and `RUN_FULL_MODEL` to their folder names:

1. Baseline: `RUN_MATRIX = ["no_groin"]`, `ENABLE_HISTORICAL_NOURISHMENT = False`,
   and a distinct `RUN_NAME_SUFFIX` (e.g. `..._no_BN`): nourishment is not part
   of the run name, so without it runs 1 and 2 overwrite each other.
2. Nourishment only: `RUN_MATRIX = ["no_groin"]`, `ENABLE_HISTORICAL_NOURISHMENT = True`.
3. Full model: `RUN_MATRIX = ["groin"]`, `ENABLE_HISTORICAL_NOURISHMENT = True`,
   with the sweep's best `GROIN_TRAPPING_RATE_M_YR` / `GROIN_DETERIORATION_FRACTION`.

Brought in line with `scripts/STYLE.md` on 2026-10-01 and proven unchanged in
behaviour; the explanations that left the code are below, word for word.

**Paths, fixed 2026-10-01.** The script used to build every path on
`PROJECT_BASE_DIR = r"/"` and the pre-move `scripts/groin/...` tree, so it
could not find its runs, its observations or its output folder on any
machine. It now finds the repo root by searching upward and reads the runs
from `output/calibration/groin_rig/`, the wet/dry table from
`../../HAT-groin-buxton-output/shoreline_position_output/`, and writes to
`three_run_comparison/` beside itself, where the existing figures already
are. Of the three runs, `..._no_groin` and `..._groin` exist in
`groin_rig/`; the baseline `HAT_1967_2018_edge_calibrated_no_BN_no_groin`
does not, so the script stops at the first panel until that run is made.

## The scripts in detail

### HAT_three_run_comparison.py

From the script's original header:

```text
HAT_three_run_comparison.py
=============================
Three-panel comparison figure: baseline (no groin, no nourishment) vs
nourishment-only vs full model (nourishment + best-fit groin from the
sensitivity sweep). Each panel shows the 1967 initial shoreline, the 2017
modeled position, and the observed 2017/2018 target -- so you can see, left
to right, how each piece of added functionality moves the model closer to
(or away from) reality. All three panels are treated identically (no
intermediate-year overlay on any of them, for consistency across panels).

X-AXIS: domain positions are shown as a simple 1-11 index for this figure
(NOT the underlying GIS D2-D12 numbering used everywhere else in the
project) -- easier to read for this specific comparison. GROIN_BOUNDARY_GIS
(5.5, the real D5/D6 interface) is converted to this same 1-11 scale
(4.5) so the groin marker lands in the correct spot either way.

SIGN CONVENTION -- corrected after visual review, verified numerically before
applying (see _flip()'s docstring and fig_three_run_comparison()): raw x_s
increases LANDWARD (Barrier3D's native direction). This script flips it so
'+' = seaward/accretion, '-' = landward/erosion -- matching FLIP_SIGN_MODEL
used in every OTHER script in this project, AND matching the standard
erosion=negative / accretion=positive convention already used throughout
this dissertation's own LRR analysis. With that flip in place, reconstructing
an observed position requires SUBTRACTING the observed wet/dry CHANGE from
the 1967 reference (planform_1967 - observed_change), matching
HAT_groin_effect_comparison.py's formula exactly -- an earlier version of
this script skipped the flip and used addition instead, which made erosion
plot as a positive number, backwards from the rest of the project.

Produces TWO versions each run: HAT_three_run_comparison.png (full title +
subtitle) and HAT_three_run_comparison_v2.png (no figure-level title, for
dropping into a paper with its own caption) -- per-panel titles ("No groin,
no nourishment", etc.) are kept in BOTH versions, since with three visually
distinct scenarios some in-figure labeling is still needed for a reader to
tell panels apart even when a caption describes the figure as a whole.

WORKFLOW -- produce the three runs BEFORE running this script:
  1. Baseline (no groin, no nourishment):
       RUN_MATRIX = ["no_groin"]; ENABLE_HISTORICAL_NOURISHMENT = False
       Give this run a DISTINCT RUN_NAME_SUFFIX (e.g. "..._no_BN") --
       nourishment on/off isn't part of the run_name, so without a distinct
       suffix this run and #2 below would silently overwrite each other.
  2. Nourishment only (no groin):
       RUN_MATRIX = ["no_groin"]; ENABLE_HISTORICAL_NOURISHMENT = True
  3. Full model (nourishment + best-fit groin):
       RUN_MATRIX = ["groin"]; ENABLE_HISTORICAL_NOURISHMENT = True
       GROIN_TRAPPING_RATE_M_YR / GROIN_DETERIORATION_FRACTION = your
       sweep's best combination

Set RUN_BASELINE / RUN_NOURISHMENT_ONLY / RUN_FULL_MODEL below to the three
resulting run folder names, then run this script directly. Produces both the
static comparison figures (titled + v2) AND an animated GIF version --
make_three_run_evolution_gif() -- showing all three panels evolving
year-by-year (1967-2017) simultaneously, same layout/colors/conventions as
the static figure, with each panel's own 1967 reference and 2018 observed
target held fixed while its modeled line moves.
```

Notes moved out of the code (_display_axis(), its docstring, as of the original):

```text
Simple 1-11 index for this figure's x-axis, instead of the underlying
GIS D2-D12 numbering used everywhere else in the project.
```

Notes moved out of the code (_flip(), its docstring, as of the original):

```text
Barrier3D's raw x_s increases LANDWARD (erosion). Flip so '+' =
seaward/accretion, '-' = landward/erosion -- matches FLIP_SIGN_MODEL=True
used in every OTHER script in this project, and matches the standard
DSAS/shoreline-change convention (erosion = negative, accretion =
positive) already used throughout this dissertation's own LRR analysis.
Verified numerically before using this (see module docstring).
```

Notes moved out of the code (load_observed_changes(), its docstring, as of the original):

```text
Observed 1967->year change (RAW sign, '+' = landward) for each year,
from the validated wet/dry change table.
```

Notes moved out of the code (make_three_run_evolution_gif(), its docstring, as of the original):

```text
Animated version of fig_three_run_comparison(): all three panels
evolve simultaneously, year-by-year (1967 through each run's true final
year), each keeping its OWN static 1967 reference and 2018 observed
target fixed while its modeled line moves -- same layout, colors, and
sign convention as the static figure, so the two are a direct visual
match (just one frozen at 2017, the other animating through it).
```

Notes moved out of the code (comment above line 102, as of the original):

```text
# NOTE: this deliberately differs from the orange/red no-groin/groin
# convention used in every OTHER figure in this project -- three categories
# to distinguish here, not a binary, so a red/amber/blue "problem -> partial
# fix -> full fix" progression reads better for THIS figure specifically.
```

Notes moved out of the code (comment above line 115, as of the original):

```text
# 1=top). Lowered from 0.9 -- at the larger legend
# font size needed for legibility, the legend box
# got physically bigger and the label started
# overlapping it at 0.9. Checked directly via
# bounding-box overlap, not just eyeballed: 0.8 was
# the first value that cleared it, 0.75 adds a
# small safety margin for real data with a
# different y-range than the synthetic test data.
```

Notes moved out of the code (comment above line 122, as of the original):

```text
# every tick, to leave room for the larger font sizes below
# without crowding (advisor's suggested fix, not just mine)
# No OCEAN_AT_BOTTOM / axis-inversion flag here, unlike the older flipped-
# convention scripts -- see the y-axis section in fig_three_run_comparison()
# for why this script's raw ('+' = landward/erosion) convention needs no
# inversion at all to get the same visual result.
```

Notes moved out of the code (comment above line 128, as of the original):

```text
# --- Font sizes -- deliberately large. This figure gets shrunk to fit a
# paper/dissertation page width (roughly 1/3 of its native 20" rendered
# width), so anything sized to "look right" at native resolution reads as
# illegibly tiny once shrunk. Advisor feedback: must be legible at 100%
# scale viewing -- tick labels, axis labels, AND legend. ---
```

Notes moved out of the code (comment above line 241, as of the original):

```text
# FLIPPED convention: '+' = seaward/accretion, '-' = landward/erosion
# (see _flip() docstring). Reconstructing the observed position now
# means SUBTRACTING the raw observed change (opposite of an earlier,
# incorrect version of this script that didn't flip and added
# instead) -- verified numerically before making this change.
```

Notes moved out of the code (comment above line 259, as of the original):

```text
# Observed final (2018, full coverage) -- reconstructed from the
# model's own 1967 reference MINUS the observed CHANGE since 1967
# (raw sign, '+' = landward -- subtracting it from the flipped
# planform correctly pushes erosion negative).
```

Notes moved out of the code (comment above line 281, as of the original):

```text
# Axis inversion restored: with the FLIPPED convention ('+' = seaward/
# accretion, '-' = landward/erosion), inverting the axis puts erosion
# (negative values) above the zero line and accretion (positive) below
# it -- matching both the standard erosion=negative convention AND the
# visual "erosion above zero" layout, simultaneously.
```

Notes moved out of the code (comment above line 290, as of the original):

```text
# _mark_groin() is called HERE, only after the y-limits are fully
# finalized -- NOT inside the loop above. With sharey=True, each panel's
# ax.get_ylim() keeps changing as LATER panels get plotted (matplotlib
# synchronizes the shared range progressively), so calling it mid-loop
# used a stale, too-narrow range and placed the label far from the top.
# Verified this directly before fixing it (not just guessed).
```

Notes moved out of the code (comment above line 373, as of the original):

```text
# Fix y-limits across the WHOLE animation (all three panels, all years,
# both static reference lines) up front, computed BEFORE inverting --
# same two-step pattern (set normal range, then reverse) already used
# in the other GIF-making functions in this project, to avoid any
# order-of-operations confusion with combining a shared range and an
# inverted axis.
```

Notes moved out of the code (comment above line 389, as of the original):

```text
# _mark_groin() must be called AFTER the shared ylim above is set, not
# inside the loop -- same ordering bug as the static figure had: calling
# it earlier uses each axis's original auto-scaled range (from its own
# plot() calls only), not the explicit all-panel range set here.
```
