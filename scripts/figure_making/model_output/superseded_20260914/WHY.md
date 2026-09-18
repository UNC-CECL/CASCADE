# superseded 2026-09-14 - the headline hindcast figure

`HAT_hindcast_final_figure.py` drew the calibrated model against the CoastSat
target, both periods, scored over D2-D89. It was the project's headline result
figure from 2026-08 until now.

## Why it was retired

It and `HAT_hindcast_final_figure_loess.py` showed the same two runs. The only
substantive difference was the scoring window - D2-D89 here, D11-D89 there -
and on 2026-09-14 the LOESS figure was changed to print BOTH windows on every
panel. At that point this script drew nothing the surviving one did not, except
three shaded bands.

Keeping two was not free. Each figure quoted the other's number in its caption,
and this one's value sat hardcoded in the LOESS script for twelve days after
the 1984 run was remade on topography v2 - the caption said 0.526 / 0.558 while
this figure printed 0.547 / 0.580. That is the specific failure two figures of
the same thing produce, and removing one removes it.

## What was NOT carried over, deliberately

The three shaded bands - D5-D7 reserved for the groin, grey outside the frozen
zone set, purple D1/D90 locked. They had already been tried on the LOESS canvas
and removed: with the geographic annotation layer in place the Buxton area
carried four overlapping fills and read as clutter. D5-D7 and D1/D90 are stated
in that figure's caption instead, and the groin line still marks D5.5. The grey
frozen-zone band is the one fact now shown nowhere, judged the least
load-bearing of the three.

If that judgement turns out wrong, this script is the record of how the bands
were drawn, and `git log` has every version of it.

## What replaced it

`../HAT_hindcast_final_figure_loess.py`, `--preset calibBE` (default) or
`--preset edgeBE`, writing to
`output/comparisons/hindcast_calibrated/hindcast_<preset>_loess_reference.png`.

The two PNGs this script last produced were DELETED on 2026-09-14, the same
day. Keeping them would have left a second set of numbers in the tree with
nothing rebuilding them when a run changed, which is the failure this
retirement was meant to end. This script regenerates them; see below.

## Running it again

It still runs, from here, unchanged - it finds the repo root by walking up for
`pyproject.toml` rather than by counting directories, so moving it did not
break its paths. It writes to the same filenames it always did, which would
land them back in the live folder; move them afterwards, or change `OUT_DIR`.
