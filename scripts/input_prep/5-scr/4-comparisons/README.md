# 4-comparisons — one source against another

Answers, not model inputs. Nothing the model reads is produced here; these
scripts take the products in `../3-rates/` and put two of them side by side.
**They read the stored tables and never recompute a rate** — if a comparison
disagrees with a rate figure, the rate product is what to fix.

One folder per question.

```
shoreline_vs_duneline/   Does the digitized dune line move with the CoastSat
                         shoreline?
dsas_vs_coastsat/        Do the two rate sources agree?
duneline_positions/      Where did the dune lines actually sit? (positions,
                         not change)
```

## shoreline_vs_duneline

Both sides are net change in metres, over the same windows, with the same sign
convention (seaward positive). That symmetry was the point of rebuilding the
dune side as an endpoint product in September 2026.

```
coastsat_vs_duneline.py
    The base comparison: 3-rates/duneline/endpoint/ against
    3-rates/coastsat/endpoint/, per window. It is also the module the other
    three import — load_chainage, the survey-date table, the beach-width
    shading — so it is the one to read first.

total_change_vs_duneline.py
    Total shoreline change (the LRR over the window it was fitted on, turned
    into metres) against the dune line's measured net change, over 1996-2010,
    2010-2024 and 1996-2024. Absorbed projected_vs_duneline.py on 2026-09-21
    when the vocabulary was settled: that script's numbers were already here
    as the *_dune_interval_m columns.

net_change_vs_duneline.py
    The same pairing over 1996-2024 and its two halves.

smoothed_loess7_vs_duneline.py
    Both curves through the model target's alongshore LOESS at 7 domains
    (3.5 km), in the two readings of the shoreline side, so the effect of the
    smoother is visible rather than assumed.
```

## dsas_vs_coastsat

Two scripts because the CoastSat side can be anchored two ways, and they
answer different questions:

```
dsas_vs_coastsat_raw.py
    Calendar windows, no smoothing: the raw per-domain mean LRR from each
    source. Separate from 6-scr-smooth/dsas_vs_coastsat/ on purpose — that
    folder exists to argue about the smoothing, and every figure in it draws
    a LOESS.

dsas_vs_coastsat_datematched.py
    The CoastSat side anchored on the shoreline SURVEY DATES instead: the mean
    position within +/-W days of each date, end minus start, at both +/-30
    days and +/-6 months. Not an OLS over a window.
```

## duneline_positions

```
duneline_positions.py
    Where the 1997, 2009 and 2023 dune lines sat — the lines standing for
    model years 1996, 2010 and 2024. Overview maps in three north-up segments,
    imagery zooms, dune-line-to-NC-12 distance, and beach width.
```

Two traps that folder's outputs record: a transect extended landward can cross
the road, and a seaward-drawn road line makes the dune-to-road distance
negative. Both are described in the data-tree README beside the figures.
