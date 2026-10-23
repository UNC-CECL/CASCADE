# 2026-10-20-spread-relaxation — the groin's step relaxes over three domains, not one

The result is the **spread groin (0, 2/3, 1/3)**: none of the updrift loss on GIS 6, two-thirds on GIS 7, one-third on
GIS 8, at strength 0.10 in both windows. Since 915690861 it is `spread_weights=(0, 2/3, 1/3)` in `cascade/groin.py`.

**Question.** Without the groin, the model's 1996 step at the groin (GIS 6 about 260 m seaward of GIS 5) relaxes
into GIS 6 alone: GIS 6 loses 71 m where the record has 9 m, and GIS 7-9 lose 0-12 m where it has 14-26 m. The
advisor asked for a three-domain northward reach and no background erosion. Does the groin fit better if the loss
GIS 6 takes across the groin face is shared over GIS 6, 7 and 8 (1/2, 1/3, 1/6)?

**Setup.** `spread_relaxation.py` patches the groin in-process (`cascade/groin.py` and the runner are not edited).
Each year, before BRIE's step, GIS 6's coming loss across the face is estimated from BRIE's own diffusion term,
f = 2 r (x_GIS5 - x_GIS6), r BRIE's diffusion number at GIS 6; GIS 6 is given back half of it and GIS 7 and 8 give
up a third and a sixth. The groin's trapping, face and leak correction are unchanged, sand is conserved, and the
spread is part of the groin's own contribution (footprint and leak correction). Calibration window 1996-2009 only,
strengths 0.06-0.19; the no-groin and runner's-groin runs of `../2026-10-20-tapered-fillet/` are the controls, on
the same settings (runner defaults after ecd6f0ff, dune growth pinned 0.005-0.01; `score` checks this).

    python spread_relaxation.py run --streams 4
    python spread_relaxation.py score     # -> spread_relaxation_scores.csv, spread_relaxation_domains.csv
    python spread_simple_gif.py --variant spread_0_23_13 --period test   # -> figures/spread_0_23_13_simple_test.gif
    python spread_simple_gif.py --variant spread_12_13_16                # -> figures/spread_12_13_16_simple_calibration.gif
    python spread_simple_gif.py --period calibration --current          # -> figures/spread_0_23_13_simple_calibration_current.gif (today's runner, as the test GIF)

Runs: `output/raw_runs/experiments/groin-wave-driven/2026-10-20-spread-relaxation/`.

## Results (1996-2009)

RMSE in metres. Gap = GIS 5 minus GIS 6 against annual CoastSat; net = end-of-window net change against the CoastSat
target (7-domain LOWESS, GIS 1-10 raw).

| layout | strength | gap | net GIS 2-89 | net GIS 2-12 | net GIS 4-9 |
|---|---|---|---|---|---|
| no groin | | 45.8 | 18.7 | 27.6 | 28.1 |
| runner (loss all in GIS 6) | 0.16 | 12.4 | 17.9 | 23.1 | 18.4 |
| spread (1/2, 1/3, 1/6 over GIS 6-8) | **0.13** | **12.1** | **17.5** | **21.1** | **12.8** |

End-of-window change, m, seaward +:

| | GIS 3 | 4 | 5 | 6 | 7 | 8 | 9 | 11 |
|---|---|---|---|---|---|---|---|---|
| observed | +12 | +24 | +9 | -9 | -26 | -16 | -14 | -43 |
| runner, 0.16 | -6 | +8 | +2 | -25 | +4 | +2 | +1 | -5 |
| spread, 0.13 | -11 | +7 | +9 | -23 | -8 | -10 | -2 | -12 |

## Findings

1. **The spread groin fits better than the runner's groin on every score**, most near the groin: GIS 4-9 RMSE
   12.8 against 18.4 m, GIS 2-12 21.1 against 23.1 m, with the gap no worse (12.1 against 12.4 m).
2. **It gives the three-domain northward effect in the record's direction.** GIS 7 and 8 now erode (-8 and -10 m
   against -26 and -16 m observed) instead of building out, and GIS 5 matches (+9 m).
3. **What it does not fix.** GIS 6 still loses 23 m against 9 m, GIS 7-9 erode about half as much as observed,
   and GIS 10-14 (18-37 m short of erosion) are beyond any groin's reach. GIS 1-3 end 4-6 m further landward than
   with the runner's groin, and GIS 11 7 m, small knock-on changes.
4. **Not yet checked:** the test window 2009-2025, and the 1/2, 1/3, 1/6 weights (fixed, not fitted).

## Fitting the weights (updrift side only)

Four more sets of shares on GIS 6, 7, 8, each at strengths 0.10-0.19 (run folders `wB`-`wE`, keys in `VARIANTS`;
the folder names are kept, the sets are named by their shares). Best strength per set (RMSE, m):

| shares on GIS 6, 7, 8 | strength | gap | net GIS 2-89 | net GIS 2-12 | net GIS 4-9 |
|---|---|---|---|---|---|
| runner (1, 0, 0) | 0.16 | 12.4 | 17.9 | 23.1 | 18.4 |
| 1/2, 1/3, 1/6 | 0.13 | 12.1 | 17.5 | 21.1 | 12.8 |
| 1/3, 1/3, 1/3 | 0.13 | 12.1 | 17.4 | 20.5 | 11.6 |
| 1/3, 1/2, 1/6 | 0.13 | 12.0 | 17.4 | 20.7 | 11.5 |
| 1/4, 1/2, 1/4 | 0.13 | 12.2 | 17.4 | 20.3 | 10.6 |
| **0, 2/3, 1/3** | **0.10** | **11.8** | **17.3** | **19.6** | **9.1** |

| end of 1996-2009, m | GIS 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|
| observed | +24 | +9 | -9 | -26 | -16 | -14 |
| 0, 2/3, 1/3 at 0.10 | +12 | +18 | -19 | -17 | -18 | -4 |
| 0, 2/3, 1/3 at 0.13 | +9 | +13 | -14 | -14 | -17 | -4 |
| 0, 2/3, 1/3 at 0.16 | +7 | +10 | -9 | -12 | -16 | -4 |

5. **The less of the loss GIS 6 takes, the better the fit; 0, 2/3, 1/3 is best on every score.** It is the advisor's
   wording taken literally (two-thirds one domain north, one-third two north): GIS 6 takes none of the loss across
   the face. At its best strength, 0.10, GIS 8 matches the record (-18 against -16), GIS 6 is 10 m short (-19
   against -9), GIS 7 erodes two-thirds as much as observed (-17 against -26) and GIS 5 gains 9 m too much. At
   0.16 GIS 6 and 8 match exactly, but the gap scores worse (18.2 m) because GIS 5 is then too high in the early
   years.
6. **A physical reading.** The groins all sit in GIS 6, so GIS 6 is the groin field itself: the sand between the
   groins is held, and the fillet loses sand from its updrift end (GIS 7-8) instead.
7. **Caveats.** 0, 2/3, 1/3 is at the edge of the grid (GIS 6 cannot take less than none without gaining), the shares were
   chosen on this window, and the test window 2009-2025 is the check.

## The spread groin (0, 2/3, 1/3) on the test window 2009-2025

`python spread_relaxation.py run --period 2009` / `score --period 2009` (-> `*_2009_2025.csv`). The spread groin at 0.06-0.16
and its own controls (no groin; the runner's groin at 0.08-0.16), groin off from 2017. The controls are run here
because the end rates were re-solved mid-study (c8092db1, 22:39: GIS 1 29.7 to 58.3 m/yr); the taper study's 2009
runs are on the old rates. `score` now also checks the end rates and the start shoreline.

| layout | strength | gap | gap 2010-16 | gap 2017-25 | net GIS 2-89 | net GIS 2-12 | net GIS 4-9 |
|---|---|---|---|---|---|---|---|
| no groin | | 23.3 | 34.7 | 5.1 | 28.5 | 52.7 | 40.3 |
| runner | 0.10 | 10.7 | 14.3 | 6.8 | 28.3 | 52.3 | 39.5 |
| runner (best gap) | 0.12 | 10.4 | 13.4 | 7.3 | 28.3 | 52.3 | 39.6 |
| **spread groin, calibration strength, unchanged** | **0.10** | 11.0 | 14.2 | 7.8 | 28.3 | **51.7** | **38.4** |
| spread groin (best gap) | 0.08 | 10.6 | 14.6 | 5.6 | 28.3 | 52.1 | 38.9 |

8. **The spread groin transfers: on the test window at its calibration strength it ties the runner's groin** on the gap
   (11.0 against 10.7 m; 14.2 against 14.3 m while the groin works, 2010-16) and is slightly better near the groin
   (GIS 4-9 38.4 against 39.5 m, GIS 2-12 51.7 against 52.3 m). It does not do worse anywhere that matters.
9. **The test window is a weak test of the groin.** The groin works only 2010-2016 and is off from the 2017 fill,
   and the window's large misses (GIS 7 +23 m against -37 m, GIS 2-4 too low by 50-110 m) are there with or
   without it: no groin scores 52.7 m on GIS 2-12 against the runner's 52.3. Before 2017 the spread groin behaves as in the
   calibration window (GIS 6 holds about 10 m better, GIS 7-8 lose 7-9 m more); after the groin is off the old
   step relaxes the same way in both, and by 2021 GIS 6 is where the runner's groin has it.

## Calibration rerun on today's runner (2026-10-21)

`python spread_relaxation.py run --current` / `score --current` (-> `*_current.csv`; runs in `current-<variant>/`).
The runner exactly as committed at 1d36cf20 (no dune-growth pin: growth 0.01-0.03; Dmaxel 12.2 m, slope 0.06; the
end rates and 1996 offset of c8092db1), with its own controls; all 8 runs share settings, end rates and start.

| layout | strength | gap | net GIS 2-89 | net GIS 2-12 | net GIS 4-9 |
|---|---|---|---|---|---|
| no groin | | 46.3 | 18.6 | 27.2 | 28.3 |
| runner (best) | 0.16 | 13.0 | 17.8 | 22.6 | 18.6 |
| **spread groin** | **0.10** | **12.2** | **17.2** | **18.4** | **8.9** |
| spread groin | 0.13 | 13.0 | 17.3 | 18.5 | 8.5 |

| end of 1996-2009, m | GIS 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|
| observed | +24 | +9 | -9 | -26 | -16 | -14 |
| runner, 0.16 | +8 | +1 | -29 | +3 | +2 | 0 |
| spread groin, 0.10 | +15 | +19 | -20 | -18 | -18 | -4 |

10. **On today's runner the spread groin still wins on every score, by a little more.** At 0.10 it beats the runner's groin at
    its best (0.16) on the gap (12.2 against 13.0 m), near the groin (8.9 against 18.6 m), on GIS 2-12 (18.4
    against 22.6 m) and island-wide (17.2 against 17.8 m). Its best strength is still 0.10, the strength that
    also transfers to the test window unchanged.

## Mirrored on the south side (`spread2`): worse, dropped

The same rule mirrored: GIS 5's coming gain across the face (BRIE's term at GIS 5) shared over GIS 5, 4, 3
(1/2, 1/3, 1/6), with the northward spread unchanged. Best at 0.06-0.10: gap 18.3-19.3 m, GIS 4-9 27.8-29.2 m,
GIS 2-12 26.2-27.0 m, all worse than the northward-only spread (12.1, 12.8, 21.1) and than the runner's groin;
worse still at higher strengths (0.19: gap 44.2 m).

| 0.10 | GIS 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| observed | +12 | +24 | +9 | -9 | -26 | -16 |
| both sides | +12 | +25 | -30 | -68 | -27 | -16 |

GIS 3-4 and 7-8 match, but the two domains at the face collapse. With GIS 5 passing half its gain on, it never
fills in, so the step at the face stays as steep as in 1996 and BRIE draws sand across it every year: GIS 6 keeps
losing, and GIS 5 loses too, since the groin's trapped sand also comes from it. The south side has to fill in next
to the groin for the step to relax at a sensible rate, so the spread stays on the north side only.
