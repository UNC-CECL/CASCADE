# 1-straight-coast — version 1 against version 6

**Question.** Before fitting a groin to Buxton, does the module behave like a groin on a coast with nothing
else going on? Tested on straight coasts at orientations from −40° to +40°.

Trimmed 2026-10-12 to the two versions that matter: the first (version 1, the source/sink dipole) and the
current candidate (version 6, the corrected drift groin). Versions 2-5 and 7, the all-version figures and
this README as it stood with all seven are in `../archive/` (see its README for why each was dropped).

**Setup.** `1-run/straight_coast_test.py` is BRIE's alongshore step alone: its own `coast_diff` table, sparse
indices and clip at zero, under option A waves (Hs 2.0, Tp 7.5, asymmetry 0.6, high-angle 0.5). No storms,
no Barrier3D, no sea level. A straight coast at angle θ0, 61 domains of 500 m, one groin at the middle face,
50 years. One step of the emulator matches a real `Brie.update()` exactly. The script still runs all seven
versions; the retired ones write to `../archive/1-straight-coast/`.

| version | what it does each year |
|---|---|
| 1, `GroinCallback` (`1-source-sink-dipole/`) | moves a fixed M (12 m) of shoreline from the downdrift cell to the updrift cell, whatever the waves do |
| 6 (`6-trapping-corrected/`) | holds back f_block (0.6) of the wave climate's net drift arriving at the groin, and gives back what BRIE's row-scaled solve creates or destroys within 5 cells of it (10 in the Buxton hindcast); BRIE untouched |

## Layout

One folder per test (since 2026-10-15). Each holds its script(s) and what they write. Tests 2 to 6 import
`1-run/straight_coast_test.py`, and some read its tables, so run it first.

```
1-run/                      straight_coast_test.py: the 50-year runs of every version
    1-source-sink-dipole/       tables/ (by year; profiles with the updrift side at positive offsets),
    6-trapping-corrected/       figures/ (each version's profiles)
    comparison/tables/          the year-50 summary of all seven and the Pelnard-Considère check
2-solve-check/              straight_coast_conservative_solve_check.py: the drift groin under a face-centred
                            (conservative) solve; tables/conservative_solve_check.csv
3-sand-balance-at-a-step/   brie_sand_balance_figure.py: does BRIE's step conserve sand at a bend?
4-v1-vs-v6/                 the side-by-side figure and maps of versions 1 and 6, and their GIFs
5-dipole-explained/         dipole_explained_figure.py: how version 1 works and why its sand did not balance
6-corrected-explained/      how version 6 works: the explained figure, its sketch GIF and the 60-year GIF
```

## Figures

Each figure's Tests / How to read / Shows caption is in its `figures/supporting/CAPTIONS.md`.

| figure | question |
|---|---|
| `1-run/<version>/figures/straight_coast_<version>_profiles.png` | Does this version hold sand updrift? (six orientations, years 1 to 50) |
| `3-sand-balance-at-a-step/figures/brie_sand_balance_at_a_step.png` | Does BRIE's alongshore step keep the island's total sand constant when nothing is added or removed? (`brie_sand_balance_figure.py`) |
| `4-v1-vs-v6/figures/straight_coast_v1_vs_v6.png` | How do the two differ? Side by side: the rule, the shoreline at −20° and +10° over 50 years, each side at year 50, and the sand budget, against orientation (`straight_coast_v1_vs_v6_figure.py`) |
| `4-v1-vs-v6/figures/straight_coast_v1_vs_v6.gif` | The same, year by year, as maps, with each panel's running sand total (`straight_coast_v1_vs_v6_gif.py`; git-ignored, rerun to make it) |
| `4-v1-vs-v6/figures/straight_coast_v1_vs_v6_map.png` | What do the two do to the coast's shape? The −20° and +10° runs as true-scale, north-up maps: the coast at its own angle, the shoreline at years 10 and 50, the sand built out and cut back, the groin and the drift (`straight_coast_v1_vs_v6_map_figure.py`; reads the tables, reruns nothing) |
| `4-v1-vs-v6/figures/straight_coast_v1_vs_v6_map.gif` | The same maps year by year over 50 years, on seven coasts from −30° to +30°, dipole above and corrected below, with each map's running sand total and a drift arrow sized to that coast's drift (`straight_coast_v1_vs_v6_map_gif.py`, about a minute; the GIF and its year-50 preview PNG are git-ignored, rerun to make them) |
| `5-dipole-explained/figures/dipole_1_in_the_model_year.png` | Where did the dipole act in a model year? Barrier3D (with its background erosion) gives each domain its shoreline change, the dipole moves 12 m of it from GIS 5 to GIS 6, BRIE moves the result alongshore; with an example of the edit (`dipole_explained_figure.py`) |
| `5-dipole-explained/figures/dipole_2_explained.png` | How did the dipole work, and why did its sand not balance? A sketch of each year's steps, the fixed 12 m it moves against the sand BRIE's step makes or loses at the bend it creates, and the running total at −20° and +10° (`dipole_explained_figure.py`; reads the by-year table, reruns nothing). The companion of `corrected_groin_explained.png` |
| `5-dipole-explained/figures/dipole_3_steps.png` | Where in a CASCADE model year did the dipole act, step by step, from a straight coast? The four steps (Barrier3D, groin hook, BRIE, management), the shoreline after each in year 1, and the shape after 50 years of repeating (+539 m of sand made) (`dipole_steps_figure.py`; `--module corrected` draws the corrected groin's version) |
| `5-dipole-explained/figures/dipole_3_steps.gif` | The same, animated: year by year from a straight coast, the first five years step by step (after Barrier3D, the groin hook, BRIE) with the active step lit, then one frame a year to year 50, with the running sand total (`dipole_steps_gif.py`; `--module corrected` makes `6-corrected-explained/figures/corrected_groin_steps.gif`; git-ignored, rerun to make them) |
| `6-corrected-explained/figures/corrected_groin_explained.png` | How does the corrected groin work? The four steps each year, the sand held back against the sand BRIE carries back, and the gap levelling off, at −20° (`corrected_groin_explained_figure.py`); animated as `corrected_groin_explained.gif` (`corrected_groin_sketch_gif.py`, git-ignored) |
| `6-corrected-explained/figures/corrected_groin_steps.png` | The same as `dipole_3_steps.png` for the corrected groin: one model year from a straight coast (it holds back 58 m in year 1, read from the waves; its leak correction is zero until a bend forms), and the balanced shape after 50 years (+4 m of sand made). Drawn by `5-dipole-explained/dipole_steps_figure.py --module corrected` |
| `6-corrected-explained/figures/corrected_groin_steps.gif` | The corrected groin's version of `dipole_3_steps.gif` (`5-dipole-explained/dipole_steps_gif.py --module corrected`; git-ignored) |
| `6-corrected-explained/figures/corrected_groin_how_it_works.gif` | Version 6 alone at −20° for 60 years (`corrected_groin_gif.py`, git-ignored) |

## Findings

1. **BRIE never moves the drift.** It moves the shoreline with a diffusivity D(θ), built from the
   derivative of its transport formula (`_coast_diff`); the transport itself (`_coast_qs`) is used only for
   inlets. Version 6 rebuilds the net drift Q_net from `_coast_qs` and BRIE's wave pdf, which reproduces
   D = −(1/depth) dQ_net/dθ to within one 1° bin (depth 19.8 m). On a grid-aligned coast 0.91 M m³/yr moves
   south; it falls to zero near +32° and reverses beyond. Version 6 reads Q_net at the angle of the face
   approaching the groin from updrift.
2. **BRIE's alongshore solve does not conserve sand at a groin.** Each cell takes its diffusivity from its
   own forward face, so where neighbouring faces differ, as they do at the bend a groin makes, sand is
   created (coast turned against the drift) or destroyed (turned with it).
3. **Version 1 is mostly that error.** It moves only 12 m a year, but the bend it makes lets the solve add
   up to +752 m (−15°) or remove up to 1,481 m (+10°) of summed shoreline in 50 years. What it builds follows
   the coast angle, not the drift: at −20° both sides end seaward, at +10° both end landward.
4. **The error is the solve, not the groin.** `2-solve-check/straight_coast_conservative_solve_check.py` reruns the drift
   groin with a face-centred solve (each face's r shared by the two cells it joins): every orientation then
   conserves to 0.0 m, with 109 m of build-out at −30°, 166 m at 0°, 22 m at +20°
   (`2-solve-check/tables/conservative_solve_check.csv`).
5. **Version 6 matches it from outside BRIE.** It stays within 27 m of zero sand change at every orientation
   and builds the updrift side by the same amount it cuts the downdrift side, at every orientation up to the
   drift reversal: 108 m at −30°, 164 m at 0°, 23 m at +20°.
6. **Why version 6's gap levels off.** What it holds back hardly changes (58 m of shoreline in year 1, 57.6 m
   in year 60 at −20°). BRIE treats the groin face as open and carries sand back across it, about 15 m a
   year at first and 38 m by year 60. The gap is 104 m after year 1, 204 m at year 25 and 215 m at year 60.
7. **On the real coast** the correction also acts on Buxton's natural bend, changing the no-groin baseline
   (`../3-real-planform/7-corrected-groin/`). The Buxton hindcast therefore applies it to the groin's own
   bends only (option B, `../../3-hindcast/2-blocking-1996-2025/2026-10-09-corrected-groin/`).

Run: `python 1-run/straight_coast_test.py [--years 50]` (about a minute), then
`python 2-solve-check/straight_coast_conservative_solve_check.py` and the figure scripts in 3 to 6.
