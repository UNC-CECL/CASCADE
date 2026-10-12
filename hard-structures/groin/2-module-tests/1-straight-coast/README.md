# 1-straight-coast — version 1 against version 6

**Question.** Before fitting a groin to Buxton, does the module behave like a groin on a coast with nothing
else going on? Tested on straight coasts at orientations from −40° to +40°.

Trimmed 2026-10-12 to the two versions that matter: the first (version 1, the source/sink dipole) and the
current candidate (version 6, the corrected drift groin). Versions 2-5 and 7, the all-version figures and
this README as it stood with all seven are in `../archive/` (see its README for why each was dropped).

**Setup.** `straight_coast_test.py` is BRIE's alongshore step alone: its own `coast_diff` table, sparse
indices and clip at zero, under option A waves (Hs 2.0, Tp 7.5, asymmetry 0.6, high-angle 0.5). No storms,
no Barrier3D, no sea level. A straight coast at angle θ0, 61 domains of 500 m, one groin at the middle face,
50 years. One step of the emulator matches a real `Brie.update()` exactly. The script still runs all seven
versions; the retired ones write to `../archive/1-straight-coast/`.

| folder | version | what it does each year |
|---|---|---|
| `1-source-sink-dipole/` | 1, `GroinCallback` | moves a fixed M (12 m) of shoreline from the downdrift cell to the updrift cell, whatever the waves do |
| `6-trapping-corrected/` | 6 | holds back f_block (0.6) of the wave climate's net drift arriving at the groin, and gives back what BRIE's row-scaled solve creates or destroys within 5 cells of it (10 in the Buxton hindcast); BRIE untouched |

Each folder holds `tables/` (by year, and profiles with the updrift side at positive offsets) and `figures/`.
`comparison/tables/` holds the year-50 summary of all seven, the Pelnard-Considère check and the
conservative-solve check.

## Figures

Each figure's Tests / How to read / Shows caption is in its `figures/supporting/CAPTIONS.md`.

| figure | question |
|---|---|
| `comparison/figures/straight_coast_v1_vs_v6.png` | How do the two differ? Side by side: the rule, the shoreline at −20° and +10° over 50 years, each side at year 50, and the sand budget, against orientation (`straight_coast_v1_vs_v6_figure.py`) |
| `comparison/figures/straight_coast_v1_vs_v6.gif` | The same, year by year, as maps, with each panel's running sand total (`straight_coast_v1_vs_v6_gif.py`; git-ignored, rerun to make it) |
| `<version>/figures/straight_coast_<version>_profiles.png` | Does this version hold sand updrift? (six orientations, years 1 to 50) |
| `6-trapping-corrected/figures/corrected_groin_explained.png` | How does the corrected groin work? The four steps each year, the sand held back against the sand BRIE carries back, and the gap levelling off, at −20° (`corrected_groin_explained_figure.py`); animated as `corrected_groin_explained.gif` (`corrected_groin_sketch_gif.py`, git-ignored) |
| `6-trapping-corrected/figures/corrected_groin_how_it_works.gif` | Version 6 alone at −20° for 60 years (`corrected_groin_gif.py`, git-ignored) |

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
4. **The error is the solve, not the groin.** `straight_coast_conservative_solve_check.py` reruns the drift
   groin with a face-centred solve (each face's r shared by the two cells it joins): every orientation then
   conserves to 0.0 m, with 109 m of build-out at −30°, 166 m at 0°, 22 m at +20°
   (`comparison/tables/conservative_solve_check.csv`).
5. **Version 6 matches it from outside BRIE.** It stays within 27 m of zero sand change at every orientation
   and builds the updrift side by the same amount it cuts the downdrift side, at every orientation up to the
   drift reversal: 108 m at −30°, 164 m at 0°, 23 m at +20°.
6. **Why version 6's gap levels off.** What it holds back hardly changes (58 m of shoreline in year 1, 57.6 m
   in year 60 at −20°). BRIE treats the groin face as open and carries sand back across it, about 15 m a
   year at first and 38 m by year 60. The gap is 104 m after year 1, 204 m at year 25 and 215 m at year 60.
7. **On the real coast** the correction also acts on Buxton's natural bend, changing the no-groin baseline
   (`../3-real-planform/corrected_groin/`). The Buxton hindcast therefore applies it to the groin's own
   bends only (option B, `../../3-hindcast/2-blocking-1996-2025/2026-10-09-corrected-groin/`).

Run: `python straight_coast_test.py [--years 50]` (about a minute), then
`python straight_coast_conservative_solve_check.py` and the figure scripts.
