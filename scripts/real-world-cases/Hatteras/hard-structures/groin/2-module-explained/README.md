# 2-module-explained — how each groin module works

Figures and GIFs that explain the groin modules step by step: where in a CASCADE model year each one acts,
what it does to the shoreline, and why its sand does or does not balance. Nothing here is a new test. The
straight-coast figures redraw runs of `../3-module-tests/2-straight-coast/1-run/straight_coast_test.py`
(BRIE's alongshore step alone, option A waves), which the scripts import from there. Split from the tests on
2026-10-15: `1-fixed-shift/` and `2-wave-driven/` were `5-dipole-explained/` and `6-corrected-explained/` in
`2-straight-coast/`. The fixed-shift groin's schematic (`groin_module_logic.png`, once `2-module-tests/figures/`) was
briefly here and is now with the fixed-shift groin's Buxton fit in `../4-hindcast/1-dipole-1967-2017/figures/`.

```
1-fixed-shift/      the fixed-shift groin: a fixed amount moved across the groin every year (once 'the
                    source/sink dipole', straight-coast version 1; the old dipole groin, removed in 563a542e)
2-wave-driven/      the wave-driven groin: a share of the drift the waves bring held back, sand balanced (once
                    'the corrected groin', straight-coast version 6; GroinManager, the runner's
                    groin since 2026-10-13)
3-brie-sand-balance/ why BRIE's alongshore step makes or loses sand at a sharp bend, from unmodified
                    BRIE on several coastline shapes: the reason the wave-driven groin needs its leak correction
```

## Figures

Each figure's Tests / How to read / Shows caption is in its `figures/supporting/CAPTIONS.md`.

| figure | question |
|---|---|
| `1-fixed-shift/figures/fixed_shift_1_in_the_model_year.png` | Where did the fixed-shift groin act in a model year? Barrier3D (with its background erosion) gives each domain its shoreline change, the fixed-shift groin moves 12 m of it from GIS 5 to GIS 6, BRIE moves the result alongshore; with an example of the edit (`fixed_shift_explained_figure.py`) |
| `1-fixed-shift/figures/fixed_shift_2_explained.png` | How did the fixed-shift groin work, and why did its sand not balance? A sketch of each year's steps, the fixed 12 m it moves against the sand BRIE's step makes or loses at the bend it creates, and the running total at −20° and +10° (`fixed_shift_explained_figure.py`; reads the by-year table, reruns nothing). The companion of `wave_driven_2_explained.png` |
| `1-fixed-shift/figures/fixed_shift_3_steps.png` | Where in a CASCADE model year did the fixed-shift groin act, step by step, from a straight coast? The four steps (Barrier3D, groin, BRIE, management), the shoreline after each in year 1, and the shape after 50 years of repeating (+539 m of sand made) (`groin_steps_figure.py`; `--module wave-driven` draws the wave-driven groin's version) |
| `1-fixed-shift/figures/fixed_shift_3_steps.gif` | The same, animated: year by year from a straight coast, the first five years step by step (after Barrier3D, the groin, BRIE) with the active step lit, then one frame a year to year 50, with the running sand total (`groin_steps_gif.py`; `--module wave-driven` makes `2-wave-driven/figures/wave_driven_3_steps.gif`; git-ignored, rerun to make them) |
| `2-wave-driven/figures/wave_driven_1_in_the_model_year.png` | Where does the wave-driven groin act in a model year, and what does it do to the numbers? The same layout as `fixed_shift_1_in_the_model_year.png`: the three steps, and the edit in year 10 at −20° split into its two parts, 57 m held back from GIS 5 to GIS 6 and a leak correction summing to −28 m that cancels the +28 m BRIE's step would make at the groin's bends (`wave_driven_in_the_model_year_figure.py`; reruns one coast) |
| `2-wave-driven/figures/wave_driven_2_explained.png` | How does the wave-driven groin work? The four steps each year, the sand held back against the sand BRIE carries back, and the gap levelling off, at −20° (`wave_driven_explained_figure.py`); animated as `wave_driven_2_explained.gif` (`wave_driven_sketch_gif.py`, git-ignored) |
| `2-wave-driven/figures/wave_driven_3_steps.png` | The same as `fixed_shift_3_steps.png` for the wave-driven groin: one model year from a straight coast (it holds back 58 m in year 1, read from the waves; its leak correction is zero until a bend forms), and the balanced shape after 50 years (+4 m of sand made). Drawn by `1-fixed-shift/groin_steps_figure.py --module wave-driven` |
| `2-wave-driven/figures/wave_driven_3_steps.gif` | The wave-driven groin's version of `fixed_shift_3_steps.gif` (`1-fixed-shift/groin_steps_gif.py --module wave-driven`; git-ignored) |
| `2-wave-driven/figures/wave_driven_60_years.gif` | The wave-driven groin alone at −20° for 60 years (`wave_driven_60_years_gif.py`, git-ignored) |

Run `python ../3-module-tests/2-straight-coast/1-run/straight_coast_test.py` first if its tables are missing,
then any script here. Each takes under a minute.
