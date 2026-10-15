# 2-module-explained — how each groin module works

Figures and GIFs that explain the groin modules step by step: where in a CASCADE model year each one acts,
what it does to the shoreline, and why its sand does or does not balance. Nothing here is a new test. The
straight-coast figures redraw runs of `../3-module-tests/2-straight-coast/1-run/straight_coast_test.py`
(BRIE's alongshore step alone, option A waves), which the scripts import from there. Split from the tests on
2026-10-15: `1-dipole/` and `2-corrected/` were `5-dipole-explained/` and `6-corrected-explained/` in
`2-straight-coast/`. The dipole's schematic (`groin_module_logic.png`, once `2-module-tests/figures/`) was
briefly here and is now with the dipole's Buxton fit in `../4-hindcast/1-dipole-1967-2017/figures/`.

```
1-dipole/           version 1, the source/sink dipole (GroinCallback, removed in 563a542e)
2-corrected/        version 6, the corrected groin (CorrectedGroinCallback, the runner's groin since 2026-10-13)
```

## Figures

Each figure's Tests / How to read / Shows caption is in its `figures/supporting/CAPTIONS.md`.

| figure | question |
|---|---|
| `1-dipole/figures/dipole_1_in_the_model_year.png` | Where did the dipole act in a model year? Barrier3D (with its background erosion) gives each domain its shoreline change, the dipole moves 12 m of it from GIS 5 to GIS 6, BRIE moves the result alongshore; with an example of the edit (`dipole_explained_figure.py`) |
| `1-dipole/figures/dipole_2_explained.png` | How did the dipole work, and why did its sand not balance? A sketch of each year's steps, the fixed 12 m it moves against the sand BRIE's step makes or loses at the bend it creates, and the running total at −20° and +10° (`dipole_explained_figure.py`; reads the by-year table, reruns nothing). The companion of `corrected_groin_explained.png` |
| `1-dipole/figures/dipole_3_steps.png` | Where in a CASCADE model year did the dipole act, step by step, from a straight coast? The four steps (Barrier3D, groin hook, BRIE, management), the shoreline after each in year 1, and the shape after 50 years of repeating (+539 m of sand made) (`dipole_steps_figure.py`; `--module corrected` draws the corrected groin's version) |
| `1-dipole/figures/dipole_3_steps.gif` | The same, animated: year by year from a straight coast, the first five years step by step (after Barrier3D, the groin hook, BRIE) with the active step lit, then one frame a year to year 50, with the running sand total (`dipole_steps_gif.py`; `--module corrected` makes `2-corrected/figures/corrected_groin_steps.gif`; git-ignored, rerun to make them) |
| `2-corrected/figures/corrected_groin_in_the_model_year.png` | Where does the corrected groin act in a model year, and what does it do to the numbers? The same layout as `dipole_1_in_the_model_year.png`: the three steps, and the edit in year 10 at −20° split into its two parts, 57 m held back from GIS 5 to GIS 6 and a leak correction summing to −28 m that cancels the +28 m BRIE's step would make at the groin's bends (`corrected_groin_in_the_model_year_figure.py`; reruns one coast) |
| `2-corrected/figures/corrected_groin_explained.png` | How does the corrected groin work? The four steps each year, the sand held back against the sand BRIE carries back, and the gap levelling off, at −20° (`corrected_groin_explained_figure.py`); animated as `corrected_groin_explained.gif` (`corrected_groin_sketch_gif.py`, git-ignored) |
| `2-corrected/figures/corrected_groin_steps.png` | The same as `dipole_3_steps.png` for the corrected groin: one model year from a straight coast (it holds back 58 m in year 1, read from the waves; its leak correction is zero until a bend forms), and the balanced shape after 50 years (+4 m of sand made). Drawn by `1-dipole/dipole_steps_figure.py --module corrected` |
| `2-corrected/figures/corrected_groin_steps.gif` | The corrected groin's version of `dipole_3_steps.gif` (`1-dipole/dipole_steps_gif.py --module corrected`; git-ignored) |
| `2-corrected/figures/corrected_groin_how_it_works.gif` | Version 6 alone at −20° for 60 years (`corrected_groin_gif.py`, git-ignored) |

Run `python ../3-module-tests/2-straight-coast/1-run/straight_coast_test.py` first if its tables are missing,
then any script here. Each takes under a minute.
