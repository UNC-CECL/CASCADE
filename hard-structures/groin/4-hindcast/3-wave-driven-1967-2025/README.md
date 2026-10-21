# 3-wave-driven-1967-2025 — the wave-driven groin fitted to Buxton

The wave-driven groin (`cascade.groin.GroinCallback`, once "corrected"; version 6 of `../../3-module-tests/`): the drift groin of `../2-blocking-1996-2025/2026-10-08-drift-groin/` plus a correction for BRIE's sand leak at the groin's own bends. The runner's only groin since 2026-10-14. Fitted on the calibration window 1996-2009 and the test window 2009-2025, and from 1967 on the groins' whole life, with full management, the solved edgeBE ends and relocations off. Scored on the GIS 5|6 gap against annual CoastSat and the wet/dry photos.

The studies before 2026-10-14 patched the groin in-process (`corrected_groin.py`); from 2026-10-14 they drive the unchanged runner.

| folder | question | status |
|---|---|---|
| `2026-10-09-corrected-groin/` | the drift groin plus the leak correction, applied to the groin's own bends only: does it fit, and does it conserve sand? Refit on the calendar-year setup 2026-10-12 | calendar refit: **0.16** for 1996-2009, **0.10** for 2009-2025, off from the 2017 fill (findings 18-22) |
| `2026-10-13-pregroin-planform/` | how strong is the groin when working? From the 1967 shoreline (before the groins), and one continuous run 1967-2025 with no reset at 1996 or 2009 | **0.25** from installation to 2009, 0.20-0.25 to 2016, off from 2017; follows the record to 17 m over 1967-2025. The groin's whole-life strength, for runs that start before the fillet exists (forecasts) |
| `2026-10-14-window-refit/` | what strength does each hindcast window need, run from its observed shoreline on today's runner? | **0.16** for 1996-2009 and **0.10** for 2009-2025 (0.25 overshoots from either start): the runner's values since 2026-10-14. `matrix_groin_effect_figure.py` draws the groin matrix against the observed net change |
| `2026-10-14-fillet-shape/` | is the model's fillet the right shape? The 1967-1978 build-up against the photos | the record's change is mostly downdrift (GIS 3-5) and its updrift build-out peaks at GIS 7; the model's sits at GIS 6 |
| `2026-10-14-face-6-7/` | where should the face be? The groins all sit in GIS 6 | a field across GIS 6 (half the sand kept in it) fits the build-up best (21.7 against 25.7 m), but the hindcast windows prefer the runner's 5\|6 face, which is kept |
| `2026-10-15-split-sample/` | run unchanged on the other window, does each window's strength still beat no groin? | the strength transfers: 0.16 run unchanged on 2009-2025 gives 14.1 m against no groin's 21.9 (skill 0.59); 0.10 on 1996-2009 gives 19.5 against 47.9 (skill 0.84) |

After the 2017 Buxton fill the groin is off; the fill's placement at GIS 6 makes its own gap (`../2-blocking-1996-2025/2026-10-05-blocking-fit-calibration/`).

Run output is under `output/raw_runs/experiments/groin-wave-driven/<study>/` (`groin/<study>/` before 2026-10-20). Run metadata written before 2026-10-15 still names these folders under `2-blocking-1996-2025/`.
