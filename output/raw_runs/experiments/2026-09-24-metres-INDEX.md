# 2026-09-24: the island offset in metres — the chain

Three studies made on one day, each answering the question the one before it
raised. Numbered in the order they ran; each folder's README or NOTE is the
record, and each links back here.

| step | folder | the question | the answer | decision |
|---|---|---|---|---|
| 1 | [`2026-09-24-metres-1-offset-units/`](2026-09-24-metres-1-offset-units/README.md) | Every calibrated run handed BRIE the island offset ÷10, though BRIE and the offset file are both metres. Can the full-scale offset be tuned (Hs, wave angles) to match? | Metres reaches RMSE 1.07 against ÷10's 1.05 once the high-angle fraction is tuned; trend removal never beats a flat line. Both are near a flat line (1.17). | **Metres, no trend removal** — the default since 2026-09-24; offset builds re-padded and renumbered v1 |
| 2 | [`2026-09-24-metres-2-wave-sensitivity/`](2026-09-24-metres-2-wave-sensitivity/README.md) | In metres, how does each wave parameter move the model, natural and under full management, in 1996–2010 and 2010–2024? | 1996–2010: ridges (larger Hs needs a larger high-angle fraction and a longer period); best natural +18% (Hs 1.25, Tp 10, asym 0.8, high-angle 0.45), best managed +17% (baseline). 2010–2024: no setting reaches a flat line (bias −1.7 to −7 m/yr against observed gain) — a window/setup problem, not the waves. | None yet: 2010–2024 needs its own look (fills, the 2021 CoastSat step, storms) |
| 3 | [`2026-09-24-metres-3-barrier3d-overwash-fix/`](2026-09-24-metres-3-barrier3d-overwash-fix/NOTE.md) | Why did 10 of step 2's runs crash with no error? | An upstream Barrier3D bug: `route_overwash` indexes row and column swapped (line 1092) — wrong cells in every run, out of bounds on narrow domains. The fix moves scores only in the third decimal and ends the crashes. | **Fix adopted**: Barrier3D branch `fix/route-overwash-axis-swap` (local) checked out; every run records its Barrier3D; the crashed cells re-run |

Before these: `2026-09-22-shoreline-offset/` tested a shoreline-derived offset
at ÷10 scale, which is what exposed the units question. After them: the ÷10
production runs were archived in `../archive/2026-09-24-pre-metres/`.

## Scripts

All in `scripts/hatteras_ms/experiments/`, named for their step:

| step | runs and scores | figures |
|---|---|---|
| 1 | `HAT_metres_1_offset_units.py` | `HAT_metres_1_offset_units_plot.py` |
| 2 | `HAT_metres_2_wave_sensitivity.py` | `HAT_metres_2_wave_sensitivity_plot.py` |
| 3 | `HAT_metres_3_overwash_fix.py` | `HAT_metres_3_overwash_fix_plot_explained.py` |

## Housekeeping

Reorganised 2026-09-24 (Hannah) from three places into this chain: step 1 was
`experiments/2026-09-24-island-offset-scale-wave-tuning/`, step 2
`sensitivity/2026-09-24-natural-waves/` (kind `sensitivity`, now `experiment`:
it is a multi-stage study, not a one-at-a-time sweep of the calibration), step
3 `experiments/2026-09-24-overwash-fix/`. Every run's recorded tag and kind was
rewritten to its new place and the index rebuilt; the Barrier3D diagnostics
and the upstream issue draft moved from step 2 to step 3; the aborted
`2026-09-23-offset-metres/` (one parameters file, no results) was removed.
The run folders of steps 1–3 are on disk only (paths past Windows' 260); each
study's README says what is committed.
