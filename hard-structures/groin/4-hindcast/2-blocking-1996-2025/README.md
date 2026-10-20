# 2-blocking-1996-2025 — the blocking groin on the calibration and test periods

> **Superseded 2026-10-14** by the wave-driven groin (`../3-wave-driven-1967-2025/`), the runner's only groin since then. Kept as the record of how the blocking pin was fitted and of the two groin versions (drift, conserving) tried on the way to the wave-driven one. The wave-driven studies lived here until 2026-10-15.

**Code these studies ran against.** The fixed-shift (the old `GroinCallback`, the fixed-shift groin) and blocking (`BlockingGroinCallback`) groins that the scripts here use or patch were removed from `cascade/groin.py` on 2026-10-14, when the wave-driven groin became the runner's only groin. The scripts ran against the code at commit `4d937b79`; check that commit out to rerun one.

The blocking groin (`cascade.groin.BlockingGroinCallback`) fitted on the calibration period 1996-2009 and checked on the test period 2009-2025, with full management, the solved edgeBE ends and relocations off. Scored on the GIS 5|6 gap against the wet/dry photos and, from 10-08, the annual CoastSat gap.

| folder | question | status |
|---|---|---|
| `2026-10-05-blocking-fit-calibration/` | which b and f fit the 1996-2009 gap; how does that groin do on 2009-2025? Also the fixed-shift groin against the blocking groin | **b 0.6, f 0.6, pinned** (replaced 2026-10-14) |
| `2026-10-08-schedule-refit/` | does the fit change if the groin fails from 1996 (the last repair) instead of at the 2004 step? | done; its schedules were carried into the wave-driven fit (`../3-wave-driven-1967-2025/2026-10-09-corrected-groin/`) |
| `2026-10-08-conserving-groin/` | does a version whose own books balance still fit? Also the three-module comparison GIF | **stopped**: the premise was wrong (see `../../3-module-tests/2-straight-coast/`) |
| `2026-10-08-drift-groin/` | option 1: does a groin that blocks the drift the waves carry (BRIE untouched) fit 1996-2009, with what b and f? | grids done: joint floor 15.8-16.7 m at b 0.08-0.10 (pin 20.9); annual and photo scores want different b; 2009-2025 test of b 0.1 f 0.4: beats the pin, not no groin (gap held too wide). With a correction for BRIE's leak at its own bends it became the wave-driven groin |

After the 2017 Buxton fill every groin run misses the test period by 37-43 m. The cause is the fill's placement at GIS 6, not the groin (10-05 README).
