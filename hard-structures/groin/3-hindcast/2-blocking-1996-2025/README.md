# 2-blocking-1996-2025 — the blocking groin on the calibration and test periods

**Code these studies ran against.** The dipole (`GroinCallback`) and blocking (`BlockingGroinCallback`) groins that the scripts here use or patch were removed from `cascade/groin.py` on 2026-10-14, when the corrected groin became the runner's only groin. The scripts ran against the code at commit `4d937b79`; check that commit out to rerun one.

The blocking groin (`cascade.groin.BlockingGroinCallback`) fitted on the calibration period 1996-2009 and checked on the test period 2009-2025, with full management, the solved edgeBE ends and relocations off. Scored on the GIS 5|6 gap against the wet/dry photos and, from 10-08, the annual CoastSat gap.

| folder | question | status |
|---|---|---|
| `2026-10-05-blocking-fit-calibration/` | which b and f fit the 1996-2009 gap; how does that groin do on 2009-2025? Also the dipole against the blocking groin | **b 0.6, f 0.6, pinned** |
| `2026-10-08-schedule-refit/` | does the fit change if the groin fails from 1996 (the last repair) instead of at the 2004 step? | done; **decision pending** (options A-C in its README) |
| `2026-10-08-conserving-groin/` | does a version whose own books balance still fit? Also the three-module comparison GIF | **stopped**: the premise was wrong (see `../../2-module-tests/1-straight-coast/`) |
| `2026-10-08-drift-groin/` | option 1: does a groin that blocks the drift the waves carry (BRIE untouched) fit 1996-2009, with what b and f? | grids done: joint floor 15.8-16.7 m at b 0.08-0.10 (pin 20.9); annual and photo scores want different b; 2009-2025 test of b 0.1 f 0.4: beats the pin, not no groin (gap held too wide) |
| `2026-10-09-corrected-groin/` | the drift groin plus module 6's leak correction, applied to the groin's own bends only: does it fit, and does it conserve sand? | calibration joint 20.4 m at b 0.2 f 0.8 (pin 20.9), island +19 m (drift +107); 2009-2025 test worse than no groin: holds the gap (−2.8 m against −50.4 m). Schedules A-C fit 1996-2009 the same (19.7-21.1 m). With the groin off from the 2017 fill (D, a hypothesis to check against the fill record), the test joint is 18.5 m, against no groin's 27.6 |
| `2026-10-13-pregroin-planform/` | how strong is the corrected groin when working? From the 1967 shoreline (before the groins), and one continuous run 1967-2025 with no reset at 1996 or 2009 | **0.25** ramped in by 1978, held to 2009, 0.20 to 2016, off from the 2017 fill; follows the record to 17 m over 1967-2025. Supersedes the calendar-year window fits (0.16, 0.10) as the groin's strength: those describe a reset fillet relaxing |

After the 2017 Buxton fill every groin run misses the test period by 37-43 m. The cause is the fill's placement at GIS 6, not the groin (10-05 README).
