# Ends solved on the full window, 1996-2025

**Question.** Which GIS 1 and GIS 90 source/sink rates make the 30-yr run (1996 through 2025) match the CoastSat LRR 1996-2025 at its two ends?

**Setup.** `scripts/hatteras_ms/HAT_full_window_1996_2025.py solve`. The unchanged runner, with the 1996 period patched in-process: end 2025, 30 model years, storms `1996_2025_storms_v3_split12_trim24`, RSLR 0.005. Shoreline v2 start (DEM-centred), full management, no groin, no relocation. Secant from the 1996-2015 ends (+3.39 / +37.60). The target is the raw domain mean at GIS 1 and LOWESS-7 at GIS 90. Probes keep no `.npz`. Solver output is in `solve.txt`, and the run logs are in `output/logs/driver/full_window_1996_2025/`.

**Answer.** It converged at step 5: **GIS 1 +144.2227, GIS 90 +68.4160 m/yr**, residuals +0.009 / +0.005 m/yr (targets +3.99 / +1.46).

| step | GIS 1 imposed | GIS 1 model | GIS 90 imposed | GIS 90 model |
|---|---|---|---|---|
| 0 | +3.39 | +1.75 | +37.60 | -0.13 |
| 1 | +24.78 | +2.93 | +52.74 | +0.35 |
| 2 | +44.07 | +3.11 | +87.11 | +2.73 |
| 3 | +138.56 | +3.95 | +68.69 | +1.50 |
| 4 | +143.92 | +3.99 | +68.12 | +1.42 |
| 5 | +144.22 | +4.00 | +68.42 | +1.46 |

**Read with care.** GIS 1 responds weakly (local gain 0.038), so a +4 m/yr target takes +144 m/yr of imposed sand. The 2010-2026 solve did the same (+172.89). GIS 90 is about double its 1996-2015 value.

**Used by.** The final edgeBE run, `raw_runs/matrix/1996_2025/edgeBE/`. Not in the site config: HATTERAS_PERIODS has no 1996-2025 window.

**Status.** current, for the full-window run only.
