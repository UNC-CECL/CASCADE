# 2026-10-20-tapered-fillet — the groin's gain spread over GIS 6, 7 and 8

**Question** (advisor's suggestion): should the fillet taper north of the groin, instead of all landing in GIS 6?
Full effect against the groin, two-thirds of it one domain north, one-third two domains north. Scaled to sum to
one, that is 1/2 of the updrift gain in GIS 6, 1/3 in GIS 7 and 1/6 in GIS 8. The downdrift loss stays in GIS 5.

**Setup.** `tapered_fillet.py` patches the taper in-process, as `../2026-10-14-face-6-7/field_windows.py` does;
`cascade/groin.py` and the runner are not edited. The taper is applied to the step and to the groin's own
contribution, so the footprint and leak correction include it. Each window runs from its observed shoreline on
the runner's defaults after commit ecd6f0ff (edgeBE, full management, relocations off, storms `v3_240`, dune
ceiling Dmaxel 12.2 m, beach slope 0.06, off from 2017), with no groin, the runner's groin and the taper at the
same strengths. Dune growth is pinned to 0.005-0.01 by the script, because the runner's value changed three times
while the study ran; `score` refuses runs whose parameter yamls differ. 26 runs, all clean; the sand-balance guard
stopped none.

    python tapered_fillet.py run --streams 4
    python tapered_fillet.py score     # -> tapered_fillet_scores.csv, tapered_fillet_domains.csv
    python tapered_fillet.py figure    # -> figures/tapered_fillet_1996_2009_2009_2025.png
    python tapered_fillet.py positions # -> figures/tapered_fillet_plan_1996_2009_2009_2025.png (1:1 plan, GIS 4-8)
    python tapered_simple_gif.py       # -> figures/tapered_simple_{calibration,test}.gif (+ last frame .png)

The figure: top, end-of-window change at GIS 1-14 against CoastSat; bottom, the groin's own effect (each layout
less the no-groin run), at each layout's best gap strength.

Runs: `output/raw_runs/experiments/groin-wave-driven/2026-10-20-tapered-fillet/`. Logs: `logs/`. The first attempt, on mixed
settings (the runner changed mid-study), is kept in `../../../../../output/raw_runs/experiments/groin-wave-driven/
2026-10-20-tapered-fillet-mixed-settings/` and `logs-mixed-settings/`; its "any groin hurts 2009-2025" result was
the settings change, not the groin (`output/raw_runs/experiments/groin-wave-driven/2026-10-20-groin-2009-diagnosis/NOTE.md`).

## Results

Gap = the GIS 5 minus GIS 6 change against annual CoastSat. Net = end-of-window net change against the CoastSat
target (7-domain LOWESS, GIS 1-10 raw). RMSE in metres; best gap strength per layout.

| window | layout | strength | gap RMSE | net RMSE GIS 2-89 | net RMSE GIS 2-12 |
|---|---|---|---|---|---|
| 1996-2009 | no groin | | 45.8 | 18.7 | 27.6 |
| 1996-2009 | runner (all in GIS 6) | 0.16 | 12.4 | **17.9** | **23.1** |
| 1996-2009 | taper | 0.19 | **12.2** | 20.0 | 34.1 |
| 2009-2025 | no groin | | 23.3 | 28.1 | 54.6 |
| 2009-2025 | runner (all in GIS 6) | 0.12 | 10.3 | **28.1** | **54.8** |
| 2009-2025 | taper | 0.12 | **9.9** | 28.3 | 56.0 |

Net change per domain at the end of the window, metres, seaward positive (`tapered_fillet_domains.csv`):

| 1996-2009 | GIS 4 | GIS 5 | GIS 6 | GIS 7 | GIS 8 |
|---|---|---|---|---|---|
| observed | +24 | +9 | -9 | -26 | -16 |
| no groin | +17 | +22 | -71 | -12 | -1 |
| runner, 0.16 | +8 | +2 | -25 | +4 | +2 |
| taper, 0.19 | -4 | -21 | -57 | +22 | +28 |

## Findings

1. **The taper ties the runner on the gap but puts the sand in the wrong places.** Its gap score is about the
   same (12.2 against 12.4 m; 9.9 against 10.3 m), but near the groin (GIS 2-12) its net change is worse than the
   runner's in both windows, and in 1996-2009 worse than no groin (34.1 against 27.6 m).
2. **Why the gap score misses it.** The gap only scores GIS 5 minus GIS 6. With the taper both sides fall
   together: GIS 6 keeps half the sand and retreats 57 m, GIS 5 loses the deficit with less sand diffusing back
   from GIS 6 and retreats 21 m, and GIS 7-8 build out 22-28 m where the record has them retreating 16-26 m.
3. **Same reason as the 2026-10-14 field test.** From 1996 the coast north of the groin erodes in the record
   (the old fillet relaxing as the groins deteriorate). Any layout that sends sand further north works against
   that. The taper sends it further than the field did, and does worse.
4. **Not tested here: the build-up.** The record's 1967-1978 build-out peaks at GIS 7
   (`../2026-10-14-fillet-shape/`), so a taper may fit a run from 1967, or a forecast of a rebuilt groin field.
   The runner's layout is kept for the hindcast windows.
5. **In 2009-2025 the groin is neutral on net change** (54.8 against 54.6 m with no groin), as on the older
   settings. Near the groin the window fits poorly with or without it: GIS 4 ends 21-24 m landward where the
   record has it 32 m seaward. On the pre-ecd6f0ff dune settings the no-groin run had GIS 4 at +24 m and GIS 2-12
   at 37.9 m, so the new dune settings, not the groin, moved it (one run each, so a pointer, not a test).
