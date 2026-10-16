# 2026-10-13-pregroin-planform

**Question.** How strong is the wave-driven groin when it is working, and does one strength history reproduce the groin's whole life, 1967-2025, without resetting to the observed shoreline at 1996 and 2009? The answer sets the strength a forecast of the 2026 rebuild should use.

**Why.** The calendar-year refit (`../2026-10-09-corrected-groin/`, findings 18-22) fitted 0.16 on 1996-2009 and 0.10 on 2009-2016. But each of those runs starts from the observed shoreline of its window, which already holds the fillet the groins built. A strength fitted that way says how fast the model relaxes a fillet it was handed, not how strong the groin is. The only years in which the groin builds its own fillet are before 1996.

## Why there is a jump at GIS 5|6 in the no-groin runs

The 1996 and 2009 start offsets step seaward 264 m and 240 m from GIS 5 to GIS 6, against 60-90 m between other domains. The step is in the CoastSat shoreline itself: re-measured, it is 63-67 m per 100 m transect across the groin field. The 1967 dune line, two years before the groins, steps 104 m. So about 100 m is the coast's natural bend toward the cape and about 150 m is the fillet the groins built. With no groin, BRIE's alongshore transport smooths that step: GIS 6 erodes and GIS 5 builds out, about 9 m net over GIS 4-7 by 2009. It is not a source or a sink.

| line | GIS 4→5 | GIS 5→6 | GIS 6→7 |
|---|---|---|---|
| 1967 dune line (before the groins) | 132 m | **104 m** | 67 m |
| 1978 dune line | 134 m | 201 m | 63 m |
| 1997 dune line | 139 m | 256 m | 98 m |
| 1996 shoreline (model start) | 145 m | **264 m** | 84 m |
| 2009 shoreline (model start) | 131 m | 240 m | 71 m |

## Setup

- **Start (`pregroin_planform.py build`).** The 1967 wet/dry shoreline (August 1967, `../../../1-observations/data/gis/wet_dry_groin.geojson`), measured on the model's 100 m transects with `step1_line_to_raw_offsets.intersect`, the code every offset uses. Re-measuring the 1996 shoreline that way reproduces its raw offsets to 0 m. The 1967 line covers GIS 2-12 only, so the start is the 1996 shoreline offset (shoreline/v4) everywhere else, and over GIS 2-12 the 1996 values plus (1967 − 1996) less that difference's straight line between GIS 2 and 12, so the two join without a step. The GIS 5→6 step is 108 m. File: `inputs/Island_Shoreline_Offsets_1996_with_1967_GIS2-12_PADDED_120.csv`.
- **Everything else** is the calibration run's: the 1996 topography, full management, the solved edgeBE ends, relocations off, the corrected groin with the tracked reference. Each child process swaps the period's offset file (and for the long runs its years, storms, sea level and fills) and drives the unchanged runner, as `corrected_groin.py` does. The runner's report, labels and metadata still show the 1996 period's; each run folder has a note saying what was patched.
- **Storms.** 1968-1983 are before the storm record, so they replay the recorded 1984-1999 (Hannah's choice, 2026-10-13); 1984-1995 come from the 1984-2004 series and 1996-2025 from the 1996-2025 series the calibration and test windows were sliced from (the two agree exactly where they overlap). The old dipole rig's 1967 storms were not reused: they come from the superseded 72 h variant.
- **Scored** on the gap across the groins (GIS 6 seaward of GIS 5) as change since 1967: the wet/dry photo dates, and annual CoastSat placed on the photos by its mean offset over the years both have (2004 left out, where the photos sit 42 m above CoastSat). Scoring the level, not only the change since 1996, is what shows whether the eras join.

## 1. The groin on the pre-groin coast, 13 years (`run`)

The calibration's storms (1997-2009) on the 1967 start, constant strength from year 0.

| strength | gap growth, year 2 | year 6 | year 13 |
|---|---|---|---|
| no groin | 0 | 1 | 1 m |
| 0.10 | 31 | 47 | 50 m |
| 0.16 | 49 | 71 | 77 m |
| 0.25 | 74 | 103 | 106 m |
| observed 1967-1978 | | | 120 m (11 years) |

1. **With no kink to smooth, no groin does nothing**, which confirms the explanation above.
2. **The groin builds a fillet that levels off** in about 6 years, where its trapping balances BRIE's smoothing. The model splits the change about evenly (GIS 5 retreats, GIS 6 builds out); the observed 1978 line is mostly downdrift retreat, spread over GIS 3-5.

## 2. The maintained life, 1967-1995 (`run --full`)

| strength | 1970 | 1978 | 1985 | 1993 | 1995 | photo RMSE (without 1995) |
|---|---|---|---|---|---|---|
| observed (photos) | 9 | 117 | 99 | 116 | 155 | |
| 0.16 | 49 | 74 | 79 | 79 | 79 | 31 m |
| 0.20 | 61 | 90 | 92 | 97 | 97 | **27 m** |
| 0.25 | 74 | 108 | 108 | 110 | 109 | 32 m |
| 0.30 | 89 | 123 | 121 | 126 | 128 | 42 m |

3. **The plateau (1985-1993, about 106 m) wants 0.25**; the overall minimum sits at 0.20 only because the model builds too fast in the first years. 1995 is the Hurricane Gordon spike, mostly storm loss on the south side.

## 3. One continuous run, 1967-2025 (`continuous_life.py`)

Fills on (Rodanthe 2014, Buxton 2017 and 2022, Avon 2022; the 1971 and 1973 Buxton fills are not in the model's list), RSLR 0.004 throughout.

4. **The window fits do not join up.** Run as one history (0.25 to 1995, then 0.16, 0.10, off), the gap drops about 25 m when the strength falls in 1996 and stays 30-40 m below the record to 2017 (`continuous_scores.csv`). Scored only as change since 1996 the same run looks good, because its shape is right; its level is not.
5. **Fitted era by era on the continuous run** (`continuous_life.py fit`, each era held while the next is fitted, scored on its level; `continuous_fit.csv`):

| era | strength | level RMSE | neighbours |
|---|---|---|---|
| 1969-1995 (scored 1978-1995) | **0.25** | 18.3 m | 0.20: 25.3; 0.30: 24.5 |
| 1996-2009 | **0.25** (the cap: no stronger than before) | 18.1 m | 0.20: 28.0; 0.15: 41.9 |
| 2010-2016 | **0.20** | 9.4 m | 0.25: 9.9; 0.15: 22.1 |
| 2017-2025 | off | 15.2 m | |

6. **Follow-ups** (`continuous_life.py followups`, `continuous_followups.csv`):
   - **The 1995 repair:** 0.30 for 1996-2009 scores 17.2 m against 18.1 m for 0.25, and 0.35 scores 26.5 m. The groin did not weaken after 1995; it may have been slightly stronger.
   - **The pre-1984 storms:** with 1968-1983 drawn at random from the recorded years instead (two seeds), 0.25 is still the minimum for 1978-1995 (18.7 m both), with 0.20 and 0.30 worse by 5-7 m.
   - **The build-up:** ramping the strength in linearly to 0.25 by 1978 brings the 1970-1977 error from 56.5 m to 7.6 m and leaves the plateau unchanged (18.8 m against 18.3 m).

**The groin's strength history, fitted on one continuous run:** a ramp to 0.25 by 1978, 0.25 to 2009 (perhaps 0.30 after the 1995 repair), 0.20 to 2016, off from the 2017 fill. It follows the record to 16.6 m over 1967-2025: 7.6 m in the build-up, 18.8 m for 1978-1995, 18.3 m for 1996-2009, 9.8 m for 2010-2016, 15.9 m after the fill (`figures/continuous_gap_1967_2025.png`).

**What it means.**
- The windowed strengths, 0.16 on 1996-2009 and 0.10 on 2009-2016 (`../2026-10-09-corrected-groin/`), are superseded as estimates of the groin's strength: they were the rate a reset fillet relaxes at. The groin held about 0.25 for some 40 years and largely stopped working only around the 2017 fill, which matches the losses reported from 2015-2016 and the record's gap, which held near 120-150 m from 1995 to 2009.
- **For a forecast of the 2026 rebuild:** about 0.25, ramped in over 5-9 years. Only the southernmost of the three groins is being rebuilt, so 0.25 is an upper bound for it.
- **For the 1996-2009 and 2009-2025 hindcast runs**, which start from the observed fillet, a strength of about 0.25 to 2009 and 0.20-0.25 to 2016 is the one consistent with the groin's whole life. The window fits remain the better fit inside each reset window.

**Caveats.** The 1996 topography from 1967; one RSLR for 58 years; pre-1984 storms are stand-ins (tested above); the 1967 start covers GIS 2-12 only; the model's fillet sits in the two domains beside the groin, where the observed one spans three or four; the no-groin run jumps 35 m at 2017 because the fill starts at GIS 6, so the post-2017 years carry the fill's footprint as much as the groin.

## Files

| file | what |
|---|---|
| `pregroin_planform.py` | `build` the 1967 start; `run` the 13-year test; `storms` and `run --full` the 1968-1995 life; `score`; `gif`, `gif-full` |
| `continuous_life.py` | `storms`, `run`, `score`, `fit`, `followups`, `figure`, `gif` for 1967-2025 |
| `inputs/` | the start offset and its table, the 1967 transect crossings, the storm files (1968-1995; 1968-2025; the two resampled variants) |
| `full_life_scores.csv`, `continuous_scores.csv`, `continuous_fit.csv`, `continuous_fit_best.txt`, `continuous_followups.csv` | the scores above |
| `figures/continuous_gap_1967_2025.png` | the gap 1967-2025: no groin, the window fits run continuously, the fitted history; the strength each uses (caption in `figures/supporting/CAPTIONS.md`) |
| `figures/*.gif` | stacked plan views, from `../2026-10-09-corrected-groin/groin_panels_gif.py`'s `render`; git-ignored, rerun to make them |

Runs are under `output/raw_runs/experiments/groin/2026-10-13-pregroin-planform/`, `2026-10-13-pregroin-full-life/`, `2026-10-13-continuous-1967-2025/` and its `-seed1`, `-seed2` variants.
