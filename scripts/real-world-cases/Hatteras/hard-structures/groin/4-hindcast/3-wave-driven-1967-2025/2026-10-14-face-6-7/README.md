# Where should the groin's face be? (2026-10-14)

**Question.** All three Buxton groins sit inside GIS 6: on transects 37, 39 and 41, its south edge, middle and north edge, about 200 m apart. On every photo from 1978 to 2016 the updrift build-out peaks at GIS 7, not GIS 6 (`../2026-10-14-fillet-shape/`). The runner's groin is one face at GIS 5|6 that puts the trapped sand in GIS 6. Does a different placement fit better, first in the build-up and then in the hindcast windows?

**Placements**, patched in-process by `field_windows.patch_groin` (`cascade/groin.py` and the runner are not edited):
- `face56`: the runner's groin as it is, the control.
- `face67`: the same groin one domain north. It takes sand from GIS 6 and gives it to GIS 7.
- `field`: domain 6 is the field. It takes sand from GIS 5 (south of the field) and gives it to GIS 7 (north of it).
- `split50` and `split33`: the field, with half or a third of the trapped sand kept in GIS 6 and the rest given to GIS 7.

    python face_experiment.py run | check | score | figure   # the 1967-1995 build-up
    python field_windows.py run --share 0.5                  # the field, half kept, in both hindcast windows
    python field_windows.py score --share 0.5                # -> field_windows_scores_share0.50.csv

## 1. The build-up, 1967-1995

**Setup.** As the "full life" runs of `../2026-10-13-pregroin-planform/`:
- the 1967 wet/dry shoreline at GIS 2-12;
- 1968-1995, with the 1968-1983 storms replaying 1984-1999;
- RSLR 0.004, no fills, the 1996 topography;
- edgeBE and full management, strength constant from 1969.

**Scoring.** Each line's change since 1967 at GIS 2-12 is compared with the 10 wet/dry photo dates from 1970 to 1995, after taking off its straight line between the end domains.

**Check.** The control reproduces the 10-13 full-life runs to 0.000 m at strengths 0 and 0.25. The 10-13 runs at 0.30 and 0.35, which used the same setup, are scored here as the control.

RMSE against the photos (m), all 10 dates. The best strength is shown for each placement; every cell is in `face_scores.csv`.

| placement | strength | shape, GIS 2-12 | downdrift, GIS 3-5 | updrift, GIS 6-9 | updrift peak |
|---|---|---|---|---|---|
| no groin | 0 | 43.0 | 51.3 | 49.3 | GIS 8 |
| face 5\|6 (runner) | 0.30 (0.25: 26.7) | 25.7 | 27.2 | 31.6 | GIS 6 |
| face 6\|7 | 0.16 | 40.6 | 41.7 | 51.0 | GIS 7 |
| field 5 to 7 | 0.16 | 25.6 | 27.3 | 31.2 | GIS 7 |
| field, a third kept in GIS 6 | 0.20 | 23.2 | 27.2 | 26.5 | GIS 7 |
| **field, half kept in GIS 6** | **0.25** (0.20: 22.6; 0.30: 22.8) | **21.7** | 26.3 | 23.9 | GIS 7 |

The observed updrift peak is at GIS 7 on most dates.

1. **The 6|7 face is wrong.** It puts the downdrift deficit in GIS 6, the field, which retreats 67-88 m where the photos show +7 to +12 m. It is worse than the control at every strength.
2. **The field puts the sand where the photos have it**: +63 m at GIS 7 by 1993 against +58 m observed. With all the sand given to GIS 7 it only ties the control (25.6 against 25.7 m), because GIS 6 then retreats 28-42 m.
3. **Keeping half the sand in GIS 6 fits best**: 21.7 m, against 25.7 m for the runner's groin. GIS 6 retreats only 6-17 m (observed +7 to +12 m) and GIS 7 stays right. The updrift error drops from 31.6 to 23.9 m.
4. **No placement fixes:**
   - the depth of the early downdrift trough (1978, GIS 3-5: -43/-129/-114 m observed);
   - the build-out's spread to GIS 8 (+44 m in 1993; about +16 m at most in the model);
   - the too-fast build (the best field run has +55 m at GIS 7 by 1978, against +28 m observed).

## 2. The hindcast windows

**Setup.** The field with half kept in GIS 6, run from each window's observed shoreline as `../2026-10-14-window-refit/` (edgeBE, full management, no relocation, off from 2017), against that refit's runs of the runner's groin. Two scores:
- the yearly gap across the groins (GIS 5 minus GIS 6) against annual CoastSat, which chose the runner's window strengths;
- net change against the CoastSat net-change target (7-domain LOWESS, GIS 1-10 raw), over the whole island and near the groin.

| window | placement | strength | gap, annual (m) | net change GIS 2-12 (m) | net change GIS 2-89 (m) |
|---|---|---|---|---|---|
| 1996-2009 | no groin | - | 47.9 | 25.5 | 19.6 |
| 1996-2009 | face 5\|6 (runner) | 0.16 | 13.4 | **21.7** | 19.0 |
| 1996-2009 | field, half kept | 0.16 / 0.19 | 12.4 / **12.3** | 27.0 / 28.7 | 19.8 / 20.1 |
| 2009-2025 | no groin | - | 21.9 | 37.5 | 23.6 |
| 2009-2025 | face 5\|6 (runner) | 0.10 | **9.7** | 38.3 | 23.8 |
| 2009-2025 | field, half kept | 0.12 | 10.7 | 39.4 | 24.0 |

Net change 1996-2009 near the groin (m, seaward +):

| | GIS 5 | 6 | 7 | 8 |
|---|---|---|---|---|
| observed | +9 | -9 | -26 | -16 |
| face 5\|6, 0.16 | +18 | -18 | +5 | +2 |
| field, half kept, 0.16 | +7 | -33 | +30 | +6 |

5. **The windows do not support the field.** It ties the runner's groin on the gap across the groins. On net change it is worse, near the groin than even no groin in 1996-2009, and about 1 m worse in 2009-2025.
6. **The reason is the period, not the placement.** In the build-up the groins were working and sand built up north of the field. From 1996 GIS 6-14 erode in the record (the old fillet relaxing as the groins deteriorate), and the field keeps feeding GIS 7: +30 m against -26 m observed. The runner's groin keeps its sand in GIS 6, which does less harm while the coast north of it erodes.
7. In 2009-2025 the two placements are within about 1 m everywhere: the groin is weak (0.10) and off from 2017. The misfit at GIS 6-9 (-35 m observed, about +27 m modelled) is not the groin's.

## 3. The build-up scored on the gap across each boundary

    python face_gap_scores.py   # -> face_gap_scores.csv

Section 1 scores the shape of the line. This scores the gap across each candidate face instead: the north domain's change since 1967 minus the south domain's, seaward positive, against the same 10 photo dates, so the 6|7 face can be scored on its own boundary. Observed: across 5|6 117, 116 and 155 m in 1978, 1993 and 1995; across 6|7 18, 44 and 33 m (`../../../1-observations/analysis/1-aerial_photo_gap/tables/photo_boundary_gaps.csv`).

RMSE of the gap (m), all 10 dates:

| variant | strength | across 5\|6 | across 6\|7 | 6\|7 gap in 1993 |
|---|---|---|---|---|
| no groin | 0 | 103.6 | **15.0** | 20 |
| face 5\|6 (runner) | 0.16 / 0.25 | 37.8 / **29.0** | 34.8 / 43.9 | −1 / −6 |
| face 6\|7 | 0.16 / 0.25 / 0.35 | 120.3 / 131.9 / 138.8 | 76.8 / 112.7 / 147.4 | 112 / 147 / 183 |
| field 5 to 7 | 0.16 / 0.25 | 44.1 / 22.1 | 62.8 / 95.3 | 99 / 133 |
| field, half kept | 0.16 / 0.20 / 0.25 | 41.3 / 31.6 / 25.5 | **23.0** / 31.1 / 41.4 | 55 / 66 / 77 |

1. **The 6|7 face fails on its own boundary too.** It moves GIS 6's sand into GIS 7, so the 6|7 gap reaches 112-183 m against 18-44 m observed, and the 5|6 gap turns negative (GIS 6 behind GIS 5) where the photos have it 116-155 m ahead. It is the worst variant on both boundaries at every strength.
2. **The runner's face misses 6|7 the other way:** with the sand kept in GIS 6, GIS 7 ends level with or behind it (−1 to −6 m in 1993 against +44 m).
3. **Only the field with half kept is close on both:** 31-41 m at 6|7 and 25-32 m at 5|6 at 0.20-0.25, though its 6|7 gap still runs high (66-77 m in 1993).
4. No groin scores best at 6|7 only because the observed 6|7 gap is small; it has no 5|6 gap at all.

## What this means

- **The runner keeps its GIS 5|6 face.** The hindcast windows, which the paper is scored on, slightly prefer it, and the main code is not changed.
- **For the forecast, 5|6 is also the right place.** Only the southernmost groin is rebuilt in 2026, on transect 37 at the GIS 5|6 boundary.
- **Recorded as a finding:** the build-up's shape fits much better when domain 6 is treated as the field, with half the trapped sand kept inside it and half against its north end (21.7 against 25.7 m). If a future run starts before the fillet exists, for example a forecast of a rebuilt three-groin field, this is the layout to use. `field_windows.patch_groin("field", 0.5)` is how to apply it.
- The windows were run only with half kept. No share would fix GIS 7 in 1996-2009, where the record erodes and every layout builds.
- The refactored `face_experiment.py` (which now uses `field_windows.patch_groin`) was not rerun; the scores above come from the runs made before the refactor.
