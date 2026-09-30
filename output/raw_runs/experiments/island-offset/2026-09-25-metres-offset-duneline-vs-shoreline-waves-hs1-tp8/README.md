# 2026-09-25 — dune line vs shoreline as the island offset, 1996–2010

Hannah, 2026-09-25: how much does setting the island's planform orientation
from the dune line, rather than the CoastSat shoreline, change the output?

| | |
|---|---|
| offsets | `2-brie-offset/1996/duneline/v1` and `1996/shoreline/v1`: both metres, the wrap-around written in the file (`HAT_ISLAND_OFFSET_SOURCE`) |
| waves | Hs 1.0 m, Tp 8 s, asymmetry 0.8, the best managed 1996–2010 setting found. **Tuned with the dune-line offset**, so the dune line has a head start at one setting; the shoreline offset gets a fair chance through the high-angle fraction, 0.3 / 0.4 / 0.45 / 0.5 / 0.55, the lever that mattered most |
| scope | natural and full management, 1996–2010; 20 runs, both offsets run fresh here on the same code and Barrier3D (route_overwash fix) |
| score | as `wave-climate/2026-09-25-wave-grid-smoothed-score/`: share of the alongshore variation explained by the model smoothed like the CoastSat target (LOWESS 10 domains, southern 10 raw), interior GIS 2–89; raw score, bias and correlation beside it |

## Layout

```
tables/all_runs.csv            every run: offset, scenario, high-angle, scores, status, run folder
figures/duneline_offset_vs_duneline_change_full_management.png
figures/shoreline_offset_vs_coastsat_total_change_full_management.png
figures/shoreline_offset_vs_coastsat_projected_change_full_management.png
figures/total_change_difference_shoreline_minus_duneline_full_management.png
                               (a) full management 1996–2010, (b) 2010–2024, high-angle 0.45;
                               redrawn 2026-09-28 in the option A study's form
                               (metres, observations LOWESS 7, no scores on the
                               figures); the two original figures (profiles_,
                               scores_vs_high_angle_) are in git history
logs/<offset>_<scenario>/<settings>.log
runs/<offset>_<scenario>/1996_2010/zeroBE/<run_name>/     on disk only
```


**2010–2024 added 2026-09-28** (full-management pair at the same settings, `run-2010`), so the figures show full management × both periods like every island-offset study: (a) 1996–2010, (b) 2010–2024. The natural runs stay in `tables/all_runs.csv`, not in the figures.

## Reproduce

```
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py run
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py run-2010
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py plot
```

## Results

Run 2026-09-25, 14:01–14:15; all 20 scored. Smoothed share of the alongshore
variation explained (raw in brackets), bias in m/yr:

| high-angle | dune line, natural | shoreline, natural | dune line, managed | shoreline, managed |
|---|---|---|---|---|
| 0.30 | −20% (−92%), −0.87 | −15% (−53%), −0.80 | +3% (−27%), −0.39 | +1% (−20%), −0.34 |
| 0.40 | +16% (−3%), −0.57 | **+23%** (+11%), −0.53 | +20% (+11%), −0.18 | **+24%** (+15%), −0.15 |
| 0.45 | +14% (+5%), −0.45 | +22% (+15%), −0.40 | **+23%** (+17%), −0.04 | +22% (+16%), −0.04 |
| 0.50 | +8% (+1%), −0.39 | +15% (+11%), −0.31 | +16% (+10%), +0.02 | +12% (+10%), +0.02 |
| 0.55 | +11% (0%), −0.36 | +19% (+10%), −0.31 | +11% (+8%), +0.04 | +5% (+7%), +0.05 |

- **Under full management the source barely matters**: at the best setting
  each scores +23–24%, bias the same; the shoreline offset peaks at
  high-angle 0.4 rather than 0.45.
- **In the natural scenario the shoreline offset scores higher at every
  high-angle** (+22% against +14% at 0.45; correlation 0.60–0.66 against
  0.44–0.59): management's dune rebuilding and fills damp the orientation's
  effect, the natural run shows it.
- The rates differ by up to 1–1.8 m/yr domain by domain, mostly at
  Tri-Village (GIS 75–90) and, natural only, GIS 15–20; elsewhere under
  0.5 m/yr.
- Neither offset produces the observed accreting peaks at GIS 18, 29 and 42:
  the orientation is not what the model is missing there.
- The waves were tuned with the dune-line offset, and the shoreline offset
  still matches or beats it.
