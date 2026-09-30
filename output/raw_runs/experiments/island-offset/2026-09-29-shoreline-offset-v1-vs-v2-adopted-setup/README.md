# 2026-09-29 — shoreline offset v1 vs v2, beside the dune line, on the adopted setup

Hannah, 2026-09-29: re-run the shoreline arm on shoreline offset **v2**, the
CoastSat mean over ±1 yr of the start DEM's lidar flights, which became CURRENT
that day. Every earlier shoreline run read **v1**, the calendar means.

The option A study (`../2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/`)
ran before three adopted changes: the per-cell dune ceilings and the beach/dune cap
fix (later on 09-28), and the split12 storm files (09-29). A v2 run made now
would differ from its v1 runs in all of those at once. Hannah chose a clean
three-way study instead, with every arm on today's setup, so **v1 → v2 is the
only difference between the two shoreline arms**.

| | |
|---|---|
| arms | `duneline` (2-brie-offset/<year>/duneline/v1), `shoreline_v1` (…/shoreline/v1, pinned with `HAT_OFFSET_VERSION_<year>_SHORELINE=v1`), `shoreline_v2` (…/shoreline/v2, CURRENT) |
| shoreline windows | v1: calendar 1995–1997 and 2009–2011. v2: 1995-10-12 → 1997-10-12 and 2008-08-17 → 2010-08-17 |
| waves | option A: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5 |
| ends | zeroBE, as in option A (the edgeBE ends were solved on the dune-line offset and would favour it) |
| scope | natural and full management × 1996–2010 and 2010–2024; relocations and groins off; 12 runs, 4–7 min each |
| setup | code defaults: adopted dune ceilings, `v3_split12_trim24` storms, beach/dune cap on added sand only; Barrier3D `hatteras/adopted` d6546c7 |
| check | each run's metadata `island_offset_version` must equal the arm's (`duneline/v1`, `shoreline/v1`, `shoreline/v2`), or scoring stops |

## Answer

**v1 → v2 is not a lever.** On every score and in every run the two shoreline
versions agree to within 0.015 of raw variance explained and 1.5 points smoothed. The
modelled total change moves by 0.8–1.5 m on average (mean absolute, interior) and at most
5.7 m in any domain. That is the size of the ~2.5 m standard error on a transect
mean. The island-wide mean shift is 0.1 m. So the choice of averaging window,
calendar or DEM-centred, does not change the result. v2 stays CURRENT because it
is dated like the DEM, not because it scores better.

**Shoreline vs dune line holds on the adopted setup, and widens slightly in
1996–2010.** Both shoreline versions are ahead of the dune line by 2–3 points raw.
Under management the gap is 7–8 points smoothed, larger than the 0–1 point in option A.
2010–2024 fails in every arm, as before (the 2021 step). On the raw score the dune
line is less bad there, by 3–4 points managed and 12–13 natural. Smoothed, the two
are level or the shoreline is ahead by up to 3 points. So the offset source remains a
small lever and is not the cause of the 2010–2024 failure.

## Scores against CoastSat (each period's own LRR, LOWESS 7, interior GIS 2–89)

| arm | scenario | period | raw explained | smoothed | bias (m/yr) | raw r |
|---|---|---|---|---|---|---|
| dune line | natural | 1996–2010 | +20.7% | +25.5% | −0.18 | 0.52 |
| shoreline v1 | natural | 1996–2010 | +23.1% | +28.9% | −0.13 | 0.51 |
| shoreline v2 | natural | 1996–2010 | **+23.4%** | **+29.2%** | −0.13 | 0.52 |
| dune line | managed | 1996–2010 | +18.1% | +16.3% | −0.03 | 0.45 |
| shoreline v1 | managed | 1996–2010 | +20.8% | +23.2% | +0.00 | 0.47 |
| shoreline v2 | managed | 1996–2010 | +20.7% | **+24.7%** | −0.00 | 0.47 |
| dune line | natural | 2010–2024 | −238% | −224% | −2.31 | 0.39 |
| shoreline v1 | natural | 2010–2024 | −251% | −221% | −2.34 | 0.37 |
| shoreline v2 | natural | 2010–2024 | −249% | −221% | −2.33 | 0.37 |
| dune line | managed | 2010–2024 | −88% | −83% | −1.46 | 0.33 |
| shoreline v1 | managed | 2010–2024 | −92% | −80% | −1.47 | 0.30 |
| shoreline v2 | managed | 2010–2024 | −91% | −80% | −1.47 | 0.30 |

Against option A (09-28, before the dune ceilings, cap fix and split storms),
1996–2010: dune line natural +19.8 → +20.7 raw and managed +18.2 → +18.1; shoreline v1
natural +22.6 → +23.1 and managed +22.4 → +20.8.

## Each offset on its own feature (metres over 14 yr, LOWESS 7)

| arm | scenario | period | graded on | explained | r | bias (m) |
|---|---|---|---|---|---|---|
| dune line | managed | 1996–2010 | dune-line net change | −36% | 0.32 | +9.6 |
| dune line | managed | 2010–2024 | dune-line net change | −64% | 0.14 | −5.5 |
| shoreline v1 | managed | 1996–2010 | total change (own LRR × 14) | +21% | 0.47 | +0.0 |
| shoreline v2 | managed | 1996–2010 | total change (own LRR × 14) | +21% | 0.47 | −0.1 |
| shoreline v1 | managed | 2010–2024 | total change (own LRR × 14) | −92% | 0.30 | −20.6 |
| shoreline v2 | managed | 2010–2024 | total change (own LRR × 14) | −91% | 0.30 | −20.5 |
| shoreline v1 | managed | 1996–2010 | projected (1996–2024 LRR × 14) | −8% | 0.43 | −7.4 |
| shoreline v2 | managed | 1996–2010 | projected (1996–2024 LRR × 14) | −9% | 0.42 | −7.4 |
| shoreline v1 | managed | 2010–2024 | projected (1996–2024 LRR × 14) | −49% | 0.24 | −6.9 |
| shoreline v2 | managed | 2010–2024 | projected (1996–2024 LRR × 14) | −49% | 0.24 | −6.9 |

Natural-run rows are in `tables/all_runs.csv`. The two features are different
targets on different estimators: each ranks an offset against its own feature,
not against the other.

## v2 minus v1, modelled total change (interior GIS 2–89)

| scenario | period | mean (m) | mean \|diff\| (m) | max \|diff\| (m) | at GIS |
|---|---|---|---|---|---|
| natural | 1996–2010 | −0.10 | 0.92 | 3.75 | 5 |
| natural | 2010–2024 | +0.12 | 1.54 | 5.57 | 52 |
| managed | 1996–2010 | −0.07 | 0.79 | 3.73 | 5 |
| managed | 2010–2024 | +0.07 | 1.35 | 5.70 | 79 |

## Figures

```
figures/total_change_difference_shoreline_v2_minus_v1.png        the answer: v2 − v1, 4 panels, same axis as the difference figures
figures/*_full_management_shoreline_v1.png                        the four house figures, shoreline v1 against the dune line
figures/*_full_management_shoreline_v2.png                        the same with shoreline v2
```

- The `duneline_offset_vs_duneline_change_*` figure appears once per suffix. The
  two copies are identical, because the dune-line arm does not depend on the shoreline version.
- **Axis note.** `total_change_difference_shoreline_minus_duneline_full_management_shoreline_v2`
  reaches −21.2 m at GIS 80 (1996–2010), past the shared −20 m. `house_figures`
  widened it to (−25, 30) with a warning, so that one figure is not on the shared
  axis. Its v1 twin stays inside (−20, 30).

## Reproduce

```
python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py run --jobs 4   # runs, then scores
python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py score
python scripts/hatteras_ms/experiments/HAT_offset_source_shoreline_v2.py plot
```

`logs/<period>/<arm>_<scenario>.log`, `logs/driver.log`; runs under
`runs/<arm>_<scenario>/<period>/zeroBE/` (on disk only, no `.npz`).
