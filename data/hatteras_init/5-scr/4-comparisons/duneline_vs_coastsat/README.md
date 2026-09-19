# Dune line vs CoastSat shoreline — methods, dates, results

Does the digitised dune line move with the satellite waterline? One folder
per hindcast window (`1984_2004/`, `1996_2010/`, `2004_2024/`, `2010_2024/`),
each holding a scatter and an alongshore figure with its tables, provenance
and captions under `supporting/`; and `alongshore_four_windows.png` at this
level, the four windows on one y axis. Written 2026-09-15 (Hannah's design,
by interview; the figures and numbers below are the state that evening).

**Since 2026-09-18 both sides are NET CHANGE between the same two dates**
(Hannah: the CoastSat-vs-dune comparison should use net position change on
both sides). Both are read from the stored products in `3-rates/`
(`duneline/endpoint` and `coastsat/endpoint`). The CoastSat LRR is no longer
drawn or scored here; it stays the model's scoring target in
`3-rates/coastsat/lrr/` and `rate_windows/coastsat/`.

```
python scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py --start-year 1984 --end-year 2004
python scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py --grid   # the stacked four
```

## The two net changes, per GIS domain, seaward positive

| side | what it is | stored product |
|---|---|---|
| **Dune line** | The end digitized line minus the start line. Each line is intersected with the 100 m transects by `duneline_to_raw_offsets.py`, one station per transect measured from the offshore datum; the ~5 transects of each 500 m domain are averaged. | `3-rates/duneline/endpoint/<window>/` |
| **CoastSat shoreline** | For each CoastSat transect, the mean chainage within **±6 months (182.6 days) of each dune-line survey date**, end minus start, so each end averages one full seasonal cycle; averaged over the ~10 transects of each domain. | `3-rates/coastsat/endpoint/<window>/` |

Both are drawn as the net change divided by the survey interval (m/yr), so
the four windows of different length share one axis. The metres are in each
window's `supporting/domain_comparison.csv`, along with the beach-width change
(shoreline minus dune).

A dune line and a waterline are different features. A gap between their
rates is beach-width change as much as disagreement, and the purpose here is
the proxy question (does the dune line track the shoreline?), not beach
width.

## Survey dates: which imagery each line is, and when it was flown

A period year reaches its line through `hat_topo_version.DUNE_LINE_FOR_YEAR`;
the line and its raw file are named for the **imagery vintage**, never the
period. The dates centre the endpoint windows and set the dune interval.

| period year | line | imagery | date | where the date comes from |
|---|---|---|---|---|
| 1984 | `duneline_1984.geojson` | USGS aerial photo (Henderson release) | **1984-09-19** | `D:\Hatteras_GIS\Aerial\1984_henderson\1984_metadata`, Calendar_Date |
| 1996 | `duneline_1997_v2.geojson` | USGS aerial photo (Henderson) — no 1996 line; the nearest island-wide survey | **1997-10-12** | `D:\Hatteras_GIS\Aerial\1997_henderson`, Calendar_Date. (A 1996 Henderson set exists, dated 1996-10-14, undigitised) |
| 2004 | `duneline_2004.geojson` | Google Earth historical capture | **2004-05-25** | Hannah, 2026-09-15; the raw_GE frames carry no date |
| 2010 | `duneline_2009.geojson` | Google Earth historical capture — no 2010 aerial imagery | **2009-05-30** | Hannah, 2026-09-15 |
| 2024 | `duneline_2023.geojson` | NOAA NGS 2023 emergency-response imagery | **not known**; assumed 2023-07-01 | `D:\Hatteras_GIS\Aerial\2023`; its metadata gives only the 2015–2023 series extent. Each affected PROVENANCE.md flags the assumption and reports a ±6 month sensitivity |

So the survey intervals are 19.68 yr (1984–2004), 11.63 yr (1996–2010, really
1997 to 2009), 19.10 yr (2004–2024, really to 2023) and 14.09 yr (2010–2024,
really 2009 to 2023). Both sides of the comparison span exactly these intervals.

Every raw dune file is built by the same shapely intersection since
2026-09-15 (the 1984 and 2004 ArcGIS exports were rebuilt that afternoon), so
a change between any two lines carries no method term. The dates live in the
script's `KNOWN_SURVEY_DATES`, keyed by vintage, and in
`2-brie-offset/dunelines/README.md`; the geojsons themselves carry none.

## What is inside each endpoint window

Median CoastSat observations per transect inside the one-year window, and
the span the observations actually cover. The **1984 window is one-sided**:
the CoastSat record begins 1984-09-21, two days after the photo, so that end
is an autumn-and-winter mean of about seven images, not a full seasonal cycle.

| window | start window | obs | end window | obs |
|---|---|---|---|---|
| 1984–2004 | 1984-09-21 to 1985-02-28 (truncated) | 7 | 2003-11-29 to 2004-11-07 | 20 |
| 1996–2010 | 1997-05-04 to 1998-03-04 | 6 | 2008-12-04 to 2009-11-05 | 9 |
| 2004–2024 | 2003-11-29 to 2004-11-07 | 20 | 2023-01-21 to 2023-12-15 | 39 |
| 2010–2024 | 2008-12-04 to 2009-11-05 | 9 | 2023-01-21 to 2023-12-15 | 39 |

No transect had an empty window in any period (906 transects checked each
time). A storm inside a window moves that end.

## Results, over 90 domains (seaward positive)

| window | dune line | CoastSat shoreline | beach width | r | slope (shoreline on dune) | RMSE (m/yr) |
|---|---|---|---|---|---|---|
| 1984–2004 | −15.4 m (−0.78 m/yr) | −19.4 m (−0.99) | −4.1 m | 0.50 | 0.49 | 1.98 |
| 1996–2010 | −16.2 m (−1.39) | −2.0 m (−0.17) | +14.2 m | 0.66 | 0.73 | 2.02 |
| 2004–2024 | −1.7 m (−0.09) | +7.8 m (+0.41) | +9.5 m | 0.37 | 0.34 | 1.89 |
| 2010–2024 | +1.4 m (+0.10) | +7.0 m (+0.50) | +5.6 m | 0.71 | 0.91 | 1.33 |

These are on the lines as re-digitized on 2026-09-18. The 1984–2004 row did
not change with the re-digitization, since neither of its lines moved.

The dune line follows the shoreline's alongshore pattern in every window.
2004–2024 is the weakest (r 0.37). In three windows of four, the dune line
lost ground relative to the waterline, so the beach widened. 1984–2004 is the
exception: there the shoreline retreated slightly more than the dune.

## Figures and conventions

Red is the dune line and dark blue the CoastSat shoreline, both as net change:
the house red/blue pair marks FEATURE here, not vintage or sign.
Village bands, the Buxton groin (solid hairline) and the Avon and Rodanthe
piers (dotted) are drawn as on every alongshore figure. Domain 1 is Cape
Point, 90 is Pea Island. Statistics are in each window's
`supporting/PROVENANCE.md`; captions in `supporting/CAPTIONS.md`.
