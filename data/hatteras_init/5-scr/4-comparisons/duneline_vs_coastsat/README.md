# Dune line vs CoastSat shoreline — methods, dates, results

Does the digitised dune line move with the satellite waterline? One folder
per hindcast window (`1984_2004/`, `1996_2010/`, `2004_2024/`, `2010_2024/`),
each holding a scatter and an alongshore figure with its tables, provenance
and captions under `supporting/`; and `alongshore_four_windows.png` at this
level, the four windows on one y axis. Written 2026-09-15 (Hannah's design,
by interview; the figures and numbers below are the state that evening).

```
python scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py --start-year 1984 --end-year 2004
python scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py --grid   # the stacked four
```

## The three rates, per GIS domain, seaward positive

| rate | what it is | source |
|---|---|---|
| **Dune line** | Two digitised lines differenced. Each line is intersected with the 100 m transects (`2-brie-offset/transects/`) by `duneline_to_raw_offsets.py`, one station per transect measured from the offshore datum line; the ~5 transects of each 500 m domain are averaged, exactly as the hindcast's end-year target loader reads the same files. End minus start, sign flipped so seaward is positive, divided by the interval between the two **survey dates**. | `2-brie-offset/raw_offsets/<vintage>_duneline_offset_raw.csv` |
| **CoastSat LRR** | The linear regression rate: for each CoastSat transect, an ordinary-least-squares slope of waterline chainage against date over the calendar window (1 Jan of the start year to 31 Dec of the end year; ~250 observations per transect). Averaged over the ~10 transects in each domain. This is the quantity the model is graded against. | `5-scr/3-rates/coastsat/lrr/<start>_<end>/transect_lrr_full.csv` |
| **CoastSat endpoint** | Built to be like-for-like with a two-survey dune line. For each transect, the mean chainage of every observation within **±6 months (182.6 days) of each survey date**, a one-year window centred on the survey so each end averages one full seasonal cycle. End mean minus start mean, divided by the survey interval. Averaged per domain the same way. | `5-scr/1-observations/coastsat_timeseries/`, computed per run into `supporting/transect_coastsat_endpoint.csv` |

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
really 2009 to 2023), while the LRR spans the calendar window in each case.

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
time). A storm inside a window moves that end; that is why the LRR is
reported beside the endpoint.

## Results, island means over 90 domains (m/yr, seaward positive)

| window | dune line | CoastSat LRR | CoastSat endpoint | r, dune vs LRR / vs endpoint | slope | endpoint vs LRR, r |
|---|---|---|---|---|---|---|
| 1984–2004 | −0.78 | −1.13 | −0.99 | 0.44 / 0.50 | 0.41 / 0.49 | 0.94 |
| 1996–2010 | −1.39 | −0.32 | −0.17 | 0.59 / 0.66 | 0.61 / 0.73 | 0.90 |
| 2004–2024 | −0.09 | +0.49 | +0.41 | 0.25 / 0.37 | 0.23 / 0.34 | 0.93 |
| 2010–2024 | +0.10 | +1.15 | +0.50 | 0.66 / 0.71 | 0.91 / 0.91 | 0.89 |

**The 1997, 2009 and 2023 lines were re-digitized on 2026-09-18**, and the
three windows that use them were rebuilt the same day. The 1984–2004 row is
unchanged (neither of its lines moved). With the old lines, dune against LRR
gave r = 0.18 (1996–2010), 0.06 (2004–2024) and 0.51 (2010–2024).

Read across: the two CoastSat estimators agree with each other in every
window (r 0.89–0.94), so none of the disagreement with the dune line is an
estimator effect. Since the re-digitization, the dune line tracks the
shoreline in every window except 2004–2024 (r 0.25). The rest of this
paragraph was written on the 09-15 lines, and its domain-by-domain reading
has NOT been re-checked against the new ones. In
every window the dune line moves landward relative to the waterline by
roughly 1 m/yr, strongest after 2004 when the shoreline advanced from Cape
Point to Avon (GIS 1–20, 28–37) while the dune kept retreating; and the dune
advanced at GIS 49–55 and 60–74 while the shoreline held. The single-domain
spikes (GIS 1 in 1984–2004 and 1996–2010, GIS 78 in 1996–2010) read as
digitising or feature questions on the individual lines, not coastal change.

## Figures and conventions

Red is the dune line, dark blue the CoastSat LRR, light blue dashed the
CoastSat endpoint: the house red/blue pair marking FEATURE here, not vintage.
Village bands, the Buxton groin (solid hairline) and the Avon and Rodanthe
piers (dotted) are drawn as on every alongshore figure. Domain 1 is Cape
Point, 90 is Pea Island. Statistics are in each window's
`supporting/PROVENANCE.md`; captions in `supporting/CAPTIONS.md`.
