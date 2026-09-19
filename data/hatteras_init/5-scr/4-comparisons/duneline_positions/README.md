# duneline_positions — where the dune line sat in 1997, 2009 and 2023

Positions, not change. These are the dune lines that stand for the model years
1996, 2010 and 2024 (`hat_topo_version.DUNE_LINE_FOR_YEAR`), shown as maps and
as alongshore profiles of where they sit relative to NC-12 and the shoreline.
Written by `scripts/input_prep/5-scr/duneline_positions/duneline_positions.py`
(2026-09-18, Hannah's design by interview). Every figure is also published to
`output/figures/shoreline/duneline_positions/`.

```
overview/duneline_positions_overview.png   the island in three north-up
                                           segments, one common scale,
                                           villages bracketed, locator inset
zooms/zoom_<site>.png                      one site per figure, over imagery
zooms/duneline_positions_zooms.png         all five sites on one sheet
context/dune_to_nc12.png                   dune line to NC-12, per domain
context/beach_width.png                    dune line to the CoastSat shoreline
supporting/                                PDFs and CAPTIONS.md under each
                                           folder; the tables below
```

**Year colours:** one ordered ramp, light grey 1997, slate 2009 and black
2023, so the order reads at a glance and red and blue stay free for
seaward and landward. On the imagery each line is edged for contrast.

## The zooms

Each zoom is a window **three model domains (1.5 km) alongshore by 950 m
across**: 650 m landward and 300 m seaward of the 2023 line. Every panel has
the same extent and scale, over the **2023 NOAA NGS orthomosaic**
(`D:\Hatteras_GIS\Aerial\2023\2023_full_aerial.tif`, read through its
overviews), so the 2023 line can be checked against the dune it traces. The
1997 and 2009 lines were not traced from this image. NC-12 is drawn white,
with the year shown by dash pattern: dotted for the 1978 alignment (1997),
short dashes for 2008 (2009), long dashes for today's (2023).

| site | window | why centred there |
|---|---|---|
| Buxton | GIS 6–8 | the domain in GIS 1–15 with the largest net dune change 1997–2023 |
| Avon | GIS 23–25 | the same rule within GIS 21–31 |
| Tri-Village | GIS 77–79 | the same rule within GIS 68–83 |
| Mirlo Beach S-curves | GIS 83–85 | the same rule within GIS 84–90 |
| Largest change elsewhere | GIS 32–34 | the domain outside every named site with the largest change |

The choices are recorded in `supporting/zoom_sites.csv`. They follow the data,
so a re-digitized line can move a window.

## The two context profiles

- **Dune line to NC-12:** along each 100 m transect, the road's distance from
  the fixed offshore datum minus the dune line's; per domain the mean.
  Positive means the road lies landward of the dune.
  - Road lines: the 1978 export for 1997 and the 2008 export for 2009 (the
    alignments the model uses), and today's NCDOT alignment for 2023
    (`4-mgmt-forcing/road_offset/raw_offset/current/`, extracted 2026-09-18).
  - Where a transect crosses the road twice, the crossing nearest the dune is
    used. At Buxton (GIS 8–9) the transects also cross the leg that turns west
    toward Frisco.
  - GIS 1–7 has no value: the transects there do not reach NC-12.
  - The jump at GIS 82–88 in 2023 is the **2022 Jug Handle bridge**, which
    moved NC-12 500–800 m into Pamlico Sound. It is real.
  - Island means: 267 m (1997), 257 m (2009), 307 m (2023). The 2023 rise is
    that bridge.
- **Beach width, dune line to CoastSat shoreline:** along each CoastSat
  transect, the shoreline position minus where the dune line crosses it.
  - The shoreline is the mean satellite position within ±6 months of the
    dune-line image date, from `3-rates/coastsat/endpoint`.
  - Transects are extended 400 m landward first, because many start seaward
    of the dune. Without that, 256–470 of the 906 missed the dune line in the
    first draw.
  - Island means: 34 m (1997), 48 m (2009), 54 m (2023). The steps (+14 m,
    +6 m) match the beach-width change in `../net_change_1996_2024/`.
  - **Caveat:** 12 transects give a negative width in 1997, all at GIS 2–5
    near Cape Point, where the 1997 line lies seaward of the satellite
    shoreline. That is a question about how the 1997 line was traced there,
    and the domain means at GIS 2–5 (1–5 m) should be read with it.

The 2023 image date is not known and is assumed to be 1 July. Beach width for
2023 inherits that assumption; the dune positions do not.

## Tables (`supporting/`)

| file | one row per |
|---|---|
| `dune_to_nc12_transects.csv` | transect × year: dune and road stations, crossings, distance |
| `dune_to_nc12_domains.csv` | domain × year: mean, min, max, n |
| `beach_width_transects.csv` | CoastSat transect × year: dune and shoreline chainage, width |
| `beach_width_domains.csv` | domain × year: mean, min, max, n |
| `zoom_sites.csv` | zoom site: window and the rule that chose it |

## Rebuilding

```
python scripts/input_prep/5-scr/duneline_positions/duneline_positions.py
python scripts/input_prep/5-scr/duneline_positions/duneline_positions.py --no-imagery
```

Run the script with the project's `.venv` Python, which has `rasterio`. The
imagery needs the D: drive; `--no-imagery` draws the zooms on a plain
background.
