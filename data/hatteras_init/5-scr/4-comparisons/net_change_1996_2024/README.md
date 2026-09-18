# net_change_1996_2024 — net shoreline change against net dune-line change

Does the dune line move with the shoreline over 1996–2024, and over each of
its halves? Both sides are **net position change in metres, seaward
positive**, measured between the same dates. The gap between them is
beach-width change. Written by
`scripts/input_prep/5-scr/net_change/net_change_1996_2024.py` (2026-09-18,
Hannah's design by interview).

```
net_change_shoreline_vs_dune.png    (a) 1997–2023, (b) 1997–2009, (c) 2009–2023
supporting/
    domain_comparison.csv           window × domain: shoreline, dune, beach-width
                                    change (shoreline − dune), same sign or not
    island_summary.csv              per window: means, r, slope, RMSE, sign counts
    net_change_shoreline_vs_dune.pdf, CAPTIONS.md
```

The figure is also published to `output/figures/shoreline/`.

## The two sides

Both are read from stored products, not computed here.

| side | product | what it is |
|---|---|---|
| shoreline | `3-rates/coastsat/endpoint/<window>/` | mean CoastSat position within ±6 months of each dune-line image date, end minus start; ~10 transects per domain |
| dune line | `3-rates/duneline/endpoint/<window>/` | the end line minus the start line along the 100 m transects; ~5 per domain |

The dates are 1997-10-12, 2009-05-30 and 2023-07-01. The 2023 flight date is
not known and is **assumed**. The lines stand in for the model years 1996,
2010 and 2024. Both sides share the 2009 date, so on each side the two
halves add up to the whole exactly; the script checks this. The per-transect
CoastSat change, with the number of positions in each end window, is in
`3-rates/coastsat/endpoint/<window>/transect_endpoint.csv`.

## Results at build time (90 domains)

| window | shoreline | dune line | beach width | r | slope (dune on shoreline) | same sign |
|---|---|---|---|---|---|---|
| 1997–2023 | +5.0 m | −14.8 m | +19.8 m | 0.81 | 0.69 | 67 / 90 |
| 1997–2009 | −2.0 m | −16.2 m | +14.2 m | 0.66 | 0.60 | 58 / 90 |
| 2009–2023 | +7.0 m | +1.4 m | +5.6 m | 0.71 | 0.56 | 71 / 90 |

The dune line tracks the shoreline's PATTERN well (r 0.66–0.81), but it moved
landward relative to it in nearly every domain. Over 1997–2023 the beach
widened by about 20 m on average, and 14 m of that came in 1997–2009, when
the shoreline barely moved and the dune line retreated 16 m.

**An island-wide offset like that is worth checking before it is read as
physics.** The 1997 line is traced from a USGS aerial photo and the 2009 line
from a Google Earth capture. A difference in how the dune edge reads on the
two kinds of imagery would appear exactly this way, as a near-uniform shift.
Not yet tested.
