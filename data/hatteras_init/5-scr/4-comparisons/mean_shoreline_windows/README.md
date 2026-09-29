# mean_shoreline_windows — one period's mean shoreline, over two windows

The CoastSat mean shoreline each period starts from, averaged over the
**3-yr calendar window** (what the shoreline island offset v1 was built from)
and the **2-yr window centred on the start DEM's lidar flights** (what v2 is
built from), differenced transect by transect. Asked for by Hannah on
2026-09-29, after the offset moved onto the DEM-centred window.
Written by `scripts/input_prep/5-scr/4-comparisons/mean_shoreline_windows/coastsat_mean_shoreline_windows.py`.

| period | 3-yr calendar | 2-yr DEM-centred | centred on |
|---|---|---|---|
| `1996/` | 1995–1997 | 1995-10-12 – 1997-10-12 | 1996 ALACE lidar, flown 1996-10-09 – 10-16 |
| `2010/` | 2009–2011 | 2008-08-17 – 2010-08-17 | 2009 USACE NCMP lidar, flown 2009-08-10 – 08-24 |

```
<period>/mean_shoreline_windows_<period>.png   (a) 2-yr minus 3-yr per transect
                                               and domain, + seaward;
                                               (b) standard error of each mean
<period>/PROVENANCE.md                         the numbers, the domains that move most
<period>/supporting/                           PDF, CAPTIONS.md, domain and
                                               transect tables
```

Both lines are read as stored in `1-observations/mean_shoreline/<window>/`;
nothing is re-averaged here. The windows share satellite passes, so no
significance is attached to a difference. The island offset is zeroed on its
own minimum, so only the alongshore-varying part of a difference reaches the
model; each PROVENANCE reports it beside the island mean.
