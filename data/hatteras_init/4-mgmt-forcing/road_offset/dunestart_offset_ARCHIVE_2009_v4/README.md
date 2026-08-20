# ARCHIVE — dune-start road offset on the **2009_v4** topography

**Frozen 2026-08-20. Do not regenerate into this folder, and do not feed it to a
model run.** The live outputs are in `../dunestart_offset/`, which is built on
`2009_v5`.

This is the complete `dunestart_offset/` tree exactly as it stood before the v5
extraction — all 18 files, including the figures and the two audit markdowns.

## Why it exists

`2009_v5` is a **different DEM**, not just re-picked windows. It reads the
gap-filled elevation set (`2009_pea_hatteras_filled`): the 2009 base with its
LiDAR holes filled from the 2014 NOAA Post-Sandy DEM, mostly landward of NC-12.
The road setback is measured from interior row 0, so it moves with any change to
the interior — and this change is far larger than v3 → v4 was.

**89 of the 90 domains have a different interior SHAPE in v5**, and they are
uniformly deeper, because land that used to be no-data is now surveyed:

| | v4 rows | v5 rows |
|---|---|---|
| GIS 11 | 62 | 157 |
| GIS 52 | 32 | 59 |
| GIS 79 | 67 | 101 |
| GIS 85 | 50 | 176 |

Only GIS 1 kept its shape (values still changed).

## What that fixed

The v4 audit in this folder diagnosed its own drowning and named the cause:

> these domains **do not drown on measured water — they drown on LiDAR coverage
> gaps.** Across the six flanking rows that fail at GIS 78/79/80 there are 106
> wet cells, of which **105 were never surveyed and 1 is genuinely measured wet.**

v5 is that diagnosis acted on, and the numbers moved the way it predicted:

| | this archive (v4) | live (v5) |
|---|---|---|
| Drowns at init, 1984 | 3 | **0** |
| Drowns at init, 2004 | 3 | **0** |
| Negative setbacks, 1984 | 6 — GIS 11, 12, 13, 84, 85, 86 | 6 — GIS **10**, 11, 12, 84, 85, 86 |
| Negative setbacks, 2004 | 0 | 0 |
| DEM | `2009_pea_hatteras` | `2009_pea_hatteras_filled` |
| Picks | `HAT_dune_search_windows_2009_v4.json` | `..._2009_v5.json` |

Because nothing drowns on v5, **no setback is moved seaward any more** — the
`MOVED_SEAWARD` correction that v4 had to apply at GIS 78–80 is simply not
needed. The v4 audit's own caveat ("any managed-vs-unmanaged result at GIS 78–80
needs that sentence attached to it") no longer applies to the live set.

## What changed in the setbacks themselves

Less than the interior change suggests, because most of the fill is well
landward of where the road sits:

| | domains changed | median Δ | range |
|---|---|---|---|
| 1984 | 21 of 82 | +0 m | −30 to +100 m |
| 2004 | 26 of 82 | +0 m | −30 to +100 m |

Largest movers, both years: **GIS 79 (400 → 500 m)** and **GIS 80 (440 → 490 m)** —
two of the three roadways the relocation logic acts on, and the same pair the
v3 → v4 mixed-version bug hit. Then GIS 23 (220 → 190), GIS 30/31/32 (±20 m),
GIS 10 (1984 20 → 0, 2004 70 → 40).

## What did NOT change, and why

The **road masks were not re-burned**, and did not need to be. They are
registered to the domain `resampled_*.tif` grids, which the gap-fill did not
touch: filled and unfilled arrays have identical shapes in all 90 domains, the
extractor ran with `REQUIRE_ROAD_MASKS = True` and did not error, and
`HAT_check_geojson_vs_mask.py` re-confirms **100 % of centreline vertices inside
the mask footprint, 0 outside**.

The **8 relocation DROWNS in `RoadSetback_audit.csv` are unchanged** (GIS 11–15
and 84, 86, 87), and the FITS/DROWNS split is 167/8 in both versions. That is
not a stale read — 24 of the 28 columns differ between the two files. Those rows
are `kind = relocation`, hypothetical prescribed-relocation scenarios, not the
initial setbacks; the initial setbacks now drown nowhere. Whether those
scenarios are reachable is a separate question this archive does not answer.

## Identifying marks

If you ever need to tell a v4 file from a v5 one:

| | this archive (v4) | live (v5) |
|---|---|---|
| `RoadOffset_dunestart_audit.md` | says `2009_v4`, picks `..._2009_v4.json` | says `2009_v5`, picks `..._2009_v5.json` |
| Setback at GIS 79 | 400 m | 500 m |
| Drowns at init | 3 per year | 0 |
| DEM in the audit header | `2009_pea_hatteras` | `2009_pea_hatteras_filled` |

See `../dunestart_offset_ARCHIVE_2009_v3/README.md` for the v3 → v4 step.
