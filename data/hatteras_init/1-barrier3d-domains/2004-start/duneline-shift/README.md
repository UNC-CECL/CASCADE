# duneline-shift (2004-start)

The 2004 dune line measured against this product's interior row 0, by
`HAT_measure_duneline_shift.py --year 2004`.

**This is a CONTROL, not an input.** Nothing consumes it. It exists to bound the
definitional offset between a digitized dune line and the extractor's row 0, by
measuring a pair that is close in date: the 2004 line against a 2009 DEM.

Read it BY REACH, not island-wide:

| reach | 2004 residual |
|---|---|
| mid-island D41–70 | **−1.4 m** |
| the 25 smallest-\|shift\| domains | −0.5 m |
| Cape Point D1–8 | +47.8 m |
| Avon–Salvo D17–40 | +42.2 m |

On the stable mid-island the line and row 0 are the same feature to within half
a cell. Where the residual is large those are the nourished and hotspot reaches,
where five years of real change sits between the 2004 line and the 2009 surface —
so the control cannot bound the definitional term there.

The island-wide median is +14.6 m, **not** zero. Quoting that number as "the
feature offset" would be wrong; the near-zero result is a mid-island result.

Superseded as the primary method by the 1984/1997 difference in
`../../1984-start/duneline-shift/`, which cancels the definitional term outright
rather than bounding it.

## Re-measured 2026-09-03

Against the product's own raster (`hat_elevation_products`) rather than the clip tree, and stamped `topo_version`. The reach structure is unchanged — mid-island still −1.4 m — so every claim above stands; only the island-wide median moved, +14.2 → +14.6 m.
