# duneline-shift

**Moved here from `1984-start/duneline-shift/` on 2026-09-03**, so the
measurement of N sits with the analysis built on it. The path is resolved by
`hat_topo_version.duneline_shift_dir(product)` — ONE definition, because eight
scripts used to build it by hand and two of them WRITE here. 2004-start keeps
the plain `2004-start/duneline-shift/` layout; that asymmetry is deliberate and
is documented in the resolver.

**This folder is an input, not a by-product.** `duneline_retreat_1984_1997.csv`
is where N comes from — `HAT_insert_seaward_rows.py` builds v4 and v5 from it,
and `HAT_report_row_insert_scope.py` reports against it. Deleting it means no
future insert version can be built.

## Removed 2026-09-03

| file | why |
|---|---|
| `superseded/dsas_shift_1984.csv` | the `--shift-source dsas` option was removed entirely. It measured the SHORELINE, not the dune line, understating dune retreat ~2.3x (GIS 85: 19.5 m against a measured 58.2 m) — and its path was already broken, pointing one level above where the file sat. `HAT_measure_dsas_shift.py` is now marked ORPHANED. |
| `superseded/duneline_retreat_1984_2004.csv` | 1984−2004 retreat scaled 20→12 yr; superseded by the unscaled 1984−1997 pair, which needs no scaling. No script referenced it. |

Sizes and reasons: `../../../archive_purge_20260903.csv`.

---

How far seaward of the model's interior row 0 the digitized dune line sits, per
domain, per profile. Written by
`scripts/input_prep/1-barrier3d-domains/HAT_measure_duneline_shift.py`.

**Sign: positive = the dune line is SEAWARD of interior row 0**, i.e. the number
of metres the island has to move seaward for row 0 to land on the digitized
line.

## Why it was measured

The 1984-start topography is a 1996 ALACE beach and foredune on a 2009 backdune,
and the NC-12 lines it is measured against are 1984-vintage. Where the island
migrated far enough between 1984 and the surveyed surface, the 1984 roadbed ends
up **seaward of the dune crest** and the setback goes negative. GIS 85 is the
only domain where it does, at −10 m, which the model CSV floors to 0.

Zero is not inert. `roadway_manager.road_relocation_checks` does
`road_setback += dune_migrated` every year, so the first year with any erosion
drives it negative and relocates the road. In the shipped
`HAT_1984_2004_calibBE_road_bdm_groin` run, the relocation year is monotone in
the input setback and not in the physics:

| GIS | model setback | first relocation | n relocations |
|---|---|---|---|
| 85 | 0 m | **year 1** | 3 |
| 86 | 0 m | **year 1** | 2 |
| 84 | 25 m | year 3 | 2 |
| 11 | 10 m | year 9 | 1 |
| 10 | 20 m | year 14 | 1 |

## Method

A domain clip is north-up, so one extractor profile is a raster row of constant
y. The dune-line geometry is intersected with that horizontal line and the
crossing's easting converted to a cross-shore cell index, then differenced
against interior row 0 for the same profile.

The frame is the extractor's own — its `c0` and per-profile `shear`, through the
inversion `cell_to_map` documents — so this is measured in the frame CASCADE
indexes. That is the difference between this and the legacy
`HAT_setback_from_lines.py` pass, which measured raw and ocean-first,
unstraightened, and is why the two disagree by a median 38 m.

A meandering line can cross one profile more than once; the crossing nearest row
0 is taken, not the first, which on the wide domains would pick a sound-side
meander.

## The control, and what it does and does not license

Run the same measurement on the 2004 pair — 2004 dune line against 2004-start,
whose DEM is 2009 — and read it **by reach**, not island-wide:

| reach | 1984 shift | 2004 control |
|---|---|---|
| mid-island D41–70 | +18.5 m | **−1.4 m** |
| the 25 domains with the smallest \|2004 shift\| | **+15.6 m** | −0.9 m |
| Cape Point D1–8 | +19.5 m | +47.8 m |
| Buxton/Avon S D9–16 | +17.9 m | +30.7 m |
| Avon–Salvo D17–40 | +12.7 m | +42.5 m |
| Rodanthe D71–83 | +47.1 m | +17.0 m |
| Pea Island D84–90 | +26.2 m | +16.9 m |

On the domains where the control says a digitized dune line and interior row 0
are the **same feature to within half a cell**, the 1984 pair still reads
+15.6 m. So the feature term is ~0 and the 1984 number is a date difference.

**What it does not license.** The island-wide 2004 median is +14.2 m, not zero.
The near-zero result is a mid-island result. Where the control is large — Cape
Point, Buxton/Avon, Avon–Salvo — that is five years of nourished and hotspot
change between the 2004 line and the 2009 surface, which is real, but it means
the control cannot bound the feature term on those reaches. `corr(1984, 2004)`
is −0.15 across the island; the two are not measuring the same thing, which is
the point.

**The 1984 side is mixed-vintage.** Row 0 is a 1996 feature only where ALACE
reached the dune; elsewhere it is 2009. So a per-domain 1984 shift spans 12 to
25 years of migration depending on the domain, and the island-wide median
(+18.9 m) is not a rate.

## An independent order-of-magnitude check

Driven by measured shoreline change, the model migrates GIS 85's dune **−70 m
over the 20-year hindcast** (3.5 m/yr). The measured 1984 shift there is
+65.9 m — the same size as two decades of its own modelled migration, which is
what a 12–25 year gap at that domain should look like.

## Files

| file | what |
|---|---|
| `duneline_shift_1984.csv` | per domain: median/p10/p90 shift, row 0 cell, dune-line cell |
| `duneline_shift_1984_profiles.csv` | per profile, for the alongshore gradient |

The 2004 control lives in `2004-start/duneline-shift/`.

## What consumes this

`HAT_insert_seaward_rows.py` reads `shift_m_median` and will not run without it —
N is a measurement, not a target. See
`1984-start/dune-topo/v1_{pad,translate,none}_measured/RUN_MANIFEST.txt`.

## Caveat carried forward

At GIS 85 the per-profile setback runs from −40 m at profile 0 to +30 m at
profile 43, near-monotone: the road is not parallel to the dune there. The
dune-line shift itself is much tighter (p10 +62, p90 +70), so the gradient is in
the road, not the dune. A single scalar shift is still fitting a median to a
ramp, and the alongshore half that was already positive gets pushed further
positive by it.
