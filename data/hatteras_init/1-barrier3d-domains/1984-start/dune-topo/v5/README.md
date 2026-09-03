# 1984-start dune-topo `v5`

`v3` (the 2026-09-02 re-pick) **plus** fabricated interior rows at the seaward
end, at every domain whose measured 1984→1997 dune-line shift rounds to at least
one cell. Built by:

```
python scripts/input_prep/1-barrier3d-domains/HAT_insert_seaward_rows.py \
    --variant pad --fill measured --shift-source date --n-rule measured \
    --domains measured --dst-version v5
```

## What problem this solves

The 1984-start topography is a 1996 ALACE beach and foredune grafted onto a 2009
backdune, but the NC-12 road lines are 1984-vintage. Between 1984 and 1996 the
island migrated landward, so at the worst domains the 1984 roadbed sits *seaward*
of the 1996 dune crest and the measured setback comes out negative:

| GIS | setback in v3 | in v5 |
|---|---|---|
| 85 | **−15 m** → floored to 0 | **+45 m** |
| 86 | **−10 m** → floored to 0 | **+20 m** |

A floored setback is not inert. `roadway_manager.road_relocation_checks` does
`road_setback += dune_migrated` every year, so the first year with any erosion
drives it negative and relocates the road. GIS 85 relocated in **model year 1**.
That is an artefact of the input, not a model result.

## What it does — and only what it does

`v5` is **purely additive**. Verified against `v3`, all 90 domains:

- every existing interior cell is bit-identical, `v5[n:] == v3`
- every dune array is copied unchanged
- no crest shaving, no landward retirement, no clamping

The island simply gets `N` cells wider at the seaward end, which moves interior
row 0 seaward and turns the setback from `(road − row0)` into
`(road − row0) + N·10`.

## How N is determined

N is a **measurement**, not a target: the per-domain median cross-shore offset
between the digitized 1984 and 1997 dune lines, measured in the extractor's own
frame by `HAT_measure_duneline_shift.py`, then `round(shift / 10)`.

The naive number — 1984 line to interior row 0 — is not usable, because it also
contains the offset between a digitized line and row 0. That is separated by
measuring the 1997 line the same way:

```
total (row0 − 1984)  =  feature (row0 − 1997)  +  date (1997 − 1984)
```

Island-wide medians: total **+23.0 m** = feature **+18.3 m** + date **+0.8 m**.
Only the **date** term is used. About 85% of the naive number was definitional.

The control is in the 2004 pair: measured identically, the 2004 dune line and
2004-start row 0 agree to +1.4 / −4.1 / −3.6 m on the stable mid-island
(GIS 40/50/60), so "digitized line vs DEM row 0" carries no material feature
offset.

## Scope — the whole island, by measurement

`v4` inserted rows only in the two historical relocation blocks (GIS 9–14,
84–87). That made *unchanged* mean two different things along the island: a
domain the measurement said needed nothing, and a domain that was never asked.
**30 domains had N ≥ 1 and were passed over by a scope decision**, several of
them a larger measured retreat than GIS 86's.

`v5` scopes by the measurement instead. **38 domains** get rows; the other 52 are
unchanged because N rounds to 0. Five of the 38 (GIS 2–6) carry no NC-12 at all —
they are included because the 1984 land missing from the 1996 survey is missing
whether or not a road ran over it.

N is **identical to v4 at all ten block domains**, so the two differ only in
scope.

## What is in the inserted cells

`--fill measured`: wherever the DEM has real dry land at the cell that a
prepended row lands on, that value is kept; only cells that were water in 1996
are invented, laid flat at the backdune platform (the median of interior rows
1–3). The real value is used as a **floor**, not taken raw — at an eroding domain
the 1996 surface is the *later and lower* one, so it bounds 1984's elevation from
below rather than estimating it. Taking it raw put a 0.57 m beach cell two rows
inside the 1984 interior at GIS 85.

Across the 38 domains, **3732 of 4900 inserted cells (76%) are measured DEM**;
the rest are fabricated platform.

## Known limitations — read these before using a result

**The interval is 13 years; the surface is 12.** The dune lines are 1984 and
**1997**, but the DEM surface at row 0 and at every inserted cell is 1996 ALACE
(verified: 100% of cells, all block domains, no 2009 or 2014). The source aerials
are 1996-10-14 and 1997-10-12 — 363 days apart — so N overstates the 1984→1996
retreat by close to one full year of it. This is **recorded, not corrected**:
scaling by 12/13 would assume retreat was steady across 13 storm-driven years,
which the record does not support. The effect is under one cell everywhere and
moves N at three domains: GIS 12 (2→1), GIS 84 (4→3), GIS 85 (6→5). Removing the
assumption entirely would mean digitizing a 1996 dune line from the 33
georeferenced 1996 frames already on disk.

**Progradation is not removal.** N is floored at 0. Where the 1997 line is
*seaward* of the 1984 one the island prograded, the 1996 survey already covers
everything the 1984 island had, and there is no missing land to restore — so
nothing is inserted and nothing is taken away. This version only ever *adds*
interior that existed in 1984 and is absent from the 1996 DEM.

**The rows are fabricated.** No survey covers land that had eroded away by 1996.
Any result that turns on island width, mean interior height, overwash flux or the
roadway relocation-room test at these domains is describing that fabrication as
much as it is describing Hatteras.

**The dune arrays are copied unchanged**, which assumes the 1984 crest height
equalled the 1996 one. Only its *position* in the model frame moves.

**The setback CSV in this folder is matched to this topography.** Reading the v1
or v3 setbacks against these arrays, or the reverse, restores exactly the
off-by-N error this version exists to remove.

## Files

| file | |
|---|---|
| `topography/`, `dunes/` | the arrays |
| `RoadSetback_1984_dunestart.csv` | setbacks matched to these arrays |
| `HAT_seaward_row_insert_audit.csv` | per-domain N, setback before/after, land rows, fabricated vs measured cells |
| `RUN_MANIFEST.txt` | build parameters and the caveats above |

## Per-domain record

| GIS | rows added | measured shift (m) | raw setback before → after (m) | % of block from DEM |
|---|---|---|---|---|
| 2 | 1 | 12.7 | - | 100% |
| 3 | 4 | 35.9 | - | 42% |
| 4 | 2 | 24.6 | - | 75% |
| 5 | 6 | 63.6 | - | 57% |
| 6 | 1 | 9.2 | - | 84% |
| 10 | 1 | 5.8 | +20 to +30 | 100% |
| 11 | 2 | 18.6 | +10 to +30 | 100% |
| 12 | 2 | 16.0 | +20 to +40 | 100% |
| 13 | 1 | 12.5 | +20 to +30 | 100% |
| 18 | 1 | 5.6 | +120 to +130 | 100% |
| 44 | 1 | 8.6 | +160 to +170 | 100% |
| 45 | 2 | 17.1 | +150 to +170 | 95% |
| 46 | 1 | 13.8 | +150 to +160 | 98% |
| 47 | 1 | 13.0 | +150 to +160 | 100% |
| 48 | 1 | 13.8 | +150 to +160 | 94% |
| 49 | 3 | 30.1 | +90 to +120 | 93% |
| 50 | 4 | 38.3 | +90 to +130 | 88% |
| 51 | 3 | 27.1 | +155 to +185 | 81% |
| 52 | 1 | 14.3 | +170 to +180 | 100% |
| 53 | 1 | 10.7 | +170 to +180 | 100% |
| 54 | 3 | 26.0 | +160 to +190 | 90% |
| 55 | 1 | 6.8 | +170 to +180 | 100% |
| 71 | 1 | 6.0 | +380 to +390 | 100% |
| 72 | 3 | 26.5 | +470 to +500 | 99% |
| 73 | 5 | 51.9 | +550 to +600 | 80% |
| 74 | 3 | 31.2 | +590 to +620 | 73% |
| 75 | 1 | 5.5 | +530 to +540 | 100% |
| 78 | 2 | 24.2 | +550 to +570 | 90% |
| 79 | 6 | 59.4 | +560 to +620 | 42% |
| 80 | 7 | 66.1 | +520 to +590 | 58% |
| 81 | 4 | 41.2 | +410 to +450 | 85% |
| 82 | 5 | 50.6 | +140 to +190 | 69% |
| 83 | 3 | 25.6 | +90 to +120 | 61% |
| 84 | 4 | 35.4 | +25 to +65 | 78% |
| 85 | 6 | 58.2 | -15 to +45 | 47% |
| 86 | 3 | 27.4 | -10 to +20 | 91% |
| 87 | 1 | 10.9 | +50 to +60 | 100% |
| 88 | 1 | 6.7 | +160 to +170 | 100% |
GIS 2–6 carry no NC-12, so they have no setback to report.

## Related figures

In `../` (`1984-start/dune-topo/`):

- `HAT_where_inserts_occur.png` — where the insert changes a domain and where it
  does not, and why; the scope change from v4
- `HAT_how_N_determined_GIS85.png`, `_GIS84.png` — the full measurement chain
- `HAT_dunelines_on_grid_GIS85.png` — both lines on the Barrier3D grid, all 90 domains
- `HAT_dunelines_on_DEM_GIS85.png` — both lines on the raw 1996-graft DEM, in EPSG:3725
