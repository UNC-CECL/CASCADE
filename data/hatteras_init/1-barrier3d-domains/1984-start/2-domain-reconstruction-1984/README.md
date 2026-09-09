# 2-domain-reconstruction-1984

**Where the 1984 seaward-row insert would touch the Barrier3D grid, and by how
many cells.** Scope only — no array is written, no elevation is fabricated.

Written by `scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/2-extent/HAT_report_row_insert_scope.py`.

## Layout (2026-09-09)

**Six steps, in the order the argument runs, the same names in the data, the
figures and the scripts** (`scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/<step>/`):
how far the dune line moved, how many rows, WHERE they go, what they contain,
the version built, what the model does. Root keeps this README,
`DUNE_TOPO_VERSION_GUIDE.md` and `figures/`. Paths are resolved by
`hat_topo_version.insert_scope_step(product, step[, sub])` and
`insert_figures_dir`; no script builds one by hand. Everything moved with
`git mv` (history follows); nothing was re-run.

| step | the question | data here | figures | scripts |
|---|---|---|---|---|
| `1-measurement/` | how far did the dune line move | `duneline-shift/` — the two crossings per profile | `1-measurement/` | `HAT_measure_duneline_shift.py`, three per-domain plotters |
| `2-extent/` | how many rows, which domains (placement-independent) | `footprint_1984_by_domain.csv` (**the table everything reads**), `footprint_1984_profiles.csv`, `HAT_footprint_1984.txt`, `row_insert_scope_by_domain.csv`, `HAT_row_insert_scope.txt` | `2-extent/island/` | `HAT_footprint_1984.py`, `HAT_report_row_insert_scope.py`, `HAT_plot_where_inserts_occur.py` |
| `3-placement/` | WHERE the rows go | `road_placement_check_1984.csv` / `.txt` (the road check); `imagery-review/` — the sheet with the verdicts and picks, and its two reports (the evidence that decides) | `3-placement/seaward/`, `behind-road/` (the two candidates), `road-check/`, `imagery-review/` | `HAT_verify_road_placement_1984.py`; `imagery-review/`: batch, window, summary |
| `4-fill/` | what the rows contain | `fill_copy_by_domain.csv`, `HAT_fill_copy_scope.txt`, `HAT_fill_options.txt` | `4-fill/` | `HAT_fill_copy_scope.py`, the fill plotters, the explainers |
| `5-build/` | the version built from 2 + 3 + 4 | a README pointing at `../dune-topo/v3/`, where the build writes | none | `HAT_build_footprint_version.py` |
| `6-result/` | what the model does with it | (none yet) | `6-result/` | the result and comparison plotters |

The extractor itself, `HAT_dune_topo_extractor.py`, stays where the resolver
finds it by path; `nodata_audit/` and `old_extractors/` are untouched.

## 2026-09-07 — the 1984 footprint, both directions (the new approach)

Written by `scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/2-extent/HAT_footprint_1984.py`. Scope
only: rows drawn BLANK, no array written, no model-facing CSV touched. Four rules
decided in interview the same day, each recorded in the script's docstring:
**symmetric** (rows added where the 1984 line lies seaward of the 1997 line,
existing rows removed where it lies landward), **paired median** (median of the
50 per-profile 1997−1984 differences, not the difference of two medians —
they disagree by one cell at 13 domains), **10 m rule** (N = trunc(shift/10):
a row only once a full cell of change is measured; the 0–9 m residual is a
column), and the **row-0 setback convention** (setback' = (road − row 0) + shift
per profile, unrounded; NOT road minus the 1984 line, which is ~20 m larger
because the digitized lines trace the toe ~19 m seaward of row 0).

| file | what |
|---|---|
| `footprint_1984_by_domain.csv` | per domain: paired shift (median, p10, p90), N, residual, rows now/after, setback now / new (row-0 convention, with spread) / derived / raw-lines, flags |
| `footprint_1984_profiles.csv` | the per-profile join the medians come from (both line crossings, row 0, road cell) |
| `HAT_footprint_1984.txt` | the report: rules, totals, the road table, the assumptions |
| `figures/3-placement/seaward/HAT_footprint_1984_grid.png` | the Barrier3D grid, added rows blank in red, removed rows hatched blue, new row 0 and NC-12 marked |
| `figures/3-placement/seaward/HAT_footprint_1984_plan.png` | plan view on the DEM with both dune lines and NC-12 1984 |
| `figures/3-placement/behind-road/HAT_row_insert_grid_behindroad.png` | the same rows placed behind NC-12 (advisor's placement): crest-to-road strip as measured, added rows behind the road as placed under the 1984 setback, removed rows from directly in front of today's road (2026-09-08); columns `insert_anchor`, `insert_row_behind_road` (the first row added or removed), `road_row_v2`, `road_row_new`, `road_cells_offset`, `rows_behind_road` in the table |
| `road_placement_check_1984.csv`, `HAT_road_placement_check_1984.txt` | **the placement check** (`HAT_verify_road_placement_1984.py`, 2026-09-08): the 1984 dune line and NC-12 on the lidar, along the profiles and perpendicular, against the row-0 setback and what v3 holds (CSV, rows, road rows); figures in `figures/3-placement/behind-road/` (island-wide) and its `rows-added/`, `rows-removed/` (per domain) |
| `fill_copy_by_domain.csv`, `HAT_fill_copy_scope.txt` | the COPY fill for the behind-road block (advisor, 2026-09-07): block = the N interior rows immediately landward of the insert point, copied cell by cell; window behind the crest where there is no road; per-domain window stats, seam and flags. Figures in `figures/4-fill/rows-added/`, `rows-removed/` (method, before/after per domain) and `4-fill/island/` |
| `../dune-topo/v3/` | **the built version** (2026-09-08): v2 + this footprint behind the road, copy fill, the 1984 setbacks (no floor). Built by `HAT_build_footprint_version.py`; test run under arm `behindroad-copy`; result in `figures/6-result/` |
| `figures/2-extent/island/HAT_footprint_1984_rows.png` | rows per domain, signed |
| `figures/2-extent/island/HAT_footprint_1984_shift.png` | the paired shift with p10–p90 and the rows kept |
| `figures/3-placement/seaward/HAT_footprint_1984_setback.png` | the road setback now and from the new row 0, per road domain |
| `figures/CAPTIONS.md` | the words under the figures (house style: no in-image titles) |
| `imagery_review_1984.csv`, `HAT_imagery_review_1984.txt` | **the aerial-imagery review** (`HAT_imagery_review_1984.py`, 2026-09-08): every changed domain and ten unchanged neighbours in the 1984 and 1997 photographs with both placements of the rows drawn on; the sheet carries the numbers and BLANK verdict columns to fill by eye. Figures in `figures/3-placement/imagery-review/rows-added/`, `rows-removed/`, `unchanged/`. `HAT_imagery_review_gui.py` is the window for filling it (picks, blink); `HAT_imagery_review_summary.py` / `HAT_imagery_review_summary.txt` the summary of the sheet as it stands |

Result: 29 domains gain 73 rows, 23 lose 47, 38 unchanged; largest +7 (GIS 80)
and −5 (GIS 63, 64, 66, 67). No new setback is negative; GIS 85/86 go from 0
(floored) to 44 and 12 m; GIS 16 falls to 1 m after losing four rows. The
earlier scope report and grid (`HAT_row_insert_scope.txt`, `row_insert_scope_by_domain.csv`,
`HAT_row_insert_grid.png`) were re-run the same day on this footprint, so they
agree with it; the add-only plan view was retired.

| file | what |
|---|---|
| `HAT_row_insert_scope.txt` | the report — symmetric since 2026-09-07, reads N from the footprint table |
| `row_insert_scope_by_domain.csv` | the same table, machine-readable, with the easting-frame cross-check |
| `figures/3-placement/seaward/HAT_row_insert_grid.png` | the Barrier3D grid as the model would hold it: added rows blank, removed rows hatched |
| ~~`figures/HAT_row_insert_plan.png`~~ | retired 2026-09-07; the plan view is `figures/3-placement/seaward/HAT_footprint_1984_plan.png` |
| `HAT_fill_options.txt` | **what to fill the inserted rows with** — the methodological options, worked at GIS 85 |
| `figures/superseded-layers/HAT_fill_options_grid_GIS85.png` | the two candidates and the control as Barrier3D domain views, with NC-12 at its actual road elevation |
| `figures/superseded-layers/HAT_fill_options_GIS85.png` | the fuller set as profiles — what the road sits on, how much is measured |
| `figures/1-measurement/rows-added/HAT_duneline_zoom_GIS83_87.png` | the two dune lines at true scale across GIS 83–87, each domain labelled with its offset and its row count |

## 2026-09-08 — the footprint against the aerial photographs (a review, not a measurement)

Hannah's colleagues asked to see the photographs before trusting where the rows
go. `scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/3-placement/imagery-review/HAT_imagery_review_1984.py` draws,
for each of the 52 changed domains and 10 unchanged neighbours (controls), the
same window in the **1984-09-19** and **1997-10-12** photographs of the USGS
Henderson release (`D:\Hatteras_GIS\Aerial`, doi 10.5066/P1CXBCDW, georeferenced
to the 2007 orthos, stated accuracy 1.2 m) and on the 1 m lidar, with the two
digitized dune lines, NC-12, interior row 0, the road rows as measured and as
placed, and **both placements of the rows**: v3's (behind NC-12 / in front of
it) and the seaward alternative (at the dune). A brightness strip along the 50
profiles sits under the panels so the sand-to-vegetation transitions read as a
curve. Nothing is picked or measured from the photographs by the script.

**The question it is built to answer** (Hannah): did the dune field actually
get narrower between 1984 and 1997, and from which side, so that the rows are
placed where the width was actually lost or gained. The verdict is written by
eye into `imagery_review_1984.csv` (`extra_width_was`: seaward_of_crest |
crest_to_road | behind_road | none; plus `dune_field_change`, `edge_moved`,
`placement_ok`, `confidence`, `notes`). Re-running the script keeps filled
verdicts. What is read off the photographs: the seaward vegetation line (the
toe the digitized lines trace) and the pavement itself, not the geojsons.

**The window.** `HAT_imagery_review_gui.py` (same folder) opens the same
domains one at a time in a Tk window: the two photographs side by side with
shared pan and zoom, every overlay a checkbox, the brightness strip under them,
the domain's numbers, and the verdict form with the vocabulary as drop-downs.
Save writes that domain's six verdict columns into `imagery_review_1984.csv`,
plus **the picks** (2026-09-09; Hannah: N is trusted, the question is WHERE the
rows go): three features per year, one click each on the photograph, the toe
(beach sand to dune vegetation), the back of the dune band (dune sand to flat
interior) and the road edge as it is in that photograph. From them the sheet
gets each feature's distance from interior row 0, the three bands per year
(dune band, back of dune to road, toe to road), their change, each feature's
shift in the footprint's sign, and, taking N x 10 m as given, where the 1984
width sat: `lost_dune_band_m`, `lost_back_to_road_m`, and the remainder
`lost_behind_road_m`. `band_suggests` names the largest share and is DERIVED;
the verdict columns stay the reviewer's. `reviewed_by` / `reviewed_at` stamp
each save. Measured by the reviewer, not by code; the batch script preserves
these columns when it re-runs. **Blink mode** puts both photographs in one panel
and flips them on Space, or on a timer. Arrow keys step through domains; the
next one is read in the background and tiles are cached under `~/.cascade/`
(outside the repo). Run the batch script first: the window reads the sheet it
writes.

**The summary.** `HAT_imagery_review_summary.py` (also the window's Summarize
button) reads the sheet as it stands and writes
`figures/3-placement/imagery-review/island/HAT_imagery_review_summary.png` and
`HAT_imagery_review_summary.txt`: the verdict per domain along the island
against the signed N; **the placement question** as stacked bars per measured
domain of where the 1984 width sat (dune band / back of dune to road / behind
the road) against N x 10 m, with the island-wide median shares by sign, which
is what the one placement rule reads off; and which feature moved. Runs on an
empty sheet and says so.

Decisions recorded in the script's docstring: 1984 + 1997 as the default pair
(1996 is in the release but judged poor; 1997 is the year the second line was
digitized from; `--years` takes any year 1978–2002 for a second pass on
uncertain domains); controls = unchanged neighbours of the changed runs, thinned
to ten; GIS 1–8 profile coordinates re-derived through the extractor's
`cell_to_map` and checked against `RoadOffset_1984_profiles.csv` to the
centimetre. The 1984 mosaic has gaps (south of GIS 63); the release's frames
under `1984_georef_TIF/Input` fill them.

## Decided 2026-09-04 — `--fill median`, built as `dune-topo/v5` (deleted 2026-09-07)

The fill is settled: keep every dry measured cell, give the water cells the
median of the block's own dry cells. Chosen under two constraints — the
extracted interior is not modified, and the simplest rule with the least
fabrication wins. `v5` carries it; `v4` (measured + floor) stays as the shipped
reference. Same footprint, same setbacks, 4765 of 4900 inserted cells from the
DEM. The argument and the one cost that was accepted with it — at GIS 85 and 86
the road sits inside the 1996 crest and the first bulldoze pass adds ~3 m to the
dune — are in `HAT_fill_options.txt`, section DECISION, and in
`../dune-topo/v5/README.md`.

## Everything for the seaward-row insert lives here

Consolidated 2026-09-03.

| | |
|---|---|
| `duneline-shift/` | **the measurement of N.** Moved here from the product root; resolved by `hat_topo_version.duneline_shift_dir`. Two scripts write into it |
| `HAT_row_insert_scope.txt` | where the insert lands and how many cells |
| `HAT_fill_options.txt` | what to fill the inserted rows with |
| `row_insert_scope_by_domain.csv` | the scope table, machine-readable |
| `figures/` | numbered by step (measurement, footprint, fill) plus the deleted layers' record and the frozen pair; `figures/README.md` maps every file |

What is deliberately NOT here: the topography itself. `dune-topo/v3` and `v4`
are the built inserts, and they stay with the other topography versions because
that is what `hat_topo_version` resolves over. This folder is the measurement
and the argument; `dune-topo/` is the arrays.

## The answer

**38 of 90 domains, 98 rows in total, largest 7 (domain 80).**

    N_metres = shift_1984 − shift_1997 = line_1997 − line_1984
    N_rows   = round(N_metres / 10), floored at 0

| N | domains |
|---:|---|
| 1 | 2, 6, 10, 13, 18, 44, 46, 47, 48, 52, 53, 55, 71, 75, 87, 88 |
| 2 | 4, 11, 12, 45, 78 |
| 3 | 49, 51, 54, 72, 74, 83, 86 |
| 4 | 3, 50, 81, 84 |
| 5 | 73, 82 |
| 6 | 5, 79, 85 |
| 7 | 80 |

N comes from `duneline-shift/duneline_retreat_1984_1997.csv`, in this folder — the same file
`HAT_insert_seaward_rows.py` reads under `--shift-source date`. Interior row 0
cancels out of that subtraction, so no assumption about where row 0 sits
survives into N, and the definitional offset between "a digitized dune line" and
"the model's row 0" cancels with it.

## This is not a proposal — it is already built

`dune-topo/v4` is the base plus exactly these rows. The report checks itself
against `v4/HAT_seaward_row_insert_audit.csv` and agrees in **all 90 domains**:
38 modified, 98 rows. It also reproduces the ten published block-scope values,
which `dune-topo/v3` still carries in its audit CSV.

So the footprint is settled. What is open is the **fill**, and whether to adopt
it at all.

## Reading the plan view

Two encodings, because one cannot carry it alone. A one-row insert is 10 m; on a
panel 3 km wide that is a third of one percent — sub-pixel. So each **domain box
is shaded by N**, which is what answers "which domains" at a glance, and the
**true-scale band** is drawn on top of it: the 1997 dune line offset seaward by
exactly N × 10 m. The band is honest about size and mostly needs a zoom to see;
it only reads at island scale where N is large, along domains 79–85.

The band is anchored on the **1997** line because the rows go in ahead of current
interior row 0, and row 0 was picked from a surface whose dune is the 1996/1997
one. It stops at N × 10 m rather than at the 1984 line itself — the gap between
the band's outer edge and the red line is the rounding in N, drawn rather than
hidden.

## The measurement at GIS 85, at true scale

`HAT_duneline_zoom_GIS83_87.png` is the five-domain span around GIS 85, drawn at
equal aspect with the same 300 m cross-shore crop as the island-wide zooms in
`0-elevation/2009-2014-1996-duneline/figures/`. Each box carries its measured
offset **and** the row count that offset becomes:

    87  +10 m → +1 row      86  +27 m → +3 rows      85  +59 m → +6 rows
    84  +36 m → +4 rows     83  +27 m → +3 rows

Across this reach the 1984 line sits consistently seaward of 1997 — the one
place on the island where the sign is stable over five neighbours *and* large.
Regenerate it with:

```
python scripts/input_prep/0-elevation/3-figures/HAT_plot_duneline_offset.py \
    --zoom 83-87 --zoom-out <path>
```

`--zoom` reads the existing offset table rather than re-measuring, so it takes
seconds, and it shares the crop rule, line styling and equal aspect with the
three standard reaches — a custom zoom cannot quietly differ from them.

## The fill — measured footprint, asserted contents

`HAT_fill_options.txt` works the choice at **GIS 85** (N = 6 rows, 60 m), the
domain the problem is named after. The premise is not negotiable: **no survey
covers ground that had eroded away by 1996**, so every option is an assertion
about 1984, not a measurement of it. The DEM does hold cells there, but they are
the 1996 surface — a later, lower landform in the same place.

What does not change between options: N, everything landward of the block, the
dune arrays, and NC-12's offset. What changes is **the ground under the road**,
which is what `bulldoze()` scrapes and what the setback is read against.

**Two candidates and a control** in the grid figure. Options are named rather
than lettered here — the two figures no longer carry the same set, so a letter
means different things depending on which one you're holding (see *The two
figures differ* below).

| fill | ground under the road | block mean | % measured |
|---|---|---|---|
| matched backdune — today's profile, moved seaward | 1.38 / 1.21 m | 1.89 m | 0% |
| **measured + median fill** — `--fill median` | 3.17 / 4.96 m | **2.17 m** | **91%** |
| raw DEM, no floor — *control* | 3.17 / 4.96 m | 1.98 m | 100% |

**Matched backdune** copies the real near-dune profile 60 m seaward, so it
carries genuine alongshore texture and puts the road on backdune. But **it
duplicates the crest**: v2's row 0 at GIS 85 *is* the 1996 crest, so the copy
leaves two interior ridges on top of Barrier3D's own dune array. It counts as 0%
measured despite being real cells — they are measurements of a *different place*.
Not a `--fill` option yet.

**Measured + median fill** keeps every dry cell exactly as measured and fills
only the cells at or below MHW, with the median of the block's own dry cells.
One guard, one constant. It is **internally consistent where the shipped rule is
not**: `measured` argues "1996 is a lower bound on 1984" and then *raises* 44% of
the block above the DEM to the platform, which asserts more than the data
supports. Only 28 of 300 cells at GIS 85 need the constant, 25 of them in the
seaward-most row, and the road never touches it. Its cost: it admits the 1996
beach ramp, so row 1 sits below row 0's fill and the block carries a dip two
cells inside the island. Added as `--fill median` on 2026-09-03 and verified on
all 38 domains.

**Raw DEM is drawn to be rejected.** Of the 300 cells in GIS 85's block, **28
(9.3%) are at or below MHW** — half the seaward-most row is water. `median`
differs from this control *only* in replacing those 28.

### The road is not part of the fill question

`bulldoze()` does `np.zeros(...) + road_ele`, so every road cell carries one
value regardless of the ground beneath. The grid figure now draws it that way —
a flat band at **0.807 m MHW**, identical in every candidate panel, because the
road elevation is an independent forcing no fill can move. What the fill still
decides is how much sand gets scraped off and handed to the dune.

**Which road-elevation file — this is a trap.** `HATTERAS_ROAD_ELEVATION_FILE`
resolves to `4-mgmt-forcing/road_elevation/RoadElevation.csv` (**0.807 m** at
GIS 85). A second file, `.../dunestart_offset/1984/RoadElevation_1984_dunestart.csv`,
reads **1.833 m** — a metre higher, because it samples along the *1984*
alignment, which at GIS 9–15 and 84–87 now lies **under the foredune** and so
returns dune, not roadbed.

### Median vs the shipped rule, all 38 domains

| | measured share of the block |
|---|---|
| shipped (measured + floor) | min 42.3%, median 93.7%, mean 86.1% |
| **new (measured + median)** | min 79.0%, median 100%, mean 98.8% |
| domains at 100% measured | 15 → **32** |

GIS 85 is close to the worst case for the shipped rule — island-wide it is
already 86% measured. **And the island barely moves:** setbacks identical in all
33 domains with a road, land rows identical in all 38, mean land elevation
differing by a median of 0.001 m and at most **−0.043 m (GIS 85)**. The two rules
disagree about *justification* far more than about the island.

### No longer considered

**flat backdune** and **measured + floor** were dropped from the grid figure on
2026-09-03. `measured + floor` is **still what is on disk** (v3, v4) — adopting
either candidate is a change *from* it, so it stays documented in the report.
Its objection: the cells it keeps at GIS 85 are the **1996 dune face**, handing
the model a road embedded in 3–5 m of ridge that had not formed there in 1984.
`taper` was dropped earlier the same day; `--fill taper` still exists in the
build script.

**Still undrawn:** an **alongshore analogue** (a neighbour's surviving backdune —
matched backdune's idea without the crest duplication) and a **mass-conservative
reconstruction**, the only option that would be derived rather than asserted.

### The two figures differ

The grid was trimmed to the two candidates plus the control; the **profile figure
still carries the fuller set**, including flat backdune and the shipped rule. So
`c` and `d` mean different things in the two figures. That is a defect, not a
design — left as-is only because trimming the profile wasn't asked for, and it is
the only place the shipped rule is still drawn beside its replacement.

## Three checks worth knowing about

**Independent frame.** `0-elevation/2009-2014-1996-duneline/` measures the same
1984-minus-1997 quantity in raw easting — no extractor `c0`, no shear, no row 0.
Identical N in 84 of 90 domains, the other 6 differ by one row, mean absolute
difference 1.77 m. That bounds the **frame**, not the lines; both readings come
from the same two geojsons.

**45 domains are negative** — the 1984 line lands *landward* of the 1997 line.
Every one is floored to N = 0 and left alone. A symmetric rule would have
implied 62 rows of **removal**, i.e. deleting surveyed interior cells on the
strength of apparent progradation. That is not done, and the asymmetry is a
choice recorded here rather than an oversight.

**The two lines may not be the same feature.** In the true-scale zooms at
`0-elevation/2009-2014-1996-duneline/figures/`, the two big excursion reaches
look like both lines straddling one ridge from opposite sides and swapping which
side they take — 1984 landward in 62–68, seaward in 78–85. That is what a
definitional difference looks like, not one-directional retreat, and it bears
directly on N. Unresolved.

## Base version

The figure and the `rows_now` column read whatever `hat_topo_version.topo_dirs`
resolves — currently **v2**, because that module reads the extractor's `VERSION`
literal *before* `dune-topo/CURRENT`, and CURRENT still says `v1`. The report
header records which was used. **N does not depend on the choice**; the drawn
interior does. Override with `--base`.

Written 2026-09-03.
