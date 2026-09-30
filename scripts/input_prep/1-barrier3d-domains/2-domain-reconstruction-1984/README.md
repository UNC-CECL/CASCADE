# 2-domain-reconstruction-1984 - the 1984 domains rebuilt from the 1996-based DEM

v2 -> v3: where the interior must grow or shrink so row 0 stands on the 1984
dune line, where the added rows go, what they contain, the build, and the
hindcast on the result. One subfolder per step, in the order the argument runs.
The stage map is in `../README.md`; the data mirror this tree under
`data/hatteras_init/1-barrier3d-domains/1984-start/2-domain-reconstruction-1984/`.

```
1-measurement/
    HAT_measure_duneline_shift.py       how far a dune line sits from interior row 0; line minus line
    HAT_plot_dunelines_on_dem.py        the lines on the raw DEM
    HAT_plot_dunelines_on_grid.py       the lines on the Barrier3D grid, all domains
    HAT_plot_how_N_is_determined.py     the chain from lines to N, one domain
2-extent/
    HAT_footprint_1984.py               the footprint table: rows added or removed, the 1984 setback
    HAT_plot_where_inserts_occur.py     the relocation blocks and the setback, from the table
    HAT_report_row_insert_scope.py      the scope report and grid figure
3-placement/
    HAT_verify_road_placement_1984.py   is NC-12 where the lines say, in three frames
    imagery-review/                     the photographs: batch figures, the GUI, the quick review, the summary
4-fill/
    HAT_fill_copy_scope.py              what the copy fill puts in the added rows
    HAT_insert_seaward_rows.py          the earlier seaward-row layers (guarded since 2026-09-07)
    HAT_plot_fill_options*.py           candidate fills (guarded)
    HAT_plot_insert_explainer*.py       how the insert is built (guarded)
5-build/
    HAT_build_footprint_version.py      the footprint as a dune-topo version (v3)
    HAT_plot_version_figures.py         a built version's figure set
6-result/
    HAT_run_version_pair.py             one scenario on two versions, identically
    HAT_compare_versions.py             v2 against v3: relocations and geometry
    HAT_plot_footprint_result.py        the first result, at the road
    HAT_version_pair_gif.py             v2 beside v3 through time
    HAT_version_pair_report.py          v2 against v3 in one report
    HAT_plot_b3d_grid.py, HAT_plot_insert_three_scales.py,
    HAT_plot_method_compare.py, HAT_plot_seaward_insert_compare.py   version comparisons
```

What not to trust: every layer built on v2 (v3-v8 of the earlier numbering)
and every run on modified topography were deleted on 2026-09-07. The scripts
that drew them are guarded and stop before drawing until a layer is rebuilt;
their figures live in `figures/superseded-layers/`.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 1-measurement/HAT_measure_duneline_shift.py

How far seaward of the model's interior row 0 does a digitized dune line sit, per domain?

From the script's original header:

```text
How far seaward of the model's interior row 0 does a digitized dune line sit?

WHY
---
The 1984-start topography is a 1996 ALACE beach on a 2009 backdune, and the
NC-12 lines it is measured against are 1984-vintage. Where the island migrated
far enough between 1984 and 1996, the 1984 roadbed ends up SEAWARD of the 1996
dune crest and the setback goes negative (GIS 85: -10 m, floored to 0, which
makes roadway_manager relocate the road in year 1). This measures the offset
that gap represents, so HAT_insert_seaward_rows.py can act on a number rather
than on a target.

METHOD -- and why it is exact rather than approximate
A domain clip is north-up, so one extractor "profile" is a raster row of
constant y. Intersect the dune-line geometry with that horizontal line, convert
the crossing's easting to a cross-shore cell index, and difference it against
interior row 0 for the same profile.

The frame comes from the extractor itself -- its own c0 and per-profile shear,
through the same inversion cell_to_map documents -- so this is measured in the
frame CASCADE indexes, not in a re-derived one. That distinction is the whole
reason the legacy RoadSetback numbers and the dunestart numbers disagree by a
median 38 m: the legacy pass measured raw and ocean-first, unstraightened.

SIGN: positive = the dune line lies SEAWARD of interior row 0, i.e. the number
of cells the island has to move seaward for row 0 to land on the digitized line.

THE CONTROL THAT MAKES THE 1984 NUMBER READABLE
Run it on the 2004 pair (2004 dune line against 2004-start, whose DEM is 2009)
and the stable mid-island comes out at +1.4 / -4.1 / -3.6 m on GIS 40/50/60 --
under half a cell. So a digitized dune line and the extractor's interior row 0
are the SAME feature to within the grid, and a large 1984 number is a date
difference rather than a definitional one.

That control matters because it settles a confound the road-offset work had
recorded as unresolvable without a same-year DEM. It is resolvable without one:
the 2004 line is close enough in date to the 2009 surface that its residual IS
the feature term, and the feature term is ~0.

Values well above zero at GIS 10/11/84/85 in the 2004 pass are NOT method error
-- those are the relocation blocks, where five years of hotspot erosion is real.
Read the mid-island domains for the method check.

INPUT   D:\Hatteras_GIS\Dunelines\duneline_<year>.geojson   (EPSG:26918)
        domain-clips-1m/domain_<N>/resampled_domain_<N>.tif (EPSG:3725)
        the extractor's own picks for the resolved version

OUTPUT  hat_topo_version.duneline_shift_dir(<product>)/duneline_shift_<year>.csv
        (1984-start: .../1984-start/2-domain-reconstruction-1984/1-measurement/duneline-shift/)
            one row per domain: median/p10/p90 shift, row 0, dune-line cell

USAGE
    python HAT_measure_duneline_shift.py                  # 1984, all domains
    python HAT_measure_duneline_shift.py --year 2004      # the control
    python HAT_measure_duneline_shift.py --domains 84,85,86
```

Notes that were in the code:

```text
The REPO copy wins. The lines used to be read straight off D:\Hatteras_GIS,
which is not version-controlled, not present on another machine, and not
something a run can record the state of. Hannah placed a curated set inside the
repo on 2026-09-02; that is now the source, and the external drive is only a
fallback so older invocations keep working.
```

```text
GEOREFERENCE AGAINST THE PRODUCT'S OWN RASTER, not the clip tree.

This was wrong until 2026-09-03. It used
1-barrier3d-domains/domain-clips-1m/domain_<N>/resampled_domain_<N>.tif
which is a DIFFERENT DEM product from the one the 1984 npy arrays were
exported from: at GIS 85 its origin is offset +0.50 m in easting and its
elevations differ (max 6.80 m against 6.19 m). Only the transform is used
here, so the cost was 0.05 cells -- N at GIS 85 went 5.87 instead of 5.82 and
rounded to 6 either way -- but mixing the extractor's frame with another
grid's transform is not a thing to leave in place.

Resolved through hat_elevation_products so the period -> product pairing is
the same single definition every other reader uses.
```

```text
map x of cross-shore cell 0 on this profile, inverting the
extractor's chain exactly as cell_to_map does
```

```text
a meandering line can cross one profile more than once; take the
crossing nearest row 0 rather than the first, which would silently
pick a sound-side meander on the wide domains
```

```text
STAMP THE TOPOGRAPHY. These numbers are measured against interior
row 0, so they are only valid for the extraction that produced it.
Without the stamp a v1-era duneline_shift_1984.csv sat unmarked
beside v3 files -- 65.9 m against the correct 75.0 m at GIS 85 --
and nothing on disk said which was which.
```

```text
1996 and 1967 are not hindcast period starts, so they have no product
of their own. They are measured in the frame of the period being
corrected -- 1984-start -- which is also the only frame in which a
difference against the 1984 line is meaningful.
```

```text
LINE MINUS LINE, in one frame. Each line is first measured against the
same interior row 0, then the two are differenced -- so row 0 drops out
algebraically and no assumption about it survives into the answer. That
is the whole point: the row-0 reference is the part that cannot be
validated, and differencing removes it rather than bounding it.
```

```text
shift = row0 - line, so (row0 - lineEARLY) - (row0 - lineLATE)
= lineLATE - lineEARLY. Negated here so the stored number reads as
RETREAT: positive = the later line is LANDWARD of the earlier one,
i.e. the dune moved landward by that many metres. Storing the raw
difference would put a negative sign on ordinary erosion and invert
every consumer that reuses this file expecting the shift_m_median
convention of the un-differenced output.
```

```text
NOT symmetric between products - 1984-start lives under
2-domain-reconstruction-1984/. hat_topo_version.duneline_shift_dir owns that.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_line()`**

```text
The dune line, dissolved and reprojected into the domain grid's CRS.

Reprojected through pyproj from the file's own declared CRS. NAD83 and
NAD83(NSRS2007) differ by centimetres here, but the transform is done rather
than assumed away -- the same rule the 1996 aerial chips follow.
```

</details>

### 1-measurement/HAT_plot_dunelines_on_dem.py

The dune lines on the raw DEM in map coordinates, and what the DEM says at each line, island-wide.

From the script's original header:

```text
The dune lines on the RAW DEM, in map coordinates, before any Barrier3D
processing -- and what the DEM says at each line, island-wide.

WHY THIS FIGURE EXISTS
    The 1984 dune line is SUBMERGED in the surveyed surface. That is not a
    problem with the measurement, it is the measurement: by 1996 the island had
    retreated past where its 1984 dune stood, so the ground at that position is
    now below MHW and no survey covers it. How far offshore the line sits is
    what sets how many interior rows have to be added.

    Drawing it on a processed Barrier3D grid cannot show this, because that grid
    starts at the water trim and is expressed in cells. This draws the geometry
    on the product's own raster, in metres, in the raster's own CRS.

    DEM: 0-elevation/2009-2014-1996 (2-resampled-10m) -- the 1996 ALACE graft,
    which is the product the 1984-start arrays are exported from. Resolved
    through hat_elevation_products, not by joining strings.

WHAT N IS, IN THESE TERMS
    N is NOT how far offshore the 1984 line is. That distance -- row 0 to the
    1984 line -- also contains the offset between a digitized line and the
    model's interior row 0, which the 1997 line measures separately:

        offshore distance (total)  =  N  +  (line vs row 0)

    Inserting the full offshore distance would put row 0 on the digitized line,
    a light/dark break at the dune toe, when row 0 is one cell landward of the
    crest. At GIS 85 that would over-insert by ~1.7 cells.

USAGE
    python HAT_plot_dunelines_on_dem.py [--domain 85]
```

Notes that were in the code:

```text
---- (a) map view ---------------------------------------------------
Elevation in the house classes, relative to MHW, so the water break in
the colour scale is the MHW contour the panel is about.
```

```text
an inset beside the map, so the bar is the map's height and not the
gridspec row's (the map is aspect-equal and shorter than its row)
```

### 1-measurement/HAT_plot_dunelines_on_grid.py

The digitized dune lines drawn on the Barrier3D grid, and the same measurement across all 90 domains.

From the script's original header:

```text
The digitized dune lines drawn ON the Barrier3D grid, and the same measurement
across all 90 domains.

WHY PUT THE LINES ON THE GRID
    N is a difference between two digitized lines, expressed in model cells.
    Every earlier figure showed that as numbers or as a profile. Drawing the
    lines on the cells they are measured in makes three things checkable by eye:

      * that the lines fall where the topography says a dune line should --
        the 1997 line should track the seaward face of the surveyed dune;
      * that interior row 0 sits LANDWARD of both, which is why the raw
        measurement carries a definitional offset at all;
      * that the band between the two lines -- the date term, which is N -- is
        a coherent alongshore feature and not per-profile noise.

    Both lines are drawn per profile, at the fractional cell where the geometry
    actually crosses that profile's raster row, so the sawtooth is real: it is
    the per-profile shear of the north-up clip, and it is present in both lines
    identically, which is why it cancels in the difference.

USAGE
    python HAT_plot_dunelines_on_grid.py [--domain 85]
```

Notes that were in the code:

```text
Resolved through hat_topo_version.duneline_shift_dir - ONE definition
of a path that eight scripts used to build by hand. Moved under
2-domain-reconstruction-1984/ on 2026-09-03.
```

```text
The bar is attached to (a) and describes only that panel; (b) and (c)
are charts.
```

```text
no legend: the four handles are those of (a), whose legend sits
directly above this panel
```

```text
Do not overstate this. The feature term is TIGHT over most of the island
(IQR +14.5 to +26.2 m) but it is not constant: it spikes to 130-145 m
around GIS 35 and 63-68, the reaches where the date term is strongly
negative -- i.e. where the shoreline prograded and the two lines are on
opposite sides of row 0. Those are the domains where the differencing
argument is weakest, and the caption says so rather than averaging
them away.
```

### 1-measurement/HAT_plot_how_N_is_determined.py

How the number of inserted rows, N, is measured: one domain, the whole chain.

From the script's original header:

```text
How the number of inserted rows, N, is measured. One domain, the whole chain.

THE QUANTITY WANTED
    How far the dune line moved between 1984 and the surveyed surface, in 10 m
    cells. That distance is how far interior row 0 has to move seaward for the
    1984 roadway to sit its true distance behind the dune.

WHY IT IS A DIFFERENCE OF TWO LINES AND NOT ONE MEASUREMENT
    The obvious measurement -- 1984 dune line against the extractor's interior
    row 0 -- confounds two things:

        row 0 - line_1984  =  (how far the island moved)          DATE
                           +  (digitized line vs the model's row 0)  FEATURE

    The feature term is not small. The 1997 line, measured the same way against
    the same row 0, sits +16.2 m seaward of it island-wide (IQR +12.8 to +21.0)
    -- a near-constant offset, which is what a definitional difference looks
    like. Island-wide it accounts for ~85% of the naive number.

    Differencing two digitized lines cancels it exactly:

        (row0 - line_1984) - (row0 - line_1997) = line_1997 - line_1984

    Row 0 drops out algebraically, so no assumption about where row 0 sits
    survives into N. And because the same person digitized the same feature from
    the same kind of imagery at both dates, the definitional term cancels too.

    N = round( median over profiles / 10 m ),  floored at 0.

USAGE
    python HAT_plot_how_N_is_determined.py [--domain 85]
```

Notes that were in the code:

```text
Resolved through hat_topo_version.duneline_shift_dir - ONE definition
of a path that eight scripts used to build by hand. Moved under
2-domain-reconstruction-1984/ on 2026-09-03.
```

```text
(a) and (b) side by side, (c) the bar below them. The arithmetic that
used to be a monospace panel (d) is in the caption.
```

### 2-extent/HAT_footprint_1984.py

Where the interior must grow or shrink for row 0 to stand on the 1984 dune line, and the 1984 setback that results.

From the script's original header:

```text
Where the Barrier3D interior has to grow or shrink for its seaward row 0 to
stand where the 1984 dune line stood - per domain, in whole 10 m cells, in
BOTH directions - and what the 1984 road setback becomes when it does.

SCOPE ONLY. No array is written, no elevation is fabricated, no model-facing
CSV is touched. The rows are drawn BLANK: this says where they land and how
many, so the fill can be argued separately against a known footprint.

DECIDED WITH HANNAH, 2026-09-07 (each one changes what the model would ingest)
    symmetric     rows are ADDED where the 1984 line lies seaward of the 1997
                  line and REMOVED where it lies landward. The earlier work
                  (layers v3-v8, deleted 2026-09-07) floored negatives to 0.
    paired median the per-domain shift is the MEDIAN OF THE 50 PAIRED
                  PER-PROFILE DIFFERENCES, line_1997 - line_1984 on the same
                  raster row. duneline_retreat_1984_1997.csv stores the
                  difference of two medians instead; the two disagree by one
                  cell at 13 domains (GIS 80: 7 vs 6). The paired form also
                  gives a p10-p90 spread, which the stored file cannot.
    10 m rule     N = trunc(shift / 10 m): a row only when a FULL cell of
                  change is measured. Understates by 0-9 m at every changed
                  domain, always toward less change; the residual is a column.
    row-0 setback the new setback keeps the model's reference (metres landward
                  of interior row 0, which is one cell behind the picked crest):
                      setback_new(p) = (road(p) - row0(p)) + (line97(p) - line84(p))
                  per profile, UNROUNDED, median per domain. NOT road minus the
                  1984 dune line: the digitized lines trace the toe, ~19 m
                  seaward of row 0 (IQR 15-28 m over the road domains), so that
                  reading is ~20 m larger everywhere and would silently change
                  the convention every run so far has used. It is kept as a
                  record column, `setback_raw84_m`, and nothing reads it.

TWO PLACEMENTS OF THE SAME ROWS (the second added 2026-09-07 evening)
    anchor = dune   rows go in at the seaward edge, between the dune and current
                    row 0, so row 0 lands on the 1984 dune line; the setback
                    becomes setback_new_m. The footprint above.
    anchor = road   Hannah's advisor: keep the strip from the crest to the road
                    AS MEASURED and put the rows BEHIND THE ROADWAY ROWS. The
                    roadway in the model is two straight rows at one setback
                    per domain - road_start = int(setback / 10) from row 0,
                    ROAD_ROWS = 2 - and the setback the model gets is the 1984
                    one (setback_new_m, no floor), so the block goes in at
                        insert_row_behind_road = int(setback_new_m / 10) + 2
                    cells landward of row 0, directly behind the road AS PLACED
                    (2026-09-08; until then it hung off today's setback and the
                    re-set road landed on it). NOT behind the GIS mask's
                    landward-most cell: that edge wanders 3-12 cells along a
                    domain (kept as `road_land_max_cell` for the record) but the
                    model never sees it. REMOVALS come out of the interior IN
                    FRONT of the road (Hannah, 2026-09-08): the |N| rows directly
                    seaward of today's roadway rows go,
                        rows int(setback_v2/10) - |N| .. int(setback_v2/10) - 1,
                    so the road's own cells and everything behind them are
                    kept as measured and the crest-to-road strip shortens at
                    its road end. The model then places NC-12 at the 1984
                    setback, which lands on the old pavement's first row or
                    the row seaward of it (`road_cells_offset`, 0 or 1: the
                    1984 setback truncates to a cell, the row count is exact).
                    Row 0 and the dune stay put; the road sits on
                    measured cells at its 1984 distance from the crest; the
                    added width is behind it. GIS 85: road rows 4-5, block 6-10. Domains with no model road (GIS 1-5, 8) use the
                    CREST ROW as the anchor instead (anchor "crest", 2026-09-08):
                    the crest is the largest alongshore-median elevation in the
                    first CREST_SEARCH_ROWS interior rows, and the block goes in
                    at crest + 1, so the crest stays at the front and the copy
                    fill takes what follows the block, as behind the road.
                    N is identical in both; only where the rows sit differs.
                    The missing ground was lost from the OCEAN side; this books
                    it on the sound side, which restores 1984 width but not the
                    1984 position of either edge.

WHAT IS ASSUMED (and cannot be checked from these files)
    * the 1984 and 1997 lines trace the SAME feature. 1997 carries metadata
      saying "light/dark elevation break"; 1984 carries none.
    * one integer per domain. The interior is rectangular, so the alongshore
      median stands for 50 profiles and the spread inside a domain is lost.
    * 1997 stands for 1996. The surface is 1996 ALACE; the line is a year later.
    * the 1984 dune crest equalled the 1996 one - the dune array is not
      re-estimated.
    * removal deletes SURVEYED rows: the |N| directly seaward of today's
      roadway rows (the road end of the crest-to-road strip). The 1996
      foredune and the cells behind the road stay.
    * only the island width moves. The ocean shoreline is the shoreline-offset
      input; rows at the dune move the bay edge.
    * the 1984 road line is the 1978 export (deliberate, recorded elsewhere).

INPUTS (all already on disk; nothing is re-measured against GIS here)
    2-domain-reconstruction-1984/1-measurement/duneline-shift/duneline_shift_{1984,1997}_profiles.csv
        per (domain, profile): the line's crossing as a cross-shore cell and
        interior row 0, both in the extractor's own c0/shear frame
    4-mgmt-forcing/road_offset/dunestart_offset/measured/1984/RoadOffset_1984_profiles.csv
        per (domain, profile): the road's seaward cell and interior row 0, same
        frame. Row 0 is asserted identical across the three files.
    4-mgmt-forcing/road_offset/dunestart_offset/measured/1984/RoadOffset_1984_domains.csv
        the setback the model currently receives (setback_model_m) and flags
    1984-start/dune-topo/<CURRENT>/topography   rows_now, and the grid figure

OUTPUTS  2-domain-reconstruction-1984/
    footprint_1984_by_domain.csv     one row per domain (the audit table)
    footprint_1984_profiles.csv      the per-profile join the medians come from
    HAT_footprint_1984.txt           the report
    figures/2-extent/HAT_footprint_1984_rows.png       rows per domain, signed
    figures/2-extent/HAT_footprint_1984_shift.png      the paired shift with spread and the rows kept
    figures/3-placement/seaward/HAT_footprint_1984_grid.png     the grid, current frame, both signs
    figures/3-placement/seaward/HAT_footprint_1984_plan.png     plan view, both lines, NC-12
    figures/3-placement/seaward/HAT_footprint_1984_setback.png  the road setback now and from the new row 0
    figures/CAPTIONS.md              the words that go under the figures

USAGE
    python HAT_footprint_1984.py            # everything
    python HAT_footprint_1984.py --no-plan  # skip the slow DEM panel
```

Notes that were in the code:

```text
The step's root holds the placement-independent figures (rows, shift); the
seaward/ subfolder the figures that assume the rows go in at the seaward
edge (the grid in the current frame, the plan view, the new setback).
```

```text
Colours: the RdBu pair the dune-line figures use. Red is 1984 / seaward /
ground ADDED; blue is 1997 / landward / ground REMOVED. Same meaning on every
panel of every figure here.
```

```text
+ = the 1984 line lies SEAWARD of the 1997 line = the island retreated
= rows to ADD. Cells grow landward, so seaward is the smaller index.
```

```text
--- the behind-the-road placement (advisor's suggestion) -----------
Anchored on the MODEL's road AS PLACED under the 1984 setback (Hannah,
2026-09-08): two straight rows at int(setback_new/10) from row 0, the
block directly behind them. Until then it hung off the road at TODAY's
setback, and once the setback moved to its 1984 value the model's road
landed on the block. The GIS mask's landward-most cell is recorded
beside it but does not place the block.
```

```text
removal: the |N| rows directly SEAWARD of today's roadway rows
(Hannah, 2026-09-08). The old pavement lands at r_v2 - |N|;
the model places the road at int(setback_new/10), 0 or 1 row
seaward of that because the 1984 setback truncates.
```

```text
what the model receives TODAY (floored, drowning-relocated), for the
before/after panel; NaN where the road is outside the managed span
```

```text
Each panel's window is the full extent of its thirty domain BOXES plus
a pad - not the strip the lines occupy, as the dune-line figure crops.
The boxes are the subject here (they carry N), so every one is shown
whole, sound side included (Hannah, 2026-09-07).
```

```text
drawn at the printed width: the panels are equal-aspect, so their height
follows from the page width and the windows' shapes
```

```text
the road AS PLACED under the 1984 setback: the pavement's
landward edge moved landward by (setback_new - measured);
the block runs landward from there
```

```text
removal: the |N| rows directly seaward of today's pavement
(cross-shore grows landward = west, so seaward is +x)
```

```text
no road: behind the crest row, i.e. landward of the 1997
line (the map proxy for row 0), for both signs
```

```text
Labels inside each box's landward (sound-side) edge: every changed
domain with its N, in its colour; every fifth unchanged domain as a
plain number, as the reference does. Inside the box rather than at
the panel edge because the boxes are staggered across the panel.
```

<details><summary>Function notes (the original docstrings)</summary>

**`fig_grid()`**

```text
Every interior, every row, in ONE frame: current interior row 0 is y = 0 on
every domain. Existing cells are drawn in the project's elevation classes,
so the dune ridge, the backbarrier flat and the sound-side marsh read as
what they are; added rows sit ABOVE 0 (between the new row 0 and the old),
blank because no fill has been chosen; removed rows are the existing rows
0..|N|-1, hatched. The black tick is where interior row 0 ends up - the
1984 dune line, one cell behind the crest - and the dark bar is NC-12 at
its measured position, which does not move.
```

**`fig_rows()`**

```text
Rows per domain, signed, on its own - the communities banded along the
axis so a domain can be placed without the map.
```

**`fig_plan()`**

```text
The footprint in plan view, in the layout of the dune-line figure
HAT_duneline_offset_lines_island_3panel.png (2026-09-07, Hannah's request):
three panels of thirty domains, each cropped to the strip the two dune lines
and NC-12 occupy, a grey hillshade backdrop with no readable elevation, no
coordinate ticks (scale bar and north arrow instead), the communities as a
bracket in the ocean margin, domain numbers on the landward edge.

Two encodings, because one cannot work alone at island scale: each domain
box is shaded by N (red added, blue removed), and the TRUE-SCALE band is
drawn on top - the 1997 dune line offset by N x 10 m seaward (add) or
landward (remove). The 1997 line is the map-space proxy for the existing
array's seaward edge; the digitized line itself sits ~19 m seaward of row 0.

anchor="road" (2026-09-07, the advisor's placement): added rows hang off
the LANDWARD edge of NC-12 as placed under the 1984 setback (the 1984
centreline offset 10 m landward, then by the shift) and run landward by
N x 10 m; removed rows hang off the SEAWARD edge of today's pavement and
run seaward by |N| x 10 m (2026-09-08: removals come out of the interior
in front of the road). Domains without a model road keep the dune anchor,
as the placement does.
```

</details>

### 2-extent/HAT_plot_where_inserts_occur.py

Where the 1984 footprint changes the domain and where it does not, zoomed on the two relocation blocks.

From the script's original header:

```text
Where the 1984 footprint changes the domain and where it does not -- zoomed on
the two relocation blocks. Both signs. TWO figures since 2026-09-07 (they
were panels of one):

    HAT_where_inserts_occur_blocks.png   (a) (b) the two relocation blocks,
        GIS 9-14 and 84-87: the measured paired shift per domain with its
        p10-p90, and the rows it becomes under the 10 m rule, +N added / -N
        removed / no change.
    HAT_where_inserts_occur_setback.png  the NC-12 setback at those domains,
        as the model receives it now and from the new row 0 (with p10-p90).
    (an island-wide panel was drawn here too, and RETIRED the same day: it
     showed the same quantity as HAT_footprint_1984_shift.png, which is the
     canonical island-wide view - Hannah's call, 2026-09-07)

REWRITTEN 2026-09-07. This used to contrast two layers (block scope v3 against
island scope v4, both add-only) through their audit CSVs. Those layers were
deleted that day and the footprint became symmetric, so the figures read ONE
table, `2-domain-reconstruction-1984/2-extent/footprint_1984_by_domain.csv`, written by
HAT_footprint_1984.py.

THREE REASONS A DOMAIN IS UNCHANGED - now only one
    Under the old block scope "unchanged" meant either "measured, no cell
    needed" or "never asked". With island-wide scope and a symmetric rule there
    is one reason left: |shift| is under a full 10 m cell.

THE INTERVAL CAVEAT, DRAWN
    The dune lines are 1984 and 1997 -- 13 years. The DEM surface at row 0 is
    1996 ALACE, so the interval wanted is 12 years. The green triangles in the
    blocks figure show N if the measurement were scaled 12/13; recorded, not
    corrected -- scaling would assume steady change across 13 storm years.

USAGE
    python HAT_plot_where_inserts_occur.py
```

### 2-extent/HAT_report_row_insert_scope.py

Which domains would have interior rows added or removed, and how many, drawn on the Barrier3D grid.

From the script's original header:

```text
Which domains would have interior rows ADDED behind the dune, which would have
rows REMOVED, and how many - drawn on the Barrier3D grid the way the model
would hold it, with a report and a per-domain table.

SCOPE ONLY. Nothing is written into a topography version, no array is
modified, no elevation is fabricated. The added rows are drawn BLANK.

REWRITTEN 2026-09-07 for the symmetric footprint. Until then this script
computed N itself as round(shift / 10) FLOORED AT ZERO and reported what the
negatives "would have implied". Hannah's decision that day made the footprint
symmetric (rows removed where the island prograded), changed the per-domain
statistic to the median of PAIRED per-profile differences, and changed the
rounding to a 10 m threshold (trunc). Those rules are applied in ONE place,
HAT_footprint_1984.py, which writes `footprint_1984_by_domain.csv`; this script
READS that table rather than re-deriving N, so the two cannot disagree. The
rules, and the assumptions behind them, are in that script's docstring.

WHAT THIS ADDS TO THE FOOTPRINT SCRIPT
    * the grid drawn AS THE MODEL WOULD HOLD IT: dune rows on top, then the
      added rows, then the existing interior pushed down the page by N; where
      rows are removed, the existing rows 0..|N|-1 are hatched between the
      dune and the rows that survive. HAT_footprint_1984_grid.png draws the
      same footprint in the CURRENT frame (row 0 fixed); this one shows the
      stack. NC-12 is drawn at its measured position in both.
    * the cross-check against the independent easting-frame measurement
      (`0-elevation/2009-2014-1996-duneline/duneline_offset_by_domain.csv`:
      raw easting in the axis-aligned box, 1 m sampling, no extractor frame,
      no row 0). Same two geojsons, so agreement bounds the FRAME, not the
      lines. The same trunc rule is applied to it.

THE PLAN VIEW moved. `HAT_row_insert_plan.png` (add-only) was retired the same
day; the plan view of the symmetric footprint is
`figures/3-placement/seaward/HAT_footprint_1984_plan.png`, drawn by HAT_footprint_1984.py.
Drawing it twice under two names would be one figure with two provenances.

OUTPUTS (1-barrier3d-domains/1984-start/2-domain-reconstruction-1984/)
    HAT_row_insert_scope.txt          the report
    row_insert_scope_by_domain.csv    per domain: signed N, shift, the easting
                                      cross-check, rows now/after
    figures/3-placement/seaward/HAT_row_insert_grid.png             the stacked grid, both signs
    figures/3-placement/behind-road/HAT_row_insert_grid_behindroad.png  the same rows behind NC-12 (--anchor road)
    figures/4-fill/HAT_fill_copy_grid_island.png   ... filled by the copy rule (--anchor road --fill copy)

USAGE
    python HAT_footprint_1984.py          # first - writes the footprint table
    python HAT_report_row_insert_scope.py
```

Notes that were in the code:

```text
Colours: the RdBu pair every 1984 dune-line figure uses. Red = 1984 line
seaward = ground ADDED; blue = 1984 line landward = ground REMOVED.
```

```text
rows 0..ins-1 stay where they are; the block goes in at `ins`; the
rest of the interior is pushed down by N (add) or stays (remove -
the caller hatches the rows that go)
```

```text
the fill rule (HAT_fill_copy_scope.py): the N rows that follow
the insert point, copied cell by cell
```

```text
the replaced rows, outlined so the fill reads as part of the
interior and still shows where it is
```

```text
NC-12 at its measured position (seaward edge, 20 m). It does not
move: with rows added it is pushed down with the interior; with
rows removed it stays where the surviving rows put it.
dune anchor: NC-12 at its MEASURED position (unfloored; at GIS 85/86
that is seaward of row 0), pushed down with the interior. Road
anchor: NC-12 where the MODEL holds it - two straight rows at
int(setback/10), the floored setback - and the block sits behind it.
```

```text
Two road symbols, drawn the same way at every road domain
(2026-09-08, Hannah: the roadway labelling was not clear).
filled  NC-12 as the model holds it in v3: two rows at the
1984 setback, int(setback_new/10) from row 0
outline NC-12 as surveyed: the pavement rows of the 1984
alignment on this surface, today's setback
Where rows are removed in front of the road the strip is still
the v2 frame (the hatched rows are on it), so the model's road
is drawn on the v2 cells it ends up on: |N| rows further down.
```

```text
Legend wording (2026-09-08, Hannah): plain, academic, no working
vocabulary - "today's setback", "v3", "as surveyed" - on the figure.
The two road symbols are the same road at two setbacks: the one measured
on the 1996 surface, and the 1984 one the model is given.
```

```text
the table and the report are the same for both placements (N does
not change); only the grid is redrawn
```

<details><summary>Function notes (the original docstrings)</summary>

**`n_cells()`**

```text
trunc(shift / 10): a row only once a FULL cell of change is measured.
The rule HAT_footprint_1984.py applies; repeated here ONLY for the
easting-frame cross-check, which that script does not carry.
```

**`_community_bar()`**

```text
The communities as a bar along the BOTTOM of the strip, names below it and
village ticks above it, from HATTERAS_ANNOTATIONS - the same object every
other island figure uses. Below rather than above (2026-09-07): above, the
names collided with the +N labels and with each other (Tri-Village's
villages).
```

**`build_strip()`**

```text
One alongshore strip as RGBA, rows DOWN the page from the dune. The dune
stays put. Added rows go in between the dune and the existing interior,
which is pushed down by N. Removed rows are the existing rows 0..|N|-1,
left in place here and hatched by the caller: what the model would hold is
the interior starting at row |N|. Existing cells carry the project's
elevation classes; white is off the array.
```

</details>

### 3-placement/HAT_verify_road_placement_1984.py

Is NC-12 placed where the 1984 road and dune lines say it was? Checked in three frames that must agree.

From the script's original header:

```text
Is NC-12 placed where the 1984 road and dune lines say it was?  (Hannah,
2026-09-08: "help me ensure that the roadway is being placed correctly and
that the offset matches reality")

The check runs in three frames and they have to agree:

  MAP        the 1984 dune line and the 1984 NC-12 centreline on the 1 m
             lidar, in map metres, no Barrier3D processing. Two distances
             per domain: ALONG the extractor's profiles (raster rows, the
             frame the model indexes) and PERPENDICULAR (frame-free, nearest
             point on the road from samples along the dune line). They differ
             only by the obliquity of the profiles to the island.
  ROW 0      the same distance expressed as the model needs it: metres
             landward of interior row 0 (one cell behind the picked 1996
             crest). The 1984 dune line is a toe, ~19 m seaward of row 0,
             so the row-0 setback is the toe-to-road distance minus that
             per-profile feature term:
                 setback_new = (road - line84) - (row0 - line97)
                             = (road - row0) + (line97 - line84)      [identity]
             It is the SAME measurement; only the reference changes.
  MODEL      what v3 receives: RoadSetback_1984_dunestart.csv must equal
             setback_new_m; the road rows are int(setback/10) and +1; the
             truncation loses 0-10 m; the array must hold rows_before + N
             rows and the road must sit on it. With rows added the road as
             placed lies on the v2 cells N rows behind today's pavement
             (the block goes in behind it); with rows removed in front of
             the road it lies on the old pavement's first row or 0-2 rows
             seaward of it (`road_cells_offset`: the row count comes from
             the shift median, the setback from the median of the sum).
             bulldoze() overwrites the road rows with the road elevation
             every year, so which measured cells lie under the pavement
             does not change the run; the distance from the crest does.

OUTPUTS  2-domain-reconstruction-1984/3-placement/road_placement_check_1984.csv     per road domain
         2-domain-reconstruction-1984/3-placement/HAT_road_placement_check_1984.txt  the report
         figures/3-placement/behind-road/
             rows-{added,removed}/HAT_road_placement_check_GIS<N>.png
                 map (1 m lidar, the two
                 dune lines, NC-12, row 0, the model's road rows and the
                 block / removed rows in map space) beside the v3 grid
             HAT_road_placement_check_island.png   every road domain: the
                 map distance along profiles against perpendicular, and
                 the setback today / the 1984 one / the one the model holds

USAGE
    python HAT_verify_road_placement_1984.py                 # 85, 63, 49, 16
    python HAT_verify_road_placement_1984.py --domains 85,84
```

Notes that were in the code:

```text
the two medians do not add: the row count follows the shift median, the
setback the median of the per-profile sum
```

```text
the model's road AS PLACED, on the v2 cells it actually covers: v3 row r
is v2 row r in front of the seam and r + |N| behind it, so with rows
removed the two road rows can straddle the seam and sit apart on the map
```

```text
window: the seaward ~700 m of the box around the road (set before the
labels, which need the edges)
```

```text
the label sits over its arrow's middle unless that is near a panel
edge, where it hangs off the arrow's inner end instead of spilling out
```

<details><summary>Function notes (the original docstrings)</summary>

**`perpendicular()`**

```text
Nearest distance from samples along the 1984 dune line (inside the
domain box) to the 1984 NC-12 centreline, minus the half width. Frame-free.
```

**`_band()`**

```text
A cross-shore band over rows r0..r1 (relative to row 0, inclusive) on
every profile of the domain, as one polygon in map space.
```

</details>

### 3-placement/imagery-review/HAT_imagery_review_1984.py

The 1984 footprint against the aerial photographs: each changed domain in the 1984 and 1997 photos, on the lidar.

From the script's original header:

```text
The 1984 footprint against the aerial photographs: for every domain where the
footprint adds or removes rows (and a set of unchanged neighbours as controls),
the same window of the island in the 1984 and 1997 photographs, on the 1 m
lidar the rows are cut from, with the two digitized dune lines, NC-12, and BOTH
candidate placements of the rows drawn on each. Under the panels, the
photograph's brightness along the domain's 50 profiles, so the sand-to-
vegetation transitions can be read as a curve without anyone picking them.

WHAT IT DECIDES - NOTHING. It is a review aid (Hannah's colleagues, 2026-09-08:
"look at the aerial imagery for the domains where we add and remove rows ...
to see how the dunes and interior have actually changed"). The judgement is
made by eye and written into `imagery_review_1984.csv`; this script fills the
numbers, draws the figures, and leaves the verdict columns blank. Re-running it
keeps whatever verdicts are already in the sheet.

THE QUESTION THE FIGURES ARE BUILT TO ANSWER (Hannah, 2026-09-08): did the dune
field actually get narrower between 1984 and 1997, and if so from which side -
so that the rows the footprint adds or removes can be placed where the width
was actually lost or gained. v3 books every change behind NC-12 (rows added)
or directly in front of it (rows removed); the seaward alternative books it at
the dune. Both are drawn on every photograph so the reviewer can say which one
the photographs support, per domain:

    extra_width_was     seaward_of_crest | crest_to_road | behind_road | none
                        where the 1984 island was wider (rows added) or
                        narrower (rows removed) than in 1997
    dune_field_change   narrower | wider | same | unclear
    edge_moved          seaward | landward | both | none
    placement_ok        yes | no | unclear      does v3's placement match
    confidence          high | medium | low
    notes               free text

WHAT IS READ OFF THE PHOTOGRAPHS (decided with Hannah, 2026-09-08)
    * the SEAWARD VEGETATION LINE - bright sand to dark vegetation, the
      clearest edge on a greyscale photograph and roughly the dune toe, the
      feature the digitized lines trace (~19 m seaward of interior row 0).
      Preferred over the wet/dry line because it is less sensitive to the
      beach state on the day: 1984-09-19 is days after Hurricane Diana and
      1997-10-12 is a normal autumn beach.
    * the ROAD CENTRELINE as visible in each photograph, so the crest-to-road
      distance can be judged per year. Both NC-12 alignments are drawn (the
      1984 line is the 1978 export, deliberately; see the road-line README),
      and the reviewer should trust the pavement in the photograph over either.

IMAGERY (D:\Hatteras_GIS\Aerial, the USGS Henderson release, doi
10.5066/P1CXBCDW): georeferenced to the Dare County 2007 orthophotos in NC
State Plane feet, ~1 ft pixels, stated horizontal accuracy 1.2 m for both 1984
and 1997 (RSS of the 2007 control, the scan resolution and the fit). 1984 is
one mosaic in UTM 18N at 0.26 m; 1997 is 32 frames, tiled per domain here, no
mosaic built. Other years in the release (1978-2002) can be asked for with
--years; 1996 is there but Hannah judged it poor (2026-09-08), so the default
pair is 1984 and 1997, the year the second dune line was digitized from.
Anything under the 1.2 m accuracy is not evidence; a Barrier3D cell is 10 m.

FRAME. Everything is drawn in map coordinates (EPSG:3725, the 1 m tiles'
frame), so ALONGSHORE_FLIP does not apply. The model rows are placed on the map
exactly as HAT_verify_road_placement_1984.py places them: map x = interior_x -
row * 10 along each profile, where (interior_x, interior_y) is interior row 0
from RoadOffset_1984_profiles.csv for the 82 road domains, and from the same
cell_to_map chain (re-run through the extractor) for GIS 1-8, which that file
does not cover. On the first road domain met, the re-run is checked against
the CSV to the centimetre, so the two sources cannot silently disagree.

BRIGHTNESS STRIP. For each profile and each 10 m cell along it (from 250 m
seaward of row 0 to the landward edge of the window), the mean pixel value of
the 10 x 10 m block; the strip is the median over the 50 profiles with the
25-75 % band, per year, each year scaled to its own 2-98 % range over the
window so that two films of different exposure can share an axis. Sand is
bright, vegetation and water dark; a step down going landward is the
vegetation line. It is a reading aid, not a measurement: nothing is picked,
nothing is written from it.

CONTROLS. Unchanged domains that border a changed one (both sides of every
add/remove run), thinned to --controls evenly along the island. Same
photographs, same reach, so the eye has a local baseline for "no change".

OUTPUTS  2-domain-reconstruction-1984/3-placement/imagery-review/   (the evidence for the placement step, 2026-09-09)
    imagery_review_1984.csv           the review sheet: numbers filled, verdicts blank
    HAT_imagery_review_1984.txt       the report: sources, rules, domain list
    ../../figures/3-placement/imagery-review/rows-{added,removed}/ , unchanged/
        HAT_imagery_review_GIS<N>.png   one per reviewed domain
    figures/CAPTIONS.md               the caption, one section

USAGE
    python HAT_imagery_review_1984.py                      # 52 changed + 10 controls
    python HAT_imagery_review_1984.py --domains 85,63,84   # a pilot
    python HAT_imagery_review_1984.py --years 1984,1997,1996
    python HAT_imagery_review_1984.py --no-controls
    python HAT_imagery_review_1984.py --resume             # after an interrupted run
```

Notes that were in the code:

```text
Written by the window (HAT_imagery_review_gui.py) from the reviewer's clicks on
each photograph: three features per year, on the profile nearest the click -
toe   the seaward vegetation line (beach sand -> dune), what the dune lines trace
back  the landward edge of the dune band (dune sand -> flat vegetated interior)
road  the seaward edge of the pavement as it is in that year's photograph
- each as a position and as metres landward of interior row 0; the bands between
them per year; their change 1997 - 1984; and where the lost / gained width sat,
given N: the dune band, the strip from the back of the dune to the road, or the
remainder, behind the road. Measured by the reviewer, not by code. `band_suggests`
is DERIVED from the picks and labelled so; the verdict columns stay the reviewer's.
```

```text
The quick review (HAT_imagery_review_quick.py, 2026-09-10): two yes/no/unclear
answers per domain - is the 1984 road offset right, is N right - kept across re-runs
like the columns above.
```

```text
1984's frames sit under 1984_georef_TIF/Input: georeferenced, 1 ft, the
input to the mosaic, and the only thing that fills the mosaic's gaps
```

```text
the mosaic is read first; the frames only fill where it has no pixels
(the 1984 mosaic has gaps, e.g. the south of GIS 63)
```

```text
the drive is external and can disappear mid-session; the batch
script stops, the window carries on from its cache
```

```text
A coarse request (an island-wide map, 2026-09-23) reads the
frame decimated -- through its overviews where it has them --
to ~res/2, instead of pulling every 0.3 m pixel into memory.
A fine request (the domain windows) reads at full resolution
as before.
```

```text
rows removed in front of the road: v3 row r maps to v2 row r
before the seam and r + |N| after it
```

<details><summary>Function notes (the original docstrings)</summary>

**`_year_files()`**

```text
The georeferenced tifs for one year, preferring a finished mosaic.

The release folders are not uniform: 1984 has a mosaic plus the frames it
was built from (under 1984_georef_TIF/Input; /zip holds archives), 1996/1997
are frames only, later years are
mosaics under other names, and the 2025 draft folder holds every year
again as frames. One rule: a `<year>_full_aerial.tif` wins; otherwise the
frames named `<year>_MMDD_*.tif` outside any Input/zip/thumbnail folder.
```

**`placements()`**

```text
Row ranges (relative to interior row 0, v2 frame) of everything drawn.

v3 (the live placement): rows added go directly behind the model's road as
placed under the 1984 setback, rows removed come out directly in front of
today's pavement; no-road domains use the crest row. The seaward
alternative: rows added between the 1984 line and row 0 (drawn as the |N|
cells seaward of row 0), rows removed as rows 0..|N|-1.
```

**`domain_window()`**

```text
The window drawn for one domain: from a little seaward of the seaward-most
dune line to WINDOW_PAD_LAND_M behind the landward-most row drawn, clipped to
the domain box; the box's full alongshore extent plus 10 m.
```

**`_edge_distance()`**

```text
Distance to the nearest no-photograph pixel, in metres, on a coarse
overview of one file (computed once per file, cached).

The frames carry a dark FRINGE inside their black border (values ~40
on 1997_1012_040163d, not 0), so "pixel > 0" cannot tell photograph
from film edge. Instead every pixel of the tile takes the frame it lies
farthest inside, the ordinary seamline rule, and a fringe a few metres
wide can never win against a frame that has real photograph there.
```

**`read()`**

```text
(H, W, 3) uint8 in the map frame at `res`, 0 where no photograph.

Where files overlap, each pixel comes from the file it lies farthest
inside (see _edge_distance). The mosaic is one file among the others,
so its gaps are filled by the frames and its own pixels win elsewhere.
```

</details>

### 3-placement/imagery-review/HAT_imagery_review_gui.py

The aerial-imagery review of the 1984 footprint as a window you work in, with measuring tools.

From the script's original header:

```text
The aerial-imagery review of the 1984 footprint as a window you work in,
instead of 62 PNGs and a spreadsheet. Same data, same overlays, same sheet as
HAT_imagery_review_1984.py - this only changes how the judgement is entered,
and lets the reviewer MEASURE the one thing the verdict rests on.

WHAT THE WINDOW SHOWS
    Left: the 1984 and 1997 photographs of one domain side by side (a third
    panel, the 1 m lidar, on request), sharing one view so pan and zoom in
    either moves both (the matplotlib toolbar below them; the mouse wheel
    zooms about the pointer). Under them the brightness strip. Every overlay
    - the two digitized dune lines, NC-12, interior row 0, the road rows,
    v3's rows, the seaward alternative, the domain box, your picks - is a
    checkbox, so the photograph can be looked at bare and the lines brought
    back.
    BLINK MODE puts both photographs in ONE panel and flips between them on
    Space (or B), or on a timer. The eye catches movement between two frames
    of the same view far better than side by side, which matters for the
    one-cell offsets most of the footprint is made of.
    Right: the domain's numbers from the footprint table, the pick buttons,
    the verdict form (the six columns of imagery_review_1984.csv, with the
    vocabulary as drop-downs), Save, Prev / Next, and Summarize.

THE PICKS (measured by the reviewer, not by code)
    Three features per year, each one click on the photograph after its button:
        toe    the seaward vegetation line, beach sand -> dune vegetation: the
               feature the digitized dune lines trace
        back   the landward edge of the dune band, hummocky dune sand -> flat
               vegetated interior
        road   the seaward edge of the pavement AS IT IS IN THAT PHOTOGRAPH
    Each click is taken on the profile nearest it and stored as a position and
    as metres landward of interior row 0 (negative = seaward). From them, per
    year, the three bands the placement question is about:
        dune_band     back - toe        the dune field
        back_to_road  road - back       the strip between the dune and the road
        toe_to_road   road - toe        the whole crest-to-road space
    and once both years are picked, their change 1997 - 1984, the shift of each
    feature (positive where the 1984 feature lay seaward, the footprint's sign),
    and - taking N as given (Hannah, 2026-09-09) - where the lost or gained
    width sat:
        lost_dune_band     = -(d_dune_band)
        lost_back_to_road  = -(d_back_to_road)
        lost_behind_road   = N x 10 m - the two above      (the remainder)
    `band_suggests` names the largest share (dune_band / back_to_road /
    behind_road / none when N is 0) and is DERIVED; the verdict is yours.
    Nothing is snapped: the number is where you clicked, on the photograph as
    georeferenced (stated accuracy 1.2 m). A domain without a model road still
    takes the road pick if a road is visible; the bands that need it stay blank
    otherwise.

WHAT IT WRITES
    On Save (Ctrl+S, or "Save & next"): the six verdict columns and the
    measured columns of the domain on screen, plus reviewed_by / reviewed_at,
    into 2-domain-reconstruction-1984/3-placement/imagery-review/imagery_review_1984.csv. Nothing else in the sheet
    is touched, and HAT_imagery_review_1984.py keeps these columns when it
    re-runs. "Summarize" runs HAT_imagery_review_summary.py on the sheet as
    it stands.

KEYS   Right / Left  next / previous domain      Ctrl+S  save
       Space or B    flip the year in blink mode   Esc     cancel a pick
       1 2 3         pick toe / back / road on the year shown (blink mode)
       (keys are ignored while the cursor is in a text box)

PERFORMANCE
    A domain takes ~5-10 s to read from the drive the first time. The arrays
    are cached under ~/.cascade/imagery_review_cache/ (outside the repo, keyed
    on domain, year, window and resolution), and the next domain in the list
    is read in the background while you look at the current one. Delete the
    cache folder to force a re-read (e.g. after changing the merge rule in
    the batch script).

USAGE
    python HAT_imagery_review_gui.py                 # every domain in the sheet
    python HAT_imagery_review_gui.py --domains 85,63 # a subset
    python HAT_imagery_review_gui.py --years 1984,1997,1996
    python HAT_imagery_review_gui.py --smoke         # open, draw one, screenshot, close
```

Notes that were in the code:

```text
where the width sat, taking N as given: lost = 1984 - 1997 (positive where
the 1984 island was wider there); the remainder is behind the road
```

### 3-placement/imagery-review/HAT_imagery_review_quick.py

The quick review: one domain at a time, the 1984 road offset drawn raw and as modeled on the photos, two questions.

From the script's original header:

```text
The basic review (Hannah's advisor, 2026-09-10): one domain at a time, the
1984 and 1997 photographs side by side, the 1984 ROAD OFFSET drawn on them
raw and as modeled, and two questions. Nothing to pick, nothing to toggle.

    offset_ok   yes | no | unclear   does the modeled 1984 road offset look
                                     right against the photographs?
    rows_ok     yes | no | unclear   does N, the rows the footprint adds or
                                     removes, look right?
    notes       free text

THE OFFSETS (decided with Hannah, 2026-09-10, revised the same afternoon)
    on the photographs, one arrow per panel in the year's colour: that year's
             dune line to the 1984 NC-12 line, both from the geojsons,
             measured along the 50 profiles (the line's crossing nearest row
             0, the road's crossing landward of it, moved HALF_M seaward to
             the pavement edge), median over the profiles. BOTH years
             reference the 1984 roadway (Hannah, 2026-09-10: "none of the
             offset ... should be in reference to the 2004 road position");
             the 2004 line is drawn for context only. What you would measure
             by hand on the map, year by year, against the road as it was.
    table    for the record in the text: `setback_raw84_m` from the footprint
             table, the same 1984 distance but to the road MASK in whole
             cells; it runs ~19 m larger than the model's reference because
             the digitized line traces the toe, a cell or two seaward of
             interior row 0.
    modeled  the setback CASCADE receives for v3 (`setback_new_m`, row-0
             convention, no floor) cut to whole cells: the road's first row is
             int(setback / 10) from interior row 0, and the pavement is the two
             rows from there. Drawn on the MODEL panel only, as the dark rows
             and an arrow from row 0 in rows.
    today    for reference only, in the text: the setback measured on the 1996
             surface (`setback_v2_m`, signed) and what the model holds for it
             (`setback_model_now_m`, floored at 0 where negative).
    Both arrows are drawn on the profile nearest the middle of the domain, so
    they start on the line they belong to; their LENGTH is the domain median,
    which is the number printed. The ~1.2 m stated accuracy of the
    georeferencing and the 10 m cell are the yardsticks.

WHAT IS ON THE PHOTOGRAPHS (settled in interview, Hannah, 2026-09-10, after
the overlays had grown to nine and "made it more confusing")
    Two verdicts, two comparisons, nothing else:
    * OFFSET. The 1984 NC-12 line and, on the middle profile, ONE white tick:
      interior row 0 measured seaward from the pavement edge - on (a) as the
      model places it for 1984 (the setback after its cut to cells, the same
      number panel (c) shows), on (b) as the DEM has it. The reviewer asks
      whether the crest visible in that year's photograph sits on the tick.
    * ROWS N. Both digitized dune lines (1984 red, 1997 blue) on both
      photographs, and one bracket between them on the middle profile
      labelled with the median shift and the rows it became. The reviewer
      asks whether each line follows the vegetation edge of its year and
      whether the retreat looks like the number.
    Plus the domain box. Removed from the photographs: interior row 0 as a
    line, the implied 1984 crest, the dune-line-to-road arrows, the row-0-to-
    road arrows, the dune search window, the 2004 road line. Their numbers
    stay in the side panel; the model-side marks stay on panel (c).

THE MODEL PANEL (upper right; the photographs stack down the left, 1984 over
1997, and the legend sits lower right - Hannah, 2026-09-10)
    The PROCESSED domain, the model input as the model holds it: the
    straightened Barrier3D grid of dune-topo/v3 (two dune rows at berm + dune
    height, then the interior, in the elevation classes of the placement
    check), cross-shore rows across with the ocean on the right and alongshore
    cells up the page, south at the bottom, so it faces the same way as the
    photographs. Nothing is mapped back through the shear: the dune is a
    straight band because that is what the model gets. With the model-side
    measurements on it: interior row 0, the road rows at the 1984 setback
    (dark), the rows the footprint inserted (red outline) or the seam it left
    (blue dashes), and the modeled offset as an arrow in rows. (The outline
    of today's pavement rows was dropped 2026-09-10: it is not part of the
    model input and read as a second road.) If dune-topo/v3 is not on disk the CURRENT version is
    drawn instead and the panel says so.

WHAT IT WRITES
    On Save: offset_ok, rows_ok, notes, reviewed_by, reviewed_at for the
    domain on screen into 2-domain-reconstruction-1984/3-placement/imagery-review/
    imagery_review_1984.csv (the sheet HAT_imagery_review_1984.py writes; that
    script keeps these columns when it re-runs). Nothing else in the sheet is
    touched, and the six older verdict columns are left as they are.
    "Save figure" writes the view on screen at 200 dpi to
        figures/3-placement/imagery-review/{rows-added,rows-removed,unchanged}/
            HAT_imagery_review_quick_GIS<N>.png
    with its caption in figures/CAPTIONS.md (one section for the set). The
    figure carries a header (domain, island section, N, the 1984 setback and
    its row), split legends for the photographs and the model input, and a
    source line (USGS release and DOI, stated accuracy, which road line each
    year is measured against).
    "Summarize" (or --summary) tallies the sheet:
        HAT_imagery_review_quick.txt                 counts, the "no" domains
        figures/3-placement/imagery-review/island/HAT_imagery_review_quick.png
            (a) N per domain along the island, coloured by rows_ok
            (b) the raw and modeled 1984 offset per domain, with offset_ok

KEYS   Right / Left  next / previous domain    Ctrl+S  save
       Space or B    flip the year in blink mode
       (keys are ignored while the cursor is in a text box)

Photographs are read through the same cache as the full window
(~/.cascade/imagery_review_cache/), so a domain seen in either is instant in
the other.

USAGE
    python HAT_imagery_review_quick.py                  # every domain in the sheet
    python HAT_imagery_review_quick.py --domains 85,63  # a subset
    python HAT_imagery_review_quick.py --summary        # tally only, no window
    python HAT_imagery_review_quick.py --smoke          # open, draw one, screenshot, close
```

Notes that were in the code:

```text
a SQUARE data window (the four panel boxes are equal squares): the rows
shown set the side; the alongshore range is centred on the domain and
padded with blank where the side exceeds 50 cells
```

```text
photographs down the left (1984 over 1997), the model input upper right,
the legend lower right; in blink mode one photograph left, the model right.
Four equal cells; every panel is a square data window with equal aspect,
so the four boxes come out the same size and line up (Hannah, 2026-09-10).
The header has its own row above; titles sit at one pad; the strip
below holds the legend when the DEM panel takes the legend's cell.
```

```text
the settled set (interview 2026-09-10): the 1984 road, both dune
lines, the domain box, one crest tick per year, one shift bracket
```

```text
---- (d) the DEM, on request: the surface row 0 was picked on, with the
diagnostic layers that belong to it (row 0 line, search window)
```

<details><summary>Function notes (the original docstrings)</summary>

**`_square_window()`**

```text
The photographs' window made SQUARE about its centre, so that the panel
boxes are equal squares (Hannah, 2026-09-10). Wraps the batch script's
domain_window; the loader reads the tiles for this window and caches them
under it, so the padding is real photograph, not blank.
```

**`line_offsets()`**

```text
The offset between that year's dune line and that year's NC-12, from the
geojsons, along the 50 profiles: per profile the line's crossing nearest
interior row 0 and the road's crossing nearest landward of it, the road
centreline moved HALF_M seaward to the pavement's edge; median over the
profiles. No cells, no mask: what the photograph's overlays show.
Returns {photo year: dict(m, n, x_line, x_road, y)} with the middle
profile's crossings for the arrow. Both years are measured to the 1984
NC-12 line (ROAD_REF): the question is the 1984 roadway position, and the
1997 arrow shows how far the dune had moved from it by then.
```

**`search_window()`**

```text
The extractor's picked dune search window for one domain, (i0, i1) as
cross-shore cells of the straightened profile array, i1 exclusive - the
same cell frame as `row0` in the profile table. None if not picked.
```

**`draw_search_window()`**

```text
The search window on the map: per profile the cells i0..i1-1 relative
to that profile's row 0 (map x = interior_x - row x 10).
```

**`draw_crest_tick()`**

```text
ONE tick per photograph on the middle profile: interior row 0 measured
seaward from the 1984 pavement edge. On the 1984 photograph it is the row
the model places for 1984 (the setback after its cut to cells, the same
number panel (c) shows); on the 1997 photograph it is the DEM's own row 0.
The reviewer judges whether the crest visible in that photograph sits on
the tick. A dotted connector from the road says what the metres refer to.
```

**`draw_shift_bracket()`**

```text
One bracket between the two dune lines on the middle profile, labelled
with the median shift and the rows it became. Ends at each line's
crossing of that profile; the LENGTH drawn is the domain median, which
is the number printed, so the bracket and the label agree.
```

**`model_r_max()`**

```text
The landward-most interior row the model panel shows: at least what the
photographs show, and past the road and the footprint's rows.
```

**`draw_model()`**

```text
The PROCESSED domain, the model input as the model holds it: the
straightened grid, dune rows a straight band, cross-shore rows across with
the ocean on the right and alongshore cells up the page (south at the
bottom, as the photographs). Nothing is mapped back through the shear.
```

</details>

### 3-placement/imagery-review/HAT_imagery_review_summary.py

What the filled imagery review says: verdicts along the island, and the toe picks against the measured shift.

From the script's original header:

```text
What the imagery review says, once the sheet is filled: the verdicts along the
island against the rows the footprint changes, and the reviewer's toe picks
against the shift the footprint measured from the digitized lines.

Reads 2-domain-reconstruction-1984/3-placement/imagery-review/imagery_review_1984.csv (written by
HAT_imagery_review_1984.py, filled by hand or through HAT_imagery_review_gui.py)
and footprint_1984_by_domain.csv. Writes nothing model-facing. Runs on a
half-filled or empty sheet and says so: unjudged domains are drawn hollow and
counted as such.

FIGURE  figures/3-placement/imagery-review/island/HAT_imagery_review_summary.png
    (a) every reviewed domain along the island: the signed rows N as a bar,
        coloured by the verdict `extra_width_was` (where the 1984 island was
        wider / narrower than in 1997): behind the road (v3's placement),
        between the crest and the road, seaward of the crest, or no real
        change; hollow where not yet judged. Controls (N = 0) as markers,
        filled where their `dune_field_change` says "same". A cross marks a
        domain whose `placement_ok` is "no".
    (b) THE PLACEMENT QUESTION, from the picks: per measured domain, where the
        1984 width sat, as stacked bars of N x 10 m - the part lost or gained
        in the dune band (toe to back of dune), in the strip from the back of
        the dune to the road, and the remainder behind the road - with the
        island-wide medians of the three shares by sign printed in the panel.
        The rule reads off the medians.
    (c) which feature moved: the shift of the toe, the back of the dune and
        the road edge (1997 - 1984, + = the 1984 feature seaward) against
        N x 10 m, per measured domain.

REPORT  2-domain-reconstruction-1984/3-placement/imagery-review/HAT_imagery_review_summary.txt
    counts by verdict for the changed domains (adds and removals apart),
    placement_ok, confidence, the controls, the domains where v3's placement
    is contradicted, and the agreement of the picks with N.

USAGE
    python HAT_imagery_review_summary.py
    python HAT_imagery_review_summary.py --sheet <other.csv> --out-dir <dir>   # a copy, for testing
```

### 4-fill/HAT_fill_copy_scope.py

What the added rows would contain under the copy fill, behind the road: scoped, drawn and audited.

From the script's original header:

```text
What the added rows would CONTAIN under the copy fill, for the behind-the-road
placement of the 1984 footprint - scoped, drawn and audited, no version written.

THE FILL (Hannah's advisor; decided with Hannah 2026-09-07)
    The block of N rows goes in directly behind the model's two roadway rows
    AS PLACED under the 1984 setback (insert_row_behind_road in
    footprint_1984_by_domain.csv = int(setback_new/10) + 2; 2026-09-08). It is filled
    with a DIRECT COPY of the N interior rows immediately landward of the
    insert point - rows r..r+N-1 of the existing interior, in order, cell by
    cell across the 50 alongshore columns - so the block fabricates no value
    and reads as the backbarrier it stands beside. The window follows N.

    * seams: the seaward junction is continuous by construction (the block's
      first row is the row that used to follow the road). The only seam is at
      the LANDWARD end, where the block's last row (a copy of r+N-1) meets the
      original row r; `seam_jump_m` is that step, alongshore mean of |dz|.
    * no-road domains (GIS 2-5): the footprint puts the block BEHIND THE
      CREST ROW (anchor "crest", insert_row = crest + 1, 2026-09-08), so the
      rule is the same as behind the road - the window is the N rows that
      follow the insert point, the crest stays at the front, nothing is
      duplicated. (Until 2026-09-08 the block went in at row 0 and the window
      skipped the crest, which left the crest row stranded between block and
      source.)
    * outliers: GIS 82/83 windows hold 10-11 m cells (Rodanthe structures in
      the DEM). Copied as measured - they are already in the interior - and
      flagged.
    * water: no window holds a cell at or below MHW, so no rule is needed and
      none is applied. If a future footprint changes that, `window_water_frac`
      is the column to watch.
    * removals need no fill. The |N| rows directly seaward of today's roadway
      rows are deleted (insert_row = int(setback_v2/10) - |N|; Hannah,
      2026-09-08: the rows come out of the interior in front of the road);
      the after-panel is read from the built version and checked.

WHAT IS ASSUMED
    * the lost 1984 ground is booked behind the road, so the backbarrier next
      to the road is the analogue for it (what was lost was ocean-side beach
      and dune);
    * the 1996/2009 backbarrier behind the road stands for 1984 backbarrier
      (no accretion or subsidence in between);
    * alongshore structure in the window appears twice cross-shore;
    * the dune array, row 0, the road rows and the model's setback are as in
      the behind-road placement - untouched.

OUTPUTS  2-domain-reconstruction-1984/
    fill_copy_by_domain.csv          per domain: N, insert row, source rows,
                                     window stats, seam jump, flags
    HAT_fill_copy_scope.txt          the report
    figures/4-fill/rows-{added,removed}/HAT_fill_copy_grid_GIS<...>.png
                                     for the example domains: the near-road
                                     interior before and after, in elevation classes
    figures/4-fill/rows-{added,removed}/HAT_fill_copy_method_GIS<...>.png
                                     the method in three stages, model frame: the
                                     domain as extracted (dune rows + interior),
                                     the N rows inserted blank behind NC-12, the
                                     copy fill with the source window and the copy
                                     drawn as an arrow

USAGE
    python HAT_fill_copy_scope.py                    # examples 80, 85, 5, 49
    python HAT_fill_copy_scope.py --domains 80,73
```

Notes that were in the code:

```text
one rule for every anchor: the window is the N rows that follow the
insert point. For a no-road domain the footprint already put that
point behind the crest row.
```

```text
If the version has been built, the after-panel IS that version's
array - read it and check it is what the rule says.
```

```text
(a) NC-12 as the model holds it today, two rows at int(setback/10);
(b) the same rows outlined as the measured pavement, and the model
road where v3's 1984 setback puts it - N rows inland, on the block
```

```text
(b), (c): the model road at its 1984 setback - the block sits
directly behind it; the measured pavement rows outlined
```

```text
the source window, in the after frame, and the copy drawn as an arrow.
The source starts at src0, which is the insert point for a road
domain but one row behind the crest for a no-road domain (GIS 2-5)
- so it is src0 + n after the insert, not ins + n.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_fig_domain_removal()`**

```text
Before and after for a removal domain: the |N| rows directly seaward of
today's roadway rows go; the road's cells and everything behind are kept.
```

**`_fig_method_removal()`**

```text
The three stages for a removal domain, model frame: as extracted, the
|N| rows identified (in front of NC-12), the rows removed and the road at
its 1984 setback.
```

**`fig_method()`**

```text
The method in three stages, one domain, all in the MODEL's frame (dune rows
on top, then the interior): (a) the Barrier3D domain as extracted, (b) the
N rows identified by the footprint inserted BLANK behind the roadway rows,
(c) the same rows filled by copying the N interior rows that follow them.
```

</details>

### 4-fill/HAT_insert_seaward_rows.py

Move the island seaward per domain by prepending interior rows, as a new dune-topo layer.

From the script's original header:

```text
Move the island seaward, per domain, by prepending fabricated interior rows.

WHY THIS EXISTS
The 1984-start topography is a 1996 ALACE beach and foredune grafted onto a
2009 backdune. The NC-12 road lines are 1984-vintage. Between 1984 and 1996 the
island migrated landward, so at the worst domains the 1984 roadbed now sits
SEAWARD of the 1996 dune crest and the measured setback comes out negative:

    GIS 85   road seaward cell 13, interior row 0 = source cell 14
             setback (13 - 14) * 10 = -10 m  ->  floored to 0 in the model CSV

A zero setback is not inert. roadway_manager.road_relocation_checks does
`road_setback += dune_migrated` every year, so the first year with any erosion
drives it negative and relocates the road. In the shipped 1984-2004 run GIS 85
relocates in YEAR 1 and twice more after that; GIS 86 the same; the relocation
year across GIS 10/11/84/85/86 is monotone in the input setback, not in the
physics. That is the artefact this script exists to remove.

WHAT IT DOES
For each selected domain it prepends N rows to the SEAWARD end of the saved
interior array, which moves interior row 0 N cells seaward and turns the road's
setback from (road - row0) into (road - row0) + N*10.

N comes from a MEASUREMENT, not a target: the per-domain cross-shore offset
between the digitized 1984 dune line and interior row 0, measured in the
extractor's own frame by HAT_measure_duneline_shift.py. Island-wide that median
is +18.9 m; at GIS 85 it is +65.9 m with a p10-p90 of +62 to +70.

The control that makes this credible is in the 2004 pair: measured the same way,
the 2004 dune line and 2004-start row 0 agree to +1.4 / -4.1 / -3.6 m on the
stable mid-island (GIS 40/50/60). So "digitized line vs DEM row 0" carries no
material feature offset, and the 1984 number is a date difference rather than a
definitional one.

WHAT IS FABRICATED, AND SAY SO
No survey covers land that was gone by 1996. The N new rows are invented. The
fill rule is explicit and recorded in the manifest; `backdune` (the default)
lays them flat at the median of interior rows 1-3, i.e. a backdune platform, NOT
at row 0's elevation -- at GIS 85 row 0 IS the 4.82 m crest, and copying it
would build a 70 m plateau at crest height.

Note what stays behind: the old row 0 becomes an interior ridge N cells inside
the island. That is a relict foredune, which is a reasonable thing for a
migrating barrier to have, but it is a consequence rather than a choice.

TWO VARIANTS, BECAUSE THE BAY SIDE IS NOT FREE
brie_coupler.offset_shoreline sets x_s per domain and the interior extends
LANDWARD from the dune, so prepending rows does not push the shoreline seaward
-- it pushes the bay edge further into the sound.

    pad         prepend N. Island gets N cells wider. Mean interior height and
                InteriorWidth_AvgTS both rise, which feeds overwash flux and the
                relocation room test.
    translate   prepend N and retire the N landward-most LAND rows per column to
                the water sentinel. Per-column land width is preserved. This is
                the barrier-migration reading: in 1984 the bay edge was also N
                cells seaward.

Neither is free and they are not equivalent; run both and compare before
picking. `none` writes the setback CSV alone and leaves the topography be.

OUTPUT
A NEW dune-topo version beside the source. v1 is never written to.

    1984-start/dune-topo/<DST_VERSION>/
        topography/   domain_<N>_topography.npy    (dam)
        dunes/        domain_<N>_dune.npy          copied unchanged
        RUN_MANIFEST.txt
        HAT_seaward_row_insert_audit.csv
        RoadSetback_1984_dunestart.csv             matched to this topography

USAGE
    python HAT_insert_seaward_rows.py                       # measured N, pad
    python HAT_insert_seaward_rows.py --variant translate
    python HAT_insert_seaward_rows.py --n-rule minimum --domains 85,86
```

Notes that were in the code:

```text
TWO measurements of the same quantity, and they disagree by a factor of ~3
at the domains this work is about (GIS 85: 65.9 m vs 19.5 m). Neither is
the truth; --shift-source names which one a build used and the manifest
records it, so no output can be read without knowing.
```

```text
THE ONE TO USE (2026-09-02). 1984 line minus 1997 line, same feature at
both ends, so the DEFINITIONAL offset between a digitized line and interior
row 0 cancels and what is left is pure date. Measured, not assumed: the
1997 line sits +16.2 m seaward of row 0 island-wide (IQR +12.8 to +21.0),
and subtracting that leaves an island-wide date term of +0.8 m -- i.e. the
ORIGINAL "duneline" source below was about 85% feature offset.
```

```text
Superseded. Kept so the earlier arms stay reproducible and so the size of
the correction stays visible rather than being quietly absorbed.
```

```text
"dsas" REMOVED 2026-09-03. It pointed at a file that had already been
moved into superseded/, so the option failed at runtime; and the
estimate measures the SHORELINE, not the dune line, understating dune
retreat by ~2.3x (GIS 85: 19.5 m against a measured 58.2 m). A route to
a known-wrong number is not worth keeping reachable.
```

```text
LAND is > 0 m MHW, not "> the water sentinel". At these domains the sound side
is MEASURED marsh and shallow bay sitting between -3 m and 0 m, so a sentinel
test calls the bay land and reports GIS 85 as 160 rows wide when its island is
37. 0 m MHW is the threshold bulldoze() itself uses (drown_threshold = 0) and
the one HAT_road_domain_views draws.
```

```text
MATCHED BACKDUNE: the existing near-dune profile copied in front of
itself, per column, so the 1984 block reproduces today's cross-shore
FORM at the 1984 position. Every value is a real measured cell of
THIS domain, shifted N cells seaward; none is a measurement at the
coordinates it lands on, so cells_from_dem is reported as 0.

matched-crest    block row k = interior row k      (k = 0 .. N-1)
Panel (b) of the 2026-09-03 fill figure, exactly:
row 0 IS the 1996 crest, so the crest appears at
the new seaward edge AND at its measured position
N cells landward. Two interior ridges by design.
matched-nocrest  block row k = interior row k + 1  (k = 0 .. N-1)
The same copy starting one row landward, so the
crest is skipped and the block is backdune only.

Built as separate rules (2026-09-04) rather than as one rule with a
switch because they are two ARMS of the fill comparison, and an arm
should be nameable from the manifest's fill_rule alone.
```

```text
Scope by the MEASUREMENT, not by where the road happens to be. Any
domain whose 1984 dune line sits >= half a cell seaward of interior
row 0 is missing 1984 land from the 1996 survey, and that is true of
domains with no NC-12 in them as well. Selecting the two relocation
blocks made "unchanged" mean two different things along the island:
30 domains had N >= 1 and were passed over by a scope decision rather
than by a measurement. This removes that discontinuity.
```

```text
`none` still credits the setback with N: it is the variant that says
"the dune-to-road distance was wrong, the island was not", so the
correction lands entirely in the CSV and the arrays are passed through.
```

```text
FLOOR THE REAL VALUE AT THE BACKDUNE PLATFORM, do not simply
take it. The measured cells here are the 1996 surface, and at
an eroding domain 1996 is the LATER and LOWER one -- so its
elevation is a lower bound on 1984's, not an estimate of it.
Taking it raw put a 0.57 m beach cell two rows inside the 1984
interior at GIS 85, which is 1996 beach standing where 1984 had
dry backdune. The floor keeps the real dune flank (1.98, 3.44)
and the road's own cell, and declines to import the beach.
```

```text
KEEP EVERY DRY MEASUREMENT, INVENT ONE NUMBER.

`measured` above justifies itself with "1996 is a lower bound
on 1984" and then RAISES 44% of the block above the DEM to
the platform, which is not a lower-bound operation - it
asserts the ground was at least backdune height. This rule
drops that second step. A dry cell is kept as measured; only
the cells with no usable measurement get a value, and that
value is the median of the block's OWN dry cells rather than
a statistic imported from interior rows 1-3.

One guard, one constant, and it never overrides a measurement
upward. At GIS 85 it takes the measured share from 47% to
91%: only 28 of 300 cells are at or below MHW, 25 of them in
the seaward-most row. NC-12 sits at rows 4-5, which are 100%
dry, so the constant does not reach the road at all.

WHAT IT ADMITS, deliberately: the 1996 beach ramp. Per-row
medians at GIS 85 run -0.00, 0.70, 1.20, 1.84, 3.17, 4.96, so
row 1 sits BELOW whatever fills row 0 and the block carries a
dip two cells inside the island. That is what the platform
floor existed to remove. It is admitted here because it is
what the measurement says, and the alternative asserts more
than the data supports.
```

```text
No dry cell anywhere in the block: nothing to take a median
of, so the backdune platform stands as the fallback and
n_real stays 0. Does not occur at any of the 38 domains, but
a silent nan-median would be worse than a stated fallback.
```

```text
Audit EVERY domain that was processed, not only the ones carrying a
road. Under --domains measured the scope is the whole island, and a
domain with no NC-12 still has rows inserted and still has to be
accountable for them.
```

<details><summary>Function notes (the original docstrings)</summary>

**`real_block()`**

```text
The DEM cells that ALREADY EXIST seaward of interior row 0, as (n, along).

build_interior fills topo[:, i] from prof_arr[i, row0[i]:], so the cells the
insert is about to cover are prof_arr[i, row0[i]-n : row0[i]] -- per profile,
because row 0 is a per-profile cut, not a horizontal one.

These are real measurements, just excluded from the interior for being
seaward of the dune pick. At GIS 85 they include cell 13 at 3.44 m, which is
the road's own seaward cell -- the cell whose elevation decides what
bulldoze() is scraping. Fabricating over it when the DEM has it measured is
a loss for nothing.

Returned in dam, NaN where the source cell is off the array. Cells at or
below water are left as NaN too: they are 1996 beach, and in 1984 that
ground was dry island. Those are the ones that genuinely have to be invented.
```

**`lower_old_crest()`**

```text
Shave the DEM's dune ridge down to the backdune platform, per column.

THE PROBLEM THIS SOLVES. Prepending N rows puts Barrier3D's dune at the 1984
dune line, but the DEM's own crest is still standing in the interior N cells
landward -- so the model starts with TWO dunes. That is the same sand counted
twice: the ridge is at the later position *because* the dune migrated there
by 1996, so in 1984 it had not formed yet. At GIS 85 it is 4.82 m, the
tallest thing in the domain, sitting on the road's own cell and shielding it
from landward.

THE COST, STATED PLAINLY. This discards a real measurement. The defence is
that it is a measurement of the wrong YEAR: the whole operation is de-aging
the surface by 12-25 years, and a 1996 crest is not a 1984 initial condition
just because it is real. The opposite choice -- keeping it -- is equally
defensible and is what --no-lower-old-crest gives you. Build both.

Walks landward from the seaward edge capping at the platform, and stops at
the first cell already at or below it, so it shaves the ridge and nothing
else. `max_reach` bounds the walk past the insert: a column whose profile
never drops back to the platform would otherwise be flattened across the
whole island, and that is a failure worth hearing about rather than
absorbing.
```

**`retire_landward_rows()`**

```text
Drown the N landward-most LAND cells of each column, into the local bay.

Per column, not per row: the island's landward edge is not a straight line,
which is the whole reason the road drown test looks at flanking rows rather
than a single width.

The retired cells take the elevation of the bay immediately landward of them
in the SAME column, not the -3 m sentinel. Stamping the sentinel would dig a
30 m trench along the sound edge of a domain whose real back-barrier is
measured marsh a few decimetres below MHW, and Barrier3D would read that as
the island having calved rather than migrated.

USES BARRIER3D'S OWN WIDTH DEFINITION, not a count of dry cells.

THIS WAS A BUG AND IT MADE `translate` BEHAVE LIKE `pad` (fixed 2026-09-02).
The first version counted every cell above 0 m MHW in a column and drowned
the landward-most n of them. Barrier3D's FindWidths (barrier3d.py:29) does
something else: it walks from row 0 and STOPS AT THE FIRST cell <= SL, so
anything beyond an interior water gap is not island at all.

On GIS 85 the two disagree in all 50 columns -- median 37.5 cells against 44,
and in the worst column 25 against 52, because the profile dips to -0.02 m at
row 26 and everything past it is sound-side marsh at 0.03-0.19 m. So the
cells being drowned were out in that marsh, which the model was never
counting; the retirement removed nothing while the prepended rows still
added. Measured effect: t=0 island width rose 263->309 m (D84), 361->402
(D85), 286->301 (D86) when it should not have moved at all.
```

</details>

### 4-fill/HAT_plot_fill_options.py

What can the added interior rows be made of? Four candidates and a control, on one domain.

From the script's original header:

```text
What can the added interior rows be MADE of? Four candidates and a control,
one domain.

THE PROBLEM
    The rows added behind the dune stand where land existed in 1984 and had
    eroded away by 1996. No survey covers that ground. The DEM does have cells
    at those coordinates, but they are the 1996 surface -- a later, lower
    landform at the same place. So every option below is a different answer to
    "what was here in 1984", and none of them is a measurement of it.

FOUR CANDIDATES AND A CONTROL
    Each is named by its RULE, in the legend, on the bar axes and in the
    caption -- the b/c/d/e tags this figure carried until 2026-09-10 read as
    panel letters beside the house style's own (a)/(b)/(c), which is a trap
    worth closing. The grid figure names the same rules the same way.

    flat backdune     flat, at the median of interior rows 1-3
    matched backdune  today's near-dune PROFILE copied to the 1984 position
                     -- "the 1984 backdune looked like the present one, just
                     further seaward". Real cells, so it carries alongshore
                     texture the flat fill cannot.
    measured + floor  the real DEM cell where it is dry land, floored at the
                     backdune platform. THIS IS THE SHIPPED RULE -- v4 at the
                     ten block domains, v5 island-wide. Same rule, new scope.
    measured + median  every dry cell kept AS MEASURED; only the cells at or
                     below MHW filled, with the median of the block's own dry
                     cells. `--fill median`. One guard, one constant, and no
                     measurement is ever raised -- 91% measured at GIS 85
                     against the floored rule's 47%.
    raw DEM          the 1996 cells as they are, no floor, no dry-land test.
                     A CONTROL, NOT A CANDIDATE -- drawn to show what the two
                     guards actually reject, in numbers.

    Dropped 2026-09-03: `taper`, a linear platform-to-row-0 ramp. Fully
    invented, and it anchored on row 0 -- which at GIS 85 IS the mis-picked
    1996 crest, so it inherited a known-bad endpoint. `--fill taper` still
    exists in HAT_insert_seaward_rows.py: removing a build capability is a
    different decision from removing a figure panel.

    `matched backdune` is NOT a --fill choice in HAT_insert_seaward_rows.py.
    It is drawn as a candidate; building it needs a new fill rule.

    Still undrawn: an alongshore analogue from a neighbouring domain, and a
    mass-conservative reconstruction. Neither is a one-line variant of the
    others, and the second is the only one that would be DERIVED rather than
    asserted.

WHAT TO LOOK FOR
    Where NC-12 lands. The road is a fixed 2-cell block at a fixed setback, so
    the only thing that changes between options is the ground under it.
    `measured + floor` puts it on 3.17 / 4.96 m because those cells are the
    1996 DUNE FACE; the flat and matched backdunes put it on backdune, where a
    road behind a dune belongs.

USAGE
    python HAT_plot_fill_options.py [--domain 85]
```

Notes that were in the code:

```text
INS_V supplies N and the post-insert setback ONLY; the figure draws
blocks it builds itself. Repointed v4 -> v5 on 2026-09-03: N is
identical at all ten block domains (verified), and v5 is the version
taken forward, so v4 no longer has to exist for this figure to build.
DELETED 2026-09-07 with every layer (only unmodified topography is kept);
the literal is kept as the name of what this drew. require_version() in
main() says so before any array is opened.
```

```text
TAGGED b/c/d/e TO MATCH THE GRID FIGURE'S PANEL LETTERS, where (a) is the
reference. The two figures are read side by side and single letters that
meant different things in each was a trap worth closing.

THE "% FROM THE DEM" IS CARRIED HERE, beside the block it describes.
It used to be recomputed further down by an `if tag == "A"` chain, which
fell through to the wrong branch the moment these tags were renamed and
silently reported 47% for every option. A number that describes a block
belongs with the block.

`matched backdune` counts as 0%: its cells ARE real measurements, but of
a different place, copied. Measured-ness here means "measured AT THESE
COORDINATES", which is the only sense that bears on whether the fill is
invented.
```

```text
(key, legend label, axis label, block, colour, % taken from the DEM)

Five house colours, one per rule; the same colour names the rule in
all three panels. "measured + floor" is the rule the layers were
built with, so it carries C["ACCENT"], the change under test.
```

```text
NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
```

### 4-fill/HAT_plot_fill_options_grid.py

The candidate interior fills as Barrier3D domain views: the grid the model would be handed under each.

From the script's original header:

```text
The candidate interior fills as BARRIER3D DOMAIN VIEWS -- the grid the model
would be handed under each, not a profile line through it.

A median profile hides the thing that decides overwash: whether the added rows
are uniform across the domain or carry alongshore structure. The flat fill has
none by construction; the matched backdune carries today's backdune texture;
the measured fill inherits the 1996 dune face and is as variable as that dune
was. The grid shows that directly.

    reference          v3, no rows added. NC-12 at its MEASURED offset.
    flat backdune      median of interior rows 1-3, per column
    matched backdune   today's near-dune profile copied to the 1984 position
    measured + floor   the real DEM cell where dry, floored at the platform.
                       THE SHIPPED RULE -- v4 at the ten block domains, v5
                       island-wide.
    measured + median  every dry cell kept as measured; only the cells at or
                       below MHW filled, with the median of the block's own
                       dry cells. `--fill median`.
    raw DEM, no floor  A CONTROL, NOT A CANDIDATE.

`taper` was drawn here until 2026-09-03 and is gone: it was fully invented and
it anchored on row 0, which at GIS 85 IS the mis-picked 1996 crest, so it
inherited a known-bad endpoint. `--fill taper` still exists in
HAT_insert_seaward_rows.py - removing a build capability is a different
decision from removing a figure panel.

`matched backdune` is NOT a `--fill` choice in HAT_insert_seaward_rows.py. It
is drawn here as a candidate; building it would need a new fill rule.

THE ROAD IS DRAWN AT ITS MEASURED OFFSET, NOT THE FLOORED ONE
    int() truncates toward zero, so a measured -15 m gives road_start = -1: the
    road lands SEAWARD of interior row 0, drawn overlapping the dune strip and
    hatched. That is not a plotting artefact, it is the failure --
    roadway_manager would evaluate xyz_interior_grid[-1:1, :], valid Python
    indexing from the LANDWARD end, and bulldoze the sound-side marsh. Flooring
    it in the figure is what made the problem invisible in earlier versions.

USAGE
    python HAT_plot_fill_options_grid.py [--domain 85] [--rows 26]
```

Notes that were in the code:

```text
INS_V supplies N and the post-insert setback ONLY; the figure draws
blocks it builds itself. Repointed v4 -> v5 on 2026-09-03: N is
identical at all ten block domains (verified), and v5 is the version
taken forward, so v4 no longer has to exist for this figure to build.
DELETED 2026-09-07 with every layer (only unmodified topography is kept);
the literal is kept as the name of what this drew. require_version() in
main() says so before any array is opened.
```

```text
THE ROAD ELEVATION CASCADE ACTUALLY USES, m MHW.

hatteras_site_config.HATTERAS_ROAD_ELEVATION_FILE resolves to THIS file.
There is a second one, road_offset/dunestart_offset/measured/1984/
RoadElevation_1984_dunestart.csv, and it must NOT be used here: it samples
along the 1984 alignment, which at the relocated domains (GIS 9-15, 84-87)
now lies UNDER the foredune, so it returns dune rather than roadbed. At
GIS 85 the two read 0.807 m and 1.833 m - a metre of difference that is
entirely the abandoned corridor being buried. See the long note beside
HATTERAS_ROAD_ELEVATION_FILE in hatteras_site_config.py.
```

```text
MATCHED BACKDUNE. The existing first n interior rows copied in front
of themselves, so the near-dune PROFILE is reproduced at the 1984
position: "the 1984 backdune looked like the present one, just
further seaward". Every value is a measured cell, so unlike the flat
fill it carries real alongshore texture - but it is counted as 0%
measured because those cells are measurements of the WRONG PLACE,
copied, not of the ground being filled.
```

```text
KEEP DRY, FILL WET WITH THE BLOCK'S OWN MEDIAN. `--fill median`.
Same dry-land test as the shipped rule, but no second step: a
measurement is never raised, and the one invented number comes from
the ground being filled rather than from interior rows 1-3.
```

```text
NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
```

<details><summary>Function notes (the original docstrings)</summary>

**`draw()`**

```text
`road_ele` paints the road block at a SINGLE elevation, which is what
the model holds: bulldoze() does `np.zeros(...) + road_ele`, so after the
first pass every road cell carries the same value regardless of the ground
that was there. Left as None for the reference panel, whose setback is
negative and whose road therefore falls outside the interior array.

The per-panel note box of quantitative results was removed 2026-09-10: the
house rule is that statistics belong in the caption, and main() now builds
one caption line per panel from the same numbers.
```

</details>

### 4-fill/HAT_plot_insert_explainer.py

How the 1984 seaward-row insert is built, on one cross-shore line at one domain, from the real arrays.

From the script's original header:

```text
How the 1984 seaward-row insert is built, on one cross-shore line at one domain
(GIS 85 by default), drawn from the real arrays.

    (a) the survey along that line, cell by cell, shaded by survey year, with
        the crest pick and the measured 1984 road position
    (b) what the extraction keeps: the crest cell becomes Barrier3D's dune
        rows, everything landward of it is the interior (row 0 first), and
        everything seaward is dropped
    (c) the insert: N rows are added at the dropped coordinates, the dune rows
        move in front of them, the old crest stays inside, and the road does
        not move -- its setback grows by 10 N
    (d) what each version writes into the N cells, over the survey values that
        are there

Alongshore MEDIANS of the domain's 50 profiles; the cross-shore axis is metres
from v2's interior row 0 (negative = seaward), the frame every version shares.

USAGE
    python HAT_plot_insert_explainer.py [--domain 85]
```

Notes that were in the code:

```text
DELETED 2026-09-07 with every layer (only unmodified topography is kept);
the literal is kept as the name of what this drew. require_version() in
main() says so before any array is opened.
```

```text
NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
```

### 4-fill/HAT_plot_insert_explainer_grid.py

The 1984 seaward-row insert explained as Barrier3D plan-view grids, on one domain.

From the script's original header:

```text
The 1984 seaward-row insert explained as Barrier3D PLAN-VIEW GRIDS (rows
cross-shore, columns alongshore, colour = elevation), on one domain.

    top row     (a) the survey around the pick, in the extractor's frame,
                    referenced to v2's interior row 0; 2009-survey cells
                    hatched; the crest pick and the measured road marked
                (b) the v2 domain as the model sees it: two dune rows then
                    the interior; nothing seaward of the crest exists
                (c) the inserted domain (v5 as the example): dune rows, the N
                    added rows outlined, the old crest a ridge inside, and
                    the road on rows 4-5 of the new interior
    bottom row  the first rows of v4-v8 side by side: same dune rows, same
                interior from old row 0 on, only the N added rows differ

Everything is drawn from the arrays on disk (v2, v4-v8 topography and dune
files, the survey-year clip) - no run output.

USAGE
    python HAT_plot_insert_explainer_grid.py [--domain 85]
```

Notes that were in the code:

```text
DELETED 2026-09-07 with every layer (only unmodified topography is kept);
the literal is kept as the name of what this drew. require_version() in
main() says so before any array is opened.
(folder on disk, the fill rule it holds, the rule as a panel title). The
version token names the LAYER, which is working vocabulary: it is carried
here and in the caption's mapping, never on the canvas (2026-09-10).
```

```text
NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
```

### 5-build/HAT_build_footprint_version.py

Build a 1984-start dune-topo version from the footprint table: rows added or removed behind the road.

From the script's original header:

```text
Build a 1984-start dune-topo VERSION from the footprint table: the symmetric
1984 footprint (rows added where the 1984 dune line lay seaward of the 1997
line, removed where it lay landward), placed BEHIND THE ROAD (or behind the
crest row where there is no model road) and filled by the COPY rule.

WHAT IS WRITTEN (dune-topo/<dst>/)
    topography/domain_<N>_topography.npy   v2's array with the block inserted
                                           (a copy of the N rows that follow the
                                           insert point) or the rows removed
    topography/domain_<N>_nodata.npy       the same row operation on the mask
    dunes/domain_<N>_dune.npy              copied unchanged: the dune stays put
    RoadSetback_1984_dunestart.csv         the 1984 setbacks: setback_new_m from
                                           the footprint table, the road measured
                                           against the 1984 dune line in the
                                           model's row-0 convention, (road - row 0)
                                           + shift per profile. Never negative, so
                                           NO FLOOR. With the rows behind the road
                                           this moves the model's road N rows
                                           inland; where rows were removed in
                                           front of the road it lands on the old
                                           pavement's first row or the row
                                           seaward of it (2026-09-08). Domains
                                           the footprint has no setback for keep
                                           v2's value.
    HAT_footprint_audit.csv                what was done to every domain
    RUN_MANIFEST.txt, README.md            provenance and the rules

THE RULES, as decided with Hannah 2026-09-07/08 (HAT_footprint_1984.py and
HAT_fill_copy_scope.py carry the argument; this script only applies them)
    N               n_cells in footprint_1984_by_domain.csv
    insert point    insert_row_behind_road: int(setback_new/10) + 2, behind
                    the model's two roadway rows AS PLACED under the 1984
                    setback; crest_row + 1 where there is no model road
                    (GIS 1-5, 8)
    add             block = z[r : r+N] copied cell by cell, inserted at r
    remove          rows r .. r+|N|-1 deleted, r = int(setback_v2/10) - |N|: the
                    |N| rows directly SEAWARD of today's roadway rows (Hannah,
                    2026-09-08: the rows come out of the interior in front of
                    the road); the road's cells and all behind them are kept
    dune array      unchanged
    setback         setback_new_m (the 1984 measurement); v2's value where absent

VERIFIED AFTER WRITING
    unchanged domains are byte-identical to v2; a changed domain has exactly
    rows_before + N rows, its rows before the insert point are identical to
    v2, and its block equals the rows that follow it (add) or its tail is
    v2's tail (remove).

USAGE
    python HAT_build_footprint_version.py --dst-version v3
    python HAT_build_footprint_version.py --dst-version v3 --overwrite
```

Notes that were in the code:

```text
--- the setback CSV: two rows, domain ids then values (the model-facing
format hatteras_site_config reads). Start from v2's so the domain list and
the format are exactly what the runner expects, then replace every value
the footprint has a 1984 setback for.
```

### 5-build/HAT_plot_version_figures.py

The figures a built dune-topo version gets, so it sits beside v1 and v2 with a figure set of its own.

From the script's original header:

```text
The figures a BUILT dune-topo version gets, so that v3 sits beside v1 and v2
with a figure set of its own (Hannah, 2026-09-09: "v3 should have figures
here as well").

WHY NOT THE EXTRACTOR'S FIGURES. v1 and v2 carry `figures/qc/` and
`figures/gis_vs_processed/`: the raw DEM profile against the extracted one,
per domain, and the island summary of windows and dune heights. Those are
figures OF AN EXTRACTION - they need the picks, the raw profiles and the
straightening frame. A built version is not extracted: it is its source
version with rows inserted or removed and the road setback re-set, so the
questions its figures answer are different - what does each domain look like
now, beside what it looked like before; where did the rows go; what changed
island-wide. Nothing here re-measures anything: every number is read from the
version's own arrays, its setback CSV and its footprint audit.

WHAT IS WRITTEN, into dune-topo/<version>/
    figures/grid/domain_NNN_grid_<version>.png    one per domain: the source
        version and this version side by side as the model holds them - the
        two dune rows on top (drawn at berm + dune height), every interior row
        down the page, elevation classes (m MHW), NC-12's two rows at each
        version's setback, and the footprint: the inserted block outlined
        (add) or the removed rows hatched in the source and the seam marked
        in the version (remove). Unchanged domains are drawn too, so the set
        is complete; their two panels are identical.
    HAT_dune_topo_summary_<version>.png           every domain on one page:
        interior rows, the road setback, mean interior elevation and mean
        dune height, source against version, with the communities banded.
        The counterpart of the extractor's summary page.
    HAT_dune_topo_island_planview_<version>_<year>_{trimmed,padded}.png
        the island in plan view at the period's dune offsets, in the
        extractor's poster style, with NC-12 drawn where the MODEL places it
        (the version's setback), not from the GIS mask - a built version's
        interior frame is no longer the mask's frame. The counterpart of the
        extractor's plan views.
    figures/README.md                             what these are and are not

USAGE
    python HAT_plot_version_figures.py                      # v3, source from its manifest
    python HAT_plot_version_figures.py --version v3 --source v2
    python HAT_plot_version_figures.py --domains 85,63      # only the grid panels of these
    python HAT_plot_version_figures.py --no-grid            # summary and plan views only
```

Notes that were in the code:

```text
What the two versions are called on the figures (no working vocabulary; the
version names themselves are in the file names and the README).
```

```text
a single-column figure per domain: the two panels side by side, the
height following the row count so a deep domain is not squashed
```

```text
one legend under the figure, two columns: the elevation classes (m MHW)
and the footprint marks; a single column is too narrow for two legends
```

```text
the villages as light bands, in the canvas's column frame: a translucent
white so they read on the ocean colour and vanish under the island
```

```text
the extractor's cmap paints masked cells in the ocean colour, which would
cover the bands; here the axes background is the ocean and the mask is clear
```

### 6-result/HAT_compare_versions.py

v2 against v3 under the same hindcast: what the 1984 reconstruction changes in what the model does.

From the script's original header:

```text
v2 against v3 under the same hindcast, side by side: what the 1984
reconstruction changes in what the model DOES. One figure and one table per
run pair, read from the two runs' saved state (the .npz) and their metadata;
nothing is re-run and nothing is re-scored.

THE PAIRS (Hannah's advisor, 2026-09-09: "a comparison between v2 and v3
under the modules' automatic behaviour, full management and calibrated")
    emergent    HAT_1984_2004_calibBE_road_bdm_groin
                full management, calibBE, groin on, the roadway and
                beach-dune modules acting on their own (no prescribed
                relocations). Both versions in output/raw_runs/version-pair/<v>/,
                run by HAT_run_version_pair.py on the same code the same day.
                (The earlier pair - v2 in the calibration tree of 2026-09-07,
                v3 in behindroad-copy of 09-08 - sat on different commits, with
                the pipeline and the live setback CSV changed between them, so
                it was not a clean pair and was re-run.)
    prescribed  HAT_1984_2004_calibBE_road_reloc_bdm_groin
                the same with the recorded 1989 (GIS 84-87) and 1999
                (GIS 9-14) relocations imposed: the control, showing what
                the v3 setbacks change when the module is not deciding.
                Both versions in output/raw_runs/version-pair/<v>/, run by
                HAT_run_version_pair.py --relocations 1.

THE FIGURES, per pair (two double-column figures since 2026-09-10, in the
house style of hat_figure_style: v2 grey C["BASE"], v3 purple C["ACCENT"],
the recorded events C["REF"]; no title sentence on the canvas, the run name
and the pair go to CAPTIONS.md beside the PNGs)
  HAT_compare_v2_v3_<pair>_relocations.png
    (a) the NC-12 setback the model starts with, per road domain, two thin
        lines with markers on a symlog axis
    (b) every year the model relocated NC-12, per domain, the recorded
        events outlined; counts in the legend. In the prescribed pair the
        1989/1999 rows are inputs, so only the OTHER relocations are the
        module's own.
    (c) the number of relocations per year, island-wide, through time
  HAT_compare_v2_v3_<pair>_geometry.png
    (a) interior width per domain at 1984 (dashed) and 2004 (solid), both
        versions: what the footprint added or removed, and what the run
        then did with it
    (b) island-mean interior width through time, (c) the difference
        v3 - v2 on its own axis, (d) island-total cumulative overwash
    (e) mean interior elevation of the land cells at 2004 per domain, both
        versions, and (f) the difference v3 - v2 on its own axis
    The island-wide shoreline skill of both runs (they are near-identical by
    construction: the shoreline offset does not read the topography) is in
    the table and the report, not on the figure.

THE TABLE  version_compare_<pair>.csv, one row per domain: initial setback,
    relocation years and count, drowned, interior rows and width at 1984
    and 2004, mean land elevation at 1984 and 2004, cumulative overwash -
    for v2, for v3, and the difference.
THE REPORT  HAT_compare_versions.txt: the skill of the four runs, the
    relocation counts and their timing against the recorded events (mean
    error over the event blocks; a domain that never relocated is censored
    and counted, not averaged), and the island-wide geometry medians.

UNITS. Barrier3D stores decametres: widths and elevations x 10 -> m;
QowTS is dam^3 per dam of shoreline per year, x 100 -> m^3/m. The buffer
is 15 domains: GIS g is index g + 14.

USAGE
    python HAT_compare_versions.py                 # both pairs, whatever exists
    python HAT_compare_versions.py --pairs emergent
```

Notes that were in the code:

```text
RESOLVED, NOT JOINED: run_layout knows where a run folder keeps each of
its files, in the new layout and the old flat one alike.
```

```text
a PRESCRIBED event is applied as a displacement of the setback and
does not raise the relocated flag: read it off the setback series
as the jump at the event year in the event's block
(displacements run 17-108 m; setbacks quantise to 10 m; an emergent
relocation the same year would carry the flag and is excluded)
```

```text
Two double-column figures per pair (2026-09-10, the house style): the road
(setback, relocations, relocations per year) and the geometry (widths,
overwash, elevation, each difference on its own axis). The run name and the
pair description go to CAPTIONS.md beside the PNGs, not onto the canvas.
```

```text
a prescribed event is an input: a tick in the run's colour across the
recorded bar where it was applied, one per version, side by side
```

<details><summary>Function notes (the original docstrings)</summary>

**`timing_score()`**

```text
Mean error of the first relocation against the recorded event year, per
block; domains that never relocate are censored, counted, not averaged.
Meaningful only when the events are NOT prescribed (scoring a run against
its own input is circular); the prescribed pair reports counts instead.
```

**`own_relocations()`**

```text
Relocations that are the module's own in the prescribed pair: everything
outside the (event year, event block) pairs.
```

</details>

### 6-result/HAT_plot_b3d_grid.py

The Barrier3D grid itself, for the extraction and the layer built from it, side by side.

From the script's original header:

```text
The Barrier3D grid itself -- the arrays CASCADE is handed -- for the
extraction and the layer built from it, side by side (--base/--insert).

WHAT IS DRAWN
    Every cell the model receives for a domain, in the model's own indexing:

        dune rows      DuneDomain, 2 rows, height above berm -> m MHW
        interior rows  InteriorDomain, row 0 seaward, increasing landward
        NC-12          the 2-cell block bulldoze() writes, at
                       road_start = int(setback / dy)

    Cross-shore runs DOWN the page (row 0 at the top, behind the dune), and
    alongshore runs across -- 50 cells, 500 m.

    The beach and shoreface are NOT in these arrays. Barrier3D carries them as
    parameters, not cells, so they are named in the caption but not drawn; the
    top of the dune strip is where the model's grid begins.

WHY BOTH VERSIONS SIDE BY SIDE
    v2 prepends N interior rows behind the dune so the 1984 roadway starts where
    it historically did. In the grid that shows as the road block moving DOWN the
    page: the dune does not move relative to the array, the ground between the
    dune and the road grows.

    At GIS 85 the difference is stark. In v1 the road occupies rows 0-1 -- on
    interior row 0, which is the dune crest itself -- because the setback
    measured -10 m and was floored to 0. In v2 it sits at rows 5-6 on backdune.

USAGE
    python HAT_plot_b3d_grid.py
    python HAT_plot_b3d_grid.py --domains 85 --rows 40
```

Notes that were in the code:

```text
Defaults, overridable per run. HARDCODING THESE IS WHAT WENT WRONG: the
script was named _v1_v2 and defaulted its output to that name, then was
repointed at v3/v4 -- so a default run would have written v3/v4 content
into a file called v1_v2 and quietly replaced a correct figure.
Default v3 -> v5 since 2026-09-03. It was v3 -> v4, which meant a bare
run wrote HAT_b3d_grid_v3_v4.png - a figure deliberately deleted as
superseded, so the default recreated the thing the cleanup removed.
DELETED 2026-09-07 with every layer (only unmodified topography is kept);
the literal is kept as the name of what this drew. require_version() in
main() says so before any array is opened.
```

```text
Named for what it DRAWS, so it cannot silently replace another pair's
figure.
```

```text
NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
```

<details><summary>Function notes (the original docstrings)</summary>

**`draw()`**

```text
One domain's grid: dune strip on top, interior below, NC-12 outlined.

The per-panel "+N rows" count came off the canvas 2026-09-10 with the
title's setback figures: the house rule puts statistics in the caption,
and main() writes one caption clause per domain from the same numbers.
```

</details>

### 6-result/HAT_plot_footprint_result.py

The first result of the 1984 footprint: the hindcast on v3 against the same run on v2, at the road.

From the script's original header:

```text
The first RESULT of the 1984 footprint: the hindcast on v3 (rows behind the
road, copy fill, 1984 setbacks) against the same run on v2, at the road.

Reads the two runs' saved model state (the roadway objects in the .npz) and
draws, for every road domain, the setback the model started with and every
year it relocated NC-12, the two versions together. The island-wide skill of
both runs is in the caption beside the figure, not on it: the figure is drawn
double-column in the house style of hat_figure_style (v2 grey C["BASE"], v3
purple C["ACCENT"], the recorded events C["REF"]).

    v2   output/raw_runs/1984_2004/calibBE/<run>              the calibration tree
    v3   output/raw_runs/behindroad-copy/1984_2004/calibBE/<run>   arm behindroad-copy

USAGE
    python HAT_plot_footprint_result.py
```

Notes that were in the code:

```text
Asked for rather than spelled: the arm moved under arms/ on 2026-09-10 and
a hand-built path missed it. find_run_dir reads either layout.
```

```text
RESOLVED, NOT JOINED: run_layout knows where a run folder keeps each of
its files, in the new layout and the old flat one alike.
```

### 6-result/HAT_plot_insert_three_scales.py

The seaward row insert at three scales: GIS 85, the two relocation blocks, and the whole island.

From the script's original header:

```text
The seaward row insert at three scales: GIS 85, the two relocation blocks, and
the whole island.

WHAT THE INSERT IS
    The 1984-start DEM is a 1996 ALACE beach on a 2009 backdune, so its dune has
    already migrated landward past the 1984 NC-12 alignment. At GIS 85 that puts
    the 1984 roadbed SEAWARD of interior row 0 -- setback -10 m, floored to 0,
    and a road that relocates in model year 1 by construction.

    The fix measures how far the dune line moved between 1984 and 1997 (both
    digitized from imagery, same feature) and prepends that many interior rows
    behind the dune, so row 0 sits at the 1984 dune position and NC-12 lands its
    true distance behind it.

WHY THREE SCALES
    Row 1 shows the mechanism on the domain the work was for.
    Row 2 shows every domain with a documented historical relocation, which is
        where the correction is actually applied.
    Row 3 shows the measured retreat for all 90 domains -- context for whether
        GIS 85 is exceptional or typical. NOTE it is MEASURED island-wide but
        APPLIED only to the ten block domains; the row 3 bars outside the shaded
        blocks are what an island-wide version WOULD insert, not what it did.

USAGE
    python HAT_plot_insert_three_scales.py
```

Notes that were in the code:

```text
Resolved through hat_topo_version.duneline_shift_dir - ONE definition
of a path that eight scripts used to build by hand. Moved under
2-domain-reconstruction-1984/ on 2026-09-03.
```

```text
DELETED 2026-09-07 with every layer (only unmodified topography is kept);
the literal is kept as the name of what this drew. require_version() in
main() says so before any array is opened.
```

```text
scope) until 2026-09-03. N is identical at the ten
BLOCK domains, so panels (a)-(c) are unchanged - but
panel (d) is NOT: it goes from 8 red bars to 38,
because v5 applies rows wherever the measurement
selects them. Every label below is derived from
VERSION and from the data, never written literally,
so switching the version cannot leave a stale caption.
```

```text
Shared vocabulary: BASE = the unmodified input, ACCENT = the change
under test, ADDED = fabricated ground. Same meanings in every figure.
```

```text
Named for the versions drawn. The old default was a fixed
"HAT_insert_three_scales.png", so repointing BASE_VERSION/VERSION from
v1/v2 to v3/v4 would have overwritten the v1/v2 figure in place.
```

```text

The monospace block of GIS 85 numbers that used to sit at gs[0, 1] came
off the canvas 2026-09-10: a statistics table is caption text under the
house rules, and it is now the last third of caption() below. Panel (a)
takes the full width in its place.
```

```text
NOTE 2026-09-08: this figure now lives in figures/superseded-layers/; the script is
guarded (no layer on disk), so nothing is written here until a layer is rebuilt.
```

### 6-result/HAT_plot_method_compare.py

One domain, two methodologies side by side in the model's frame: the v2 extraction and the v3 footprint.

From the script's original header:

```text
One domain, two methodologies, side by side in the model's frame:

    v2   (default base) the re-pick extraction v3 is built on. The 1984 road
         was measured against interior row 0, came out NEGATIVE (seaward of
         row 0, in the dune), and was FLOORED to 0 - so the model placed NC-12
         on rows 0-1, at the dune, and relocated it in year 1. No rows added.
         Same pick set as v3, so the panels differ ONLY by the rows and the
         setback. (--base v1 shows the original extraction instead, which
         also differs by the 2026-09-02 re-pick.)
    v3   the re-pick base (v2) + the symmetric 1984 footprint placed directly
         behind the road AS PLACED under its 1984 setback (road against the
         1984 dune line, row-0 convention; no floor) and filled by copying the
         rows that follow. The road sits on measured cells; the block is behind it.

Each panel: the two dune rows (berm + dune height) on top, the interior below
in elevation classes, a metres axis on the right; NC-12 as the model places it
(dark band); the measured 1984 road position (outlined); and the inserted rows
outlined in the accent colour. The interior depth, the retreat the block stands
for and what each version's hindcast did with the road (v1: arm
pea1989basenoreloc; v3: arm behindroad-copy; both calibBE, full management,
prescribed relocations off) go to the CAPTIONS.md beside the figure, not onto
the canvas; the figure is drawn double-column in the house style of
hat_figure_style.

USAGE
    python HAT_plot_method_compare.py                 # GIS 85
    python HAT_plot_method_compare.py --domain 86 --rows 40
```

Notes that were in the code:

```text
v1 named dunestart_offset_ARCHIVE_1984start_v1/, which became a dated
superseded folder; resolved since 2026-09-18.
```

### 6-result/HAT_plot_seaward_insert_compare.py

Two 1984-start versions drawn against each other in a common frame (default: the two extractions, v1 and v2).

From the script's original header:

```text
Two 1984-start versions drawn against each other in a common frame.

    Default pair since 2026-09-07 -- the two EXTRACTIONS, the only versions kept:
    v1   the original pick set (2026-08-27). GIS 85 setback -10 m, floored to 0.
    v2   the re-pick with NC-12 visible (2026-09-02); what CURRENT names.

    The script was written for the seaward-row insert (an extraction against a
    layer with N rows inserted behind the dune, N measured as the 1984-1997
    dune-line difference). The layers v3-v8 were deleted 2026-09-07 -- only
    unmodified topography is kept -- so no insert pair exists on disk.
    Any pair that does can be compared with --versions "a:label;b:label".

Land width is drawn with BARRIER3D's definition (stop at the first cell at or
below sea level), not a count of dry cells, so the panel agrees with what the
model actually computes. v2 preserves it exactly.

Everything is drawn in a COMMON frame: distance landward of v1's interior row 0.
A variant whose row 0 has moved seaward therefore starts at negative x, and the
fabricated ground is the part left of zero. Plotting each variant from its own
row 0 would hide exactly the thing being compared.

USAGE
    python HAT_plot_seaward_insert_compare.py
    python HAT_plot_seaward_insert_compare.py --domains 85,86
```

Notes that were in the code:

```text
Drawn thick to thin, because the versions coincide over most of the profile
and equal linewidths would show only the last one drawn.
The default pair is v1 vs v2, the two extractions (2026-09-07). It was v2 vs
v4 (base vs the measured+floor layer) until the layers were deleted that
day, and v1 vs v2 in the pre-re-pick numbering before 2026-09-03. The three
width variants this script was originally written for -- v1_pad_measured,
v1_translate_measured, v1_none_measured -- were DELETED on 2026-09-02: they
predated the island-width fix, so all three behaved as `pad`, and no run was
ever built from them. Pass --versions to compare anything else.
```

```text
The origin shift is read off the ARRAY, not guessed from the
version name: a variant that inserted rows is taller than v1 by
exactly n. A name-suffix test silently plotted later variants
unshifted -- superimposing them on v1 and hiding the comparison.
```

```text
THE ROAD DOES NOT MOVE. Its position in this common frame is
v1's own raw setback; what every variant changes is where row 0
sits relative to it. Drawing one block per variant, as an
earlier version of this figure did, draws the opposite claim.
```

<details><summary>Function notes (the original docstrings)</summary>

**`b3d_width()`**

```text
Island cells per alongshore column, BARRIER3D'S definition.

Reproduces FindWidths (Barrier3D/barrier3d/barrier3d.py:29): walk landward
from row 0 and stop at the FIRST cell at or below sea level. Anything past
an interior water gap is not island.

This panel used to count every dry cell in the column instead, which is a
different number -- on GIS 85 it disagrees in all 50 columns, median 44
against 37.5. That is the same mistake that made `translate` silently behave
like `pad`, and drawing it here would have shown v2's width changing when the
model sees it as identical to v1's.
```

**`baseline_setbacks()`**

```text
v1's own raw setbacks, from the road-offset measurement.

v1 has no insert audit -- it is the thing the inserts are measured against --
so its baseline number has to come from the file that produced it, or the
comparison table prints nan in the row the reader most needs.
```

</details>

### 6-result/HAT_run_version_pair.py

Run one hindcast scenario on two dune-topo versions identically, so they differ only in topography and setbacks.

From the script's original header:

```text
Run ONE hindcast scenario on two dune-topo versions, identically, so the two
runs differ in nothing but the topography and its setback CSV. The pair the
advisor asked for (2026-09-09): v2 (the extraction) against v3 (the 1984
reconstruction), full management, calibBE, groin on, under the modules'
automatic behaviour - and the same pair with the recorded 1989/1999
relocations prescribed, as the control.

WHY A SCRIPT (the same reason HAT_run_row_insert_set.py is one). Two pieces
of global state select a version and both must be put back: the forcing-tree
setback CSV, which hatteras_site_config.py hardcodes, is copied per version
from dune-topo/<version>/ and restored in `finally`; the topography version
goes through HAT_TOPO_VERSION_1984_START, which outranks CURRENT and dies
with the subprocess. hat_run.yaml is ignored (HAT_IGNORE_SETTINGS=1).

WHERE THE RUNS LAND
    output/raw_runs/version-pair/<version>/1984_2004/calibBE/<run_name>/
    via HAT_RUN_KIND=version HAT_RUN_TAG="version-pair/<version>", so the
    runs file under raw_runs/versions/version-pair/<version>/ (2026-09-16).
    The run name is the same for both versions by design; the tag tells them
    apart, on disk and in the `kind`/`tag` columns of run_index.csv.

USAGE
    python HAT_run_version_pair.py --relocations 1          # the prescribed control
    python HAT_run_version_pair.py                          # emergent (the modules decide)
    python HAT_run_version_pair.py --versions v2,v3 --dry-run
```

### 6-result/HAT_version_pair_gif.py

v2 beside v3 through time: the relocation animations, with the two panels as the two dune-topo versions.

From the script's original header:

```text
v2 beside v3 through time: the animations of the relocation comparison, but
with the two PANELS being the two dune-topo versions under ONE run scenario,
instead of the two relocation arms under one version. This is the view that
shows what the inserted and removed cells did (Hannah, 2026-09-09: "I want
to see how the inserted cells affected things").

WHAT IS DRAWN, per scenario (emergent: the modules decide; prescribed: the
recorded 1989/1999 relocations imposed) and per alongshore window:
    road_topography_<window>.gif   Barrier3D's own interior grids painted
        year by year, NC-12 on them, v2 left and v3 right, one colour scale
        and one year clock. Where v3 added rows the island is wider behind
        the road from year 0; where it removed them, narrower in front.
    road_relocation_<window>.gif   the dune line and the road as lines,
        landward-positive from each domain's year-0 dune line, a star where
        the module relocated, a ring where a prescribed move was applied.
One folder per place - the whole island, the two event blocks, and the two
reaches where the footprint is largest, Pea Island (GIS 78-87, rows added) and
the Avon-Tri-Village removals (GIS 62-68) - with `topography.gif` and
`dune-and-road.gif` in each.

Everything is read from the runs' saved state through the comparison
script's own loaders; nothing is re-run. The makers are the ones the
relocation comparison uses (cascade_pipeline.plotting.road_relocation_gif);
only the panel labels and the pairing differ.

WHERE  output/comparisons/relocation/versions/v2_vs_v3/<scenario>/<place>/
       (its own folder: it is a cross-version comparison, not a set of one
       version, so it does not belong under 1984_2004/v2/ or v3/)

USAGE
    python HAT_version_pair_gif.py                       # both scenarios
    python HAT_version_pair_gif.py --scenarios emergent
```

Notes that were in the code:

```text
the relocation comparison moved into hatteras_ms/experiments/ on 2026-09-13
(cfd0b475); this import broke silently until 2026-09-17
```

```text
ONE FOLDER PER PLACE, two files in each (Hannah, 2026-09-09: "organize the
figures better"): a reader opens the reach they care about and finds both
views of it side by side. The relocation-comparison windows first, then the
two reaches where the footprint is largest.
(folder, title, topography window, line window)
```

```text
both panels are the SAME scenario: neither run carried a prescribed
move in the emergent pair, both did in the prescribed one
```

### 6-result/HAT_version_pair_report.py

v2 against v3 in one report: the relocation comparison, with the two versions side by side in every section.

From the script's original header:

```text
v2 against v3 in ONE report: the relocation comparison's report.txt, but with
the two dune-topo versions side by side in every section instead of one
version per file. The companion of HAT_version_pair_gif.py, which does the
same for the animations (Hannah, 2026-09-10: "a txt report directly comparing
the two versions, similar to v3/calibBE_groin/report.txt").

WHAT IS READ. The two per-version sets that HAT_relocation_comparison.py
wrote, v2/calibBE_groin/tables/ and v3/calibBE_groin/tables/, the run
metadata of the four runs behind them, and the v3 footprint audit
(dune-topo/v3/HAT_footprint_audit.csv: which domains got rows and how many).
Nothing is re-run and nothing is re-scored: every number here is one of
theirs, or a difference of two of theirs. A section whose numbers disagree
with the per-version reports means one of the three is stale.

WHAT IS WRITTEN, all under output/comparisons/relocation/versions/v2_vs_v3/
    report.txt          the console output of this run, with the provenance
                        of the four runs and the two table sets above it
    tables/*.csv        the side-by-side tables the report prints, in full
                        (the report prints the historical domains and the
                        domains that differ; the CSVs hold every road domain)

USAGE
    python HAT_version_pair_report.py
    python HAT_version_pair_report.py --preset calibBE --set calibBE_groin
```

Notes that were in the code:

```text
the relocation comparison moved into hatteras_ms/experiments/ on 2026-09-13
(cfd0b475); this import broke silently until 2026-09-17
```

<details><summary>Function notes (the original docstrings)</summary>

**`read_set()`**

```text
The seven CSVs of one per-version comparison set, plus the 'generated'
stamp of the report beside them.
```

**`read_audit()`**

```text
v3's footprint audit: one row per GIS domain, n_cells signed (+ added,
- removed, 0 unchanged).
```

**`_tab()`**

```text
A ruled table: pandas' own per-column formatting, columns separated by
' | ' and a rule under the header (Hannah, 2026-09-10: lines between the
columns so the tables read more easily). Right-aligned like to_string.
```

</details>
