# figure_making/island — the island as the model starts it

The site and the model's starting surface: where Hatteras is, how the 90
domains tile it, the t=0 elevation of every domain, and one domain as a grid.

```
study_area_figures.py           the generic site figures: study area, domain framework,
                                one domain, forcing timelines, dune lines (1-site/,
                                2-observations/, 3-model-inputs/)
initialization_figures.py       the t=0 elevation surface, per start year and treatment
                                (3-model-inputs/1-domains/initial_island/)
visualize_domain_topography.py  one domain's topography as a heatmap, from a run's .npz
```

What not to trust: `visualize_domain_topography.py` names a run
(`HAT_1978_1997_natural`) that has been deleted, so it cannot run until
`CASCADE_OUTPUT_FILE` points at a current run's `.npz`. `study_area_figures.py`
fetches satellite tiles unless run with `--vector`.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### initialization_figures.py

The island as the model starts it: the t=0 elevation surface of every domain, one figure per start year.

From the script's original header:

```text
The island as the model starts it: one page figure per orientation year in
YEARS, showing the t=0 elevation surface every domain is initialised from.

    python scripts/figure_making/island/initialization_figures.py

Writes to output/figures/3-model-inputs/1-domains/initial_island/<year>/<scheme>/ (ORGANIZATION.md
rule 1), PNG at the top of the scheme folder and PDF + CAPTIONS.md under its
supporting/, through the house `save()`. Three figures per folder:

    island_<year>.png                  DETRENDED, three panels: (a) the BRIE
                                       shoreline offset, (b) the 90 real
                                       domains, (c) the same with the 15
                                       buffer domains at each end
    island_<year>_absolute.png         ABSOLUTE, one panel: every domain at
                                       its true offset, 90 real domains
    island_<year>_absolute_buffers.png ABSOLUTE, one panel: the same with the
                                       buffer domains

...for each of the four hindcast starts (1984, 1996, 2004, 2010) in each of
the two elevation treatments (classes, terrain): 24 figures from two
topography loads. The views are the same data placed different ways, and each
caption points at its companions. SINCE 2026-10-01 TERRAIN ONLY: the house
style switched domain topography from classes to terrain, so the classes/
treatment was dropped (12 figures, <year>/terrain/).

ONE FOLDER PER START, THEN PER TREATMENT
    Two dozen figures in one folder read as two dozen unrelated images; the
    pairing that matters is the views of ONE start, in ONE treatment. So
    <year>/<scheme>/, each with its own CAPTIONS.md, and the file names are
    identical across schemes -- only the path says which treatment it is
    (Hannah, 2026-09-17). The names keep their year, because a figure pulled
    out of its folder must still say which start it is.

WHY THE ABSOLUTE VIEW IS TWO FILES, NOT TWO PANELS
    The real-domain map and the with-buffers map were panels (a) and (b) of
    one figure until 2026-09-17. Each is a map of the whole reach in its own
    right and each is used on its own, so as a pair they cost the page two
    half-height panels to say one thing twice. Split, each gets the full
    column width and its own cross-shore window -- the real-domain map no
    longer reserves the extra kilometre the buffer offsets need. The price is
    that the two no longer share metres per inch, so each caption states its
    span and its exaggeration.

WHICH SURFACE EACH YEAR GETS
    Not one. 1984 and 1996 read the 1984-start extraction (2009+2014 with the
    1996 ALACE graft); 2004 and 2010 read 2004-start (2009+2014 alone).
    YEAR_PRODUCT in hat_topo_version pairs them and product_for_year() is the
    accessor. Until 2026-09-17 this script resolved topo_dirs() ONCE at import,
    so it drew the 2004-start island under a 1984 label -- the exact failure
    product_for_year()'s docstring was written to prevent.

WHY ONE FIGURE DRAWS THE ISLAND LEVEL, NOT AS A DIAGONAL
    Each domain sits at its own BRIE shoreline offset, and those offsets span
    about 6 km across the reach. Drawn in absolute cross-shore space the 2 km
    island became a thin diagonal ribbon crossing an 8 km canvas that was
    three-quarters empty water, and no cross-shore detail survived at page
    width. So the map panels are DETRENDED: every domain is drawn against its
    own frame origin, which is what its elevation array actually holds, and
    the offset that would have displaced it is drawn as panel (a) above, on
    the same alongshore axis. Nothing is hidden -- the offset becomes a number
    you can read off an axis instead of a slope you have to estimate
    (Hannah, 2026-09-17).

    The absolute view is kept alongside it, restyled, because a single canvas
    IS the initial condition and the detrended panels are a rearrangement of
    it (Hannah, 2026-09-17). It carries about 9 km of cross-shore, so it takes
    its own exaggeration -- see VERTICAL_EXAGGERATION_ABSOLUTE.

    Panels (b) and (c) share one alongshore scale in metres per inch, so the
    real span in (b) sits directly above its own position in (c) and the
    buffer domains are visibly the part that sticks out.

THE CROSS-SHORE IS EXAGGERATED, AND SAYS SO
    A 60 km reach and a 2 km island cannot share a scale on a 190 mm page:
    at 1:1 the map panels would be 0.22 in tall. VERTICAL_EXAGGERATION and
    VERTICAL_EXAGGERATION_ABSOLUTE set the factor for their figure, are the
    only place a panel height comes from, and both captions state them. The
    old poster stretched the cross-shore by an unstated `fig_w * aspect * 1.8`
    capped at 7.5 in.

HOUSE STYLE
    Elevation is drawn in CLASSES, not a ramp (`elevation_cmap()`): the back
    barrier sits a few decimetres below MHW and the dune is metres above it,
    so a continuous ramp renders the whole island as one flat tone. Until
    2026-09-17 this script used `plt.cm.terrain` under a `FuncNorm`, which
    also gave the figure TWO colours for water -- the -3 m sentinel fill came
    out terrain navy, the uncovered canvas came out house blue -- so the
    padded back-barrier frame read as deep ocean. One class, one colour now.

    The terrain treatment kept beside it obeys the SAME water rule: every cell
    below 0 m MHW is masked and painted C["WATER"], and the ramp is truncated
    to the 0-4 m land range it actually draws from, so its colourbar cannot
    advertise a blue no cell uses. The two schemes differ only in how they
    colour LAND (Hannah, 2026-09-17).

    The figure is sized by `figsize()` at the width it will be printed, its
    caption lives in CAPTIONS.md rather than on the canvas, and the local
    rcParams block that overrode the house ink with '#1a1a2e' is gone.
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5), so this block is
independent of whatever this script calls its own repository variable.
```

```text
Derived from this file's location (scripts/figure_making/island/) so the
script runs on any checkout without editing a hardcoded path.
```

```text
Orientation years to render: every hindcast period start. Offsets come from
2-brie-offset/<year>/, the same PADDED_120 files the hindcast run script
initializes from, and the topography from that year's product (below).
```

```text
WHICH EXTRACTION -- resolved, not pinned, the same way the hindcast runner and
the groin sweep worker resolve it. This was hardcoded '2009_v2', a directory
that has since been moved into 2009-dune-topo/incorrect/, so the poster either
failed to load or drew arrays the run does not use. topo_dirs() reads VERSION
out of HAT_dune_topo_extractor.py, so the figure always shows the surface the
model is actually initialised from.
```

```text
Paths come from what topo_dirs() RETURNED, not re-joined from parts - the
tree went period-first on 2026-08-25 and re-joining would have kept pointing
at a folder that no longer exists. topo_dirs() with no product resolves
2004-start, which is the surface this poster showed before the restructure.
```

```text
array_name() is the single definition of these filenames - the same one
the extractor writes with. Nothing here spells a name.
```

```text
...AND THE PRODUCT IS RESOLVED PER YEAR, NOT ONCE. The four starts do not
share one surface: 1984 and 1996 read 1984-start (the 2009+2014 mosaic with
the 1996 ALACE graft), 2004 and 2010 read 2004-start (2009+2014 alone).
YEAR_PRODUCT in hat_topo_version is the single definition of that pairing,
and product_for_year() its accessor -- whose docstring names precisely the
bug this script had until 2026-09-17: "a loop body that resolves topo_dirs()
once, outside the loop, and silently gives every year the same interiors."
It did, so the 1984 figure was drawing the 2004-start island.
```

```text
The extractor works in a 200-row (2000 m) cross-shore frame but writes the
topography trimmed to the interior, so the back-barrier water rows are absent
from the .npy files. Refill them with the extractor's water sentinel so every
domain spans the full frame. Both values come from RUN_MANIFEST.txt
(TOPO_ROWS, SENTINEL_WATER_M / WATER_CLAMP_M) in the version folder.
```

```text
HOW FAR THE CROSS-SHORE IS STRETCHED. The alongshore scale is fixed by the
page (120 domains across a 190 mm column); this is the only other free
number in the layout, and the caption states it. 8 puts the map panels at
about 1.75 in, enough to read a dune line against a back barrier.
```

```text
The absolute-placement figure carries about 9 km of cross-shore against the
detrended figure's 2, so it cannot use the same factor and still fit a page:
8x would make one of its panels 7.9 in tall. 3x puts them near 3 in.
```

```text
Compositing (padding, unit conversion, canvas assembly) is shared with the QC
notebook via cascade_pipeline.plotting.init_planview; only the page-figure
styling below is local to this script.
```

```text
build_canvas leaves max(offset) + topo_rows + 5 rows; with no offsets
that is five rows of NaN above the frame.
```

```text
One alongshore scale, in metres per inch, shared by both map panels and the
offset panel. Everything else in the layout is derived from it, so the panels
cannot disagree about where a domain is.
```

```text
GIS-domain coordinates of the padded array: index 0 is the 15th buffer south
of GIS 1, so domain i sits at i - NUM_BUFFER_DOMAINS + FIRST_FILE_NUMBER.
```

```text
WATER IS MASKED, NOT COLOURED. Every cell below 0 m MHW -- real
back-barrier, and the -3 m sentinel the extractor's trimmed rows are
refilled with -- goes to `set_bad`, which both schemes set to
C["WATER"]. Masking rather than letting each colormap render its own
low end is what makes the two treatments agree about water; it also
matches the class scale, whose first bound is [-99, 0), so nothing
changes in the class figures. Strictly BELOW: a cell at exactly 0.0 m
falls in the [0, 0.5) land class, and it still does here.
```

```text
The house rule is CLASSES (hat_figure_style: a hard break at 0 m, one colour
for water, so the only distinction that matters -- which cells are land -- is
the one you see first). The terrain ramp the figure used until 2026-09-17 is
kept beside it, at Hannah's request, because a continuous ramp shows the
smooth cross-shore gradient a class boundary hides, and because it is what
the extractor's own QC view still draws.

They are NOT mixed in a folder: a reader who sees both treatments of one
island in one place has to work out that they are the same data. Each gets
its own <year>/<scheme>/ with its own supporting/ and CAPTIONS.md, so the
file names are identical and only the path says which is which
(Hannah, 2026-09-17).
```

```text
ONE WATER RULE, BOTH SCHEMES (Hannah, 2026-09-17). Water is every cell below
0 m MHW and it is C["WATER"], full stop -- the class rule, applied to the
ramp as well. `_map_panel` masks those cells and `set_bad` paints them, so
both treatments take the identical path and cannot drift apart.

Terrain used to draw its own bottom (navy) for everything under its -1 m
floor, which on this canvas is the -3 m sentinel that refills each domain's
trimmed rows. That put TWO blues in one image -- navy for the padded
back-barrier frame inside a domain, pale blue for canvas no domain covers --
and the navy read as deep ocean. Now there is one.

Because water no longer comes off the ramp, the ramp must not claim it: the
terrain colormap is TRUNCATED to the land part it actually uses, its old
sea-level position upward, under a plain 0-4 m Normalize. Land colours are
unchanged -- the old FuncNorm sent elevation e >= 0 to colormap position
SEA_LEVEL_POS + (1 - SEA_LEVEL_POS) * e / ELEV_MAX_M, which is exactly what
the truncated map under a linear norm does -- but the colourbar no longer
shows a blue band no cell is drawn from.
```

```text
THE TERRAIN SCHEME KEEPS TERRAIN'S OWN NAVY (Hannah, 2026-09-17). It is
plt.cm.terrain(0.0), the floor the ramp used to render the -3 m sentinel at,
and it is what the figure looked like before this rewrite. The pale house
water was tried first and sits too close in value to the saturated green the
ramp starts at for the shoreline to read.

What is NOT restored is the old two-blue split. Back then navy came off the
ramp (the padded back-barrier frame inside a domain) while the uncovered
canvas came off set_bad in house blue, so one image had two colours for
water and the navy read as deep ocean. Water is masked under one rule now,
so this navy is every cell below 0 m MHW and nothing else. The class figures
keep C["WATER"] exactly; the schemes share the rule, not the shade.
```

```text
WHAT IS DRAWN ON TOP OF THE WATER FOLLOWS IT. The dashed lines bracketing the
real span cross open water for most of their length, so their colour is a
property of the water they sit on, not of the figure: near-black INK on the
pale class blue, white on the terrain navy. Picking one ink for both would
lose the line in whichever scheme disagreed (Hannah, 2026-09-17).
```

```text
bounds[0] is the water class, which no cell reaches any more; the
ticks start at the 0 m break either way.
```

```text
The reach as one map: every domain at its true BRIE cross-shore offset, so
the seaward bend from Pea Island to Cape Point is geometry rather than a
curve on a separate axis. This is the view the figure had before 2026-09-17,
kept because a single canvas IS the initial condition and a detrended panel
is a rearrangement of it -- but drawn to the same rules as everything else.
```

```text
The two buffer ramps sit at opposite ends of the cross-shore window
-- the south one high, the north one low -- so each tag goes to the
corner its own ramp leaves empty rather than onto the island.
```

<details><summary>Function notes (the original docstrings)</summary>

**`dune_offset_file()`**

```text
Padded BRIE dune-offset CSV for one orientation year, THROUGH THE
VERSION THE RUN READS.

Every start's offsets were put under `<year>/v<n>/` with a CURRENT marker
on 2026-09-15, so the flat `<year>/Island_Dune_Offsets_...csv` this used to
read no longer exists and the poster could not be redrawn. The version is
resolved the way the runner resolves it, so a re-version moves this with
it rather than breaking it again.

Since 2026-09-18 that is literally true: hat_topo_version.offset_file is
the function the runner's config calls. The copy that stood here ignored
HAT_OFFSET_VERSION_<year> and fell back to the newest v<n> where the
runner raises.
```

**`elevation_file_paths()`**

```text
The 120 padded-order array paths for one product's topography dir.

The 15 domains at each end are the SAME sampled array repeated, which is
why the buffers read as a regular comb in the absolute figure.
```

**`product_for()`**

```text
(grids, product, version, interior row range) for one start year.

Loaded once per PRODUCT, not once per year: 1984 and 1996 share a surface,
as do 2004 and 2010, so four figures cost two loads.
```

**`detrended_canvas()`**

```text
The composited surface with every domain on its OWN frame origin.

Passing zero offsets to the shared compositor is exactly the detrending:
`build_canvas` places domain i at row `offset_cells[i]`, so zeros put each
one where its elevation array starts. The offsets themselves are drawn in
panel (a) instead of being spent on 6 km of empty canvas.
```

**`absolute_canvas()`**

```text
The composited surface with every domain at its TRUE BRIE offset.

The initial condition as one map, which is what the figure showed before
2026-09-17 and what `island_<year>_absolute.png` shows again: the reach
bends seaward by about 6 km from Pea Island to Cape Point, and here that
bend is geometry rather than a curve on a separate axis.
```

**`_map_panel()`**

```text
One elevation strip, seaward at the bottom.

The image extent comes from the canvas's OWN row count, so a cell always
lands at its true cross-shore distance; `top_km` then sets the window.
Two panels drawn from canvases of different heights -- which is what the
real-only and with-buffers canvases are, the buffer offsets running further
seaward than any real domain -- therefore still share one y scale.
```

**`scheme_colours()`**

```text
(cmap, norm, colourbar ticks, water colour) for one elevation treatment.

`_map_panel` masks every cell below sea level and `set_bad` paints it, so
both schemes take the identical path to water; only the shade differs.
```

**`year_dir()`**

```text
`output/figures/3-model-inputs/1-domains/initial_island/<year>/<scheme>/`. `save()` creates it,
and the `supporting/` inside it, so nothing here makes a directory.
```

**`figure_absolute()`**

```text
Draw and save ONE absolute-placement map: a single panel, its own file.

The real-domain map and the with-buffers map were panels (a) and (b) of one
figure until 2026-09-17. They are separate figures now (Hannah): each is a
map of the whole reach in its own right, each is used on its own, and as a
pair they cost the page two half-height panels to say one thing twice.

Standing alone, each also gets the full column width and its OWN
cross-shore window, so neither carries the other's dead space -- the
real-domain map no longer reserves the extra kilometre the buffer offsets
need. The price is that the two no longer share metres per inch, so the
caption states the span and the exaggeration of each.
```

</details>

### study_area_figures.py

The generic site figures: where Hatteras is, how the 90 domains tile it, one domain as a grid, the forcing.

From the script's original header:

```text
The generic site figures: where Hatteras Island is, how the 90 model domains
tile it, and what one domain looks like as a Barrier3D grid. For manuscripts
and talks, drawn from the model's own inputs so a map and a run cannot
disagree.

    python scripts/figure_making/island/study_area_figures.py [--vector] [--only NAME] [--talk]

Writes into the numbered output/figures/ layout (1-site/, 2-observations/,
3-model-inputs/; see SUBJECT below), PNG at the top and PDF + CAPTIONS.md
under supporting/; --talk writes the same paths under output/figures/talk/:

    study_area.png        the reach on satellite imagery with the 90 domain
                          boxes, NC-12, the villages and structures, and a
                          regional inset
    domain_framework.png  the domains by role (boundary, scored interior,
                          community zones, buffers) over the island outline,
                          with a regional inset
    domain_framework_vertical.png
                          the same, NORTH UP as a portrait page: the reach
                          running up the page, ocean to the east
    domain_metrics.png    the island width and highest cell the 2004-start
                          extraction gives each domain
    site_overview.png     study_area and domain_framework as one two-panel
                          figure on one frame, for a manuscript
    domain_grid.png       one domain (GIS 45): its box on imagery, the 10 m
                          elevation array the model reads, the road cells, and
                          the mean cross-shore profile
    forcing_timeline_1984.png, forcing_timeline_1996.png
                          one figure per hindcast CHAIN (1984-2004 + 2004-2024,
                          and 1996-2010 + 2010-2024): its windows, storms per
                          year and the largest runup, Duck sea level with the
                          window trends, and the road and nourishment events by
                          domain and year. Each is drawn on its own span, so
                          the two do not share an axis
    observed_rates.png    CoastSat LRR per domain for the four windows
    dune_lines.png        the digitised dune line 1984-2023: per-domain
                          movement since 1984, and two detail maps
    domain_schematic.png  the cross-section and plan of one domain as the
                          coupled model builds it: shoreface, berm, dune rows,
                          interior, road setback, bay
    management_footprint.png  both NC-12 alignments, the relocations, the
                          fills and the bridge on the reach
    reach_elevation.png   both extractions at 10 m, all 90 domains

THE FRAME
    The island is drawn ROTATED a quarter turn so the reach runs left to
    right, GIS 1 (Cape Point) at the left and GIS 90 (Pea Island) at the
    right -- the same direction as every alongshore chart in the project
    (DOMAIN_AXIS_LABEL). A north-up map of a 60 km strip trending NNW is a
    thin diagonal that wastes the page. The turn is exactly 90 degrees, not
    the fitted reach axis (82.4 degrees), because the domain boxes are
    axis-aligned in UTM: at 90 degrees every box is a level rectangle on the
    page and the reach steps in y where the coast bends (Hannah, 2026-09-17).
    The north arrow points right, and is not drawn with `_north_arrow()`,
    which assumes north is up.

WHICH DOMAIN POLYGONS
    5-scr/2-transect-frame/transect_domains/HAT_domains.json: the 90 boxes, 2000 m cross-shore
    by 500 m alongshore, axis-aligned in UTM 18N, that every per-domain DEM
    clip and elevation array was cut from (the repository copy of
    D:/Hatteras_GIS/domains.geojson). NOT the retired map_elements/archive/
    domains_1000m_20251014/HAT_domains.shp: an older 1000 x 500 m set whose ID runs about
    nine domains south of the model's, found 2026-09-17 when its box for
    "45" did not contain the array for domain 45.

LAYERS AND THEIR OWNERS
    domain boxes                      5-scr/2-transect-frame/transect_domains/HAT_domains.json
    island outline                    map_elements/hatteras_outline/
    NC-12 centrelines                 hat_topo_version.road_line_file(1978|2008)
    domain elevation arrays           hat_topo_version.npy_dirs("2004-start")
    road cell masks                   hat_topo_version.road_mask_file(2008, gis)
    villages, piers, groin, zones     hatteras_site_config
    imagery                           Esri World Imagery through contextily,
                                      cached under the user's temp dir;
                                      --vector draws the outline instead and
                                      needs no network
    locator coastline                 map_elements/natural_earth/
                                      (Natural Earth 10 m states, clipped)
```

Notes that were in the code:

```text
WHICH FOLDER EACH FIGURE BELONGS TO. output/figures/ is organised by the
question a figure answers, in the paper's order (the numbered layout,
Hannah 2026-09-29; before that by subject, 2026-09-17): (subject, *parts)
as figure_dir() takes them.
```

```text
THE VECTOR BASE MAP. Water a pale blue-grey, land an ivory, edges a mid
grey: the conventional quiet base of a journal map, on which the black
road and the grey role fills carry the figure (Hannah, 2026-09-17: "more
professional and academic"). C["WATER"] stays for elevation classes, where
it means a cell at or below sea level.
```

```text
white type on imagery needs the opposite of the house halo: a dark stroke,
so a number over pale sand is as legible as one over dark water
```

```text
The boxes are axis-aligned in UTM with their alongshore side due
north, so the frame is an exact quarter turn: north to the right.
A reach axis FITTED through the centroids (82.4 deg) tilted every
box by 7.6 deg on the page (Hannah, 2026-09-17: "perfectly
horizontal"). The reach then steps in y as the coast bends, which
is the boxes' real stagger.
```

```text
Salvo, Waves and Rodanthe are single domains 5 and 6 apart -- under 3 km,
which is narrower than the names are wide on a reach-length panel, so
they printed as "Salvo WavesRodanthe". They are measured and pushed
apart, with the leader running from the moved name back to its own
domains (found 2026-09-17 when the villages were put back on the
management map).
```

```text
How far out the names actually reach, so whatever is placed beyond them
can be told to stay clear.
```

```text
x_rot = ox + (northing - oy); y_rot = oy - (easting - ox)   (theta = +90 deg)
x_rot = ox - (northing - oy); y_rot = oy + (easting - ox)   (theta = -90 deg)
```

```text
THE NAMES SIT INSIDE ONE GRATICULE CELL EACH, not across a line: the
cells are 2 deg of longitude, about half an inch here, so a name wider
than that is set on two lines (Hannah, 2026-09-17: no label text over a
border or a line).
```

```text
only the two names a reader needs: a locator carrying five at 6 pt is a
thicket at this size. Under about an inch wide even two collide, so the
state name goes and the box keeps the meaning (the vertical framework's
locator is 0.82 in).
```

```text
AT 8 PT A SMALL LOCATOR HOLDS TWO NAMES, not three: one line of
"Atlantic Ocean" is 2.5 deg wide against graticule cells of 2 deg, and
the main map names the ocean anyway. It returns on a wide locator,
clear of the reach box (which reaches 35.1 N) and of the lines.
```

```text
the latitude labels go on whichever side has open water beside them:
on the portrait map the locator's right edge is a kilometre from the
domain numbers, so they go left, over the sound
```

```text
the regional inset in the upper right, above the village row: one inch
square, the sound padded to make the room
```

```text
THE NAMES SIT ON THE OCEAN SIDE and the domain numbers on the sound
side (Hannah, 2026-09-17), which frees the whole north-west corner for
the locator: the island runs 5 km east as it goes north, so at the top
of the frame the sound is at its widest.

THE REACH SITS WELL EAST OF CENTRE, 9 km of window west of it. The
window's WIDTH is fixed by equal scale (it is the height times the
panel's shape), so the only way to widen the sound is to draw the map
shorter: `panel_in` came down from 7.8 to 6.5 in to buy those 3 km,
which costs the reach a seventh of its scale and buys the locator a
fifth more side. The east margin is what the names need and no more.
```

```text
land, clipped to the window (the outline file also holds Ocracoke and
the mainland, 20 km west of anything this figure shows)
```

```text
every tenth domain numbered on the SOUND side, the villages named on
the OCEAN side beside their own domains
```

```text
the ends and the water bodies
the end names: inside the frame with a kilometre to spare, and beside
the reach rather than on it -- placed on the axis they sat on NC-12,
which runs on past both ends of the modelled domains
```

```text
north is up here, so the ordinary arrow applies; a labelled UTM frame
replaces the scale bar
```

```text
the locator in the UPPER LEFT, the open sound at the north end
(Hannah, 2026-09-17). It clears the northernmost village name by about
a tenth of an inch, and the (a) letter goes ABOVE the frame rather than
in the corner the locator now holds.
as wide as the sound is at the north end of the frame, less the room
its own latitude labels need on the right and the domain numbers take
at 455 km easting. Every tenth of an inch more reaches a tenth of an
inch further south, where the island is further west and the sound
narrower, so this is close to the ceiling for a locator in this corner.
```

```text
THE LEGEND UNDER THE LOCATOR (Hannah, 2026-09-17), in the open sound:
the locator ends a fifth of the way down the panel and nothing else is
drawn below it until the sound's own name. Anchored to the locator's
foot rather than to a corner, so the two read as one block.
```

```text
THE TWO HINDCAST CHAINS, ONE FIGURE EACH. A start year runs a PAIR of windows
that tile at their shared year, and the 1984 pair and the 1996 pair are two
alternative readings of the same forty years, not four panels of one record:
putting all four on one frame invited a reader to compare windows that are
never run together (Hannah, 2026-09-17). Each chain is drawn on its own span
so the years it actually covers get the full width of the page, which is why
the two figures are NOT on a common axis.
```

```text
Cut to the years the chain COVERS (its first window's start to its last
window's end), not to the axis, which carries a margin year at each end:
clipping by the axis alone left a half-drawn storm bar in the margin, and
a record outside the span would still set the y limits of its panel.
```

```text
an early event's label would run off the left of the axis, so
it goes to the right of its bar instead
```

```text
the numbers on the OCEAN side, clear of the panel letter in the
other corner, and only where the whole of one clears the frame
```

```text
the arrow sits between two domain numbers, not beside one: the
numbers fall at the middle of each 500 m box, so the gap at
mid-panel is the one place in the water that is always free
```

```text
(b) the same domain in plan, on the same cross-shore axis: dune rows
then interior, alongshore up the page
```

```text
ROW ASSIGNMENT: ONE row for as long as the labels fit on it side by
side, and a second only when they no longer do. Stacking eagerly -- a
new row whenever two labels wanted the same place on the reach -- looks
tidy per label and reads badly, because the outer label's leader then
has to cross the inner label to reach its band, which is what the 1989
leader was doing to the bridge label. A reach 60 km long has room to
separate three labels sideways; use it.
```

```text
One row's worth of clear space is the tallest label plus half a line of
air, so two rows can never print into each other whatever the wording.
```

```text
The first row clears both the offset asked for and anything already
occupying the space (the village row), plus this row's own half height.
```

```text
One window for both panels. Equal aspect and a shared width mean a
shared height too, so the padding here is the worst case of the two:
village names soundward, one row of annotation seaward.
```

```text
THE PERIOD STARTS THIS FIGURE SPEAKS TO, and the digitised line each one
reads. Looked up rather than typed: there are only two NC-12 lines on
disk, 1978 and 2008, and ROAD_LINE_FOR_YEAR is what decides which start
gets which. Change PERIOD_STARTS and the key follows.

THE KEY NAMES THE START YEARS, not the tracing vintages (Hannah,
2026-09-17: "it is the same for the corresponding start years"). For
2010 that is exact -- both relocations predate the 2008 trace and the
road does not move between 2008 and 2010. For 1996 it holds everywhere
EXCEPT GIS 84-87, because the 1989 Pea Island relocation falls inside
the 1978-1996 gap; derived/1996/PROVENANCE.md is explicit that a 1996
road is the 1984 road with that one event applied, and that no line of
1996 vintage exists. The caption carries the exception, because the
same stretch is labelled "relocated 1989" right beside it.
```

```text
Opaque enough to be unmistakable over both the land tan and the water
blue, light enough that the island outline and the domain edges survive
under it. ONE VALUE PER FAMILY, each passed to both the tinted boxes and
that family's legend swatch, so the key and the map cannot disagree.
The road tint is lighter because C["ROAD"] is near-black where
C["ADDED"] is a mid orange, so equal alpha does not read as equal
weight, and the road tint carries the bridge hatch on top of it at
GIS 84-87. Matched by APPEARANCE, not by arithmetic.
```

```text
UPPER RIGHT (Hannah, 2026-09-17). The letter keeps its bold weight,
so this is two artists, and the letter's x is set from the MEASURED
width of the name: the two names differ in length and a fixed gap
would leave one pair tight and the other loose.
```

```text
---- (a) BEACH NOURISHMENT -------------------------------------------
NC-12 IS NOT DRAWN HERE. It is panel (b)'s subject, and the same blue
line appearing in (a) as unlabelled context makes a reader ask whether
it means something in (a) too (Hannah, 2026-09-17). What that costs is
the reason the Buxton and Rodanthe fills run past their villages --
they follow the road corridor -- and that is a sentence, so it lives in
the caption and the rules table rather than in a line the key cannot
label without ambiguity. The domains are numbered and the villages
named, so the footprints still read.
```

```text
TERSE, because a label's width is what makes it collide: where and
how much. The project's full name and agency are in the rules table.
```

```text
both vintages, the earlier in red, so a relocation is visible as the gap
The early line carries every relocation that has already happened by its
period start, so it IS that start's alignment rather than the raw trace.
Consequence worth knowing: the red and blue lines now coincide at
GIS 84-87, because the road did not move there between 1996 and 2010 --
the 1989 event is before both. The only divergence left is GIS 9-14,
the 1999 relocation, which is the one that does fall between them.
```

```text
The bridge span CONTAINS the 1989 relocation, so it is a hatch laid over
the tint rather than a second fill: GIS 82, 83 and 88 read as hatch
alone, 84-87 as hatch over the tint, and 9-14 as tint alone.
```

```text
---- one legend for both panels --------------------------------------
Each swatch is a domain box marked the way the panels mark one, so the
key restates the drawing rather than the colours. INK at 0.35 is
draw_reach's own box edge.
```

```text
the shear: the trend of easting against alongshore cell, fitted on the
box centres, removed column by column
```

<details><summary>Function notes (the original docstrings)</summary>

**`fig_path()`**

```text
`output/figures/<numbered subject>/.../<name>.png`, the folder created.
Under `talk/` in talk mode, mirroring the same path.
```

**`Frame()`**

```text
The rotation that lays the reach out left to right, and which way is
seaward in it.
```

**`north_arrow_rotated()`**

```text
A north arrow pointing to true north in a rotated, equal-aspect frame;
`length` is a fraction of the axes height. (x, y) is the corner the arrow
keeps clear of: the arrow is placed so that neither end crosses it.
```

**`tiles()`**

```text
Tile mosaic for a UTM window: (img, extent(left, right, bottom, top))
in `t_crs`, warped from Web Mercator.
```

**`letter_corner()`**

```text
The panel letter inside a frame but clear of it. `_letter_inside` puts
it at (0.03, 0.985), where its white backing lands on the two spines of
the corner and reads as a nick in the frame. Every lettered map panel in
this script uses this instead (Hannah, 2026-09-17: no label text on top
of a border or a line).
```

**`credit_figure()`**

```text
The tile attribution once for the whole figure, in the bottom margin.
On a panel narrower than about an inch the in-panel credit is wider than
the free water and its backing lands on the frame.
```

**`ReachArt()`**

```text
What draw_reach drew that a caller needs back: the road colour it chose
for this basemap, and how far the village names reach from the soundward
edge of the boxes, so annotations placed beyond them can clear them.
```

**`draw_reach()`**

```text
The rotated reach: imagery or the outline, the domain boxes, NC-12,
the ends, the structures. Villages are named in ONE ROW on the sound
side, each with a leader down to its domains, so they read as a labelled
axis rather than scattered text (Hannah, 2026-09-17 evening). Returns
the road colour used.
```

**`village_row()`**

```text
The five village names on one line along the sound side of the reach,
a leader from each down to the sound edge of its domains where the gap
is more than `leader_gap_m`.
```

**`utm_frame()`**

```text
UTM 18N coordinates on a quarter-turned map: northing runs along the
x axis and easting up the y axis (decreasing upward when the ocean is at
the bottom). Ticks at round kilometres, labelled in km. A labelled frame
replaces the scale bar; the north arrow stays because the frame is
turned.
```

**`framework_legend()`**

```text
One column INSIDE the map, in a corner of open water, with the faint
white backing the house style allows for a legend inside the axes
(Hannah, 2026-09-17). The ocean corner is the one that never competes
with the domains; a legend under the map cost the figure an inch of
height for nothing.
```

**`framework_handles()`**

```text
SHORT labels. What each class MEANS -- scored against CoastSat, carries
the alongshore boundary condition, roadway management off -- is in the
caption; spelling it in the legend made the longest entry five sixths of
a single-column panel, which is why the legend used to sit under the map
instead of in it.
```

**`regional_inset()`**

```text
The south-eastern US coast with the reach marked: a vector locator
from the Natural Earth 10 m states layer, in plain longitude and
latitude with a labelled graticule, so it needs no scale bar or north
arrow (the house rule for a labelled frame). Until 2026-09-17 evening
this was a grey tile basemap whose own small labels fought ours.
```

**`reach_figure()`**

```text
A double-column figure whose height is set by the window's aspect, so
an equal-aspect map panel fills it.
```

**`domain_metrics()`**

```text
Per domain, from the extraction arrays (m NAVD88, -10 water): the
median across the alongshore rows of the land width, and of each row's
highest cell. The width is capped by the 2000 m box.
```

**`fig_domain_framework_vertical()`**

```text
The domain framework NORTH UP, as a portrait figure: the reach runs up
the page, GIS 1 at the bottom and GIS 90 at the top, the ocean to the
right. Every other map in this script is turned a quarter turn, because a
45 km strip trending NNW wastes a landscape page; a portrait page is the
one shape that fits it unturned, and unturned is the orientation a reader
can carry to any other map of the Outer Banks (Hannah, 2026-09-17).

Nothing is rotated here, so the domain boxes are level by construction
(they are axis-aligned in UTM) and north is up: `_north_arrow` applies,
and the UTM frame reads the ordinary way round, easting along the bottom
and northing up the side. The window is sized from `panel_in` so the map
keeps equal scale in both directions.
```

**`fig_site_overview()`**

```text
The study-area imagery and the domain framework as one two-panel
figure on one frame, for a manuscript that should carry one map of the
reach rather than two.
```

**`storm_record()`**

```text
One row per storm, 1984-2023, with its calendar year: the 1984-2004
and 2004-2024 series tile at 2004.
```

**`line_position()`**

```text
Mean easting of each vintage's line inside each domain box, m. The
boxes are axis-aligned and the coast trends 8 degrees from north, so a
difference of eastings between vintages is the cross-shore movement to
within one per cent; positive is seaward.
```

**`processed_domain()`**

```text
The extractor's output for one domain: interior (rows cross-shore,
ocean first; columns alongshore) and the dune row, both in decametres.
```

**`half_is_lower()`**

```text
True when `half` ("sea" or "sound") is the LOWER half of a domain box
on the page. One expression, shared by the map and the legend, so a
flipped reach cannot leave the key describing the other side.
```

**`HalfBox()`**

```text
A legend entry that IS a domain box with one half marked, drawn the way
`tint_domains` marks it. Handled by HalfBoxHandler.
```

**`HalfBoxHandler()`**

```text
Draws a HalfBox: the domain box's outline at draw_reach's own weight,
with the marked half shaded inside it.
```

**`relocate_line()`**

```text
A rotated road line with a per-domain LANDWARD displacement applied.

The digitised NC-12 lines are 1978 and 2008 exports, so a line read by a
later period start is only that period's alignment if the relocations
between the two dates are applied to it. The 1989 Pea Island event falls
inside the 1978-1996 gap, which left the red line at its pre-1989 position
at GIS 84-87 while the panel labelled that same stretch "relocated 1989"
(Hannah, 2026-09-17).

The shift is per DOMAIN and therefore STEPS at a domain boundary rather
than tapering. That is the forcing, not a drafting choice: CASCADE moves
the road by a whole number of cells per domain, and a smooth taper here
would draw something the model does not do.

Args:
    rotated: the road GeoSeries already in the frame.
    displacement_m: {gis: metres landward}, as HATTERAS_ROAD_EVENTS holds it.
```

**`tint_domains()`**

```text
Tint the seaward or soundward HALF of each selected domain box.

WHY HALVES. Every management category is a statement about whole domains,
so each is drawn on the domains themselves rather than on a bar alongside
them (Hannah, 2026-09-17). At the north end they land on the SAME
domains -- the Rodanthe fill is GIS 84–89, the 1989 relocation 84–87, the
bridge 82–88 -- so one tint per whole box would stack three colours on one
rectangle and none of them would be legible. Splitting the box gives each
family its own half, and an overlap reads as both halves being marked
rather than as a fourth colour nobody can name. The split is also where
the things are: the fill goes on the ocean side, the road sits landward.

The boxes are axis-aligned rectangles in the rotated frame (see Frame),
so a half is exactly half of the geometry's bounds.

Args:
    rdom: the rotated domain GeoSeries.
    sel: boolean mask over it.
    half: "sea", "sound", or "all" for the whole box. Since the figure
        was split into one panel per family (2026-09-17) each family has
        the domain to itself again, so "all" is what the management
        panels use; the halves remain for any panel that must carry two
        families at once.
```

**`measure_m()`**

```text
Half-width and half-height of each string, in DATA metres.

Label placement has to compare text against reach distances, and a
character count cannot: it is a different physical width at every font
size. The text is rendered on a throwaway artist and converted through the
(equal-aspect) data transform instead.
```

**`separate_x()`**

```text
Push labels apart along the reach until none overlaps, in place.

`placed` is a list of dicts carrying "x" (wanted position) and "hw" (half
width), in any order. Left-to-right push first, then a right-to-left pull
back inside the window, so a label cannot be shoved off the canvas by the
one before it, and a margin keeps any of them off the frame.
```

**`label_lanes()`**

```text
Annotation labels beside the reach, on as few straight rows as they fit.

WHY THIS EXISTS. Every label here used to carry its own hand-tuned
perpendicular offset (3300, 4500, 6500 ...), each tuned once against one
rendering. Two labels whose bands are close in the reach then overlap,
and nothing in the figure prevents it -- the 1989 relocation and the 2022
bridge annotations were printing on top of each other, and the Buxton
fill label was running under the domain boxes.

Labels in a row share one perpendicular offset, so the annotations read
as a labelled axis rather than as scattered text -- the convention
village_row() already uses for the village names (Hannah, 2026-09-17
evening). A label goes in the first row where it clears the labels
already there; within a row the labels are then pushed apart along the
reach until none can touch, and a leader runs from each to the band it
names. Widths are MEASURED from the rendered text, so this holds at any
font size and for any wording.

Args:
    items: [(gis_mid, anchor_xy, text)] -- the reach position the label
        belongs to, the point on its band to lead to, and the text.
    side: +1 to place seaward of the reach, -1 soundward.
    y_edge: the domain boxes' edge on that side, in frame coordinates.
    first_offset_m: perpendicular distance from `y_edge` to the first row.
    clear_m: a distance from `y_edge` already occupied -- the village
        names, say. The first row is pushed outside it rather than
        printing on top of it.
    clear_pad_m: the air left between that and the first row.
    row_step_m: distance between rows. Default: the tallest label plus a
        half line, MEASURED -- a literal step here is what let the 1989
        and bridge labels sit one line apart and touch.
    pad_m: the clear space kept between two labels in a row.
    leader_gap_m: draw a leader only when label and band are further
        apart than this.

Returns:
    The outermost y the labels reach, so the caller can keep the legend,
    the scale bar and the north arrow clear of them.
```

**`fig_management_footprint()`**

```text
Two panels on one reach: (a) what was added to the beach, (b) what was
done to the road.

ONE MAP CARRIED BOTH until 2026-09-17, and the two families competed for
the same domains and the same label space. At the north end the Rodanthe
fill (GIS 84-89), the 1989 relocation (84-87) and the bridge span (82-88)
all cover the same ground, so a whole-box tint per family would have
stacked three colours on one rectangle; that is what forced each family
onto half a box, and it still left six annotations and five village names
competing for two label rows. Split into panels, each family has the WHOLE
domain again -- which is what the model applies it to -- and each panel
carries two or three labels. The panels share ONE window, so the same
domain is the same place on the page and a reader can still see that the
Rodanthe fill and the bridge cover the same ground (Hannah, 2026-09-17).

The two panels are deliberately identical in structure: villages named
soundward, annotations seaward, so the eye learns the layout once.
```

**`reach_mosaic()`**

```text
Domains lo..hi of one extraction as one (cross-shore, alongshore)
mosaic, m NAVD88, water -10, ocean at the top, cropped to the rows that
hold land plus a margin.

Each box has its own easting, so stacking the arrays edge to edge (as
until 2026-09-17) drew the coast as a sawtooth. The boxes are first laid
at their true eastings, which makes the coast continuous; then the
LINEAR alongshore trend of the box eastings across the strip is fitted
and every alongshore column is shifted by it -- the shear the extractor's
own `straighten` applies -- so the coast runs level and only its residual
bend remains. Cross-shore distances within a column are unchanged. Drawn
at true eastings the bend made each strip 3-4 km tall and the figure
overran the page. Returns the mosaic and its cross-shore extent
(m, in the sheared frame: top, bottom).
```

**`talk_mode()`**

```text
A projector reads pale fills as white and thin lines as nothing, and a
slide is seen from further away than a page. So: the water barely off
the white of the slide (the land still ivory), type one point larger,
hairlines heavier, and the output kept apart under talk/.
```

**`anchor_at()`**

```text
Where a label's leader stops: the seaward edge of the tinted block,
as a SPAN rather than a midpoint.

Anchoring every leader at its block's centre made the two north-end
labels converge into a V, because the bridge span (GIS 82-88) and the
1989 relocation (84-87) have almost the same centre. Given the span,
label_lanes attaches each leader to the nearest part of the block it
names, so the two arrive at opposite ends and stay apart (Hannah,
2026-09-17).
```

**`panel()`**

```text
The shared basemap. The scale bar and the north arrow are drawn
once, on (b): the panels are the same map at the same scale, and a
second set would say they might not be.
```

</details>

### visualize_domain_topography.py

One CASCADE domain's topography as a cross-shore vs alongshore heatmap, dune and interior labelled.

From the script's original header:

```text
Visualize CASCADE Domain Topography
Creates a cross-shore vs alongshore elevation heatmap for a single CASCADE domain,
with labeled dune domain and interior domain regions.

Usage:
    python visualize_domain_topography.py

Date: January 2025
```

Notes that were in the code:

```text
Common options: 300 (nearshore focus), 500 (moderate), 1000 (wide view)
Set CROP_CROSSSHORE = False to use CASCADE's evolved island width
```

```text
Get domain topography at specified time step
DomainTS is in units of dam (decameters), so multiply by 10 to get meters
```

```text
Get dune domain at specified time step
DuneDomain is also in dam, convert to meters and add berm elevation
```

```text
Combine dune and interior domains
Stack dune domain (first few rows) with interior domain
```

```text
Calculate figure dimensions to maintain proper aspect ratio
aspect_ratio = physical_width / physical_height
```

```text
Add region labels
Dune domain label (near bottom)
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_elevation_data()`**

```text
Load elevation data for a specific domain from CASCADE comparison.

Parameters:
-----------
filepath : str
    Path to NPZ file containing CASCADE comparison
domain_idx : int
    Domain index to extract (0-indexed)
time_step : int
    Time step to extract. Default -1 gets the final time step.
max_crossshore_m : float or None
    Maximum cross-shore distance to include (m). If None, includes full domain.
dy : float
    Cross-shore resolution (m) for calculating crop index
    
Returns:
--------
elevation : ndarray
    2D array of elevation values (cross-shore x alongshore) in meters MSL
```

**`create_topography_plot()`**

```text
Create a heatmap visualization of domain topography.

Parameters:
-----------
elevation : ndarray
    2D array of elevation values (cross-shore x alongshore)
dy : float
    Cross-shore resolution (m)
dx : float
    Alongshore resolution (m)
dune_width : int
    Number of cross-shore cells in dune domain
domain_num : int
    Domain number for title
z_lim : float
    Maximum elevation for colorbar
    
Returns:
--------
fig, ax : matplotlib figure and axis objects
```

**`print_domain_statistics()`**

```text
Print useful statistics about the domain topography.

Parameters:
-----------
elevation : ndarray
    2D array of elevation values
dy : float
    Cross-shore resolution (m)
dx : float
    Alongshore resolution (m)
dune_width : int
    Number of cross-shore cells in dune domain
```

</details>
