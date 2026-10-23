# road_relocation

How far NC-12 moved between two digitised vintages, per Barrier3D domain.

Port of `from_roya/road_relocation_dis.py` (was `scripts/input_prep/roya_files/` until 2026-09-18) (Pea Island, 1992→1996) onto
Hatteras. Same measurement — sample the old road inside each domain, measure
each sample to the whole new road — with the Hatteras road files, the
90-polygon domain file, and a signed direction added.

```
HAT_road_relocation_distance.py
```

```
python HAT_road_relocation_distance.py
HAT_RELOC_FROM=1978 HAT_RELOC_TO=2008 python HAT_road_relocation_distance.py
```

Writes to `data/hatteras_init/4-mgmt-forcing/road_relocation/<from>_<to>/`:
a per-domain CSV, the sample points as GeoJSON, and three figures (below).
`<from>` and `<to>` are LINE vintages -- the lines on disk are the 1978 and
2008 exports, filed under those years since 2026-09-15, so the folder is
`1978_2008/`. Which hindcast period reads which line is
`hat_topo_version.ROAD_LINE_FOR_YEAR` (1984 → 1978, 2004 → 2008).

**This is not a forcing.** `road_setback` comes from `../road_offset/`,
measured against interior row 0 of the model grid. This is a GIS-frame
observation of two digitised lines, and nothing reads it. Its use is to say
whether an observed relocation happened and where, so a *modelled* relocation
can be checked against one.

## The zeros are not measurements

**70.5% of the 1978 line's vertices are identical to a 2008 vertex to the
millimetre.** The 2008 line was digitised by editing a copy of the 1978 one,
and only the stretches that visibly moved were re-drawn. Every shared vertex
measures exactly 0.000 m by construction.

So a 0 m domain here means *nobody edited that stretch* — not that the road
held still. The script measures the shared fraction itself and prints it
before any result.

A second artefact sits behind the first: some stretches *were* re-drawn, but
only **re-traced** — the new centreline wanders a metre or two and never
leaves the old road's own footprint. That is one road digitised twice, not a
road that moved.

Every domain therefore carries a `classification`:

| | |
|---|---|
| `no_edit` | `coincident_fraction ≥ 0.90` — the line was copied through unedited |
| `redigitized` | largest displacement anywhere `< REDIGITIZE_MAX_M` (5 m) |
| `relocated` | the road actually moved — **the only measured subset** |

Nothing is dropped from the CSV; only `relocated` domains are summarised and
drawn. **Never average across the whole table** — it averages real movement
against copy and re-tracing alike, which is why no island-wide mean is printed.

### Why 5 m

About half a road width: NC-12 is ~8 m wide, and two digitisings of one road
off different photos disagree by that much from georeferencing alone. The
Hatteras data splits cleanly there and the exact value does not matter —
re-traced domains top out at **3.25 m**, real relocations start at **12.2 m**,
so anything from ~3.5 to ~12 m gives the same answer.

Direction corroborates it independently. The re-traced domains are internally
coherent (`sign_agreement = 1.00`) but flip sign arbitrarily between
neighbours — 16–18 seaward, 19–22 landward, 24–27 seaward. That is what a
smoothed line does. A road that moved does not change its mind every few
hundred metres.

## What the columns mean

| Column | |
|---|---|
| `mean_relocation_m` | unsigned nearest-distance, old road → new road. Magnitude only, and a *lower bound* on cross-shore movement: where the lines cross obliquely the nearest point is diagonally alongshore |
| `mean_signed_landward_m` | same displacements, signed **+ landward / − seaward**, from which side of the northward-oriented old road the new road falls on (`OCEAN_ON_RIGHT`) |
| `coincident_fraction` | share of samples on unedited shared geometry |
| `sign_agreement` | share of *moved* samples agreeing on direction; `ROADS_CROSS` below 0.90 means the domain mean is averaging both ways |
| `oblique_fraction` | share of samples where the road runs >60° off north; over half and the domain is flagged `OBLIQUE_SIGN` |
| `classification` | `no_edit` / `redigitized` / `relocated` — see above |

The sign has one blind spot. Around the Cape Point bend the island turns
east–west, the ocean stops being on the right of a northward road, and the
landward/seaward convention breaks. Those domains are flagged `OBLIQUE_SIGN` —
on the 1978→2008 pair that is **domain 8 and nothing else**. Magnitude is
unaffected; it never used the tangent.

## The three figures

They answer three questions at three incompatible scales, so they are three
files rather than one page.

| File | Question |
|---|---|
| `..._alongshore.png` | **Where** along the island did the road move? |
| `..._sites.png` | **What** did the move look like? |
| `..._domain_map.png` | **Which** domains carry a relocation? |

Since 2026-09-10 the three are drawn in the house style of the dune-line
figures (`apply_style()` in `0-elevation/3-figures/HAT_plot_duneline_offset.py`:
Arial, panel letters, RdBu poles with the earlier vintage red and the later
blue, north arrow and scale bar on tickless maps, frameless legends outside the
axes) and carry no title sentences or statistics lines. That text, with the
numbers filled from the CSV, is in `CAPTIONS.md` beside the figures, written by
the same run. The `sites` figure labels each domain with its measured
displacement and the value CASCADE is forced with, the measurement rounded to
the nearest 10 m cell (`hatteras_site_config.py`, ROUNDED TO WHOLE CELLS).

The island is 8 km wide and 45 km long — aspect 0.18 — so any true-scale map of
the whole thing is a hairline in a column of white space. `sites` sidesteps
that by only ever drawing ~2 km; `domain_map` rotates the island 90° clockwise
so it runs south (left) → north (right), matching the domain axis of
`alongshore`. Rotation preserves distance, so its scale bar is still honest.

**alongshore** — signed bars per domain, with a pale bar behind for the largest
displacement anywhere in that domain (the gap between them is how much of the
domain moved). Artefact domains shaded, town spans along the foot so a domain
number means a place, each site bracketed and named.

**sites** — one true-scale zoom per relocation site. Sites are found from the
data (contiguous runs of `relocated`), not named in the script, so this still
works on a different pair of vintages. **All panels share one scale**: Buxton
spans 8 domains and Rodanthe 4, and framing each on its own extent would render
them at different metres-per-inch, making the shorter site look like the bigger
relocation. The old road is a dashed spine drawn **on top of** its own sample
points — underneath, the points hide it and the two alignments no longer
visibly pull apart — plus an OCEAN label, because the signed metric depends
entirely on which side that is.

**domain_map** — every domain in its real place, filled by what it carries:
white where the road never reaches, light grey never-edited, darker grey
re-traced, and the relocated ones outlined in black and filled by distance on
the same colour scale as `sites`.

## Result, 1978 → 2008 (the 1984 and 2004 periods' roads)

Of 83 road-carrying domains: **56 never edited, 15 re-traced, 12 relocated.**

The 12 are two coherent stretches, **every one of them landward**, with
nothing at all in between:

| Domains | Where | Largest domain mean |
|---|---|---|
| 8–15 | Buxton | **77 m landward** (domain 11) |
| 84–87 | Rodanthe | **109 m landward** (domain 85) |

Across those 12: mean 45.1 m, median 40.1 m, 0 seaward. Domains 8 and 15 are
flagged `ROADS_CROSS` — they are the tapered ends of the Buxton relocation,
where the new alignment rejoins the old and the domain mean mixes both
directions. Domain 8 also carries `OBLIQUE_SIGN` (Cape Point bend), so read its
magnitude and ignore its direction.

Both stretches sit where the hindcast already expects trouble; Rodanthe is the
S-curves. Treat the magnitudes as the movement of a *digitised centreline*
between a 1978 and a 2008 photo, which is the interval those two files really
span — the 1984/2004 labels are the hindcast periods they stand in for.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_road_relocation_distance.py

How far NC-12 moved between two digitised vintages, per Barrier3D domain.

From the script's original header:

```text
How far NC-12 moved between two digitised vintages, per Barrier3D domain.

Hatteras port of from_roya/road_relocation_dis.py (was roya_files/ until 2026-09-18). Same measurement: sample
the OLD road inside each domain, measure each sample point to the WHOLE new
road, and summarise per domain. Everything Hatteras-specific -- the road
vintages, the 90-polygon domain file, the CRS chain -- is in CONFIG below.

    python HAT_road_relocation_distance.py
    HAT_RELOC_FROM=1978 HAT_RELOC_TO=2008 python HAT_road_relocation_distance.py

WHAT THE NUMBER IS
`mean_relocation_m` is an UNSIGNED nearest-distance from the old road to the
new one. It says how far the road moved, not which way, and because it takes
the *nearest* point on the new road it is a lower bound on cross-shore
movement: where the two roads cross obliquely the nearest point is diagonally
alongshore, not straight across.

So this script also reports a SIGNED column, `mean_signed_landward_m`, built
from the same displacement vectors:

  landward = the new road is LEFT of the old road heading south -> north
           = west, i.e. away from the ocean         (see OCEAN_ON_RIGHT)

Positive is landward, negative is seaward. Read the signed column when you want
direction and the unsigned one when you want magnitude; `sign_agreement` tells
you what fraction of a domain's samples agreed on the direction, and anything
below ~0.9 means the two roads cross inside that domain and the domain mean is
averaging two directions.

That convention has one place it breaks. Around the Cape Point bend the island
turns east-west, the ocean is no longer on the right of a northward road, and
the sign becomes meaningless -- so domains where the road runs more than 60
degrees off north are flagged OBLIQUE_SIGN. On the 1978-2008 pair that is
domain 8 and nothing else. Magnitude is unaffected: it never used the tangent.

READ THIS BEFORE READING THE ZEROS
The two Hatteras road files SHARE MOST OF THEIR GEOMETRY. 558 of the 1978
line's 791 vertices are identical to a 2008 vertex to the millimetre: the 2008
line was digitised by editing a copy of the 1978 one, and only the stretches
that visibly moved were re-drawn. Roughly 72% of the old line is therefore
exactly 0.000 m from the new line by construction.

A 0.00 m domain in this table means NOBODY EDITED THAT STRETCH. It is not
evidence the road held still. Real movement can only be claimed where the
lines actually diverge, so the script measures the shared fraction itself,
reports it per domain as `coincident_fraction`, and flags those domains
NO_EDIT rather than letting them read as a measurement.

A second artefact sits behind the first: some stretches WERE re-drawn, but
only re-traced -- the new centreline wanders a metre or two and never leaves
the old road's own footprint. That is one road digitised twice, not a road
that moved, so domains whose largest displacement anywhere stays under
REDIGITIZE_MAX_M are flagged REDIGITIZED.

Every domain therefore carries a `classification`:

    no_edit      the line was copied through unedited
    redigitized  re-traced within the road's own width
    relocated    the road actually moved  <- the only measured subset

Nothing is dropped from the CSV, but only `relocated` domains are summarised
and drawn. Any island-wide average over the full table is meaningless -- it
averages real movement against copy and re-tracing alike.

WHAT THIS IS, AND IS NOT
Not a setback. `road_setback` comes from ../road_offset/, measured against
interior row 0 of the model grid; this is a GIS-frame observation of how far
apart two digitised lines lie.

It IS a CASCADE forcing as of 2026-08-20. `HATTERAS_ROAD_EVENTS` in
scripts/site_layer/hatteras_site_config.py reads `mean_signed_landward_m` out of the
1978_2008 CSV for the domains its two historical relocation events move -- GIS
84-87 (1989) and GIS 9-15 (1999). It replaced eleven hand-entered literals
attributed to a 1978->1997 ArcGIS measurement whose 1997 line is not in the
repo. So re-running this script with those vintages CHANGES WHAT THE MODEL IS
FORCED WITH; it is no longer a diagnostic that can be regenerated freely. The
config refuses any domain this script classifies `no_edit` or `redigitized`,
so a re-run that reclassifies one of those eleven raises at import rather than
quietly forcing a number the lines do not support.

A NOTE ON THE VINTAGES
nc12_1978.geojson and nc12_2008.geojson were digitised off 1978 and 2008
imagery -- the nearest usable coverage to the two hindcast period starts. The
labels are the periods they stand in for, not the photo dates, so the interval
measured here is really ~30 years, not 20.

REQUIREMENTS
  geopandas, shapely, numpy, pandas, matplotlib
```

Notes that were in the code:

```text
The line files and the output folder resolve through hat_topo_version.py
(2026-09-18); scripts/ goes on the path here because the file's own
sys.path setup comes later.
```

```text
Which pair of road vintages to compare. Override from the shell to compare a
different pair without editing the file.
LINE vintages, not period starts (2026-09-15): the lines are 1978 and 2008
exports, filed under those years, and the output folder is named by them --
road_relocation/1978_2008/. Which period reads which line is
hat_topo_version.ROAD_LINE_FOR_YEAR.
```

```text
The same 90-polygon domain file the shoreline and groin work uses. It left
the scripts tree with the observed shoreline data (commit 17a0334f) and
lives under data/ now; the old scripts/input_prep/5-scr/CoastSat/ path was
still here until 2026-09-15.
Resolved through hat_observed_rates.py since 2026-09-18.
```

```text
Three figures, each answering one question, rather than one page trying to
answer all three at incompatible scales.
```

```text
Domains to analyse. None means every domain the old road touches.
Barrier3D runs GIS 9-90; the road file reaches 8-90.
```

```text
NAD83(2011) / UTM Zone 18N. The road files are EPSG:2264 (NC State Plane,
US survey feet) and the domains EPSG:3725, so everything is reprojected here
and every distance below is metres.
```

```text
NC-12 runs south -> north up Hatteras with the Atlantic to the east, so the
ocean is on the RIGHT of the direction of travel and landward is LEFT. This
is what turns an unsigned distance into a signed one; flip it for a site
where the ocean sits on the other hand.
```

```text
Town and village spans, used only to name the relocation sites on the
figure. Taken from hatteras_site_config rather than restated here, so a
domain number means the same place in this figure as in every other one.
```

```text
The domains the model's two historical events actually move. The site
figure labels each of these with the displacement it is forced with, so
a reader sees the number that reaches CASCADE, not only the colour; the
domains a site spans but the events leave alone (GIS 8, 15) say so.
```

```text
The forcing is the measurement rounded to the nearest cell (see
hatteras_site_config, ROUNDED TO WHOLE CELLS); label both.
```

```text
Below this fraction of samples agreeing on direction, the two roads cross
inside the domain and the signed mean is mixing landward with seaward.
```

```text
A sample closer than this to the new road is sitting on SHARED geometry --
a vertex the digitiser never edited -- not on a measured non-movement. Well
below any plausible digitising precision, so nothing real is discarded.
These samples carry no direction and are excluded from the signed statistics.
```

```text
Above this coincident fraction, the domain is unedited copy rather than a
measurement, and is flagged NO_EDIT.
```

```text
A domain whose LARGEST displacement anywhere still falls below this was
re-traced, not relocated: the new centreline never leaves the old road's own
footprint. NC-12 is ~8 m wide, so 5 m is about half a road width -- two
digitisings of one road off different photos disagree by that much from
georeferencing alone.

The Hatteras data splits cleanly here and the exact value does not matter:
the re-traced domains top out at 3.3 m and the real relocations start at
12.2 m, so anything from ~3.5 to ~12 m gives the same answer. Corroborated
by direction -- re-traced domains are internally coherent but flip sign
arbitrarily between neighbours (16-18 seaward, 19-22 landward, 24-27
seaward), which is what a smoothed line does and a moved road does not.
```

```text
Northward component of the road's local tangent below which the landward /
seaward convention stops being safe: 0.5 is a bearing more than 60 degrees
off north. Around the Cape Point bend the island turns east-west and the
ocean is no longer on the right of a northward road, so those domains are
flagged OBLIQUE_SIGN and their sign should not be trusted. Magnitude is
unaffected -- it never depended on the tangent.
```

```text
Orient northward so "left" means the same thing everywhere, whichever
way the digitised line happens to run.
```

```text
Each vintage may be one feature or many; combine into a single geometry so
the nearest-point search sees the whole road.
```

```text
Run before anything is measured, because it decides how the whole table
should be read. Two lines that share vertices were not digitised
independently, and every shared vertex contributes a 0.000 m "relocation"
that is an artefact of the editing workflow.
```

```text
Intersect only the OLD road with the current domain. This determines
which portion of the old road belongs to this domain.
```

```text
Samples sitting on shared, unedited geometry carry no direction: their
sign is numerical noise on a sub-millimetre displacement. Keep them out
of every signed statistic, and count them instead.
```

```text
Where the road runs east-west rather than north-south -- around the
Cape Point bend -- "ocean on the right of northward travel" stops
holding, and the sign of the displacement cannot be trusted. The
northward component of the local tangent is what detects it.
```

```text
Three-way classification. Only `relocated` domains carry a measurement
of the road moving; the other two are artefacts of how the lines were
drawn, and are kept in the table but out of every summary.
```

```text
Three figures, because they work at three incompatible scales and one page
cannot serve all of them:

1. alongshore   WHERE along the island the road moved      (domain axis)
2. sites        WHAT the move looked like                  (~2 km, true scale)
3. domain map   WHICH domains carry a relocation           (45 km, true scale)

Hatteras is 8 km wide and 45 km long -- aspect 0.18 -- so any true-scale map
of the whole island is a hairline in a column of white space. Figure 3 gets
round that by rotating the island to run left-right; figure 2 sidesteps it by
only ever drawing ~2 km at a time.

HOUSE STYLE (Hannah, 2026-09-10: "more professional and academic styled").
The dune-line figures' style block is reused, not restated: Arial 8-10 pt,
thin dark-grey axes, the ColorBrewer RdBu poles (earlier vintage red, later
blue), panel letters, north arrow and scale bar on tickless maps, frameless
legends outside the axes, and NO title sentences, statistics lines or
footnotes on the canvas. That text is in CAPTIONS.md beside the figures,
written at the end of this section with the numbers filled from the table,
so it cannot go stale against the picture without the CSV going stale too.
```

```text
Legend wording: no working vocabulary (Hannah, 2026-09-08). A reader outside
the project has to be able to take each entry literally.
```

```text
One colour scale shared by figures 2 and 3, so a colour means the same
distance in both. Viridis: sequential, print- and CVD-safe.
```

```text
Shade the domains that carry no measurement of the road moving. Both stay
lighter than the data; re-traced is the darker of the two because it is the
one a reader is likelier to mistake for a small relocation.
```

```text
Pale bar = the largest displacement anywhere in the domain, solid bar = the
domain mean. The gap between them is how much of the domain actually moved.
Landward carries the later vintage's blue, seaward the earlier one's red,
so the sign reads off the colour before the axis is consulted.
```

```text
The villages as the house bands behind the axis, named once, so a domain
number means a place. `town_bands` reads the same spans from the site config.
```

```text
A STRIP, not a wash: this panel already shades two data classes
("centreline unchanged", "re-digitised") full height in grey, and a
third full-height grey for the villages cannot be told from them.
```

```text
ONE SCALE ACROSS ALL PANELS. Buxton spans 8 domains and Rodanthe 4, so
framing each on its own extent renders them at different metres-per-inch and
the shorter site looks like the bigger relocation. Every panel gets the same
vertical span (the tallest site, padded); its width follows its own site, so
a centimetre means the same distance in both and no panel is mostly white.
```

```text
Drawn at the printed width: the panels share what the colour bar and the
margins leave of a double column, and their common height follows.
```

```text
One domain of context on each side, so the road is seen entering and
leaving the relocation rather than starting at the panel edge.
```

```text
The old road goes ON TOP of its own sample points, as a thin dashed
spine. Drawn underneath it is invisible -- the samples sit exactly on it
-- and the reader loses the one thing this panel exists to show: the two
alignments pulling apart.
```

```text
Domain number and, under it, the measured displacement beside the value
CASCADE is forced with (that mean rounded to the nearest cell; see
hatteras_site_config, ROUNDED TO WHOLE CELLS). Anchored low in the left
of each box: the road runs up the right side of every site box, and
crosses mid-height only in the Buxton loop (GIS 8), so the lower-left
corner is clear in all of them. A domain the site spans but no event
moves says so, rather than being read as a prescribed move from its
colour.
```

```text
A box whose foot sits in the panel's bottom band shares it with the
scale bar; that one is labelled from its top-left corner instead.
```

```text
The whole signed metric rests on which side the ocean is, so the map has
to say. North-up, so the ocean side of a northward road is the right.
```

```text
Every domain in the file, in its real place, coloured by what it carries.

Drawn north-up this is 8 km across and 45 km tall and nothing is legible, so
the island is rotated 90 degrees clockwise to run south (left) -> north
(right). That puts the ocean at the bottom and matches the domain axis of
figure 1, so the two figures read the same way round. Rotation preserves
distance, so the scale bar is still honest.
```

```text
Carry the classification onto the domain polygons. Domains the road never
reaches get their own category rather than being lumped in with "no edit":
there is nothing there to have edited.
```

```text
The highlight: relocated domains filled by how far the road moved, so the
figure says which domains AND how much in one read.
```

```text
Label every tenth domain along the strip, plus the ENDS of each relocation
run: at the printed width, numbering every relocated domain runs the labels
into each other, and the bracket above already spans the run.
```

```text
Name each relocation site above the strip it belongs to. The statistics
are in the caption.
```

```text
Which way is which, now that the map is rotated off north: the ends named,
the ocean named, and a north arrow pointing along the strip.
```

```text
One caption per figure, numbers filled from the table this run wrote, so a
document takes its caption from here and the picture stays a picture.
```

<details><summary>Function notes (the original docstrings)</summary>

**`sample_line_geometry()`**

```text
Generate regularly spaced sample points along line geometry.

Returns (point, tangent) pairs. The tangent is the local direction of the
sampled line, oriented northward, and is what gives the distance a sign.
```

**`local_tangent()`**

```text
Unit direction of `line` at `distance` along it, oriented so it points
north. Returns None where the tangent cannot be resolved.
```

**`signed_relocation()`**

```text
Distance from `point` to `target_geometry`, signed positive landward.

The sign is the side of the old road the new road falls on: the z of the
cross product of the northward tangent with the displacement vector is
positive to the left, and left is landward when the ocean is on the right.
```

**`contiguous_runs()`**

```text
Group sorted domain numbers into runs, allowing gaps up to `max_gap` so a
single unmeasured domain does not split one relocation into two sites.
```

**`site_label()`**

```text
Name a run of domains after the place it overlaps.

Falling back to the nearest village matters here: the northern relocation
sits at domains 84-87, past the end of every span in the site config, and
would otherwise go unnamed on the figure. It is the stretch north of
Rodanthe -- so say that, rather than claim it IS Rodanthe.
```

**`add_scale_bar()`**

```text
A capped bar in data units, so it scales with the panel and cannot
disagree with it. White halo so it reads on any backdrop.
```

</details>
