# 2-brie-offset - where each domain sits cross-shore at model year zero

The island offset BRIE places the domains by: a digitised dune line (or the
CoastSat shoreline) intersected with the 100 m transects, turned into a
relative offset per domain and padded to the model's 120 domains. Builds are
versioned under `data/hatteras_init/2-brie-offset/<year>/<source>/`, and the
runner reads the one named by `CURRENT`.

```
1-produce/
    build_island_offset.py        the driver: a dune line in, a padded model input out
    duneline_to_raw_offsets.py    dune line -> raw per-transect offsets (the ArcGIS step, in shapely)
    island_offset_hybrid.py       raw offsets -> relative offset per domain, padded to 120
    colleague_old_version/        a colleague's original scripts, kept as they were (not restyled)
2-figures/
    compare_offset_sources.py     dune line against shoreline for one start year
    HAT_compare_offset_versions.py  two builds of one start year
    HAT_plot_offset_profile_1to1.py the padded profile at true scale
```

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 1-produce/build_island_offset.py

A digitised dune line in, a model input out: the padded 120-domain island offset for one start year.

From the script's original header:

```text
build_island_offset.py -- a digitised dune line in, a model input out
One command from a geojson under 2-brie-offset/dunelines/ to
the padded 120-domain file a hindcast start reads, following the same steps
that were run by hand until 2026-09-15 (Hannah: "one driver that chains the
two steps"):

    1. duneline_to_raw_offsets.py   line x 100 m transects -> per-transect
                                    stations from the offshore datum, written
                                    to raw_offsets/<vintage>_duneline_offset_raw.csv
    2. island_offset_hybrid.py      first row per transect, domain mean, zeroed
                                    on the minimum, padded to 120 with the
                                    model's smooth wrap-around -> <year>/v<n>/
    3. HAT_compare_offset_versions  the new build against the previous one,
                                    in both frames (model, and fixed datum)
    4. CURRENT <- v<n>              unless --no-current
    5. PROVENANCE.md                written from the run, not by hand

The two scripts are unchanged and still run on their own.

VINTAGE VS PERIOD YEAR
    The geojson is named for the IMAGERY vintage of the line (duneline_1997);
    --year is the hindcast start it serves (1996). The pairing is
    hat_topo_version.DUNE_LINE_FOR_YEAR, the only place it is spelled, and
    this driver refuses a pair the table does not hold rather than guess.

VERSIONS
    Every build is a version: <year>/v1/, v2/, ... and a CURRENT file naming
    the one the runner reads (hatteras_site_config._island_offset_file).
    A version is never overwritten; ask for the next one. Each version folder
    keeps a copy of the raw file it was built from, so it can be rebuilt with
    island_offset_hybrid.py --raw-file after the vintage's raw has moved on.

USAGE
    python build_island_offset.py --duneline duneline_1984.geojson --year 1984
    python build_island_offset.py --duneline duneline_1997.geojson --year 1996 \
        --compare-with v1
    python build_island_offset.py --duneline duneline_2009.geojson --year 2010   # once in the table

    --version vN       name the version (default: the next free one)
    --compare-with X   a version folder under <year>/ to compare against
                       (default: the previous version, else a superseded_*
                       folder if one exists, else no comparison)
    --validate-against a GIS export of the same line, passed to step 1
    --no-current       build but leave CURRENT as it is

EXTENDED GEOMETRY (2026-09-16, the Pea Island extension experiment)
    python build_island_offset.py --duneline duneline_1997.geojson --year 1996 \
        --geometry n115

    Runs step 1 with --extension (the transects beyond GIS 1-90, numbered on
    the continued grid, to raw_offsets/ext/) and step 2 with --geometry (the
    surveyed raw plus the extension, zeroed on the surveyed minimum, padded,
    to <year>/ext/<geometry>/). Not a version: CURRENT is untouched, no
    comparison is drawn, and the provenance lands in the ext folder.
```

Notes that were in the code:

```text
geojson properties copied into the provenance when present; imagery_date is
the one every new line should carry (the 1984 and 2004 dates had to be
recovered from a metadata file and from memory).
```

```text
<year>/duneline/ since 2026-09-22: this driver only ever builds from a
dune line, so it writes into that source's folder, not the start's.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_echo()`**

```text
Print a step's output on whatever console this is (Windows cp1252
included) without dying on a character it cannot show.
```

**`build_extension()`**

```text
The extended-reach build: steps 1 and 2 in their extension modes,
then a provenance file in <year>/ext/<geometry>/.
```

</details>

### 1-produce/duneline_to_raw_offsets.py

Dune line -> raw per-transect offsets: the ArcGIS step, redone in shapely.

From the script's original header:

```text
Dune line -> raw per-transect offsets (the step ArcGIS used to do)

Every file in data/hatteras_init/2-brie-offset/raw_offsets/ was, until
2026-09-15, an ArcGIS export: buffer the digitised dune line by 1.5 m, intersect
the buffer with 1 m points generated along the 100 m transects, and export the
attribute table. The distance that matters is ORIG_LEN, the station of the
point along its transect measured from the transect's start on the offshore
datum line. island_offset_hybrid.py then keeps one point per transect, averages
the ~5 transects in each 500 m domain, and pads for CASCADE.

This script does the same intersection in shapely, so a re-digitised line can
be turned into a raw file from the repo alone, with no GIS session and no
external drive. One row per transect, holding the exact intersection station.

VALIDATED 2026-09-15 against the ArcGIS export of the v1 1997 line
(raw_offsets/1997_duneline_offset_raw.csv): all 450 transects in GIS 1-90
match, mean difference -1.01 m, sd 0.32 m, worst 1.9 m. The constant metre is
the GIS convention, not geometry: of the ~3 one-metre stations inside the
1.5 m buffer the export lists the LANDWARD-most first, and the downstream
scripts take the first row per transect. The exact intersection sits ~1 m
seaward of it. This cancels in island_offset_hybrid.py (each year is zeroed
on its own minimum) and cancels in an end-year difference of two files built
by THIS script; a difference between a GIS-built and a shapely-built file
carries the metre. Pass --validate-against to reproduce those numbers.

WHICH CROSSING when a transect meets the line more than once: the landward-most
(largest station), which is what the GIS first-row convention returned. The
count is written to n_crossings so those transects can be found.

USAGE
    python duneline_to_raw_offsets.py --duneline duneline_1997.geojson         --out 1997_v2_duneline_offset_raw.csv         --validate-against 1997_duneline_offset_raw.csv

    --duneline   a file under 2-brie-offset/dunelines/, or a path
    --out        a file name under 2-brie-offset/raw_offsets/, or a path

EXTENSION MODE (2026-09-16, the Pea Island extension experiment)
    python duneline_to_raw_offsets.py --duneline duneline_1997.geojson --extension

    The 172 transects the surveyed polygon join left without a domain -- Pea
    Island north of GIS 90, and the last kilometre south of GIS 1 -- are given
    one by the SAME line-intersects-polygon join onto Hannah's whole-island
    polygons (hat_extension_domains.join_lines), and intersected the same
    way. A transect no polygon covers (the kilometre south of GIS 1, and the
    slivers between polygons) is dropped, as the surveyed join dropped it. Written
    to raw_offsets/ext/<vintage>_duneline_offset_raw_ext.csv, the same columns
    as the surveyed file, and the transect-to-domain table once to
    transects/transects_100m_ext.csv. The surveyed file is not touched.
```

Notes that were in the code:

```text
The 100 m transects, 10 km long, each starting on the offshore datum line
(x = 460198 in EPSG:3725) and running west across the island. domain_id is
the ArcGIS spatial join onto the 500 m domain polygons; 450 of the 622
transects fall in GIS 1-90, five per domain. Copied 2026-09-15 from
hard-structures/groin/HAT-groin-gis-analysis/gis_data/, see transects/README.md.
```

```text
Metadata columns copied from the dune-line feature when it carries them
(the 1997 lines do; 1984 and 1967 carry none).
```

```text
The layer is an ArcGIS join export: every column is prefixed with the
table it came from ("Transects_100m.LineID"). Strip to the leaf name and
keep the first of any duplicates (OBJECTID and Shape_Length appear twice).
```

```text
A transect runs due west from the datum line, so either end's
northing is the transect's.
```

```text
The same line-intersects-polygon join that placed the surveyed
transects, onto Hannah's whole-island polygons (2026-09-16). A
transect no polygon covers is dropped, as the surveyed join
dropped those in the slivers between polygons.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_transects()`**

```text
The 100 m transects with a domain number each.

Surveyed mode: the 450 the ArcGIS polygon join placed in GIS 1-90.
Extension mode: the 172 it left unplaced, numbered by their northing on
their whole-island polygon (hat_extension_domains.join_lines); the table is
also written to transects/transects_100m_ext.csv so the numbering is on
disk beside the layer it extends.
```

</details>

### 1-produce/island_offset_hybrid.py

Raw per-transect offsets -> the relative island offset per domain, padded to 120 domains for the model.

From the script's original header:

```text
Hatteras CASCADE Dune Offset Pipeline

This script:
1. Reads a raw feature-to-baseline intersection CSV (one vintage).
2. Calculates the relative offset per domain (metres, baseline = minimum).
3. Pads the result for CASCADE with the smooth wrap-around the model uses
   (cascade_pipeline.hindcast.pad_offset_ring): BRIE's domain is periodic,
   so the buffers carry the shoreline from GIS 90 back round to GIS 1 along
   a cubic Hermite matched to the island's end slopes. The padded file is
   therefore exactly what offset_mode "metres" hands Cascade.
4. Saves a diagnostic figure: the padded profile, and the shoreline angle
   BRIE reads between neighbouring domains against its ~42 degree limit.

UNITS: metres throughout, from the raw file's ORIG_LEN (EPSG:3725) to the
padded file. Nothing here converts to decametres.

PADDING HISTORY: until 2026-09-24 (every v1 build) the buffers were a local
slope segment plus a linear bridge, clipped at 0. The runner never used them
in metres mode -- it replaced them with this closure -- so the file and its
diagnostic showed a buffer the model did not see; those builds are now
<start>/<source>/superseded_20260924_pre-metres/v1. The current v1 (built 2026-09-24) writes the
closure itself (Hannah: "option (a)").
```

Notes that were in the code:

```text
The island_offset/ tree was renamed 2-brie-offset/ (raw_offsets/ plus
hindcast_<year>/), which left every path here dead. Anchored on the repo
root and on YEAR, so either hindcast start can be produced (2026-09-10).
```

```text
A start year is admissible here when hat_topo_version.DUNE_LINE_FOR_YEAR
pairs it with a dune-line vintage (1996 reads the 1997 line). Since
2026-09-15; before that the raw file had to exist under the period's own
name, which for 1996 meant a copy of the 1997 file.
```

```text
VERSIONED OUTPUT (2026-09-15). A start year can hold more than one build
when its dune line is re-digitised: 1996/v1/ is the build from the v1 1997
line (ArcGIS intersection), 1996/v2/ from duneline_1997_v2 (shapely
intersection, duneline_to_raw_offsets.py). Which one the runner reads is the
CURRENT file in <year>/, resolved by hatteras_site_config.island_offset_file.
Without --version the files land flat in <year>/, as 1984 and 2004 still are.
```

```text
raw_offsets/<vintage>_duneline_offset_raw.csv is the file the END-YEAR
TARGET loader reads too, so it always holds the CURRENT build of that
vintage. To rebuild an older version from its own raw file (each version
folder keeps a copy), name that file here.
```

```text
EXTENDED GEOMETRY (2026-09-16, the Pea Island extension experiment). A
named reach from hat_extension_domains: the surveyed raw for GIS 1-90 plus
raw_offsets/ext/<vintage>_duneline_offset_raw_ext.csv for the domains
beyond, zeroed on the SAME minimum as the surveyed build (checked against
<year>/CURRENT), padded with the same buffers, written to
<year>/ext/<geometry>/. Not a version: CURRENT is untouched.
```

```text
WHICH FEATURE THE OFFSET IS MEASURED FROM (2026-09-22). Every build until
then came from a digitised DUNE line. "shoreline" builds from the CoastSat
window mean instead (scripts/input_prep/5-scr/1-observations/mean_shoreline/),
which is a different FEATURE, not a newer reading of the same one -- so it
is a separate source with its own v1, never a v2 of the dune build.
The dune source keeps the flat <year>/v<n>/ layout it has always had, so
nothing the runner resolves moves; a non-default source nests one level
deeper, <year>/<source>/v<n>/. Everything after this point is identical:
same domain mean, same zeroing on the build's own minimum, same padding.
```

```text
A source names its own raw file. The shoreline's is named for the AVERAGING
WINDOW, not a vintage year (shoreline_raw_file_for_year), because a window
mean is what a satellite shoreline has instead of a survey date.
```

```text
ext/ sits under the SOURCE too (2026-09-22): an extended geometry is built
from the same feature as the surveyed reach it extends, and it is checked
against that source's CURRENT a few hundred lines below.
```

```text
What to call the feature in a title, an axis and a warning, so a
shoreline build is not labelled "Dune" on its own diagnostic figure.
```

```text
No build version in the title: the folder carries it, and a title naming
it went stale when the builds were renumbered v2 -> v1 (2026-09-28).
```

```text
The extension must not move the surveyed reach: same zero, same
values as the build the matrix runs read. A different minimum here
would shift every GIS 1-90 offset and the experiment would no
longer be about the buffer.
Through the resolver since 2026-09-22, when the builds moved under
<year>/<source>/: this joined <year>/ and CURRENT by hand and would
have looked for the surveyed build one level too high.
```

<details><summary>Function notes (the original docstrings)</summary>

**`shoreline_angles_deg()`**

```text
The angle BRIE reads between each padded domain and the next, wrapping
from the last back to the first (brie.py: atan2(diff(x_s), dy)).
```

**`plot_buffer_diagnostic()`**

```text
The padded profile and the shoreline angle BRIE reads, one figure.

(a) offset along the padded domains, buffers shaded, GIS numbering on the
real reach; (b) the angle between neighbouring domains, with the ~42
degree limit past which BRIE's shoreline goes anti-diffusive.
```

</details>

### 2-figures/compare_offset_sources.py

Two sources for one start year: the island offset from the dune line against the one from the shoreline.

From the script's original header:

```text
Two SOURCES for one start year: the dune line against the shoreline
Compares the island offset built from the digitised DUNE line with the one
built from the CoastSat SHORELINE, for one hindcast start. Written 2026-09-22,
the day the shoreline source was added.

NOT THE SAME QUESTION AS HAT_compare_offset_versions.py
    That script compares two VERSIONS of one source: the same feature
    re-digitised, where the difference is a correction and the interesting
    number is how many domains moved. This one compares two SOURCES: two
    DIFFERENT FEATURES on the island, where the difference is the beach
    between them and is supposed to be there. Different question, different
    framing, different caption -- so a sibling script rather than a flag.

THE TRAP THIS FIGURE EXISTS TO SHOW
    Correlate the two profiles and you get r = 1.0000, which looks like the
    two features agreeing about the island. They do not. Both are dominated by
    the same ~6.2 km of cape curvature, against which the 17 m standard
    deviation of their difference is 0.3%. So the figure states the
    correlation nowhere and shows the gap instead.

WHICH FRAME IS DRAWN, AND WHY IT HAD TO CHANGE
    Both profiles are drawn as the STATION FROM THE SHARED OFFSHORE DATUM, not
    as the min-zeroed offset the model is handed.

    The figure was drawn in the model frame until 2026-09-22, when Hannah
    asked the question that breaks it: "shouldn't the dune always be behind
    the shoreline?" It should, and on the ground it is -- the shoreline is
    seaward of the dune line in 90 of 90 domains, by 17-97 m. But the two
    model files are each zeroed on their OWN most seaward domain, and those
    minima are 45.6 m apart (dune 1953.198, shoreline 1907.607), so

        model_diff = -beach_width + 45.6 m

    and the dune line comes out apparently seaward wherever the beach is
    narrower than 45.6 m. That was 69 of 90 domains -- exactly the 69 the old
    figure shaded "dune line seaward", a physically impossible claim generated
    entirely by the zeroing.

    The datum frame has no such constant, and it is the SAME SHAPE (a model
    offset is this minus the build's own minimum). So the profiles are
    unchanged, the band between them is the beach, and it is on the correct
    side of the island everywhere. What the model reads is still one
    subtraction away, and the CSV holds both frames.

HOW IT IS DRAWN, AND WHAT THAT COSTS
    Six consecutive sections on a 2 x 3 grid, each a VERTICAL strip:
    alongshore up the page, cross-shore across it, south at the bottom, the
    axis inverted so the ocean is on the right. The panels read as the island.
    Why a GRID and not a row of six: see SECTIONS below -- columns take width
    from the very axis the gap is measured on.

    These are the ABSOLUTE profiles, and the beach is a small fraction of the
    cross-shore range they span, so the band is thin. Removing a smooth trend
    from both sources would open it up (14-43% of a panel rather than 2-3%),
    and was built and then taken back out (Hannah, 2026-09-22: "I actually
    don't like the detrend, I want to see the original shoreline shape"). The
    shape is the point; each panel therefore states its own beach range as a
    number.

OUTPUT   2-brie-offset/<year>/shoreline/<v>/comparisons/<a>_vs_<b>/
    (filed with the shoreline build it was drawn against since 2026-09-29;
    until then <year>/comparisons/<a>_vs_<b>/, which could not say which
    shoreline version it held once there were two)
    offset_<year>_<a>_vs_<b>.csv        per domain, both frames, columns named
                                        for the SOURCE not for a version
    offset_<year>_<a>_vs_<b>.png/.pdf   the two profiles as vertical strips of
                                        island, caption in CAPTIONS.md
    README.md                           what the folder is

    hat_topo_version.offset_source_comparison_dir resolves it.

USAGE
    python compare_offset_sources.py --year 1996
    python compare_offset_sources.py --year 1996 --a duneline --b shoreline
    python compare_offset_sources.py --year 1996 --shoreline-version v2

    Each source defaults to its CURRENT build; --duneline-version and
    --shoreline-version name one instead (2026-09-29, so a build can be
    compared before it becomes CURRENT).
```

Notes that were in the code:

```text
HOW EACH SOURCE IS NAMED AND DRAWN.

The house RdBu pair (C_1984 red / C_1997 blue) means EARLIER and LATER
vintage, and these two are not a time order -- drawing them red and blue
would say the dune line came first. Blue is already the CoastSat colour in
this style sheet, so the shoreline keeps it and the dune line takes the
accent purple.
```

```text
The island in sixths, drawn as vertical strips on a 2 x 3 GRID (2026-09-22).

THE GRID IS THE POINT, not the section count. These panels are columns, so
every extra one alongside takes width away from the offset axis -- the axis
the gap is measured on -- faster than the tighter zoom gives back. Measured,
as the widest gap actually rendered on paper at 300 dpi:

sections, across   panel width   gap on paper
3                  2.07 in       13 px
4                  1.52 in       10 px
6                  0.96 in        7 px   <- MORE panels, LESS gap
10                  0.52 in        6 px

Splitting into ROWS instead gives the width back, so the zoom is kept and
the offset axis is not squeezed:

6 as 2 x 3          2.07 in       15 px
9 as 3 x 3          2.07 in       23 px

6 on a 2 x 3 grid is the chosen point: half again the separation of 4 across,
15 domains a panel, and panels still tall enough (3.9 in) to read as strips
of coast.
```

```text
Stations grow LANDWARD from the offshore datum, so a - b is positive
where b lies seaward of a. With a=duneline and b=shoreline that is the
beach width.
```

```text
---- figure ---------------------------------------------------------- #
FOUR SECTIONS, DRAWN VERTICALLY, ABSOLUTE OFFSETS (2026-09-22, Hannah:
"I want to see the original shoreline shape ... should you do vertical
instead like the actual island shape").

Alongshore runs UP the page and the offset across it, so the panels read
as four consecutive strips of Hatteras with south at the bottom. The x
axis is INVERTED: an offset grows landward (the 100 m transects run due
west from the offshore datum), so inverting it puts the ocean on the
right, where it is.

No detrending: these are the profiles as the model reads them. That
costs legibility and it is worth being honest about how much -- a
quarter of the island still spans 800-2200 m of offset, so the widest
gap between the two sources is 1.3-3.3% of a panel's width. The filled
band is what carries it; in the flatter sections it is a sliver.
```

```text
DRAWN IN THE FIXED-DATUM FRAME, not the model frame (2026-09-22, Hannah:
"shouldn't the dune always be behind the shoreline?").

It should, and on the ground it is: the shoreline is seaward of the dune
line in 90 of 90 domains. But the MODEL files are each zeroed on their
own most seaward domain, and those two minima are 45.6 m apart, so
differencing them puts the dune line apparently seaward wherever the
beach is narrower than 45.6 m -- which is 69 of 90 domains. The earlier
version of this figure drew exactly that and labelled it "dune line
seaward", a physically impossible claim produced entirely by the zeroing.

The station from the shared offshore datum carries no such constant, and
it is the SAME SHAPE: a model offset is this minus the build's own
minimum. So nothing about the profiles is lost, the band between them is
the beach, and it is on the correct side everywhere.
```

```text
One band, one colour, one direction: the station grows LANDWARD from
the datum and the shoreline's is always the smaller, so the band is
always the beach. It needs no second colour for a sign that cannot
change.
```

```text
A domain is a count, so its ticks are whole numbers; the default
locator offered 42.5 and 45.0.
```

```text
The beach as a number, inside the panel: at this scale the eye
cannot measure the band, and the title has no room for it.
```

```text
ABOVE the panels, not below: at the bottom the legend and the shared
offset label are both "outside lower centre" and constrained_layout
stacks them on top of each other.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_town_bands_alongshore_y()`**

```text
Village spans, for a panel whose ALONGSHORE axis is the vertical one.

hat_figure_style.town_bands draws the same thing against a horizontal
alongshore axis and has no orientation switch, so this is its axhspan
twin. The spans themselves still come from the one owner,
hatteras_site_config.HATTERAS_ANNOTATIONS -- only the axis differs.
```

**`_window_gap_months()`**

```text
Months between the centre of the shoreline averaging window and the
dune line's survey date, or None if either is unavailable.

Computed, never typed: it was typed once, as "about nine months", and it
is 15.4 (found 2026-09-22 when the figure asked for it). The survey date
has one owner, duneline_endpoint.survey_date, so this asks it.
```

**`_shoreline_raw()`**

```text
The raw file a shoreline BUILD was made from: each version folder keeps
a copy (<window>_shoreline_offset_raw.csv), and its name is the window.
```

**`_shoreline_window()`**

```text
(first, last) day of the build's averaging window, read from its raw
file's name: 1995_1997 (calendar years) or 1995-10-12_1997-10-12.
```

**`_vintage_label()`**

```text
What the source IS for this start year, spelled for a legend: a dune
line carries an imagery vintage, a shoreline carries a window.
```

**`_unpadded()`**

```text
The 90-domain file the model would read from one source's build
(CURRENT unless a version is named), zeroed on that build's own most
seaward domain.
```

**`_raw_domain_means()`**

```text
Per-domain mean station from the shared offshore datum. One row per
transect first, so a domain with more transects does not weight twice.
```

</details>

### 2-figures/HAT_compare_offset_versions.py

Two builds of one start year: where did the island offset move?

From the script's original header:

```text
Two builds of one start year: where did the island offset move?

Compares two versions under data/hatteras_init/2-brie-offset/<year>/ (the
unpadded 90-domain files each version's island_offset_hybrid.py run wrote) and,
when both versions have a raw per-transect file, the ABSOLUTE distances behind
them. The unpadded files are each zeroed on their own minimum, so their
difference is the change the MODEL sees; the raw files share the offshore
datum, so their difference is where the dune line was actually moved.

Written 2026-09-15 for 1996 v1 (the ArcGIS intersection of duneline_1997) vs
v2 (the shapely intersection of duneline_1997_v2, local corrections only).
The raw comparison uses the v1 line re-intersected by the SAME shapely script,
so the 1 m station convention of the GIS export (see
duneline_to_raw_offsets.py) does not appear as a change.

Outputs, in <year>/<b>/:
    offset_<year>_<a>_vs_<b>.csv    per domain: both versions, both frames
    offset_<year>_<a>_vs_<b>.png/.pdf   two panels, caption in CAPTIONS.md

USAGE
    python HAT_compare_offset_versions.py --year 1996 --a v1 --b v2         --raw-a <path to the v1 line's shapely raw> --raw-b 1997_v2_duneline_offset_raw.csv

    # a superseded build: --a is its folder, --label-a what the outputs call it
    python HAT_compare_offset_versions.py --year 1996         --a superseded_20260919_pre-redigitized/v2 --label-a superseded_v2 --b v1         --raw-a <its raw> --raw-b <v1's raw>

--label-a/--label-b (2026-09-23) name a build in the file stem, the columns,
the legend and the caption. They default to --a/--b; they exist because a
superseded build's folder is a path, and a slash cannot go in a file name.
```

Notes that were in the code:

```text
Under the SOURCE's folder since 2026-09-22 -- this joined <year>/ and the
version by hand, which after the split would have written the comparison
into a <year>/v2/ that no longer exists.
```

```text
What panel (b) is a picture OF depends on what changed between the raw
files: the line, the method that measured it, or (unrecorded) unknown.
Two shapely raws differ ONLY where the line does -- the intersection is
deterministic against the same transects -- so that case is a line
change even when both name the same file (edited in place).
```

```text
Two CoastSat window means: nothing was digitised, the averaging
window changed (2026-09-29, the DEM-centred shoreline v2).
```

<details><summary>Function notes (the original docstrings)</summary>

**`_raw_provenance()`**

```text
(method, line file, built date) behind one raw file, read from the file.

duneline_to_raw_offsets.py stamps `built_by`, `built` and `duneline_file`
on every row; the ArcGIS exports it replaced (raw_offsets/superseded_
20260915_gis_exports/) carry none, so an export's line is unknown unless
the caller names it (--line-a/--line-b). Until 2026-09-23 panel (b) was
always titled "the re-digitised line", which was wrong for 1984 and 2004:
there the line is the same geojson and only the intersection method
changed.

A FILE NAME IS NOT A LINE. duneline_2009.geojson was re-digitised in place
on 2026-09-18, so both 2010 raws name it and differ by up to 66 m. Hence
the caller's rule: two shapely raws are compared as a line change whatever
the names say, because the intersection is deterministic.
```

</details>

### 2-figures/HAT_plot_offset_profile_1to1.py

The padded island offset profile at true 1:1 scale, so the buffer's shape reads as it is.

From the script's original header:

```text
The padded offset profile at true 1:1 scale
The buffer diagnostic that island_offset_hybrid.py draws stretches 7.5 km of
cross-shore offset across a panel that spans 72 km alongshore, so the island
reads as a deep V. This draws the same padded profile with equal axes: 1 m
alongshore is 1 m cross-shore, the way a map draws it (Hannah, 2026-09-16).

    python HAT_plot_offset_profile_1to1.py --year 1996 --geometry n115
    python HAT_plot_offset_profile_1to1.py --year 1996              # the CURRENT surveyed build
    python HAT_plot_offset_profile_1to1.py --file 1996/v1/Island_Dune_Offsets_1996_PADDED_120.csv
    python HAT_plot_offset_profile_1to1.py --all                    # every padded build on disk

Writes <build dir>/<the build's own stem>_buffer_diagnostic_1to1.png beside
the original diagnostic (exaggerated axes, kept), with the PDF and caption
under supporting/. Named `_v2` until 2026-09-23, which read as a build number
inside a v<n>/ folder.
```

Notes that were in the code:

```text
Through hat_topo_version since 2026-09-22, when every build moved under
<year>/<source>/. This read CURRENT out of <year>/ itself and would now
find no CURRENT there, silently falling back to the year folder.
```

```text
Named from the file it was drawn FROM, not from the dune stem
(2026-09-22): this was hardcoded, so the shoreline build's diagnostic
landed in 1996/shoreline/v1/ calling itself Island_Dune_Offsets.
```

<details><summary>Function notes (the original docstrings)</summary>

**`geometry_of()`**

```text
(geometry, first, last) from a padded file's length; None if no
geometry pads to that many domains.
```

**`every_padded_file()`**

```text
Every padded build under 2-brie-offset, superseded ones included.

One glob per SOURCE (2026-09-22). This matched only the dune stem, so
--all quietly skipped every shoreline build; the stems come from
hat_topo_version so a new source cannot be forgotten here again.
```

</details>
