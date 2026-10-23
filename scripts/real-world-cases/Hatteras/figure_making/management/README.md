# figure_making/management — NC-12 and nourishment

The management the hindcast applies, drawn from `hatteras_site_config`: a
table of the rules and a timeline of every event. Output goes to
`output/figures/3-model-inputs/4-management/`.

```
management_rules_table_figure.py  the rules as a table (manuscript and slide versions)
management_timeline_figure.py     every management event as a domain-by-year timeline
```


## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### Deleted 2026-10-01

- `management_investigation_plot.py` — natural vs roadway vs nourishment from
  saved runs; its `RUN_PATHS` named `HAT_1984_2004_*_Hs2p5` runs that no
  longer exist anywhere under `output/raw_runs/`.
- `superseded_20260918/diagnose_road_drowning.py` (with its WHY.md) — retired
  2026-09-18, the last caller of the old `domain_{n}_topography_<year>.npy` names.

To recover, by the path it was last committed under:

```
git log --diff-filter=D --oneline -- scripts/figure_making/management/<path>
git show <commit>^:scripts/figure_making/management/<path>
```

### management_rules_table_figure.py

The management rules the hindcast applies, as a table: a manuscript version and a slide version.

From the script's original header:

```text
The management rules the hindcast applies, as a table.

TWO OUTPUTS, ONE SOURCE
    rules_table.png         the manuscript table: 190 mm, booktabs rules, no
                            fills, no title and no note on the canvas (house
                            style -- that text is in CAPTIONS.md beside it).
    rules_table_slide.png   the same rows for a projector: wider, larger type,
                            a title and the note ON the canvas, because a slide
                            has no caption to carry them. This is the ONE
                            deliberate departure from the "nothing on the
                            canvas" rule in figure_making/STYLE.md.

WHAT CHANGED 2026-09-17, and why each change was needed
  * THE HEADER WAS PRINTING A FILE PATH. The column label was 'Model\nDomains';
    the 2026-09-14 path-anchoring pass read "\nDomains" as a path fragment and
    rewrote the literal to `str(_PATH_REPO / "nDomains")`, so every rendering
    since carried C:\Users\...\CASCADE\nDomains across the header row.
  * THE NUMBERS WERE STALE. Three rows were typed by hand in an earlier draft
    and never followed the config: Rodanthe read D85-88 / 1,620,000 cy against
    the configured 84-89 / 1,600,000, Avon read D23-26 against 21-28 (the
    config's own comment names 23-26 as the superseded footprint, corrected
    2026-08-22), and the no-road reach read D1-6 against
    HATTERAS_FIRST_ROAD_DOMAIN = 9. Every domain span, year and volume in this
    figure is now READ FROM hatteras_site_config, which is what the caption
    always claimed. Only the prose in the Description column is editorial.
  * TYPE AND WIDTH. It was set in DejaVu Serif on a 15 in canvas, so its 9 pt
    body reduced to about 4.5 pt in a two-column manuscript -- the exact
    failure hat_figure_style.figsize() exists to prevent. It is Arial on 190 mm
    now, and the type is the size it will be printed at.
  * THE DECORATION. A black title bar, a grey header bar, zebra-striped fills
    and boxed column rules are web-table conventions; a journal table is three
    horizontal rules and white space. The only colour left is the accent tick
    beside each section heading, in the SAME two colours the other management
    figures use for fill (C["ADDED"]) and NC-12 (C["ROAD"]).
  * WRAPPING IS MEASURED, not guessed. Line breaks came from textwrap at a
    hand-tuned 80 characters, which is a different physical width at every font
    size; text is now wrapped against the measured width of the column, and row
    heights follow the number of lines that produces.
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5), so this block is
independent of whatever this script calls its own repository variable.
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
The Description column is the only editorial text in the table: what the
project was and where it was placed, keyed by the project name in the config
so a renamed or added project shows up with its own config note rather than
silently inheriting someone else's sentence.
```

```text
The permanent settlement footprints. Inside a village the road is a
street network that is maintained, not relocated, so the rule is off.
```

```text
Location | Domains | Years | Description. The first three are sized to
their own widest entry plus a gutter, so the Description column gets
everything that is left rather than a share fixed by hand.
```

```text
Lay the rows out: wrap every cell, then give the row the height its
tallest cell needs.
```

```text
SECTION HEADING. An accent tick in the colour this rule wears in the
other management figures, then the name and a one-line gloss.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_measurer()`**

```text
width_in(text, fontsize, bold, italic) -> rendered width in inches.

A scratch canvas, so the real figure can be sized from the lines the text
actually takes rather than from a character count that means a different
width at every font size.
```

**`put()`**

```text
Top-aligned in the cell, which is how a wrapped journal table sets
a row whose columns have different line counts.
```

</details>

### management_timeline_figure.py

The management of NC-12 and the beach as a timeline, domain against year: the runs and the whole record.

From the script's original header:

```text
The management of NC-12 and the beach as a timeline: domain against year.

TWO OUTPUTS, ONE BUILDER
    timeline_1996_2024.png   THE RUNS. Starts at the first period start, so
                             every mark is an event some run applies.
    timeline_1984_2024.png   THE RECORD. The whole management history, with
                             1984-1996 drawn as the stretch before the runs.
                             The 1989 Pea Island relocation only appears here.

WHAT CHANGED 2026-09-17, and why
  * THE PERIODS WERE STALE. Y0/PBREAK/Y1 were typed as 1984/2004/2024, the
    pair the project ran first. hatteras_site_config defines four periods now
    and the ones in use are 1996->2010 and 2010->2024 (Hannah, 2026-09-17:
    "1996 and 2010 are my main starting years now"). The span, the break and
    both period labels are READ FROM HATTERAS_PERIODS -- change PERIOD_STARTS
    and both figures follow.
  * THE RED WAS NOT THE HOUSE RED. The village bands carried a one-off salmon
    (#e6b39a, the settlement tint of the site figures) and it sat right beside
    C_1984_FILL on the Period 1 bar, so the period pair read as two oranges
    rather than as the red/blue vintage pair every other figure uses (Hannah,
    2026-09-17: "make the red match the red we have been using instead of this
    orange"). The period bars now carry the vintage pair with their saturated
    edge and text -- the same red as the 1996 NC-12 alignment on the map --
    and the ZONES moved to neutral greys, because on this figure colour is
    reserved for management and period. The map keeps the settlement tint,
    where nothing else is red.
  * ONE VISUAL LANGUAGE WITH THE MAP. The two figures drew the same three
    management families in different encodings -- here an orange bar, a black
    bar and a blue block; there an orange tint, a grey tint and a hatch. Both
    use the map's now, so a reader learns the key once.
  * THE CALLOUT BOXES ARE GONE. Rounded white boxes with coloured borders and
    curved arrows are a slide idiom, and six of them, each hand-placed, were
    most of the ink. Labels are plain text at a MEASURED, collision-checked
    position with a hairline leader -- what the reach figures use.
  * THE LEGEND LOST ITS BOX, ITS TITLE AND HALF ITS ENTRIES. The zone swatches
    repeated the zone names printed up the right-hand side, and the period
    swatches repeated the labelled period bar.
  * THE PANEL WAS THREE QUARTERS EMPTY at 5.6 in of height.
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5), so this block is
independent of whatever this script calls its own repository variable.
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
The period starts in use. Each one's end comes from HATTERAS_PERIODS, so the
bar below the axis states the run windows rather than a memory of them.
```

```text
COLOUR IS FOR MANAGEMENT AND PERIOD ON THIS FIGURE. The zones are greys so
that nothing competes with the vintage red/blue; the site map keeps the warm
settlement tint, where nothing else is red.
```

```text
The map's encoding, so a reader learns one key: fill orange, relocation the
NC-12 ink, the bridge a hatch.
```

```text
EVERY BOX BUT ITS OWN. A label sits beside the bar it names, so the
clearance pad around it overlaps that bar by design; testing against
it let the bar veto its own label, and the Buxton fill label -- boxed
in by Avon above and the axis below -- was dropped entirely.
```

```text
NEVER DROP A LABEL. If nothing is clear, take the first position that
at least fits on the panel: a label that overlaps is a flaw a reader
can see and work around, one that is missing is a figure that lies.
```

```text
Explicit margins: the locator strip hangs off the left of the axes and
the zone names off the right, and with a tight bbox both were growing
the saved figure past the 190 mm double column -- which is how 8 pt type
becomes 7 pt on the page. Reserving the room here keeps the content
inside the width figsize() asked for.
```

```text
PERIOD BAR. The vintage pair, saturated on the edge and in the text, so
Period 1 is the same red as the 1996 NC-12 alignment on the map.
```

```text
no years on this one: the axis and the dashed break carry them,
and the segment is too narrow at 40 years to hold them
```

<details><summary>Function notes (the original docstrings)</summary>

**`events_in()`**

```text
(year0, year1, gis_lo, gis_hi, kind, label) for the window, from the
config. An event before `y0` is dropped; on the runs window that is what
removes the 1989 Pea Island relocation, which precedes both periods.
```

**`_clear()`**

```text
True if the segment misses every rectangle in `rects`.

Sampled rather than clipped analytically: at 64 steps the spacing is far
finer than the narrowest event box, and the intent is legibility, not a
proof.
```

**`place()`**

```text
Label each event without printing on a bar, on another label, OR
dragging its leader across one.

A label starts beside its own bar and, if that position is taken, tries the
next candidate offset; a hairline leader is drawn whenever it ends up away
from its anchor. Positions come from the RENDERED size of the text, so this
holds if the wording or the window changes -- which the six hand-placed
callout boxes it replaces did not.

THE LEADER IS CHECKED TOO (Hannah, 2026-09-17). Testing only the label box
let a position pass whose leader then ran straight down through a coloured
box on its way to the bar below: the Buxton fill label sat above Avon, so
its leader crossed the Avon fill. A candidate is now rejected unless the
line it would draw also misses every box but its own.
```

**`build()`**

```text
One timeline. `y0` is the first year on the axis; `pre_run` draws the
stretch before the first period start as unmodelled context.
```

</details>
