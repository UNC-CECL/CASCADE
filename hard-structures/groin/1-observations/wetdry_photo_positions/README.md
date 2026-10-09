# wetdry_photo_positions — shoreline position per domain from the surveys

```
HAT_geometric_distance_sanity_check.py   the per-domain distances -> the three CSVs here, two check figures -> figures/
HAT_groin_two_shorelines_figure.py       the groin in one picture, from Change_from_wetdry -> figures/
Change_from_wetdry_1967_D2_D12.csv       24 dated wet/dry surveys, GIS 2-12, change since 1967,
                                         landward-positive: the observed gap every fit reads
Change_from_duneline_1967_D2_D12.csv     the same from the dune lines
geometric_distances_all_shorelines.csv   raw distance to the datum line, every shoreline
figures/                                 groin_two_shorelines, sanity_check_spatial_map, sanity_check_distance_profile;
                                         PDFs and CAPTIONS.md in supporting/
```

The CSVs stay at this level because the hindcast fits and the sweep config read them by path. The
geojsons the sanity check reads are the 1967 rig's inputs
(`../../3-hindcast/1-dipole-1967-2017/inputs/groin_init/island_offset/input/`). Until 2026-10-08 the
script sat in the rig's `input_prep/shoreline_position/` and its products in the rig's
`shoreline_position_output/`; until 2026-10-09 the two sanity-check PNGs sat loose here and the
two-shorelines figure had its own `1-observations/figures/` folder.

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning and choices behind it are here. Moved out of the scripts on 2026-10-01, when they were brought in line with `scripts/STYLE.md` (code unchanged, proven with `style_equivalence_check.py`; path fixes listed per script).

### HAT_geometric_distance_sanity_check.py

A general tool: per-domain mean distance to the offshore datum line for any
mix of shoreline geojsons (dune lines, wet/dry lines), with two figures so the
numbers can be checked by eye before they are trusted. It replaced
`HAT_duneline_geometric_distance.py` and `HAT_wetdry_geometric_distance.py`.

**Labels.** Single-feature files with no `year` property (e.g.
`duneline_1967.geojson`) are labelled from the filename; multi-feature dated
files (`wet_dry_shorelines_groin.geojson`) as `wetdry_<year>`, with the month
appended when a year repeats. Each label carries a sort key (year plus a month
fraction) so figures colour chronologically whatever the label text.

**Distance.** Each shoreline is clipped to each domain polygon and sampled at
`N_SAMPLES` (25) points per segment; the domain value is the mean distance to the
datum line.

**Checks.** `VALIDATION_KEY` (duneline_1967) is compared, min-subtracted,
against `VALIDATION_KNOWN_GOOD` within `VALIDATION_TOLERANCE_M` (15 m), the
values from the earlier dune-line script. `SEAWARD_CHECK_PAIRS` prints whether
the first line is seaward (smaller distance) of the second at every shared
domain: wet/dry should be seaward of the dune line everywhere.

**Outputs.** The raw-distance CSV, one change-from-reference CSV per
`REFERENCE_KEYS` entry (the format `HAT_groin_effect_comparison.py`'s
observed-change logic reads; `[]` writes the raw CSV only), the spatial map
(real coordinates, zoomed to the domains: the figure most likely to catch a
wrong CRS, wrong file or reversed geometry) and the distance profile, against
`FIGURE_REFERENCE_KEY`.

**Path fix (2026-10-01).** `INPUT_DIR` was the dead literal
`/hard-structures/groin\hindcast_groin_test\input_prep\shoreline_position\input`
and `OUTPUT_DIR` `/scripts/groin/HAT-buxton-hindcast-groin-test/input_prep/shoreline_position/shoreline_position_output`.
Both now resolve from the repo root: the inputs to
`groin_init/island_offset/input/` (where all seven geojsons are) and the output
to `1-observations/wetdry_photo_positions/` (where its products are).

From the script's original header:

```text
HAT_geometric_distance_sanity_check.py
========================================
General-purpose tool: computes per-domain distance-to-datum for ANY set of
shoreline geojsons (dune lines, wet/dry lines, or a mix -- doesn't matter
which), using the offshore datum line + domain polygons, and produces two
diagnostic figures so the result can be checked BY EYE, not just trusted
numerically.

Consolidates and generalizes HAT_duneline_geometric_distance.py and
HAT_wetdry_geometric_distance.py into one reusable tool. Handles both file
styles automatically:
  - one feature per file, no 'year' property (e.g. duneline_1967.geojson)
    -> labeled from the filename
  - many dated features in one file (e.g. wet_dry_shorelines_groin.geojson)
    -> labeled from each feature's own 'year' (+ month, if a year repeats)

FIGURES
-------
1. Spatial map (real coordinates): domain polygons, the datum line, and
   every shoreline overlaid, colored by year (older=blue, newer=red via
   colormap) -- lets you SEE whether shorelines actually fall inside the
   expected domains and stay roughly parallel to the datum line. This is
   the figure most likely to catch a wrong CRS, wrong file, or reversed
   geometry immediately, before ever looking at a number.
2. Distance profile: distance-to-datum (or change-from-REFERENCE_KEY, if
   set) per domain, one line per shoreline, same color scheme as the map.

Also runs the same validation/sanity checks already established:
  - if VALIDATION_KNOWN_GOOD is set, checks a named shoreline's min-
    subtracted values against it (as done for 1967 dune line).
  - prints a pairwise seaward/landward check between any two named
    shorelines you list in SEAWARD_CHECK_PAIRS (e.g. wet/dry vs dune line
    at the same year -- wet/dry should be seaward everywhere).

Author: Hannah A. Henry, UNC CECL
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_shorelines()`**

```text
Load every shoreline feature in a file, auto-labeling each one.

- Multi-feature files with a 'year' property (e.g. wet/dry): labeled
  'wetdry_<year>' or 'wetdry_<year>_<month>' if the year repeats.
- Single-feature files with no 'year' property (e.g. duneline_YYYY.geojson):
  labeled from the filename stem.

Returns dict {label: (geometry, sort_key)}; sort_key is a float year
(with a small month-based fraction added for duplicate years) so
figures can be colored/ordered chronologically regardless of label text.
```

**`fig_spatial_map()`**

```text
Real-coordinate map: domain polygons, datum line, every shoreline
colored chronologically. The figure most likely to catch a CRS/file/
geometry mistake by eye before trusting any number.
```

**`fig_distance_profile()`**

```text
Distance-to-datum (or change-from-reference) per domain, one line
per shoreline, colored chronologically to match fig_spatial_map.
```

**`save_change_table()`**

```text
Change-from-reference_key CSV for every shoreline, ready to feed
HAT_groin_effect_comparison.py's observed-change logic directly. Mirrors
HAT_duneline_geometric_distance.py / HAT_wetdry_geometric_distance.py's
combined-change-table format.
```

</details>

<details><summary>Notes that were comments in the code</summary>

- Change-from-reference CSVs are saved for each key in this list (one CSV per key), ready to feed straight into HAT_groin_effect_comparison.py's observed- change logic. Add both a dune-line and a wet/dry reference if you want to compare either baseline. Leave empty to skip (raw-distance CSV only).
- Which reference (if any) the *figure's* second panel shows change relative to -- independent from REFERENCE_KEYS above (that controls saved CSVs).
- Optional validation against a known-good min-subtracted series (see HAT_duneline_geometric_distance.py). Set KEY to None to skip.
- Optional pairwise seaward/landward checks: (should_be_seaward, reference). Prints whether the first is seaward of (smaller distance than) the second at every domain both have data for.
- Zoom to the domain polygons' own extent, not the full (much longer) shoreline/datum-line length -- that's the point of this figure.

</details>

### HAT_groin_two_shorelines_figure.py

The groin in one picture: the shorelines at GIS 5 and GIS 6 from the 24 dated wet/dry surveys in `Change_from_wetdry_1967_D2_D12.csv`, and the gap between them. Writes `figures/groin_two_shorelines.png`, its PDF and caption into `figures/supporting/`. Nothing here runs the model. Brought in line with `scripts/STYLE.md` on 2026-10-01. Until 2026-10-09 it sat in its own `1-observations/figures/` folder; it moved here, beside the table it reads, to match the other two analyses (output path only changed; regenerated, same gap values).


The gap between the sheltered (GIS 6) and unprotected (GIS 5) shorelines is
the groin's effect, and needs no definition to read. Plotted seaward-positive;
the timeline figure keeps the source's landward-positive sign.

From the script's original header:

```text
The groin in one picture: two shorelines, and the gap between them.

WHY THIS EXISTS ALONGSIDE THE TIMELINE FIGURE
    The fillet is a DERIVED quantity -- a difference between two domains -- so
    it needs a sentence of explanation before anyone can read it, and it invites
    the wrong reading ("150 m of new beach") when it actually means "150 m of
    relative position".

    Plotting the two shorelines themselves needs no explanation. Both retreat.
    One retreats far less. The gap between the lines IS the groin's effect, and
    its width over time is the whole story:

        widening gap  ->  the groin is trapping
        closing gap   ->  the groin has stopped

    Everything the timeline figure says is visible here without a definition.

WHAT IT IS FOR
    Deciding how to apply the module. The groin's job in the model is to hold
    the updrift line above the downdrift one. Reading the gap tells you directly
    what the module has to do in each hindcast period, and where it cannot help.

DATA
    Shoreline change since a fixed 1967 datum at the two domains flanking the
    structure -- D5 downdrift, D6 updrift -- from 24 dated wet/dry surveys
    (`Change_from_wetdry_1967_D2_D12.csv`, GIS analysis).

    Plotted SEAWARD-POSITIVE so the lines fall as the shoreline retreats, which
    is how a reader expects an eroding coast to look. The source table is
    landward-positive, so it is negated here; the timeline figure keeps the
    source convention, which is why its curve rises where these fall.

STYLE, 2026-09-11
    Drawn under the project house style (`scripts/site_layer/hat_figure_style.py`): a
    190 mm column, so the 9 pt type on the canvas is 9 pt on the page rather
    than the 4 pt a 13 in canvas reduced to. The sheltered side is the ACCENT
    purple and the unprotected side BASE grey -- the RdBu red/blue pair is
    reserved for the 1984/1997 vintages, and these two lines are places, not
    vintages. Nothing on the canvas that belongs in a caption: the title
    sentence and the two footnote paragraphs are now in CAPTIONS.md beside the
    image, and the two gap widths they quoted as hand-written constants are
    measured off the plotted series.

Usage:
    python HAT_groin_two_shorelines_figure.py

Writes groin_two_shorelines.png (and .pdf) beside this file.

Author: Hannah A. Henry, UNC CECL
```

Notes moved out of the code (two_shorelines(), its docstring, as of the original):

```text
{year: (updrift, downdrift)} change since 1967, SEAWARD-positive.
```

Notes moved out of the code (endpoint_rate(), its docstring, as of the original):

```text
Rate of `series` between the first and last survey inside [start, end].

Endpoints, not a regression: the two rates quoted for the hindcast periods
have always been endpoint differences, and the caption says so.
```

Notes moved out of the code (period_strip(), its docstring, as of the original):

```text
The hindcast windows as a band against the top edge, named once.

A strip rather than a full-height wash: the panel already carries a fill
between the two shorelines, and two full-height greys on one panel cannot
be told apart.
```

Notes moved out of the code (comment above line 79, as of the original):

```text
# The sheltered side is the thing under test; the unprotected side is the
# baseline it is read against, and the band between them is the effect.
```

Notes moved out of the code (comment above line 174, as of the original):

```text
# The gap at its widest and at the last survey, measured off the series
# rather than written in: both were hand-set constants until 2026-09-11,
# and the widest is NOT where the constants put it. They read 2004, the end
# of period 1, at 150 m; the largest gap in the record is the 1995 survey
# at 155 m, a single-survey spike driven by the downdrift side. The label
# goes where the data is and the caption carries both numbers.
```

