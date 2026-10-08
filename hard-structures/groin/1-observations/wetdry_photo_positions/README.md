# wetdry_photo_positions — shoreline position per domain from the surveys

One script, `HAT_geometric_distance_sanity_check.py`, and its products beside it:
`Change_from_wetdry_1967_D2_D12.csv` (24 dated wet/dry surveys, GIS 2-12, change since 1967,
landward-positive: the observed gap every fit reads), `Change_from_duneline_1967_D2_D12.csv`,
`geometric_distances_all_shorelines.csv` and the two by-eye figures. The geojsons it reads are the
1967 rig's inputs (`../../3-hindcast/1-dipole-1967-2017/inputs/groin_init/island_offset/input/`).
Until 2026-10-08 the script sat in the rig's `input_prep/shoreline_position/` and its products in
the rig's `shoreline_position_output/`.

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
