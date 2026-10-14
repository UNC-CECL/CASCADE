# aerial_photo_gap — the groin's effect from the aerial photos alone

**The photo record, on its own.** 24 dated aerial-photo wet/dry lines, 1967-2023, measured as each GIS
domain's distance to a fixed offshore datum line and taken as change since 1967. Two things come out of
it: the per-domain position table every hindcast fit reads (`Change_from_wetdry_1967_D2_D12.csv`), and
the gap across the groins from the photos alone. CoastSat on its own is `../2-coastsat_gap/`; the two
records side by side are `../3-photos_vs_coastsat/`. Until 2026-10-14 this folder was
`1-wetdry_photo_positions/`, and the photo gap was drawn by `HAT_groin_two_shorelines_figure.py` (below).

```
HAT_geometric_distance_sanity_check.py   the per-domain distances -> the three CSVs here, two check figures -> figures/
photo_gap.py                             the gap across the groins from the photo table -> tables/
photo_gap_figures.py                     figures from tables/ -> figures/photo_gap_<n>_<what>.png
Change_from_wetdry_1967_D2_D12.csv       24 dated wet/dry surveys, GIS 2-12, change since 1967,
                                         landward-positive: the observed gap every fit reads
Change_from_duneline_1967_D2_D12.csv     the same from the dune lines
geometric_distances_all_shorelines.csv   raw distance to the datum line, every shoreline
tables/                                  photo_gap_annual.csv, photo_era_rates.csv
figures/                                 photo_gap_1 to _5 in the order of the analysis; the two
                                         sanity_check_* figures; PDFs and CAPTIONS.md in supporting/
```

```
python HAT_geometric_distance_sanity_check.py   # only when a shoreline geojson changes
python photo_gap.py
python photo_gap_figures.py                     # the GIF takes a few minutes
```

The CSVs stay at this level because the hindcast fits and the sweep config read them by path. The
geojsons the sanity check reads are the 1967 rig's inputs
(`../../../3-hindcast/1-dipole-1967-2017/inputs/groin_init/island_offset/input/`). Until 2026-10-08 the
script sat in the rig's `input_prep/shoreline_position/` and its products in the rig's
`shoreline_position_output/`; until 2026-10-09 the two sanity-check PNGs sat loose here and the
two-shorelines figure had its own `1-observations/figures/` folder.

## The photo gap (2026-10-14)

How a survey becomes a per-domain position is figure 2: the mean due-east distance from the domain's piece of line
to a fixed north-south datum line ~7.4 km offshore, as change since 1967. Each side is its domain's row in the change table, negated to seaward positive; the gap is north (GIS 6,
updrift) minus south (GIS 5, downdrift), so a rising gap is the groin holding the north side. The two
1972 surveys are averaged. Era rates are least-squares slopes over the surveys inside each era.

1. **Both sides retreated; the gap is the north side retreating less.** By 1996 the north side was 49 m
   landward of 1967 and the south side 183 m, a 135 m gap. *Figure 4.*
2. **The gap built fast, then slowly.** +11.2 ± 1.4 m/yr over 1967-1978 (5 surveys), then +1.6 ± 1.3 m/yr
   to the 155 m peak in 1995 (8 surveys). *Figure 5.*
3. **After the last repair it narrowed slowly.** −1.4 ± 0.6 m/yr over 1996-2016 (6 surveys), then it fell
   to 63-94 m around the 2017 fill and the reported structural losses (2015-2021); 74 m in 2023.
4. **The photos are sparse where it matters.** One survey between 1997 and 2008 (2004, at 150 m), so the
   photos cannot date the turn on their own. When it happened is read from CoastSat
   (`../2-coastsat_gap/`), and where the photos and CoastSat differ is in `../3-photos_vs_coastsat/`, whose figures 2-3 are figure 4
   here with CoastSat's annual record added in the same frame.

The figures follow the analysis step by step: the raw lines, how a line becomes a number, the
whole table, the two domains either side of the groins, then the rates.

| # | figure | question | answer |
|---|---|---|---|
| 1 | `photo_gap_1_wetdry_lines` | what was digitized? | 24 dated wet/dry lines, 1967-2023, over GIS 2-12 and close up at the groins; the lines fan out most just south of the groins |
| 2 | `photo_gap_2_measuring_position` | how does one survey become a position per domain? | clip the line at each domain, sample 25 points, take each point's due-east distance to the north-south datum line ~7.4 km offshore, average, subtract 1967. Because the datum runs north-south, it is an east-west distance: the shoreline at GIS 5-6 runs about 21° from north, so a change reads about 7% larger than the shore-normal change, on both sides about alike |
| 3 | `photo_gap_3_positions_by_domain` | where did each domain's shoreline go, survey by survey? | the position table drawn whole: everything retreats, GIS 4-5 just south of the groins far more than GIS 6-8 north of them |
| 4 | `photo_gap_4_two_sides` | where was each side, and how wide was the gap? | GIS 6 minus GIS 5: the gap grew to 155 m in 1995, was 135 m in 1996 and 74 m in 2023. Animated as `photo_gap_4_two_sides.gif` (git-ignored, rerun to make it) |
| 5 | `photo_gap_5_stage_rates` | how fast did the gap change in each stage? | +11.2 m/yr building (1967-78), +1.6 to the peak (1978-95), −1.4 after the last repair (1996-2016) |
| 6 | `photo_gap_6_placement` | where along GIS 4-8 did the shoreline build out and retreat, and how big is the gap across GIS 5/6 and 6/7? | in 100 m strips the shoreline rises steadily across the field, from the deepest retreat in GIS 4-5, through zero between the groins (about 170 m north of the southern groin), to the most building-out (+20-30 m by 1989-95) at the 6/7 boundary, falling off through GIS 7. The 5/6 gap is 97-155 m over 1978-97, the 6/7 gap 14-44 m. Added 2026-10-14 |

Until 2026-10-14 figures 4 and 5 were `photo_gap_1_two_sides` and `photo_gap_2_stage_rates`. The two
`sanity_check_*` figures, written by the position script, are its by-eye checks and sit outside this
sequence.

| table | contents |
|---|---|
| `photo_gap_annual.csv` | each survey year: north and south position since 1967 (seaward positive), the gap, surveys averaged |
| `photo_era_rates.csv` | gap rate ± SE per stage, with the survey count |
| `photo_boundary_gaps.csv` | each survey year: the gap across GIS 5/6 and across 6/7 (north domain minus south, since 1967) |
| `photo_alongshore_change.csv` | each survey: change since 1967 in 100 m strips of northing across GIS 4-8, the position table's method per strip |

## Where should the model groin sit: GIS 5/6 or 6/7? (2026-10-14)

The groin field fills almost all of GIS 6: the southern groin is about 75 m inside the 5/6 boundary, the
northern groins end 25-35 m short of 6/7. So the two candidate boundaries are the two ends of one field.

From the photos alone (figure 6): the starved side begins at the southern groin and the fillet sits against
the northern groins, at 6/7. Domain means hide this: GIS 6 averages a ramp from about −50 m to +30 m
(1989), so its mean held near zero, the 5/6 gap is large and the 6/7 gap small. A single boundary groin can
place either the starved side (5/6) or the fillet (6/7) at a boundary, not both. The CoastSat rates and the
two records together are in `../2-coastsat_gap/` figure 9 and `../3-photos_vs_coastsat/` figure 6.

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
to `1-observations/analysis/1-aerial_photo_gap/` (where its products are).

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

### photo_gap.py and photo_gap_figures.py

Added 2026-10-14, when the observations were split by source. `photo_gap.py` carries the photo half of
what `../2-coastsat_gap/coastsat_gap_across_groins.py` used to compute (the photo gap and its era rates);
the GIF in `photo_gap_figures.py` is the one that was `gap_across_groins_0_what_the_gap_is.gif`.

### HAT_groin_two_shorelines_figure.py (deleted 2026-10-14)

The groin in one picture, the two shorelines at GIS 5 and 6 and the gap between them. Replaced by
`photo_gap_figures.py` figure 4, which draws the same series. Its caption was out of date: it marked the
last repair at 1996 (the record says 1995, see `../2-coastsat_gap/README.md`), used the old 1984-2004 and
2004-2024 hindcast periods, and quoted a module result from an earlier sweep. The script and its notes
are in git at 563a542e (`hard-structures/groin/1-observations/analysis/1-wetdry_photo_positions/`).
