# comparison - the groin's effect against the observed shoreline

| script | what it does | writes |
|---|---|---|
| `HAT_groin_effect_comparison.py` | a no-groin and a groin 1967-2017 run against the observed wet/dry shoreline at each checkpoint year | `COMPARISON_OUTPUT_DIR`: an overview, one figure per checkpoint (each with a `_v2`), a shoreline GIF |

The runs come from `../HAT_groin_hindcast_1967_2017.py` (or the three-way
runner in `"single"` mode); the observed table is
`../../../../../1-observations/wetdry_photo_positions/Change_from_wetdry_1967_D2_D12.csv`.

## The scripts in detail

Moved out of the script on 2026-10-01, when it was brought in line with
`scripts/STYLE.md`. Passages from its own header are kept word for word.

### HAT_groin_effect_comparison.py

A standalone publication-figure script. It compares a `no_groin` and a
`groin` 1967-2017 hindcast run against the real observed shoreline at
MULTIPLE checkpoint years (1997 mid-run, 2017 endpoint), isolating the
groin's (and its deterioration's) effect on modelled shoreline position over
time.

Generalised from the original single-endpoint (1967->1997) version to an
arbitrary list of checkpoint years, so the 1997 mid-run point is a genuine
validation target -- not just the final year -- which matters specifically
for the 1995->2003 deterioration ramp the extended run is testing.

**Observed change methodology -- important.** Observed change comes from the
validated WET/DRY change table (`Change_from_wetdry_1967_D2_D12.csv`, built by
`HAT_geometric_distance_sanity_check.py`) -- NOT the dune-line raw CSVs used in
earlier versions of this script. Two decisions behind that switch, both
validated empirically, not just asserted:

1. CASCADE's own x_s is defined (`barrier3d.py`) as x_t + LShoreface -- a
   shoreface-toe-based "shoreline" conceptually much closer to a water-line /
   MHW-type contour than to the dune/vegetation line. Confirmed: the 1967
   wet/dry line sits seaward of the 1967 dune line at every domain, exactly as
   physically expected.
2. The change table is internally self-consistent -- wetdry_1967 is the
   baseline, every other wetdry_YYYY is differenced against that SAME fixed
   reference throughout, never re-baselined per year (re-baselining to each
   year's own minimum domain was shown to silently erase real signal -- D12's
   own raw distance moved ~33 m between 1967 and 2017, which a per-year
   re-zeroing would hide).

Only the wetdry_* columns of that table are used -- the table also carries
duneline_* columns (referenced to the same wetdry_1967 baseline) for
diagnostic purposes, but mixing one into the observed series would
reintroduce the exact feature mismatch this switch fixed. It must be the
WETDRY-referenced table specifically, not the dune-line-referenced one.

Each year's observed change is already computed as raw_distance[year] -
raw_distance[1967], sharing ONE fixed absolute baseline throughout. It is
returned in the raw sign (+ = landward/erosion), which is what `_flip()`
expects downstream; a missing file or column gives None for that year
(omitted from the plots; the run continues).

**One shared reference frame, multiple checkpoints.**

- the 1967 shoreline (the model's own dune-line-based initial planform -- a
  DIFFERENT feature than the wet/dry observed series, but that is fine: every
  comparison here is CHANGE from each series' own 1967 reference, never
  absolute position across the two);
- for each year in `CHECKPOINT_YEARS` (e.g. 1997, 2017): OBSERVED (the 1967
  reference plus the observed wet/dry 1967->year change), modelled no groin
  (the model's own row at that calendar year), and modelled with groin.

Both runs share an identical 1967 initial condition, so either run's row 0
gives the same reference. The observed target is the 1967 reference minus the
observed change: the change is "+ = landward", and subtracting it flips it to
"+ = seaward", matching the model's `_flip()` convention, before it is added
onto the 1967 reference.

`_year_to_row()` maps a calendar year to its matrix row (row 0 = START_YEAR,
row t = START_YEAR + t), and returns None with a printed warning if the run
does not reach that year -- a guard against the END_YEAR-is-exclusive
convention silently pulling the wrong row.

`CHECKPOINT_YEARS` are the observed validation targets -- add or remove years;
everything else (figures, GIF markers, observed loading) follows, as long as a
matching `wetdry_<year>` column exists in the table. Each checkpoint year gets
its own alpha / marker / linestyle (`CHECKPOINT_STYLE`), so 1997 and 2017 stay
distinguishable in the overview while colour still encodes series TYPE
(observed black, no groin orange, groin red) consistently.

**The figures.** The overview puts ALL checkpoint years on one frame -- busier
than the per-year figure by design (up to 3 + 3 x len(CHECKPOINT_YEARS)
series). The per-checkpoint figure isolates the groin's effect AT ONE YEAR
against the real observed target.

`show_title=True` gives the full figure with title, subtitle and the grey
italic source footer (for standalone use, e.g. slides). `show_title=False` is
the "v2": no title, subtitle or footer -- just the plot -- for figures whose
caption sits in the surrounding text (e.g. a dissertation figure with a
numbered caption); pass a shorter height for it, e.g. (11, 5), since it needs
no room for the title block. Both share identical axis label sizes
(`AXIS_LABEL_FONTSIZE`), so they drop into the same document or deck without
looking mismatched.

The top x-axis is alongshore distance in km from D2 (the left edge of the
window), at the project's 500 m per domain. The per-checkpoint legend's left
edge is anchored at the 3 km mark on that axis rather than flush against the
right border, to keep it clear of the with-groin / no-groin lines further
right.

The groin label sits near the top of the FINAL displayed axis -- after the
later `ax.set_ylim(ax.get_ylim()[::-1])` inversion when `OCEAN_AT_BOTTOM` --
offset to the right of the dashed line so it never overlaps the line or the
x-axis. When `OCEAN_AT_BOTTOM`, the eventual top of the chart is the current
(pre-inversion) bottom of the y-range, so "near the top" means a small
fraction there.

**The GIF.** Shoreline position evolving year by year for the no-groin
(orange) and with-groin (red) runs, on the same 1967-anchored frame as the
static figures. The 1967 shoreline and every checkpoint's OBSERVED target are
thin static dashed references that stay put; only the two modelled lines
move. The year label per frame is exact integers (row t = START_YEAR + t) --
the previous `np.linspace(START_YEAR, END_YEAR, nt)` approach did not land on
clean years whenever nt - 1 did not divide the span evenly. The y-limits are
fixed for the whole animation up front from the full range of every series,
static and animated, so the axis never rescales frame to frame -- a rescaling
axis is the most common way animated plots end up jittery or misleading.

**Usage.** Set `RUN_NO_GROIN` / `RUN_GROIN` to the two run folder names and
`COMPARISON_SUBFOLDER` to a label for this comparison (so different
comparisons never overwrite each other), then run. It saves, under
`COMPARISON_OUTPUT_DIR`:

    groin_effect_overview.png        (all checkpoints, one figure)
    groin_effect_overview_v2.png     (same, no title/footer)
    groin_effect_{year}.png          (per checkpoint year)
    groin_effect_{year}_v2.png       (same, no title/footer)
    shoreline_evolution.gif

The footer still says "Obs: digitized dune-line offsets", from before the
switch to the wet/dry table.

**Paths.** As of 2026-10-01 the script's paths were fixed after the restyle
(see the folder history in git): the repo root is found by searching upward,
the wet/dry table is read from `1-observations/wetdry_photo_positions/`,
and figures go to `hard-structures/groin/3-hindcast/1-dipole-1967-2017/results/comparison/<COMPARISON_SUBFOLDER>/`.
Before that, `PROJECT_BASE_DIR` was `r"/"` and both paths named the
pre-reorganisation `scripts/groin/...` tree, so the script could not find its
table or write its figures. `RUN_DATA_DIR` still reads runs from
`output/raw_runs/`; the rig runner now writes to `output/calibration/groin_rig/`.
