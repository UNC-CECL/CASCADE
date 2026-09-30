# figure_making/model_output — figures drawn from finished runs

Cross-run and per-run figures: the hindcast result, the scenario grid, the GIS 11
relocation question, the plan-view animation, and the tool that redraws a
run's own figures. The manuscript copies go to `output/figures/5-results/`;
working copies to `output/comparisons/`. `superseded_20260914/` and
`superseded_20260918/` hold retired versions (see their WHY.md).

```
hindcast_final_figure_lowess.py  THE hindcast result: both periods against the LOWESS curve
scenario_grid.py                 every scenario x preset x period on one page
gis11_relocation_drown_figure.py why GIS 11 drowns at a 30 m relocation setback
planview_evolution_gif.py        animated plan view of one run: elevation and the road
rerender_run_figures.py          redraw a finished run's figures without re-running it
```

What not to trust: each figure reads the runs named in its CONFIG; a figure is
only as current as those runs. `rerender_run_figures.py` refuses to redraw if
the runner's figure conventions have moved (`check_conventions`).

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### gis11_relocation_drown_figure.py

Why GIS 11 drowns at a 30 m standard relocation setback and not at 77 m: the wet-cell criterion is a cliff.

From the script's original header:

```text
Why GIS 11 drowns at a 30 m standard relocation setback, and not at 77 m.

THE QUESTION THIS ANSWERED, AND WHAT WAS DECIDED
    Introducing a flat 30 m relocation target (`relocation_setback_m`) drowned
    NC-12 at GIS 11 in all eight 1984-2004 reloc-arm runs. Nothing drowned
    under the previous per-domain measured targets. The figure exists to settle
    whether 97 m is a physically implausible place to put a road, or whether it
    is merely one cell past a threshold.

    It is the latter, and that is the point of the third panel: the criterion
    takes whole-cell values only, and cells 7 and 8 are 0% wet against cell 9
    at 24%. There is no gradual approach to failure to read a margin off.

    DECIDED 2026-09-01: the standard is 20 m, and the whole matrix runs at it.
    A 20 m target lands the 1999 event at 87 m, which is cell 8, and all eight
    drownings go away -- confirmed by re-running the twelve reloc arms, where
    relocation counts moved only 26 -> 28 and no new drowning appeared
    anywhere. But 87 m clears the threshold by ONE CELL, so 20 m is not robust
    to different forcing; it is the dry side of a cliff, not a margin.

    THE UNDERLYING COUPLING WAS NOT FIXED. `_apply_relocation` still adds a
    surveyed displacement to a modelled position, so any future change to the
    emergent rule will move where the historical events land. Anchoring it to
    an absolute setback (initial measured setback + cumulative displacement)
    would decouple the two and let the standard be chosen on its own merits.
    That remains open.

THE CHAIN, WHICH IS NOT WHAT IT LOOKS LIKE
    The standard did NOT push the road into the bay directly. It raised GIS
    11's emergent relocation from 10 m to 30 m in 1993, so the setback stood at
    20 m rather than 0 m when the historical 1999 event fired. That event is
    stored as a DISPLACEMENT, and `_apply_relocation` adds it to whatever the
    model's current setback is:

        measured target:   0 m + 77 m  =  77 m
        30 m standard:    20 m + 77 m  =  97 m

    So a change to the emergent rule moved where a PRESCRIBED historical
    relocation lands. The measured displacements were surveyed against the real
    road; adding them to a modelled road inherits the model's drift.

WHAT DROWNS IT
    `bulldoze` drowns a road when more than 20% of the cells in the row
    BORDERING the road are at or below 0 m MHW. At GIS 11 that criterion has a
    cliff between 80 m and 90 m of setback -- 0% wet against 24% wet, one cell
    apart. 77 m clears it; 97 m does not.

    The road occupies `int(setback / 10)` and the cell behind it, and the
    bordering row is the one after that, so the axis is drawn in the same whole
    cells the model indexes in rather than in smooth metres.

Usage:
    python gis11_relocation_drown_figure.py [--out PATH]

Reads output/comparisons/relocation/standard_setback/GIS11_profiles.npz, the
per-domain extract taken before the superseded runs were deleted.
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). This file drew in
matplotlib's defaults until 2026-09-17 -- it never called apply_style().
```

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
The manuscript copy, with the other figures by subject (2026-09-18): the
same panels, the headline and its two lines moved to supporting/CAPTIONS.md.
Written only when --out is not given.
```

```text
ON since the layout was redrawn for the house-style column (2026-09-18).
Before that the labels, title and legend overlapped, and it was held off so
a broken figure could not land in output/figures/.
```

```text
The two candidate positions are 20 m apart on a 200 m axis, so a
label above each band overlaps its neighbour. Rotated inside the band
each label sits on the thing it names and cannot collide; at the
190 mm column only the distance fits, and the legend says which
target each colour is (2026-09-18).
```

```text
A 0% bar draws nothing, which reads as "not evaluated" rather than
"evaluated and dry". Every setback gets a visible stub.
```

```text
Three callouts within 20 m of each other on a 130 m axis. Stacked
horizontally they ran off the axis into panel (b) at the 190 mm column
(2026-09-18); rotated, each sits in its own bar's column above the limit
line, where the three dry cells leave the panel empty.
```

```text
SNAP TO THE CELL, not to the metre. bulldoze indexes the road at
int(setback / 10), so 77 m and 70 m are the SAME road position and
the same bar. Pointing the callout at its raw metre value drops it
between two bars and invites the reader to interpolate a criterion
that only ever takes whole-cell values.
```

```text
HOUSE STYLE (2026-09-18): apply_style() only, no local rcParams; the
headline and its two lines are the caption, not the canvas.
```

```text
300 dpi, the house savefig resolution. It was 170, which is a screen
export: the geometry was already right at 190 mm, the pixels were not.
```

<details><summary>Function notes (the original docstrings)</summary>

**`wet_fraction()`**

```text
Fraction of the bordering row at or below the drown threshold.

Mirrors `bulldoze`: the road starts at `int(setback / 10)`, runs
`ROAD_CELLS` cells, and the row checked is the one after it.

Args:
    grid: (cross_shore, alongshore) interior elevations in m MHW.
    setback_m: Road setback in metres.

Returns:
    (fraction_wet, bordering_row_index, bordering_row) or (nan, idx, None)
    when the bordering row is past the end of the interior.
```

**`draw_profile()`**

```text
One cross-shore profile with the candidate road positions on it.

Args:
    axis: Axes to draw on.
    grid: (cross_shore, alongshore) interior elevations in m MHW.
    setbacks: (setback_m, label, colour) tuples to mark.
    year_label: Calendar year for the panel title.
```

**`draw_cliff()`**

```text
Wet fraction of the bordering row against setback, with the 20% limit.

This is the panel that answers the question: the criterion is not a slope,
it is a step between two adjacent cells, so 77 m and 97 m sit either side
of it and 87 m -- where a 20 m standard would land -- sits on the safe side
by one cell.

Args:
    axis: Axes to draw on.
    grid: (cross_shore, alongshore) interior elevations in m MHW.
```

</details>

### hindcast_final_figure_lowess.py

The hindcast result figure: both periods' full-management runs against the LOWESS reference curve.

From the script's original header:

```text
The calibrated hindcast against the LOWESS reference curve — presentation figure.

THE hindcast result figure. It had a companion, HAT_hindcast_final_figure.py,
which drew the same two runs scored over D2-D89; on 2026-09-14 this figure took
on both scoring windows and the companion was retired to
model_output/superseded_20260914/. The three features below were what distinguished
the two, and are now simply what this figure has.

    SHARED Y AXIS
        Both periods on identical limits, so the panels can be read against
        each other rather than only against their own observations. That
        comparison is the point: period 2 sits almost entirely above period 1,
        which is the post-Isabel recovery and the nourishment era showing up as
        a whole-island shift in the rate, not as a local feature.

    THE LOWESS CURVE ONLY, NOT THE SPLICED TARGET
        The calibration target is not one curve: GIS 1-10 are raw per-domain
        means and D11 north is the 7-domain LOWESS. That splice is right for
        calibrating -- the raw means keep the short-wavelength signal the
        source/sink field has to answer for -- but it makes an awkward figure,
        because the eye reads a change of estimator as a change of coast.
        Here the LOWESS is drawn throughout.

        D1-D10 IS DRAWN DASHED. The project excludes the LOWESS there by
        convention (`skip_southern_domains = 10`), because the smoother is
        poorly constrained at the end of its range and Cape Point's
        attachment-detachment cycle is exactly the short-wavelength signal a
        3.5 km smoother destroys. Dashing it shows the data without implying it
        carries the same weight.

    BOTH SCORING WINDOWS, PRINTED PER PANEL
        D11-D89 is where the LOWESS curve and the calibration target are the
        SAME numbers, so that statistic is the one for the line actually drawn.
        D2-D89 is the project's canonical skill window -- rmse_interior_m_yr in
        run_index.csv, and what the groin fit and the source/sink convergence
        were ranked on -- and additionally takes in D2-D10, where the target is
        the raw spliced domain mean and the model is at its worst (RMSE ~1.1
        against ~0.45 north of it). Both are correct for what they measure and
        mixing them is the error, which is why both are on the panel and
        labelled. They are COMPUTED here, both of them: the D2-D89 pair used to
        be a literal quoted from the companion figure, and it sat twelve days
        stale after the 1984 run was remade on topography v2.

ALSO ADDED FOR PRESENTATION
    An observational spread band (+/- 1 SD of transect rates within each
    domain, from `unc_m_yr`'s parent transect table). A domain is a 500 m
    average over ~10 transects, and showing that spread makes clear how much
    real alongshore variability the domain mean hides -- and therefore which
    model-observation gaps are meaningful and which sit inside the noise.

    Place names along the top, from PHYSICAL_ZONES, so an audience can locate
    features without a separate map.

Usage:
    python hindcast_final_figure_lowess.py [--preset edgeBE|zeroBE]

ON THE 1996 -> 2010 -> 2024 CHAIN since 2026-09-18, nogroin arm, edgeBE by
default. calibBE is kept in PRESETS but is not solved on this chain.

Writes output/comparisons/hindcast_calibrated/hindcast_<preset>_lowess_reference.png
and, with PUBLISH, output/figures/5-results/hindcast_<preset>.png with its
caption in supporting/CAPTIONS.md. Both are the same house-style figure; the
title and note that were drawn on the canvas are the caption (2026-09-18).
```

Notes that were in the code:

```text
The manuscript copy goes with the other figures, by subject
(output/figures/README.md); the presentation version stays in OUT_DIR.
```

```text
ON since the layout was redrawn for the house-style column (2026-09-18).
Before that the labels, title and legend overlapped, and it was held off so
a broken figure could not land in output/figures/.
```

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. Applied at
MODULE level even though this file imports matplotlib inside a function --
apply_style() sets rcParams, so it only has to run before the figure is
built, and it was never called at all until 2026-09-17.
```

```text
THE CANONICAL CHAIN, 1996 -> 2010 -> 2024 (moved from 1984/2004 on
2026-09-18). Ends come from HATTERAS_PERIODS. The full-management arm is
road_bdm in 1996-2010, which schedules no fill, and road_bdm_nourish after.
GROIN: the matrix on this chain is nogroin only -- the groin is still being
fitted in the sweep -- so the figure draws the nogroin arm and says so.
```

```text
Non-overlapping display spans. PHYSICAL_ZONES overlaps at D9-D10 (Cape Point
and Buxton-Avon both claim them) and assign_physical_zone resolves that by
first match; a label strip has to pick one, so it picks the same one.
```

```text
A domain is 500 m holding ~10 transects. Plotted at the domain integer
they stack into a vertical column, which reads as one uncertain value
rather than as an alongshore gradient. Spread them across the domain in
transect order instead. The index is the LAST numeric field of the id --
the first is the CoastSat site number and grabbing it sorts by site.
```

```text
The matrix runs carry the offset token since the metres fix (2026-09-24);
without it this found only the archived ÷10 runs.
```

```text
Resolved rather than joined by hand: a hand-built path has no slot for
the arm component and reads an arm-scoped run as missing. The retired
companion carried a twin of this function, in
model_output/superseded_20260914/HAT_hindcast_final_figure.py.
```

```text
THE PRESET OWNS EVERY WORD THAT NAMES THE CONFIGURATION. Added 2026-09-14
alongside the edgeBE companion. --preset already existed here, but the
filename, the title and the footnote were hardcoded to calibBE, so any other
preset overwrote the calibrated figure with a page still calling itself
calibrated. `stem` must stay in step with the same table in
the retired companion's table (model_output/superseded_20260914/), which named
the same stems -- kept aligned so its output and this figure's still sort
together in the folder.
```

```text
Label heights are tuned PER FIGURE, not in the shared config, because
they depend on this figure's y range and on what is drawn where. The
defaults (groin 0.68, piers 0.76) put "Buxton Groin" through the D1-D10
transect scatter and "Avon Pier" through both shoreline curves.

Buxton Groin  0.78  below the "Buxton" town label (0.92 put the two on
top of each other at the 09-18 column width) and
still above the scatter in (a); the
"Buxton" town label is centred at D7.5 so its box
clears D5.5 horizontally.
Avon Pier     0.85  the "Avon" town span is ALSO centred on D26, and it
sits at 0.90, so the pier has to hang below it.
Rodanthe Pier 0.72  the "Rodanthe" village line is at D80 with its label
at 0.84; dropping the pier clears that, and D79 is
deeply erosional in both periods so the mid-panel is
empty there.
```

```text
WHAT THE COMPANION FIGURE REPORTS, COMPUTED, NOT REMEMBERED. The
footnote quotes the D2-D89 score so a reader can see why it differs
from this figure's D11-D89 one. It used to be a literal, and on
2026-09-14 it was found to be twelve days stale: written 09-02, it
still held the pre-v2 numbers after the 1984 run was remade on
topography v2 on 09-07. Same target, same model, same window as
the runner's own interior metric -- so it cannot drift from it.
```

```text
The D1-D10 transect scatter is drawn too, and at Cape Point it runs past
both curves; a limit set on the curves alone clipped it (09-18). It only
widens the axis as far as it reaches, without the label headroom the
curves get -- the place names sit well north of D1-D10.
```

```text
Asymmetric: the top needs room for the place-name strip, the bottom does
not, and a symmetric pad on a shared axis wastes a band of panel (b).
```

```text
HOUSE STYLE (2026-09-18). This block drew for a 15 in presentation page
-- 11-15 pt type, a suptitle and a paragraph of note on the canvas --
and when figsize() put it on the 190 mm column on 09-17 the ylabels
clipped and the note ran under the legend. Now: house type sizes, the
letter and period at the left above each panel, the two scoring windows
at the right above it (they used to sit in a box over the curves), the
place names once on (a), and the title and note in CAPTIONS.md.
```

```text
Communities, village centres, piers, groins and shoal zones, from the
shared layer every other cascade_pipeline figure uses -- so this
figure cannot disagree with the annotated run plots about where the
Buxton groin or the Avon pier is. Named on (a) only; (b) has the
same bands, and the same eight names twice is clutter.
```

```text
The D5-D7 (groin-reserved) and D1/D90 (locked) spans were dropped
from this figure: with the geographic annotation layer in place the
Buxton area carried four overlapping fills and read as clutter. Both
facts are stated in the caption instead.
```

```text
Individual transects over D1-D10. The LOWESS is dashed there because
the smoother is unreliable at the end of its range; the scatter is
the actual evidence, and it shows the Cape Point spread the smooth
curve cannot represent.
```

```text
BOTH SCORING WINDOWS, ON THE FIGURE. Added 2026-09-14. This figure
scores D11-D89, the span where the LOWESS curve IS the calibration
target; the project's canonical skill column (rmse_interior_m_yr,
recorded for every run in run_index.csv) is D2-D89, which also takes
in D2-D10, where the target is the raw spliced mean and the model is
at its worst. Printing only the narrow one reads as a better model
rather than a shorter ruler, and printing it in a separate figure is
what let the two drift apart.
```

```text
One key for both panels, in reading order: the observations, then the
model and its misfit per period, then the place layer.
```

```text
Figure-level, below both panels. Too many entries to sit inside an axis
without covering something.
```

<details><summary>Function notes (the original docstrings)</summary>

**`lowess_and_spread()`**

```text
(lowess_by_domain, sd_by_domain) for one period.

The LOWESS is the widest window's smoothed curve across ALL domains --
build_target_table would splice raw means over D1-D10, which is what this
figure is deliberately not doing.
```

</details>

### planview_evolution_gif.py

Animated plan view of the island through a run: elevation, the road, and every relocation.

From the script's original header:

```text
Animated plan view of the island through a run: elevation, not shoreline.

WHAT WAS MISSING, AND WHY THIS EXISTS
    Every matrix run already writes four GIFs and all four animate the
    SHORELINE -- a one-dimensional cross-shore position per domain, plotted
    against domain number. None of them shows elevation. The plan-view canvas
    that `init_planview` builds, the one the initialization figures use, was
    only ever drawn at t = 0.

    So there was no way to watch the barrier itself evolve: where the interior
    lowers, where overwash reaches, where the island narrows. This draws that
    canvas once per model year.

THE HISTORY IS ALREADY ON DISK -- no re-run is needed
    Barrier3D keeps `DomainTS`, one interior grid per year, and `x_s_TS`, the
    shoreline position per year, for all 120 padded domains. Both survive in
    the run's `.npz`, which is why those files are ~300 MB rather than the
    ~20 KB the shoreline matrix costs. This reads them and plots; it does not
    re-run anything.

TWO THINGS THAT HAVE TO BE RIGHT
    RAGGED GRIDS. `DomainTS[t]` is NOT a fixed shape -- Barrier3D trims water
    rows, so one domain-year is (174, 50) and another (200, 50), and the count
    changes as the island evolves. Each grid is put back on the full frame with
    `pad_cross_shore`, exactly as `load_domain_grids` does for the static
    figure, so the two are directly comparable.

    A FIXED FRAME. `build_canvas` sizes the canvas from the offsets it is
    given, so a per-year canvas is a per-year height and the animation would
    breathe. Offsets are taken against ONE reference for the whole run, and the
    axes get one y limit for every frame, so what moves in the GIF is the
    island and not the camera.

UNITS
    `DomainTS` is in decameters and `x_s_TS` likewise; the plan-view config
    carries `dam_to_m = 10.0` and a 10 m cell, so a decameter of shoreline
    movement is exactly one canvas row. Elevations are converted to meters to
    match the shared colorbar.

ROAD PLACEMENT IS VERIFIED AGAINST THE MODEL, NOT ASSERTED
    `bulldoze` indexes the roadway as `road_start = int(road_setback / 10)`
    rows into the interior grid and flattens those rows to `road_ele`. So the
    road it draws is checkable: the bulldozed rows are exactly constant
    alongshore in `DomainTS`, and in the 1984-2004 calibBE groin run
    `floor(setback / 10 m)` lands on that constant row in every domain-year
    tested. The overlay uses the same `road_rows()` the static figure uses, so
    the map, the static figure and the model all index the road identically.

    What the map does NOT show is the dune line the setback is measured from:
    `DomainTS` is the INTERIOR only, and Barrier3D keeps `DuneDomain`
    separately, seaward of interior row 0. A setback of 0 m therefore draws on
    the interior's seaward edge, which is one dune width (20 m) landward of the
    dune crest.

Usage:
    python planview_evolution_gif.py <run_directory> [--fps 3] [--out PATH]
```

Notes that were in the code:

```text
HOUSE STYLE, TYPEFACE ONLY: this script writes ANIMATION frames, and the
printed-width rule exists so 9 pt type is 9 pt on a page. A frame is never
printed, so its figsize is the frame size and is left alone (2026-09-17).
```

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
Defined the same way HAT_hindcast_1984_2024.py:255 defines it. It is not
exported by the site config, and the period table's paths are relative to it.
```

```text
The relocated colour from road_planview's own style, so a relocation marker
here and a relocated road bar in the static figure read as the same thing.
```

```text
A road the model has stopped managing is still a road on the ground, so it is
drawn -- greyed, at its last managed position, and never counted again.
```

```text
THE CANVAS FRAME IS THE ISLAND OFFSET, NOT x_s. This was wrong until
2026-08-31 and the error was visible: x_s is Barrier3D's own shoreline
coordinate, and because every domain is a separate Barrier3D it varies by
only ~0.55 km across the island. The REAL alongshore geometry lives in
the BRIE island-offset file and spans ~6.3 km. Using x_s as the frame
compressed the island's diagonal about elevenfold, so the road sat far
from the dune line it is measured against and the figure disagreed with
HAT_road_island_planview_1984.png, which is built from the offsets.

The offset is the STATIC frame; x_s supplies only its CHANGE, so the
island is placed where the initialization figures place it and still
migrates through the run.
```

```text
THE ROAD MOVES TOO, and it is the reason to watch this rather than the
shoreline GIFs: a setback is measured from the dune line, so a road that
never relocates still closes on the ocean as the barrier retreats. The
roadway manager keeps _road_setback_TS per domain per year; where a
domain carries no road the series is absent and its entry stays 0, which
road_rows() renders as NaN rather than as a road at the dune line.
```

```text
WHERE THE ROAD ACTUALLY IS, from this run rather than from a site-wide
constant. cascade._roadway_management_module is a per-domain mask of the
domains the roadway manager actually manages -- 55 of 90 in the 1984
calibBE run, matching that run's own road_management_summary.csv exactly.

HATTERAS_FIRST/LAST_ROAD_DOMAIN (9 and 90) is the REACH, not the road:
NC-12 is present 9-20, 32-67 and 84-90, with real gaps at 21-31 and
68-83. Drawing the whole reach put road through both gaps. Using the
run's own mask also means a scenario with roadway management off draws no
road at all, which is correct and which a constant cannot express.
```

```text
WHEN THE MODEL STOPS MANAGING A ROAD it returns from RoadwayManager.update
BEFORE writing that year's time series, so _road_setback_TS stays 0 for
every remaining year. Read literally that draws an abandoned road pinned
to the dune line for the rest of the run -- exactly where a drowned road
is not. _road_ele_TS is the reliable per-year flag: the manager stops the
moment the road elevation would go below 0 m MHW, and leaves 0 behind.
```

```text
An abandoned road keeps its last managed setback so it can be drawn
greyed where it actually sits, rather than snapping to the dune line.
```

```text
RELOCATIONS ARE EVENTS, not a state: _road_relocated_TS is a 0/1 flag
per domain per year, raised in the year the roadway manager moves the
road. Kept per year rather than accumulated so a frame shows what
happened THAT year and the running total can be built from it.
```

```text
A setback of ZERO means the road sits ON the dune line, not that there is no
road. road_planview.road_rows() cannot tell those apart -- it returns NaN for
any setback <= 0 -- so a road whose setback decayed to zero vanished from the
animation and reappeared if it later relocated seaward.

That is not rare and it is mostly NOT relocation: in the 1984-2004 calibBE
groin run, 31 of 90 domains cross zero during the run and only five ever
relocate (GIS 10, 11, 84, 85, 86). The rest are the dune line simply catching
up with a road that never moved -- GIS 9 decays 40 -> 30 -> 20 -> 10 -> 0 m
and stays there. Blanking it draws the road as absent exactly when it is most
exposed, which is backwards.

The road's real extent is a site fact, not something to infer from a setback:
hatteras_site_config says GIS 9-90 carry NC-12 and "Domains 1-8 (Cape Point)
have no road in the modelled span". So presence comes from that, and a zero
setback inside the road reach is nudged just above zero -- floor(eps / 10 m)
is 0, so the bar lands exactly on the dune line, which is where the road is.
This keeps road_rows() as the single implementation of the geometry rather
than reimplementing it here where it could drift.
```

```text
No legend here: the map's legend already names both bar colours, and a
second copy under the panel only competes with the year axis.
```

```text
The run name carries the period, so the calendar year is recoverable
without trusting an attribute that may not exist.
```

```text
Domains that have relocated at least once BY a given year, so the map can
keep showing where the road has already had to move rather than flashing
it for a single frame and losing it.
```

```text
THE CROSS-SHORE AXIS IS STRETCHED. Cells are square, so the distortion is
just the pixel aspect of the axes box; stating it means the reader can
take the shapes at face value instead of guessing at them.
```

```text
Static furniture, drawn once. The run name is provenance, not a title,
so it goes in the subtitle at reading weight rather than in bold above
the panel where it competed with the map.
```

```text
THE ROAD. Placed by the same road_rows() the static figure uses --
offset + floor(setback / cell) -- so the two agree by construction
rather than by a second implementation that could drift. A road the
model has given up on is drawn separately, greyed, at the last
position it was managed in, so abandonment reads as abandonment
rather than as a road that has snapped onto the shoreline.
```

```text
RELOCATIONS, marked in the year they happen and kept afterwards as a
faint tick. A relocation is the one thing in this animation the model
DECIDES rather than suffers, so it is drawn as an event marker above
the road rather than as a change of road colour, which would be
invisible on a 1-cell bar. An event that cannot move the road -- a
prescribed relocation setback of 0 m -- gets a hollow marker, because
counting those as retreat is the easiest way to overstate this figure.
```

```text
THE COUNTER, in the empty ocean at lower left where it competes with
nothing. Running and final totals both, because "how often did the
road have to move" is the question this animation exists to answer,
and a per-year flash never answers it.
```

```text
A SCALE BAR, because the two axes are at different scales and a
reader measuring the island off the tick labels alone will get the
alongshore distance wrong.
```

```text
Cells are an implementation detail; kilometres are what a reader
measures the island in.
```

```text
The colorbar is drawn once, outside the frame loop: plot_canvas would
otherwise add a new one on every frame and shrink the axes each time.
```

```text
One legend for the road, at figure level so a frame redraw cannot drop
it and it never lands on the barrier.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_period_start()`**

```text
The hindcast period this run belongs to, for picking its offset file.

Read from the run NAME first: HATTERAS_PERIODS is keyed on the period start
year, the run name carries it, and a cascade attribute may not.

Args:
    run_dir: The run directory, whose name carries the period.
    cascade: The loaded Cascade, used only as a fallback.

Returns:
    A key present in HATTERAS_PERIODS.

Raises:
    SystemExit: If no period can be identified.
```

**`RunHistory()`**

```text
Everything one run contributes to the animation, already per-year.

Attributes:
    grids: grids[t] is a padded-order list of metre-valued interior arrays
        on the full cross-shore frame.
    offsets: offsets[t] is the matching per-domain canvas row origin.
    setbacks: setbacks[t] is the per-domain drawable setback in metres,
        zeroed where no road should be drawn that year.
    managed: managed[t] is the per-domain mask of roads the model is still
        managing that year.
    relocations: relocations[t] is the per-domain relocation flag.
    rebuilds: rebuilds[t] is the per-domain dune-rebuild flag.
    prescribed_setback_m: per-domain relocation setback the run was given.
        Where this is 0 a "relocation" cannot move the road landward at
        all, so the event is a rebuild in place, not a retreat.
    road_domains: per-domain mask of domains carrying a road at all.
    start_year: calendar year of frame 0, or None.
```

**`load_history()`**

```text
Per-year grids, shoreline offsets and roadway state from a run's .npz.

Args:
    run_dir: A matrix run directory holding exactly one .npz.

Returns:
    A RunHistory.

Raises:
    SystemExit: If no .npz is present, or it holds no domain history.
```

**`_drawable_setbacks()`**

```text
Setbacks with roadless domains blanked and on-dune roads kept visible.

Args:
    setbacks_m: Padded per-domain setbacks in metres, one per domain.
    has_road: Per-domain boolean mask of the domains carrying a road.

Returns:
    A float array: 0.0 where there is no road, so road_rows() blanks it,
    and at least _ZERO_SETBACK_EPS_M where there is one so a road whose
    setback has decayed to zero still draws, on the dune line.
```

**`_managed_at()`**

```text
Whether the model was still managing this road in this model year.

The roadway manager writes a positive road elevation every year it runs and
returns without writing once it gives the road up, so a zero elevation is
the abandonment flag. Year 0 is the initial state, before any update.
```

**`build_relocation_tally()`**

```text
Splits relocation events by whether they win the road any clearance.

A relocation resets the setback to the domain's PRESCRIBED relocation
setback, and CASCADE has no separate parameter for that: cascade_groin.py
re-assigns `road_relocation_setback = road_setback`, the domain's STARTING
setback, every year. So a domain that starts with the road on the dune line
can only ever be relocated back onto the dune line.

That is the case at GIS 85 and 86, whose measured 1984 setback floored to
0 m. Each event still drags the road one 10 m cell landward -- it rides the
dune toe -- but it ends the year with the same zero clearance it began
with, so the next cell of retreat re-fires it. The counts are exact: GIS 85
retreats 72.7 m (7.3 cells) and relocates 7 times, GIS 86 retreats 59.9 m
(6.0 cells) and relocates 6 times. One relocation per cell.

A domain with clearance behaves completely differently. GIS 9 holds 40 m,
absorbs 38 m of retreat over the whole run and never triggers at all.

So the split is NOT moved against not-moved -- every event moves the road.
It is relocated-with-clearance against re-pinned-at-the-dune-line, and
reporting one total for both would say the road retreated eighteen times
when it bought itself room five times.

Args:
    history: A RunHistory.

Returns:
    (moves_by_year, pinned_by_year, cumulative_moves, cumulative_pinned)
    where the first two are per-domain boolean arrays per year and the last
    two are integer running totals per year.
```

**`draw_timeline()`**

```text
Draws the static relocation timeline strip below the map.

A plan view has no time axis, so animating one leaves the reader with no
sense of WHEN anything happened -- only of what is on screen now. The strip
carries the whole run at once and the cursor says where in it this frame
sits, which is what turns a loop into a record.

Args:
    axis: Axes for the strip.
    moves_by_year: Per-year per-domain masks of relocations that move.
    pinned_by_year: Per-year per-domain masks of 0 m rebuilds in place.
    start_year: Calendar year of frame 0, or None.

Returns:
    The cursor Line2D, to be moved each frame.
```

</details>

### rerender_run_figures.py

Redraw a finished run's figures without re-running the model, from its saved shoreline matrix.

From the script's original header:

```text
Redraw a finished run's figures WITHOUT re-running the model.

WHY THIS EXISTS
    The per-run figures are built during a run, so a change to the plotting
    package only reaches a run that is executed again. After the 2026-09-10
    restyle that left 164 run folders holding figures in the previous look --
    two LOWESS curves, the wave height and the SLR rate in the title, a 22 in
    canvas. Re-running them would be 5-11 hours of model time and would rewrite
    240 MB of archive per run to change a picture.

    Every input those figures need is already saved beside them:

        <run>_shoreline_matrix.npy    (annual_states, padded_domains), metres
        <run>_run_metadata.json       period, wave height, BE state, run name

    so this reads those two, recomputes the plotted rate with the same
    function the run used, rebuilds the CoastSat series from the same CSVs,
    and calls the same plotting entry points. The .npz is never opened.

WHAT IT WILL AND WILL NOT REPRODUCE
    The two rate PNGs are exact: `compute_lrr` is deterministic on the saved
    matrix, and the CoastSat side is read from files on disk.

    The GIFs are exact EXCEPT for the roadway-relocation markers. Those come
    off each RoadwayManager's `_road_relocated_TS`, which lives only in the
    .npz, so `--gifs` draws no relocation markers unless `--open-npz` is
    given. A run whose GIFs carry markers is therefore left alone by default:
    the script detects them from the run's road-management table and SKIPS
    the GIFs for that run, rather than quietly dropping the markers. `--open-npz`
    loads the archive for those runs and keeps them.

WHAT IT NEVER TOUCHES
    The .npz, the .npy matrix, every CSV and TXT in the run folder, and
    output/raw_runs/run_index.csv. It only overwrites image files, and only
    the ones it can rebuild.

USAGE
    python rerender_run_figures.py --dry-run
    python rerender_run_figures.py --arm matrix/1984_2004/calibBE
    python rerender_run_figures.py --match "*calibBE*groin" --gifs
    python rerender_run_figures.py --run-dir output/raw_runs/.../HAT_...
    python rerender_run_figures.py --arm matrix --ylim=-10,10 --ylim-real=-7.5,7.5
    python rerender_run_figures.py --arm sensitivity --lowess-only
```

Notes that were in the code:

```text
These four MUST match section 8/9 of HAT_hindcast_1984_2024.py. They are
restated rather than imported because importing that module runs a hindcast.
The assertion in `check_conventions` catches them drifting apart.
```

```text
The two windows the runner added on 2026-09-11 (section 8 of
HAT_hindcast_1984_2024.py). Missing here until 2026-09-27, so a 1996 or
2010 run was redrawn with NO CoastSat curve: build_coastsat_series found
no active dataset and the overlay drew nothing. Keep in step with the
runner's list.
```

```text
The full-record rate the PROJECTED position change is built from (the
advisor's target, 2026-09-19): the 1996-2024 LRR carried onto a run window.
```

```text
Resolved, not joined: this OVERWRITES the figure the run already has,
so it has to land wherever that figure currently lives -- the new
figures/ subfolder, or the old flat name if the run has not moved.
```

```text
REFRESH ONLY. A run whose GIFs were disabled has none on disk, and
drawing four for it now would add files the run never produced and
change what the tree contains. Skip it.
```

<details><summary>Function notes (the original docstrings)</summary>

**`position_change_jobs()`**

```text
(reference, scaled cs_series, observed legend, caption phrase) for the
two observed references a run's position change is drawn against.

Named by the window the rate was FITTED on (the 09-21 vocabulary): the
run window's own LRR x span is TOTAL change; the 1996-2024 LRR x span is
PROJECTED change. Both are LRR x the run's span in years.
```

**`check_conventions()`**

```text
Fail loudly if the hindcast's figure conventions have moved.

A re-render that silently used a different estimator or LOWESS window than
the run would put two incompatible curves in one folder, which is exactly
the failure this script exists to clean up.
```

**`has_relocation_markers()`**

```text
Whether this run's GIFs carry roadway-relocation markers.

Read from the run's road-management table, which the run writes beside
the figures; the per-year flags themselves live only in the .npz.
Resolved rather than joined: the table is `tables/road_management.csv`
in the current layout and `road_management_summary.csv` in the old one.
```

</details>

### scenario_grid.py

Every management scenario, every solved preset, both periods, on one page against the CoastSat target.

From the script's original header:

```text
Every scenario, every preset, both periods, on one page.

THE FIGURE
    Periods as rows and source/sink presets as columns, left to right in
    order of increasing correction. Each panel draws one line per management
    scenario against the section 8 CoastSat target.

    NEITHER DIMENSION IS FIXED. The rows come from PERIOD_STARTS through
    HATTERAS_PERIODS -- the canonical chain is 1996-2010 and 2010-2024 since
    2026-09-17 -- and a preset only gets a column if be_rates() has it solved
    for EVERY period drawn. calibBE is solved for 1984 and 2004 only, so on
    the current chain it is dropped rather than drawn as a column that can
    never be filled.

    A cell with no run on disk draws nothing, so the script PRINTS the
    missing (period, preset, scenario) combinations: a sparse grid should
    read as runs not yet done, not as a result.

    So the two contrasts read on different axes: scanning ACROSS a row shows
    what the source/sink term does, and the spread WITHIN a panel shows what
    management does. The target is the same heavy black line in every panel,
    which is what makes the across-row read a skill comparison rather than
    just a shape comparison.

WHAT IS ON THE Y AXIS
    Shoreline change rate, m/yr, (+) seaward -- read from each run's own
    `*_shoreline_change_rate.csv`. That file is written by the pipeline from
    the same array section 12 scores, so this figure and the reported skill
    numbers cannot disagree about what a run did.

SHARED SCALES
    x is shared down each column: GIS domain 1-90, the whole island.

    y is shared across ALL SIX PANELS by default, so a change rate has the
    same height everywhere on the page and the two periods can be compared
    directly by eye. That is the whole point of the figure: if the axes
    differed, a period-2 line that looked steeper than a period-1 line might
    only be a different scale, and every amplitude read would need a glance
    at the tick labels first.

    The cost is small here, and was measured rather than assumed. Period 1
    spans 6.40 m/yr across every run and the target, period 2 spans 8.23, and
    the two together span 8.33 -- so on a common axis period 1 still occupies
    77% of the height. There is no meaningful squashing to trade away.
    `--y-per-row` restores an independent range per period for the case where
    one period's detail has to be read closely.

COLOUR
    A single-hue sequential ramp ordered by management intensity: natural
    (lightest) through to full_management (darkest). The ordering is in the
    colour, so "more management" reads as "darker" without consulting the
    legend, and a single hue stays legible under the common colour-vision
    deficiencies. The observed target is black and heavier than any model
    line, so it never competes with a scenario for attention.

    Relocation arms are the SAME colour as their non-reloc twin, dashed. They
    sit almost exactly on top of it -- the two differ in the fifth decimal of
    mean bias -- and drawing them as a distinct colour would imply a
    separation that is not there. Dashed-over-solid shows the overlap
    honestly, and any real divergence would immediately stand out.

Usage:
    python scripts/figure_making/model_output/scenario_grid.py
    python scripts/figure_making/model_output/scenario_grid.py --no-reloc --out FIG.png
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). This file drew in
matplotlib's defaults until 2026-09-17 -- it never called apply_style().
```

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
HAT_run_all lives in scripts/hatteras_ms/, named outright since this file
moved to scripts/figure_making/model_output/ (2026-09-18); it used to be
found as this file's grandparent.
```

```text
THE RUN DRIVER'S OWN GUARD. Which (period, scenario) pairs are distinct
runs is decided by HAT_run_all.scenario_applies -- full_no_fill only exists
where a fill is actually scheduled -- and this figure asks it rather than
keeping a second copy of the rule that could disagree (2026-09-17).
Importing is safe: HAT_run_all does its work under a __main__ guard.
```

```text
The manuscript copy, with the other figures by subject (2026-09-18). Written
only from a default run: an --out or any flagged variant is a working figure
and must not overwrite it.
```

```text
THE CANONICAL CHAIN, 1996 -> 2010 -> 2024 (Hannah, 2026-09-17). Ends come
from HATTERAS_PERIODS, so changing PERIOD_STARTS moves the whole figure.
```

```text
ONLY PRESETS SOLVED FOR EVERY PERIOD DRAWN. calibBE is solved for 1984 and
2004 only, so on the new chain it would be a column that can never be
filled -- not a gap in the runs but a preset that does not exist for those
windows. Asking be_rates() is what decides, so this cannot go stale.
```

```text
Section 8's settings, matching the runner, so the target drawn here is the
curve the runs were scored against rather than a second opinion.
```

```text
WHICH MODEL COLUMN THESE PANELS DRAW. The observed curve on every panel
is a CoastSat LRR -- a per-transect OLS slope through the period -- so
the model side is read from lrr_m_yr, the run's matching OLS slope
through its annual states, rather than change_rate_m_yr, which is a net
displacement over a span. Set to "change_rate_m_yr" to redraw a
pre-2026-08-22 version of this figure.
```

```text
Starts at 0.35, not 0.0: the pale end of any sequential map disappears
against white, and the lightest scenario still has to be readable.
```

```text
Resolved rather than joined: runs forced off the calibration wave
climate sit under an arm component this join had no slot for. The
default arm is the calibration one, which is what this grid draws.
```

```text
lrr_m_yr where the run has it, change_rate_m_yr otherwise.
The target on these axes is a CoastSat LRR, so the model side
has to be one too; a run written before the column existed
still plots, and says so.
```

```text
The end year comes from HATTERAS_PERIODS, not from start + 20: the older
pair happened to be 20-year windows, the canonical chain is 1996-2010
and 2010-2024, both 14 (2026-09-17).
```

```text
bands and lines only: six panels cannot each carry the
eight annotation names at this width
```

```text
Applied after every panel is drawn, because the shared case needs the
pooled range of the whole figure and cannot be set row by row.
```

```text
The strip the legend needs depends on how many scenarios actually had
runs, so it is measured rather than fixed: with only full_management
on disk the key is one row, not three.
```

```text
The title and the y-axis note used to be drawn here. Both are caption
material under figure_making/STYLE.md, and on a 190 mm six-panel grid they
were also the two widest things on the page (2026-09-17).
WHAT IS MISSING, SAID OUT LOUD. A cell with no run draws nothing, and a
near-empty grid looks like a result rather than an absence of runs
(2026-09-17: the move to 1996/2010 left most scenarios unrun).
```

```text
A cell only counts if the driver says it is a DISTINCT run: a scenario
that collapses onto another in this period is not a gap, and reporting
it as one would leave the figure permanently claiming missing work.
```

```text
NO TIGHT BBOX: it trims to content, and with the labels hanging outside
the axes this saved at 14.20 in wide however figsize() was set -- nearly
double the column. tight_layout above already reserves the margins.
```

<details><summary>Function notes (the original docstrings)</summary>

**`classify()`**

```text
Maps a run directory name to (scenario, relocations).

Reads the switch tokens rather than matching whole names, because the
token set is not the same in both periods -- period 2 carries a
nourishment token that period 1 has no reason to.

Returns:
    (scenario_key, reloc_bool), or (None, None) for a run this figure
    does not draw (anything with the groin attached).
```

**`load_runs()`**

```text
Every nogroin run for one period/preset, keyed by (scenario, reloc).

Returns:
    {(scenario, reloc): Series indexed by GIS domain}, empty if the
    preset directory does not exist.
```

</details>
