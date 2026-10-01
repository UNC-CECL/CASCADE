# sensitivity_analysis - varying one forcing at a time

Drives the hindcast repeatedly with one parameter moved off its calibrated
value, and reads the result.

```
hindcast_sensitivity.py          the sweep driver: one hindcast run per cell
plot_sensitivity.py              the sweep's figures: skill per axis, alongshore, road outcomes
plot_sensitivity_pair.py         two cells of one axis, stacked, when their curves overlap
plot_sensitivity_zoom.py         one axis over a value range, beside the full set
plot_hs_experiment.py            the wave-height experiment specifically (source/sink by zone)
```

The driver's manifest goes to `output/calibration/sensitivity/`
(`sensitivity_<start>.jsonl`); the figures go beside the cells, under
`output/raw_runs/sensitivity/figures/<start>_<end>_<preset>/` (since
2026-09-28). `output/calibration/sensitivity/figures/README.md` is one of the
three decision records in the output tree.

**A cell that moves a forcing earns a name token**, so it lands in its own
directory. Without one it would derive the same name as the matrix run it is
being compared against, and the last to finish would wear the production name.
Cells file under `output/raw_runs/sensitivity/<axis>/<period>/<preset>/` (the
driver sets `HAT_RUN_KIND=sensitivity`; the axis is read off the token), and
the index keys a cell on (run_name, kind, tag) so it and its matrix baseline
are two rows. No model state is written for a cell.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### hindcast_sensitivity.py

Run the hindcast once per sweep cell, with one setting moved off its calibrated value.

From the script's original header:

```text
Parameter sensitivity for the Hatteras hindcast, either period.

WHY THIS REPLACED THE SCRIPTS THAT WERE HERE
    `HAT_waveSensitivity_1984_2004.py` and
    `HAT_waveheight_Sensitivity_1984_2004.py` both built the model by hand.
    Between them they had: a SyntaxError that made one of them unrunnable,
    input paths that no longer resolved (`RoadSetback_1984.csv` under
    `raw_offset/`, `storms_1984_2004_base.npy` under `hindcast_storms/` --
    neither exists; the only copies live under `old_method_offset/` and
    `old/testing_storms/`), a direct `Cascade(...)` call against the STOCK
    class rather than the sandbox one every hindcast run uses, and a comment
    pinning them to "HAT_hindcast_1984_2024_old version.py". Anything they
    produced would have described a different model on legacy forcing.

    This drives the hindcast instead of reimplementing it. Each cell is one
    ordinary run of `HAT_hindcast_1984_2024.py` with one setting moved through
    the environment, exactly the way `HAT_run_all.py` drives the matrix. There
    is no second copy of the model setup here, so the sweep inherits the period
    table, the measured setbacks, the groin and the LRR estimator
    automatically, and cannot drift from them.

    Nothing in the notebook or the headless .py knows this file exists. Those
    two were touched once, on 2026-09-01, to make the swept values settings
    rather than literals and to give an off-default value a run-name token;
    that is the whole of their involvement, and both changes stand on their own.

BOTH PERIODS COME FREE
    The period is `--start-year`, and everything period-specific is already in
    HATTERAS_PERIODS. Nothing in this file knows what year it is.

CELLS CANNOT COLLIDE WITH THE MATRIX, BUT THE TWO AXES DO IT DIFFERENTLY
    Without some separator every cell would derive the SAME name as the matrix
    run beside it, and the last one to finish would be left wearing the
    production name -- the failure `output/calibration/groin/README.md` documents for
    the rig sweep. Two mechanisms prevent it, and which one applies depends on
    the axis:

      * `relocation_setback` earns a NAME token (`rset40`) from
        `cascade_pipeline.hindcast`, so the cell is a separate directory beside
        the matrix run and a separate row in run_index.csv.

      * The four WAVE axes earn one too (`waveHs1p2`), again since 2026-09-16.
        Between 09-01 and 09-16 the runner filed a wave cell by a forcing ARM
        under the matrix run's own name instead, which fanned one sweep out
        into twelve top-level folders; the purpose layout put the token back.

    Every cell is filed under raw_runs/sensitivity/<axis>/<period>/<preset>/
    by run_registry: this driver sets HAT_RUN_KIND=sensitivity and the runner
    derives the axis from the name's trailing token. The index key is
    (run_name, kind, tag), so a cell and its baseline -- which differ only in
    the token -- are two rows.

WHAT IS SWEPT
    Five axes, one at a time, each around its calibration value. Only the swept
    parameter moves, so every cell differs from its baseline in exactly one
    thing. Four are the wave climate; the fifth is the relocation target, which
    is not wave physics but is the newest and least settled number in the
    model -- `hat_run.yaml` records that 20 m clears the GIS 11 drowning
    threshold by ONE 10 m cell.

Usage:
    python hindcast_sensitivity.py --start-year 1984 --param wave_height
    python hindcast_sensitivity.py --start-year 1984 --param all --dry-run
    python hindcast_sensitivity.py --start-year 2004 --param wave_height \
        --values 1.5,2.0,2.5
```

Notes that were in the code:

```text
The "each domain relocates to its own measured offset" case. The environment
carries strings, and an EMPTY string reads as unset in HAT_hindcast_config
(`if raw != ""`), so this case has to be spelled: "" would silently leave the
cell at the default and the result would be indistinguishable from one.
```

```text
The swept axes. `setting` is the HAT_hindcast_config field, which fixes both
the environment name (HAT_ + upper case) and the default the name token is
measured against -- so a value equal to the default produces an untokened
name, which would be the matrix run. Those cells are skipped, not run.

`arm` is a per-axis callable (start_year, end_year) -> environment overrides,
for an axis that would be inert or misnamed in the default arm. A callable
rather than a dict because whether the override is right DEPENDS ON THE
PERIOD -- see relocation_arm.
```

```text
The four wave axes are centred on option A (Hs 2.0, Tp 7.5, asymmetry
0.6, high-angle 0.5; the defaults since 2026-09-27) and were set with
Hannah on 2026-09-28. Each keeps one cell past where the 09-24..09-27
wave studies (raw_runs/experiments/wave-climate/) saw the model break, so
the record shows the edge rather than stopping short of it:
Hs 0.75 and Tp 12  -- the barrier drowned at Hs 0.75 / Tp 10 and at
Tp 12 with Hs 1-1.25 (both at other settings)
asymmetry 0.5      -- below it the net drift reverses
high-angle 0.55    -- past ~0.5 most of the coast turns anti-diffusive
The USGS hindcast mean Hs for this coast is ~1.2-1.3 m; 2.0 is a fitted
value, and 2.0-2.5 was a flat optimum in the 09-27 Hs check.
```

```text
NOT wave physics, and swept for a different reason. 20 m was decided on
2026-09-01 because 30 m drowned NC-12 at GIS 11 in all eight 1984-2004
reloc arms, and hat_run.yaml records that 20 clears that threshold by ONE
10 m cell. A number chosen with one cell of margin needs its margin
measured rather than asserted. 0 is the GIS 85/86 ratchet in its purest
form; MEASURED is the pre-2026-08-31 behaviour as it actually was.

RUN WITH THE HISTORICAL EVENTS ON. The drowning being asked about is a
prescribed 1999 DISPLACEMENT landing on top of the emergent target, so
with the events off this axis would measure something else and report a
reassuring flat line. Period 1 only: no event falls in 2004-2024, and
there the axis moves emergent relocations alone.
```

```text
The child prints a few non-ASCII characters (arrows, en dashes);
without this its stdout is cp1252 on Windows and the capture below
raised UnicodeDecodeError on every cell (2026-09-16).
```

```text
The axis's own arm first, then the swept value. In this order, so an axis
can move the arm it is measured in but cannot overwrite the one value
that defines the cell.
```

```text
The run prints where it landed; echoing it makes the sweep log a map
from parameter value to directory without a second convention.
```

<details><summary>Function notes (the original docstrings)</summary>

**`relocation_arm()`**

```text
The arm the relocation-target axis has to be measured in, for a period.

Period 1 holds the 1989 and 1999 events, and the drowning this axis asks
about is a prescribed DISPLACEMENT landing on top of the emergent target --
with the events off, the axis would measure something else and report a
reassuring flat line. So it forces them on there.

Period 2 holds no RelocationEvent (only the 2022 Jug Handle BridgeEvent),
so forcing them on would add a `reloc` token claiming history the period
does not have, and would name the cells for a baseline the matrix has never
run. There the axis moves EMERGENT relocations alone, which is a real but
different question, and the ordinary arm is the right one.

Read from HATTERAS_ROAD_EVENTS rather than by testing the year, so adding
an event to that table brings this with it.

Args:
    start_year: Period start.
    end_year: Period end.

Returns:
    Environment overrides for the arm, possibly empty.
```

**`normalise()`**

```text
One swept value as the run will see it: a float, or None for `measured`.

Comparison against the default has to happen on this, not on what was
typed. 20, "20" and 20.0 are one cell and MEASURED and None are one cell,
but `20 == "20"` and `"measured" == None` are both False, and either would
run a duplicate of the baseline under a name that hides which one it is.
```

**`build_environment()`**

```text
The environment one sweep cell runs under.

Mirrors `HAT_run_all.run_once`, including HAT_IGNORE_SETTINGS: without it
the settings this does not set would come from whatever experiment was last
left in `hat_run.yaml`, and every cell would silently inherit it. Stray
HAT_* variables in the calling shell are dropped for the same reason.

Args:
    start_year: Period start, a key of HATTERAS_PERIODS.
    sweep: One entry of SWEEPS.
    value: The value for this cell.
    args: Parsed CLI arguments.

Returns:
    The environment dict for subprocess.
```

**`cells_for()`**

```text
The values of a sweep that are worth running.

A value equal to the calibration default produces no name token, so the run
would derive -- and collide with -- the matrix run's own directory. That
cell is not a sensitivity result anyway: it IS the baseline, and reading it
from there is both free and more honest than re-running it under a name
that hides which one it is.

Args:
    sweep: One entry of SWEEPS.

Returns:
    (runnable_values, skipped_default) where skipped_default is the
    calibration value if it appears in the sweep, else None.
```

**`run_cell()`**

```text
Runs one sweep cell and returns its outcome.

Args:
    start_year: Period start year.
    sweep: One entry of SWEEPS.
    value: The parameter value for this cell.
    args: Parsed CLI arguments.

Returns:
    A dict recording the cell, suitable for the manifest.
```

</details>

### plot_hs_experiment.py

Where a higher Hs changes the source/sink correction, zone by zone (Hs 2.5 vs 3.0).

From the script's original header:

```text
Where a higher Hs changes the source/sink correction, zone by zone.

WHAT THIS ANSWERS
    The totals say the required correction field shrinks ~6% at Hs 3.0. They
    also hide the finding: the reaches move in OPPOSITE directions, and the
    ones that improve are not the ones carrying most of the correction. A
    single number for the island would report a modest win and conceal that
    the mid-island got worse.

WHY DIVERGING, AND WHY ORDERED SOUTH TO NORTH
    The quantity is a signed change either side of "no difference", which is a
    polarity encoding: two hues with a neutral midpoint, never a sequential
    ramp. And the zones are laid out by their domain range rather than sorted
    by value, so the panel reads as an alongshore profile -- which is what
    exposes that the improvement is concentrated at one end of the island.

INPUT
    The two pass-0 calibrations under output/calibration/hs/, produced by
    be_zone_residual_fit.py with HAT_BE_OUTPUT_DIR redirected. Nothing
    here reads or writes the production calibration.
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). This file drew in
matplotlib's defaults until 2026-09-17 -- it never called apply_style().
```

```text
THE ENCODING IS NOT FIXED, so it is detected rather than assumed. The
analysis writes through whatever encoding stdout has: run at a Windows
console it emits cp1252 (the case be_apply_fit_to_config.py's RATES_ENCODING documents),
run with stdout redirected to a file it emits UTF-8. Assuming either one
turns the en-dash in "Buxton-Avon Transition" into mojibake in the zone
labels -- a replacement character one way, "a-EUR-quote" the other.
```

```text
A diverging pair with a neutral midpoint: less correction needed is the good
direction and gets the cool hue, more correction the warm one. Deliberately
NOT the sensitivity figures' sequential ramp -- that encodes magnitude along
one hue, and this quantity has a sign.
```

```text
The en-dash in "Buxton-Avon Transition" is ALREADY LOST in the file: the
analysis wrote a literal U+FFFD because its own output encoding could not
represent the character. No read encoding recovers it, so it is repaired
here for display. Fixing it at the source would mean the analysis writing
UTF-8 explicitly, which is a change to a script this experiment is
deliberately not modifying.
```

```text
NOT `first`/`last`: those are DataFrame methods, so a column of either
name is reachable by [] but NOT by attribute -- row.first silently
returns a bound method, which is how the zone labels came out as
"<bound method NDFrame.first of zone ...>" on the first draft.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load()`**

```text
One pass-0 calibration, indexed by domain.

Tries each encoding in turn and takes the first that decodes cleanly. A
wrong guess does not raise -- cp1252 decodes any byte -- so UTF-8 is tried
first and only a genuine decode failure falls through to it.
```

</details>

### plot_sensitivity.py

Figures for the hindcast sensitivity sweep: skill per axis, alongshore rates per cell, road outcomes.

From the script's original header:

```text
Figures for the hindcast parameter sensitivity sweep.

WHY THIS REPLACED plot_sensitivity_vs_coastsat.py
    That script read a "timestamped session folder" produced by the two
    hand-built sensitivity scripts that have since been deleted, and it needed
    SESSION_DIR and two CoastSat CSV paths retyped at the top of the file before
    every use. Neither the folder layout nor those scripts exist any more.

    This reads the run registry instead. `hindcast_sensitivity.py` writes a
    manifest line per cell recording the parameter, its value and the directory
    the run landed in; `run_index.csv` carries the skill metrics; each run
    directory carries its own per-domain rates. Nothing is retyped and nothing
    is passed between the two scripts except the manifest.

THE OBSERVED LAYER IS THE TARGET, NOTHING ELSE
    The target table comes from `build_target_table`, the production scoring
    path, so the curve drawn is exactly what the RMSE is computed against. It
    is a HYBRID: GIS 1..skip_southern_domains are raw per-domain means (LOWESS
    is suppressed near Oregon Inlet, where boundary effects dominate and
    smoothing would hide the gradient), and the rest is the 7-domain LOWESS.
    The alongshore panels draw that one curve plus the transect dots over
    D1..skip only, where the target is not a LOWESS. The full transect scatter
    and the unsmoothed domain-mean line were dropped 2026-09-29: over D11-90
    they duplicated the target at higher noise and buried the model curves.

INTERIOR METRICS, NOT ISLAND-WIDE
    Every skill number here is the interior one (GIS 2-89). The two end domains
    carry imposed background-erosion rates under edgeBE and calibBE alike, so
    island-wide skill partly scores the boundary condition rather than the model.

Usage:
    python plot_sensitivity.py --start-year 1984 --preset edgeBE
    python plot_sensitivity.py --start-year 2004 --preset edgeBE
    python plot_sensitivity.py --start-year 1984 --circularity
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). Missed by the first
sweep because its figsize is computed, not a literal (2026-09-17).
```

```text
The base point is the matrix run's value, read from field_default. The
2026-09-28 sweep is centred on option A (Hs 2.0 / Tp 7.5 / asym 0.6 /
ahf 0.5), the defaults since 2026-09-27. Between those two dates the wave
fields were pinned here to the /10-offset climate (2.5 / 8 / 0.7 / 0.1),
whose sweep was archived on 2026-09-24.
```

```text
Reading order, most informative first. Alphabetical put the relocation axis
-- the one that barely moves shoreline skill -- in the leftmost panel and the
first file, which is the opposite of how these should be read. Used for both
the panel columns and the numeric filename prefixes, so the figure and a
directory listing agree.
```

```text
Beside the cells they describe (Hannah, 2026-09-28). Until then they went to
OUT_ROOT/figures, where the pre-2026-09-24 sweeps' directories still are.
```

```text
Where the runner reads them (COASTSAT_BASE_DIR in HAT_hindcast_1984_2024.py).
The scripts/input_prep/5-scr/CoastSat path this held until 2026-09-16 no
longer exists; the loader only WARNED, so every observed layer was empty.
Resolved through hat_observed_rates.py (2026-09-18), not typed.
```

```text
Must match section 8.1 of HAT_hindcast_1984_2024.py. These are the LowessConfig
defaults, so the two agree by construction rather than by copying -- but the
run metadata records the target it actually used, and check_target_matches()
below asserts against that rather than trusting this line.
```

```text
The model curves are a SEQUENTIAL ramp: the swept values are ordered
magnitudes, not categories, so one hue light-to-dark is the correct encoding
and a categorical cycle would imply the values are unrelated. Warm on purpose
-- the observed layer owns the blues (#5BA3C9 scatter, #6BAED6 / #08519C
LOWESS), so a blue model ramp would collide with the thing it is measured
against.
```

```text
The current setting, on the alongshore panels only. It needs a hue that is in
NEITHER family already on those axes -- the warm model ramp or the blue
observed layer -- because it is the one curve a reader goes looking for.
Black was too close to the navy target and, worse, the legend never said
which value it was: the calibration cell is skipped by the sweep (it IS the
baseline), so it never appears on the colourbar either.
```

```text
The manifest stores the run directory as an ABSOLUTE path, written
when the sweep ran, so any later move of output/raw_runs makes it
stale -- the 2026-09-10 tree change moved every sweep run under
sweeps/<family>/. The run NAME is stable, so re-resolve from it and
fall back to what was recorded.
```

```text
A 2026-09-01..16 cell filed by arm under the matrix run's own
name. Those were re-run with the token and deleted; one still
in a manifest is skipped rather than drawn against itself.
```

```text
RESOLVED THROUGH THE REGISTRY, whatever the manifest recorded: the
tree has moved twice (09-10 sweeps/, 09-16 sensitivity/) and the
run NAME is the stable handle.
```

```text
The index key, (run_name, kind, tag): a cell is a sensitivity row and
its baseline a matrix row, so nothing below looks a run up by name alone.
```

```text
`measured` (None) is not a point on the numeric axis; it sorts last and is
drawn as its own category rather than being given a fake number.
```

```text
RESOLVED, NOT JOINED: the rate CSV is tables/shoreline_change_rate.csv in
the new run layout and {run}_shoreline_change_rate.csv in the old one.
```

```text
The full-record rate PROJECTED change is built from (2026-09-19 advisor
target): the 1996-2024 LRR carried onto the run window.
```

```text
RESOLVED, NOT JOINED BY HAND, AND IN THE ARM THE ROW CLAIMS. This
join had no slot for the arm component, so a run filed under one
resolved to a path that does not exist -- and the `continue` below
then skipped the check in silence, which is the opposite of what a
check is for. find_run_dir raises instead, naming the arms the run
IS under. Passing row.arm rather than letting it default is what
makes this figure read the run it is scoring: the index carries the
arm for exactly this, and as of the 2026-09-01 backfill every row's
arm is the directory it is actually in.
```

```text
A run predating the metadata field records no target and cannot be
checked. That is the only thing this skip is allowed to mean now.
```

```text
The two windows the runner added 2026-09-11 (COASTSAT_DATASETS in
HAT_hindcast_1984_2024.py). Without them a 1996 or 2010 sweep has no
active dataset and the plotter raises.
```

```text
One place, applied once at import. Set here rather than per-figure so every
panel in the directory carries the same type sizes and rules -- a set of
figures read side by side is the unit that has to look consistent, not any
one of them.
```

```text
Two spines, not four. A box around the data adds a rule the reader has to
look past on every panel and encodes nothing.
```

```text
The calibration value is a POINT ON THE CURVE, not just a reference
line. Without it the series jumps straight across its own reference
-- the Hs cells run 2.0 then 3.0 and the segment between them hides
the 2.5 the whole sweep is measured against.
```

```text
An axis the model barely responds to autoscales to its own noise
and reads as a cliff. Floor the span on the run's own RMSE so a
flat response LOOKS flat, at whatever scale that run lives at.
```

```text
Without this the calibration diamond is half-clipped by the
spine on any axis whose calibration value is an endpoint of the
swept range -- which is three of the five.
```

```text
`measured` has no numeric x. It is reported in the caption strip
rather than dropped, and rather than floated inside the axes where it
can land on the curve.
```

```text
The observed layer is ONE curve: the scoring target, which is the
7-domain LOWESS over D11-90 and the raw domain means over D1-10. The
full transect scatter and the unsmoothed domain-mean zigzag were drawn
too (Hannah, 2026-09-29) and buried the model curves under a second,
noisier blue layer; the transect dots stay only over D1-10, where the
target is not a LOWESS and the reader needs to see what it is made of.
```

```text
The line runs through D1-{skip} too, as the raw domain means the target
uses there: without it the southern dots were hard to read as a curve
(Hannah, 2026-09-29, after trying it dots-only).
```

```text
NOT the cell's sibling: since 2026-09-10 a sweep cell lives in
<preset>/sweeps/<family>/ and its baseline is a scenario run one level
up, so sibling arithmetic pointed at sweeps/<family>/<base>. The
registry row carries the period and preset, so ask for the directory.
```

```text
Named with its VALUE, not just "calibration run". On the Hs panel this is
the line the reader is looking for -- where the current setting sits among
the alternatives -- and "calibration run" does not answer that.
```

```text
Solid, not dashed: with a white halo the dash gaps let the ramp curves
show through as orange/red flecks, and the line read as broken (Hannah,
2026-09-29). Weight and the thin halo carry its identity instead.
```

```text
Reference curves only -- the swept values are on the colourbar. Below the
axes, so it cannot cover data at any y-limit.
```

```text
Discrete colourbar: one swatch per cell, because the swept values are not
evenly spaced and the colours are assigned by rank. A continuous bar would
imply a linear mapping from value to colour that is not what was drawn.
```

```text
The blocked-relocation panel is usually all zeros -- real information,
but it does not need half the figure to say it.
```

```text
Further out than the default: this figure's y-labels are long and
the short lower panel brings them closer to the letter.
```

```text
Two entities, fixed colours: the preset owns the colour, so adding a third
preset later cannot repaint these two.
```

```text
PANEL (b) EXISTS BECAUSE OF THE MAGNITUDE GAP. calibBE fits 48 domains to
the target and sits around 0.5-1.2; edgeBE fits two and sits around
1.2-3.0. On one shared axis the calibBE basin is squashed into the bottom
third and the reader cannot compare the two SHAPES -- which is the entire
question. (b) re-plots each curve as a rise above its own minimum, so the
location and width of each basin can be read directly.
```

```text
The calibration Hs is a POINT ON THIS CURVE. Without it calibBE draws
a straight 2.0-to-3.0 segment straight over the top of its own
minimum at 2.5 -- hiding a feature the panel exists to show.
```

```text
Annotated at the TOP of panel (a), which is empty at the calibration Hs --
both curves are near their minima there, well below. At the foot it
grazed the calibBE diamond; moving the legend out from under the title
freed this corner up.
```

```text
Below the axes, where it cannot cover a curve at any y-limit. The diamond
is explained once in the caption instead of doubling every series entry.
```

```text
One directory per (period, preset). Flat output put 22 files with long
names in one listing, where the only way to find the pair you wanted to
compare was to read every name.
```

<details><summary>Function notes (the original docstrings)</summary>

**`period_component()`**

```text
The "<start>_<end>" path component of a period, from the period table.

Not start + 20: the 1996 window ends at 2010, and a sweep run for it files
under 1996_2010 exactly as its scenario run does.
```

**`load_cells()`**

```text
Every completed cell for one period and preset, newest wins.

The manifest is append-only, so re-running a cell adds a second line rather
than replacing the first. Keeping the LAST line per (sweep, value) is what
makes a re-run authoritative -- and it has to be, because a cell re-run
after the resize_interior_domain fix supersedes the one that failed.

Args:
    start_year: Period start year.
    preset: Source/sink preset the cells ran under.

Returns:
    DataFrame with sweep, setting, value, sort_key, run_dir, run_name.
```

**`baseline_name()`**

```text
The matrix run a cell is a departure from.

Every cell is its baseline plus one trailing token, so the baseline is
the name with that token removed. Derived rather than rebuilt from switches
because the token is the only difference by construction -- see the run-name
tokens in cascade_pipeline.hindcast.

Args:
    run_name: A sensitivity cell's directory name.

Returns:
    The baseline run name.
```

**`load_index()`**

```text
run_index.csv, indexed by (run_name, kind, tag).

KEYED ON KIND AND TAG TOO, NOT DEDUPLICATED. A run name describes the SCENARIO
and an arm describes the FORCING, so one name legitimately appears once per
arm -- and the index is keyed on (run_name, Hs_m, arm) for that reason.
This function previously collapsed those with
`drop_duplicates(keep="last")`, which does not pick the calibration run: it
picks whichever row sorts last. On 2026-09-01 that made
`HAT_1984_2004_edgeBE_road_bdm_groin` resolve to the `waveHs3_probe` arm --
a Newton probe at Hs = 3.0 -- and that row is the BASELINE every edgeBE
wave cell is drawn against, so each panel measured its sweep against a
different experiment. Silently: the row exists and reads cleanly.

Since 2026-09-16 the index is keyed on (run_name, kind, tag): a cell is
a `sensitivity` row tagged with its axis and its baseline a `matrix` row,
so the two are distinct rows even though their names differ by one token.
load_run_index translates a file from before that date.

Returns:
    DataFrame indexed by (run_name, kind, tag), one row per key. Numeric
    columns are parsed; the registry reads the file as text.

Raises:
    ValueError: If a key is duplicated, which the rebuild makes impossible
        and so means a fault in the index rather than something to pick
        a winner from.
```

**`model_rates()`**

```text
Per-domain modelled LRR for one run.

Args:
    run_dir: The run's directory.
    run_name: The run's name, which names the run folder and prefixes the
        files still at its root.

Returns:
    DataFrame with gis_domain and lrr_m_yr.
```

**`model_position_change()`**

```text
Per-domain modelled shoreline position change, end minus start (m).

From the saved shoreline matrix, the way rerender_run_figures.py
--position-change draws each run's own figure, so the two agree.
```

**`position_layers()`**

```text
The observed layer in metres: (cs_series, target, legend label).

`reference` "total" is the window's own LRR x span; "projected" is the
1996-2024 LRR x span (the 09-21 vocabulary: named by the FIT window).
Scaling after the LOWESS is exact -- see scale_coastsat_series.
```

**`check_target_matches()`**

```text
Assert every run was scored against the target this figure draws.

The skill numbers plotted here come from run_index.csv, and the curve comes
from build_target_table with TARGET_WINDOW. If a run was scored against a
different window, the panel would show a curve that is not the one its RMSE
refers to. Each run records its target in run_metadata.json, so this is
checkable rather than assumed.
```

**`coastsat_layers()`**

```text
The CoastSat series and the spliced target, built the hindcast's way.

Args:
    start_year: Period start year, which selects the active dataset.

Returns:
    (cs_series, target_frame).
```

**`panel_label()`**

```text
Puts a lower-case panel letter outside the axes, top left.

Outside, so it can never land on a data point. Journals want the letter on
every panel of a multi-panel figure and want it in a fixed place -- reading
a ten-panel grid means finding (f) without hunting for it.
```

**`plot_skill_overview()`**

```text
Interior RMSE and bias against the swept value, one column per axis.

Two ROWS rather than two y-scales on one axes: RMSE and bias are different
measures on different scales, and overlaying them on twin axes is the one
chart form that reliably misleads, because the crossing point of the two
curves is an artifact of the scaling choice.
```

**`plot_alongshore()`**

```text
Modelled LRR per domain for every cell of one axis, over the observed layer.

`position` None draws rates. Otherwise it is (reference, observed legend
label): the model curves are end-minus-start position change (m) and
cs_series / target must already be in metres (position_layers).

THE LEGEND IS NOT IN THE DATA AREA. A 16-entry legend box placed inside
these axes sat directly on the target curve and the southern transect
scatter, which is the densest and most important corner of the panel. The
swept values become a discrete colourbar on the right -- correct anyway,
because they are an ordered ramp rather than unrelated categories -- and
the four reference curves get a horizontal legend beneath the axes.
```

**`plot_relocation_outcomes()`**

```text
Road outcomes against the relocation target.

A bar chart, not a line: the question is a count at each candidate value,
and the interesting feature is a THRESHOLD between two adjacent values, not
a slope. `measured` sits beside the numeric bars as its own category because
it is a different policy, not a larger number.
```

**`plot_circularity()`**

```text
The Hs skill curve under both presets, which is what makes the circularity visible.

calibBE fits a per-domain background-erosion rate against the CoastSat LRR
at the calibration Hs, so its RMSE minimum says as much about the fit as
about the wave climate. edgeBE imposes a rate on the two end domains only,
so its interior is generated by the physics. Two series, so both are named
in a legend and neither is identified by colour alone.
```

</details>

### plot_sensitivity_pair.py

Two settings of one sweep axis, stacked, each against the observed target.

From the script's original header:

```text
Two settings of one sweep axis, stacked, each against the observed target.

The alongshore panels in plot_sensitivity.py overlay every cell; when two
neighbouring values differ only in a few domains, the lines sit on top of each
other. This stacks two runs (top, bottom) on shared axes so the difference
reads by eye. A value equal to the calibration default is the matrix baseline.

Usage:
    python plot_sensitivity_pair.py --start-year 1996         --sweep wave_angle_high_fraction --top 0.5 --bottom 0.51
```

### Deleted 2026-10-01

- `plot_sensitivity_vs_coastsat.py` — replaced by `plot_sensitivity.py`. It
  read a session folder written by the deleted
  `HAT_waveSensitivity_1984_2004.py` and CoastSat CSVs under
  `scripts/input_prep/CoastSat/`, neither of which exists.
- `superseded_20260902/` — two markdown guides from the earlier sweep
  (domain-by-domain analysis; verification and interpretation) and their WHY.md.

To recover, by the path it was last committed under:

```
git log --diff-filter=D --oneline -- scripts/sensitivity_analysis/<path>
git show <commit>^:scripts/sensitivity_analysis/<path>
```

### plot_sensitivity_zoom.py

A zoomed view of one sweep axis over a value range, beside the full set.

From the script's original header:

```text
A zoomed view of one sweep axis over a value range, beside the full set.

plot_sensitivity.py draws every cell of an axis on one panel. A fine sweep
added later (e.g. high-angle fraction 0.51-0.54, 2026-09-29) would crowd that
panel, so this draws only the cells inside [--lo, --hi] with the plotter's own
functions and writes them to a sub-folder. The standard 01-06 set is untouched.

Usage:
    python plot_sensitivity_zoom.py --start-year 1996         --sweep wave_angle_high_fraction --lo 0.5 --hi 0.55
```
