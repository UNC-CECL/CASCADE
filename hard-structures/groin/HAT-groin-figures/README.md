# HAT-groin-figures - the groin's observed effect and how the module represents it

Three figures, each a standalone script, each writing its PNG and PDF beside
itself and its caption into `CAPTIONS.md` through `hat_figure_style.caption()`.
`GROIN_PLAN.md` (one folder up) cites all three. Nothing here runs the model.

| script | figure | reads |
|---|---|---|
| `HAT_groin_two_shorelines_figure.py` | `groin_two_shorelines.png` - start here: the two shorelines either side of the groin and the gap between them | the wet/dry change table |
| `HAT_groin_timeline_figure.py` | `groin_timeline_and_hindcast.png` - the fillet's history against the module's trapping schedule | the wet/dry change table |
| `HAT_groin_module_logic_figure.py` | `groin_module_logic.png` - the module's arithmetic, and why it cannot close the gap | the wet/dry change table, and the 1984-2024 sweep cells in `output/calibration/groin/fullperiod_1984_2024/` |

The wet/dry change table is
`../HAT-groin-buxton-output/shoreline_position_output/Change_from_wetdry_1967_D2_D12.csv`
(24 dated surveys, GIS 2-12, change since 1967, landward-positive).

    python HAT_groin_two_shorelines_figure.py
    python HAT_groin_timeline_figure.py
    python HAT_groin_module_logic_figure.py

Brought in line with `scripts/STYLE.md` on 2026-10-01: header and author
block, one-line comments, CONFIG block, functions in the order `main()` calls
them. Proven unchanged in behaviour with `style_equivalence_check.py`. The
explanations that left the code are below, word for word.

## The scripts in detail

### HAT_groin_two_shorelines_figure.py

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

### HAT_groin_timeline_figure.py

The fillet (GIS 5 minus GIS 6) through time, in two opposite regimes, against
the trapping rate `GroinCallback` applies on its absolute calendar (install
1969, ramp from the 1996 repair to the 2003 storm, then M x f).

From the script's original header:

```text
The Buxton groin's observed history, and how it maps onto the hindcast.

One figure answering three questions that kept getting tangled:

    WHEN did the groin affect the shoreline?
        Continuously from 1969, but in two opposite regimes. It TRAPPED through
        2004 and has been RELEASING since. The turning point coincides with the
        2003 storm damage, not with anything in the model.

    IS THERE A DIFFERENT SIGNAL IN EACH HINDCAST PERIOD?
        Yes, and they have opposite sign. Period 1 (1984-2004) is still
        accumulating; period 2 (2004-2024) is draining. Both rates are computed
        here and quoted in the caption.

    HOW DOES THE MODULE REPRESENT THAT?
        `GroinCallback` carries an absolute calendar schedule -- install 1969,
        deterioration onset 1996 (last repair), linear ramp to 2003 (storm),
        then hold at M*f. The lower panel draws that schedule against the
        observations, which is the clearest way to see that the module's
        built-in timeline already matches the measured history, and where it
        cannot follow.

WHAT THE MODULE CANNOT DO, DRAWN RATHER THAN FOOTNOTED
    Trapping is bounded at >= 0, so the groin can stop adding sand but cannot
    actively drain the fillet. Period 2's observed release is therefore outside
    what the parameterisation can produce at any (M, f), and the lower panel
    hatches that region so the limitation is visible next to the data.

DATA
    Fillet is x_s[D5] - x_s[D6] against a fixed 1967 datum, from
    `Change_from_wetdry_1967_D2_D12.csv` -- 24 dated wet/dry surveys produced by
    the GIS analysis in HAT-groin-gis-analysis. Landward-positive, so a rising
    curve means the updrift side is holding while the downdrift side retreats,
    which is what a groin builds.

STYLE, 2026-09-11
    Drawn under the project house style (`scripts/site_layer/hat_figure_style.py`). Three
    things changed beyond type and colour. The canvas is a 190 mm column
    instead of 14 in, so its type survives a page. The phase washes are gone:
    four overlapping shades on one panel (three phases plus two hindcast
    windows) could not be told apart, so the phases are named with their rates
    in the caption and the hindcast windows are a strip against the top edge.
    And the title sentences and the footnote paragraph are in CAPTIONS.md
    beside the image, with every number in them computed here rather than
    written in.

Usage:
    python HAT_groin_timeline_figure.py

Writes groin_timeline_and_hindcast.png (and .pdf) beside this file.

Author: Hannah A. Henry, UNC CECL
```

Notes moved out of the code (observed_fillet(), its docstring, as of the original):

```text
{year: fillet in metres} against the fixed 1967 datum.
```

Notes moved out of the code (retreat_at(), its docstring, as of the original):

```text
(updrift, downdrift) retreat since 1967, positive landward, one survey.

The caption's "175 m against 25 m" pair: read from the table rather than
written into the prose. The column is matched on its survey year, not
built by format -- the real names carry a trailing unit ("_m") and two of
them a month ("1972_july"), so a formatted name misses.
```

Notes moved out of the code (effective_trapping(), its docstring, as of the original):

```text
M_eff for a calendar year -- mirrors GroinCallback._effective_trapping_rate.
```

Notes moved out of the code (phase(), its docstring, as of the original):

```text
(first survey, last survey, change, rate) inside [start, end].
```

Notes moved out of the code (period_strip(), its docstring, as of the original):

```text
The hindcast windows as a band against the top edge, named once.

A strip rather than a full-height wash, for the same reason the phase
washes were dropped: several overlapping greys on one panel read as one.
The twin of this helper is in HAT_groin_two_shorelines_figure.py.
```

Notes moved out of the code (comment above line 79, as of the original):

```text
# The structure's documented history. These are the values GroinCallback is
# configured with, not fitted quantities.
```

Notes moved out of the code (comment above line 85, as of the original):

```text
# The chosen parameters. M is an EFFECTIVE, grid-specific rate -- see the
# accompanying GROIN_PLAN.md for why it is not a sediment flux.
```

Notes moved out of the code (comment above line 93, as of the original):

```text
# The schedule is the thing under test; the observations are the record it is
# read against. The RdBu vintage pair is not used here -- nothing on this
# figure is a 1984/1997 vintage.
```

Notes moved out of the code (comment above line 191, as of the original):

```text
# Clear of the curve: at 0.025 the install label ran through the 1970
# survey, which reads as a label for the marker.
```

Notes moved out of the code (comment above line 214, as of the original):

```text
# The floor the module cannot go below, hatched rather than washed so it
# cannot be mistaken for another shaded period.
```

### HAT_groin_module_logic_figure.py

The module's four lines of arithmetic as a schematic, beside the observed gap
and two sweep cells (groin off, and M 60 f 0.5 standing in for the chosen
f 0.6, which that sweep's grid never ran).

From the script's original header:

```text
How the module's logic produces (and fails to produce) the observed gap.

The two-shoreline figure shows WHAT happened. This one shows WHY the module can
follow part of it and not the rest, by putting the mechanics next to the
consequence.

THE MECHANICS (left panel)
    `GroinCallback` is four lines of arithmetic, applied once a year just
    before BRIE's alongshore solve:

        M_eff        = M, tapering to M*f between 1996 and 2003
        dx_updrift   = -M_eff      seaward advance at D6
        dx_downdrift = +M_eff      landward retreat at D5

    So each year the module pushes the two sides APART by 2*M_eff, and BRIE's
    diffusion immediately starts spreading that dipole back out. The modelled
    gap is the balance of those two.

WHY THAT MATTERS (right panel)
    M_eff is bounded at >= 0. There is no value of M or f that makes the module
    pull the two sides together. It can widen the gap, or -- at f = 0 -- stop
    widening it and let diffusion slowly close it. It can never actively close
    it.

    The observations do close it after 2004. Diffusion alone manages a small
    fraction of that; the fraction is computed here from the groin-off run and
    quoted in the caption. So the post-2004 narrowing is outside the
    parameterisation, and no choice of (M, f) reaches it. That is the single
    fact that determines how the module should be applied.

DATA
    Observed gap from the wet/dry surveys; modelled gap from the continuous
    1984-2024 sweep cells, which carry a full annual trajectory each.

STYLE, 2026-09-11
    Drawn under the project house style (`scripts/site_layer/hat_figure_style.py`): a
    190 mm column rather than a 15 in canvas, the sheltered side in ACCENT
    purple and the unprotected side in BASE grey to match the two-shoreline
    figure, and the schematic's two shouted sentences and the three footnote
    paragraphs moved into CAPTIONS.md beside the image.

    The "about a tenth" that footnote asserted is now measured off the
    groin-off run each time the figure is drawn, and it does NOT survive the
    measurement as it was stated. Diffusion alone closes the gap by about a
    sixth of the observed NET change across the window, which is where that
    number came from; its post-2004 RATE is about a sixty-sixth of the observed
    closure rate, which is what the sentence claimed. The caption carries both,
    because they disagree by a factor of ten and the rate is the one the
    argument rests on.

Usage:
    python HAT_groin_module_logic_figure.py

Writes groin_module_logic.png (and .pdf) beside this file.

Author: Hannah A. Henry, UNC CECL
```

Notes moved out of the code (observed_gap(), its docstring, as of the original):

```text
{year: gap in metres} = D5 retreat minus D6 retreat, from 1967 datum.
```

Notes moved out of the code (modelled_gap(), its docstring, as of the original):

```text
Modelled gap per year for one sweep cell, referenced to its own year 0.
```

Notes moved out of the code (rate(), its docstring, as of the original):

```text
Endpoint rate of `series` over the samples inside [start, end].
```

Notes moved out of the code (draw_mechanics(), its docstring, as of the original):

```text
Schematic: what the module adds to the two domains each year.

Laid out for a HALF of a 190 mm column, which is about a quarter of the
canvas this panel had before 2026-09-11. Every label was rewritten short
and moved OUTWARD from the structure: at the printed width the old
full-width text ran through the groin and through its own arrows, and a
label touching an arrow reads as a label FOR that arrow.
```

Notes moved out of the code (comment above line 148, as of the original):

```text
# the two shorelines, each labelled at its OUTER end and on the side away
# from the arrows
```

Notes moved out of the code (comment above line 160, as of the original):

```text
# What the module applies each year.
# Each label sits beyond its own arrowhead and on the SAME SIDE of the
# groin as its arrow, on one line: beside the shaft it landed at the same
# height as the shoreline label next to it, and on the far side it ran
# straight through the structure.
```

Notes moved out of the code (comment above line 206, as of the original):

```text
# The chosen pair is (60, 0.6), but the continuous-window sweep's
# f grid runs 0.1/0.3/0.5/0.7/0.9, so f=0.6 was never run here.
# M60_f0.50 is its nearest neighbour and the label says so.
```

Notes moved out of the code (comment above line 257, as of the original):

```text
# BOTH comparisons, because they disagree by a factor of ten and the
# footnote this caption replaces quoted only the flattering one. It
# said diffusion manages "about a tenth" of the observed closure; that
# holds for the net change across the window ({net_share}), not for the
# post-2004 rate ({rate_share}), which is what the sentence claimed.
```
