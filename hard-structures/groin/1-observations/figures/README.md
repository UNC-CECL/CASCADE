# figures

The groin in one picture: the shorelines at GIS 5 and GIS 6 from 24 dated wet/dry surveys, and the gap between them. Reads `../wetdry_photo_positions/Change_from_wetdry_1967_D2_D12.csv`. Nothing here runs the model.

    python HAT_groin_two_shorelines_figure.py

Writes `groin_two_shorelines.png` beside the script, and its PDF and caption into `supporting/`. Brought in line with `scripts/STYLE.md` on 2026-10-01; the notes that left the code are below, word for word.

## HAT_groin_two_shorelines_figure.py

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
