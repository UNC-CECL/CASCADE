# figures

The fillet's history against the dipole module's trapping schedule (install 1969, ramp, then M x f). Reads the wet/dry change table. Superseded with the dipole.

    python HAT_groin_timeline_figure.py

Writes `groin_timeline_and_hindcast.png` beside the script, and its PDF and caption into `supporting/`. Brought in line with `scripts/STYLE.md` on 2026-10-01; the notes that left the code are below, word for word.

## HAT_groin_timeline_figure.py

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
    the GIS analysis in 1-observations/shoreline_rates_by_era. Landward-positive, so a rising
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
