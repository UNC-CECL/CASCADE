# figures

The dipole module's arithmetic as a schematic, beside the observed gap and two dipole sweep cells. Reads the wet/dry change table and the 1984-2024 dipole sweep cells, which moved to `D:\CASCADE_offload\output\calibration\groin\fullperiod_1984_2024\` on 2026-10-07, so it only redraws with that drive attached.

    python HAT_groin_module_logic_figure.py

Writes `groin_module_logic.png` beside the script, and its PDF and caption into `supporting/`. Brought in line with `scripts/STYLE.md` on 2026-10-01; the notes that left the code are below, word for word.

## HAT_groin_module_logic_figure.py

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
