# figure_making/shoreline — shoreline change, observed and modelled

The observed shoreline record and its rates (CoastSat, the target; DSAS, the
independent check), the chainage animations of the raw CoastSat record, and
a poster tool.

```
plot_coastsat_calibration_periods.py   THE observed-rate figure: CoastSat per run period
plot_coastsat_poster.py                the same, LOWESS-smoothed (7 domains)
dsas/plot_dsas_calibration_periods.py  the DSAS check, same style, under supporting/
dsas/dsas_from_gis.py                  format the GIS-computed DSAS domain statistics
dsas/dsas_rate_verification.py         check the DSAS rates against known values
dsas/dsas_shoreline_analysis.py        DSAS rates, all in one: table and four figures
chainage/shoreline_chainage_*.py       40 years of CoastSat chainage: decade panels, GIFs
rodanthe_erosion_example_poster.py     Rodanthe erosion trends, for a poster
```

What not to trust:
- `dsas/dsas_shoreline_analysis.py` writes its CSV and four PNGs into the
  folder it is run from, not under `output/` (ORGANIZATION.md rule 1).
- The three older DSAS scripts run straight through with no `main()`.
- `rodanthe_erosion_example_poster.py` had its header stranded below an
  inserted import (from the 2026-09-22 rename); it is restored as the header,
  and its original text is below.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### chainage/shoreline_chainage_alldata_evolution.py

CoastSat shoreline chainage over 40 years, every observation: decade panels and an animated GIF.

From the script's original header:

```text
Visualises CoastSat shoreline chainage across Hatteras Island over 40 years.
Uses ALL individual CoastSat observation dates for the GIF (not seasonal medians),
giving the highest temporal resolution view of shoreline change.

Outputs
1.  Four 10-year panel figures (one PNG each):
        1984–1994, 1994–2004, 2004–2014, 2014–2024
    Shared 1984 baseline and shared y-axis across all panels for direct comparison.
    Lines coloured light→dark as years pass within each decade.

2.  A decadal GIF cycling slowly through the four panel PNGs.

3.  A high-temporal-resolution GIF with one frame per CoastSat observation date
    (~1200+ frames), showing individual acquisition snapshots across the full record.

Usage
    Edit the CONFIG section, then run:
        python shoreline_chainage_alldata_evolution.py

Dependencies
    pip install pandas numpy matplotlib scipy imageio pillow tqdm
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). This file drew in
matplotlib's defaults until 2026-09-17 -- it never called apply_style().
```

```text
Typeface only: this script writes ANIMATION frames, and the printed-width
rule does not apply to something that is never printed. Its figsize is
the frame size and is left as it is.
```

```text
The per-transect CoastSat timeseries live in the data tree, not under
scripts/; the old value was a driveless path that never existed (2026-09-10).
```

```text
"input_preperation" is the pre-2026 folder name; the lookup now lives
under data/hatteras_init/5-scr/2-transect-frame/transect_domains/ (moved out of the
scripts tree 2026-09-12; hat_observed_rates.py resolves it).
```

```text
Products go to output/, not beside the script (2026-09-13). These were
absolute paths into a home directory, so they resolved on one machine
and dumped 1212 files into the code tree. Rule 1 and rule 5 of
ORGANIZATION.md.
```

<details><summary>Function notes (the original docstrings)</summary>

**`make_gif_cmap()`**

```text
Custom diverging colormap: dark crimson→red→grey→blue→deep navy.
Stays saturated near zero so small deviations are visible on grey background.
Negative = erosion = red, Positive = accretion = blue.
```

**`save_gif_pil()`**

```text
Assemble frames into a GIF using PIL with correct duration handling.
PIL has a 2× duration bug; passing duration_s * 500 corrects it.
```

</details>

### chainage/shoreline_chainage_annual_evolution.py

CoastSat shoreline chainage over 40 years, annual medians: decade panels and a seasonal GIF.

From the script's original header:

```text
Visualises CoastSat shoreline chainage across Hatteras Island over 40 years.

Outputs
1.  Four 10-year panel figures (one PNG each):
        1984–1994, 1994–2004, 2004–2014, 2014–2024
    Each shows annual median chainage deviation from the period-start baseline
    for every individual transect.  Lines coloured light→dark as years pass.

2.  A GIF animating spring (Apr–May) and fall (Oct–Nov) seasonal median
    shoreline profiles year by year, 1984–2024.

Usage
    Edit the CONFIG section, then run:
        python shoreline_chainage_alldata_evolution.py

Dependencies
    pip install pandas numpy matplotlib scipy imageio tqdm
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). This file drew in
matplotlib's defaults until 2026-09-17 -- it never called apply_style().
```

```text
Typeface only: this script writes ANIMATION frames, and the printed-width
rule does not apply to something that is never printed. Its figsize is
the frame size and is left as it is.
```

```text
The per-transect CoastSat timeseries live in the data tree, not under
scripts/; the old value was a driveless path that never existed (2026-09-10).
```

```text
"input_preperation" is the pre-2026 folder name; the lookup now lives
under data/hatteras_init/5-scr/2-transect-frame/transect_domains/ (moved out of the
scripts tree 2026-09-12; hat_observed_rates.py resolves it).
```

```text
Products go to output/, not beside the script (2026-09-13). These were
absolute paths into a home directory, so they resolved on one machine
and dumped 1212 files into the code tree. Rule 1 and rule 5 of
ORGANIZATION.md.
```

```text
── Bottom x-axis: replace generic label with domain tick marks ───────────
Place domain number ticks on the primary axis
```

```text
── Compute 1984 baseline ONCE — reused by all panels and GIF ────────────
Use full-year 1984 median (not spring-only) for panels so early Landsat
years with sparse seasonal coverage still get a valid baseline.
```

```text
── Pre-compute all annual deviations from 1984 baseline across full record
so we can set a shared y-axis range before drawing any panel.
```

```text
Use PIL directly — imageio duration is unreliable across versions/viewers.
PIL duration is always milliseconds, consistent across all GIF viewers.
```

```text
Custom diverging colormap: saturated red (erosion/negative) → mid-grey
→ saturated blue (accretion/positive). Stays visible near zero unlike
RdBu which passes through white. Negative = erosion = red, Positive = accretion = blue.
```

<details><summary>Function notes (the original docstrings)</summary>

**`annotate_axes_publication()`**

```text
Apply the full publication annotation suite to ax.
ylim must be the (ymin, ymax) already set on ax so labels are placed correctly.
```

</details>

### dsas/dsas_from_gis.py

Format the DSAS domain statistics already calculated in GIS for the Python scripts.

From the script's original header:

```text
Format DSAS shoreline change data for Python scripts
Uses domain-level statistics already calculated in GIS

This script formats pre-calculated domain means from GIS for use in Python analysis.
No transect-level filtering or aggregation needed - just clean formatting.
```

Notes that were in the code:

```text
Anchored 2026-09-22. These were drive-rooted into a home directory, naming
data/hatteras_init/shoreline_change/ -- a folder that has not existed since
the DSAS tables moved into the 5-scr tree. Both FILENAMES survived the move
unchanged, so this is a resolved relocation, not a guess: the resolver owns
the folder (rule 6) and the repo is found by searching upward (rule 5).
```

### dsas/dsas_rate_verification.py

Check the calculated DSAS shoreline change rates against known values.

From the script's original header:

```text
Verification: Compare Calculated Rates to Known Values
Checks if your calculated shoreline change rates make sense
```

Notes that were in the code:

```text
The DSAS tables moved into the data tree 2026-09-13 (rule 1). These
were read by bare filename, so they only resolved when you happened
to run from this folder.
```

```text
HOUSE STYLE. Like the poster script, this imported record_caption at the
BOTTOM and never called apply_style(), so it was not in the project
typeface at all (2026-09-17).
```

### dsas/dsas_shoreline_analysis.py

DSAS shoreline change rates, all in one: domain table and four figures.

From the script's original header:

```text
Complete Shoreline Change Rate Analysis - All-in-One
Analyzes DSAS data and creates publication-quality visualizations
Run this once and get everything!
```

Notes that were in the code:

```text
The DSAS tables moved into the data tree 2026-09-13 (rule 1). These
were read by bare filename, so they only resolved when you happened
to run from this folder.
```

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5). This file drew in
matplotlib's defaults until 2026-09-17 -- it never called apply_style().
```

### dsas/plot_dsas_calibration_periods.py

Observed DSAS shoreline change rate per domain for the run periods: the independent check on CoastSat.

From the script's original header:

```text
Observed shoreline change rate per domain from the DSAS transect record --
the INDEPENDENT CHECK on the CoastSat figure, not the target.

CoastSat is what the model is graded against (see
coastsat_calibration_periods.png, and hat_observed_rates), so this one lives
under supporting/ and is drawn in exactly the same style, so the two can be
laid side by side and only the data differs (Hannah, 2026-09-17).

TWO THINGS THIS FIGURE CANNOT MATCH EXACTLY, and both are stated in the
caption rather than smoothed over:

  THE WINDOWS. DSAS has five digitised shorelines -- 1978, 1987, 1997, 2009,
  2019 -- so it cannot be cut to the run periods. The nearest pairs are used:
  1997-2009 stands in for 1996-2010 (one year in at each end) and 2009-2019
  for 2010-2024 (one year late, five years short). Until 2026-09-17 it drew
  1978-1997 and 1997-2019, which straddle the run periods rather than
  approximating them.

  THE ESTIMATOR. This is an END-POINT RATE: the first and last shoreline of
  the window, over the elapsed years. CoastSat is an OLS slope through every
  transect observation in the window, which is what the model is scored with
  (see hat_observed_rates and the LRR note in HAT_hindcast_methods). With two
  shorelines an OLS slope IS the end-point rate, so the two agree in form
  here; they would not if a third DSAS vintage fell inside a window.
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5), so this block is
independent of whatever this script calls its own repository variable.
```

```text
The run periods this is checking, and the DSAS pair that stands in for each.
Both are stated so a reader never has to infer the offset.
```

```text
THE LIMITS GO FIRST. town_bands() skips any span outside the current
view and clamps a label to the visible part of its span, so calling it
before set_xlim silently dropped Buxton (GIS 7-8): the axes had
autoscaled to the shoal spans and 6.5-8.5 fell outside them. Band and
label both vanished, with no error (2026-09-17).
```

```text
a second row, clear of the village names at 0.985: Avon Shoals
spans Avon and Wimble Shoals spans Tri-Village, so the two sets
of labels overlap in x and must differ in y
```

```text
NAMED, like the shoals and the structures: if a band is worth
drawing it is worth naming (Hannah, 2026-09-17). town_bands puts
these at the top of the panel, and structures() already knows to
tuck its own labels under them.
```

```text
A domain is NaN where no transect in it carries BOTH of the window's
shorelines; matplotlib breaks the line and the fill there rather than
bridging a gap that has no data behind it.
```

```text
A CHECK, SO IT LIVES UNDER supporting/. The caption is keyed to the
figure's name in the folder's one CAPTIONS.md, beside the primary's.
```

### plot_coastsat_calibration_periods.py

Observed CoastSat shoreline change rate per domain, one curve per run period: the primary shoreline figure.

From the script's original header:

```text
Observed CoastSat shoreline-change rate per domain, one curve per run period.

This is the PRIMARY shoreline figure. CoastSat is the target the model is
graded against, so the DSAS version of the same plot is the independent check
and lives under supporting/ (Hannah, 2026-09-17).

WHAT CHANGED 2026-09-17, and why
  * THE PERIODS WERE THE OLD PAIR. It read the 1984-2004 and 2004-2024 LRR
    products; the canonical chain is 1996 -> 2010 -> 2024 now. That is NOT a
    relabel: the windows have their own LRR products and the rates genuinely
    differ (GIS 1 is -4.16 m/yr over 1984-2004 and +3.23 over 1996-2010), so
    the curves move with the labels. PERIOD_STARTS drives both, and the ends
    come from HATTERAS_PERIODS.
  * THE ANNOTATIONS WERE ITS OWN. Village spans, groin and pier lines were
    re-declared here as literals -- a fourth copy of what the site config
    already holds -- and drawn in this file's own style. They come from
    town_bands() and structures() now, so they match the reach figures
    exactly and cannot disagree with the config.
  * BOTH SHOAL ZONES. Only Wimble was drawn; HATTERAS_ANNOTATIONS.shoal_zones
    has Avon Shoals (GIS 24-39) as well, and both are real (Hannah,
    2026-09-17: "do both").
  * THE PLACE NAMES CAME OFF THE CANVAS. Buxton / Avon / Tri-Village / Salvo /
    Waves / Rodanthe were printed in the data area; the house style puts that
    naming in the caption, and the bands still show where they are.
  * THE DIRECTION MARKERS WENT INTO THE AXIS LABEL. "Accretion" and "Erosion"
    with arrows sat in the panel; the axis says "+ seaward" now, which is the
    same statement in the place the house style keeps it.
  * LEGEND OUT OF THE PANEL, frameless, below -- as the management figures.
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5), so this block is
independent of whatever this script calls its own repository variable.
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
The canonical chain. Each end comes from the site config, and the LRR
product for a window lives under <start>_<end>/, so naming the starts names
the data as well as the labels.
```

```text
SHOAL ZONES, both of them, under everything. They are the one warm
accent on this panel; the periods own the red/blue.
THE LIMITS GO FIRST. town_bands() skips any span outside the current
view and clamps a label to the visible part of its span, so calling it
before set_xlim silently dropped Buxton (GIS 7-8): the axes had
autoscaled to the shoal spans and 6.5-8.5 fell outside them. Band and
label both vanished, with no error (2026-09-17).
```

```text
a second row, clear of the village names at 0.985: Avon Shoals
spans Avon and Wimble Shoals spans Tri-Village, so the two sets
of labels overlap in x and must differ in y
```

```text
village spans, from the site config, drawn as every alongshore figure
draws them; unlabelled, because the caption names them
NAMED, like the shoals and the structures: if a band is worth
drawing it is worth naming (Hannah, 2026-09-17). town_bands puts
these at the top of the panel, and structures() already knows to
tuck its own labels under them.
```

```text
structures() measures text against the settled layout, so it goes after
the data and after anything that resizes the axes
```

### plot_coastsat_poster.py

The LOWESS-smoothed companion to the CoastSat calibration-periods figure (7 domains, 3.5 km).

From the script's original header:

```text
The LOWESS-SMOOTHED companion to coastsat_calibration_periods.png: the same
CoastSat rates over the same two run periods, smoothed over a 7-domain (3.5 km)
window so the alongshore pattern reads without the domain-to-domain scatter.

Drawn in exactly the style of the primary, so the two can be laid side by side
and only the smoothing differs.

WHAT CHANGED 2026-09-17
  * THE PERIODS ARE THE CANONICAL CHAIN, 1996 -> 2010 -> 2024, read from
    HATTERAS_PERIODS, and each curve comes from that window's own LRR product.
    It drew 1984-2004 / 2004-2024 before.
  * IT WAS 13 INCHES WIDE, a poster size, saved with a tight bbox so the
    labels rather than figsize() decided the width. At that size its 12 pt
    bold axis labels reduce to about 5 pt on a page. It is 190 mm now, with
    explicit margins and no tight bbox. Render it at a poster width
    deliberately if a poster copy is wanted.
  * THE HOUSE STYLE WAS NEVER APPLIED. The file imported record_caption at the
    BOTTOM and never called apply_style(), so it was not in the project
    typeface at all despite sitting beside figures that are.
  * OFF THE CANVAS: a two-line bold title, the S/N end labels, the place names
    (Buxton / Avon / Tri-Village / Salvo / Waves / Rodanthe), the
    Accretion/Erosion markers and a footnote paragraph. All caption material
    under figure_making/STYLE.md, and the caption carries it now.
  * THE COLOURS ARE THE VINTAGE PAIR. It used a dark blue and a brown of its
    own; the earlier period is C_1984 red and the later C_1997 blue here, as
    in every other figure that draws two periods.
  * THE LEGEND IS OUT OF THE PANEL. It was boxed and sat on the data.
```

Notes that were in the code:

```text
HOUSE STYLE: one typeface and one palette across every figure in this
project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
found by searching upward (ORGANIZATION.md rule 5), so this block is
independent of whatever this script calls its own repository variable.
```

```text
THE LIMITS GO FIRST. town_bands() skips any span outside the current
view and clamps a label to the visible part of its span, so calling it
before set_xlim silently dropped Buxton (GIS 7-8): the axes had
autoscaled to the shoal spans and 6.5-8.5 fell outside them. Band and
label both vanished, with no error (2026-09-17).
```

```text
a second row, clear of the village names at 0.985: Avon Shoals
spans Avon and Wimble Shoals spans Tri-Village, so the two sets
of labels overlap in x and must differ in y
```

```text
NAMED, like the shoals and the structures: if a band is worth
drawing it is worth naming (Hannah, 2026-09-17). town_bands puts
these at the top of the panel, and structures() already knows to
tuck its own labels under them.
```

### Deleted 2026-10-01

- `plot_shoreline_from_npz.py` — yearly shoreline GIFs and the rate figure
  from a saved run; it read `HAT_2004_2024_base_newbufferv3.npz` by absolute
  path (a run deleted long before) and the retired 1984-2004 / 2004-2024
  windows.

To recover, by the path it was last committed under:

```
git log --diff-filter=D --oneline -- scripts/figure_making/shoreline/<path>
git show <commit>^:scripts/figure_making/shoreline/<path>
```

### rodanthe_erosion_example_poster.py

CoastSat shoreline erosion around Rodanthe: annual positions and LRR trends per domain, for a poster.

From the script's original header:

```text
CoastSat Shoreline Erosion Trends — Rodanthe Area
Hatteras Island, NC

Produces a single square publication-quality figure for poster use showing:
  - Annual mean shoreline position per domain (thin lines) — shows the raw signal
  - LRR trend line per domain (bold, full-period) — shows the erosion rate

Two messages conveyed:
  1. Rodanthe has been eroding consistently over the full record
  2. Different parts of Rodanthe erode at different rates

Data source: CoastSat satellite-derived shorelines

Usage
  1. Edit the CONFIG section below.
  2. Run: python plot_coastsat_rodanthe_erosion.py
```

Notes that were in the code:

```text
Anchored 2026-09-14. These named a home directory or a tree renamed
twice over, so none resolved. Rule 5 of ORGANIZATION.md.
```

```text
A poster figure, so it goes to output/figures/ with the others (2026-09-18);
it used to be written in among the rate fits, 5-scr/coastsat_lrr/rodanthe_plots
(those older copies are in 5-scr/archive/rodanthe_plots/).
```
