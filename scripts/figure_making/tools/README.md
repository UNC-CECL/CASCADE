# figure_making/tools — helpers for the figures, not figures themselves

Scripts that keep the published figures reproducible and findable, and two
small viewers. `regenerate_all_figures.py` redraws everything under
`output/figures/` and then runs `figure_index.py`, which writes
`output/figures/README.md`.

```
regenerate_all_figures.py  redraw every published figure from its producer, in layout order
figure_index.py            write output/figures/README.md: layout, captions, producers, dates
clip_nc_coast.py           rebuild the clipped NC coastline map layer from the D: drive
color_picker.py            preview a hex-to-hex colour gradient (needs rgb_gradient)
npy_view.py                show one saved .npy elevation array
```

What not to trust: `npy_view.py` points at
`data/hatteras_init/topography/2009_FIXED/`, which no longer exists; edit the
path before use. `clip_nc_coast.py` needs the D: GIS drive.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### clip_nc_coast.py

Rebuild the clipped NC 1:80k coastline the figures read, from the source on the D: GIS drive.

From the script's original header:

```text
Rebuild data/hatteras_init/map_elements/nc_coast_80k/ from the NC 1:80k
coastline on the D: GIS drive.

    python scripts/figure_making/tools/clip_nc_coast.py

WHY
    The overwash maps drew the coast straight off
    D:/Hatteras_GIS/Outlines/nc_80k/nc_80k.shp -- 6 MB for the whole state --
    so they could not be made without the drive. They only ever use the land
    in a window a few km around the domain boxes, so the repository keeps that
    window: every polygon within NC_COAST_PAD_M (25 km) of the boxes, clipped,
    as a geojson (0.75 MB) in EPSG:4326 like its source.

    Re-run it only if the source changes; the output is what the figures read.
```

### color_picker.py

Preview a linear colour gradient between two hex colours, with each step's hex code.

### figure_index.py

Write output/figures/README.md: the figure layout, every figure with what it shows and what drew it.

From the script's original header:

```text
Writes `output/figures/README.md`, the map of the figures: the layout, a
"where do I find..." table, and one table per folder listing every figure,
what it shows, the script that draws it and the day it was drawn.

    python scripts/figure_making/tools/figure_index.py

Re-run after adding or redrawing a figure; `regenerate_all_figures.py` runs it
at the end.

WHY IT IS GENERATED
    A hand-kept index of a folder ~25 scripts write into is stale the day after
    it is written, and a wrong index is worse than none. Everything in the
    table is recorded beside the figures:
      * "shows" is the first sentence of the figure's entry in its folder's
        `supporting/CAPTIONS.md`;
      * "drawn by" comes from `supporting/producers.json`, which
        `regenerate_all_figures.py` writes by noting which files each producer
        touched; for a figure redrawn by hand since, it falls back to searching
        the scripts tree for the figure's file name;
      * "drawn" is the PNG's modification date, so a stale figure shows.
    A dash in "shows" or "drawn by" is a real finding: nobody wrote down what
    the figure shows, or nothing can reproduce it.

THE LAYOUT (Hannah, 2026-09-29): numbered in the paper's order, see LAYOUT.
    The folder names come from hat_figure_style.FIGURE_SUBJECTS / INPUT_STEPS.
```

Notes that were in the code:

```text
Every folder that holds figures, in reading order, with one line saying what
it answers. A folder on disk that is missing here is listed under "other" so
it cannot hide.
```

```text
"Where do I find ...": the questions people actually ask, pointed at a path
under output/figures/. Each path is checked; a missing one is flagged.
```

<details><summary>Function notes (the original docstrings)</summary>

**`source_captions()`**

```text
{file name: caption} from every CAPTIONS.md under data/ and
output/comparisons/. A figure PUBLISHED as a copy (the mean-shoreline
images, copied from data/hatteras_init/5-scr/) has its caption beside the
original, not beside the copy; this finds it by file name.
```

**`searched_producers()`**

```text
{figure stem: script}, by searching the scripts tree for the file name.
The fallback for a figure redrawn by hand after the last full regeneration.
```

</details>

### npy_view.py

Show one saved .npy elevation array as an image.

### regenerate_all_figures.py

Redraw every figure published under output/figures/ from its own script, then rewrite the index.

From the script's original header:

```text
Redraws every figure published under output/figures/ from its own script, in
the numbered layout's order, then rewrites the index (output/figures/README.md).

    python scripts/figure_making/tools/regenerate_all_figures.py            # everything
    python scripts/figure_making/tools/regenerate_all_figures.py --only 4-model-mechanics
    python scripts/figure_making/tools/regenerate_all_figures.py --list

WHY IT EXISTS
    output/figures/ is gitignored and fed by ~25 scripts, so "are the figures up
    to date?" had no answer short of remembering which script drew what. On
    2026-09-29, 44 of 182 figures were still from 09-17, before option-A waves
    and the 09-28 adoption. This is the answer now: one list of producers, run
    in order, each one's log kept, failures reported at the end rather than
    stopping the rest. Add a producer here when a script starts publishing to
    output/figures/.

    The logs go to output/logs/scratch/figures_<timestamp>/, one per step.
```
