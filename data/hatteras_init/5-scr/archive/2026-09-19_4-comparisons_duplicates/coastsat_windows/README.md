# coastsat_windows — the observed CoastSat rate, window by window

Every file here is written by
`scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_windows.py`. Nothing is
edited by hand, and a re-run regenerates the whole folder.

```
lrr_four_windows.png        START HERE: the four model windows as a 2 x 2,
                            the 1984-start period left, the 1996-start right
supporting/
    CAPTIONS.md             caption for the 2 x 2
    lrr_four_windows.pdf
    lrr_windows_wide.csv    the four windows' domain mean and std side by side
    y_bounds.txt            the shared y bound (±8 m/yr) and how it was set

1984_2004/  1996_2010/  2004_2024/  2010_2024/
    lrr_<start>_<end>.png   one window, full width, on the shared ±8 axis
    supporting/             its PDF and caption

1996_2024/                  THE LONG WINDOW: context, not a model target
    lrr_1996_2024_halves.png   (a) 1996-2024 filled by sign; (b) 1996-2010
                               and 2010-2024 as lines; fills and shoals marked
    supporting/
        lrr_1996_2024_halves.csv   the three domain means side by side
        lrr_1996_2024_halves.pdf, CAPTIONS.md
```

## Which figures compare directly

The four window figures and the 2 x 2 share one y axis, **±8 m/yr**, so
any panel can be compared with any other. The `1996_2024/` figure does
**not** use it. Its axis is the tightest whole metre that holds every line
(±7), and its caption says so. That figure is also published to
`output/figures/shoreline/`.

## Rebuilding

```
python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_windows.py                     # 2 x 2 + four windows
python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_lrr_windows.py --overlay 1996_2024 # 1996_2024/
```

Set `PYTHONIOENCODING=utf-8` on Windows if the console chokes on the
en-dashes. The rates behind these figures are in
`3-rates/coastsat/lrr/<start>_<end>/domain_lrr_summary.csv`.

## History

Until 2026-09-18 every figure sat flat in this folder with one shared
`supporting/`, so the ±8 window figures and the ±7 long-window figure lay
side by side with nothing to tell them apart. The layout now matches
`coastsat_vs_duneline/`: one folder per window, the grid at the top.
