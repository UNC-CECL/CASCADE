# 3-rates — the model targets

The rate fits the model is graded against. These are **model inputs**: the
hindcast runner reads the LRR table on every run, and the calibrated
source/sink preset is fitted against it. Everything else here is a reading of
the same record taken a different way.

Grouped by source, then by what kind of number it is.

```
coastsat/
    lrr/           the OLS rate fit — the thing runs are actually graded on
    endpoint/      net change between +/-6-month means at the dune-line dates
    5yr_bins/      the OLS in successive 5-year bins: WHEN did it change?
    window_convergence/  the OLS on NESTED window families pinned at each
                   end: WHICH WINDOWS recover the long-term rate?
    total_change/  the rate as a distance, and the projections
    extension/     the same fit beyond the 90 surveyed domains
duneline/
    endpoint       end line minus start line
rates_figures.py   one house-style figure per window, for all of the above
```

## A window is an interval, and the ends are not equivalent

Folders are named `<start>_<end>` because a rate fit spans an interval (rule 2
of `ORGANIZATION.md`). The five windows are **not** peers: `1996 -> 2010 ->
2024` is the main chain, `1984_2004` and `2004_2024` are the older one, and
`1996_2024` is context that nothing is graded against.
`data/hatteras_init/5-scr/WINDOWS.md` says which is which, generated from
`hat_observed_rates.WINDOW_ROLE` by `../tools/windows_index.py` so the
table cannot drift from the code.

## Why the figure script is not inside a product folder

`rates_figures.py` draws the figure for *every* product here, beside its
tables. Putting it in any one product folder would misfile it four times over.
Run it after whichever product you rebuilt. It replaced the autoscaled
quick-looks the LRR fit used to draw itself, which are archived under
`5-scr/archive/coastsat_lrr_quicklooks_20260918/`.

Two figures stay at a window folder's top level with `lrr_<w>.png` —
`smoothing_windows_<w>.png`, from
`coastsat/lrr/coastsat_lrr_smoothing_windows.py`. By contrast
`coastsat_lrr_transect_zoom.py` writes into a `transects/` subfolder, because
that family grows a file per reach and per variant.

## The scripts

```
coastsat/lrr/coastsat_domain_lrr.py
    The fit. One window per run; reads the lookup from ../2-transect-frame/
    and the chainage from 1-observations/coastsat_timeseries/. Tables only
    since 2026-09-18 — the figure is drawn by rates_figures.py.
    (Was coastsat_domain_lrr_fixed.py until 2026-09-22; the "_fixed" had no
    unfixed sibling left to distinguish it from.)

coastsat/lrr/coastsat_lrr_windows.py
    Every window side by side on ONE y axis. Each window's own figure is
    autoscaled, so four of them together draw a 2 m/yr swing as tall as a
    7 m/yr one and the eye reads the wrong story. Also the module the other
    figure scripts import for the shoals, fills and structure drawing.

coastsat/lrr/coastsat_lrr_smoothing_windows.py
    The unsmoothed field plus the LOESS at 3, 5 and 10 domains (1.5, 2.5 and
    5.0 km) on one axis, darkest being the window runs are graded at.

coastsat/lrr/coastsat_lrr_transect_zoom.py
    One window at transect resolution over a short reach, transects on the
    x axis instead of domains.

coastsat/endpoint/coastsat_endpoint.py
    The mean CoastSat position within +/-6 months of each dune-line survey
    date, end minus start, in metres and m/yr. Centred on the dune dates so
    that a gap between the shoreline and the dune line is a real difference
    and not an artefact of comparing different moments.

coastsat/5yr_bins/coastsat_5yr_bins.py         the table
coastsat/5yr_bins/coastsat_5yr_bins_figure.py  the figure
    Bin edges snap to whole years and bins under 3.75 yr are dropped, so a
    short tail bin cannot manufacture a large rate.

coastsat/window_convergence/coastsat_window_convergence.py
    Which windows recover the long-term rate? The same OLS on two families of
    NESTED windows, one pinned at each end -- forward walks the END out from
    1996, backward walks the START back from 2024 -- so the pair brackets the
    answer rather than giving one side of it. Both converge on the 1996-2024
    rate BECAUSE that rate is the longest window of each family, so what the
    sweep gives is the window at which the curve stops leaving a tolerance,
    not a match test. Three tolerances are scored side by side rather than one
    being chosen -- but the HEADLINE, what the figures draw and the READMEs
    lead with, is CI OVERLAP: a window passes when its own 95% interval reaches
    the reference's. It is the only one that is not a number somebody picked,
    and it is the only one shaped right -- a funnel, wide where the record is
    short, narrowing onto the reference -- because it asks a five-year fit to
    be indistinguishable from the long-term rate rather than as precise as it.
    It barely moves the answer (27 -> 26 years forward, 25 -> 22 backward), and
    that is the finding: the windows disagree because the shoreline changed,
    not because the fits were noisy. Three scales: `sites/` draws eight transects in full,
    `all_transects/` runs all ~906 and profiles them alongshore, and
    `domain_means/` aggregates that sweep to the unit the model is actually
    graded on -- a groupby, not a refit, because the target is the MEAN of
    transect slopes. Averaging turned out to buy nothing: the domain mean
    settles at the same median window as a single transect, which says the
    window-to-window disagreement is shoreline behaviour and not per-transect
    noise. It draws its own figures -- per-unit panel families, not the
    alongshore profile rates_figures.py covers.

    WHY the windows disagree is answered next door, in
    1-observations/detrended_position/: one island-wide excursion, which on its own
    biases the fitted rate by -0.37 m/yr over 1996-2010 and +0.88 m/yr over
    2010-2024.

coastsat/total_change/coastsat_total_change.py
    A rate turned into a distance, named for the window it was FITTED on:
      --product total      LRR(W) x the years of W      -> total_change/
      --product projected  the 1996-2024 LRR onto a window it was not fitted
                           on                           -> projected/
    The distinction is the whole point of the split; see the docstring.

coastsat/extension/coastsat_extension_lrr.py
    The same fit for the coast beyond GIS 90, for the Pea Island extension
    experiment. Numbers those transects through hat_extension_domains rather
    than the 90-polygon lookup, which stops at 90.

duneline/duneline_endpoint.py
    End dune line minus start line, per 100 m transect and per domain.
    Endpoint and not LRR on purpose: the quantity is net change, and an OLS
    through the intermediate lines would answer a different question. It
    replaced a dune-line LRR product on 2026-09-18.
```
