# 5-scr — observed shoreline change

> **Lost?** [`FIGURES.md`](../../../FIGURES.md) is the one-page index of which figure answers which question. [`WINDOWS.md`](WINDOWS.md) says which window is which — two chains, and one context window nothing is graded against.

Where the observations live: the satellite shoreline record, the rate fits the
model is graded against, and the layers that tie transects to Barrier3D
domains.

## Layout

The folders are grouped by what they are for, numbered in the order the work
runs (2026-09-18):

```
1-observations/              measured or digitized; nothing here is fitted by us
    coastsat_timeseries/     raw per-transect chainage, one folder per CoastSat site
    dsas_1978_2019/          DSAS rates: a different source, different transects
    shoreline_inventory/     study-area and reference shorelines
2-transect-frame/
    transect_domains/        the transect-to-domain lookup, the transect layer,
                             the domain polygons, and the verification set
3-rates/                     THE MODEL TARGETS LIVE HERE: tables plus one
                             house-style figure per window (rates_figures.py).
                             Grouped by source since 2026-09-18; README.md
                             there indexes every product
    coastsat/
        lrr/<window>/        the OLS rate fits: 1984_2004 1996_2010
                             2004_2024 2010_2024, and 1996_2024 (CONTEXT
                             only, not a model target, below)
        endpoint/<window>/   NET CHANGE between ±6-month means at the
                             dune-line dates, m and m/yr (2026-09-18)
        5yr_bins/<window>/   the OLS in successive 5-year bins, 1996_2010
                             2010_2024 1996_2024 (rebuilt 2026-09-18)
    duneline/
        endpoint/<window>/   NET CHANGE between the two DUNE LINES bounding
                             each window, m and m/yr (2026-09-18; replaced
                             duneline_lrr/, an OLS through every line in the
                             window -- Hannah: "we are tracking net change")
4-comparisons/                 one folder per QUESTION (reorganized 2026-09-19)
    shoreline_vs_duneline/   does the dune line move with the shoreline? Every
                             alongshore figure in metres with the beach-width gap
        net_change/<window>/     observed net change on both sides, per window
                                 (1984_2004 ... 1996_2024), stacked figure, chains/
        projected/1996_2024/     the shoreline as the 1996-2024 LRR projected over
                                 the dune-line interval
    duneline_positions/      WHERE the 1997/2009/2023 dune lines sat: an
                             island overview, imagery zooms, and the dune
                             line's distance to NC-12 and to the shoreline
archive/                     kept, not for use
    coastsat_lrr_superseded_20260810/   retired windows (1978-1997, 1997-2019);
                                        6-scr-smooth's DSAS comparison still reads them
    coastsat_5yr_bins/                  the three earlier 5-year-bin runs (the
                                        1984-2004 / 2004-2024 framing)
    coastsat_lrr_quicklooks_20260918/   the autoscaled bar and scatter PNGs that
                                        used to sit beside each LRR fit
    duneline_lrr_retired_20260918/      the dune-line OLS, retired for endpoint
    2026-09-19_4-comparisons_duplicates/ coastsat_windows/ and duneline_windows/:
                                        their figures duplicated 3-rates (window
                                        and chain figures); the 2 x 2 moved to
                                        3-rates/coastsat/lrr/lrr_four_windows
    rodanthe_plots/                     poster figures; the script writes to
                                        output/figures/shoreline/ now
```

Until 2026-09-18 all of these sat side by side at the top of `5-scr/`. Four
were renamed in the move:

| now | was |
|---|---|
| `1-observations/dsas_1978_2019/` | `scr-dsas-1978-2019/` |
| `3-rates/coastsat/5yr_bins/` | `coastsat_timeseries_lrr/` |
| `4-comparisons/coastsat_windows/` | `coastsat_lrr_windows/` (archived 2026-09-19) |
| `4-comparisons/trajectory_patterns/` | `shoreline_change_patterns/` (output deleted 2026-09-19) |

On 2026-09-19 `4-comparisons/` went from eight folders to two (Hannah):
`duneline_vs_coastsat/`, `net_change_1996_2024/` and `projected_vs_duneline/`
merged into `shoreline_vs_duneline/` (and on 2026-09-21 `projected_vs_duneline`
was absorbed again, into `total_change/`, when the vocabulary was settled); `coastsat_windows/` and
`duneline_windows/` archived as duplicates of `3-rates/`;
`two_period_comparison/` (03-31) and `trajectory_patterns/` (06-09) deleted as
stale, on the old 1984/2004 periods (their scripts still exist and would
regenerate them).

### The geojsons are not in git

The four geojsons here (the 80 MB `CoastSat_transect_layer.geojson`,
`nc_shorelines.geojson` at 20 MB, `cascade_area.geojson` and
`wet_dry_groin.geojson`) stay on disk but off GitHub (2026-09-18;
`.gitignore`). **A fresh clone does not have them.** Anything that reads the
transect layer or the reference shorelines needs a copy from someone who has
one. The lookup CSVs, `HAT_domains.json` and the fitted rates are tracked,
so the model runs without them.

## A window, not a year

`3-rates/coastsat/lrr/1996_2010/` is named for an **interval**, because a rate
fit spans one. A survey — a dune line, a road alignment — is a moment, and is
named for its year instead. That split is the naming rule the whole init tree
follows.

## Do not hardcode these paths

Resolve them through **`scripts/site_layer/hat_observed_rates.py`**:

```python
import sys; sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_observed_rates import lrr_csv, transect_lookup, windows

lrr_csv(1996, 2010)     # raises, naming the windows on disk, if absent
```

Every folder above has a constant there (`COASTSAT_LRR_ROOT`, `DSAS_ROOT`,
`TRANSECT_LAYER`, `DOMAIN_BOXES`, ...). Since 2026-09-18 every script reads
5-scr through the resolver, so this regrouping changed that one file and no
script. Keep it that way: a new script that types `"5-scr" / ...` is the
thing that makes the next move expensive.

`transect_lrr_full.csv` is a **model input**: section 8 of the hindcast runner
reads it on every run, and the calibrated source/sink preset is fitted against
it. It is not a figure product.

## Moved here 2026-09-12

All of this lived under `scripts/input_prep/5-scr/CoastSat/` until then, which
put a model input in the code tree — so anyone archiving `data/hatteras_init/`
shipped a model that could not run.

The **producers stayed** in `scripts/input_prep/5-scr/`:

| producer (`scripts/input_prep/5-scr/`) | writes to (`data/hatteras_init/5-scr/`) |
|---|---|
| `CoastSat/coastsat_domain_mapping.py` | `2-transect-frame/transect_domains/` |
| `CoastSat/coastsat_domain_lrr_fixed.py` | `3-rates/coastsat/lrr/<window>/` |
| `CoastSat/coastsat_extension_lrr.py` | `3-rates/coastsat/lrr/<window>/ext/` |
| `CoastSat_timeseries/coastsat_lrr_5year_bins.py` | `3-rates/coastsat/5yr_bins/` |
| `duneline_endpoint/duneline_endpoint.py` | `3-rates/duneline/endpoint/<window>/` |
| `coastsat_endpoint/coastsat_endpoint.py` | `3-rates/coastsat/endpoint/<window>/` |
| `coastsat_total_change/coastsat_total_change.py` | `3-rates/coastsat/total_change/<window>/`; `--product projected` -> `3-rates/coastsat/projected/<window>/` (was `coastsat_lrr_projected/`) |
| `net_change/net_change_1996_2024.py` | `4-comparisons/shoreline_vs_duneline/net_change/chains/` |
| `total_change_vs_duneline/total_change_vs_duneline.py` | `4-comparisons/shoreline_vs_duneline/total_change/` (was `lrr_net_change/`; absorbed `projected/` 2026-09-21) |
| `duneline_positions/duneline_positions.py` | `4-comparisons/duneline_positions/` |
| `CoastSat/coastsat_lrr_windows.py` | `3-rates/coastsat/lrr/lrr_four_windows` (the 2 x 2 only, since 2026-09-19) |
| `duneline_vs_coastsat/duneline_vs_coastsat.py` | `4-comparisons/shoreline_vs_duneline/net_change/<window>/` |
| `shoreline_change_patterns/` | `4-comparisons/trajectory_patterns/` (not current: old periods) |
| `CoastSat/coastsat_two_period_comparison.py` | `4-comparisons/two_period_comparison/` (not current: old periods) |
| `shoreline_inventory/HAT_shoreline_inventory.py` | `1-observations/shoreline_inventory/` |

## Rebuilding a window

```
python scripts/input_prep/5-scr/CoastSat/coastsat_domain_lrr_fixed.py \
    --start-year 1996 --end-year 2010
```

Writes `transect_lrr_full.csv`, `domain_lrr_summary.csv` and two figures into
`3-rates/coastsat/lrr/<start>_<end>/`. Re-running an existing window
overwrites it in place; the fit is deterministic, so that is reproducible
rather than destructive.

## The long window, 1996-2024 (context, not a target)

`3-rates/coastsat/lrr/1996_2024/` is the long-term rate over the whole
canonical chain, built 2026-09-18 with the same script and schema as the
other windows (`--start-year 1996 --end-year 2024`), plus the Pea Island
extension under `ext/` (`coastsat_extension_lrr.py`, same arguments).
**No run is graded against it.** The model chain is still 1996→2010 and
2010→2024, graded window by window. `lrr_csv(1996, 2024)` resolves it like
any window, so a script that takes a window argument can read it. Nothing
iterates the windows on disk, so it cannot be picked up as a target by
accident.

Nourishment is left in: Rodanthe 2014 sits mid-window, and Buxton and Avon
2022 sit two years from its end. A fill is a step, which a single slope fits
poorly, so read those domains as "includes placed sand". Nothing is masked.

Its figures are `3-rates/coastsat/lrr/1996_2024/lrr_1996_2024.png` and the
chain figure `3-rates/coastsat/lrr/chains/lrr_chain_1996_2010_2024.png`. Until
2026-09-19 it was `4-comparisons/coastsat_windows/1996_2024/lrr_1996_2024_halves.png`
(`coastsat_lrr_windows.py --overlay 1996_2024`, now retired and archived), two stacked panels: (a) the
long window filled by sign, with the model-input fill footprints as bars
above it; (b) 1996-2010 (grey) and 2010-2024 (black) as plain lines. The
shoal zones (Avon, Wimble) are faint amber-hatched, outlined boxes behind
the data in both panels, only to show where they are. The y
axis is the tightest whole metre that holds every line (±7 m/yr), NOT the
shared ±8 of the single-window figures; the caption says so. The three
means are side by side in `supporting/lrr_1996_2024_halves.csv`. The figure is also published to
`output/figures/shoreline/`.

## Comparing windows

`3-rates/coastsat/lrr/lrr_four_windows.png` is the one place the four windows
are drawn against each other (in `4-comparisons/coastsat_windows/` until
2026-09-19; since then `coastsat_lrr_windows.py` draws only this 2 x 2, the
per-window figures being `rates_figures.py`'s). The text below describes it as it was. `scripts/input_prep/5-scr/CoastSat/coastsat_lrr_windows.py`
reads each `domain_lrr_summary.csv` through the resolver, pins the y axis at
the largest |mean| over all of them plus 1 m, rounded up to the metre
(written to `supporting/y_bounds.txt`), and writes `lrr_four_windows`, a
2 x 2 with the 1984-start period in the left column and the 1996-start period
in the right, plus one figure per window in its own `<start>_<end>/` folder
(the `duneline_vs_coastsat/` layout, 2026-09-18; they sat flat before). See
that folder's README. The per-window `domain_lrr_bar.png`
autoscales, so it is not the figure to compare across windows.

## Dune line against the shoreline

`4-comparisons/shoreline_vs_duneline/README.md` is the map (since 2026-09-19);
`4-comparisons/shoreline_vs_duneline/net_change/README.md` (was
`4-comparisons/duneline_vs_coastsat/README.md`) is the methods report: the
three rates, the imagery date behind each line, what is inside each endpoint
window, and the results table. `4-comparisons/shoreline_vs_duneline/net_change/<start>_<end>/`
differences two digitized dune lines (`2-brie-offset/raw_offsets/`, read the
way the hindcast loader reads them) and puts the per-domain rate against the
CoastSat shoreline two ways: the LRR already in
`3-rates/coastsat/lrr/<start>_<end>/`, and an endpoint rate from the mean
CoastSat position in a one-year window about each survey date. Seaward
positive throughout. Built by
`scripts/input_prep/5-scr/duneline_vs_coastsat/duneline_vs_coastsat.py`;
the survey dates are in each folder's `supporting/PROVENANCE.md`, with
the tables, PDFs and captions (a figure folder shows figures). `--grid`
redraws every window on disk as `alongshore_four_windows.png` on one y
axis: four full-width panels stacked, the 1984-start pair above the
1996-start pair, each as wide as a single-window figure (Hannah, 2026-09-15:
the 2 x 2 was too small to read). `--layout grid` gives the 2 x 2 by period
that `3-rates/coastsat/lrr/lrr_four_windows` uses. The 1984 line is
the 1984-09-19 photo and the 2004 line the 2004-05-25 Google Earth capture;
neither file carries its date, the script's `KNOWN_SURVEY_DATES` does.
