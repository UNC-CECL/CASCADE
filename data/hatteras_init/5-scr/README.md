# 5-scr — observed shoreline change

Where the observations live: the satellite shoreline record, the rate fits the
model is graded against, and the layers that tie transects to Barrier3D
domains.

```
coastsat_timeseries/        raw per-transect chainage, one folder per CoastSat site
transect_domains/           the transect-to-domain lookup, the transect layer,
                            the domain polygons, and the verification set
coastsat_lrr/               rate fits, ONE FOLDER PER WINDOW
    1984_2004/  1996_2010/  2004_2024/  2010_2024/
    old_time_periods/       retired windows, kept for comparison, not for use
    two_period_comparison/  rodanthe_plots/  old_dsas_comparisons/  custom/
coastsat_timeseries_lrr/    the 5-year-bin fits
shoreline_inventory/        study-area and reference shorelines
shoreline_change_patterns/  trajectory classification output
scr-dsas-1978-2019/         DSAS rates — a different source, different transects
```

## A window, not a year

`coastsat_lrr/1996_2010/` is named for an **interval**, because a rate fit
spans one. A survey — a dune line, a road alignment — is a moment, and is named
for its year instead. That split is the naming rule the whole init tree
follows.

## Do not hardcode these paths

Resolve them through **`scripts/hat_observed_rates.py`**:

```python
import sys; sys.path.insert(0, str(REPO / "scripts"))
from hat_observed_rates import lrr_csv, transect_lookup, windows

lrr_csv(1996, 2010)     # raises, naming the windows on disk, if absent
```

`transect_lrr_full.csv` is a **model input**: section 8 of the hindcast runner
reads it on every run, and the calibrated source/sink preset is fitted against
it. It is not a figure product.

## Moved here 2026-09-12

All of this lived under `scripts/input_prep/5-scr/CoastSat/` until then, which
put a model input in the code tree — so anyone archiving `data/hatteras_init/`
shipped a model that could not run. Twenty-two files built that path by hand;
that is why there is a resolver now.

The **producers stayed** in `scripts/input_prep/5-scr/`, one folder per
product, and each output folder here is named for the producer folder it came
from so the pairing is readable in both directions.

## Rebuilding a window

```
python scripts/input_prep/5-scr/CoastSat/coastsat_domain_lrr_fixed.py \
    --start-year 1996 --end-year 2010
```

Writes `transect_lrr_full.csv`, `domain_lrr_summary.csv` and two figures into
`coastsat_lrr/<start>_<end>/`. Re-running an existing window overwrites it in
place; the fit is deterministic, so that is reproducible rather than
destructive.
