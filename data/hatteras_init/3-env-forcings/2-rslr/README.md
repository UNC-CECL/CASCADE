# rslr — relative sea level rise at the Duck gauge

The rate the model raises sea level by, one value per hindcast window, and the
record it was fitted on. Everything here is written by
`scripts/input_prep/3-env-forcings/2-rslr/duck_rslr_analysis.py`; nothing is
edited by hand except this file.

```
record/
    duck_8651370_meantrend.csv      the NOAA download, untouched
fits/
    duck_rslr_rates.csv             ONE ROW PER WINDOW: the number the model uses
    duck_rslr_timeseries_<w>.csv    the monthly record inside window <w>, with
                                    the fitted trend and the residual
figures/
    duck_rslr_full_record.png/.pdf  the record with the windows and trends
    duck_rslr_windows.png/.pdf      one panel per window
    duck_rslr_residuals.png/.pdf    residuals from each trend
    CAPTIONS.md                     the figure captions (nothing is on the canvas)
```

## The record

NOAA CO-OPS station 8651370, Duck, NC: monthly mean sea level with the average
seasonal cycle removed, in metres relative to the station's most recent MSL
datum, 1978 to 2026. Downloaded 2026-05-04 (the file's date) from the CO-OPS
sea level trends product for the station; the exact URL was not recorded at
the time. The file is kept exactly as it came, four header lines and all, and
the script parses around them.

The record is MSL-datum, not MLLW. Before 2026-09-15 the figures' axis label
said MLLW; the numbers were never affected, only the label.

## The fit

An ordinary least squares line through the monthly values inside each window,
both end years inclusive, so 1984-2004 is fitted on 252 months. The 95%
confidence interval on the slope is t × SE with n − 2 degrees of freedom.

| window | slope (m/yr) | ± 95% CI | n | in the site config |
|---|---|---|---|---|
| 1984-2004 | 0.00391 | 0.00134 | 252 | 0.004 |
| 2004-2024 | 0.00639 | 0.00132 | 249 | 0.006 |
| 1996-2010 | 0.00402 | 0.00225 | 180 | 0.004 |
| 2010-2024 | 0.00651 | 0.00217 | 177 | 0.007 |

`fits/duck_rslr_rates.csv` is the same table at full precision with R², p and
the intercept. Its last column, `config_m_yr`, is the slope rounded to 0.001,
which is how `scripts/site_layer/hatteras_site_config.py` carries it in
`HATTERAS_PERIODS[<start>]["sea_level_rise_rate"]`.

**The config does not read this file.** The four literals are typed there by
hand, and this table is the record of where they came from. Rounding is the
only difference, and it is worth knowing about: 2004-2024 rounds down and
2010-2024 rounds up, so the two look 0.001 m/yr apart in the config and
0.0001 m/yr apart in the gauge record.

## Two pairs of windows

1984-2004 and 2004-2024 partition the record; 1996-2010 and 2010-2024 tile
the later part of it. The pairs overlap each other, which is why the
full-record figure has two panels: each pair is drawn on its own, the earlier
window red and the later blue, the house convention for two vintages.

## Reorganised 2026-09-15

Fifteen files sat flat in this folder: the record, four series, and ten
figures in a presentation palette with the statistics printed on the canvas.
The fitted rates themselves were not written anywhere; they lived in those
annotations and in the config literals. That day the folder was split into
record/fits/figures, the rates table was added, and the figures were redrawn
under the house style (`scripts/site_layer/hat_figure_style.py`) with the four
per-window close-ups and four residual plots folded into one panel figure
each. The per-window series CSVs are unchanged.
