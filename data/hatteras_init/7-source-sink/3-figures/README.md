# 3-figures — the source/sink calibration, in pictures

Five figures, grouped by the question each answers rather than by the script
that drew it. Each folder carries its own `CAPTIONS.md` beside the figures, as
everywhere else in the data tree.

**One folder per calibration pair** (since 2026-09-18). The five below are in
`1984_2004__2004_2024/`, the pair that was fitted, iterated and exported.
`1996_2010__2010_2024/` has `1-field/` only: that pair's fit writes its two
field figures, but the method and limits figures are drawn for the default pair
alone. Before 2026-09-18 the default pair's folders sat directly in
`3-figures/`.

| folder | figure | the question it answers |
|---|---|---|
| `1-field/` | `fig_be_rates` | What is the calibrated field? The per-domain source/sink rate for both periods, and the three forecast scenarios built from them. |
| | `fig_be_diagnostic` | How well does the model match the record? Observed versus modelled shoreline rate, and the residual between them, per period. |
| `2-method/` | `fig_be_zones_and_corrections` | Which domains were eligible at all, and how much of each rate came from the one-shot solve versus the iteration. |
| | `fig_be_convergence` | Did the fixed-point solve converge, and was the zone set fixed before it ran rather than grown to fit? |
| `3-limits/` | `fig_groin_reserved_residual` | Why the largest residual in the hindcast, at D5-D7, is deliberately left uncorrected. |

## Why `3-limits` holds one figure

Because that figure is the one most often misread. At convergence the
calibration leaves its biggest misfit at D6 — about 2.3 m/yr in period 1 and
−2.9 in period 2, roughly twice the next worst domain — and on its face that
looks like a failure. It is not: D5-D7 are the Buxton groin's footprint, and
the residual there is the groin module's shortfall, not a background-erosion
term. Letting the source/sink absorb it would close the same gap twice and make
the M/f fit unfalsifiable.

Filed beside the method figures it reads as another piece of supporting
evidence. On its own it reads as what it is — a stated boundary on what this
calibration claims to do.

## Both formats

Every figure is written as PNG and PDF with the same stem. Note that `*.png` is
gitignored repository-wide, so **only the PDFs are tracked** — a figure with no
PDF exists on one machine only.

## Regenerating

From `scripts/input_prep/7-source-sink`:

```
HAT_BE_BASE_PRESET=calibBE python 2-calibrate/be_zone_residual_fit.py   # 1-field/
python 3-figures/plot_be_convergence.py                                 # 2-method/
python 3-figures/plot_be_zones.py                                       # 2-method/
python 3-figures/plot_groin_reserved_residual.py                        # 3-limits/
```

The two `1-field/` figures come out of the calibration script itself; the other
three are standalone and read the live config, so they cannot drift from the
field they document. Stage 4 refuses to export figures older than the newest
calibBE run, so regenerate these before re-exporting.

Style: `scripts/site_layer/hat_figure_style.py`. No in-image titles or footnotes; the words
are in each folder's `CAPTIONS.md`.
