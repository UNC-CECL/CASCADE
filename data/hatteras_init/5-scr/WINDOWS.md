# Which window is which

*Written 2026-09-22 by `scripts/input_prep/5-scr/tools/windows_index.py` from `site_layer.hat_observed_rates.WINDOW_ROLE`. Edit the dict, not this file.*

Most products under `3-rates/` and `4-comparisons/` are filed one folder per window, and the folders are named only by their years. They are **not** peers.

| window | role | |
|---|---|---|
| `1996_2010` | Calibration period | **current chain** — the model is fitted here |
| `2010_2024` | Test period (held out) | **current chain** — held out; nothing is fitted to it |
| `1996_2024` | Full period (context, not graded) | **current chain** — spans the whole current chain; CONTEXT ONLY, no run is graded against it |
| `1984_2004` | Calibration period, legacy chain | legacy chain — the 1984-start chain, superseded as the main chain 2026-09 |
| `2004_2024` | Test period, legacy chain | legacy chain — the 1984-start chain, superseded as the main chain 2026-09 |

The **current chain** is 1996 → 2010 → 2024: the model is fitted on the first half and the second is held out. `1996_2024` spans both and exists for context only — no run is graded against it, and a number read from it is not a model result.

The **legacy chain** is 1984 → 2004 → 2024. Both of its periods are still live in `hatteras_site_config.HATTERAS_PERIODS` and runs on them still exist, but it was superseded as the main chain in September 2026. A figure in a `1984_2004/` or `2004_2024/` folder is not the current answer to anything unless you meant to ask about that chain.

## Which product covers which window

Read off the disk, so it is what is actually there:

| product | 1996–2010 | 2010–2024 | 1996–2024 | 1984–2004 | 2004–2024 |
|---|---|---|---|---|---|
| `3-rates/coastsat/lrr` | yes | yes | yes | yes | yes |
| `3-rates/coastsat/endpoint` | yes | yes | yes | yes | yes |
| `3-rates/coastsat/total_change` | yes | yes | yes | — | — |
| `3-rates/coastsat/projected` | yes | yes | — | — | — |
| `3-rates/coastsat/5yr_bins` | yes | yes | yes | — | — |
| `3-rates/duneline/endpoint` | yes | yes | yes | yes | yes |
| `4-comp/coastsat_endpoint_vs_duneline_endpoint` | yes | yes | yes | yes | yes |
| `4-comp/coastsat_total_change_vs_duneline_endpoint` | yes | yes | yes | — | — |
| `4-comp/coastsat_projected_vs_duneline_endpoint` | yes | yes | — | — | — |

The gaps are deliberate, not missing work:

- `projected/` has no `1996_2024`: there the rate window IS the change window, so the answer is `total_change/1996_2024` and building it twice under two names is the confusion the 2026-09-21 rename removed. See `3-rates/README.md` for the total / projected / observed vocabulary.
- `projected/`, `total_change/` and `5yr_bins/` cover the current chain only. They were built after it became the main chain. The two `*_projected_vs_duneline_endpoint` windows are the two halves: over the full period the rate window IS the change window, so that case is the `total_change` product.
- `lrr/`, both `endpoint/` products and `coastsat_endpoint_vs_duneline_endpoint/` cover all five, because they predate the switch.
