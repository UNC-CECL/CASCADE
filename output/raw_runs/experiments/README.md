# Experiments, by theme

One question each, grouped by the investigation it belongs to (reorganised 2026-09-25, Hannah). Each study folder is `<date>-<what it tested>`; its README or NOTE is the record. Each theme's README has a table of what every study asked, what it found and whether it is still current. A study's folder path under `experiments/` is also its runs' tag in `run_index.csv`.

## Start here: the wave climate to use (settled 2026-09-27)

| | 1996–2010 | 2010–2024 | end rates, GIS 1 / GIS 90 (m/yr) |
|---|---|---|---|
| **A, same waves both windows (default)** | Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5 | same | 1996: +4.8394 / +17.545 · 2010: +18.8 / +24.535 |
| B, one change between windows | as A | **Hs 2.5**, rest as A | 2010: +8.0 / +40.399 (1996 as A) |

Both scenarios (natural and full management) use the same row. Raw share of
alongshore variation explained (option A): 1996–2010 managed +20%, natural +24%.
2010–2024 is not fitted by any setting, because of the 2021 CoastSat step
(managed −135%, natural −699%; option B improves this to −122% and −538%).

- The decision and its figures: [`wave-climate/2026-09-27-wave-recommendation/README.md`](wave-climate/2026-09-27-wave-recommendation/README.md)
- The end rates as a file: [`end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/tables/ends.json`](end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/tables/ends.json).
  Use `ends_m_yr`, which holds the Hs-2 values. The Hs-1 and Hs-2.5 solves are kept under `history`.
  The top-level `reference_waves` key records the first solve (Hs 1) and is not the current setting.
  Pass the ends as `HAT_BE_OVERRIDE="1=<GIS 1>,90=<GIS 90>"` with the edgeBE preset.
- The report, covering the same waves across windows vs period-specific waves and every test: https://claude.ai/artifact/L3ezkgxG7aQNDDPMYazDjG (private)
- **In the code since 2026-09-27:** option A is the default. The ends are in `hatteras_site_config.HATTERAS_BE_EDGE_ONLY` and the waves are the `HAT_hindcast_config` / `hat_run.yaml` defaults, so a plain edgeBE run uses them with no overrides. Option B is recorded as `HATTERAS_WAVE_OPTION_B` in the site config, with the command that runs it (`HAT_HS=2.5` plus `HAT_BE_OVERRIDE`).

## [`island-offset/`](island-offset/README.md)

How the island's planform (the BRIE island offset) is set: its units, and whether it comes from the dune line or the CoastSat shoreline.

- `2026-09-22-shoreline-offset-at-div10-scale` (superseded)
- `2026-09-24-metres-1-offset-units` (current)
- `2026-09-25-offset-source-duneline-vs-shoreline` (current)

## [`wave-climate/`](wave-climate/README.md)

Tuning the four wave parameters (Hs, Tp, asymmetry, high-angle fraction) against the CoastSat alongshore rates.

- `2026-09-24-metres-2-wave-sensitivity` (record)
- `2026-09-25-wave-grid-smoothed-score` (record)
- `2026-09-26-wave-shortlist-ends-solved` (record)
- `2026-09-27-wave-grid-fixed-ends` (current)
- `2026-09-27-wave-recommendation` (current)

## [`end-domain-boundaries/`](end-domain-boundaries/README.md)

The source/sink rates locked at the two end domains (GIS 1 and 90): what they must carry against each target.

- `2026-09-16-end-domains-solved-on-duneline` (superseded)
- `2026-09-18-end-domains-solved-on-redigitized-duneline` (superseded)
- `2026-09-19-end-domains-2010-recheck` (superseded)
- `2026-09-19-end-domains-solved-on-lrr-1996-2024` (superseded)
- `2026-09-27-ends-resolved-metres-offset` (current)
- `2026-09-27-ends-solved-on-duneline-option-a` (current)
- `2026-09-27-ends-solved-on-lrr-1996-2024-option-a` (current)

## [`topography-and-domains/`](topography-and-domains/README.md)

What the model domains contain: dune footprints, rows added or removed, and extending the reach onto Pea Island.

- `2026-09-02-pea-island-row-insert-control` (record)
- `2026-09-08-dune-footprint-behind-road` (record)
- `2026-09-16-pea-island-domain-extension` (record)

## [`code-checks/`](code-checks/README.md)

Whether a code or model change moves the results: re-runs against stored runs, and the Barrier3D route_overwash fix.

- `2026-09-14-calibrated-pair-rerun-current-code` (record)
- `2026-09-14-relocation-arm-rerun-new-code` (record)
- `2026-09-14-relocation-rounding-probes` (record)
- `2026-09-14-site-config-split-check` (record)
- `2026-09-24-metres-3-barrier3d-overwash-fix` (current)

## Chains across themes

- [`2026-09-24-metres-INDEX.md`](2026-09-24-metres-INDEX.md): the 2026-09-24 metres work, steps 1-3 (offset units, wave sensitivity, the Barrier3D fix), and its 2026-09-25 follow-ups.

Renamed on 2026-09-25 (old name → new): 
- `2026-09-22-shoreline-offset` → `island-offset/2026-09-22-shoreline-offset-at-div10-scale`
- `2026-09-16-dune-edgesolve` → `end-domain-boundaries/2026-09-16-end-domains-solved-on-duneline`
- `2026-09-18-dune-edgesolve` → `end-domain-boundaries/2026-09-18-end-domains-solved-on-redigitized-duneline`
- `2026-09-19-edgesolve-2010` → `end-domain-boundaries/2026-09-19-end-domains-2010-recheck`
- `2026-09-19-edgesolve-lrr1996_2024` → `end-domain-boundaries/2026-09-19-end-domains-solved-on-lrr-1996-2024`
- `2026-09-02-pea1989` → `topography-and-domains/2026-09-02-pea-island-row-insert-control`
- `2026-09-08-behindroad-copy` → `topography-and-domains/2026-09-08-dune-footprint-behind-road`
- `2026-09-16-peaisland-ext` → `topography-and-domains/2026-09-16-pea-island-domain-extension`
- `2026-09-14-currency` → `code-checks/2026-09-14-calibrated-pair-rerun-current-code`
- `2026-09-14-paramsplit` → `code-checks/2026-09-14-site-config-split-check`
- `2026-09-14-probe` → `code-checks/2026-09-14-relocation-rounding-probes`
- `2026-09-14-recode` → `code-checks/2026-09-14-relocation-arm-rerun-new-code`
