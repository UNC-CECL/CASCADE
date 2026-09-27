# Experiments, by theme

One question each, grouped by the investigation it belongs to (reorganised 2026-09-25, Hannah). Each study folder is `<date>-<what it tested>`; its README or NOTE is the record. A study's folder path under `experiments/` is also its runs' tag in `run_index.csv`.

## [`island-offset/`](island-offset/README.md)

How the island's planform (the BRIE island offset) is set: its units, and whether it comes from the dune line or the CoastSat shoreline.

- `2026-09-22-shoreline-offset-at-div10-scale`
- `2026-09-24-metres-1-offset-units`
- `2026-09-25-offset-source-duneline-vs-shoreline`

## [`wave-climate/`](wave-climate/README.md)

Tuning the four wave parameters (Hs, Tp, asymmetry, high-angle fraction) against the CoastSat alongshore rates.

- `2026-09-24-metres-2-wave-sensitivity`
- `2026-09-25-wave-grid-smoothed-score`
- `2026-09-26-wave-shortlist-ends-solved`
- `2026-09-27-wave-grid-fixed-ends`
- `2026-09-27-wave-recommendation`

## [`end-domain-boundaries/`](end-domain-boundaries/README.md)

The source/sink rates locked at the two end domains (GIS 1 and 90): what they must carry against each target.

- `2026-09-16-end-domains-solved-on-duneline`
- `2026-09-18-end-domains-solved-on-redigitized-duneline`
- `2026-09-19-end-domains-2010-recheck`
- `2026-09-19-end-domains-solved-on-lrr-1996-2024`
- `2026-09-27-ends-resolved-metres-offset`

## [`topography-and-domains/`](topography-and-domains/README.md)

What the model domains contain: dune footprints, rows added or removed, and extending the reach onto Pea Island.

- `2026-09-02-pea-island-row-insert-control`
- `2026-09-08-dune-footprint-behind-road`
- `2026-09-16-pea-island-domain-extension`

## [`code-checks/`](code-checks/README.md)

Whether a code or model change moves the results: re-runs against stored runs, and the Barrier3D route_overwash fix.

- `2026-09-14-calibrated-pair-rerun-current-code`
- `2026-09-14-relocation-arm-rerun-new-code`
- `2026-09-14-relocation-rounding-probes`
- `2026-09-14-site-config-split-check`
- `2026-09-24-metres-3-barrier3d-overwash-fix`

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
