# wave-climate

Tuning the four wave parameters (Hs, Tp, asymmetry, high-angle fraction) against the CoastSat alongshore rates.

## Start here: the wave climate to use (settled 2026-09-27)

| | 1996–2010 | 2010–2024 | end rates, GIS 1 / GIS 90 (m/yr) |
|---|---|---|---|
| **A, same waves both windows (default)** | Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5 | same | 1996: +4.8394 / +17.545 · 2010: +18.8 / +24.535 |
| B, one change between windows | as A | **Hs 2.5**, rest as A | 2010: +8.0 / +40.399 (1996 as A) |

Both scenarios (natural and full management) use the same row. Raw share of
alongshore variation explained (option A): 1996–2010 managed +20%, natural +24%.
2010–2024 is not fitted by any setting, because of the 2021 CoastSat step
(managed −135%, natural −699%; option B improves this to −122% and −538%).

- The decision and its figures (**read first**): [`2026-09-27-wave-recommendation/README.md`](2026-09-27-wave-recommendation/README.md)
- The end rates as a file: [`../end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/tables/ends.json`](../end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/tables/ends.json).
  Use `ends_m_yr`, which holds the Hs-2 values. The Hs-1 and Hs-2.5 solves are kept under `history`.
  The top-level `reference_waves` key records the first solve (Hs 1) and is not the current setting.
  Pass the ends as `HAT_BE_OVERRIDE="1=<GIS 1>,90=<GIS 90>"` with the edgeBE preset.
- The report, covering the same waves across windows vs period-specific waves and every test: https://claude.ai/artifact/L3ezkgxG7aQNDDPMYazDjG (private)
- **In the code since 2026-09-27:** option A is the default. The ends are in `hatteras_site_config.HATTERAS_BE_EDGE_ONLY` and the waves are the `HAT_hindcast_config` / `hat_run.yaml` defaults, so a plain edgeBE run uses them with no overrides. Option B is recorded as `HATTERAS_WAVE_OPTION_B` in the site config, with the command that runs it (`HAT_HS=2.5` plus `HAT_BE_OVERRIDE`).

Related, outside this theme: the end rates are in [`../end-domain-boundaries/2026-09-27-ends-resolved-metres-offset`](../end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/README.md). The dune line vs shoreline offset test, run at Hs 1 rather than the recommended waves, is in [`../island-offset/2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8`](../island-offset/2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8/README.md).

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-09-24-metres-2-wave-sensitivity`](2026-09-24-metres-2-wave-sensitivity/README.md) | One parameter at a time around Hs 1 / Tp 8 / asym 0.8 / high-angle 0.45, zeroBE. | Hs and high-angle trade off along a ridge; 1996–2010 fits, 2010–2024 does not. | record; superseded by the four-parameter grids |
| [`2026-09-25-wave-grid-smoothed-score`](2026-09-25-wave-grid-smoothed-score/README.md) | All four at once (coarse 108 cells, then refined), zeroBE, scored on the smoothed model. | The zeroBE ranking on the smoothed score; its top 10 per window and scenario became the 09-26 shortlist. | record; its runs are the zeroBE baseline for later studies |
| [`2026-09-26-wave-shortlist-ends-solved`](2026-09-26-wave-shortlist-ends-solved/README.md) | The top settings with the ends solved separately for each. | On the smoothed score the best waves barely move once their ends are solved (1996–2010: Hs 1.25–1.5, asym 0.7, high-angle 0.45); 2010–2024 stays below a flat line. | record |
| [`2026-09-27-wave-grid-fixed-ends`](2026-09-27-wave-grid-fixed-ends/README.md) | The grid with the ends fixed (first the Hs-1 ends, then targeted reruns on the Hs-2 ends), plus one-parameter changes for 2010–2024. | Hs 2 / Tp 7.5 / asym 0.6 / high-angle 0.5 wins on the raw score in both windows; only an Hs change helps 2010–2024. | **current** — the main evidence; read the last section |
| [`2026-09-27-wave-recommendation`](2026-09-27-wave-recommendation/README.md) | Which settings to use, same or per period? | Options A (same) and B (Hs 2.5 in 2010–2024), with figures and the Hs 2.25–3.0 check. | **current — read this first** |
| [`2026-09-28-hs2p5-check-adopted-setup`](2026-09-28-hs2p5-check-adopted-setup/README.md) | Does Hs 2.5 beat Hs 2.0 on the adopted setup (overwash fixes, trimmed storms, dune-cap fix), each on its own ends? | No: tied on the raw score (1996 22.3% vs 22.5%, 2010 −74% vs −77%), for +14 m/yr more at GIS 90. Option A stands. | **current** |

**Status** — **current**: its answer is in use now. **superseded**: a later study
re-asked it; follow the pointer. **record**: a finished check or a result from an
earlier set-up (÷10 offset, Hs 2.5 calibration), kept so the number can be traced.

Back to [the map](../README.md).
