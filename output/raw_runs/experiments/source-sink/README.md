# source-sink

The per-domain background erosion (BE) field: what each domain needs so the calibration run ends where CoastSat says the shoreline went, and whether that carries into the test.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-10-05-be-domain-solve-1996_2009`](2026-10-05-be-domain-solve-1996_2009/README.md) | What BE does each of the 90 domains need on the 1996–2009 calibration (BE set 1)? | 10 passes, interior RMSE 0.49 m; field −4.85 to +12.43 m/yr (`fields/be_field_step10.csv`). Smoothing the field makes the fit worse (0.49 → 4.33 m). | **current**; the `domainBE` preset (1996 and 2009) |
| [`2026-10-06-test-target-and-the-2021-step`](2026-10-06-test-target-and-the-2021-step/README.md) | How much of the 2009–2025 test misfit is the island-wide 2020→2021 CoastSat step? | The bias is the step: domainBE bias −18.2 → −1.1 m (no step), −1.4 (pre-step window). RMSE does not improve: the pattern doesn't transfer. | record; set 2 not derived |

**Status** — **current**: its answer is in use now. **record**: a finished check, kept so the number can be traced.

Back to [the map](../README.md).
