# end-domain-boundaries

The source/sink rates locked at the two end domains (GIS 1 and 90): what they must carry against each target.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-09-16-end-domains-solved-on-duneline`](2026-09-16-end-domains-solved-on-duneline/NOTE.md) | What would the ends carry if solved on the dune line instead of CoastSat? | Ends shrink 2–3×; GIS 1 flips sign depending on the reading. CoastSat stays the target. | superseded by the 09-18 re-solve (re-digitized lines) |
| [`2026-09-18-end-domains-solved-on-redigitized-duneline`](2026-09-18-end-domains-solved-on-redigitized-duneline/NOTE.md) | The same, on the 09-18 re-digitized dune lines. | Same reading; the dune-solved runs are a model set in the model-vs-observed figures. | superseded by `2026-09-27-ends-solved-on-duneline-option-a` |
| [`2026-09-19-end-domains-2010-recheck`](2026-09-19-end-domains-2010-recheck/NOTE.md) | Do the config's 2010 ends (+72.6 / +31.3) still close after the re-digitizing? | Yes, to within 0.02 m/yr. | superseded by `2026-09-27-ends-resolved-metres-offset` |
| [`2026-09-19-end-domains-solved-on-lrr-1996-2024`](2026-09-19-end-domains-solved-on-lrr-1996-2024/NOTE.md) | The ends against the full 1996–2024 LRR instead of each window's own. | 1996: +28.5 / +24.5; 2010: +37.1 / +25.4. | superseded by `2026-09-27-ends-solved-on-lrr-1996-2024-option-a` |
| [`2026-09-27-ends-resolved-metres-offset`](2026-09-27-ends-resolved-metres-offset/README.md) | The ends under the metres offset, at the wave climate in use. | At Hs 2 / Tp 7.5 / asym 0.6 / high-angle 0.5: 1996 +4.8394 / +17.545; 2010 +18.8 / +24.535 (`tables/ends.json`). | **current**; in the site config since 2026-09-27 |
| [`2026-09-27-ends-solved-on-duneline-option-a`](2026-09-27-ends-solved-on-duneline-option-a/NOTE.md) | The ends solved on the dune line under option A (for target_comparison and model_vs_observed). | mean3: 1996 −3.0 / +7.6, 2010 +3.4 / +15.0. Interior RMSE barely moves (1.05 → 1.08, 2.25 → 2.28 vs CoastSat). | **current** |
| [`2026-09-27-ends-solved-on-lrr-1996-2024-option-a`](2026-09-27-ends-solved-on-lrr-1996-2024-option-a/NOTE.md) | The ends solved on the 1996–2024 CoastSat LRR under option A (target_comparison/projected). | 1996 +4.5 / +27.5, 2010 +4.8 / +20.5. | **current** |

**Status** — **current**: its answer is in use now. **superseded**: a later study
re-asked it; follow the pointer. **record**: a finished check or a result from an
earlier set-up (÷10 offset, Hs 2.5 calibration), kept so the number can be traced.

Back to [the map](../README.md).
