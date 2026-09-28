# 2026-09-28 — the edgeBE matrix on the LOESS-10 option A ends

The 11 edgeBE no-groin matrix runs made on 2026-09-27 (commit 2ea0d140) at the option A waves, with the end rates solved against the **LOESS-10** CoastSat target:

| period | GIS 1 | GIS 90 |
|---|---|---|
| 1996–2010 | +4.8394 | +17.545 |
| 2010–2024 | +18.8 | +24.535 |

They were moved here intact on 2026-09-28, before the matrix was re-run. That day the runner's target went to LOESS 7, and the ends were re-solved and adopted: 1996 +4.8394 / +18.2545, 2010 +18.8657 / +24.2358 (`end-domain-boundaries/2026-09-28-ends-resolved-loess7/`). Hannah: "rerun the matrix with the new ends".

The zeroBE matrix runs were **not** moved. They impose no end rates, so the new ends do not change them, and they stay in `matrix/`.

- `manifests/driver_manifest.jsonl` is the driver manifest as it stood. Its edgeBE rows carry the old source/sink digest, so the re-run did not skip them.
- `manifests/run_index_snapshot.csv` is the run index before the move.

Kept for comparison. Do not use for analysis.
