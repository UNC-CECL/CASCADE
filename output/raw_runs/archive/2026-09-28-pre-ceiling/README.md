# 2026-09-28 — the matrix before the dune ceiling and storm adoption

All 22 matrix runs (edgeBE and zeroBE, 1996–2010 and 2010–2024) as they stood on 2026-09-28. They were made on:

| | |
|---|---|
| Barrier3D | `fix/route-overwash-axis-swap` @ 49fd069 |
| dune ceiling | `Dmaxel` never set, so the 3.4 m NAVD88 default (3.04 m MHW) applied to every dune |
| storms | `v3_72`: events over 72 h dropped |
| end rates | the LOESS-7 option A ends: 1996 +4.8394 / +18.2545, 2010 +18.8657 / +24.2358 |

**Moved here intact on 2026-09-28, before the matrix was re-run on the adopted setup** (Hannah: "include the overwash fixes, keep option A, go ahead"). The adopted setup is:

- **Barrier3D** `hatteras/adopted` (the three overwash fixes and per-cell dune ceilings from the starting dunes; `HATTERAS_FIXES.md` on that branch).
- **Storms** `v3_trim24`: every event kept, trimmed to 24 h around its peak.
- **End rates** re-solved on that setup.

The evidence is in `experiments/storms-and-overwash/` (six studies) and `experiments/code-checks/2026-09-28-*`. The default setup reproduces those experiments exactly (`code-checks/2026-09-28-adoption-end-to-end/`).

**Why every run moved, not only edgeBE.** The dunes and storms change every scenario, including zeroBE.

- `manifests/driver_manifest.jsonl` is the driver manifest as it stood. Its keys carry neither the Barrier3D version nor the storm file, so left in place it would have made `HAT_run_all.py` skip the re-run.
- `manifests/run_index_snapshot.csv` is the run index before the move.

To reproduce one of these runs, check out `fix/route-overwash-axis-swap` in the Barrier3D repository and use the parameter template and storm variant (`v3_72`) of commit 37704b2d or earlier.

Kept for comparison. Do not use for analysis.

## The wave sensitivity sweep on this matrix (added 2026-09-28, 19:15)

`sensitivity/` holds the 48 one-at-a-time wave cells run 16:24–18:57 around
the edgeBE full-management runs in `matrix/` (Hs, Tp, asymmetry, high-angle;
both windows), with their figures (`sensitivity/figures/<window>_edgeBE/`),
driver logs (`sensitivity/logs/`) and manifests
(`manifests/sensitivity_{1996,2010}.jsonl`). They ran on the same setup as the
matrix here, so they moved with it; `sensitivity/README.md` has the values and
what they showed. The sweep is re-run on the adopted setup into
`raw_runs/sensitivity/`.

**On disk only.** A cell's name appears twice in its file paths, and under
this archive the longest reaches 288 characters, past what git on Windows
accepts without `core.longpaths`. The 153 of their files that had been staged
were unstaged before the move. The figures were drawn against `matrix/` here
before it moved; `plot_sensitivity.py` will not redraw them now, because it
looks for its baselines among the live matrix runs.
