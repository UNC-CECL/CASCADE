# 2026-09-29 — matrix runs before the split12 storms

The 22 no-groin matrix runs (edgeBE and zeroBE, 1996–2010 and 2010–2024), moved here on 2026-09-29 before they were re-run.

**What they ran on:**
- Barrier3D `hatteras/adopted` (per-cell dune ceilings, overwash fixes)
- the beach/dune cap limited to added sand
- option A waves
- ends 1996 +4.3509/+19.0935 and 2010 +8.0/+21.2582
- storms **`v3_trim24`**

Hannah, 2026-09-29: "use split12 as the storm series going forward", then "re-run the matrix and re-solve the ends". The storm series is now `v3_split12_trim24`, in which grouped events are split at ≥12 h below the berm, so Fran 1996 and Jose 2017 are back. See `data/hatteras_init/3-env-forcings/3-storms/PROVENANCE.md` and `experiments/storms-and-overwash/2026-09-29-event-splitting/`.

`manifests/driver_manifest.jsonl` is the driver manifest (43 matrix_nogroin entries). Its keys do not carry the storm series, so left in place the driver would have skipped the re-runs.

`matrix/figures/` was not moved; it is regenerated from the new runs.

Kept for comparison (trim24 against split12 on the full matrix). Do not use for analysis.
