# 2026-09-28 (evening) — managed runs before the dune-cap fix

The 12 matrix runs with the beach/dune manager on (full management with and
without relocation, full management without fill, beach/dune only; edgeBE and
zeroBE; 1996-2010 and 2010-2024), moved here on 2026-09-28 before they were re-run.

They ran on the adopted setup (Barrier3D `hatteras/adopted`, `v3_trim24` storms,
ends 1996 +4.3509/+19.0935, 2010 +8.0/+22.4937). CASCADE's `filter_overwash`
still clipped the WHOLE dune cell to 4 m above the berm every year. That deleted
natural village dune above the cap: 24 domains in 1996-2010 and 17 in 2010-2024,
about 150,000 m3 of starting dune each time.

Hannah, 2026-09-28: "go with option 1, clip only the bulldozed sand". Since then
the cap limits only the sand the manager adds (`cascade/beach_dune_manager.py`,
`DUNE_CAP_APPLIES_TO = "added sand only"`, recorded in run metadata as
`bdm_dune_cap_applies_to`). The analysis is in `output/comparisons/adoption_2026-09-28/`.

`manifests/driver_manifest.jsonl` is the driver manifest before their entries
were removed. Its keys do not carry the CASCADE change, so left in place the
driver would have skipped the re-runs. The natural and roadway-only runs did
not change and stay in `matrix/`.

Kept for comparison. Do not use for analysis.

## `matrix_2010_edgeBE_old_gis90_end/`

After the fix, the 2010-2024 GIS 90 end moved from +22.4937 to +21.2582 m/yr
(experiments/end-domain-boundaries/2026-09-28-ends-resolved-dunecap/). Every
2010 edgeBE run is therefore re-run. Three runs on the old end are kept here:
- natural and roadway-only, which the cap fix itself does not touch
- the full-management run made WITH the fix at the old end, which was the
  residual check that showed GIS 90 at +0.185

1996 held its ends, so its natural and roadway-only runs stay in `matrix/`.
