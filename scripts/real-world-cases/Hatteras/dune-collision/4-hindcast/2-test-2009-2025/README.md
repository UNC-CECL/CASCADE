# 2-test-2009-2025 — does the 1996-2009 calibration hold in 2009-2025?

`HAT_collision_test_2009.py`. Nothing refit: the hindcast's own defaults with the module as calibrated
on 1996-2009. Cs 1.24e-4 and a 45 m beach limit to the crest until 2026-11-03; the dry beach, Cs 1.02e-4 and
25 m since (`*_dry` variants). The script header lists the variants (control, collision, nocape and its growth,
Hs, beach and road variants).
Runs, scores and notes: `output/raw_runs/dune-storm-response/2026-10-29-collision-test-2009-2025/`
(`NOTE.md`, `tables/scores.csv`).

Added 2026-11-03 to 2026-11-04 (each in the script with its reason):

- `control_rd198`: the module off with the 198 m road rule, the pair of `nocape_rd198_dry`.
- Rodanthe village (GIS 77-83): `control_rd198_ow30`, `_ow59` (30%, 59% of overwash bulldozed onto the dune),
  `control_rd198_upkeep` (dunes rebuilt to their 2009 height, the sand off the beach).
- `*_vf40`, `_vf65`, `_vf90`: every village returning that share of its overwash to the shoreface. 90% at
  Rodanthe alone became the default (2026-11-04); `RODANTHE_40` holds it at 40% for the steps below.
- The model-change steps: `step1_base`, `step2_module`, `step3_road` (per-cell toes) and `step2_module_dt`,
  `step3_road_dt`, `nocape_rd198_dry_toeavg` (one toe and one beach factor per domain from the 2004, 2005 and
  2009 lidar, `TOE_DOMAIN_AVG`; nothing after the start). `nocape_rd198_dry_toeavg_cw` is the same with the
  per-cell beach widths (2011-2020 survey medians). `nocape_rd198_dry_ts`: toes and widths from the 2009 start
  profile.

Added 2026-11-04 (evening), Larson alone (`LARSON_ALONE`: the beach factor 1 everywhere, Cs 1.9e-5), the default
since: `nocape_rd198_dry_toeavg_la` (the final setup), `step2_module_la`, `step3_road_la` (the steps), and the Cs
sweep `la_cs10` ... `la_cs80` (shown, not used to choose Cs). Since `hat_run.yaml` now sends no beach widths, the
script gives the beach-factor variants the toe table's widths itself, so they rerun as before.
