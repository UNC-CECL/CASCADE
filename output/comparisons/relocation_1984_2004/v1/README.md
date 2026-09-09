# v1 — superseded

The six comparison sets (three source/sink presets, groin off and on) and the
dune-position check were made on 2026-09-01 from runs on **dune-topo v1**, the
original pick set of 2026-08-27. v1 was replaced as CURRENT by v2 (the re-pick
base) on 2026-09-04 and is no longer what a run loads by default. These are the
numbers `scripts/hatteras_ms/RELOCATION_COMPARISON_RESULTS.md` quotes and the
sets the 2026-08-27 archive was compared against; kept as that record.
Each `report.txt` names the two runs it read, with their version and commit.

**`calibBE_groin/` cannot be regenerated** (2026-09-09): its arm A,
`output/raw_runs/1984_2004/calibBE/HAT_1984_2004_calibBE_road_bdm_groin`, was
re-run on v2 on 2026-09-07 while its `reloc` twin stayed on v1, so the
calibration tree no longer holds a v1 pair for it and the script refuses a
cross-version pair. Its files are the 2026-09-01 render (no relocation tracker,
road drawn as a line). The other five sets were re-rendered 2026-09-09 with the
tracker and the road as CASCADE's two rows.
