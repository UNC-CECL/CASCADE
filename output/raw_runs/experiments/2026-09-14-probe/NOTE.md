# 2026-09-14-probe

**Question.** Probes around the 2026-09-14 change that rounds a prescribed
relocation displacement to whole 10 m cells (see the note in
`hatteras_site_config.py` near `round_to_cell`), all on 1984-start v1:
`natural` and `roadonly` (zeroBE, no relocations), `v1setbacks` (roadway-only
on the v1-era setback CSV), `paired` (calibBE with relocations, rounded) and
`unrounded` (the same without rounding: RMSE 0.523102 against 0.522873).

**Answer lives in.** The recode comparison beside this folder
(`2026-09-14-recode`) and `tools/HAT_compare_rerun.py`.

**Runs deletable?** Yes; each is a single probe whose numbers are in the index.
