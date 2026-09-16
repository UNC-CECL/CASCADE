# 2026-09-02-pea1989

**Question.** The control for the GIS 84-86 seaward-row insert (the Pea Island
1989 relocation): the 1984-2004 calibBE full-management groin-on run on
1984-start **v1** with the setback CSV as shipped, once with the prescribed
relocations (`base`) and once without (`basenoreloc`). Every insert arm it was
the control for was deleted on 2026-09-07 (unmodified topography only).

**Made by.** `scripts/hatteras_ms/experiments/HAT_run_crest_experiment.py`.

**Answer lives in.** `output/experiments/pea1989_crest/` (frozen) and
`data/hatteras_init/1-barrier3d-domains/LINEAGE.md`.

**Runs deletable?** Yes once the crest figures are no longer needed; they are
on v1 and reproducible from the script. Model state was kept because the
comparison reads roadway objects.
