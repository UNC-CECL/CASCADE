# 2026-09-19-edgesolve-2010

**Question.** The 2010-2024 end-domain solve (`../2026-09-16-edgesolve-2010/`,
now in `archive/2026-09-19-pre-redigitized-sens-exp/`) re-asked after the 2009
dune line was re-digitized: do the locked end values in
`HATTERAS_BE_EDGE_ONLY[2010]` (+72.6 at GIS 1, +31.3 at GIS 90) still close the
CoastSat misfit on offsets 2010/v1?

**Method.** `be_dune_edgesolve_loop.py --target coastsat --windows 2010`
(the `--target coastsat` mode was added for this): the same Newton solve as
09-16 (model lrr_m_yr against target_lrr_m_yr, GIS 1 raw, GIS 90 LOESS-10),
bracketed by the re-run 09-18 matrix zeroBE and edgeBE full-management runs.
full_management, nourishment on, no groin, relocations off, Hs 2.5.

**Answer.** Yes. At step 0 the matrix edgeBE run already sat at -0.018 / +0.014
m/yr residual; one probe (+72.8 / +31.2) closed it to +0.004 / +0.002
(`loop_log.csv`). The 0.1-0.2 m/yr differences are within noise, so the config
values stand unchanged (no re-run of the matrix needed).

**Layout.** `coastsat/step1/2010_2024/edgeBE/<run>/`, `logs/`, `loop_log.csv`.
