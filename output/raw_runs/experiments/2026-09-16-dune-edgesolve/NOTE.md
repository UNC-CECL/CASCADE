# 2026-09-16-dune-edgesolve

**Question.** What do the two locked end domains carry if they are solved
against the DUNE LINE instead of CoastSat? Every edgeBE pair in the config
closes the model's OLS shoreline rate on the CoastSat LRR target at GIS 1
and 90. Here the same Newton solve is run with the observation replaced by
the digitised dune line: end vintage minus start vintage per domain over
the survey interval, seaward positive, the quantity
`output/comparisons/duneline_vs_modeled_windows/` draws. Hannah, 2026-09-16,
by interview:

- all four windows;
- two readings of the target at an end domain, `raw` (the domain's own
  value) and `mean3` (the mean of it and its two inward neighbours,
  GIS 1-3 / 88-90), because at GIS 1 in 1996-2010 they differ by 5 m/yr;
- the ENDPOINT estimator on the model side (`change_rate_m_yr`), a dune
  line being two surveys; the CoastSat solves use the OLS slope;
- an experiment only: no preset, no config change. Both comparison figure
  sets are re-drawn with the dune-solved runs beside the CoastSat-solved
  ones.

**Brackets.** Each solve takes its first secant through two runs that
already impose 0 and the CoastSat-solved pair: for 1996-2010 and 2010-2024
the current matrix zeroBE and edgeBE road_bdm arms (09-15 and 09-16, the
code of today); for 1984-2004 and 2004-2024 the matrix pairs predate the
LRR refit and the 20 m relocation standard, so both are re-run here under
`brackets/`, on today's code and today's resolved versions (1984-start v2,
offsets 1984/v1; 2004-start v1, offsets 2004/v1).

**Layout.**

```
brackets/<window>/{zeroBE,edgeBE}/<run>/   fresh brackets for 1984 and 2004
raw/step<k>/<window>/edgeBE/<run>/         probes, target read raw
mean3/step<k>/<window>/edgeBE/<run>/       probes, target read as a 3-domain mean
logs/                                      one log per run
RESULTS.md                                 the solved pairs and the skill tables
```

Solve arithmetic: `HAT_be_edge_domain_solve.py --target duneline
--dune-smooth {raw,mean3} --estimator endpoint`, one `--kind` and `--tag`
per run.
