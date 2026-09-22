# 2026-09-16-dune-edgesolve - findings

Eight solves, all converged (residual under 0.02 m/yr at both ends, mostly
under 0.01), 24 runs including the four fresh brackets. Numbers in
RESULTS.md; the runs each solve stands on are in solved.csv. Figures:
`output/comparisons/observed_vs_modeled_windows/dune_ends_<raw|mean3>/` (the
CoastSat observation) and
`output/comparisons/duneline_vs_modeled_windows/sensitivity/dune-solved-ends/<raw|mean3>/` (the
dune-line observation), same layout as their CoastSat-solved parents.

## 1. The ends move a lot; the sign at GIS 1 is not robust

Against the dune line the end values are a third to a half of the CoastSat
ones, and GIS 1 changes sign in three of four windows:

| window | CoastSat | dune raw | dune mean3 |
|---|---|---|---|
| 1984-2004 | -42.6 / +13.0 | +31.8 / +9.5 | -1.6 / +6.2 |
| 1996-2010 | +32.2 / +10.0 | -2.7 / +2.1 | -28.6 / -3.1 |
| 2004-2024 | +50.3 / +46.7 | -2.7 / +18.1 | -6.2 / +15.9 |
| 2010-2024 | +72.6 / +31.3 | +27.8 / +12.3 | +27.2 / +10.8 |

GIS 90 is stable across the two readings (within 3 m/yr) and always smaller
than CoastSat's, by a factor of 1.5 to 3. GIS 1 is not: in 1984-2004 the
raw reading asks for +31.8 and the three-domain mean for -1.6, in 1996-2010
for -2.7 against -28.6. That is one domain's worth of dune line (GIS 1 alone
is +4.05 m/yr in 1984-2004 while GIS 2-3 are near zero; GIS 2 is -10.5 m/yr
in 1996-2010) deciding the sign of a boundary term, because the gain is
about a tenth and every metre a year of target costs ten of imposed rate.
The CoastSat target at GIS 1 is also a raw domain mean, but of ten
satellite transects through ~250 dates; the dune line is five transects
through two.

## 2. The interior hardly notices

Moving the ends by 30 to 70 m/yr moves the interior score by 0.2 to 0.6 m/yr
in bias and up to 0.6 in RMSE. The footprint of an end value is the ten or
so domains BRIE diffuses it over, and the solve was scored on GIS 2-89.

Against CoastSat every dune-solved run is worse than the CoastSat-solved
one it started from, as it must be: RMSE rises 0.2 to 0.6 in three windows
and is flat in 1996-2010 raw (1.11 vs 1.14). Against the dune line the
dune-solved runs are better in RMSE in every window (by 0.2 to 0.6) but not
uniformly in bias: 2004-2024 and 2010-2024 go from near zero to -0.6 and
-0.7, because closing the ends on a dune line that barely moved pulls a
model that already erodes too fast further landward across the reach.

## 3. What this says

- The end values are boundary-artefact absorbers, and the target they
  absorb towards is a choice. The two observations disagree with each other
  at the ends by more than the model disagrees with either in the interior,
  so the choice sets the ends without settling the interior.
- The dune line is the noisier target for this purpose. A single-domain
  reading at GIS 1 flips the sign of the solve between windows and between
  readings; a three-domain mean tames 1984 and destabilises 1996. There is
  no smoothing that is right for both, which is the same conclusion the
  CoastSat protocol reached when it kept D1-D10 raw.
- Keeping CoastSat as the end target is the defensible default: it is the
  observation the interior is scored on, its LRR is the same estimator the
  model side uses, and its transect density makes a raw end-domain value
  mean something. The dune-solved pairs are a sensitivity, not a
  replacement; if a dune-line-based end is ever wanted, mean3 at GIS 90
  and a reading wider than one domain at GIS 1 is the place to start.

## Provenance

Brackets for 1984 and 2004 re-run here on today's code (1984-start v2,
offsets 1984/v1; 2004-start v1, offsets 2004/v1) because the matrix pairs
predate the LRR refit and the 20 m relocation standard; their CoastSat
residuals reproduce the matrix rows. 1996 and 2010 use the matrix pairs of
09-15 and 09-16. Every probe ran full_management, relocations off, groin
off, Hs 2.5, edgeBE with HAT_BE_OVERRIDE at both ends. Solve script:
`be_edge_domain_solve.py --target duneline --dune-smooth {raw,mean3}
--estimator endpoint`; results: `be_dune_edgesolve_results.py`.
