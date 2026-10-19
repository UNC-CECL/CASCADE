# 2026-10-19 dune crest: DEM vs model year 0 and year 3

Asked by the advisor: is the dune threshold an elevation or a height above the berm, what are the
real dune crests along the island, and does CASCADE build the dunes up early (so too little overwash)?
Read-only on the current matrix runs (edgeBE, managed = road_bdm_groin / road_bdm_nourish_groin;
natural in the tables). Producer: scripts/hatteras_ms/experiments/HAT_dune_crest_dem_vs_model.py.

- The ceiling is given as an ELEVATION (Dmaxel 5.5 m NAVD88). Barrier3D turns it into a height above
  the berm: Dmax = Dmaxel - BermEl = 5.5 - 1.7 = 3.8 m. The dune files and DuneDomain are height above
  a fixed berm of 1.7 m NAVD88 (not a measured berm).
- Model year 0 equals the DEM dune files (max difference 2 mm / 6 mm).
- By year 3 the low cells have grown toward the ceiling and the tall ones have been cut down to it.
  P10 height above berm: 1996 2.23 -> 3.61 m; 2009 1.35 -> 3.01 m (managed; natural 1.77).
- Cells below the year's highest run-up, model vs DEM held fixed: 1999 (Dennis) 0.2% vs 6.6%;
  2012 6.8% vs 21.7% (2010 24% for both). So yes, the dunes outgrow the real low spots within 2-3 years.
