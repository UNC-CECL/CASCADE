# diagnose_road_drowning.py -- retired 2026-09-18

Drew the INITIAL cross-shore profile of chosen domains with the road and the
cells RoadwayManager's drowning check inspects at t=0, to explain why roads
drowned on the first time step.

Retired rather than repaired (Hannah, 2026-09-18):

* It reads topography from `topography/2009/2009_v2/domain_<N>_topography_2009.npy`,
  a layout that went period-first on 2026-08-25 (arrays now resolve through
  `hat_topo_version.topo_dirs(product)` / `array_path()`, with no year in the
  name), and the buffer from `hatteras_init/buffer/` rather than
  `1-barrier3d-domains/buffer/`.
* Its road state came from `output/raw_runs/HAT_1984_2004_base/`, a run
  folder from before the purpose-based raw_runs layout.
* The problem it diagnosed -- roads drowning at t=0 in GIS 12, 14, 15, 50,
  51, 86 -- went away with the gap-filled DEM (v5): 3 a year to 0.

To look at a t=0 profile again, start from `topo_dirs(product)` for the
arrays and `road_setback_file(year)` for the road, not from this.

Kept, not deleted: rule 4 of ORGANIZATION.md.
