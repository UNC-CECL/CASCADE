# Retired 2026-09-18

* **plot_overwash.py** -- the first overwash plot of a hindcast run, read from
  an absolute path into a run folder (`HAT_1984_2004_basestorms_Hs2p0`) that
  no longer exists. Its NPZ loader lives on, copied, in
  `../compare_overwash_figures.py` and `../compare_overwash_observed.py`, which
  resolve their run through `cascade_pipeline.run_registry`.
* **overwash_pea_early.py** (was `file_from_roya/overwash_pea.py`) -- an
  earlier copy of Roya's Pea Island overwash script. Pea Island is a different
  site; the maintained copy, which reads BermEl and MHW from the YAML, is
  `scripts/other_ms/pea_island_ms/overwash_pea.py`. Renamed on retirement so
  the two cannot be mistaken for each other.

Kept, not deleted: rule 4 of ORGANIZATION.md.
