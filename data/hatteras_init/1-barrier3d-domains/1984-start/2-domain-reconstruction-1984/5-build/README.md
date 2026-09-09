# 5-build — the version built from the footprint

Nothing is written here. `HAT_build_footprint_version.py` (scripts
`2-domain-reconstruction-1984/5-build/`) takes the rows of `../2-extent/`, the placement
of `../3-placement/` and the fill of `../4-fill/` and writes a dune-topo
version, with its own README, run manifest and footprint audit:

    ../../dune-topo/v3/     the behind-road placement, copy fill, 1984 setbacks (2026-09-08)

A version built with a different placement, once the imagery review has
decided one, goes in `../../dune-topo/v4/` and is recorded here the same way.

`HAT_plot_version_figures.py` (same scripts folder) gives a built version its
figure set, written into the version folder beside the extractor's for v1 and
v2: the per-domain grid (source beside version), the summary page and the plan
views. Run it after every build: `python HAT_plot_version_figures.py --version v3`.
