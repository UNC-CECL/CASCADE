# map_elements - the shared map layers

Not a stage, so not numbered: the layers any figure or map draws on,
whichever stage it belongs to. Resolve them through
`scripts/site_layer/hat_map_layers.py`; do not type the paths.

```
hatteras_outline/     HAT_island_outline.shp (+ sidecars). Hatteras Island
                      only: the NC 1:80k land polygons it overlaps. Checked
                      2026-09-18 against the source: 88.85 km2 either way,
                      zero difference.                          ISLAND_OUTLINE
nc_coast_80k/         nc_80k_hatteras_window.geojson. NC 1:80k coastline
                      (D:/Hatteras_GIS/Outlines/nc_80k/, EPSG:4326), every
                      polygon within 25 km of the domain boxes, clipped: the
                      sound shores, Ocracoke and Pea Island that a map window
                      reaches past the island. Rebuilt by
                      scripts/figure_making/tools/clip_nc_coast.py.  NC_COAST
natural_earth/        ne_10m_states_southeast_us.geojson, the locator-map
                      states; see its README.                   NE_STATES
archive/
    domains_1000m_20251014/   the retired domain polygons; see WHY.md
```

## The domain boxes are not here, on purpose

The model's domains are
`5-scr/2-transect-frame/transect_domains/HAT_domains.json` (90 boxes,
2000 x 500 m, EPSG:3725). `hat_map_layers.DOMAIN_BOXES` points there, so a
figure script needs one import and there is still only one copy.
`D:/Hatteras_GIS/domains.geojson` is the same file (geometry and every
attribute identical, checked 2026-09-18); the figure and analysis scripts
that read it off the drive read the repository copy now. The 0-elevation
producers still name the drive, because they need it for the DEMs anyway.

## Moved here 2026-09-18

These sat in `data/hatteras_init/9-figures/map_elements/`. `9-figures` was
split by kind: the map layers are data and came here; the written house
style is documentation for whoever writes a figure script and went to
`scripts/figure_making/STYLE.md`; the rendered style sheet is a figure and
went to `output/figures/style/`. Both of those are written by
`write_style_sheet()` in `scripts/site_layer/hat_figure_style.py`.

## Not in git

The shapefile family (`*.shp`, `*.dbf`, `*.shx`, `*.prj`, ...) is ignored
repository-wide, and the geojsons here are untracked, so a clone has none of
these layers. `nc_coast_80k/` can be rebuilt from the D: drive; the outline
and the Natural Earth clip cannot be rebuilt from the repository.
