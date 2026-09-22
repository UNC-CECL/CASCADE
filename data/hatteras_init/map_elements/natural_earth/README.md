# natural_earth - the locator-map coastline

`ne_10m_states_southeast_us.geojson`: Natural Earth 1:10m Admin 1 states and
provinces (with lakes), public domain, https://www.naturalearthdata.com/,
downloaded 2026-09-17 and clipped to 84.5-73 W, 31-40 N: the thirteen states
the regional inset of the site figures can show. Columns `name`, `postal`.
Read by `regional_inset()` in scripts/figure_making/island/study_area_figures.py.
Kept as a geojson because the shapefile family is gitignored and a locator
map that cannot be redrawn from the repository is not reproducible.
