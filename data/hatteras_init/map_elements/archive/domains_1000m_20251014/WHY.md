# Retired 2026-09-18: the 1000 m domain polygons

`HAT_domains.shp` (+ sidecars, dated 2025-10-14): 92 boxes, 1000 x 500 m,
EPSG:26918, columns `ID`, `LabelYN`, `LRR_1978_2`.

They are not the model's domains. Those are 90 boxes, 2000 x 500 m, in
`5-scr/2-transect-frame/transect_domains/HAT_domains.json`, and this set's IDs
run about nine domains south of them. That came to light on 2026-09-17, when
this set's box for "45" did not contain the elevation array for domain 45.
Nothing reads these any more; they are kept because figures drawn before
2026-09-17 may have used them.
