# `nc_shorelines.geojson` before the 1997 re-stamp — retired 2026-09-22

The NC Coastal Management shoreline file as it stood **before two 1997
features were re-dated**.

Both carried the placeholder stamp `1/1/1997`; they are the Hatteras-area
features, and were re-stamped to the flight date `9/27/1997`. The other 21
features in the file kept `1/1/1997` on purpose — they are different surveys
elsewhere in North Carolina, and asserting a Hatteras flight date statewide
would replace one wrong date with another. Which two to change was decided
spatially, by intersecting each 1997 feature with `cascade_area.geojson`.

The full reasoning, including the EPSG:3725 export artefact that had to be
read as UTM 18N for that intersection, is in the parent `PROVENANCE.md`.

**This copy is irreplaceable.** These geojsons are gitignored, so there is no
version of the pre-edit file in git history — this folder is the only record
of what the dates were before. Do not delete it to save space.

Rule 4 of `ORGANIZATION.md`: retirement is a dated folder with a note, and
deletes nothing.
