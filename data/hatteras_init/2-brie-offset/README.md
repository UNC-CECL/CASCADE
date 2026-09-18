# 2-brie-offset - where each domain starts, cross-shore

The BRIE island offset: one distance per domain, setting where that domain's
barrier sits at model year zero.

```
dunelines/        the digitised dune lines, duneline_<vintage>[_v<n>].geojson,
                  named for the IMAGERY year (here since 2026-09-15; they were
                  under 1-barrier3d-domains/raw-duneline-geojson/)
transects/        the 100 m transect layer every offset is measured along
raw_offsets/      one CSV per dune-line VINTAGE, per transect; a period finds
                  its vintage through hat_topo_version.DUNE_LINE_FOR_YEAR
<year>/           one period start: a PROVENANCE.md index, a CURRENT file,
<year>/v<n>/      and one build per version -- every start is versioned since
                  2026-09-15 (1984 v1, 1996 v2, 2004 v1, 2010 v1 from the
                  2009 line)
<year>/superseded_*/   the flat build a start had before it was versioned
```

Each build holds the padded 120-domain file the model reads, the unpadded 90,
a buffer diagnostic figure (and, since 2026-09-16, `_buffer_diagnostic_v2.png`: the
same padded profile at 1:1 scale, drawn by `HAT_plot_offset_profile_1to1.py --all`,
because the original stretches 7.5 km of offset across 72 km alongshore and
reads as a deep V; the slope is the coast's bearing, not its curvature), a copy
of the raw file it was built from, the
comparison with the previous build, and a PROVENANCE.md written by the
driver. `hatteras_site_config._island_offset_file` resolves the path: env
`HAT_OFFSET_VERSION_<year>` first, then CURRENT, then the only version
present. Run metadata records `island_offset_version`, so a run made on a
flat build (before 2026-09-15) shows it as unset.

## Do not type these paths

Resolve them through `scripts/site_layer/hat_topo_version.py`:
`offset_file(year, "padded" | "input" | "unpadded")` for a start's CURRENT build
(`offset_version()` holds the choice: `HAT_OFFSET_VERSION_<year>`, then
`CURRENT`, then the only `v<n>`; `hatteras_site_config` calls it too),
`dune_raw_file_for_year()` / `duneline_geojson()` for the lines, and
`RAW_OFFSET_DIR`, `DUNELINE_DIR`, `TRANSECT_FILE_100M` for the rest. Until
2026-09-18 about fifteen scripts typed these, and four read layouts that no
longer existed: the flat per-start build, `hindcast_<year>/` folders, and a
`dunelines/` under `1-barrier3d-domains/`. Two more took "the first sorted
match" across every build, which for 1984 and 2004 is the superseded flat one.

## Two things that bite

**The padded files are each zeroed on their own most seaward domain**, so
differencing two of them subtracts a constant and flips the sign of the mean.
To compare two years, difference the RAW files: they share a fixed offshore
datum.

**A period is not always served by a line of its own year.** The 1996 start
reads the 1997 line, the nearest island-wide survey; the pairing is
`DUNE_LINE_FOR_YEAR` in `scripts/site_layer/hat_topo_version.py` and nowhere else. See
`raw_offsets/PROVENANCE.md`, which also records that the 2017 file is a
Buxton-only clip and cannot stand in for anything island-wide.

## From a dune line to a model input

One command (2026-09-15):

```
python scripts/input_prep/2-brie-offset/1-produce/build_island_offset.py \
    --duneline duneline_1984.geojson --year 1984          # -> 1984/v<next>/, CURRENT
```

which runs, in order, the three scripts beside it:

```
duneline_to_raw_offsets.py --duneline <geojson> --out <raw csv>   line x transects -> raw_offsets/<vintage>_...
island_offset_hybrid.py --year <y> --version vN [--raw-file ...]   raw -> <year>/vN/
HAT_compare_offset_versions.py --year <y> --a <prev> --b vN        two builds side by side
```

then writes `CURRENT` and the version's `PROVENANCE.md` from the run. The
driver refuses a geojson whose vintage the table does not pair with `--year`,
and never overwrites a version. Until 2026-09-15 the first step was an ArcGIS
session; `raw_offsets/PROVENANCE.md` records which files came from which,
and the 1 m convention that separated them.
