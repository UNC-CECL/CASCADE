# 4-mgmt-forcing - what people do to the island

NC-12 and the nourishment record: the forcing CASCADE's management modules
spend.

```
road_offset/      where the road sits at a period's START, per domain
road_elevation/   how high it is - ONE set, no year
road_relocation/  the measured displacement between the two digitised lines
nourishment/      FIGURES ONLY: when and where the fills were placed, the
                  site-config projects against the record, and the m^3/m each
                  domain receives (scripts/input_prep/4-mgmt-forcings/
                  beach_nourishment.py). Not an input; the model reads the
                  list in hatteras_site_config
Hatteras_Management_Timelines.xlsx   the source record behind the nourishment
                                     projects in hatteras_site_config
```

## Do not type these paths

Every path in this tree resolves through `scripts/site_layer/hat_topo_version.py`:
the per-period helpers (`road_setback_file`, `road_line_file`,
`road_mask_file`, `road_relocation_file`, `legacy_setback_file`) and the roots
beside them (`ROAD_ELEVATION_FILE`, `MGMT_RECORD_XLSX`, `NOURISHMENT_DIR`, ...).
About thirty scripts typed them until 2026-09-18, and four were reading a
folder that had been renamed a week earlier without knowing it.

## Three organising rules, not one (settled 2026-09-15)

A hindcast can start in 1984, 1996, 2004 or 2010. Not everything here is per
start year, and it pays to be explicit about which is which:

* **State at year zero is per START YEAR.** The road setback is where NC-12
  sits relative to interior row 0 when the run begins, so every start year
  owns one file, `road_offset/dunestart_offset/{measured,derived}/<year>/`.
  Two are measurements (1984, 2004); two are built from those (1996 = 1984
  plus the 1989 Pea Island relocation; 2010 = 2004 unchanged).
* **History is ONE timeline, and the run window selects from it.** The
  relocation events and the nourishment projects are single lists in
  `hatteras_site_config.py` (`HATTERAS_ROAD_EVENTS`,
  `HATTERAS_NOURISHMENT_PROJECTS`); a run fires whatever falls between its
  start and end year. Nothing in this tree is a per-period copy of that record,
  and none should be added.
* **Elevation has no year at all.** One file, for every period; see below.

So a new start year needs a setback file plus entries in
`scripts/site_layer/hat_topo_version.py` (`ROAD_LINE_FOR_YEAR`, `ROAD_SETBACK_KIND`,
`YEAR_PRODUCT`) and `hatteras_site_config.py` (`HATTERAS_PERIODS`), and
nothing else here.

## Two axes, so two kinds of integer

The digitised NC-12 lines are **1978 and 2008 exports that stand in for the
1984 and 2004 starts**. Until 2026-09-15 they were filed under the start
years, so `raw_offset/1984` meant "the 1978 line" while
`dunestart_offset/1984` meant "the 1984 start". They are now filed under the
line's true vintage:

```
road_offset/raw_offset/1978/, 2008/     the lines (nc12_1978, nc12_2008)
road_offset/raster/1978/, 2008/         their masks on the domain grids
road_relocation/1978_2008/              the displacement between them
```

The one place a start year is paired with a line is
`scripts/site_layer/hat_topo_version.py:ROAD_LINE_FOR_YEAR` -- 1984 and 1996 read the
1978 line, 2004 and 2010 the 2008 line -- and the helpers beside it
(`road_line_file`, `road_mask_file`, `road_setback_file`) refuse a start year
where a line vintage is expected. A folder named 1984 or 2004 in this tree is
therefore always a PERIOD; 1978 and 2008 are always LINES.

**Road elevation carries no year on purpose.** There are two DEMs and under the
road they differ by a median 0.22 m, but that difference is an uncorrected
survey offset rather than a roadbed, so it is kept out of the forcing.

**A setback is metres landward of interior row 0**, so it belongs to the
extraction it was measured against. Spending one version's setbacks on another
version's arrays measures from a row that moved. Every folder under
`road_offset/dunestart_offset/` carries a `PROVENANCE.md` saying which line,
which extraction, and -- for `derived/` -- which measured file is behind it.
