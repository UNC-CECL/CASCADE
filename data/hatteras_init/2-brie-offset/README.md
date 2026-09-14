# 2-brie-offset - where each domain starts, cross-shore

The BRIE island offset: one distance per domain, setting where that domain's
barrier sits at model year zero.

```
raw_offsets/      one CSV per dune-line survey, per transect, from GIS
hindcast_<year>/  the CASCADE-facing files for one period start
```

Each `hindcast_<year>/` holds the padded 120-domain file the model reads, the
unpadded 90, and a buffer diagnostic figure.

## Two things that bite

**The padded files are each zeroed on their own most seaward domain**, so
differencing two of them subtracts a constant and flips the sign of the mean.
To compare two years, difference the RAW files: they share a fixed offshore
datum.

**A file named for a period is not always a survey of that year.** The 1996
file is built from the 1997 line, the nearest island-wide survey. See
`raw_offsets/PROVENANCE.md`, which also records that the 2017 file is a
Buxton-only clip and cannot stand in for anything island-wide.

Built by `scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year`.
