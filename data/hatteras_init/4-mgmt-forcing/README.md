# 4-mgmt-forcing - what people do to the island

NC-12 and the nourishment record: the forcing CASCADE's management modules
spend.

```
road_offset/      where the road sits, per domain, per period
road_elevation/   how high it is - ONE set, no year
road_relocation/  the measured displacement of each historical relocation
Hatteras_Management_Timelines.xlsx   the source record behind the nourishment
                                     projects in hatteras_site_config
```

**Road elevation carries no year on purpose.** There are two DEMs and under the
road they differ by a median 0.22 m, but that difference is an uncorrected
survey offset rather than a roadbed, so it is kept out of the forcing.

**A setback is metres landward of interior row 0**, so it belongs to the
extraction it was measured against. Spending one version's setbacks on another
version's arrays measures from a row that moved. Folders under
`road_offset/dunestart_offset/` are named for the PERIOD, and the derived ones
carry a `PROVENANCE.md` saying which survey is really behind them.
