# duneline_windows — net dune-line change, a long window and its halves

Written by `scripts/input_prep/5-scr/duneline_windows/duneline_windows.py`
(2026-09-18). It is the dune-line counterpart of
`../coastsat_windows/1996_2024/`, drawn with the same code so the two sit
side by side.

```
1996_2024/
    duneline_change_1997_2023_halves.png   (a) 1997→2023 filled by sign;
                                           (b) 1997→2009 grey, 2009→2023 black
    supporting/
        duneline_change_1997_2023_halves.csv   per domain: the three line
                                               positions and the three changes
        duneline_change_1997_2023_halves.pdf, CAPTIONS.md
```

The figure is also published to `output/figures/shoreline/`.

## What the number is

**Net position change in metres. Not a rate, and not a fit.** For each line,
the per-transect station is its distance from the fixed offshore datum
(`2-brie-offset/raw_offsets/<vintage>_duneline_offset_raw.csv`, read the way
the hindcast loader reads it). The stations are averaged over the ~5
transects of each 500 m domain. The change is the start position minus the
end position, so **seaward is positive**.

- **No survey dates are involved**, and the 2023 imagery date is unknown.
- **The halves add up to the whole exactly**, because they share the 2009
  line. The script asserts this. The whole is always the sum of its parts,
  unlike the CoastSat LRR, where a jump at the join can put the long-term
  rate outside both halves.

## Which lines

A period year reads its line through `hat_topo_version.DUNE_LINE_FOR_YEAR`:
1996 → 1997, 2010 → 2009, 2024 → 2023. The figure is labelled with the real
line years, and the caption says what they stand in for. All three lines
were re-digitized on 2026-09-18; see `2-brie-offset/dunelines/README.md`.

Results at build time, over the 90 domains:

| change | island mean | domains landward | range |
|---|---|---|---|
| 1997→2023 | −14.8 m | 61 / 90 | −90 to +55 m |
| 1997→2009 | −16.2 m | 73 / 90 | −69 to +29 m |
| 2009→2023 | +1.4 m | 37 / 90 | −45 to +38 m |

## Rebuilding

```
python scripts/input_prep/5-scr/duneline_windows/duneline_windows.py
python scripts/input_prep/5-scr/duneline_windows/duneline_windows.py --years 1996 2010 2024
```
