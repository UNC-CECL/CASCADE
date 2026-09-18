# 3-rates/duneline/endpoint/1996_2024 - provenance

Written 2026-09-18 17:37 by scripts/input_prep/5-scr/duneline_endpoint/duneline_endpoint.py.

Net change between the two dune lines that bound the window, per 100 m transect, seaward positive; per domain the mean over its transects. `change_m` needs no date; `rate_m_yr` is `change_m` over the survey interval.

| line | raw stations | survey date |
|---|---|---|
| 1997 | `1997_duneline_offset_raw.csv` (written 2026-09-18 15:01) | 1997-10-12 |
| 2023 | `2023_duneline_offset_raw.csv` (written 2026-09-18 15:02) | 2023-07-01 **ASSUMED** |

Interval 25.72 yr (1996 and 2024 read the 1997 and 2023 lines through DUNE_LINE_FOR_YEAR). The 2023 flight date is not known and is centred on 1 July; `rate_m_yr` inherits that assumption, `change_m` does not.

## Island summary

450 transects in 90 domains. Domain mean change -14.8 m (-0.58 m/yr); 61 of 90 domains landward; range -90.4 to +54.6 m.
