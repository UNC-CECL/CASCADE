# 3-rates/duneline/endpoint/2004_2024 - provenance

Written 2026-09-18 17:37 by scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py.

Net change between the two dune lines that bound the window, per 100 m transect, seaward positive; per domain the mean over its transects. `change_m` needs no date; `rate_m_yr` is `change_m` over the survey interval.

| line | raw stations | survey date |
|---|---|---|
| 2004 | `2004_duneline_offset_raw.csv` (written 2026-09-15 17:25) | 2004-05-25 |
| 2023 | `2023_duneline_offset_raw.csv` (written 2026-09-18 15:02) | 2023-07-01 **ASSUMED** |

Interval 19.10 yr (2004 and 2024 read the 2004 and 2023 lines through DUNE_LINE_FOR_YEAR). The 2023 flight date is not known and is centred on 1 July; `rate_m_yr` inherits that assumption, `change_m` does not.

## Island summary

450 transects in 90 domains. Domain mean change -1.7 m (-0.09 m/yr); 42 of 90 domains landward; range -77.7 to +61.0 m.
