# 3-rates/duneline/endpoint/2010_2024 - provenance

Written 2026-09-18 17:37 by scripts/input_prep/5-scr/duneline_endpoint/duneline_endpoint.py.

Net change between the two dune lines that bound the window, per 100 m transect, seaward positive; per domain the mean over its transects. `change_m` needs no date; `rate_m_yr` is `change_m` over the survey interval.

| line | raw stations | survey date |
|---|---|---|
| 2009 | `2009_duneline_offset_raw.csv` (written 2026-09-18 15:01) | 2009-05-30 |
| 2023 | `2023_duneline_offset_raw.csv` (written 2026-09-18 15:02) | 2023-07-01 **ASSUMED** |

Interval 14.09 yr (2010 and 2024 read the 2009 and 2023 lines through DUNE_LINE_FOR_YEAR). The 2023 flight date is not known and is centred on 1 July; `rate_m_yr` inherits that assumption, `change_m` does not.

## Island summary

450 transects in 90 domains. Domain mean change +1.4 m (+0.10 m/yr); 37 of 90 domains landward; range -44.7 to +37.6 m.
