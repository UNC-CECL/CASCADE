# lrr_smoothed — the lrr field at two LOWESS widths

The test of 5 against 7 domains (2.5 km against 3.5 km) of alongshore smoothing,
on the candidate windows 1996–2015 and 2010–2026 and their reference 1996–2026
(Hannah, 2026-10-02). 7 domains is the group's range and the width every run is
graded at; 5 is the narrower alternative.

```
<window>/
    lrr_lowess5_vs_lowess7_<window>.png      transects, raw domain means, both curves
    supporting/
        lrr_lowess5_vs_lowess7_<window>.pdf
        lrr_lowess5_vs_lowess7_<window>.csv  domain_number, mean_lrr, lowess5, lowess7, lowess7_minus_lowess5, unc95_lowess7, ref_lowess7_1996_2026, ref_unc95_lowess7_1996_2026
```

The rates are the window's own lrr tables (`../lrr/<window>/transect_lrr_full.csv`),
smoothed with the same `spliced_lowess_series` the scoring target uses, so GIS 1–10
keep their raw means and the two curves are identical there by construction. Every
figure shares the y axis of the windows' lrr figures.

On the two candidate windows the 1996–2026 rate is drawn as a dashed grey line,
smoothed at 7 domains so it is compared like with like (`--reference`). North of
GIS 10 the 7-domain curve sits 0.49 m/yr below it for 1996–2015 (RMS 0.83) and
0.89 m/yr above it for 2010–2026 (RMS 1.22).

The 7-domain curve and the reference carry a 95% band: each transect's fit
half-width (`unc_m_yr`), smoothed the same way as the rate. These are OLS
intervals, which assume independent residuals; a shoreline record is serially
correlated (seasons, storms), so the bands are too narrow and "bands apart"
overstates how distinguishable two curves are.

`2010_2026/lrr_lowess7_2010_2020_and_2010_2026_vs_1996_2026.png` is the 2021-step
check: the same start with the record cut at 2020. North of GIS 10, 2010–2020 sits
0.25 m/yr below the long-term rate and 2010–2026 0.89 above it; 2010–2026 sits
1.14 m/yr above 2010–2020 (RMS 1.61). Most of 2010–2026's excess over the
long-term rate comes from 2021–2025.

Every figure of the candidate windows (here and the lrr window figures) shares one
y axis, ±12 m/yr, taken over 1996–2015, 2010–2026, 1996–2026 and 2010–2020
(`coastsat_lrr_windows.candidate_half`).

The CoastSat record ends 13 January 2026, so a window ending 2026 holds two weeks
of that year.

Producer: `scripts/input_prep/5-scr/3-rates/coastsat/lrr_smoothed/coastsat_lrr_smoothed.py`
(`--windows`, `--widths`).

| window | sd of lowess7 − lowess5 (GIS 11–90) | largest gap | domains > 0.5 m/yr |
|---|---|---|---|
| 1996_2015 | 0.265 m/yr | 1.04 m/yr at GIS 25 | 6 of 80 |
| 2010_2026 | 0.189 m/yr | 0.59 m/yr at GIS 66 | 2 of 80 |
| 1996_2026 | 0.143 m/yr | 0.53 m/yr at GIS 29 | 1 of 80 |
