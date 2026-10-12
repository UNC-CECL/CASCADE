# target-end-window

Options for the start and end lines of each period's target: calibration 1996–2009 and test 2009–2025. Each run ends on 1 Jan of its end year. The current end lines are CoastSat means centred on a DEM date (or its anniversary), so they fall about 7 months after the model's end. For calibration, that means the line includes storms the model never sees. Each study tries other end lines against the same start and scores the existing runs. No new runs.

Until 2026-10-09 this topic was `calibration-end-window/`.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-10-07-end-window-options`](2026-10-07-end-window-options/README.md) | How much does each target move if its end line is centred elsewhere? | **Calibration:** every alternative is 3.1–3.6 m seaward (LOWESS-7), with the same shape (r ≥ 0.988); the whole year 2009 has the current line's mean date and is still +3.5 m, so the gap is the late-2009/2010 storms. **Test:** the alternatives move the target 1–3 m, against a 10–18 m model shortfall; the whole year 2025 equals the current line. | diagnostic; the calendar 2008 ±6-month calibration end was adopted 2026-10-09 |
| [`2026-10-09-6mo-starts`](2026-10-09-6mo-starts/README.md) | Does a ±6-month start shoreline (offset and target start) change the runs, against ±1 yr? | Barely: model net change moves RMS 0.8–1.6 m, bias about 1 m, RMSE under 0.5 m. The ±6-month lines need a 5-position minimum to cover the island. | diagnostic; ±1 yr kept |
| [`2026-10-10-test-end-calendar`](2026-10-10-test-end-calendar/README.md) | Which calendar-year window should replace the DEM-anniversary test end? | Calendar 2025 ±6 months (calibration's rule, off-centre because CoastSat ends 2026-01-13) or calendar year 2025 only: either moves bias < 0.3 m and RMSE about 1 m. | diagnostic; choice pending |
| [`2026-10-10-end-options-calendar`](2026-10-10-end-options-calendar/README.md) | The selected end (last modelled year ±6 months) against four calendar-based ends, both periods, no DEM-dated ends | All within about 1.5 m bias and 2 m RMSE. The selected window is the only one that falls the same way in both periods (about 63% of its data in the last modelled year, mean date 7–8 months before the model end). | diagnostic; replaces the 10-07 option set |

Back to [the map](../README.md).
