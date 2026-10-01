# storms - the groin study's storm files

```
HAT_resample_grointest_storms.py   real 1984-2004 storms resampled into 30 years -> ../../groin_init/storms/1967_1997/
HAT_build_1967_2017_storms.py      resample + two real periods stitched, 1967-2017 -> ../../groin_init/storms/1967_2017/
```

Run the resample first; the build reads it from `groin_init/storms/input_storms/`.
Map: `../../README.md`.

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning and choices behind it are here. Moved out of the scripts on 2026-10-01, when they were brought in line with `scripts/STYLE.md` (code unchanged, proven with `style_equivalence_check.py`; path fixes listed per script).

### HAT_resample_grointest_storms.py

The storm file for the 1967-1997 groin test, made by resampling the existing,
validated 1984-2004 series into 30 model years. No new storm physics.

**Why resample.** The groin acts through BRIE's wave-driven alongshore
transport, set by the wave climate parameters, not the storm file. Storms drive
Barrier3D's cross-shore overwash and dune response, very nearly orthogonal to
the groin test. So the test needs a record of the right length with realistic
event statistics, not a 1967-1997 storm climate. Storm timing is therefore not
historical: right for a function test, not for a forced hindcast.

**Modes.** `"bootstrap_years"` (used) copies all storms of a randomly chosen
source year into each target year, keeping within-year clustering and realistic
annual counts. `"bootstrap_events"` draws storms i.i.d. to a Poisson annual
count (`EVENTS_PER_YEAR`, default the source mean), breaking within-year
correlation. Seed 1967.

**Durations** are truncated at `MAX_STORM_DUR` (72 h), as the storm creator's
`max_storm_dur` did then, and storms under `MIN_STORM_DUR` (8 h) dropped.

**Format** (as `historical_storm_creation_v3_HAT.py`): columns time (1-based
model year), Rhigh and Rlow (dam MHW, (NAVD88 - MHW)/10), period (s), duration
(h). Saved as .npy (CASCADE input) and .csv (inspection).

**Path fix (2026-10-01).** `SAVE_DIR` was the dead literal
`/scripts/groin/HAT-buxton-hindcast-groin-test/groin_init/storms/1967_1997`; it now
resolves from the repo root to `groin_init/storms/1967_1997/`, where its
`1967_1997_grointest_storms` products are.

From the script's original header:

```text
HAT_resample_grointest_storms.py
============================
Build the storm file for the 1967-1997 groin TEST run by RESAMPLING an existing
CASCADE storm series into a new N-year window. No new storm physics -- this is a
reformat/retime of real events so the record covers the run length.

Why resample instead of create new storms
------------------------------------------
The groin acts through BRIE's wave-driven alongshore transport, which is set by
the WAVE CLIMATE parameters (wave_height, period, asymmetry), NOT by the storm
file. The storm file drives Barrier3D overwash/dune response (cross-shore), which
is very nearly orthogonal to the groin test. So we do not need a 1967-1997 storm
climate -- we need a record of the right LENGTH with realistic event statistics.
Resampling the real 1984-2004 events preserves their Rhigh/Rlow/period/duration
distributions and per-year storm counts while retiming them into N years.

Storm timing is therefore NOT historical for 1967-1997. State this in the run
docstring; it is correct for a function test, not a forced hindcast.

Format (matches historical_storm_creation_v3_HAT.py comparison)
-----------------------------------------------------------
Columns: time, Rhigh, Rlow, period, duration
  time     : 1-based model year (1 .. N)
  Rhigh    : max run-up  [decameters MHW]  (NAVD88 - MHW)/10, as in the creator
  Rlow     : min run-up  [decameters MHW]
  period   : wave period at peak TWL [s]
  duration : storm duration [hrs]
Saved as both .npy (CASCADE input) and .csv (inspection).
```

<details><summary>Function notes (the original docstrings)</summary>

**`_by_year()`**

```text
Split source storms into {source_year: rows[:,1:] (Rhigh,Rlow,period,dur)}.
```

</details>

<details><summary>Notes that were comments in the code</summary>

- Anchored 2026-09-14: absolute into a home directory, or into a tree renamed since. Rule 5 of ORGANIZATION.md.
- How to build each target year's storms: "bootstrap_years" : for each target year, copy ALL storms from a randomly chosen source year (preserves within-year clustering and realistic annual counts). Recommended. "bootstrap_events": draw individual storms i.i.d. to hit a target annual count (breaks within-year correlation; use only if you want a smoother, less clustered series).
- Duration cap (match your creator's max_storm_dur convention: TRUNCATE, do not drop). Set None to leave source durations untouched.

</details>

### HAT_build_1967_2017_storms.py

Stitches the 1967-2017 storm file for the extended groin-deterioration
hindcast from three sources, each on its own 1-based model-time index
(Barrier3D matches storms on `self._StormSeries[:, 0] == self._time_index`,
which starts at 1 at each run's first modelled year: not a calendar year, not
0-based).

| source | what | kept (own index) | shifted to |
|---|---|---|---|
| `1967_1997_grointest_storms.npy` | the 30-yr bootstrap resample above (7.03 storms/yr) | 1-17 (1967-1983) | 1-17 |
| `1984_2004_storms_v3_72.npy` | real storms: Stockdon (2006) runup on WIS waves + NOAA Duck tide gauge 8651370, 21 years | 1-21 (1984-2004) | 18-38 |
| `2004_2024_storms_v3_72.npy` | same method, 21 years | 2-14 (2005-2017) | 39-51 |

**The 2004 overlap.** Period 1's last year and Period 2's first are both
calendar 2004 and byte-for-byte identical (9 storms); `verify_no_duplicate_storms`
checks this before Period 2's copy is dropped. Result: one continuous index,
1967 = 1 through 2017 = 51 (17 resampled + 21 + 13 real years), checked for gaps
and extras before saving.

**END_YEAR convention.** The groin hindcast computes
RUN_YEARS = END_YEAR - START_YEAR, exclusive of END_YEAR (the original 1967-1997
run, END_YEAR = 1997, only simulated through 1996). To simulate all 51 years
here, through 2017, set END_YEAR = 2018, not 2017. The reminder the script
prints still names `HAT_groin_hindcast_1967_1997.py`; that script was deleted on
2026-10-01 (`git show <commit>^:hard-structures/groin/HAT-groin-buxton-output/1967_1997_run/HAT_groin_hindcast_1967_1997.py`,
commit from `git log --diff-filter=D --oneline -- <that path>`), and the same
convention holds in `../../../HAT-groin-buxton-output/1967_2017_run/`. The old
header's pointer to that script's "Section 2c" for the nourishment citation was
already dangling: it had no Section 2c.

**Path fix (2026-10-01).** `INPUT_DIR` and `OUTPUT_PATH` were dead literals under
`/scripts/groin/HAT-buxton-hindcast-groin-test/groin_init/storms/`; both now
resolve from the repo root to `groin_init/storms/input_storms/` (all three
sources are there) and `groin_init/storms/1967_2017/1967_2017_groin_storms.npy`.

From the script's original header:

```text
HAT_build_1967_2017_storms.py
==============================
Builds the combined 1967-2017 storm file for the extended groin-deterioration
mini hindcast, by stitching together three source files that each use their
own independent 1-based model-time-step index (confirmed against
barrier3d.py: storms are matched via `self._StormSeries[:, 0] == self._time_index`,
where `_time_index` starts at 1 at each run's own first modeled year -- NOT a
calendar year, and NOT 0-based).

SOURCES (see HAT_groin_hindcast_1967_1997.py Section 2c for the nourishment
citation note; storm sourcing is separate from that):
  1967_1997_grointest_storms.npy   - 30-yr bootstrap resample of the 1984-2004
                                      real storm distribution (7.03 storms/yr),
                                      built for the original 1967-1997 test.
                                      Only years 1967-1983 (its own time 1-17)
                                      are reused here -- real storms take over
                                      from 1984 onward, so the resampled years
                                      18-30 (1984-1996) in that file are no
                                      longer needed.
  1984_2004_storms_v3_72.npy       - real storms, Stockdon (2006) runup applied
                                      to WIS wave data + NOAA Duck tide gauge
                                      (station 8651370). Period 1, 21 model
                                      years (1984-2004 inclusive), own time 1-21.
  2004_2024_storms_v3_72.npy       - same methodology, Period 2, 21 model years
                                      (2004-2024 inclusive), own time 1-21.

THE 2004 OVERLAP: Period 1 (time=21) and Period 2 (time=1) both represent the
same real calendar year 2004, confirmed byte-for-byte identical (same 9 storms,
same Rhigh/Rlow/period/duration -- verified below in verify_no_duplicate_storms).
Period 1's copy is kept; Period 2 contributes only its 2005-2017 years (its own
time 2-14) to avoid double-counting 2004's storms.

RESULT: one continuous 1-based time index, 1967=1 through 2017=51 (51 model
years total: 17 resampled + 21 real Period-1 + 13 real Period-2).

IMPORTANT -- END_YEAR convention: HAT_groin_hindcast_1967_1997.py computes
RUN_YEARS = END_YEAR - START_YEAR, which is EXCLUSIVE of END_YEAR (the
original 1967-1997 run, with END_YEAR=1997, only ever simulated through 1996).
To actually simulate through 2017 inclusive (all 51 years in this file), set
END_YEAR = 2018 in that script's Section 3, not 2017.
```

<details><summary>Function notes (the original docstrings)</summary>

**`verify_no_duplicate_storms()`**

```text
Confirm 2004 (Period 1's time==21, Period 2's time==1) really is the
same real storm record before we rely on dropping one copy of it.
```

</details>

<details><summary>Notes that were comments in the code</summary>

- Shifts applied to each segment's time column so the combined file is one continuous 1..51 index with no gaps or overlaps.

</details>
