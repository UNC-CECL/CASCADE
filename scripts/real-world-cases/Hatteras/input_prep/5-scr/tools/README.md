# tools — checks and indexes

Neither of these makes a product. They check one, or describe one. That is why
they are here rather than under a numbered stage: filing them beside a product
would suggest they are part of building it.

```
windows_index.py
    Regenerates data/hatteras_init/5-scr/WINDOWS.md from
    hat_observed_rates.WINDOW_ROLE, so the table cannot drift from the code.

    Five window folders sit as peers under most 3-rates and 4-comparisons
    products — 1984_2004, 1996_2010, 1996_2024, 2004_2024, 2010_2024 — and
    they are not equivalent: two chains, and one context window nothing is
    graded against. A folder named only by its years says none of that, so
    1984_2004 reads as a current result. This renders the roles, plus the
    coverage grid showing which product actually has which window.

        python scripts/input_prep/5-scr/tools/windows_index.py

coastsat_rates_check.py
    Internal-consistency checks on one CoastSat rate window, before it is
    trusted as a model target:
      1. NaN audit        how many transects have no LRR, and why
      2. Domain mean      does the stored domain mean equal the mean of the
                          transect LRRs in the same CSV?
      3. Transect count   does n_transects match the actual row count?
      4. Per-domain       a side-by-side table, to spot a domain that looks off

        python tools/coastsat_rates_check.py                          # 2010-2024
        python tools/coastsat_rates_check.py --start-year 1996 --end-year 2010

    Reads finished products and writes ONE file, the per-domain comparison:
    3-rates/coastsat/lrr/<window>/checks/rates_check_<window>.csv. A `checks/`
    subfolder, so a verification artefact is never mistaken for a product.
```

### It was broken until 2026-09-22, in three ways

Found by running it rather than reading it -- an earlier version of this
README described it from its docstring and said it "reads finished products
and writes nothing, so it is safe to run at any time". All three were wrong.

1. **It died on a cp1252 console** before finishing the first check, on the
   U+2713 it prints as a pass mark. It now carries the same stdout guard its
   siblings do (`export_be_calibration.py`), reconfiguring to UTF-8 rather
   than ASCII-ifying, so the next mark someone types cannot reintroduce it.
2. **Its output path was a dead literal**, drive-rooted into
   `scripts/input_preperation/CoastSat_verification/` -- note the old spelling;
   that tree has not existed since the rename. It resolves through
   `hat_observed_rates` now, and creates the folder it writes to.
3. **Its window was hardcoded to 2004_2024**, which is not on the canonical
   1996 -> 2010 -> 2024 chain. The window is now an argument, defaulting to
   2010-2024.

Verified on both current windows: all four checks pass on each.

Neither script is imported by anything. Both are run by path.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### coastsat_rates_check.py

Check the CoastSat rate outputs are internally consistent before they reach CASCADE.

From the script's original header:

```text
CoastSat Pipeline Data Verification
Runs targeted checks on the outputs of the CoastSat → CASCADE pipeline
to confirm that transect LRR values, domain averages, and transect counts
are internally consistent before these values enter CASCADE.

Checks performed
  1. NaN audit        — how many transects have missing LRR, and why
  2. Domain mean match — do the pre-computed domain means equal the manual mean
                         of transect LRRs in the same CSV?
  3. Transect count   — does n_transects in the summary match the actual count
                         in the transect file?
  4. Per-domain detail — prints a side-by-side table for every domain so you
                         can spot any domain that looks off

Usage
Edit the CONFIG section to point at your files, then run:
    python coastsat_rates_check.py

Outputs
  Console report (always)
  verification_report.csv  — detailed per-domain comparison table
```

Notes that were in the code:

```text
The window to check. It was hardcoded to 2004_2024, which is not on the
canonical 1996 -> 2010 -> 2024 chain, so the default now follows the chain
and either end can be overridden.
```

```text
Where to save the per-domain comparison table. Resolved through
hat_observed_rates, beside the window it checks, in a `checks/` subfolder so
a verification artefact is never mistaken for a product. Until 2026-09-22
this was a drive-rooted literal naming scripts/input_preperation/ -- a tree
renamed long ago -- so the script could not finish even when it ran.
```

```text
Tolerance for floating-point comparison in Check 2 (m/yr)
Differences smaller than this are treated as matching
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_never_die_on_a_print()`**

```text
Stop a console encoding from killing a finished check.

This file prints U+2713 and U+2717 as pass/fail marks, which a Windows
cp1252 console cannot encode, so `print` raises UnicodeEncodeError -- it
died on the FIRST check, before reporting anything. Reconfigure rather
than ASCII-ify, so the next mark someone types cannot reintroduce it.
The same guard its siblings carry (export_be_calibration.py).
```

**`check_nan_lrr()`**

```text
How many transects are missing LRR values?
A large number here usually means transect IDs didn't match CSV filenames.
```

**`check_domain_means()`**

```text
Does the mean_lrr in domain_lrr_summary.csv equal the manual mean
of lrr_m_yr values for each domain in transect_lrr_full.csv?

A mismatch means something changed between when the summary was computed
and the current transect file, or a different set of transects was used.
```

**`check_transect_counts()`**

```text
Does n_transects in domain_lrr_summary.csv match the actual count
of rows per domain in transect_lrr_full.csv?

Includes both valid and NaN-LRR transects (n_transects = total in domain).
```

**`build_detail_table()`**

```text
Build a side-by-side per-domain table combining:
  - Values from domain_lrr_summary.csv
  - Manually computed values from transect_lrr_full.csv
```

</details>

### windows_index.py

Write data/hatteras_init/5-scr/WINDOWS.md from the window definitions in the code.

From the script's original header:

```text
Write data/hatteras_init/5-scr/WINDOWS.md from the definitions in
site_layer.hat_observed_rates, so the table cannot drift from the code.

WHY THIS EXISTS.  Five window folders sit as peers under most 3-rates and
4-comparisons products -- 1984_2004, 1996_2010, 1996_2024, 2004_2024,
2010_2024 -- and they are NOT equivalent: two chains, and one context window
nothing is graded against. A folder named only by its years says none of
that, so 1984_2004 reads as a current result (Hannah, 2026-09-21). The roles
live in `hat_observed_rates.WINDOW_ROLE`; this script renders them, plus the
coverage grid showing which product actually has which window, which is the
other question a reader cannot answer by looking.

USAGE
    python scripts/input_prep/5-scr/tools/windows_index.py
```
