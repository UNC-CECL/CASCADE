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
    ** DOES NOT RUN AS OF 2026-09-22. ** See below.

    Internal-consistency checks on the CoastSat rate outputs, before they are
    trusted as model targets:
      1. NaN audit        how many transects have no LRR, and why
      2. Domain mean      does the stored domain mean equal the mean of the
                          transect LRRs in the same CSV?
      3. Transect count   does n_transects match the actual row count?
      4. Per-domain       a side-by-side table, to spot a domain that looks off
```

### coastsat_rates_check.py is broken, in two ways

Found 2026-09-22 by running it. An earlier version of this README said it read
finished products, wrote nothing and was safe to run at any time. All three
were wrong.

1. **It dies on a cp1252 console** before finishing the first check. It prints
   a U+2713 check mark, and a Windows console cannot encode that:
   `UnicodeEncodeError: 'charmap' codec can't encode character '✓'`. The
   fix is the one its siblings already carry -- reconfigure stdout to UTF-8 at
   import, as `export_be_calibration.py` and the 6-scr-smooth scripts do.

2. **It writes, to a tree deleted long ago.** `OUTPUT_REPORT_CSV` is a
   drive-rooted literal naming `scripts/input_preperation/CoastSat_verification/`
   -- note the old spelling; that folder has not existed since the tree was
   renamed. `detail.to_csv(OUTPUT_REPORT_CSV)` at the end would fail even if
   the encoding did not stop it first. It is one of the 20 unresolvable paths
   `hat_layout_check.py` reports under rule 5.

Neither is hard to fix, and neither has been, because nothing runs this. The
checks it performs are worth having; if you want it back, the work is the
stdout reconfigure plus routing the report through
`site_layer/hat_observed_rates.py` instead of a typed path.

Neither script is imported by anything. Both are run by path.
