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
