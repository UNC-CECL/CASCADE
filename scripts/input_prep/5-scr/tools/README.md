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
    Internal-consistency checks on the CoastSat rate outputs, before they are
    trusted as model targets:
      1. NaN audit        how many transects have no LRR, and why
      2. Domain mean      does the stored domain mean equal the mean of the
                          transect LRRs in the same CSV?
      3. Transect count   does n_transects match the actual row count?
      4. Per-domain       a side-by-side table, to spot a domain that looks off

    It reads finished products and writes nothing, so it is safe to run at
    any time.
```

Neither script is imported by anything. Both are run by path.
