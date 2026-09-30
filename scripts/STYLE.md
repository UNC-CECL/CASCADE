# How a script is written

A script should read top to bottom in a minute. It says **what** it does; the
**why** lives in the folder's `README.md`. Where files go is
[`ORGANIZATION.md`](../ORGANIZATION.md); how figures look is
[`figure_making/STYLE.md`](figure_making/STYLE.md), and how a figure script is
put together is [`figure_making/GUIDE.md`](figure_making/GUIDE.md).

The reference implementation is `input_prep/5-scr/template/`: copy its shape.

This applies to new scripts and to any script being rewritten. Older scripts
with long explanatory docstrings are not wrong; bring them in line when you
next rewrite them, not in a sweep.

---

## 1. Layout, in this order

```python
"""
Step 2: shoreline change rate (OLS slope, LRR) per transect and per zone.

    python shoreline_rates_template.py --start-year 1996 --end-year 2024

Window = 1 Jan start year through 31 Dec end year. Input: one CSV per
transect; the lookup from step 1. Needs pandas, numpy, scipy, matplotlib.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
...

# --- CONFIG ------------------------------------------------------------------
START_YEAR, END_YEAR = 1996, 2024
MIN_OBS      = 10
MAX_P_VALUE  = 1.0       # off; screening on p drops stable transects (README)
# -----------------------------------------------------------------------------


# Read one transect CSV: date + position, by column order
def load_timeseries(path: Path) -> pd.DataFrame:
    ...


# Run: window, fit, screen, average by zone, write
def main() -> None:
    ...
    # Fit every transect inside the calendar-year window
    ...
    # Write both tables
    ...


if __name__ == "__main__":
    main()
```

1. **Docstring**: one line saying what the script does (with `Step N:` when
   order matters), the command to run it, one short paragraph of inputs and
   dependencies, then the author block.
2. **Imports**, standard library first.
3. **CONFIG**: every value a user might change, between the two ruled lines.
   Nothing below it should need editing.
4. **Functions**, in the order `main()` calls them.
5. **`main()`**, then the `if __name__ == "__main__":` guard.

## 2. Comments: one short line, and only that

- **One line above each function** saying what it does: `# OLS slope of
  position vs time, with 95% CI`.
- **One line above each block in `main()`** naming the step: `# Write both
  tables`.
- **An end-of-line comment** only where a value could be misread: a unit,
  a switch that is off on purpose, a pointer to the README.
- **No banners** except the CONFIG rules. No paragraphs of reasoning, and no
  history ("changed on 09-14 because..."). The reasoning goes in the folder's
  README, and the history goes in git.

If a comment needs more than one line, it belongs in the README.

## 3. Run order is in the folder names

When scripts must run in sequence, number their folders: `1-zone-join/`,
`2-shoreline-change/`. Scripts that can run in either order share a number.
Folders use `N-name`; file names do not carry numbers (see ORGANIZATION.md
rule 2 and each tree's naming section).

## 4. The author block

Every script Hannah writes ends its docstring with:

```
Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: YYYY-MM-DD
```

- **Version** is the date the script last changed in substance. Bump it when
  you change what the script does or produces, not for a typo.
- A script whose header is a `#` comment block carries the same four lines as
  `# ` comments.
- A notebook carries it under the title in the first markdown cell.
- **Never on someone else's work.** Code from colleagues or upstream CASCADE
  (`other_ms/`, `from_lexi/`, `from_roya/`, `colleague_old_version/`,
  `cascade/`) keeps its original authorship.
- **A script adapted from someone else's** gets the block with one line above
  it crediting the source, by path and author:

  ```
  Adapted from: from_lexi/historical_storm_creation_v3.ipynb, by Lexi
  Author:  Hannah A. Henry, Coastal Environmental Change Lab,
  ...
  ```

  Borrowing an idea is not adapting: a script that only took inspiration says
  so in its README, not in this line.

## 5. Habits

**Paths come from the repo, never typed.** Find the root by searching upward,
and take data locations from the `site_layer/` resolvers
(`hat_observed_rates`, `hat_topo_version`, `hat_env_forcings`, ...). See
ORGANIZATION.md rules 5 and 6.

```python
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
```

**Units are in the name.** Every variable and column carries its unit:
`_m`, `_m_yr`, `_dam`, `_yr`. State the sign convention once (seaward
positive). What each quantity is and where it is converted:
[`UNITS.md`](../UNITS.md).

**Fail loudly, never guess.** A missing input, CRS, column or window stops the
script with a message naming what *is* available. No silent fallbacks, no
default that quietly stands in for real data.

**Report, don't correct.** When two things disagree, write the mismatch out as
its own column or note; don't adjust either side to fit. For example, the
endpoint template reports `interval_obs_yr` beside the calendar interval
instead of swapping one for the other.

**Arguments for what varies, CONFIG for what doesn't.** A window, period or
version that changes run to run is a command-line argument. Constants live in
CONFIG. Runner settings come from the `HAT_` environment variables
(`HAT_hindcast_config`).

**Outputs say where they came from.** The window or version is in the output
filename (`rates_per_zone_1996_2010.csv`). A product folder carries a
`PROVENANCE.md`: the script that wrote it, when, and from which inputs.

**Print UTF-8 on Windows.** A script that prints anything beyond ASCII
(✓, →, ±) reconfigures its output first, or it crashes on the cp1252 console:

```python
for stream in (sys.stdout, sys.stderr):
    stream.reconfigure(encoding="utf-8", errors="replace")
```

**Experiments don't edit main code.** An experiment drives the unchanged runner
(swapping inputs in-process or through a worktree on `PYTHONPATH`) and answers
one question. Adopting its answer is a separate, deliberate change.

**Check against a known answer.** A script that makes a model target or a
shared template gets a test with an answer you already know, such as the
synthetic coast in `input_prep/5-scr/template/README.md`. Run it after every
edit.
