"""
write_windows_md.py
==============================================================================
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
    python scripts/input_prep/5-scr/tools/write_windows_md.py
==============================================================================
"""

from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))

from site_layer.hat_observed_rates import (  # noqa: E402
    COMPARISONS, CURRENT_CHAIN, LEGACY_CHAIN, RATES, WINDOW_NOTE, WINDOW_ROLE,
    is_current_chain,
)

WINDOWS = [(1996, 2010), (2010, 2024), (1996, 2024), (1984, 2004), (2004, 2024)]
# The products that are filed one folder per window. Label -> path.
PRODUCTS = {
    "3-rates/coastsat/lrr": RATES / "coastsat" / "lrr",
    "3-rates/coastsat/endpoint": RATES / "coastsat" / "endpoint",
    "3-rates/coastsat/total_change": RATES / "coastsat" / "total_change",
    "3-rates/coastsat/projected": RATES / "coastsat" / "projected",
    "3-rates/coastsat/5yr_bins": RATES / "coastsat" / "5yr_bins",
    "3-rates/duneline/endpoint": RATES / "duneline" / "endpoint",
    "4-comp/coastsat_endpoint_vs_duneline_endpoint":
        COMPARISONS / "shoreline_vs_duneline" / "coastsat_endpoint_vs_duneline_endpoint",
    "4-comp/coastsat_total_change_vs_duneline_endpoint":
        COMPARISONS / "shoreline_vs_duneline" / "coastsat_total_change_vs_duneline_endpoint",
    "4-comp/coastsat_projected_vs_duneline_endpoint":
        COMPARISONS / "shoreline_vs_duneline" / "coastsat_projected_vs_duneline_endpoint",
}


def coverage_grid() -> list[str]:
    """Which product has which window, read off the disk, not asserted."""
    head = "| product | " + " | ".join(f"{s}–{e}" for s, e in WINDOWS) + " |"
    rows = [head, "|" + "---|" * (len(WINDOWS) + 1)]
    for label, root in PRODUCTS.items():
        cells = ["yes" if (root / f"{s}_{e}").is_dir() else "—" for s, e in WINDOWS]
        rows.append(f"| `{label}` | " + " | ".join(cells) + " |")
    return rows


def main() -> int:
    out = _REPO / "data" / "hatteras_init" / "5-scr" / "WINDOWS.md"
    cur = " → ".join(str(y) for y in CURRENT_CHAIN)
    leg = " → ".join(str(y) for y in LEGACY_CHAIN)
    lines = [
        "# Which window is which",
        "",
        f"*Written {dt.datetime.now():%Y-%m-%d} by "
        "`scripts/input_prep/5-scr/tools/write_windows_md.py` from "
        "`site_layer.hat_observed_rates.WINDOW_ROLE`. Edit the dict, not this "
        "file.*",
        "",
        "Most products under `3-rates/` and `4-comparisons/` are filed one "
        "folder per window, and the folders are named only by their years. "
        "They are **not** peers.",
        "",
        "| window | role | |",
        "|---|---|---|",
    ]
    for w in WINDOWS:
        mark = "**current chain**" if is_current_chain(w) else "legacy chain"
        lines.append(f"| `{w[0]}_{w[1]}` | {WINDOW_ROLE[w]} | {mark} — "
                     f"{WINDOW_NOTE[w]} |")
    lines += [
        "",
        f"The **current chain** is {cur}: the model is fitted on the first "
        f"half and the second is held out. `1996_2024` spans both and exists "
        "for context only — no run is graded against it, and a number read "
        "from it is not a model result.",
        "",
        f"The **legacy chain** is {leg}. Both of its periods are still live in "
        "`hatteras_site_config.HATTERAS_PERIODS` and runs on them still exist, "
        "but it was superseded as the main chain in September 2026. A figure "
        "in a `1984_2004/` or `2004_2024/` folder is not the current answer to "
        "anything unless you meant to ask about that chain.",
        "",
        "## Which product covers which window",
        "",
        "Read off the disk, so it is what is actually there:",
        "",
        *coverage_grid(),
        "",
        "The gaps are deliberate, not missing work:",
        "",
        "- `projected/` has no `1996_2024`: there the rate window IS the "
        "change window, so the answer is `total_change/1996_2024` and building "
        "it twice under two names is the confusion the 2026-09-21 rename "
        "removed. See `3-rates/README.md` for the total / projected / observed "
        "vocabulary.",
        "- `projected/`, `total_change/` and `5yr_bins/` cover the current "
        "chain only. They were built after it became the main chain. The two "
        "`*_projected_vs_duneline_endpoint` windows are the two halves: over "
        "the full period the rate window IS the change window, so that case "
        "is the `total_change` product.",
        "- `lrr/`, both `endpoint/` products and "
        "`coastsat_endpoint_vs_duneline_endpoint/` cover all five, because "
        "they predate the switch.",
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out.relative_to(_REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
