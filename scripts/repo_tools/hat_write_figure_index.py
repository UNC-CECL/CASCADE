"""
hat_write_figure_index.py
==============================================================================
Write FIGURES.md at the repo root: one page that answers "which figure do I
open for X" across the three trees that hold shoreline-change work --
`data/hatteras_init/5-scr/3-rates/`, `.../4-comparisons/` and
`output/comparisons/`.

WHY THIS EXISTS (Hannah, 2026-09-21).  Each of those trees has a good README
of its own and none of them spans the others, so answering a question meant
already knowing which tree it lived in. Three things in particular were easy
to get wrong and are stated here once:

  * the estimator            rate (m/yr) vs distance (m), and which window
                             the rate was FITTED on -- see the vocabulary in
                             3-rates/README.md
  * the window               five window folders sit as peers across two
                             chains; see WINDOWS.md
  * the units                model_vs_observed/ is in m/yr; every other
                             comparison tree is in METRES. Two figures of the
                             same comparison in different units is the single
                             most confusable pair in the project.

HOW IT STAYS HONEST.  The map below is editorial -- a person decides which
question a figure answers -- but every path in it is CHECKED against the disk
when this runs, and a missing one is reported and marked in the output rather
than silently written. So the index cannot quietly rot the way a hand-kept
list does. Run it after adding or renaming a figure.

USAGE
    python scripts/repo_tools/hat_write_figure_index.py
    python scripts/repo_tools/hat_write_figure_index.py --check   # no write
==============================================================================
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())

SCR = "data/hatteras_init/5-scr"
RATES = f"{SCR}/3-rates"
COMP = f"{SCR}/4-comparisons"
OUT = "output/comparisons"

# question -> [(what you get, path, note)]. A path ending in "/" is a folder.
# <w> is substituted per window where a row covers several.
SECTIONS: list[tuple[str, str, list[tuple[str, str, str]]]] = [
    (
        "What did the shoreline do?",
        "CoastSat satellite waterline, observations only. No model anywhere in "
        "these.",
        [
            ("The rate, m/yr",
             f"{RATES}/coastsat/lrr/<w>/lrr_<w>.png",
             "OLS slope through every satellite position in the window. **This "
             "is the model's scoring target.**"),
            ("That rate as a distance, m",
             f"{RATES}/coastsat/total_change/<w>/total_change_<w>.png",
             "the window's own rate x its own years, against the observed "
             "change. TOTAL change — nothing extrapolated."),
            ("The long-term rate applied to a half, m",
             f"{RATES}/coastsat/projected/<w>/projected_<w>.png",
             "the 1996–2024 rate x 14 yr, on a window it was NOT fitted on. "
             "PROJECTED. 1996_2010 and 2010_2024 only."),
            ("Two snapshots differenced, m and m/yr",
             f"{RATES}/coastsat/endpoint/<w>/coastsat_endpoint_<w>.png",
             "mean position ±6 months about each dune-line date. No rate fit."),
            ("The rate in 5-year bins",
             f"{RATES}/coastsat/5yr_bins/<w>/lrr_5yr_bins_<w>.png",
             "is the trend steady inside the window?"),
            ("How much the alongshore smoothing changes it",
             f"{RATES}/coastsat/total_change/<w>/smoothed/",
             "and `projected/<w>/smoothed/`. Read the bias, not r — a smoother "
             "inflates r on both sides."),
        ],
    ),
    (
        "What did the dune line do?",
        "The digitized dune line, observations only.",
        [
            ("Net change between the two lines",
             f"{RATES}/duneline/endpoint/<w>/duneline_endpoint_<w>.png",
             "end line minus start line. Measured, never fitted."),
            ("Where each line actually sat",
             f"{COMP}/duneline_positions/",
             "maps, imagery zooms, distance to NC-12, beach width."),
        ],
    ),
    (
        "Did the dune line move with the shoreline?",
        "Both observations on one panel, the gap between them shaded as "
        "beach-width change. **Metres.**",
        [
            ("Shoreline as two snapshots",
             f"{COMP}/shoreline_vs_duneline/coastsat_endpoint_vs_duneline_endpoint/"
             "<w>/coastsat_endpoint_vs_duneline_<w>_alongshore.png",
             "observed vs observed; `..._scatter.png` beside it."),
            ("Shoreline as its OWN window's trend",
             f"{COMP}/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/"
             "<w>/coastsat_total_change_vs_duneline_<w>_two_panel.png",
             "also `_shaded_gap` and `_overlay`; the dune-interval value is a "
             "column in `domain_comparison.csv`."),
            ("Shoreline as the LONG-TERM trend carried onto a half",
             f"{COMP}/shoreline_vs_duneline/coastsat_projected_vs_duneline_endpoint/"
             "<w>/coastsat_projected_vs_duneline_<w>_two_panel.png",
             "the 1996–2024 LRR × 14 yr against the dune line measured over "
             "that half. 1996_2010 and 2010_2024 only."),
            ("One long-term prediction vs two dune-line outcomes",
             f"{COMP}/shoreline_vs_duneline/coastsat_projected_vs_duneline_endpoint/"
             "all_windows_stacked/"
             "coastsat_projected_vs_duneline_1996_2010_2024_halves_overlay.png",
             "the two halves stacked. The shoreline side is IDENTICAL in both "
             "panels, so every difference between them is the dune line's. Read "
             "beside the `coastsat_total_change_...` sheet of the same name."),
            ("The same two halves, each on its OWN rate",
             f"{COMP}/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/"
             "all_windows_stacked/"
             "coastsat_total_change_vs_duneline_1996_2010_2024_halves_overlay.png",
             "the counterpart of the row above. The pair separates what the "
             "long-term trend PREDICTS from what it was FITTED on."),
            ("Does smoothing change any of it?",
             f"{COMP}/shoreline_vs_duneline/smoothed_loess7/",
             "both sheets again with BOTH curves LOESS-smoothed at 7 domains "
             "(3.5 km). Read the beach width, not r — a symmetric smoother "
             "inflates r on both sides."),
            ("The whole period above its two halves",
             f"{COMP}/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/"
             "all_windows_stacked/"
             "coastsat_total_change_vs_duneline_1996_2010_2024_stacked.png",
             "three panels, shoreline and dune line as lines with the gap "
             "shaded. A different question from the two-panel sheets above."),
            ("Did they change pace together?",
             f"{COMP}/shoreline_vs_duneline/coastsat_total_change_vs_duneline_endpoint/"
             "change_between_periods/"
             "coastsat_total_change_vs_duneline_change_between_periods.png",
             "second half minus first half, on both sides. Only for the "
             "total product: the projected one uses the same rate in both "
             "halves, so its difference is zero by construction."),
        ],
    ),
    (
        "How does the model compare to the observations?",
        "**These are in m/yr, not metres** — the one tree that is. Sliced by "
        "which observation, not by estimator; the cross-reference below says "
        "which estimator each folder uses.",
        [
            ("vs the CoastSat shoreline",
             f"{OUT}/model_vs_observed/vs_shoreline/domain_means/"
             "model_vs_shoreline_means_<w>.png",
             "`_grid.png` puts all four windows on one sheet."),
            ("vs the shoreline, as graded",
             f"{OUT}/model_vs_observed/vs_shoreline/smoothed/"
             "model_vs_shoreline_smoothed_<w>.png",
             "the LOESS form the runner actually scores."),
            ("vs the dune line",
             f"{OUT}/model_vs_observed/vs_duneline/endpoint_net_change/"
             "model_vs_duneline_netchange_<w>.png",
             "and `net_change_smoothed/` beside it."),
            ("vs both at once",
             f"{OUT}/model_vs_observed/vs_shoreline_and_duneline/"
             "model_vs_shoreline_and_duneline_<w>.png",
             "each solve drawn in its own target's estimator."),
            ("Does the answer survive changing the estimator or the solve?",
             f"{OUT}/model_vs_observed/sensitivity/",
             "`ends-swapped`, `dune-raw-solve`, `mixed-estimator`. The arm is "
             "in every filename."),
        ],
    ),
    (
        "Which target should the model be graded on?",
        "CoastSat against the dune line, with the runs, as **net change in "
        "metres** over each 14-yr window.",
        [
            ("Start here",
             f"{OUT}/target_comparison/projected/paired/"
             "target_and_own_run_projected_<w>.png",
             "each target with its own run. **`projected/` is the target in "
             "use**: the 1996–2024 LRR x 14 yr."),
            ("The same on each window's own rate",
             f"{OUT}/target_comparison/total_change/paired/"
             "target_and_own_run_total_change_<w>.png",
             "kept for the record — what the runner grades against."),
            ("Neither target fitted anywhere",
             f"{OUT}/target_comparison/projected/ends_unsolved/",
             "zeroBE: no source/sink term in any domain, so all 90 are the "
             "model's own response."),
            ("Does the grading window matter?",
             f"{OUT}/target_comparison/smoothing_scale/"
             "projected_vs_model_<w>.png",
             "**No** — and r was never the number to read. See its "
             "PROVENANCE.md."),
            ("Both targets AND the model, smoothed",
             f"{OUT}/target_comparison/smoothed_loess7_with_cascade/",
             "the two smoothed sheets with the zeroBE run over them in dark "
             "green. Nothing in that run was fitted to either target, so all "
             "90 domains are the model's own response."),
            ("The numbers",
             f"{OUT}/target_comparison/projected/tables/skill.csv",
             "bias, RMSE and r per window, model set and target."),
        ],
    ),
]

# Item 4 of the 2026-09-21 tidy: model_vs_observed slices by OBSERVATION while
# every other tree slices by ESTIMATOR. Rather than rename it -- its axis
# genuinely is a different question -- say which estimator each folder uses.
CROSSREF = [
    ("vs_shoreline/domain_means", "CoastSat", "OLS rate (`lrr_m_yr`)",
     "raw domain means"),
    ("vs_shoreline/smoothed", "CoastSat", "OLS rate (`lrr_m_yr`)",
     "spliced LOESS — the form the runner grades"),
    ("vs_duneline/endpoint_net_change", "dune line",
     "endpoint rate (`change_rate_m_yr`)", "raw domain means"),
    ("vs_duneline/net_change_smoothed", "dune line",
     "endpoint rate (`change_rate_m_yr`)", "spliced LOESS"),
    ("vs_shoreline_and_duneline", "both", "endpoint rate, both solves",
     "raw domain means"),
    ("sensitivity/mixed-estimator", "dune line",
     "**OLS** rate against an **endpoint** observation",
     "deliberately mismatched, as a sensitivity"),
]


def resolve(path: str) -> list[tuple[str, bool]]:
    """Expand <w> over the windows that exist, and check each path."""
    if "<w>" not in path:
        return [(path, (REPO / path).exists())]
    out = []
    for s, e in [(1996, 2010), (2010, 2024), (1996, 2024), (1984, 2004), (2004, 2024)]:
        p = path.replace("<w>", f"{s}_{e}")
        if (REPO / p).exists():
            out.append((p, True))
    return out or [(path, False)]


def build() -> tuple[str, list[str]]:
    missing: list[str] = []
    L = [
        "# Which figure do I open?",
        "",
        f"*Written {dt.datetime.now():%Y-%m-%d} by "
        "`scripts/repo_tools/hat_write_figure_index.py`, which checks every "
        "path below against the disk. Re-run it after adding or renaming a "
        "figure.*",
        "",
        "Three things to fix before reading any of these:",
        "",
        "1. **Which window.** Five window folders sit as peers across two "
        f"chains and one is context only — [`{SCR}/WINDOWS.md`]({SCR}/WINDOWS.md).",
        "2. **Which estimator.** A rate turned into a distance is named by the "
        "window it was *fitted* on: **total** = same window, **projected** = "
        f"carried onto another, **observed** = no rate — [`{RATES}/README.md`]"
        f"({RATES}/README.md).",
        "3. **Which units.** `output/comparisons/model_vs_observed/` is in "
        "**m/yr**. Every other comparison tree is in **metres**. The same "
        "comparison exists in both, and they are not interchangeable.",
        "",
        "Every figure carries its quantity, window and method in its title, "
        "and a full caption in `supporting/CAPTIONS.md` beside it.",
        "",
    ]
    for title, blurb, rows in SECTIONS:
        L += [f"## {title}", "", blurb, "", "| for | open | note |", "|---|---|---|"]
        for what, path, note in rows:
            hits = resolve(path)
            ok = all(o for _, o in hits)
            if not ok:
                missing.append(path)
            shown = path if len(hits) > 1 or "<w>" in path else hits[0][0]
            windows = ""
            if "<w>" in path:
                windows = " <br>*windows:* " + ", ".join(
                    p.split("/")[-2] if p.split("/")[-1].endswith("/") else
                    [w for w in ("1996_2010", "2010_2024", "1996_2024",
                                 "1984_2004", "2004_2024") if w in p][0]
                    for p, _ in hits)
            mark = "" if ok else " **← MISSING**"
            L.append(f"| {what} | `{shown}`{mark}{windows} | {note} |")
        L.append("")
    L += [
        "## Cross-reference: what `model_vs_observed/` actually plots",
        "",
        "That tree is sliced by which **observation** the model is held "
        "against, while `3-rates/` and `4-comparisons/` are sliced by "
        "**estimator**. Same underlying observations, different question, so "
        "the folder names do not line up. This is the translation:",
        "",
        "| folder | observation | model estimator | reading |",
        "|---|---|---|---|",
    ]
    for folder, obs_, est, reading in CROSSREF:
        L.append(f"| `{folder}` | {obs_} | {est} | {reading} |")
    L += [
        "",
        "All of it in **m/yr**. To see the same comparison as a distance in "
        "metres, use `output/comparisons/target_comparison/`.",
        "",
        "## Where the rest lives",
        "",
        "| | |",
        "|---|---|",
        "| Finished figures for the paper | `output/figures/<subject>/` |",
        "| How the repo is laid out | [`ORGANIZATION.md`](ORGANIZATION.md) |",
        "| Figure house style | `scripts/site_layer/hat_figure_style.py`, "
        "`figure_making/STYLE.md`, `output/figures/style/` |",
        "| Model runs | `output/raw_runs/`, indexed in `run_index.csv` |",
        "",
    ]
    return "\n".join(L), missing


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--check", action="store_true",
                    help="report missing paths without writing the file")
    a = ap.parse_args(argv)
    text, missing = build()
    for m in missing:
        print(f"MISSING  {m}")
    if a.check:
        print(f"{len(missing)} missing path(s)")
        return 1 if missing else 0
    out = REPO / "FIGURES.md"
    out.write_text(text, encoding="utf-8")
    print(f"wrote {out.name}  ({len(missing)} missing path(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
