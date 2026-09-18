#!/usr/bin/env python3
"""Export the converged source/sink field to data/hatteras_init/7-source-sink.

WHY THIS EXISTS
    `hatteras_site_config.py` is the source of truth -- it is what the runs
    import -- but it is a Python module in the scripts tree, which is a poor
    place to look for "what were the calibrated values". The data directory
    already held per-period copies, and they had drifted: written 2026-06-15,
    they carry GIS 1 = -40 against the current -41.8, zeros across D2-D11 that
    are now +1.4 to +2.6, and the 2004 file was truncated mid-dict. The config
    itself carries a comment warning readers not to trust them.

    So this regenerates them FROM the config rather than alongside it, and adds
    the context a bare dict cannot carry: which domains were eligible to be
    corrected at all, which were withheld and why, and how much of each final
    value came from the one-shot solve versus the iteration.

WHAT IT WRITES (paths from site_layer/hat_source_sink.py since 2026-09-18)
    4-export/be_rates_<period>.py        the dict, same shape as the files
                                         it replaces
    4-export/be_calibration_domains.csv  one row per domain: zone,
                                         eligibility, the pass-0 and final
                                         rate, what the iteration added,
                                         and the residual still standing
    README.md                            at the top of 7-source-sink/:
                                         provenance, and the caveats that
                                         matter

    It READS the default pair's 2-calibrate/1984_2004__2004_2024/ and
    3-figures/1984_2004__2004_2024/. The figures are written straight there
    by stage 3, so the data directory is self-contained -- someone handed
    just this folder can see what was done, not only what came out.

    convergence_history.json is NOT copied up. It lives once, in the pair's
    2-calibrate/ folder, beside the calibration that wrote it -- a second
    byte-identical copy at the top level gave one fact two owners, and the
    two would drift the first time a pass was re-run without re-exporting.

    The superseded 2026-06-15 files are MOVED to archive/superseded_<date>/
    rather than deleted -- they are what earlier runs were built against, so
    they are history, not clutter.

Usage:
    python scripts/input_prep/7-source-sink/4-export/HAT_export_be_calibration.py [--check]

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import re
import shutil
import sys


def _never_die_on_a_print():
    """Stop a console encoding from killing a finished export.

    This file prints en- and em-dashes, which a Windows cp1252 console cannot
    encode, so `print` raises UnicodeEncodeError. Its sibling
    HAT_be_zone_residual_fit.py lost a completed calibration to exactly that
    on 2026-08-28 -- the crash landed between computing the numbers and
    writing them out. Reconfigure rather than ASCII-ify, so the next dash
    someone types cannot reintroduce it.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            try:
                stream.reconfigure(errors="replace")
            except Exception:
                pass


_never_die_on_a_print()
from datetime import date

import pandas as pd

_HERE = pathlib.Path(__file__).resolve()
PROJECT_BASE_DIR = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_BASE_DIR / "scripts"))

from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402

from site_layer import hat_source_sink as _be  # noqa: E402

# Resolved by hat_source_sink.py since 2026-09-18. The README stays at the
# top of 7-source-sink/; the exported field goes to 4-export/, the step that
# writes it; the calibration and figures read are the DEFAULT pair's folders.
DATA_DIR = _be.BE_ROOT
EXPORT_DIR = _be.EXPORT_DIR
FIGURE_DIR = _be.figures_dir()
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
# The calibration products moved into the data tree 2026-09-12, so this
# reads them from there rather than from beside the script that made them.
CALIB_OUT = _be.calibrate_dir()
CONFIG = PROJECT_BASE_DIR / "scripts" / "site_layer" / "hatteras_site_config.py"

# WHERE THE FIGURES COME FROM (corrected 2026-09-12). This used to read them
# out of the calibration OUTPUT directory and copy them into 3-figures/ --
# but HAT_be_zone_residual_fit.py writes its figures straight there, and
# the copies in the output directory were older. So the export overwrote the
# CURRENT figures with SUPERSEDED ones, quietly, every time it ran.
#
# 3-figures/ is the record now, and the staleness check below reads the same
# place. The superseded copies were moved under superseded_20260825/.
FIGURE_SOURCE_IS_THE_RECORD = True

# Paths RELATIVE TO the pair's figure folder (FIGURE_DIR,
# 3-figures/1984_2004__2004_2024/), which is why each carries its subfolder.
# The figures were grouped on 2026-09-14 by the question each answers:
# 1-field is what the calibration produced and what it was fitted against,
# 2-method is the evidence that the way it was produced holds up, and
# 3-limits is what it deliberately does not do. A flat folder gave those
# three equal weight, and the limits figure is the one most often read as a
# calibration failure rather than a stated boundary.
FIGURES = (
    ("2-method/fig_be_zones_and_corrections.png",
     "which domains qualified, and the correction each received"),
    ("2-method/fig_be_convergence.png",
     "the iteration sequence and the frozen zone set"),
    ("1-field/fig_be_diagnostic.png",
     "observed vs modelled rate, and the residual, per period"),
    ("1-field/fig_be_rates.png",
     "the BE field, hindcast and forecast scenarios"),
    ("3-limits/fig_groin_reserved_residual.png",
     "why the D5-D7 residual is left uncorrected"),
)

# The masked-iteration lineage. The unmasked attempt (backups 214732, 220039)
# is deliberately NOT here: it was abandoned for correcting outside the
# geomorphological zone set, and mixing its values into the record would
# misrepresent what was shipped.
# The one-shot solve, for the be_pass0_* / iteration_added_* split. The apply
# step backs up BEFORE it writes, so the pass-0 field is the backup taken before
# PASS 1 -- the file before pass 0 holds the superseded field, not this lineage.
# Re-pointed 2026-09-14 with the corrected iteration, whose backups were kept;
# the previous target was the retired lineage's pass-0 file and was never
# committed, which is why these two columns exported empty.
# This pointed at scripts/, where the backup has not been since it was filed
# under 2-calibrate/prebe/, so the next export would have written the pass-0
# columns empty. One definition now, shared with HAT_plot_be_zones.py
# (2026-09-18).
PASS0_BACKUP = _be.PASS0_BACKUP

_ROW = re.compile(r"^\s*(\d+):\s*([+-]?\d+\.?\d*),\s*(?:#\s*(.*))?$", re.M)
NL = chr(10)


def rates_from(path, period):
    """{gis: (rate, label)} for one period out of a site-config file."""
    text = pathlib.Path(path).read_text(encoding="utf-8")
    block = text.split("HATTERAS_BE_RATES_CALIBRATED")[1]
    segment = block.split(f"{period}:")[1]
    segment = segment[:segment.find("},")]
    return {int(m.group(1)): (float(m.group(2)), (m.group(3) or "").strip())
            for m in _ROW.finditer(segment)}


def stale_against_runs(paths):
    """Which of `paths` predate the newest calibBE run output.

    The exporter COPIES whatever is on disk; it does not regenerate. So a
    figure left over from before the last hindcast would be published into the
    data directory looking exactly as authoritative as a current one. This is
    the check that catches that -- it is the failure mode that actually
    happened on 2026-08-25, when fig_be_convergence.png was three minutes older
    than the run that had just been rebuilt.
    """
    # Depth-agnostic on purpose. This was a fixed-depth
    # "*/calibBE/*/*_shoreline_change_rate.csv" glob, which reaches
    # <period>/calibBE/<run>/ but NOT the <arm>/<period>/calibBE/<run>/ a run
    # forced off the calibration wave climate is filed under. A missed run can
    # only ever make `newest` OLDER, so the failure is the staleness check
    # passing a figure it should have caught -- silently, and in the direction
    # that publishes the stale file.
    #
    # Runs are found by their metadata file, which stays at the run folder's
    # root under its full name; the rate CSV itself has moved into tables/ and
    # dropped the prefix, so globbing for it would miss a migrated run -- again
    # in the direction that publishes the stale file. run_layout.resolve reads
    # either layout.
    runs = []
    for meta in RAW_RUNS.rglob("*_run_metadata.json"):
        if meta.parent.parent.name != "calibBE":
            continue
        run_name = meta.name[: -len("_run_metadata.json")]
        csv_path = resolve_run_file(meta.parent, "rate_csv", run_name)
        if csv_path.is_file():
            runs.append(csv_path)
    if not runs:
        return [], None
    newest = max(r.stat().st_mtime for r in runs)
    return [p for p in paths if p.exists() and p.stat().st_mtime < newest], newest


def analysis_module():
    spec = importlib.util.spec_from_file_location(
        "_loess",
        _HERE.parent.parent / "2-calibrate" / "HAT_be_zone_residual_fit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def render_dict(period, rates, module):
    """The dict file, in the shape of the files it replaces."""
    lines = [
        f"# Converged source/sink field, {period}-{period + 20}.",
        "#",
        "# GENERATED by scripts/input_prep/7-source-sink/4-export/HAT_export_be_calibration.py",
        "# from hatteras_site_config.py, which is the source of truth. Do not",
        "# edit this file: edit the config, or re-run the calibration.",
        "#",
        "# Frozen-zone masked iteration, converged "
        f"{date.today().isoformat()}. Domains outside the frozen zone set carry",
        "# 0.0 because they were never eligible for correction, NOT because",
        "# their residual was zero -- see be_calibration_domains.csv.",
        "",
        "DOMAIN_BE_RATES_CALIBRATED = {",
    ]
    for gis in range(1, 91):
        rate, label = rates.get(gis, (0.0, ""))
        line = "    %3d: %+.1f," % (gis, rate)
        notes = []
        if label:
            notes.append(label)
        if gis in module.GROIN_RESERVED_DOMAINS:
            notes.append("RESERVED for the groin module")
        if notes:
            line += "  # " + " | ".join(notes)
        lines.append(line)
    lines.append("}")
    return "\n".join(lines) + "\n"


def build_table(module, final, pass0, metrics):
    """One row per domain: eligibility, rates, what iteration added, residual."""
    rows = []
    for gis in range(1, 91):
        row = {"domain": gis,
               "physical_zone": module.assign_physical_zone(gis)}
        if gis in module.LOCKED_DOMAINS:
            status = "locked (solved separately)"
        elif gis in module.GROIN_RESERVED_DOMAINS:
            status = "reserved for groin module"
        elif (gis in module.FROZEN_ZONE_DOMAINS[1984]
              or gis in module.FROZEN_ZONE_DOMAINS[2004]):
            status = "correctable"
        else:
            status = "withheld (outside frozen zone set)"
        row["status"] = status
        for period, tag in ((1984, "1984_2004"), (2004, "2004_2024")):
            row[f"eligible_{tag}"] = gis in module.FROZEN_ZONE_DOMAINS[period]
            f = final[period].get(gis, (0.0, ""))[0]
            row[f"be_final_{tag}"] = round(f, 3)
            if pass0 is None:
                row[f"be_pass0_{tag}"] = ""
                row[f"iteration_added_{tag}"] = ""
            else:
                p = pass0[period].get(gis, (0.0, ""))[0]
                row[f"be_pass0_{tag}"] = round(p, 3)
                row[f"iteration_added_{tag}"] = round(f - p, 3)
        if metrics is not None and gis in metrics.index:
            row["residual_1984_2004"] = round(
                float(metrics.loc[gis, "raw_residual_p1"]), 3)
            row["residual_2004_2024"] = round(
                float(metrics.loc[gis, "raw_residual_p2"]), 3)
        rows.append(row)
    return pd.DataFrame(rows)


def caveats_block():
    """What this export could NOT carry, and why.

    Generated rather than hand-written into the README, because the README is
    rewritten on every export and a hand-added note would vanish on the next
    run -- which is how the file came to disagree with the config in the first
    place.
    """
    lines = []
    if not PASS0_BACKUP.exists():
        lines.append(
            f"* **`be_pass0_*` and `iteration_added_*` are empty.** They split "
            f"each final rate into the one-shot solve and what the iteration "
            f"added, and that split can only come from `{PASS0_BACKUP.name}`, "
            f"the field written before the first pass. That file is not on "
            f"disk and was never committed, so it cannot be recovered. The "
            f"FINAL values are unaffected -- they come from the config.")
    stale, newest = stale_against_runs(
        [FIGURE_DIR / name for name, _ in FIGURES])
    if stale:
        from datetime import datetime
        when = datetime.fromtimestamp(newest).strftime("%Y-%m-%d")
        names = ", ".join(f"`{pathlib.Path(x).name}`" for x in stale)
        lines.append(
            f"* **{names} predate(s) the newest calibrated run ({when}).** "
            f"`fig_be_zones_and_corrections.png` cannot be redrawn for the "
            f"same reason the pass-0 columns are empty: its lower panels need "
            f"that lost backup. It still shows the masked iteration that "
            f"produced these values, which has not been re-run -- but it is "
            f"older than the runs and is marked here rather than passed off "
            f"as current.")
    if not lines:
        return ""
    body = (chr(10) + "## What this export could not carry" + chr(10)
            + chr(10))
    return body + (chr(10) + chr(10)).join(lines) + chr(10)


def readme(module, table, history):
    figure_rows = "; ".join(f"`{name}` - {what}"
                            for name, what in FIGURES)
    caveats = caveats_block()
    p1 = history["passes"]["1984_2004"]
    p2 = history["passes"]["2004_2024"]
    # Computed, not typed: these read 42/57 for the lineage retired on
    # 2026-09-14 and were still saying so after the field had been re-derived.
    # A number stated in prose beside the table it contradicts is worse than no
    # number at all.
    _edge = history["baselines"]["edgeBE"]
    closed_p1 = 100.0 * (_edge["1984_2004"] - p1[0]["rmse"]) / _edge["1984_2004"]
    closed_p2 = 100.0 * (_edge["2004_2024"] - p2[0]["rmse"]) / _edge["2004_2024"]
    n_corr = int((table["status"] == "correctable").sum())
    n_held = int((table["status"] == "withheld (outside frozen zone set)").sum())
    return f"""# Source/sink (background erosion) calibration — converged field

Generated by `scripts/input_prep/7-source-sink/4-export/HAT_export_be_calibration.py` from
`scripts/site_layer/hatteras_site_config.py`, which is what the model actually imports.
**Do not edit these files** — edit the config, or re-run the calibration and
re-export.

## What is here

| file | what it is |
|---|---|
| `4-export/be_rates_1984_2004.py` | the converged field, period 1 |
| `4-export/be_rates_2004_2024.py` | the converged field, period 2 |
| `4-export/be_calibration_domains.csv` | per domain: zone, eligibility, pass-0 and final rate, what the iteration added, residual still standing |
| `2-calibrate/1984_2004__2004_2024/` | what the fit itself wrote: `DOMAIN_BE_RATES*.txt`, `be_zone_metrics.csv`, `cascade_base_lrr.csv`, and `convergence_history.json` — every pass, both baselines, and the abandoned unmasked attempt |
| `2-calibrate/1996_2010__2010_2024/` | the same for the 1996/2010 pair, which is fitted but not exported |
| `2-calibrate/prebe/` | the config as it stood before each apply pass; the pass-0 split reads one of these |
| `3-figures/1984_2004__2004_2024/` | {figure_rows} |
| `3-figures/1996_2010__2010_2024/` | `1-field/` for the 1996/2010 pair |
| `archive/` | the `superseded_*` folders: the 2026-06-15 files this replaces, kept because earlier runs were built against them, and the retired 08-24 lineage |

Every period pair has its own folder under `2-calibrate/` and `3-figures/`,
named `<p1start>_<p1end>__<p2start>_<p2end>`; the numbers match the script
steps in `scripts/input_prep/7-source-sink/`. Resolve these paths through
`scripts/site_layer/hat_source_sink.py`, never by typing them (2026-09-18,
when the default pair stopped writing to the unlabelled root).
There is no `1-` folder because step 1, `1-prepare/`, writes into each run's
folder under `output/raw_runs/`, not here.

The two `be_rates_*.py` are **generated data, not code** — a header comment
and one dict literal, nothing executable. Nothing imports them. They are
Python because the field is a Python dict in the config, so this form pastes
straight back in and carries each domain's zone label as a trailing comment.
`be_calibration_domains.csv` is the machine-readable form of the same thing.

## How the field was produced

The calibration measures the residual of a base run against the CoastSat LRR
target and imposes it as background erosion. That assumes giving a domain
X m/yr moves its shoreline rate by X m/yr — and it does not, because BRIE
diffuses an imposed rate alongshore. One pass closes only {closed_p1:.0f}% (period 1) and
{closed_p2:.0f}% (period 2) of the misfit.

So the solve is **iterated**: each pass re-measures what the current field
leaves and adds it, which converges without needing to know the surviving
fraction. That fraction is not a constant — roughly 0.8–1.2 for a contiguous
block of corrections, roughly 0.1 for one alternating at the grid scale — so
dividing by it instead would have amplified narrow features tenfold into rates
indefensible as sediment fluxes.

| period | edgeBE | pass 0 | converged | passes |
|---|---|---|---|---|
| 1984–2004 | {history['baselines']['edgeBE']['1984_2004']:.4f} | {p1[0]['rmse']:.4f} | **{p1[-1]['rmse']:.4f}** | {len(p1) - 1} |
| 2004–2024 | {history['baselines']['edgeBE']['2004_2024']:.4f} | {p2[0]['rmse']:.4f} | **{p2[-1]['rmse']:.4f}** | {len(p2) - 1} |

RMSE of `lrr_m_yr` against the CoastSat target, D2–D89, on the
full-management groin-on run of each period. Stopped when a pass bought under
5%.

## Zone membership is fixed, and that is the point

Zone identification is the scientific step — it says a stretch of coast has a
real sediment-budget deficit and names the process. Magnitude is arithmetic.
Iterating both lets the arithmetic rewrite the science: re-deriving zones each
pass let progressively less coherent features cross the significance threshold
as the real ones were satisfied, and because adding background erosion at a
domain pushes sediment into its neighbours, later passes began correcting the
spillover of earlier ones.

So zones were identified **once**, from the pass-0 residual under the ordinary
significance and width rules, and then held fixed (`FROZEN_ZONE_DOMAINS`).

- **{n_corr} domains correctable** — inside the frozen zone set
- **{n_held} domains withheld** — outside it, and left at 0.0 *however large
  their residual*, which is honest unexplained variance rather than a fitted
  constant
- **{len(module.GROIN_RESERVED_DOMAINS)} domains reserved** (D5–D7) — the
  Buxton groin's footprint
- **2 domains locked** (D1, D90) — solved separately by buffer-cell
  reproduction

An unmasked run of the same iteration scored better
({history['_abandoned_unmasked']['1984_2004']:.4f} /
{history['_abandoned_unmasked']['2004_2024']:.4f}) and was abandoned. The gap
is the fit available only by correcting outside justifiable zones.
Those two numbers were measured on 2026-08-24, against the lineage retired in
`archive/superseded_20260914/`, and have not been re-measured: the unmasked variant has
no zone set by construction, so the zone-set correction does not apply to it.
They size the declined fit; they are not a current score.

## Caveats

- **A 0.0 does not mean "no residual".** It means the domain was never
  eligible. The CSV carries the residual for every domain so the two can be
  told apart.
- **D5–D7 are reserved, not fitted.** The residual there is the groin module's
  shortfall — too little fillet built in period 1, no release in period 2 —
  and absorbing it into the sediment budget would double-count against the
  M/f fit.
- **D1 and D90 are not sediment budgets.** They are boundary absorbers, and
  carry rates about ten times the interior because only ~10% of an imposed
  edge rate survives diffusion.
- **This covers the two CALIBRATED periods only.** 1996-2010 and 2010-2024 are
  wired in `HATTERAS_PERIODS` but carry no interior fit, so there is nothing to
  export for them. 1996 has solved end domains; see `HATTERAS_BE_EDGE_ONLY`.
{caveats}"""


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true",
                        help="report what would be written, write nothing")
    parser.add_argument("--allow-stale", action="store_true",
                        help="export even when a figure predates the newest "
                             "calibBE run. Only for the case where you know a "
                             "figure does not depend on run output.")
    args = parser.parse_args()

    module = analysis_module()
    final = {p: rates_from(CONFIG, p) for p in (1984, 2004)}
    # THE PASS-0 FIELD IS OPTIONAL NOW (2026-09-12). It used to be required,
    # which made this exporter unrunnable: the backup it names was never
    # committed and is not on disk, so the two columns it feeds cannot be
    # reconstructed. Refusing to export at all meant the VALUES stayed stale
    # for the sake of a provenance column -- and stale values are the failure
    # this file exists to prevent. Absent, the pass-0 and iteration-added
    # columns are written empty and the README says why.
    pass0 = ({p: rates_from(PASS0_BACKUP, p) for p in (1984, 2004)}
             if PASS0_BACKUP.exists() else None)
    if pass0 is None:
        print(f"  NOTE: {PASS0_BACKUP.name} is absent, so be_pass0_* and "
              f"iteration_added_* are left empty.")

    metrics_path = CALIB_OUT / "be_zone_metrics.csv"
    metrics = (pd.read_csv(metrics_path).set_index("domain")
               if metrics_path.exists() else None)
    history = json.loads((CALIB_OUT / "convergence_history.json")
                         .read_text(encoding="utf-8"))

    table = build_table(module, final, pass0, metrics)
    counts = table["status"].value_counts()
    print("  domain status")
    for name, n in counts.items():
        print(f"    {name:<38} {n:>3}")
    if pass0 is not None:
        for period, tag in ((1984, "1984_2004"), (2004, "2004_2024")):
            moved = table[table[f"iteration_added_{tag}"].abs() > 1e-9]
            print(f"  {tag}: iteration changed {len(moved)} domains, "
                  f"mean |added| "
                  f"{moved[f'iteration_added_{tag}'].abs().mean():.3f}, "
                  f"max {moved[f'iteration_added_{tag}'].abs().max():.2f} m/yr")

    # BEFORE anything is written, and before --check returns, so a dry run
    # reports the same refusal a real one would. Checking it later left a
    # half-finished export on disk: values written, figures not.
    stale, newest = (([], None) if args.allow_stale
                     else stale_against_runs([FIGURE_DIR / name
                                              for name, _ in FIGURES]))
    if stale:
        from datetime import datetime
        when = datetime.fromtimestamp(newest).strftime("%Y-%m-%d %H:%M")
        message = [
            "",
            "REFUSING TO EXPORT: "
            f"{len(stale)} figure(s) predate the newest calibBE run ({when}):",
        ]
        message += [f"    {path.name}" for path in stale]
        message += [
            "",
            "  They would be published as if current. Regenerate first:",
            "    HAT_BE_BASE_PRESET=calibBE python "
            "2-calibrate/HAT_be_zone_residual_fit.py",
            "    python 3-figures/HAT_plot_be_convergence.py",
            "    python 3-figures/HAT_plot_be_zones.py",
            "    python 3-figures/HAT_plot_groin_reserved_residual.py",
            "  then re-run this. Use --allow-stale only if you mean to ship "
            "figures older than the runs.",
        ]
        raise SystemExit(NL.join(message))

    if args.check:
        print("\n  --check: nothing written")
        return 0

    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    superseded = _be.ARCHIVE / f"superseded_{date.today():%Y%m%d}"
    for name in ("1984_2004_values", "2004_2024_values"):
        old = DATA_DIR / name
        if old.exists():
            superseded.mkdir(parents=True, exist_ok=True)
            shutil.move(str(old), str(superseded / name))
            print(f"  moved {name} -> archive/{superseded.name}/")

    for period, tag in ((1984, "1984_2004"), (2004, "2004_2024")):
        path = EXPORT_DIR / f"be_rates_{tag}.py"
        path.write_text(render_dict(period, final[period], module),
                        encoding="utf-8")
        print(f"  wrote {path.name}")

    # 3-figures/ IS the record, so there is nothing to copy into it -- the
    # analysis scripts write here directly. This block now only reports what
    # is present, which is what the README claims.
    figure_dir = FIGURE_DIR
    figure_dir.mkdir(parents=True, exist_ok=True)
    missing = [name for name, _ in FIGURES
               if not (figure_dir / name).exists()]
    print(f"  3-figures/{_be.DEFAULT_PAIR_TAG}/ holds {len(FIGURES) - len(missing)} of "
          f"{len(FIGURES)} record figures")
    if missing:
        # Named rather than skipped silently: a figure absent from the record
        # is indistinguishable from one that was never made.
        print(f"    MISSING, not copied: {', '.join(missing)}")

    table.to_csv(EXPORT_DIR / "be_calibration_domains.csv", index=False)
    print("  wrote be_calibration_domains.csv")

    # convergence_history.json is NOT copied up (dropped 2026-09-14). It lives
    # once, in 2-calibrate/, beside the calibration that wrote it. The copy
    # that used to sit here was byte-identical the day it was made and would
    # have diverged the first time a pass ran without a re-export -- at which
    # point a reader has two files with one name and no way to tell which the
    # figures were drawn from. The README points at the one copy instead.
    (DATA_DIR / "README.md").write_text(readme(module, table, history),
                                        encoding="utf-8")
    print("  wrote README.md")
    print(f"\n  -> {DATA_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
