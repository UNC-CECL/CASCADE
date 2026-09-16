#!/usr/bin/env python3
r"""
HAT_run_supersession_report.py
==============================================================================
Which runs under output/raw_runs are candidates for retirement, and why.

READ-ONLY. It moves nothing and deletes nothing. Hannah's instruction on
2026-09-10 was "flag candidates, retire nothing yet", and a script that can
only report cannot be misread as one that acts.

WHAT "SUPERSEDED" MEANS HERE, AND WHAT IT DOES NOT
    A run is flagged when something it was built on has since moved, which
    makes it incomparable with a run made today. That is not the same as
    wrong: a run is a faithful record of the inputs it had, and the two
    1984-2004 archives were kept for exactly that reason. The judgement about
    whether an incomparable run is still worth keeping is the reader's.

THE CHECKS
    topography   the run's dune-topo version against the product's CURRENT.
                 The 2026-09-03 re-pick moved 1984-start from v1 to v2, so a
                 v1 run measures a different island from a run made today.
    calibration  the source/sink `values_digest` within one preset AND ONE
                 PERIOD. The field is calibrated per period, so a preset having
                 a different digest in each period is by design; two digests
                 inside one period would mean it was recalibrated and only some
                 runs remade. (Grouping on the preset alone reported the
                 by-design difference as 55 superseded runs on 2026-09-10.)
    duplicates   one run name in more than one place. Usually legitimate --
                 an arm is a different forcing of the same scenario -- so
                 these are listed for orientation, not flagged.

WHERE   output/raw_runs/SUPERSEDED_CANDIDATES.md

USAGE
    python HAT_run_supersession_report.py
    python HAT_run_supersession_report.py --print
==============================================================================
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
RAW_RUNS = REPO / "output" / "raw_runs"
DOMAIN_ROOT = REPO / "data" / "hatteras_init" / "1-barrier3d-domains"
REPORT = RAW_RUNS / "SUPERSEDED_CANDIDATES.md"
# An arm component that names a dune-topo version, e.g. the v3 of
# arms/version-pair/v3.
_VERSION_TOKEN = re.compile(r"v\d+")


def current_versions() -> dict:
    """The dune-topo version each product's CURRENT marker names."""
    out = {}
    for product in sorted(p.name for p in DOMAIN_ROOT.iterdir() if p.is_dir()):
        marker = DOMAIN_ROOT / product / "dune-topo" / "CURRENT"
        if marker.is_file():
            out[product] = marker.read_text(encoding="utf-8").strip()
    return out


def versions_on_disk(product: str) -> list:
    root = DOMAIN_ROOT / product / "dune-topo"
    if not root.is_dir():
        return []
    return sorted(p.name for p in root.iterdir() if p.is_dir() and p.name.startswith("v"))


def expected_version(run, current):
    """The version this run SHOULD be on, which is not always CURRENT.

    An arm can name a version -- `version-pair/v3` holds v3 against v2, and
    `behindroad-copy` was built on the v3 footprint layer. Those runs are on
    that version deliberately, so measuring them against CURRENT reports a
    deliberate choice as drift, and acting on it would destroy the comparison
    the arm exists for (flagged 2026-09-11).
    """
    for part in run["rel"].split("/"):
        if _VERSION_TOKEN.fullmatch(part):
            return part
    return current.get(run["product"])


def _is_stale(run, current):
    want = expected_version(run, current)
    return want is not None and run["topo"] != want


def load_runs() -> list:
    runs = []
    for meta in sorted(RAW_RUNS.rglob("*_run_metadata.json")):
        try:
            d = json.load(open(meta, encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        ident = d.get("identity", {})
        period = d.get("period", {})
        source = d.get("source/sink", {})
        runs.append(dict(
            rel=meta.parent.relative_to(RAW_RUNS).as_posix(),
            name=ident.get("run_name") or meta.parent.name,
            timestamp=ident.get("timestamp") or "",
            topo=ident.get("topo_dune_version"),
            product=ident.get("topo_product"),
            commit=(ident.get("git_commit") or "")[:8],
            dirty=bool(ident.get("git_dirty")),
            preset=source.get("preset"),
            digest=source.get("values_digest"),
            start_year=period.get("start_year"),
            period=f"{period.get('start_year')}_{period.get('end_year')}",
            state_bytes=_state_bytes(meta.parent),
        ))
    return runs


# ARMS THAT NAME A VERSION are judged against that version, not CURRENT: they
# exist to compare two islands, so their state is kept whatever CURRENT says.
VERSION_ARMS = ("version-pair", "behindroad-copy")


def _state_bytes(run_dir: Path) -> int:
    """Size of this run's model state, or 0 if it was not saved.

    The .npz is the ONLY artifact that lets a deep-dive figure be re-derived
    without re-running. Everything else -- the rate table, the shoreline
    matrix, the metadata -- is written regardless and is small.
    """
    return sum(f.stat().st_size for f in run_dir.glob("*.npz"))


def state_verdict(run, current) -> str:
    """Does the keep rule keep this run's model state?

    THE RULE (Hannah, 2026-09-14): keep the state of runs on the CURRENT
    topography, and of arms that name a version deliberately. Drop the rest.

    Why it is defensible: a run on an island a re-pick behind cannot be
    compared against a run made today, so re-plotting it deeply is answering a
    question nobody can ask. Its tables and metadata survive either way.

    Why it is not automatic: a run is one to three minutes of compute, but
    re-running reproduces it at TODAY'S code, not the code it was made with.
    That is what the state is really insuring against, and it is why nothing
    here deletes anything.
    """
    if not run["state_bytes"]:
        return "none"
    if any(arm in run["rel"] for arm in VERSION_ARMS):
        return "keep (version arm)"
    if run["topo"] == current.get(run["product"]):
        return "keep (current topography)"
    return "would free"



def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--print", action="store_true", dest="echo",
                    help="also write the report to stdout")
    args = ap.parse_args()

    runs = load_runs()
    if not runs:
        raise SystemExit(f"no runs found under {RAW_RUNS}")
    current = current_versions()

    # --- topography ------------------------------------------------------
    stale_topo = [r for r in runs if _is_stale(r, current)]
    topo_groups = defaultdict(list)
    for r in stale_topo:
        topo_groups[(r["product"], r["topo"])].append(r)

    # --- calibration ------------------------------------------------------
    # WITHIN A PERIOD. The background-erosion field is calibrated per period,
    # so one preset legitimately carries a different digest for 1984-2004 than
    # for 2004-2024. Grouping on the preset alone reported that by-design
    # difference as 55 superseded runs on 2026-09-10, which was wrong; the
    # split that would matter is two digests for one preset in ONE period.
    by_preset = defaultdict(Counter)
    newest = {}
    for r in runs:
        key = (r["preset"], r["period"])
        by_preset[key][r["digest"]] += 1
        if r["timestamp"] > newest.get((key, r["digest"]), ""):
            newest[(key, r["digest"])] = r["timestamp"]
    split_presets = {k: c for k, c in by_preset.items() if len(c) > 1}

    # --- duplicates -------------------------------------------------------
    by_name = defaultdict(list)
    for r in runs:
        by_name[r["name"]].append(r)
    duplicates = {n: rs for n, rs in by_name.items() if len(rs) > 1}

    L = []
    w = L.append
    w("# Runs that may be superseded")
    w("")
    w(f"Written {datetime.now():%Y-%m-%d %H:%M} by `HAT_run_supersession_report.py`, "
      f"over {len(runs)} runs under `output/raw_runs`. **Nothing here has been "
      f"moved or deleted.** A flag means a run is no longer comparable with one "
      f"made today, not that it is wrong: a run is a faithful record of the "
      f"inputs it had.")
    w("")

    # 1
    w("## 1. Built on a topography that is no longer CURRENT")
    w("")
    w("The dune-topo version each product now points at:")
    w("")
    w("| product | CURRENT | versions on disk |")
    w("|---|---|---|")
    for product, cur in sorted(current.items()):
        w(f"| `{product}` | {cur} | {', '.join(versions_on_disk(product)) or 'none'} |")
    w("")
    if not stale_topo:
        w("Every run is on its product's CURRENT version.")
    else:
        w(f"**{len(stale_topo)} of {len(runs)} runs** are on another version:")
        w("")
        w("| product | run version | runs | newest | CURRENT |")
        w("|---|---|---|---|---|")
        for (product, topo), rs in sorted(topo_groups.items()):
            w(f"| `{product}` | {topo} | {len(rs)} | "
              f"{max(r['timestamp'] for r in rs)[:10]} | {current[product]} |")
        w("")
        w("A v1 1984-start run measures a different island from one made today: "
          "the 2026-09-03 re-pick changed the dune search windows, and with them "
          "the interiors and the road setbacks. Comparing a v1 run against a v2 "
          "run attributes the pick difference to whatever the figure is about.")
        w("")
        w("**A run that names a version is judged against that version, not "
          "CURRENT.** `versions/version-pair/v3` holds v3 against v2 and "
          "`experiments/2026-09-08-behindroad-copy` was built on the v3 footprint layer; both are "
          "on v3 deliberately. Re-running them on CURRENT would destroy the "
          "comparison they exist for, so they are not listed above.")
    w("")

    # 2
    w("## 2. Forced with a different background-erosion field")
    w("")
    w("The field is calibrated PER PERIOD, so one preset carrying a different "
      "digest for 1984-2004 than for 2004-2024 is by design. What would matter "
      "is two digests for one preset within ONE period: that means the field "
      "was recalibrated and only some runs were remade.")
    w("")
    w("| preset | period | digest | runs |")
    w("|---|---|---|---|")
    for (preset, period), counter in sorted(by_preset.items()):
        for digest, n in sorted(counter.items()):
            w(f"| `{preset}` | {period} | `{digest}` | {n} |")
    w("")
    if not split_presets:
        w("**No split.** Every preset carries exactly one digest within each "
          "period, so no run is forced with a superseded field.")
    else:
        for (preset, period), counter in sorted(split_presets.items()):
            ranked = sorted(counter.items(),
                            key=lambda kv: newest.get(((preset, period), kv[0]), ""),
                            reverse=True)
            w(f"**`{preset}` in {period} carries {len(counter)} digests** -- "
              f"the older one(s) are superseded:")
            w("")
            for rank, (digest, n) in enumerate(ranked):
                note = "current" if rank == 0 else "superseded"
                w(f"- `{digest}`, {n} runs, newest "
                  f"{newest.get(((preset, period), digest), '')[:10]} ({note})")
            w("")
    w("")

    # 3
    w("## 3. One name in more than one place")
    w("")
    w("Usually legitimate: an arm is a different forcing of the same scenario, "
      "so the name is meant to repeat. Listed for orientation.")
    w("")
    if not duplicates:
        w("Every run name is unique.")
    for name, rs in sorted(duplicates.items()):
        w(f"**`{name}`**")
        w("")
        w("| when | topo | where |")
        w("|---|---|---|")
        for r in sorted(rs, key=lambda r: r["timestamp"]):
            w(f"| {r['timestamp'][:16]} | {r['topo']} | `{r['rel']}` |")
        w("")

    # 4
    # ---- model state ----------------------------------------------------
    w("## 4. Model state on disk")
    w("")
    by_verdict = {}
    for r in runs:
        v = state_verdict(r, current)
        if v == "none":
            continue
        by_verdict.setdefault(v, []).append(r)
    total = sum(r["state_bytes"] for rs in by_verdict.values() for r in rs)
    w(f"**{sum(len(v) for v in by_verdict.values())} runs carry a model state "
      f"file, {total / 1e9:.1f} GB in total.** The `.npz` is the only artifact "
      f"that lets a deep-dive figure be re-derived without re-running; the rate "
      f"table, the shoreline matrix and the metadata are written regardless.")
    w("")
    w("**The keep rule (2026-09-14):** keep the state of runs on the CURRENT "
      "topography, and of arms that name a version deliberately. A run on an "
      "island a re-pick behind cannot be compared against one made today, so "
      "re-plotting it deeply answers a question nobody can ask.")
    w("")
    w("| verdict | runs | GB |")
    w("|---|---|---|")
    for verdict in sorted(by_verdict):
        rs = by_verdict[verdict]
        w(f"| {verdict} | {len(rs)} | {sum(r['state_bytes'] for r in rs)/1e9:.1f} |")
    w("")
    free = by_verdict.get("would free", [])
    if free:
        w(f"### The {len(free)} the rule would free")
        w("")
        w("Listed so the decision is reviewable. **Nothing deletes these** -- "
          "they are git-ignored, so removal cannot be undone, and it stays a "
          "deliberate act taken after reading this.")
        w("")
        for r in sorted(free, key=lambda r: -r["state_bytes"]):
            w(f"* `{r['rel']}` — {r['state_bytes']/1e6:.0f} MB, "
              f"{r['product']}/{r['topo']}")
        w("")

    w("## 5. Provenance")
    w("")
    dirty = sum(1 for r in runs if r["dirty"])
    w(f"- **{dirty} of {len(runs)}** runs were made from a dirty working tree, "
      f"so their commit alone does not reproduce them.")
    commits = Counter(r["commit"] for r in runs if r["commit"])
    w(f"- {len(commits)} distinct commits across the tree; the most used is "
      f"`{commits.most_common(1)[0][0]}` with {commits.most_common(1)[0][1]} runs.")
    dates = Counter(r["timestamp"][:10] for r in runs if r["timestamp"])
    w(f"- Run dates: " + ", ".join(f"{d} ({n})" for d, n in sorted(dates.items())) + ".")
    w("")
    w("## What this report does not decide")
    w("")
    w("Whether an incomparable run is worth keeping. Two 1984-2004 archives "
      "were deliberately kept when they were superseded, because a superseded "
      "run is still the only record of what those inputs produced. Retiring "
      "anything is a separate, deliberate step.")
    w("")

    text = "\n".join(L)
    REPORT.write_text(text, encoding="utf-8")
    if args.echo:
        sys.stdout.write(text + "\n")
    sys.stdout.write(f"\nwrote {REPORT}\n")
    sys.stdout.write(f"  {len(stale_topo)} run(s) on a non-CURRENT topography\n")
    sys.stdout.write(f"  {len(split_presets)} preset(s) carrying more than one "
                     f"background-erosion field\n")
    sys.stdout.write(f"  {len(duplicates)} run name(s) in more than one place\n")


if __name__ == "__main__":
    main()
