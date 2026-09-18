#!/usr/bin/env python3
r"""
HAT_list_runs.py
==============================================================================
What is under output/raw_runs, grouped by whatever makes runs comparable.

WHY A LISTING RATHER THAN A FOLDER
    98 of 164 runs sit on a topography that is no longer CURRENT, and opening
    a preset folder does not say which. The obvious fix is to nest by
    topography version, and it is the wrong one: 2004-start has only one
    version, so 63 runs would gain a level that says nothing; a sweep run is
    already five levels down; and topography is not the only axis that decides
    whether two runs are comparable, so nesting one of them just moves the
    question.

    Every run records all of it -- topo_product, topo_dune_version,
    be_values_digest -- in its metadata and in run_index.csv. So the question
    "which runs are on v1" is a query, and this is the query.

    `--stamp` writes the same facts as a one-line BUILT_ON.txt inside each run
    folder, so a folder opened on its own also answers it.

USAGE
    python HAT_list_runs.py                    # by topography (the default)
    python HAT_list_runs.py --by erosion
    python HAT_list_runs.py --by preset --detail
    python HAT_list_runs.py --only-stale       # just what is not CURRENT
    python HAT_list_runs.py --stamp            # write BUILT_ON.txt per run
==============================================================================
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
RAW_RUNS = REPO / "output" / "raw_runs"
import sys as _b3dsys
from pathlib import Path as _B3DP
_b3dsys.path.insert(0, str(next(_q for _q in _B3DP(__file__).resolve().parents
                                if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_topo_version as _b3d  # noqa: E402
DOMAIN_ROOT = _b3d.DOMAIN_ROOT
STAMP_NAME = "BUILT_ON.txt"


def current_versions() -> dict:
    out = {}
    if not DOMAIN_ROOT.is_dir():
        return out
    for product in sorted(p.name for p in DOMAIN_ROOT.iterdir() if p.is_dir()):
        marker = DOMAIN_ROOT / product / "dune-topo" / "CURRENT"
        if marker.is_file():
            out[product] = marker.read_text(encoding="utf-8").strip()
    return out


def load_runs() -> list:
    runs = []
    for meta in sorted(RAW_RUNS.rglob("*_run_metadata.json")):
        try:
            d = json.load(open(meta, encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        ident, period, source = (d.get("identity", {}), d.get("period", {}),
                                 d.get("source/sink", {}))
        runs.append(dict(
            dir=meta.parent,
            rel=meta.parent.relative_to(RAW_RUNS).as_posix(),
            name=ident.get("run_name") or meta.parent.name,
            timestamp=ident.get("timestamp") or "",
            product=ident.get("topo_product"),
            topo=ident.get("topo_dune_version"),
            commit=(ident.get("git_commit") or "")[:8],
            dirty=bool(ident.get("git_dirty")),
            preset=source.get("preset"),
            digest=source.get("values_digest"),
            period=f"{period.get('start_year')}_{period.get('end_year')}",
        ))
    return runs


def newest_digest_per_preset(runs) -> dict:
    """The digest the most recent run of each (preset, period) carries.

    PER PERIOD, not per preset: the background-erosion field is calibrated for
    one period at a time, so one preset legitimately has a different digest in
    1984-2004 than in 2004-2024. Keying on the preset alone called that
    by-design difference "superseded".
    """
    newest = {}
    for r in runs:
        key = (r["preset"], r["period"])
        if r["timestamp"] > newest.get(key, ("", ""))[0]:
            newest[key] = (r["timestamp"], r["digest"])
    return {k: d for k, (_, d) in newest.items()}


def group_key(run, by, current, newest):
    if by == "topo":
        cur = current.get(run["product"])
        tag = ("CURRENT" if run["topo"] == cur
               else f"not CURRENT, which is {cur}" if cur else "no CURRENT marker")
        return f"{run['product']} / {run['topo']}   ({tag})"
    if by == "erosion":
        key = (run["preset"], run["period"])
        tag = "current" if run["digest"] == newest.get(key) else "superseded"
        return f"{run['preset']} / {run['period']} / {run['digest']}   ({tag})"
    if by == "preset":
        return f"{run['period']} / {run['preset']}"
    return run["period"]


def parent_bucket(rel: str) -> str:
    """The folder a run sits in, collapsed to what says its purpose.

    Purpose layout (2026-09-16): matrix/<period>/<preset>, sensitivity/<axis>,
    experiments/<tag>, versions/<tag>. The 09-10 layout's sweeps/<family>
    collapses to 'sweeps' as before.
    """
    parts = rel.split("/")[:-1]
    if parts and parts[0] in ("sensitivity", "experiments", "versions", "archive"):
        keep = 2
        # a two-level tag: <set>/<member>
        if len(parts) > 2 and not parts[2][:4].isdigit():
            keep = 3
        parts = parts[:keep]
    if "sweeps" in parts:
        parts = parts[: parts.index("sweeps") + 1]
    return "/".join(parts) or "."


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--by", default="topo",
                    choices=("topo", "erosion", "preset", "period"),
                    help="what to group by (default: topo)")
    ap.add_argument("--detail", action="store_true", help="list every run name")
    ap.add_argument("--only-stale", action="store_true",
                    help="only runs not on their product's CURRENT topography")
    ap.add_argument("--stamp", action="store_true",
                    help=f"write {STAMP_NAME} into each run folder and exit")
    args = ap.parse_args()

    runs = load_runs()
    if not runs:
        raise SystemExit(f"no runs found under {RAW_RUNS}")
    current = current_versions()
    newest = newest_digest_per_preset(runs)

    if args.stamp:
        for r in runs:
            cur = current.get(r["product"])
            state = ("the CURRENT topography" if r["topo"] == cur
                     else f"NOT CURRENT -- {r['product']} is now {cur}")
            field = ("the current background-erosion field for this preset "
                     "and period"
                     if r["digest"] == newest.get((r["preset"], r["period"]))
                     else "an OLDER background-erosion field than the newest "
                          "run of this preset and period")
            (r["dir"] / STAMP_NAME).write_text(
                f"{r['name']}\n"
                f"  run        {r['timestamp']}\n"
                f"  topography {r['product']}/{r['topo']}  ({state})\n"
                f"  source/sink {r['preset']}, digest {r['digest']}\n"
                f"             ({field})\n"
                f"  commit     {r['commit']}"
                f"{'  DIRTY TREE' if r['dirty'] else ''}\n"
                f"\n"
                f"Written by HAT_list_runs.py --stamp. Regenerate it after the\n"
                f"CURRENT marker or the calibration moves; it describes how this\n"
                f"run compares with today's, not what it contains.\n",
                encoding="utf-8")
        sys.stdout.write(f"wrote {STAMP_NAME} into {len(runs)} run folder(s)\n")
        return

    if args.only_stale:
        runs = [r for r in runs
                if r["product"] in current and r["topo"] != current[r["product"]]]
        if not runs:
            sys.stdout.write("every run is on its product's CURRENT topography\n")
            return

    groups = defaultdict(list)
    for r in runs:
        groups[group_key(r, args.by, current, newest)].append(r)

    total = 0
    for key in sorted(groups):
        members = groups[key]
        total += len(members)
        sys.stdout.write(f"\n{key}   {len(members)} run(s)\n")
        for bucket, n in sorted(Counter(parent_bucket(r["rel"]) for r in members).items()):
            sys.stdout.write(f"    {bucket:<48} {n}\n")
        if args.detail:
            for r in sorted(members, key=lambda r: r["rel"]):
                sys.stdout.write(f"        {r['name']}\n")
    sys.stdout.write(f"\n{total} run(s)\n")


if __name__ == "__main__":
    main()
