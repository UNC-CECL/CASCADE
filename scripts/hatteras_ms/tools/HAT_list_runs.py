#!/usr/bin/env python3
"""
What is under output/raw_runs, grouped by whatever makes runs comparable.

    python scripts/hatteras_ms/tools/HAT_list_runs.py
    python scripts/hatteras_ms/tools/HAT_list_runs.py --by erosion
    python scripts/hatteras_ms/tools/HAT_list_runs.py --only-stale

Grouped by topography (default), erosion preset or others; --stamp writes
BUILT_ON.txt into each run. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path


# Walk up until a directory holds data/hatteras_init
def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = REPO / "output" / "raw_runs"
import sys as _b3dsys
from pathlib import Path as _B3DP
_b3dsys.path.insert(0, str(next(_q for _q in _B3DP(__file__).resolve().parents
                                if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_topo_version as _b3d  # noqa: E402
DOMAIN_ROOT = _b3d.DOMAIN_ROOT
STAMP_NAME = "BUILT_ON.txt"
# -----------------------------------------------------------------------------


# The CURRENT dune-topo version of each product
def current_versions() -> dict:
    out = {}
    if not DOMAIN_ROOT.is_dir():
        return out
    for product in sorted(p.name for p in DOMAIN_ROOT.iterdir() if p.is_dir()):
        marker = DOMAIN_ROOT / product / "dune-topo" / "CURRENT"
        if marker.is_file():
            out[product] = marker.read_text(encoding="utf-8").strip()
    return out


# Every run's metadata
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


# The digest the most recent run of each (preset, period) carries
def newest_digest_per_preset(runs) -> dict:
    newest = {}
    for r in runs:
        key = (r["preset"], r["period"])
        if r["timestamp"] > newest.get(key, ("", ""))[0]:
            newest[key] = (r["timestamp"], r["digest"])
    return {k: d for k, (_, d) in newest.items()}


# The group a run falls in, for the chosen grouping
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


# The folder a run sits in, collapsed to what says its purpose
def parent_bucket(rel: str) -> str:
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


# Run: list the runs, grouped
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
