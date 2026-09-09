#!/usr/bin/env python3
r"""
HAT_version_pair_gif.py
==============================================================================
v2 beside v3 through time: the animations of the relocation comparison, but
with the two PANELS being the two dune-topo versions under ONE run scenario,
instead of the two relocation arms under one version. This is the view that
shows what the inserted and removed cells did (Hannah, 2026-09-09: "I want
to see how the inserted cells affected things").

WHAT IS DRAWN, per scenario (emergent: the modules decide; prescribed: the
recorded 1989/1999 relocations imposed) and per alongshore window:
    road_topography_<window>.gif   Barrier3D's own interior grids painted
        year by year, NC-12 on them, v2 left and v3 right, one colour scale
        and one year clock. Where v3 added rows the island is wider behind
        the road from year 0; where it removed them, narrower in front.
    road_relocation_<window>.gif   the dune line and the road as lines,
        landward-positive from each domain's year-0 dune line, a star where
        the module relocated, a ring where a prescribed move was applied.
One folder per place - the whole island, the two event blocks, and the two
reaches where the footprint is largest, Pea Island (GIS 78-87, rows added) and
the Avon-Tri-Village removals (GIS 62-68) - with `topography.gif` and
`dune-and-road.gif` in each.

Everything is read from the runs' saved state through the comparison
script's own loaders; nothing is re-run. The makers are the ones the
relocation comparison uses (cascade_pipeline.plotting.road_relocation_gif);
only the panel labels and the pairing differ.

WHERE  output/comparisons/relocation_1984_2004/v2_vs_v3/<scenario>/<place>/
       (its own folder: it is a cross-version comparison, not a set of one
       version, so it does not belong under v2/ or v3/)

USAGE
    python HAT_version_pair_gif.py                       # both scenarios
    python HAT_version_pair_gif.py --scenarios emergent
==============================================================================
"""
from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "hatteras_ms"))
import HAT_relocation_comparison as RC  # noqa: E402   loaders, windows, gif config
from cascade_pipeline.run_info import RunInfo  # noqa: E402
from cascade_pipeline.plotting.road_relocation_gif import (  # noqa: E402
    make_road_relocation_gif, make_topography_gif)

RAW = REPO / "output" / "raw_runs" / "version-pair"
OUT = REPO / "output" / "comparisons" / "relocation_1984_2004" / "v2_vs_v3"
SCENARIOS = {
    "emergent": ("HAT_1984_2004_calibBE_road_bdm_groin",
                 "full management, calibBE, groin; the modules decide on their own"),
    "prescribed": ("HAT_1984_2004_calibBE_road_reloc_bdm_groin",
                   "the same with the recorded 1989 and 1999 relocations prescribed"),
}
LABEL = {"v2": "v2 — the extraction (1996 dune, 2009 interior; today's setbacks)",
         "v3": "v3 — the 1984 reconstruction (rows behind / in front of NC-12; 1984 setbacks)"}
# ONE FOLDER PER PLACE, two files in each (Hannah, 2026-09-09: "organize the
# figures better"): a reader opens the reach they care about and finds both
# views of it side by side. The relocation-comparison windows first, then the
# two reaches where the footprint is largest.
#   (folder, title, topography window, line window)
WINDOWS = (
    ("1-island", "the whole island",
     RC.TOPO_WINDOWS[0][1:], RC.GIF_WINDOWS[0][1:]),
    ("2-event-1999_GIS9-14", "the 1999 relocation block (GIS 9-14)",
     RC.TOPO_WINDOWS[1][1:], RC.GIF_WINDOWS[1][1:]),
    ("3-event-1989_GIS84-87", "the 1989 relocation block (GIS 84-87)",
     RC.TOPO_WINDOWS[2][1:], RC.GIF_WINDOWS[2][1:]),
    ("4-rows-added_PeaIsland_GIS78-87", "Pea Island, where the footprint adds the most rows (GIS 78-87)",
     (76, 90), (76, 90)),
    ("5-rows-removed_AvonTriVillage_GIS62-68", "Avon to Tri-Village, where it removes the most (GIS 62-68)",
     (58, 72), (58, 72)),
)
FILES = {"topography": "topography.gif", "lines": "dune-and-road.gif"}


def write_readme(out_dir: Path, key: str, name: str, what: str, runs: dict) -> None:
    lines = [f"# v2 beside v3 \u2014 {key}", "", f"{name}: {what}.", "",
             f"Written {datetime.now():%Y-%m-%d %H:%M} by `HAT_version_pair_gif.py` from",
             f"`{runs['v2'].relative_to(REPO)}` and `{runs['v3'].relative_to(REPO)}`.",
             "Left panel v2 (the extraction), right panel v3 (the 1984 reconstruction), one year clock.",
             "", "One folder per place, two views in each:", "",
             "| folder | place | `topography.gif` | `dune-and-road.gif` |", "|---|---|---|---|"]
    for folder, title, (tlo, thi), (llo, lhi) in WINDOWS:
        lines.append(f"| `{folder}/` | {title} | Barrier3D's interior grids year by year, NC-12 on them, "
                     f"GIS {tlo}-{thi} | the dune line and the road as lines, landward-positive from each "
                     f"domain's year-0 dune line, a star where the module relocated, a ring where a prescribed "
                     f"move was applied, GIS {llo}-{lhi} |")
    (out_dir / "README.md").write_text(chr(10).join(lines) + chr(10), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="v2 beside v3 through time")
    ap.add_argument("--scenarios", default=",".join(SCENARIOS))
    ap.add_argument("--readme-only", action="store_true", help="rewrite the folder READMEs, render nothing")
    a = ap.parse_args()
    span = (RC.HATTERAS_FIRST_ROAD_DOMAIN, RC.HATTERAS_LAST_ROAD_DOMAIN)
    targets = RC.historical_targets(RC.START_YEAR, RC.END_YEAR)
    for key in [k.strip() for k in a.scenarios.split(",") if k.strip()]:
        name, what = SCENARIOS[key]
        runs = {v: RAW / v / "1984_2004" / "calibBE" / name for v in ("v2", "v3")}
        missing = [v for v, d in runs.items() if not d.is_dir()]
        if missing:
            print(f"{key}: no run for {missing} under {RAW} - skipped")
            continue
        out_dir = OUT / key
        out_dir.mkdir(parents=True, exist_ok=True)
        if a.readme_only:
            write_readme(out_dir, key, name, what, runs)
            print(f"{key}: README rewritten")
            continue
        print(f"{key}: {name}\n  loading v2 and v3 ...")
        casc = {v: RC.load_cascade(str(d)) for v, d in runs.items()}
        series = {v: RC.road_series(casc[v], RC.HATTERAS_DOMAINS, *span) for v in casc}
        shore = {v: RC.load_shoreline_matrix(str(runs[v])) for v in runs}
        info = {v: RunInfo(run_name=name, run_dir=str(runs[v]), start_year=RC.START_YEAR, end_year=RC.END_YEAR)
                for v in runs}
        back = {v: RC.back_barrier_matrix(casc[v]) for v in casc}
        written = []
        for folder, title, (tlo, thi), (llo, lhi) in WINDOWS:
            wdir = out_dir / folder
            wdir.mkdir(parents=True, exist_ok=True)
            if all(shore[v] is not None for v in shore):
                r = make_road_relocation_gif(
                    (shore["v2"], info["v2"]), (shore["v3"], info["v3"]), series["v2"], series["v3"],
                    llo, lhi, str(wdir / FILES["lines"]), back_a=back["v2"], back_b=back["v3"],
                    event_years=targets, gif_config=RC.GIF_CONFIG, label_a=LABEL["v2"], label_b=LABEL["v3"],
                    title=f"NC-12 and the dune line, v2 beside v3 \u2014 {title}")
                if r:
                    written.append(Path(r))
            r = make_topography_gif(
                casc["v2"], casc["v3"], series["v2"], series["v3"], tlo, thi, str(wdir / FILES["topography"]),
                RC.START_YEAR, event_years=targets, gif_config=RC.GIF_CONFIG,
                label_a=LABEL["v2"], label_b=LABEL["v3"],
                title=f"Hatteras topography and NC-12, v2 beside v3 \u2014 {title}",
                planform_note=RC.PLANFORM_NOTE)
            if r:
                written.append(Path(r))
        write_readme(out_dir, key, name, what, runs)
        print(f"  {len(written)} animations -> {out_dir}")


if __name__ == "__main__":
    main()
