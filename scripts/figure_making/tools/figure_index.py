"""
figure_index.py
==============================================================================
Writes `output/figures/README.md`, the map of the figures: the layout, a
"where do I find..." table, and one table per folder listing every figure,
what it shows, the script that draws it and the day it was drawn.

    python scripts/figure_making/tools/figure_index.py

Re-run after adding or redrawing a figure; `regenerate_all_figures.py` runs it
at the end.

WHY IT IS GENERATED
    A hand-kept index of a folder ~25 scripts write into is stale the day after
    it is written, and a wrong index is worse than none. Everything in the
    table is recorded beside the figures:
      * "shows" is the first sentence of the figure's entry in its folder's
        `supporting/CAPTIONS.md`;
      * "drawn by" comes from `supporting/producers.json`, which
        `regenerate_all_figures.py` writes by noting which files each producer
        touched; for a figure redrawn by hand since, it falls back to searching
        the scripts tree for the figure's file name;
      * "drawn" is the PNG's modification date, so a stale figure shows.
    A dash in "shows" or "drawn by" is a real finding: nobody wrote down what
    the figure shows, or nothing can reproduce it.

THE LAYOUT (Hannah, 2026-09-29): numbered in the paper's order, see LAYOUT.
    The folder names come from hat_figure_style.FIGURE_SUBJECTS / INPUT_STEPS.
==============================================================================
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_figure_style as _hs  # noqa: E402

FIGURES = _hs.FIGURES_ROOT
SCRIPTS = REPO / "scripts"
PRODUCERS = FIGURES / "supporting" / "producers.json"

# Every folder that holds figures, in reading order, with one line saying what
# it answers. A folder on disk that is missing here is listed under "other" so
# it cannot hide.
LAYOUT = {
    "1-site": "Where the reach is, how the 90 model domains tile it, and what one domain is.",
    "2-observations": "What was measured: CoastSat shoreline, the digitised dune lines, the mean shoreline on imagery.",
    "2-observations/shoreline": "CoastSat shoreline-change rates by domain and period.",
    "2-observations/duneline": "The digitised dune lines: where they sit, how they moved, beach width, distance to NC-12.",
    "2-observations/shoreline_vs_duneline": "CoastSat shoreline change against dune-line change, side by side.",
    "2-observations/mean_shoreline": "The CoastSat mean shoreline over each start DEM's window, drawn on aerial imagery.",
    "2-observations/mean_shoreline/line_and_band": "The mean line and its spread.",
    "2-observations/mean_shoreline/with_domains": "The same with the model domain boxes.",
    "2-observations/mean_shoreline/with_positions": "The same with every satellite position in the window.",
    "3-model-inputs": "How each model input is built from its source data, one folder per data/hatteras_init/ step.",
    "3-model-inputs/0-elevation": "Which survey supplies each domain's topography, and the resample to 10 m.",
    "3-model-inputs/1-domains": "How one domain's Barrier3D arrays are cut from the DEM.",
    "3-model-inputs/1-domains/initial_island": "The whole island as the model starts it, per start year.",
    "3-model-inputs/2-brie-offset": "How the BRIE shoreline offset (the island's planform) is built.",
    "3-model-inputs/3-forcing": "Storms and sea level: how the storm series is built and the forcing record per chain.",
    "3-model-inputs/4-management": "NC-12 setbacks, the management timeline and rules, and where the reach is managed.",
    "3-model-inputs/5-observed-target": "How the observed shoreline-change target is built from CoastSat.",
    "3-model-inputs/7-source-sink": "The source/sink (background erosion) field and the end-domain solve.",
    "4-model-mechanics": "How the models work, drawn from a finished Hatteras run.",
    "4-model-mechanics/barrier3d": "Barrier3D: one domain through a storm year and over a window.",
    "4-model-mechanics/brie": "BRIE: alongshore diffusion, the wave asymmetry, and the domain order.",
    "4-model-mechanics/cascade": "CASCADE: the coupling loop, the island as held, the shoreline split, the management modules.",
    "4-model-mechanics/storm_routing": "How Barrier3D routes overwash, replayed storm by storm.",
    "5-results": "What the model gives: hindcasts against CoastSat, every management scenario, worked examples.",
    "style": "The house style sheet every figure is drawn under.",
}

# "Where do I find ...": the questions people actually ask, pointed at a path
# under output/figures/. Each path is checked; a missing one is flagged.
FAQ = [
    ("Where is the study area / a map of the domains?", "1-site/study_area.png"),
    ("Which way do the domain numbers run?", "4-model-mechanics/brie/brie_domain_orientation.png"),
    ("What does one domain look like as the model reads it?", "1-site/domain_grid.png"),
    ("What are the observed shoreline-change rates?", "2-observations/shoreline/observed_rates.png"),
    ("How do the two calibration periods compare (observed)?", "2-observations/shoreline/coastsat_calibration_periods.png"),
    ("How has the dune line moved?", "2-observations/duneline/duneline_positions_overview.png"),
    ("Shoreline change vs dune-line change?", "2-observations/shoreline_vs_duneline/coastsat_endpoint_vs_duneline_1996_2010_2024_stacked.png"),
    ("What island does the model start from?", "3-model-inputs/1-domains/initial_island/1996/classes/island_1996.png"),
    ("How is the island planform (BRIE offset) built?", "3-model-inputs/2-brie-offset/offset_build_1996.png"),
    ("Which storms reach the model?", "3-model-inputs/3-forcing/storm_events_by_duration.png"),
    ("When and where was the reach managed?", "3-model-inputs/4-management/timeline_1996_2024.png"),
    ("How is the observed target built?", "3-model-inputs/5-observed-target/observed_target_1996_2010.png"),
    ("How does CASCADE couple Barrier3D and BRIE?", "4-model-mechanics/cascade/cascade_coupling_loop.png"),
    ("What does BRIE's wave asymmetry do?", "4-model-mechanics/brie/brie_asymmetry_explained.png"),
    ("How does overwash get routed in a storm?", "4-model-mechanics/storm_routing/storm_routing_hours.png"),
    ("How well does the hindcast match CoastSat?", "5-results/hindcast_edgeBE.png"),
    ("What does each management scenario do?", "5-results/scenario_grid.png"),
]
SKIP_PARTS = {"supporting", "talk"}


def first_sentence(text: str) -> str:
    """The caption's first sentence, which is written to stand alone."""
    text = " ".join(text.split())
    m = re.search(r"(.+?\.)(?:\s|$)", text)
    s = (m.group(1) if m else text)
    return s if len(s) <= 220 else s[:217].rstrip() + "…"


def captions(folder: Path) -> dict[str, str]:
    """{figure file name: caption} from the folder's supporting/CAPTIONS.md."""
    md = folder / "supporting" / "CAPTIONS.md"
    if not md.is_file():
        return {}
    out = {}
    for m in re.finditer(r"^\*\*`([^`]+)`\.\*\*\s*(.*?)(?=\n\*\*`|\Z)",
                         md.read_text(encoding="utf-8"), re.S | re.M):
        out[m.group(1)] = first_sentence(m.group(2))
    return out


def source_captions() -> dict[str, str]:
    """{file name: caption} from every CAPTIONS.md under data/ and
    output/comparisons/. A figure PUBLISHED as a copy (the mean-shoreline
    images, copied from data/hatteras_init/5-scr/) has its caption beside the
    original, not beside the copy; this finds it by file name."""
    out: dict[str, str] = {}
    for root in (REPO / "data", REPO / "output" / "comparisons"):
        for md in root.rglob("CAPTIONS.md"):
            if "archive" in md.parts or any(p.startswith("superseded") for p in md.parts):
                continue
            for name, text in captions(md.parent.parent if md.parent.name == "supporting"
                                       else md.parent).items():
                out.setdefault(name, text)
    return out


def recorded_producers() -> dict[str, str]:
    """{path under output/figures: script} as regenerate_all_figures.py saw it."""
    if PRODUCERS.is_file():
        return json.loads(PRODUCERS.read_text(encoding="utf-8"))
    return {}


def searched_producers() -> dict[str, str]:
    """{figure stem: script}, by searching the scripts tree for the file name.
    The fallback for a figure redrawn by hand after the last full regeneration."""
    hits: dict[str, str] = {}
    for py in SCRIPTS.rglob("*.py"):
        if "__pycache__" in py.parts or "superseded" in str(py):
            continue
        text = py.read_text(encoding="utf-8", errors="ignore")
        rel = py.relative_to(REPO).as_posix()
        for m in re.finditer(r'["\']([a-z0-9_]+)\.png["\']', text):
            hits.setdefault(m.group(1), rel)
    return hits


def figure_dirs() -> list[Path]:
    """Every folder under output/figures holding a PNG, talk/ and supporting/ aside."""
    dirs = {p.parent for p in FIGURES.rglob("*.png")}
    return sorted(d for d in dirs
                  if not SKIP_PARTS.intersection(d.relative_to(FIGURES).parts))


def tree_lines() -> list[str]:
    """The layout as an indented tree with figure counts."""
    lines = ["```", "output/figures/"]
    for key, blurb in LAYOUT.items():
        d = FIGURES / key
        if not d.is_dir():
            continue
        n = sum(1 for p in d.rglob("*.png") if "supporting" not in p.parts)
        depth = key.count("/")
        name = key.split("/")[-1] + "/"
        short = blurb.split(":")[0].split(".")[0]
        lines.append(f"{'  ' * (depth + 1)}{name:<{30 - 2 * depth}} {n:>3}  {short}")
    talk = FIGURES / "talk"
    if talk.is_dir():
        n = sum(1 for p in talk.rglob("*.png") if "supporting" not in p.parts)
        lines.append(f"  {'talk/':<30} {n:>3}  projector versions, same paths as above")
    lines.append("```")
    return lines


def main() -> Path:
    recorded = recorded_producers()
    searched = searched_producers()
    elsewhere = source_captions()
    lines = [
        "# output/figures — the map",
        "",
        "Every figure this project publishes for a manuscript, poster or talk, in the",
        "paper's order: where → what was observed → what the model is fed → how it works",
        "→ what it gives. **Generated** by `scripts/figure_making/tools/figure_index.py`;",
        "do not edit by hand.",
        "",
        "**Bring everything up to date:** `python scripts/figure_making/tools/regenerate_all_figures.py`",
        "(one folder: `--only 5-results`). Each figure's PDF and caption sit in the",
        "`supporting/` folder beside it; `talk/` holds projector versions at the same paths.",
        "Retired figures go to `output/archive/<date>_<what>/`, never back in here.",
        "",
    ]
    lines += tree_lines() + [""]

    lines += ["## Where do I find …", "", "| question | figure |", "|---|---|"]
    for q, rel in FAQ:
        mark = "" if (FIGURES / rel).is_file() else " **(missing)**"
        lines.append(f"| {q} | [`{rel}`]({rel}){mark} |")
    lines.append("")

    known = set(LAYOUT)
    dirs = figure_dirs()
    order = {k: i for i, k in enumerate(LAYOUT)}
    dirs.sort(key=lambda d: (order.get(d.relative_to(FIGURES).as_posix(), 999),
                             d.relative_to(FIGURES).as_posix()))
    total = 0
    for d in dirs:
        key = d.relative_to(FIGURES).as_posix()
        # initial_island/<year>/<scheme> fold into their parent's heading
        head = next((k for k in sorted(known, key=len, reverse=True)
                     if key == k or key.startswith(k + "/")), None)
        caps = captions(d)
        pngs = sorted(d.glob("*.png"))
        if not pngs:
            continue
        title = key if head == key else f"{key}"
        blurb = LAYOUT.get(key) or (LAYOUT.get(head, "") if head else "**Not in the layout** — add it to figure_index.LAYOUT.")
        lines += [f"### `{title}/`", "", blurb, "",
                  "| figure | shows | drawn by | drawn |", "|---|---|---|---|"]
        for png in pngs:
            rel = png.relative_to(FIGURES).as_posix()
            script = recorded.get(rel) or searched.get(png.stem, "—")
            when = datetime.fromtimestamp(png.stat().st_mtime).strftime("%Y-%m-%d")
            lines.append(f"| [`{png.name}`]({rel}) | {caps.get(png.name) or elsewhere.get(png.name, '—')} | "
                         f"`{script}` | {when} |")
            total += 1
        lines.append("")

    out = FIGURES / "README.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out} ({total} figures listed, talk/ aside)")
    return out


if __name__ == "__main__":
    main()
