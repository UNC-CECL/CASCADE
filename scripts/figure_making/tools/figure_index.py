"""
figure_index.py
==============================================================================
Writes `output/figures/README.md`: one table per subject folder, listing every
figure, what it shows and the script that draws it.

    python scripts/figure_making/tools/figure_index.py

WHY IT IS GENERATED
    A hand-kept index of a folder that eleven scripts write into is stale the
    day after it is written, and a wrong index is worse than none. Everything
    in the table is already recorded beside the figures: the "shows" column is
    the first sentence of the figure's entry in `supporting/CAPTIONS.md`, and
    the "drawn by" column is found by searching the scripts tree for the
    figure's own file name. Re-run it after adding a figure.

    A figure with no caption entry, or none that any script names, is listed
    with a dash and is a real finding: it means nobody wrote down what it
    shows, or it is an orphan nothing can reproduce.

THE LAYOUT IT DOCUMENTS
    output/figures/ is organised by SUBJECT, not by the script that drew the
    figure -- the same rule scripts/figure_making/ follows (ORGANIZATION.md
    rule 1). `talk/` mirrors the subjects with the projector versions, and a
    retired figure goes to a dated `superseded_*` folder with a WHY.md
    (rule 4), never deleted.
==============================================================================
"""
from __future__ import annotations

import re
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
import sys as _hssys
from pathlib import Path as _HSP
_hssys.path.insert(0, str(next(_q for _q in _HSP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_figure_style as _hs  # noqa: E402
FIGURES = _hs.FIGURES_ROOT
SCRIPTS = REPO / "scripts"

# the order the subjects read in, and one line saying what each holds
SUBJECTS = {
    "site": "Where the reach is, how the 90 domains tile it, and what one domain is",
    "forcing": "The record the model eats: storms, sea level, and the management timeline",
    "management": "NC-12 and beach nourishment: where, when, and under what rules",
    "shoreline": "Shoreline change, observed and modelled",
    "initialization": "The island as the model starts it",
}
SKIP_DIRS = {"supporting", "talk"}


def first_sentence(text: str) -> str:
    """The caption's first sentence, which is written to stand alone."""
    text = " ".join(text.split())
    m = re.search(r"(.+?\.)(?:\s|$)", text)
    return (m.group(1) if m else text)[:300]


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


def drawn_by() -> dict[str, str]:
    """{figure stem: script path relative to the repository}, by searching the
    scripts tree for the figure's file name or stem."""
    hits: dict[str, str] = {}
    for py in SCRIPTS.rglob("*.py"):
        if "__pycache__" in py.parts or "superseded" in str(py):
            continue
        try:
            text = py.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        rel = py.relative_to(REPO).as_posix()
        for m in re.finditer(r'["\']([a-z0-9_]+)\.png["\']', text):
            hits.setdefault(m.group(1), rel)
        # the site-figure script names its figures in a SUBJECT table, not in
        # a path literal
        for m in re.finditer(r'^\s*"([a-z0-9_]+)":\s*"(?:site|forcing|management|shoreline)"',
                             text, re.M):
            hits.setdefault(m.group(1), rel)
    return hits


def by_folder() -> dict[str, str]:
    """{subject: script} for a script that builds its file name at run time
    (an f-string), where searching for the name cannot find it. Such a script
    still names the FOLDER it writes into, which is enough to attribute the
    figures in it."""
    out: dict[str, str] = {}
    for py in SCRIPTS.rglob("*.py"):
        if "__pycache__" in py.parts or "superseded" in str(py):
            continue
        try:
            text = py.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for subject in SUBJECTS:
            if f'"figures" / "{subject}"' in text or f'"figures", "{subject}"' in text:
                out.setdefault(subject, py.relative_to(REPO).as_posix())
    return out


def main() -> Path:
    scripts_by_stem = drawn_by()
    scripts_by_folder = by_folder()
    lines = [
        "# figures — by subject",
        "",
        "Every figure this project draws for a manuscript, a poster or a talk.",
        "**Generated** by `scripts/figure_making/tools/figure_index.py`; re-run it",
        "after adding a figure rather than editing this file.",
        "",
        "One folder per subject, not per script (ORGANIZATION.md rule 1). A figure's",
        "PNG sits at the top of its subject folder; its PDF and its caption go under",
        "`supporting/` (scripts/figure_making/STYLE.md). `talk/` mirrors the subjects with",
        "projector versions: water almost white, a point more type, heavier lines.",
        "A retired figure goes to a dated `superseded_*` folder with a note, never",
        "to the bin (rule 4).",
        "",
    ]
    for subject, blurb in SUBJECTS.items():
        folder = FIGURES / subject
        if not folder.is_dir():
            continue
        caps = captions(folder)
        pngs = sorted(p for p in folder.glob("*.png"))
        lines += [f"## {subject}", "", blurb, "",
                  "| figure | shows | drawn by |", "|---|---|---|"]
        for png in pngs:
            shows = caps.get(png.name, "—")
            script = scripts_by_stem.get(png.stem) or scripts_by_folder.get(subject, "—")
            lines.append(f"| `{png.name}` | {shows} | `{script}` |")
        lines.append("")

    talk = FIGURES / "talk"
    if talk.is_dir():
        n = len(list(talk.rglob("*.png")))
        lines += ["## talk", "",
                  f"The same {n} figures for a projector, under `talk/<subject>/`. Drawn by the",
                  "same scripts with `--talk`.", ""]
    for sup in sorted(FIGURES.glob("superseded_*")):
        n = len(list(sup.rglob("*.png")))
        lines += [f"## {sup.name}", "",
                  f"{n} retired figures; see `{sup.name}/WHY.md`.", ""]

    out = FIGURES / "README.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out} ({sum(1 for l in lines if l.startswith('| `'))} figures listed)")
    return out


if __name__ == "__main__":
    main()
