"""
Create an empty research project laid out the way README.md describes.

    python project_layout_template.py my_project
    python project_layout_template.py my_project --stages observations model-inputs runs

Makes scripts/, data/ and output/ with matching numbered stages, a README.md in
every folder, and scripts/lib/paths.py as the one place locations are decided.
Never overwrites a file that already exists. Needs only the standard library.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
DEFAULT_STAGES = ["observations", "processing", "analysis"]
# -----------------------------------------------------------------------------

README_STUB = """# {title}

{purpose}

## What is here

## What made it

## What not to trust
"""

PATHS_PY = '''"""
Where everything lives -- the one place. Scripts import from here and never
type a path, so moving a folder is one edit.
"""

from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
DATA = ROOT / "data"
OUTPUT = ROOT / "output"
FIGURES = OUTPUT / "figures"
LOGS = OUTPUT / "logs"


# A stage folder under data/, by name or number; unknown names fail loudly
def stage(name: str) -> Path:
    for d in sorted(DATA.iterdir()):
        if d.is_dir() and (d.name == name or d.name.split("-", 1)[-1] == name
                           or d.name.split("-", 1)[0] == str(name)):
            return d
    known = ", ".join(d.name for d in sorted(DATA.iterdir()) if d.is_dir())
    raise FileNotFoundError(f"no data stage {name!r}; known: {known}")
'''

PYPROJECT = """[project]
name = "{name}"
version = "0.1.0"
"""


# Write a file only if it isn't there yet
def write_new(path: Path, text: str, made: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_text(text, encoding="utf-8")
        made.append(path)


# The tree: (folder, purpose) pairs, stages numbered and mirrored
def layout(stages: list[str]) -> list[tuple[str, str]]:
    numbered = [f"{i}-{s}" for i, s in enumerate(stages, 1)]
    rows = [
        (".", "The project. Code in scripts/, what it reads in data/, what it "
              "makes in output/."),
        ("scripts", "All code. Folders and this README only at this level."),
        ("scripts/lib", "Modules the scripts share; produces nothing."),
        ("data", "Everything the analysis reads, by stage. Could be archived "
                 "and handed to someone on its own."),
        ("output", "Everything the code produces. Safe to delete and rebuild."),
        ("output/figures", "Finished figures: PNGs at the top, PDFs and "
                           "CAPTIONS.md under supporting/."),
        ("output/logs", "Run logs."),
        ("output/archive", "Retired products, as YYYY-MM-DD_<what>/."),
    ]
    for n in numbered:
        rows.append((f"scripts/{n}", f"Code for stage {n}; its data is data/{n}/."))
        rows.append((f"data/{n}", f"Data for stage {n}; made by scripts/{n}/."))
    return rows


# Run: build the tree, the resolver and the root marker
def main() -> None:
    ap = argparse.ArgumentParser(description="Create an empty project layout.")
    ap.add_argument("root", type=Path)
    ap.add_argument("--stages", nargs="+", default=DEFAULT_STAGES)
    a = ap.parse_args()

    made: list[Path] = []
    for folder, purpose in layout(a.stages):
        d = a.root / folder
        title = a.root.name if folder == "." else folder
        write_new(d / "README.md", README_STUB.format(title=title, purpose=purpose), made)
    write_new(a.root / "scripts/lib/paths.py", PATHS_PY, made)
    write_new(a.root / "pyproject.toml", PYPROJECT.format(name=a.root.name), made)

    print(f"{len(made)} new files under {a.root.resolve()}")
    for p in sorted(a.root.rglob("*")):
        if p.is_dir():
            print("  " + p.relative_to(a.root).as_posix() + "/")


if __name__ == "__main__":
    main()
