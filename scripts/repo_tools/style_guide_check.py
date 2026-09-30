"""
Report where scripts depart from scripts/STYLE.md: header, author block, comments, CONFIG, banners.

    python scripts/repo_tools/style_guide_check.py scripts/analyze_output
    python scripts/repo_tools/style_guide_check.py path/to/one_script.py

Advisory: prints every departure per file and exits 0. Standard library only.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import ast
import re
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
MAX_HEADER_LINES = 12          # docstring lines above the author block
SKIP_DIRS = re.compile(r"(__pycache__|supersed|archive|legacy|other_ms|colleague_old_version|from_lexi|from_roya)", re.I)
AUTHOR_LINES = ("Author:  Hannah A. Henry", "Contact: hahenry@unc.edu", "Version: ")
BANNER = re.compile(r"^\s*#\s*(={5,}|-{5,}|─{3,}|═{3,})")
CONFIG_RULE = re.compile(r"^# --- CONFIG -+$|^# -{20,}$")
CODE_LIKE = re.compile(r"^\s*#\s*(\w+\s*=|dict\(|\)|\]|\}|[\w.]+\()")
# -----------------------------------------------------------------------------


# Every departure from the guide in one file
def check(path: Path) -> list[str]:
    src = path.read_text(encoding="utf-8")
    lines = src.splitlines()
    tree = ast.parse(src)
    out = []

    doc = ast.get_docstring(tree, clean=False) or ""
    if not doc:
        out.append("no module docstring")
    else:
        body = doc.strip("\n").split("\n")
        if not all(any(l.startswith(a) for l in body) for a in AUTHOR_LINES):
            out.append("author block missing or incomplete")
        above = next((i for i, l in enumerate(body) if l.startswith(("Author:", "Adapted from:"))), len(body))
        if above > MAX_HEADER_LINES:
            out.append(f"header is {above} lines above the author block (max {MAX_HEADER_LINES})")

    defs = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))]
    for d in defs:
        first = min([d.lineno] + [x.lineno for x in d.decorator_list])
        if first < 2 or not lines[first - 2].lstrip().startswith("#"):
            out.append(f"line {first}: {d.name}() has no one-line comment above it")
        ds = ast.get_docstring(d, clean=False)
        if ds and len(d.body) > 1:
            out.append(f"line {d.lineno}: {d.name}() still has a docstring")

    upper = [n for n in tree.body if isinstance(n, ast.Assign) and any(
        isinstance(t, ast.Name) and t.id.isupper() for t in n.targets)]
    if upper and not any(l.startswith("# --- CONFIG") for l in lines):
        out.append(f"{len(upper)} module constants but no CONFIG block")

    run = []
    for i, l in enumerate(lines + [""], 1):
        if BANNER.match(l) and not CONFIG_RULE.match(l.strip()):
            out.append(f"line {i}: banner")
        if l.lstrip().startswith("#") and not CODE_LIKE.match(l) and not CONFIG_RULE.match(l.strip()):
            run.append(i)
            continue
        if len(run) > 1:
            out.append(f"lines {run[0]}-{run[-1]}: {len(run)}-line comment block")
        run = []
    return out


# Run: check every script under the given paths and summarise
def main() -> None:
    ap = argparse.ArgumentParser(description="Check scripts against scripts/STYLE.md.")
    ap.add_argument("paths", nargs="+", type=Path)
    a = ap.parse_args()
    files = [f for p in a.paths for f in ([p] if p.is_file() else sorted(p.rglob("*.py")))
             if not SKIP_DIRS.search(f.as_posix())]
    clean = 0
    for f in files:
        problems = check(f)
        print(f"{'OK  ' if not problems else 'NOTE'}  {f}")
        for p in problems:
            print(f"        {p}")
        clean += not problems
    print(f"{len(files)} files, {clean} follow the guide, {len(files) - clean} with notes")


if __name__ == "__main__":
    main()
