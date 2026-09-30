"""
Check that the explanation removed from each restyled script landed in its folder README.

    python scripts/repo_tools/style_readme_coverage.py scripts/analyze_output --ref bab312a3

For every script: the docstring and comment prose it had at --ref but no longer
has, reduced to content words, and the share of them found in the script's
README section (the `### <path>` heading naming it) or still in the script.
Lists the missing words to read. Advisory, exits 0. Standard library only.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import ast
import io
import re
import subprocess
import tokenize
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
MIN_WORD = 4
OTHERS = r"colleague_old_version|other_ms|from_lexi|from_roya"   # colleagues' code, never restyled
SKIP = {"this", "that", "with", "from", "have", "been", "were", "which", "when", "than", "then",
        "they", "them", "their", "there", "here", "into", "onto", "also", "only", "each", "every",
        "over", "under", "does", "done", "used", "uses", "same", "what", "where", "while", "would",
        "could", "should", "because", "before", "after", "about", "above", "below", "these", "those",
        "since", "until", "still", "just", "more", "most", "less", "such", "other", "both", "either",
        "none", "never", "always", "once", "twice", "must", "will", "shall", "being", "very"}
# -----------------------------------------------------------------------------


# The prose of a source: docstrings and comments, one entry per line
def prose(src: str) -> list[str]:
    out = []
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            d = ast.get_docstring(node, clean=True)
            if d:
                out += d.splitlines()
    for node in ast.parse(src).body[1:]:          # a header stranded below an inserted line
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) \
                and isinstance(node.value.value, str):
            out += node.value.value.splitlines()
    for tok in tokenize.generate_tokens(io.StringIO(src).readline):
        if tok.type == tokenize.COMMENT:
            out.append(tok.string.lstrip("#").strip())
    return [" ".join(l.split()) for l in out if l.strip()]


# Content words of some text, lower case
def words(text: str) -> set[str]:
    return {w for w in re.findall(r"[a-z][a-z0-9_]+", text.lower())
            if len(w) >= MIN_WORD and w not in SKIP}


# The README section for one script: from its ### heading to the next
def section(readme: str, rel: str) -> str:
    parts = re.split(r"(?m)^### ", readme)
    return "\n".join(p for p in parts if rel in p.split("\n", 1)[0])


# The file's text at a git ref, or None
def at_ref(ref: str, path: Path) -> str | None:
    r = subprocess.run(["git", "show", f"{ref}:{path.resolve().relative_to(REPO).as_posix()}"],
                       cwd=REPO, capture_output=True)
    return r.stdout.decode("utf-8") if r.returncode == 0 else None


# Run: per script, removed prose vs its README section
def main() -> None:
    ap = argparse.ArgumentParser(description="Did the removed explanations reach the README?")
    ap.add_argument("folder", type=Path)
    ap.add_argument("--ref", default="HEAD")
    a = ap.parse_args()

    readme = (a.folder / "README.md").read_text(encoding="utf-8")
    for f in sorted(a.folder.rglob("*.py")):
        if "__pycache__" in f.parts or re.search(r"supersed|archive|" + OTHERS, f.as_posix(), re.I):
            continue
        old = at_ref(a.ref, f)
        if old is None:
            continue
        new = f.read_text(encoding="utf-8")
        try:
            kept = set(prose(new))
            removed = " ".join(l for l in prose(old) if l not in kept)
        except (SyntaxError, tokenize.TokenError) as e:
            print(f"  ??   {f}: does not parse ({type(e).__name__}), skipped")
            continue
        gone = words(removed)
        if not gone:
            print(f"  --   {f}: nothing removed")
            continue
        rel = f.relative_to(a.folder).as_posix()
        sec = section(readme, rel)
        found = words(sec) | words(new)
        missing = sorted(gone - found)
        share = 100.0 * (1 - len(missing) / len(gone))
        flag = "OK  " if share >= 90 else "READ"
        print(f"{flag} {share:5.1f}%  {f}  ({len(gone)} words removed"
              f"{', no README section' if not sec else ''})")
        if missing:
            print(f"        not found: {', '.join(missing[:25])}{' ...' if len(missing) > 25 else ''}")


if __name__ == "__main__":
    main()
