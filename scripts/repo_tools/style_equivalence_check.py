"""
Prove a restyled script does the same thing as before: same code, only moved.

    python style_equivalence_check.py scripts/analyze_output/compare_runs/rate_windows.py
    python style_equivalence_check.py --ref bab312a3 scripts/analyze_output

Compares each file against its version at --ref (default HEAD). Passes only if,
with docstrings stripped: every function and class is identical; module-level
statements are the same multiset; functions and constant assignments may move,
but no name is read before it is bound and every other statement keeps its
order. Exits 1 on any failure. Standard library only.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())


# Remove every docstring, and bare module-level strings (no-ops), so trimming one is not a code change
def strip_docstrings(tree: ast.AST) -> ast.AST:
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                node.body = body[1:] or [ast.Pass()]
    tree.body = [s for s in tree.body if not (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant)
                                              and isinstance(s.value.value, str))]
    return tree


OTHERS = r"colleague_old_version|other_ms|from_lexi|from_roya"   # colleagues' code: inventory only
# -----------------------------------------------------------------------------


# One statement as comparable text, positions ignored
def dump(node: ast.AST) -> str:
    return ast.dump(node, include_attributes=False)


# Names a module-level statement binds, and names it reads when it runs
def binds(stmt: ast.stmt) -> set[str]:
    if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        return {stmt.name}
    if isinstance(stmt, (ast.Import, ast.ImportFrom)):
        return {(a.asname or a.name).split(".")[0] for a in stmt.names}
    out = set()
    for n in ast.walk(stmt):
        if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store):
            out.add(n.id)
    return out


# Names a statement reads when it runs (not inside nested function bodies)
def reads(stmt: ast.stmt) -> set[str]:
    if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        parts = list(stmt.decorator_list)
        if isinstance(stmt, ast.ClassDef):
            parts += stmt.bases + [k.value for k in stmt.keywords]
        else:
            parts += stmt.args.defaults + [d for d in stmt.args.kw_defaults if d]
        nodes = [n for p in parts for n in ast.walk(p)]
    else:
        nodes = []
        stack = [stmt]
        while stack:                      # don't descend into nested function bodies
            n = stack.pop()
            nodes.append(n)
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                continue
            stack.extend(ast.iter_child_nodes(n))
    return {n.id for n in nodes if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}


# May this statement move? Definitions and plain assignments may; nothing else
def movable(stmt: ast.stmt) -> bool:
    return isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef,
                             ast.Assign, ast.AnnAssign))


# Every way the new source could behave differently from the old
def compare(old_src: str, new_src: str) -> list[str]:
    problems = []
    old = strip_docstrings(ast.parse(old_src)).body
    new = strip_docstrings(ast.parse(new_src)).body

    if Counter(map(dump, old)) != Counter(map(dump, new)):
        gone = Counter(map(dump, old)) - Counter(map(dump, new))
        added = Counter(map(dump, new)) - Counter(map(dump, old))
        for label, c, src in (("removed or changed", gone, old), ("added or changed", added, new)):
            for s in src:
                if c[dump(s)]:
                    c[dump(s)] -= 1
                    problems.append(f"{label}: line {s.lineno}: {ast.unparse(s).splitlines()[0][:90]}")
        return problems

    anchored_old = [dump(s) for s in old if not movable(s)]
    anchored_new = [dump(s) for s in new if not movable(s)]
    if anchored_old != anchored_new:
        problems.append("statements that cannot move (imports, calls, if/for/with) changed order")

    module_names = set().union(*(binds(s) for s in old)) if old else set()
    for name in module_names:                              # rebindings keep their order
        seq_old = [dump(s) for s in old if name in binds(s)]
        seq_new = [dump(s) for s in new if name in binds(s)]
        if seq_old != seq_new:
            problems.append(f"'{name}' is bound in a different order")

    # Same names resolved as before; a call sees at least what it saw before
    def seen(stmts):
        out, bound = {}, set()
        for s in stmts:
            out.setdefault(dump(s), []).append((frozenset(reads(s) & bound), frozenset(bound)))
            bound |= binds(s)
        return out

    seen_old, seen_new = seen(old), seen(new)
    reach = reach_factory(old, module_names)
    for s in new:
        d = dump(s)
        (r_old, all_old), (r_new, all_new) = seen_old[d][0], seen_new[d][0]
        seen_old[d].pop(0), seen_new[d].pop(0)
        if r_old != r_new:
            problems.append(f"line {s.lineno}: resolves {sorted(r_old ^ r_new)} differently")
        calls = not isinstance(s, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) \
            and any(isinstance(n, ast.Call) for n in ast.walk(s))
        need = reach(s) & all_old
        if calls and not all_new >= need:
            problems.append(f"line {s.lineno}: a call now runs before {sorted(need - all_new)} exist")
    return problems


# Module globals a statement can reach: direct reads, local functions it calls, inputs of results it reads
def reach_factory(stmts, module_names):
    defs = {s.name: s for s in stmts
            if isinstance(s, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}
    body_reads = {n: {x.id for x in ast.walk(d) if isinstance(x, ast.Name)
                      and isinstance(x.ctx, ast.Load)} & module_names for n, d in defs.items()}
    carried: dict[str, set[str]] = {}
    for s in stmts:
        if not isinstance(s, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            got = _close(reads(s) & module_names, defs, body_reads, carried)
            for b in binds(s):
                carried[b] = got

    def reach(stmt):
        return _close(reads(stmt) & module_names, defs, body_reads, carried)
    return reach


# Transitive closure of names through function bodies and carried results
def _close(names, defs, body_reads, carried):
    seen, todo = set(), list(names)
    while todo:
        n = todo.pop()
        if n in seen:
            continue
        seen.add(n)
        todo += list(body_reads.get(n, ())) + list(carried.get(n, ()))
    return seen


# The file's text at a git ref, or None if it did not exist there
def at_ref(ref: str, path: Path) -> str | None:
    rel = path.resolve().relative_to(REPO).as_posix()
    r = subprocess.run(["git", "show", f"{ref}:{rel}"], cwd=REPO, capture_output=True)
    return r.stdout.decode("utf-8") if r.returncode == 0 else None


# Run: compare each file with its version at --ref, and check none went missing
def main() -> None:
    ap = argparse.ArgumentParser(description="Prove restyled scripts are unchanged in behaviour.")
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--ref", default="HEAD")
    a = ap.parse_args()

    files = [f for p in a.paths for f in ([p] if p.is_file() else sorted(p.rglob("*.py")))
             if "__pycache__" not in f.parts and not re.search(OTHERS, f.as_posix())]
    failed = 0

    # Inventory: every script at --ref must still exist at the same path
    for p in a.paths:
        if p.is_dir():
            rel = p.resolve().relative_to(REPO).as_posix()
            listed = subprocess.run(["git", "ls-tree", "-r", "--name-only", a.ref, "--", rel],
                                    cwd=REPO, capture_output=True, text=True).stdout.splitlines()
            for f in listed:
                if f.endswith(".py") and not (REPO / f).exists():
                    print(f"GONE  {f}")
                    failed += 1
    for f in files:
        old = at_ref(a.ref, f)
        if old is None:
            print(f"NEW   {f}")
            continue
        new = f.read_text(encoding="utf-8")
        problems = compare(old, new)
        doc_note = "  (uses __doc__: --help text may change)" if "__doc__" in new else ""
        print(f"{'FAIL' if problems else 'OK  '}  {f}{doc_note}")
        for p in problems:
            print(f"        {p}")
        failed += bool(problems)
    print(f"{len(files)} files, {failed} failed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
