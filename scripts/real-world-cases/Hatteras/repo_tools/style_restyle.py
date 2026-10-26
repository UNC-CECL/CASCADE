"""
Reassemble a script into the scripts/STYLE.md layout: header, imports, CONFIG, functions, guard.

    from style_restyle import restyle, current_order, strip_function_docstrings
    restyle(path, header, order, comments, config_comments, replace=[...])

A library for the style sweep, driven from a per-script plan. It only moves
code and rewrites comments and docstrings; always prove the result with
style_equivalence_check.py before committing it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import ast
import re
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
NL = chr(10)
RULE = "# " + "-" * 77
CONFIG_HEAD = "# --- CONFIG " + "-" * 66
# -----------------------------------------------------------------------------


# Contiguous comment lines directly above line `start` (1-based), after `prev_end`
def _gap_comment(lines, start, prev_end):
    out, i = [], start - 2
    while i >= prev_end and lines[i].lstrip().startswith("#") and not re.match(r"#\s*[-=]{5,}", lines[i].strip()):
        out.insert(0, lines[i])
        i -= 1
    return out


# Top-level def/class names in their current order, to start a plan from
def current_order(path):
    body = ast.parse(Path(path).read_text(encoding="utf-8")).body
    return [s.name for s in body if isinstance(s, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))]


# Rewrite `path`: header, kept preamble, CONFIG, defs in `order`, `late` assignments, `replace` pairs
def restyle(path, header, order, comments, config_comments=None, late=(), replace=()):
    p = Path(path)
    src = p.read_text(encoding="utf-8")
    eol = "\r\n" if "\r\n" in src else "\n"
    src = src.replace("\r\n", "\n")
    lines = src.split("\n")
    body = ast.parse(src).body
    shebang = lines[0] + "\n" if lines[0].startswith("#!") else ""

    # Imports, sys.path set-up, private globals and guarded set-up belong to the preamble
    def is_pre(s):
        if isinstance(s, (ast.Import, ast.ImportFrom)):
            return True
        if isinstance(s, ast.Expr) and isinstance(s.value, ast.Call):
            return True
        if isinstance(s, ast.Assign) and len(s.targets) == 1 and isinstance(s.targets[0], ast.Name) \
                and s.targets[0].id.startswith("_") and s.targets[0].id not in comments:
            return True
        if isinstance(s, (ast.If, ast.Try)) and not (isinstance(s, ast.If) and "__name__" in ast.unparse(s.test)):
            return True
        return False

    i = 1 if (body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant)) else 0
    doc_end = body[0].end_lineno if i else 0
    j = i
    while j < len(body) and is_pre(body[j]):
        j += 1
    pre_end = body[j - 1].end_lineno if j > i else doc_end
    preamble = NL.join(lines[doc_end:pre_end]).strip(NL)

    # A statement's source, decorators included
    def seg(s):
        first = min([s.lineno] + [d.lineno for d in getattr(s, "decorator_list", [])])
        return NL.join(lines[first - 1:s.end_lineno]), first

    # Names a statement binds
    def binds(s):
        if isinstance(s, (ast.Import, ast.ImportFrom)):
            return {(a.asname or a.name).split(".")[0] for a in s.names}
        return {n.id for n in ast.walk(s) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store)}

    module_bound = set().union(*(binds(s) for s in body)) if body else set()
    ready = set().union(*(binds(s) for s in body[:j])) if j else set()
    lifting = [True]

    # Lift into the preamble only while order is kept and the inputs already exist there
    def liftable(s):
        if not is_pre(s):
            return False
        loads = {n.id for n in ast.walk(s) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
        if lifting[0] and (loads & module_bound) <= ready:
            return True
        lifting[0] = False
        return False

    defs, config, late_out, guard, pre_extra = {}, [], [], None, []
    prev_end = pre_end
    for s in body[j:]:
        text, first = seg(s)
        above = _gap_comment(lines, first, prev_end)
        prev_end = s.end_lineno
        if isinstance(s, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            defs[s.name] = _drop_docstring(s, text)
        elif isinstance(s, ast.If) and "__name__" in ast.unparse(s.test):
            guard = text
        elif liftable(s):
            pre_extra.append(NL.join(([above[-1]] if len(above) == 1 else []) + [text]))
            ready |= binds(s)
        else:
            names = [t.id for t in (s.targets if isinstance(s, ast.Assign) else [s.target])
                     if isinstance(t, ast.Name)] if isinstance(s, (ast.Assign, ast.AnnAssign)) else []
            key = names[0] if names else None
            if config_comments is not None and key in config_comments:
                above = [f"# {config_comments[key]}"] if config_comments[key] else []
            elif len(above) > 1:
                above = []
            (late_out if key in late else config).append(NL.join(above + [text]))

    missing, extra = set(defs) - set(order), set(order) - set(defs)
    if missing or extra:
        raise SystemExit(f"order mismatch: missing {sorted(missing)}, unknown {sorted(extra)}")

    parts = [shebang + '"""' + NL + header.strip(NL) + NL + '"""', NL.join([preamble] + pre_extra)]
    if config:
        parts.append(CONFIG_HEAD + NL + NL.join(config) + NL + RULE)
    for name in order:
        parts.append(f"# {comments[name]}{NL}{defs[name]}" if comments.get(name) else defs[name])
    parts += late_out
    if guard:
        parts.append(guard)
    out = (NL * 3).join(x for x in parts if x) + NL
    out = out.replace('"""' + NL * 3, '"""' + NL, 1)
    for a, b in replace:
        if out.count(a) != 1:
            raise SystemExit(f"replace target found {out.count(a)} times: {a[:60]!r}")
        out = out.replace(a, b)
    p.write_text(out.replace(NL, eol), encoding="utf-8")


# A def's source with its docstring removed (kept if it is the whole body)
def _drop_docstring(node, text):
    b = node.body
    if b and isinstance(b[0], ast.Expr) and isinstance(b[0].value, ast.Constant) \
            and isinstance(b[0].value.value, str) and len(b) > 1:
        tl = text.split(NL)
        first = min([node.lineno] + [d.lineno for d in getattr(node, "decorator_list", [])])
        off = node.lineno - first
        a, z = b[0].lineno - node.lineno, b[0].end_lineno - node.lineno
        del tl[a + off:z + off + 1]
        text = NL.join(tl)
    return text


# Remove every def/class docstring whose body has other statements; returns how many
def strip_function_docstrings(path):
    p = Path(path)
    src = p.read_text(encoding="utf-8")
    eol = "\r\n" if "\r\n" in src else "\n"
    lines = src.replace("\r\n", NL).split(NL)
    kill = []
    for node in ast.walk(ast.parse(NL.join(lines))):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            b = node.body
            if len(b) > 1 and isinstance(b[0], ast.Expr) and isinstance(b[0].value, ast.Constant) \
                    and isinstance(b[0].value.value, str):
                kill.append((b[0].lineno, b[0].end_lineno))
    for a, z in sorted(kill, reverse=True):
        del lines[a - 1:z]
    p.write_text(eol.join(lines), encoding="utf-8")
    return len(kill)
