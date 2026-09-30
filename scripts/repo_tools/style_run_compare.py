"""
Run scripts before and after a restyle in a scratch copy, and compare what they wrote.

    python style_run_compare.py --root C:/Users/hanna/sw_verify/root \
        --before <scripts dir at the old commit> --after <restyled scripts dir> \
        --spec runs.json --out C:/Users/hanna/sw_verify/results

runs.json is a list of {"label", "script", "args": [...], "env": {...}, "set": {...}},
the script path relative to scripts/. Each run executes in --root with the old
then the new scripts/ in place; every file written under data/ and output/ is
kept and compared. "set" replaces top-level assignments in memory (e.g. an
interactive MODE), the same way on both sides; the file is never edited.
Never touches the real repository. Needs numpy and Pillow.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# --- CONFIG ------------------------------------------------------------------
WATCHED = ("data", "output")
TIMEOUT_S = 3 * 3600
SKIP_CONTENT = {".pdf"}  # checked for presence only; the PNG carries the same drawing
TEXT_SUFFIXES = {".csv", ".txt", ".md", ".json", ".yaml", ".yml", ".tsv", ".log", ".html"}
VOLATILE = [r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}(:\d{2}(\.\d+)?)?", r"\b\d+(\.\d+)? ?s\b",
            r"\b\d+(\.\d+)? ?(sec|seconds|min)\b"]
ENV = {"SOURCE_DATE_EPOCH": "0", "MPLBACKEND": "Agg", "PYTHONHASHSEED": "0",
       "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"}
# -----------------------------------------------------------------------------

# Run a script as __main__ with some top-level assignments replaced (argv: path, json)
LAUNCHER = """
import ast, json, os, sys
path, sets = sys.argv[1], json.loads(sys.argv[2])
sys.argv = [path] + sys.argv[3:]
tree = ast.parse(open(path, encoding="utf-8").read())
hit = set()
for n in tree.body:
    if isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name) \
            and n.targets[0].id in sets:
        n.value = ast.copy_location(ast.parse(repr(sets[n.targets[0].id]), mode="eval").body, n.value)
        hit.add(n.targets[0].id)
missing = set(sets) - hit
if missing:
    raise SystemExit(f"set: no top-level assignment to {sorted(missing)}")
sys.path.insert(0, os.path.dirname(path))
exec(compile(tree, path, "exec"), {"__name__": "__main__", "__file__": path})
"""


# Size and mtime of every file under the watched folders
def snapshot(root: Path) -> dict[str, tuple[int, int]]:
    out = {}
    for top in WATCHED:
        for dirpath, _, names in os.walk(root / top):
            for n in names:
                p = Path(dirpath) / n
                st = p.stat()
                out[p.relative_to(root).as_posix()] = (st.st_size, st.st_mtime_ns)
    return out


# Put one side's scripts/ in place by renaming, which is instant
def activate(root: Path, side: str) -> None:
    live = root / "scripts"
    if live.exists():
        _rename(live, root / f"_scripts_{(root / '.active').read_text().strip()}")
    _rename(root / f"_scripts_{side}", live)
    (root / ".active").write_text(side)


# Rename, retrying briefly: Windows can hold a folder open for a moment after a run
def _rename(src: Path, dst: Path, tries: int = 20) -> None:
    for i in range(tries):
        try:
            src.rename(dst)
            return
        except PermissionError:
            if i == tries - 1:
                raise
            time.sleep(3)


# One side: run, keep stdout and every file written
def run_side(root: Path, spec: dict, side: str, keep: Path) -> dict:
    activate(root, side)
    before = snapshot(root)
    env = {**os.environ, **ENV, **spec.get("env", {})}
    t0 = time.time()
    script = str(root / "scripts" / spec["script"])
    cmd = ([sys.executable, "-c", LAUNCHER, script, json.dumps(spec["set"])] if spec.get("set")
           else [sys.executable, script])
    p = subprocess.run([*cmd, *spec.get("args", [])],
                       cwd=root / "scripts" / Path(spec["script"]).parent, env=env,
                       capture_output=True, text=True, encoding="utf-8", errors="replace",
                       timeout=TIMEOUT_S)
    after = snapshot(root)
    written = sorted(k for k, v in after.items() if before.get(k) != v)
    keep.mkdir(parents=True, exist_ok=True)
    for rel in written:
        dst = keep / "files" / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / rel, dst)
    (keep / "stdout.txt").write_text(p.stdout, encoding="utf-8")
    (keep / "stderr.txt").write_text(p.stderr, encoding="utf-8")
    return {"code": p.returncode, "seconds": round(time.time() - t0, 1), "written": written,
            "created": [k for k in written if k not in before]}


# Text with the run root and timestamps masked
def normalise(text: str, root: Path) -> str:
    text = text.replace(str(root), "<ROOT>").replace(root.as_posix(), "<ROOT>")
    for pat in VOLATILE:
        text = re.sub(pat, "<T>", text)
    return text


# Same content? Text ignoring timestamps, arrays by value, images by pixel
def same(a: Path, b: Path, root: Path) -> str | None:
    suf = a.suffix.lower()
    if suf in SKIP_CONTENT or a.read_bytes() == b.read_bytes():
        return None
    if suf in TEXT_SUFFIXES:
        ta = normalise(a.read_text(encoding="utf-8", errors="replace"), root)
        tb = normalise(b.read_text(encoding="utf-8", errors="replace"), root)
        return None if ta == tb else "text differs"
    if suf in (".npy", ".npz"):
        la, lb = np.load(a, allow_pickle=True), np.load(b, allow_pickle=True)
        if suf == ".npy":
            return None if np.array_equal(la, lb, equal_nan=True) else "array differs"
        if set(la.files) != set(lb.files):
            return "npz keys differ"
        bad = [k for k in la.files if not np.array_equal(la[k], lb[k], equal_nan=True)]
        return f"arrays differ: {bad[:5]}" if bad else None
    if suf in (".png", ".gif", ".jpg", ".jpeg", ".tif", ".tiff"):
        from PIL import Image, ImageSequence
        ia, ib = Image.open(a), Image.open(b)
        fa = [np.asarray(f.convert("RGBA")) for f in ImageSequence.Iterator(ia)]
        fb = [np.asarray(f.convert("RGBA")) for f in ImageSequence.Iterator(ib)]
        ok = len(fa) == len(fb) and all(x.shape == y.shape and np.array_equal(x, y) for x, y in zip(fa, fb))
        return None if ok else "pixels differ"
    return "bytes differ"


# Every difference between the two sides of one run
def compare(root: Path, label_dir: Path, rb: dict, ra: dict) -> list[str]:
    problems = []
    if rb["code"] != ra["code"]:
        problems.append(f"exit code {rb['code']} -> {ra['code']}")
    if set(rb["written"]) != set(ra["written"]):
        problems.append(f"files written differ: only before {sorted(set(rb['written']) - set(ra['written']))[:5]}, "
                        f"only after {sorted(set(ra['written']) - set(rb['written']))[:5]}")
    for rel in sorted(set(rb["written"]) & set(ra["written"])):
        why = same(label_dir / "before/files" / rel, label_dir / "after/files" / rel, root)
        if why:
            problems.append(f"{rel}: {why}")
    for stream in ("stdout.txt", "stderr.txt"):
        ta = normalise((label_dir / "before" / stream).read_text(encoding="utf-8"), root)
        tb = normalise((label_dir / "after" / stream).read_text(encoding="utf-8"), root)
        if ta != tb:
            problems.append(f"{stream} differs")
    return problems


# Remove the folders the before side's files left empty, so the after side finds none (a run guard counts them)
def prune_empty_dirs(root: Path, created) -> None:
    dirs = {(root / rel).parent for rel in created}
    for d in sorted(dirs, key=lambda p: len(p.parts), reverse=True):
        while d != root and d.is_dir() and not any(d.iterdir()):
            d.rmdir()
            d = d.parent


# Run: install both script versions, then run and compare every spec entry
def main() -> None:
    ap = argparse.ArgumentParser(description="Compare script outputs before and after a restyle.")
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--before", type=Path, required=True)
    ap.add_argument("--after", type=Path, required=True)
    ap.add_argument("--spec", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", nargs="*", help="labels to run")
    a = ap.parse_args()
    a.root, a.before, a.after, a.out = (x.resolve() for x in (a.root, a.before, a.after, a.out))

    # Install both sides as renamable folders inside the root
    for side, src in (("before", a.before), ("after", a.after)):
        dst = a.root / f"_scripts_{side}"
        if (a.root / ".active").exists() and (a.root / ".active").read_text().strip() == side:
            dst = a.root / "scripts"
        shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__"))

    specs = json.loads(a.spec.read_text(encoding="utf-8"))
    report = []
    for spec in specs:
        if a.only and spec["label"] not in a.only:
            continue
        d = a.out / spec["label"]
        shutil.rmtree(d, ignore_errors=True)
        print(f"-- {spec['label']}: {spec['script']} {' '.join(spec.get('args', []))}", flush=True)
        rb = run_side(a.root, spec, "before", d / "before")
        for rel in rb["created"]:               # so the after side writes them fresh, not over them
            (a.root / rel).unlink(missing_ok=True)
        prune_empty_dirs(a.root, rb["created"])
        ra = run_side(a.root, spec, "after", d / "after")
        problems = compare(a.root, d, rb, ra)
        status = "SAME" if not problems else "DIFFERS"
        if rb["code"] != 0:                     # a crash on both sides proves nothing
            status = "CANNOT RUN"
            err = (d / "before" / "stderr.txt").read_text(encoding="utf-8").strip().splitlines()
            problems.insert(0, "before side exits {}: {}".format(rb["code"], err[-1] if err else ""))
        print(f"   {status}  exit {rb['code']}/{ra['code']}  {len(rb['written'])} files  "
              f"{rb['seconds']}s/{ra['seconds']}s", flush=True)
        for p in problems[:15]:
            print(f"     {p}", flush=True)
        report.append({"label": spec["label"], "script": spec["script"], "status": status,
                       "exit": [rb["code"], ra["code"]], "files": len(rb["written"]),
                       "problems": problems})
    (a.out / "report.json").write_text(json.dumps(report, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
