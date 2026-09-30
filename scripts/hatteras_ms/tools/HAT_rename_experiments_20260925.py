"""Rename the older experiments and group every experiment by theme (2026-09-25).

Hannah, 2026-09-25: "rename the old experiments so their names are clearer
about what they tested", then "organize the experiment folders ... by what
the theme or investigation is". Names and themes approved the same day:

    experiments/<theme>/<date>-<what it tested>/

What this does, while nothing is running:
  1. removes the parameters-only folders a paused sweep leaves behind (a run
     killed at start-up), which would block its resume
  2. moves each study to <theme>/<new name>
  3. rewrites the old path/name inside every text file of every study
     (run-metadata tags, NOTE/README/FINDINGS, logs and CSVs naming runs),
     fixing the relative links between studies (one level deeper now) and
     notes the old name at the top of the renamed studies' NOTE.md
  4. rewrites the old names in every tracked text file outside the archive
     (scripts -- a study's folder is also its run tag -- comparison READMEs
     and captions) and in the memory notes
  5. writes experiments/README.md (the map) and a README per theme
  6. rebuilds the run index and checks every run is under its new name
The archive and retired_runs.csv are history and are left as they were.
cascade_pipeline.run_registry.check_tag allows 4 tag levels since this change.

    python scripts/hatteras_ms/tools/HAT_rename_experiments_20260925.py [--dry-run]

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
EXP = ROOT / "output" / "raw_runs" / "experiments"
MEMORY = Path.home() / ".claude" / "projects" / "C--Users-hanna-PycharmProjects-CASCADE" / "memory"

THEMES = {
    "island-offset": "How the island's planform (the BRIE island offset) is set: its units, "
                     "and whether it comes from the dune line or the CoastSat shoreline.",
    "wave-climate": "Tuning the four wave parameters (Hs, Tp, asymmetry, high-angle fraction) "
                    "against the CoastSat alongshore rates.",
    "end-domain-boundaries": "The source/sink rates locked at the two end domains (GIS 1 and "
                             "90): what they must carry against each target.",
    "topography-and-domains": "What the model domains contain: dune footprints, rows added "
                              "or removed, and extending the reach onto Pea Island.",
    "code-checks": "Whether a code or model change moves the results: re-runs against stored "
                   "runs, and the Barrier3D route_overwash fix.",
}
# old folder name -> (theme, new name)
MOVES = {
    "2026-09-22-shoreline-offset": ("island-offset", "2026-09-22-shoreline-offset-at-div10-scale"),
    "2026-09-24-metres-1-offset-units": ("island-offset", "2026-09-24-metres-1-offset-units"),
    "2026-09-25-offset-source-duneline-vs-shoreline":
        ("island-offset", "2026-09-25-offset-source-duneline-vs-shoreline"),
    "2026-09-24-metres-2-wave-sensitivity": ("wave-climate", "2026-09-24-metres-2-wave-sensitivity"),
    "2026-09-25-wave-grid-smoothed-score": ("wave-climate", "2026-09-25-wave-grid-smoothed-score"),
    "2026-09-16-dune-edgesolve": ("end-domain-boundaries", "2026-09-16-end-domains-solved-on-duneline"),
    "2026-09-18-dune-edgesolve":
        ("end-domain-boundaries", "2026-09-18-end-domains-solved-on-redigitized-duneline"),
    "2026-09-19-edgesolve-2010": ("end-domain-boundaries", "2026-09-19-end-domains-2010-recheck"),
    "2026-09-19-edgesolve-lrr1996_2024":
        ("end-domain-boundaries", "2026-09-19-end-domains-solved-on-lrr-1996-2024"),
    "2026-09-02-pea1989": ("topography-and-domains", "2026-09-02-pea-island-row-insert-control"),
    "2026-09-08-behindroad-copy": ("topography-and-domains", "2026-09-08-dune-footprint-behind-road"),
    "2026-09-16-peaisland-ext": ("topography-and-domains", "2026-09-16-pea-island-domain-extension"),
    "2026-09-14-currency": ("code-checks", "2026-09-14-calibrated-pair-rerun-current-code"),
    "2026-09-14-paramsplit": ("code-checks", "2026-09-14-site-config-split-check"),
    "2026-09-14-probe": ("code-checks", "2026-09-14-relocation-rounding-probes"),
    "2026-09-14-recode": ("code-checks", "2026-09-14-relocation-arm-rerun-new-code"),
    "2026-09-24-metres-3-barrier3d-overwash-fix":
        ("code-checks", "2026-09-24-metres-3-barrier3d-overwash-fix"),
}
CHAIN_INDEX = "2026-09-24-metres-INDEX.md"      # stays at experiments/, links updated
PAUSED_SWEEP = "2026-09-25-wave-grid-smoothed-score"
TEXT = {".md", ".py", ".txt", ".json", ".yaml", ".yml", ".csv", ".log", ".ipynb", ".jsonl"}
DRY = "--dry-run" in sys.argv


def new_path(old):
    theme, name = MOVES[old]
    return f"{theme}/{name}"


def replacements(inside_study):
    """Ordered (old, new) string pairs. Relative links first: from inside a
    study, ../<other> becomes ../../<theme>/<other> and the chain INDEX moves
    one level up; then every remaining bare name gets its theme."""
    pairs = []
    if inside_study:
        pairs.append((f"../{CHAIN_INDEX}", f"../../{CHAIN_INDEX}"))
        for old in MOVES:
            pairs.append((f"../{old}", f"../../{new_path(old)}"))
    for old in sorted(MOVES, key=len, reverse=True):
        pairs.append((old, new_path(old)))
    return pairs


def rewrite(path, pairs):
    try:
        s0 = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return 0
    s = s0
    for a, b in pairs:
        if a in s:
            # never prefix a name that already carries its theme
            s = s.replace(b, "\0").replace(a, b).replace("\0", b)
    if s != s0:
        if not DRY:
            path.write_text(s, encoding="utf-8")
        return 1
    return 0


def runs(d):
    return len(list(d.rglob("*_run_metadata.json")))


def partial_run_dirs():
    out = []
    # 2010_2024 only: the pause stops the sweep as its 2010-2024 cells start.
    # A DROWNED run also leaves only its parameters file, and those (1996-2010,
    # Hs 0.75 / Tp 10) are results that must stay.
    for d in (EXP / PAUSED_SWEEP).glob("runs/*/2010_2024/*/*"):
        files = [p for p in d.iterdir()] if d.is_dir() else []
        if len(files) == 1 and files[0].name.endswith("-parameters.yaml"):
            out.append(d)
    return out


def write_readmes():
    lines = ["# Experiments, by theme", "",
             "One question each, grouped by the investigation it belongs to "
             "(reorganised 2026-09-25, Hannah). Each study folder is "
             "`<date>-<what it tested>`; its README or NOTE is the record. A study's "
             "folder path under `experiments/` is also its runs' tag in `run_index.csv`.", ""]
    for theme, what in THEMES.items():
        lines += [f"## [`{theme}/`]({theme}/README.md)", "", what, ""]
        studies = sorted(n for o, (t, n) in MOVES.items() if t == theme)
        lines += [f"- `{n}`" for n in studies] + [""]
        tl = [f"# {theme}", "", what, "", "| study | record |", "|---|---|"]
        for n in studies:
            rec = next((f for f in ("README.md", "NOTE.md", "FINDINGS.md")
                        if (EXP / theme / n / f).is_file()), None)
            tl.append(f"| `{n}` | " + (f"[{rec}]({n}/{rec})" if rec else "") + " |")
        tl += ["", "Back to [the map](../README.md)."]
        if not DRY:
            (EXP / theme / "README.md").write_text("\n".join(tl) + "\n", encoding="utf-8")
    lines += ["## Chains across themes", "",
              f"- [`{CHAIN_INDEX}`]({CHAIN_INDEX}): the 2026-09-24 metres work, "
              "steps 1-3 (offset units, wave sensitivity, the Barrier3D fix), and its "
              "2026-09-25 follow-ups.", "",
              "Renamed on 2026-09-25 (old name → new): "] + [
        f"- `{o}` → `{new_path(o)}`" for o, (t, n) in MOVES.items() if o != n]
    if not DRY:
        (EXP / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    # Resumable (2026-09-25: the first pass stopped at a folder held open by
    # a shell's working directory): a study already at its new place is
    # skipped; the text rewrite cannot double a prefix, and the rename note
    # is written once.
    done = {o for o in MOVES if not (EXP / o).exists() and (EXP / new_path(o)).is_dir()}
    for old in MOVES:
        if old not in done:
            assert (EXP / old).is_dir(), f"missing {old}"
            assert not (EXP / new_path(old)).exists(), f"{new_path(old)} already exists"
    if done:
        print(f"already moved: {sorted(done)}")
    stray = sorted(p.name for p in EXP.iterdir()
                   if p.is_dir() and p.name not in MOVES and p.name not in THEMES)
    assert not stray, f"experiments not assigned a theme: {stray}"
    partial = partial_run_dirs()
    before = {old: runs(EXP / (new_path(old) if old in done else old)) for old in MOVES}
    print(f"{len(MOVES)} studies, {sum(before.values())} runs; "
          f"{len(partial)} parameters-only folders from the paused sweep to remove")
    tracked = subprocess.run(
        ["git", "-C", str(ROOT), "grep", "-l", "-e", "\\|".join(MOVES), "--", ".",
         ":!output/raw_runs/experiments", ":!output/raw_runs/archive", ":!*retired_runs.csv",
         ":!*run_index.csv", ":!scripts/hatteras_ms/tools/HAT_rename_experiments_20260925.py"],
        capture_output=True, text=True).stdout.split()
    # the new drivers are not tracked yet
    extra = [str(p.relative_to(ROOT)).replace("\\", "/")
             for p in (ROOT / "scripts" / "hatteras_ms" / "experiments").glob("*.py")
             if any(o in p.read_text(encoding="utf-8") for o in MOVES)]
    tracked = sorted(set(tracked) | set(extra))
    print(f"files outside the experiments to update: {len(tracked)}")
    for f in tracked:
        print("   ", f)
    mem = [p for p in MEMORY.glob("*.md") if any(o in p.read_text(encoding="utf-8") for o in MOVES)]
    print(f"memory notes to update: {len(mem)}")
    if DRY:
        return 0
    import shutil
    for d in partial:
        shutil.rmtree(d)
    for t in THEMES:
        (EXP / t).mkdir(exist_ok=True)
    n_in = 0
    for old in MOVES:
        dst = EXP / new_path(old)
        if old in done:
            continue
        (EXP / old).rename(dst)
        for p in dst.rglob("*"):
            if p.is_file() and p.suffix.lower() in TEXT:
                n_in += rewrite(p, replacements(inside_study=True))
        theme, name = MOVES[old]
        note = dst / "NOTE.md"
        if old != name and note.is_file() and "Renamed 2026-09-25" not in note.read_text(encoding="utf-8"):
            s = note.read_text(encoding="utf-8")
            first, _, rest = s.partition("\n")
            note.write_text(f"{first}\n\n*Renamed 2026-09-25 from `{old}` and filed under "
                            f"`{theme}/` (Hannah: names say what was tested, grouped by "
                            f"theme).*\n{rest}", encoding="utf-8")
    n_in += rewrite(EXP / CHAIN_INDEX, replacements(inside_study=False))
    n_out = sum(rewrite(ROOT / f, replacements(inside_study=False)) for f in tracked)
    n_mem = sum(rewrite(p, replacements(inside_study=False)) for p in mem)
    write_readmes()
    after = {old: runs(EXP / new_path(old)) for old in MOVES}
    assert before == after, (before, after)
    print(f"rewritten: {n_in} files in the studies, {n_out} outside, {n_mem} memory notes")
    sys.path.insert(0, str(ROOT / "scripts"))
    from cascade_pipeline.run_registry import rebuild_run_index
    rebuild_run_index(ROOT / "output" / "raw_runs")
    import pandas as pd
    idx = pd.read_csv(ROOT / "output" / "raw_runs" / "run_index.csv")
    exp = idx[idx["kind"] == "experiment"]
    first = exp["tag"].astype(str).str.split("/").str[0]
    assert first.isin(list(THEMES)).all(), sorted(set(first) - set(THEMES))
    print(f"index rebuilt: {len(exp)} experiment runs, every one under a theme")
    return 0


if __name__ == "__main__":
    sys.exit(main())
