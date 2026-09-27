"""Archive the 2010-2024 runs of the fixed-ends wave grid made on the Hs-1 ends (2026-09-27).

Hannah, 2026-09-27: "re-solve the 2010 ends at Hs 2 and rerun". The 2010-2024
half of wave-climate/2026-09-27-wave-grid-fixed-ends was run with GIS 1 at
+137.581 m/yr, solved at Hs 1, which overshoots at the Hs 2 settings the grid
prefers. Those runs (coarse, refine, cross; both scenarios) move to
output/raw_runs/archive/2026-09-27-fixed-ends-2010-ends-solved-at-hs1/ with
their logs, so the rerun neither skips them nor mixes with them. The
1996-2010 runs stay. Tags in the run metadata are rewritten to the archive
path; the run index is rebuilt by the rerun's scoring.

    python scripts/hatteras_ms/tools/HAT_archive_fixed_ends_2010_20260927.py [--dry-run]
"""
import shutil
import sys
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
RAW = ROOT / "output" / "raw_runs"
STUDY_TAG = "wave-climate/2026-09-27-wave-grid-fixed-ends"
STUDY = RAW / "experiments" / STUDY_TAG
PERIOD = next((a for a in sys.argv[1:] if a.isdigit()), "2010")
WIN = f"{PERIOD}_{int(PERIOD) + 14}"
ARCH_TAG = f"2026-09-27-fixed-ends-{PERIOD}-ends-solved-at-hs1"
ARCH = RAW / "archive" / ARCH_TAG
DRY = "--dry-run" in sys.argv


def nruns(d):
    return len(list(d.rglob("*_run_metadata.json")))


moves = []
for kind in ("runs", "logs"):
    for src in sorted((STUDY / kind).glob(f"*/{WIN}")):
        moves.append((src, ARCH / kind / src.parent.name / WIN))
before = sum(nruns(s) for s, _ in moves if s.parts[-3] == "runs")
print(f"{len(moves)} folders to move, {before} runs")
for s, d in moves:
    print("  ", s.relative_to(RAW), "->", d.relative_to(RAW))
if DRY:
    sys.exit(0)
assert not ARCH.exists(), ARCH
for s, d in moves:
    d.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(s), str(d))
n = 0
for f in list((ARCH / "runs").rglob("*_run_metadata.json")) + list((ARCH / "runs").rglob("*_run_metadata.txt")):
    t0 = t = f.read_text(encoding="utf-8")
    t = t.replace(f"{STUDY_TAG}/runs/", f"{ARCH_TAG}/runs/")
    t = t.replace('"run_kind": "experiment"', '"run_kind": "archive"').replace('"kind": "experiment"', '"kind": "archive"')
    if t != t0:
        f.write_text(t, encoding="utf-8")
        n += 1
after = nruns(ARCH / "runs")
assert after == before, (before, after)
OLD_ENDS = {"2010": "GIS 1 +137.581 (accepted, not converged), GIS 90 +14.731",
            "1996": "GIS 1 +7.0245, GIS 90 +9.8283"}[PERIOD]
MISS = {"2010": "GIS 1 overshoots (model about +10 against +6.9 m/yr observed; several hundred "
                "metres of position change at GIS 1)",
        "1996": "the ends miss by about +2.4 (GIS 1) and -1.9 m/yr (GIS 90) in the managed run"}[PERIOD]
(ARCH / "README.md").write_text(f"""# {ARCH_TAG}

The {WIN} runs of `experiments/{STUDY_TAG}/` (coarse, refine, cross; natural
and full management) made with the {WIN} ends solved at Hs 1 / Tp 8 /
asym 0.8 / high-angle 0.45: {OLD_ENDS} m/yr. At the adopted waves (Hs 2 / Tp 7.5
/ asym 0.6 / high-angle 0.5) {MISS}. Superseded 2026-09-27 (Hannah: re-solve the
ends at the adopted waves, for both windows); the targeted rerun is in the study
folder. {after} runs; their logs under `logs/`. Kept intact; the study README
records what they showed.
""", encoding="utf-8")
print(f"moved {after} runs; {n} metadata files retagged")
