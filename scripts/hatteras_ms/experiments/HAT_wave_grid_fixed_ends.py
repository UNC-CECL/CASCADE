"""The four-parameter wave grid with the end domains FIXED (2026-09-27).

Hannah, 2026-09-26/27: sweep the wave climate with the end source/sink
terms held fixed, the ends first re-solved under the metres offset. The
fixed pair per window comes from
end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/tables/ends.json
(solved at the step-2 baseline waves, full management, against CoastSat) and
is imposed in every run -- natural and managed alike -- through the edgeBE
preset and HAT_BE_OVERRIDE.

Everything else is HAT_wave_grid_smoothed_score's, unchanged (the same coarse
grid, refine, cross-runs and smoothed score), so this study and the zeroBE
grid (wave-climate/2026-09-25-wave-grid-smoothed-score) compare one to one.
This module only points that driver at its own folder, preset and ends.

WHERE: output/raw_runs/experiments/wave-climate/2026-09-27-wave-grid-fixed-ends/

    python scripts/hatteras_ms/experiments/HAT_wave_grid_fixed_ends.py run all --jobs 8
    python scripts/hatteras_ms/experiments/HAT_wave_grid_fixed_ends.py score

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_wave_grid_smoothed_score as grid  # noqa: E402

ENDS_FILE = (grid.RAW_RUNS / "experiments" / "end-domain-boundaries"
             / "2026-09-27-ends-resolved-metres-offset" / "tables" / "ends.json")
ENDS = {int(p): v for p, v in json.loads(ENDS_FILE.read_text(encoding="utf-8"))["ends_m_yr"].items()}

grid.TAG = "wave-climate/2026-09-27-wave-grid-fixed-ends"
grid.STUDY_DIR = grid.RAW_RUNS / "experiments" / grid.TAG
grid.TABLES_DIR = grid.STUDY_DIR / "tables"
grid.LOGS_DIR = grid.STUDY_DIR / "logs"
grid.PRESET = "edgeBE"
_base_env = grid.run_env


def run_env(phase, scenario, period, s):
    e = _base_env(phase, scenario, period, s)
    e["HAT_BE_OVERRIDE"] = f"1={ENDS[period]['1']},90={ENDS[period]['90']}"
    return e


grid.run_env = run_env


def build_targeted(period, archive_tag=None, n=15):
    """The targeted list for one window (2026-09-27): the top n on the raw score
    in that window from this study's runs on the superseded ends (read from
    archive_tag if already archived, else from tables/all_runs.csv) and from the
    zeroBE grid, the other window's top 5 here, and the adopted candidates."""
    import json
    import pandas as pd
    K = list(grid.KEYS)
    tg = grid.common.coastsat_target(period)
    if archive_tag:
        rows = []
        for md in (grid.RAW_RUNS / "archive" / archive_tag / "runs").glob(f"*/{grid.window(period)}/edgeBE/*/*_run_metadata.json"):
            w = json.loads(md.read_text(encoding="utf-8"))["wave climate"]
            rows.append(dict(scenario=md.parts[-5].split("_", 1)[1], hs=float(w["wave_height_m"]),
                             wave_period_s=float(w["wave_period_s"]), wave_asymmetry=float(w["wave_asymmetry"]),
                             wave_angle_high_fraction=float(w["wave_angle_high_frac"]),
                             raw_variance_explained=grid.score_run(md.parent, tg)["raw_variance_explained"]))
        old = pd.DataFrame(rows)
    else:
        old = pd.read_csv(grid.TABLES_DIR / "all_runs.csv")
        old = old[(old.status == "scored") & (old.period_start == period)]
    z = pd.read_csv(grid.RAW_RUNS / "experiments" / "wave-climate" / "2026-09-25-wave-grid-smoothed-score" / "tables" / "all_runs.csv")
    z = z[(z.status == "scored") & (z.period_start == period)]
    here = pd.read_csv(grid.TABLES_DIR / "all_runs.csv")
    other = here[(here.status == "scored") & (here.period_start != period)]
    cells = set()
    for sc in grid.SCENARIOS:
        for src, m in ((old, n), (z, n), (other, 5)):
            x = src[src.scenario == sc].drop_duplicates(K).nlargest(m, "raw_variance_explained")
            cells |= {(sc, *[float(r[k]) for k in K]) for _, r in x.iterrows()}
        for c in ([(2.0, 7.5, 0.6, 0.5)] if sc == "full_management" else [(2.0, 10.0, 0.6, 0.5), (2.0, 10.0, 0.5, 0.5)]):
            cells.add((sc, *c))
    t = pd.DataFrame(sorted(cells), columns=["scenario", *K])
    t.to_csv(grid.TABLES_DIR / f"targeted_{period}_cells.csv", index=False)
    print(f"targeted {period}: {len(t)} cells ({t.groupby('scenario').size().to_dict()})", flush=True)
    return t


def run_targeted(jobs=8, period=2010):
    """The 2010-2024 rerun, TARGETED (2026-09-27, Hannah: "is there not a faster
    version"): after the 2010 ends were re-solved at Hs 2 (GIS 1 +137.6 -> +18.8),
    only the settings that can change a decision are re-run, as phase "target":
    the top 15 on the raw score in 2010-2024 from the superseded fixed-ends runs
    (archive/2026-09-27-fixed-ends-2010-ends-solved-at-hs1) and from the zeroBE
    grid, the 1996-2010 top 5 (for the one-setting and one-parameter picks),
    and the adopted candidates. The list is tables/targeted_2010_cells.csv.
    Coarse cells already run on the new ends are kept."""
    import pandas as pd
    f = grid.TABLES_DIR / f"targeted_{period}_cells.csv"
    if not f.is_file() and period == 2010:
        f = grid.TABLES_DIR / "targeted_2010_cells.csv"
    t = pd.read_csv(f)
    cells = [("target", r.scenario, period, {k: float(r[k]) for k in grid.KEYS}) for _, r in t.iterrows()]
    grid.check_barrier3d()
    grid.common.keep_awake()
    grid.run_cells(cells, jobs)
    grid.cmd_score()


if __name__ == "__main__":
    if sys.argv[1:2] == ["build-targeted"]:
        build_targeted(int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else None)
        sys.exit(0)
    if sys.argv[1:2] == ["targeted"]:
        per = int(sys.argv[2]) if len(sys.argv) > 2 else 2010
        print(f"fixed ends (m/yr): {per}: GIS1 {ENDS[per]['1']:+.3f}, GIS90 {ENDS[per]['90']:+.3f}", flush=True)
        sys.exit(run_targeted(period=per))
    print(f"fixed ends (m/yr): " + ", ".join(f"{p}: GIS1 {v['1']:+.3f}, GIS90 {v['90']:+.3f}"
                                            for p, v in ENDS.items()), flush=True)
    sys.exit(grid.main())
