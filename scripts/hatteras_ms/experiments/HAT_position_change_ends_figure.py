"""Figure 3 redrawn on the position-change ends (2026-09-28).

Hannah: "redraw the figure with the position-change ends when done". The solve
(HAT_resolve_ends_on_position_change.py) runs full management only; this runs
the natural scenario at the same ends (option A waves), then draws
HAT_wave_recommendation_figures.fig3 on the two runs per window. The LRR-ends
figure is left as it is, for comparison.

    python scripts/hatteras_ms/experiments/HAT_position_change_ends_figure.py

Natural runs land under <study>/runs/natural_final/; the figure goes to the
study's figures/ folder.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_resolve_ends_on_position_change as S  # noqa: E402
import HAT_wave_recommendation_figures as W  # noqa: E402

E, grid = S.E, S.R.grid
STUDY_DIR = grid.RAW_RUNS / "experiments" / S.TAG
E.TAG, E.STUDY_DIR = S.TAG, STUDY_DIR
E.TABLES_DIR, E.LOGS_DIR = STUDY_DIR / "tables", STUDY_DIR / "logs"
WAVES = pd.Series(W.REC)


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    j = json.loads((E.TABLES_DIR / "ends.json").read_text(encoding="utf-8"))
    ends = {int(p): {1: v["1"], 90: v["90"]} for p, v in j["ends_m_yr"].items()}
    rep = pd.read_csv(E.TABLES_DIR / f"ends_1996_2010_{E.label(WAVES)}.csv")
    runs = {}
    for _, r in rep.iterrows():
        p = int(r.period[:4])
        runs[(p, "full_management")] = STUDY_DIR / r.run_dir
    grid.check_barrier3d()
    jobs = [("natural", "final", p, WAVES, ends[p]) for p in ends]
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(E.launch, jobs))
    for p in ends:
        d = E.find_run("natural", "final", p, WAVES)
        if d is None:
            raise SystemExit(f"natural {p}: no run ({E.step2.stop_reason(E.log_path('natural', 'final', p, WAVES))})")
        runs[(p, "natural")] = d
    out = STUDY_DIR / "figures"
    out.mkdir(parents=True, exist_ok=True)
    W.apply_style()
    with plt.rc_context(W.p2.SCREEN_RC):
        png = W.fig3(runs=runs, ends={p: (e[1], e[90]) for p, e in ends.items()},
                     png=out / "recommended_vs_coastsat_position_change_ends.png",
                     solved_on="the CoastSat end-minus-start position change (GIS 90 LOWESS 7)")
    print(png)
    return 0


if __name__ == "__main__":
    sys.exit(main())
