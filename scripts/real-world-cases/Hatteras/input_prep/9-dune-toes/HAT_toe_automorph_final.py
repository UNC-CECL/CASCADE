"""
The Automorph dune toes (Itzkin et al. 2020, adapted for Hatteras: HAT_toe_method_automorph.py) as checked by hand in
HAT_toe_automorph_review.py: one toe per transect and one per domain.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_automorph_final.py

Steps:
1. The transect toe: the toe you accepted or moved in the review (automorph_review_fresh_1996.json).
2. The screen of the adapted method, which the review window does not apply: at least TOE_MIN_M_NAVD88, and at least
   TOE_BELOW_CREST_M below the crest. A toe you moved by hand is kept regardless.
3. The domain toe: the median of the kept transects, when at least MIN_KEPT; otherwise the mean of the nearest domains
   each side that have one (source "filled"), as HAT_toe_method_line.py.

Writes, in data/hatteras_init/9-dune-toes/1-transects-1996/1-toes/ (D:): automorph_transect_toes_1996.csv and
automorph_domain_toes_1996.csv (same gis and toe_m_navd88 columns as domain_toes_1996.csv, with that table's toe and
the current toe beside it). Not read by the model yet.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-09
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from HAT_toe_common import BERM_M, MHW, PICKS, TOES, TOE_BELOW_CREST_M, TOE_MIN_M_NAVD88, YEAR  # noqa: E402
from HAT_toe_method_line import _medians  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
REVIEW = PICKS / f"automorph_review_fresh_{YEAR}.json"
# -----------------------------------------------------------------------------


def transect_toes():
    t = pd.DataFrame(json.loads(REVIEW.read_text(encoding="utf-8"))).T.rename_axis("transect").reset_index()
    t["transect"], t["gis"] = t.transect.astype(int), t.gis.astype(int)
    for c in ("crest_x", "crest_z", "toe_x", "toe_z"):
        t[c] = pd.to_numeric(t[c], errors="coerce")
    checked, moved = t.toe_status.isin(["accepted", "moved"]), t.toe_status == "moved"
    t["passes_screen"] = (t.toe_z >= TOE_MIN_M_NAVD88) & (t.crest_z - t.toe_z >= TOE_BELOW_CREST_M)
    t["kept"] = checked & (t.passes_screen | moved)
    t["toe_m_navd88"] = t.toe_z.where(t.kept)
    return t[["transect", "gis", "crest_status", "toe_status", "crest_x", "crest_z", "toe_x", "toe_z", "passes_screen",
              "kept", "toe_m_navd88"]].sort_values(["gis", "transect"])


def main():
    t = transect_toes()
    d = pd.DataFrame({"gis": range(1, 91)})
    d["toe_m_navd88"], d["n_kept"], d["source"] = _medians(t)
    d["toe_m_mhw"] = d.toe_m_navd88 - MHW
    old = pd.read_csv(TOES / f"domain_toes_{YEAR}.csv").set_index("gis")
    d["line_toe_m_navd88"] = d.gis.map(old.toe_m_navd88)
    d["current_toe_m_navd88"] = d.gis.map(old.current_toe_m_navd88)
    TOES.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(TOES / f"automorph_transect_toes_{YEAR}.csv", index=False)
    d.round(3).to_csv(TOES / f"automorph_domain_toes_{YEAR}.csv", index=False)

    checked = t.toe_status.isin(["accepted", "moved"])
    print(f"{int(t.kept.sum())} of {int(checked.sum())} checked toes kept; dropped by the screen: "
          f"{t[checked & ~t.kept].transect.tolist()}")
    print(f"domains: {int((d.source == 'median').sum())} medians, filled {d[d.source == 'filled'].gis.tolist()}")
    print(f"domain toe: mean {d.toe_m_navd88.mean():.2f}, range {d.toe_m_navd88.min():.2f}-{d.toe_m_navd88.max():.2f} m "
          f"NAVD88; below the {BERM_M:g} m berm: {d[d.toe_m_navd88 < BERM_M].gis.tolist()}")
    for name, col in (("current", "current_toe_m_navd88"), ("line method", "line_toe_m_navd88")):
        dz = (d.toe_m_navd88 - d[col]).dropna()
        print(f"vs {name}: mean {dz.mean():+.2f} m, RMSE {np.sqrt((dz ** 2).mean()):.2f} m, more than 0.5 m apart "
              f"{int((dz.abs() > 0.5).sum())}, largest {dz.abs().max():.2f} m (GIS {d.gis[dz.abs().idxmax()]})")
    print(f"wrote {TOES}")


if __name__ == "__main__":
    main()
