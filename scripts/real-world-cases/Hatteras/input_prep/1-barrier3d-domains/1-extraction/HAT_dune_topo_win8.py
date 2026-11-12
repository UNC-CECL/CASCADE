"""
The 1996-start dune and interior arrays with the earlier fixed dune search (2026-11-11, test): each profile's dune is
the highest cell within 8 cells (80 m) landward of its own beach start (the first cell more than 0.5 m above MHW), as
the Benton-parity extractor did (CASCADE PeaIsland branch, scripts/dune_topo_extractor.py, DUNE_WINDOW_PX = 8),
instead of v3's hand-picked window per domain. Everything else is HAT_dune_topo_extractor.py's: the same DEM
profiles, straightening, berm, minimum dune height and interior built behind each profile's own dune.

    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/1-extraction/HAT_dune_topo_win8.py

Writes dune-topo/v3_win8/ beside v3 in the 1996-start product (D:). CURRENT stays v3; a run picks this set with
HAT_TOPO_VERSION_1996_START=v3_win8 (HAT_dune_collision.py WIN8_START). Also writes v3_win8/start_crest_compare.csv:
each domain's median starting crest, hand window against 8 cells, and road_setback_add_m, how far seaward of v3's
the 8-cell dune sits (median over profiles), which the runs add to the road setback so NC-12 stays in place.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-11
"""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_spec = importlib.util.spec_from_file_location("hat_extractor", _HERE.parent / "HAT_dune_topo_extractor.py")
X = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(X)

# --- CONFIG ------------------------------------------------------------------
WINDOW_PX = 8                                        # Benton-parity DUNE_WINDOW_PX
VERSION = "v3_win8"
OUT_DIR = X.DUNE_TOPO_ROOT / VERSION
# -----------------------------------------------------------------------------


def find_dunes_win8(prof_arr, start_beach, i0=None, i1=None):
    """Highest cell in [beach start, beach start + WINDOW_PX) per profile; i0/i1 ignored."""
    n_along = prof_arr.shape[0]
    elev, loc = np.full(n_along, np.nan), np.full(n_along, -1, dtype=int)
    for i, s in enumerate(start_beach):
        if s < 0:
            continue
        w = prof_arr[i, s:s + WINDOW_PX]
        ok = w > X.SENTINEL_WATER_M + 1e-9
        if ok.any():
            k = int(np.argmax(np.where(ok, w, -np.inf)))
            elev[i], loc[i] = float(w[k]), int(s + k)
    return elev, loc


def main():
    topo_dir, dune_dir = OUT_DIR / "topography", OUT_DIR / "dunes"
    topo_dir.mkdir(parents=True, exist_ok=True)
    dune_dir.mkdir(parents=True, exist_ok=True)
    windows = json.loads(Path(X.WINDOW_JSON).read_text(encoding="utf-8"))
    hand_find = X.find_dunes
    rows = []
    for g in range(1, 91):
        stem = f"domain_{g}"
        with contextlib.redirect_stdout(io.StringIO()):
            dom = X.load_profiles(Path(X.LOAD_PATH) / f"{stem}.npy")
        z, sb = dom["z"], dom["start_beach"]
        w = windows[stem]
        hand, hand_loc = hand_find(z, sb, int(w["i0"]), int(w["i1"]))
        X.find_dunes = find_dunes_win8                  # extract_domain calls X.find_dunes
        try:
            res = X.extract_domain(stem, z, sb, 0, 0, topo_dir, dune_dir, shear=dom["shear"],
                                   obliquity_deg=dom["obliquity_deg"])
        finally:
            X.find_dunes = hand_find
        rows.append(dict(gis=g, hand_crest_m_navd88=float(np.nanmedian(hand)) + X.MHW_M,
                         win8_crest_m_navd88=float(np.nanmedian(res["dune_elev"])) + X.MHW_M,
                         interior_rows_win8=res["interior_rows"],
                         # how far seaward the 8-cell dune sits of v3's, median over profiles, whole cells: the road
                         # setback (measured from v3's dune) grows by this so NC-12 stays where it is
                         road_setback_add_m=X.CELL_SIZE_M * float(np.median(
                             (hand_loc - res["dune_loc"])[(hand_loc >= 0) & (res["dune_loc"] >= 0)]))))
    t = pd.DataFrame(rows)
    t["win8_minus_hand_m"] = t.win8_crest_m_navd88 - t.hand_crest_m_navd88
    t.round(3).to_csv(OUT_DIR / "start_crest_compare.csv", index=False)
    (OUT_DIR / "README.md").write_text(
        "# v3_win8 (test, 2026-11-11)\n\nv3's profiles with the earlier fixed dune search: the highest cell within 8 "
        "cells (80 m) of each profile's beach start, not v3's hand-picked windows. Made by "
        "`input_prep/1-barrier3d-domains/1-extraction/HAT_dune_topo_win8.py`. Not CURRENT; runs select it with "
        "HAT_TOPO_VERSION_1996_START=v3_win8.\n", encoding="utf-8")
    print(t.round(2).to_string(index=False))
    print(f"wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
