"""
A second dune row for the 1996 start (2026-11-11, test): the seaward (front) row stays v3's, the highest cell in each
profile's hand-picked dune window; the landward (back) row becomes the median elevation of the cells in that same
window (clipped to the profile's beach start, as find_dunes clips it), so the dune is no longer two equal rows at the
crest height. Windows, interior and road are v3's.

    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/1-extraction/HAT_dune_topo_maxmed.py

Barrier3D copies a 50-value dune file into both rows (DuneParamMultipleRows unset), so the back row is not written as
a dune file: it goes to dune-topo/v3_maxmed/dunes_back/domain_<N>_dune_back.npy (50 values, dam above the berm, in
the order of v3's dune files, floored at v3's minimum dune height) and HAT_dune_collision.py (MAXMED_START) puts it in
DuneDomain[0, :, 1] after the model is built, checking that row 0 matches v3's file first. Also writes
v3_maxmed/back_row_compare.csv: each domain's median front and back row and window width.

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
VERSION = "v3_maxmed"
OUT_DIR = X.DUNE_TOPO_ROOT / VERSION
V3_DUNES = X.DUNE_TOPO_ROOT / "v3" / "dunes"
# -----------------------------------------------------------------------------


def window_median(prof_arr, start_beach, i0, i1):
    """Median elevation of each profile's window cells, the window clipped to its beach start as find_dunes does."""
    n_along, n_cross = prof_arr.shape
    med = np.full(n_along, np.nan)
    for i in range(n_along):
        a = max(i0, int(start_beach[i])) if (X.CLIP_WINDOW_TO_BEACH and start_beach[i] >= 0) else i0
        w = prof_arr[i, a:min(i1, n_cross)]
        w = w[w > X.SENTINEL_WATER_M + 1e-9]
        if w.size:
            med[i] = float(np.median(w))
    return med


def main():
    back_dir = OUT_DIR / "dunes_back"
    back_dir.mkdir(parents=True, exist_ok=True)
    windows = json.loads(Path(X.WINDOW_JSON).read_text(encoding="utf-8"))
    berm_mhw = X.BERM_ELEV_NAVD_M - X.MHW_M
    rows = []
    for g in range(1, 91):
        stem = f"domain_{g}"
        with contextlib.redirect_stdout(io.StringIO()):
            dom = X.load_profiles(Path(X.LOAD_PATH) / f"{stem}.npy")
        z, sb = dom["z"], dom["start_beach"]
        i0, i1 = int(windows[stem]["i0"]), int(windows[stem]["i1"])
        crest, _ = X.find_dunes(z, sb, i0, i1)
        front = np.full(X.ALONG_COLS, np.nan)
        front[:len(crest)] = np.where(np.isfinite(crest), np.maximum(crest - berm_mhw, X.MIN_DUNE_H_M),
                                      X.MIN_DUNE_H_M)[:X.ALONG_COLS]
        # The extractor's order must be the dune file's: rebuild v3's front row and compare before writing
        v3 = np.load(V3_DUNES / f"{stem}_dune.npy") * 10
        if not np.allclose(front, v3, atol=1e-6):
            raise SystemExit(f"{stem}: front row does not reproduce v3's dune file (max diff "
                             f"{np.nanmax(np.abs(front - v3)):.3f} m); the back row would be misaligned")
        med = window_median(z, sb, i0, i1)[:X.ALONG_COLS]
        back = np.where(np.isfinite(med), np.maximum(med - berm_mhw, X.MIN_DUNE_H_M), X.MIN_DUNE_H_M)
        back = np.minimum(back, front)                    # never above the front row's crest
        np.save(back_dir / f"{stem}_dune_back.npy", back * 0.1)
        rows.append(dict(gis=g, window_cells=i1 - i0, front_m_above_berm=float(np.median(front)),
                         back_m_above_berm=float(np.median(back)),
                         back_over_front=float(np.median(back / front))))
    t = pd.DataFrame(rows)
    t.round(3).to_csv(OUT_DIR / "back_row_compare.csv", index=False)
    (OUT_DIR / "README.md").write_text(
        "# v3_maxmed (test, 2026-11-11)\n\nv3's dunes with a distinct back row: the front (seaward) row v3's window "
        "maximum, the back row the median of the same window's cells, in `dunes_back/` (50 values, dam above the "
        "berm). Made by `input_prep/1-barrier3d-domains/1-extraction/HAT_dune_topo_maxmed.py`; HAT_dune_collision.py "
        "(MAXMED_START) sets DuneDomain[0, :, 1] from it after the build. Not a dune-topo version: runs keep v3's "
        "files.\n", encoding="utf-8")
    pd.set_option("display.width", 160)
    print(t.round(2).describe().to_string())
    print(f"wrote {OUT_DIR}")


if __name__ == "__main__":
    main()
