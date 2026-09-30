r"""
HAT_excess_overwash_diagnosis.py -- why does the model overwash more than the imagery shows?
==============================================================================
THE QUESTION (Hannah, 2026-09-28). Every storm series, the committed one
included, overwashes far more domains than the observed record in most image
windows (storms-and-overwash/2026-09-28-storm-length-selection, stage 1).

THREE EXPLANATIONS, EACH TESTED DIRECTLY (read-only: existing runs and inputs)
    1  volume   the extra overwash is real but too small to see in imagery
    2  dunes    the model's dune row is lower than the real foredune, either
                from the extraction (a clipped search window leaves the true
                crest in interior row 0, behind a lower dune row) or from the
                dunes changing during the run
    3  storms   Rhigh (Stockdon R2 on WIS Hs, slope 0.06) clears the dunes by
                too much

    `hidden`  per domain, the start-of-run gap between the dune row's crest
              and the highest of the dune row + first N interior rows: the
              foredune height the dune row does not carry
    `cells`   per image x domain cell (managed runs, drop72 and trim24): the
              storms credited with its overwash, their Rhigh, the pre-storm
              crest (after growth, as Barrier3D tests it), the margin, the
              volume, and the domain's hidden crest
    `summary` how false alarms and hits differ on each of those

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-excess-overwash-diagnosis/
==============================================================================

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))
import HAT_storm_length_selection as S  # noqa: E402

TAG = "storms-and-overwash/2026-09-28-excess-overwash-diagnosis"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
N_ROWS = 10          # interior rows (100 m) behind the dune row searched for a higher crest
DAM = 10.0
VARIANTS = ("drop72", "trim24")


def run_dir(w, v):
    if v == "drop72":
        return S.MATRIX / S.wtag(w) / "edgeBE" / S.CONTROLS[(w[0], "full_management")]
    return S.run_path(v, "full_management", w)


def hidden():
    """Start-of-run dune-row crest vs the highest ground in the dune row plus
    the first N_ROWS interior rows, per domain (m MHW, alongshore medians)."""
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rows = []
    for w in S.WINDOWS:
        c = S.load_state(run_dir(w, "drop72"))
        for gis in range(1, 91):
            b = c.barrier3d[DOM.gis_to_pad(gis)]
            dune = (np.asarray(b.DuneDomain[0]).max(axis=1) + b.BermEl) * DAM         # (50,) m MHW
            inter = np.asarray(b.DomainTS[0])[:N_ROWS] * DAM                           # (N, 50) m MHW
            ridge = np.maximum(dune, inter.max(axis=0))
            rows.append(dict(window=S.wtag(w), gis=gis, dune_crest_m=float(np.median(dune)),
                             ridge_crest_m=float(np.median(ridge)),
                             hidden_m=float(np.median(ridge - dune)),
                             hidden_share=float(np.mean(ridge - dune > 0.25)),
                             dune_min_m=float(dune.min())))
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "hidden_crest.csv", index=False)
    for w, g in t.groupby("window"):
        print(f"{w}: median dune-row crest {g.dune_crest_m.median():.2f} m MHW, median ridge within "
              f"{N_ROWS * 10} m {g.ridge_crest_m.median():.2f} m; domains with >0.5 m hidden: "
              f"{(g.hidden_m > 0.5).sum()}, >1 m: {(g.hidden_m > 1).sum()}")
    return t


def storm_table(c, summ):
    """Per (pad, storm): Rhigh, pre-storm crest (after growth), margin, and the
    overwash share credited to it (S.storm_shares)."""
    shares = S.storm_shares(c, summ)
    rows = []
    for (p, i), vol in shares.items():
        b = c.barrier3d[p]
        t = int(summ.at[i, "time"])
        dd = np.array(b.DuneDomain, dtype=float, copy=True)
        dd, _, _ = type(b).DuneGrowth(b, dd, t)
        crest = dd[t].max(axis=1)
        crest[crest < b._DuneRestart] = b._DuneRestart
        el = (crest + b._BermEl) * DAM
        rh = summ.at[i, "Rhigh"] * DAM
        rows.append(dict(pad_index=p, storm=i, t=t, rhigh_m=rh, crest_min_m=float(el.min()),
                         crest_mean_m=float(el.mean()), margin_min_m=rh - float(el.min()),
                         margin_mean_m=rh - float(el.mean()), cells_over=int((el < rh).sum()),
                         volume_m3_per_m=vol))
    return pd.DataFrame(rows)


def cells():
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    ovm = S.overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    hid = pd.read_csv(EXP_DIR / "tables" / "hidden_crest.csv")
    out = []
    for w in S.WINDOWS:
        for v in VARIANTS:
            c = S.load_state(run_dir(w, v))
            summ = pd.read_csv(S.summary_csv(w, v))
            st = storm_table(c, summ)
            dated = pd.to_datetime(summ.EndTime) - ovm.GRACE
            for _, im in ovm.images_in(w, obs).iterrows():
                o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
                in_win = set(summ.index[(dated > im["from"]) & (dated <= im.date)])
                maxrh = summ.loc[list(in_win), "Rhigh"].max() * DAM if in_win else np.nan
                for gis in range(1, 91):
                    ob = o.get(gis, np.nan)
                    if np.isnan(ob):
                        continue
                    p = DOM.gis_to_pad(gis)
                    s = st[(st.pad_index == p) & st.storm.isin(list(in_win))]
                    h = hid[(hid.window == S.wtag(w)) & (hid.gis == gis)].iloc[0]
                    out.append(dict(window=S.wtag(w), variant=v, obs_id=im.Obs_ID, image=im.date.date(), gis=gis,
                                    observed=int(ob), model=int(s.volume_m3_per_m.sum() > 0),
                                    volume_m3_per_m=float(s.volume_m3_per_m.sum()),
                                    n_storms=len(s), max_margin_min_m=float(s.margin_min_m.max()) if len(s) else np.nan,
                                    max_margin_mean_m=float(s.margin_mean_m.max()) if len(s) else np.nan,
                                    crest_min_m=float(s.crest_min_m.min()) if len(s) else np.nan,
                                    window_max_rhigh_m=maxrh, hidden_m=h.hidden_m,
                                    start_dune_crest_m=h.dune_crest_m, start_ridge_m=h.ridge_crest_m))
            print(f"  cells {S.wtag(w)} {v}", flush=True)
    t = pd.DataFrame(out)
    t.to_csv(EXP_DIR / "tables" / "cells.csv", index=False)
    return t


def summary():
    t = pd.read_csv(EXP_DIR / "tables" / "cells.csv")
    t["kind"] = np.select([(t.observed == 1) & (t.model == 1), (t.observed == 0) & (t.model == 1),
                           (t.observed == 1) & (t.model == 0)], ["hit", "false alarm", "miss"], "correct no")
    cols = ["volume_m3_per_m", "max_margin_min_m", "max_margin_mean_m", "crest_min_m", "hidden_m",
            "start_dune_crest_m", "window_max_rhigh_m"]
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(t.groupby(["window", "variant", "kind"])[cols].median().round(2).to_string())
        # false-alarm rate by hidden-crest class, per window
        t["hidden_class"] = pd.cut(t.hidden_m, [-1, 0.25, 0.5, 1.0, 10], labels=["<0.25", "0.25-0.5", "0.5-1", ">1"])
        neg = t[t.observed == 0]
        print("\nfalse-alarm rate by hidden foredune height (observed-absent cells):")
        print(neg.groupby(["window", "variant", "hidden_class"], observed=True).model.agg(["mean", "size"]).round(2).to_string())
    t.to_csv(EXP_DIR / "tables" / "cells.csv", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["hidden", "cells", "summary", "all"])
    a = ap.parse_args()
    if a.action in ("hidden", "all"):
        hidden()
    if a.action in ("cells", "all"):
        cells()
    if a.action in ("summary", "all"):
        summary()


if __name__ == "__main__":
    main()
