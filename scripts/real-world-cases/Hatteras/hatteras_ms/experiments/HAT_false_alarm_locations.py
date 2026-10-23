"""
Where does slower dune growth add its false alarms: on dunes the surveys show were low, or elsewhere?

    python scripts/hatteras_ms/experiments/HAT_false_alarm_locations.py

Reads two 1996-2009 runs of the low growth-rate study (no growth, 5.0 m, and the
control, 5.5 m with 0.55-0.95), scores each image and domain as overwash_vs_model
does, and sets every model overwash against the surveyed dune: the lowest crest
of the domain in the last survey flown before the image (1996-start/survey-crests,
1996 dune windows), compared with the highest storm run-up between the image and
the one before it. A false alarm on a dune the surveys put below that run-up is a
place reality should have overtopped and the imagery says did not (the storm,
overwash or imagery side); one on a dune the surveys put above it is a place the
model's dune is too low (the dune side). Read-only.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_growth_low_rates as L  # noqa: E402
import HAT_dune_crest_vs_surveys as V  # noqa: E402
G, U = L.G, L.U
import overwash_vs_model as ovm  # noqa: E402
from site_layer import hat_overwash as ow  # noqa: E402
from site_layer import hat_env_forcings as env  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-overwash/2026-10-19-false-alarm-locations"
EXP_DIR = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
RUNS = {"no growth (5.0 m)": G.name(5.0, (0.0, 0.0)), "control (5.5 m, 0.55-0.95)": G.name(*G.CONTROL)}
START, FIRST = 1996, 1997          # the window, and the year DuneDomain[0] is dated 1 January of
MHW = 0.36
# -----------------------------------------------------------------------------


# Every image x domain of the window for one run: observed, modelled (overwash_vs_model's rule)
def cells(variant):
    c = U.load(G.run_dir(START, variant))
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    orig = ovm.load_qow
    ovm.load_qow = lambda w, arm: (q, variant)
    try:
        obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
        d = ovm.compare((START, U.PERIODS[START][0]), "managed", obs, ovm.HEADLINE_THR).dropna(subset=["observed"])
    finally:
        ovm.load_qow = orig
    mod = V.model_crests(START, variant).reshape(-1, 90, 50)          # (states, domains, cells), m NAVD88
    return d, mod


# Highest run-up between each image and the one before, m NAVD88, and the year of that storm
def runup(d):
    s = pd.read_csv(env.storm_summary_file(*ovm.model_years((START, U.PERIODS[START][0]))),
                    parse_dates=["StartTime", "EndTime"])
    s["dated"] = s.EndTime - ovm.GRACE
    out = {}
    for (since, image), _ in d.groupby(["since", "image"]):
        w = s[(s.dated > pd.Timestamp(since)) & (s.dated <= pd.Timestamp(image))]
        out[(since, image)] = (w.Rhigh.max() * 10 + MHW, int(w.loc[w.Rhigh.idxmax(), "calendar_year"])) \
            if len(w) else (np.nan, np.nan)
    return out


# The lowest surveyed crest of each domain in the last survey flown before a date, m NAVD88
def surveyed_low(sv, date):
    flown = sv.groupby("survey").flown.first().str[:7].map(lambda x: pd.Timestamp(x + "-01"))
    before = flown[flown <= pd.Timestamp(date)]
    if before.empty:
        return None, None
    year = before.idxmax()
    g = sv[sv.survey == year]
    return g.groupby("gis").crest.min(), year


def classify(d, mod, ru, sv):
    rows = []
    for _, r in d.iterrows():
        peak, yr = ru[(r.since, r.image)]
        low, year = surveyed_low(sv, r.since)
        if low is None or np.isnan(peak):
            continue
        k = max(int(yr) - FIRST, 0)                                    # the model state before that storm
        rows.append(dict(image=r.image, since=r.since, gis=r.gis, observed=int(r.observed), model=int(r.model),
                         runup_m_navd88=round(peak, 2), survey_used=year,
                         surveyed_low_crest=round(float(low.get(r.gis, np.nan)), 2),
                         model_low_crest=round(float(mod[k, r.gis - 1].min()), 2)))
    t = pd.DataFrame(rows)
    t["surveyed_dune_below_runup"] = t.surveyed_low_crest < t.runup_m_navd88
    t["outcome"] = np.select([(t.observed == 1) & (t.model == 1), (t.observed == 1) & (t.model == 0),
                              (t.observed == 0) & (t.model == 1)], ["hit", "miss", "false_alarm"], "neither")
    return t


def main():
    sv = V.surveys(START)
    rows, summ = [], []
    for label, variant in RUNS.items():
        d, mod = cells(variant)
        t = classify(d, mod, runup(d), sv).assign(run=label)
        rows.append(t)
        for outcome in ("hit", "miss", "false_alarm", "neither"):
            o = t[t.outcome == outcome]
            summ.append(dict(run=label, outcome=outcome, cells=len(o),
                             on_surveyed_dune_below_runup=int(o.surveyed_dune_below_runup.sum()),
                             on_surveyed_dune_above_runup=int((~o.surveyed_dune_below_runup).sum()),
                             median_surveyed_low_crest=o.surveyed_low_crest.median(),
                             median_model_low_crest=o.model_low_crest.median(),
                             median_runup=o.runup_m_navd88.median()))
    t, s = pd.concat(rows, ignore_index=True), pd.DataFrame(summ).round(2)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "cells_with_surveyed_dunes.csv", index=False)
    s.to_csv(EXP_DIR / "tables" / "outcomes_by_surveyed_dune.csv", index=False)

    # The false alarms no growth adds over the control, by image and domain
    key = ["image", "gis"]
    fa = {lab: set(map(tuple, t[(t.run == lab) & (t.outcome == "false_alarm")][key].astype(str).values)) for lab in RUNS}
    a, b = list(RUNS)
    added = t[(t.run == a) & t[key].astype(str).apply(tuple, axis=1).isin(fa[a] - fa[b])]
    by_img = added.groupby("image").agg(added=("gis", "size"), runup=("runup_m_navd88", "first"),
                                        below=("surveyed_dune_below_runup", "sum"),
                                        median_surveyed_low=("surveyed_low_crest", "median"),
                                        median_model_low=("model_low_crest", "median"))
    by_img.to_csv(EXP_DIR / "tables" / "added_false_alarms_by_image.csv")
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(s.to_string(index=False))
        print(f"\nfalse alarms added by no growth: {len(added)}; on a surveyed dune below the run-up: "
              f"{int(added.surveyed_dune_below_runup.sum())}")
        print(by_img.to_string())
        print("\nadded false alarms by domain (top):", added.gis.value_counts().head(12).to_dict())


if __name__ == "__main__":
    main()
