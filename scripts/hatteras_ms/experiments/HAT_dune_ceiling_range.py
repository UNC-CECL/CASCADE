"""
Across a range of island-wide dune ceilings, how much overwash does the model make, and is it where the record puts it?

    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_range.py run
    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_range.py score

Sets Dmaxel in-process to 4.5-6.0 m NAVD88 (5.5 is the parameter file as
adopted, run unchanged); managed edgeBE runs of both periods. Scores overall
overwash skill and, per GIS domain, observed against modelled overwash
frequency, by zone, with an alongshore figure. Details:
scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-overwash-and-dunes/2026-10-13-dune-ceiling-range"
CEILINGS = {"c450": 4.5, "c500": 5.0, "c525": 5.25, "c550": None, "c575": 5.75, "c600": 6.0,   # m NAVD88
            "c840": 8.4, "c1210": 12.1,   # 2009 lidar crest P95 and maximum: "the DEM maximum" as a ceiling
            "c539": 5.39, "c573": 5.73}   # each period's own starting-DEM median crest (2009, 1996)
ADOPTED_M = 5.5               # what c550 runs with, from the parameter file
ZONES = {"Cape Point": (1, 6), "communities": None, "road": None}
COMMUNITY_GIS = set(range(7, 9)) | set(range(21, 32)) | set(range(68, 84))
# -----------------------------------------------------------------------------


def zone(gis):
    return "Cape Point" if gis <= 6 else ("communities" if gis in COMMUNITY_GIS else "road")


def ceiling(v):
    return ADOPTED_M if CEILINGS[v] is None else CEILINGS[v]


# Child process: set the island-wide ceiling, then run the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    value = CEILINGS[variant]
    if value is not None:
        orig = H.build_cascade

        def build_cascade(*a, **k):
            c = orig(*a, **k)
            U.set_uniform(c, value)
            print(f"EXPERIMENT dune ceiling = {value:.2f} m NAVD88 (parameter file's 5.5 overridden)", flush=True)
            return c

        H.build_cascade = build_cascade
    sys.path.insert(0, str(U.HINDCAST.parent))
    sys.argv = [str(U.HINDCAST)]
    runpy.run_path(str(U.HINDCAST), run_name="__main__")


def adopt():
    U.TAG, U.EXP_DIR, U._HERE = TAG, U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG, _HERE
    U.VARIANTS = dict(CEILINGS)


# Observed and modelled overwash per assessed image and domain, overwash_vs_model's rule
def overwash_cells(c, obs, dom, storm_years):
    import overwash_vs_model as ovm
    from site_layer import hat_env_forcings as env
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    s = pd.read_csv(env.storm_summary_file(*storm_years), parse_dates=["StartTime", "EndTime"])
    s["dated"] = s.EndTime - ovm.GRACE
    largest = s.loc[s.groupby("time").Rhigh.idxmax()].set_index("time").dated
    imgs = ovm.images_in((storm_years[0], storm_years[1] + 1), obs)
    rows = []
    for _, im in imgs.iterrows():
        o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
        for gis in range(1, 91):
            p = dom.gis_to_pad(gis)
            yrs = [t for t in range(1, q.shape[0]) if q[t, p] > ovm.HEADLINE_THR]
            m = any(im["from"] < largest[t] <= im.date for t in yrs if t in largest)
            rows.append(dict(obs_id=im.Obs_ID, gis=gis, observed=o.get(gis, np.nan), model=int(m)))
    d = pd.DataFrame(rows).dropna(subset=["observed"])
    o, m = d.observed.astype(int), d.model.astype(int)
    d["outcome"] = np.select([(o == 1) & (m == 1), (o == 1) & (m == 0), (o == 0) & (m == 1)],
                             ["hit", "miss", "false_alarm"], "neither")
    d["zone"] = d.gis.map(zone)
    return d


def score():
    import storms_vs_overwash as SV
    from site_layer import hat_overwash as ow
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    rows, per_dom = [], []
    for start, (end, storm_years) in U.PERIODS.items():
        for v in CEILINGS:
            rd = U.run_dir(start, v)
            if rd is None:
                print(f"  missing {start} {v}")
                continue
            c = U.load(rd)
            d = overwash_cells(c, obs, DOM, storm_years)
            g = d.groupby("gis").agg(observed=("observed", "mean"), modelled=("model", "mean"), images=("obs_id", "size"))
            g = g.reset_index().assign(period=f"{start}_{end}", ceiling_m_navd88=ceiling(v))
            per_dom.append(g)
            sk = SV.skill(d)
            row = dict(period=f"{start}_{end}", ceiling_m_navd88=ceiling(v),
                       skill=round(sk["skill"], 3), hit_rate=round(sk["hit_rate"], 3),
                       false_alarm_rate=round(sk["false_alarm_rate"], 3),
                       hits=sk["hits"], misses=sk["misses"], false_alarms=sk["false_alarms"],
                       observed_share=round(float(d.observed.mean()), 3), modelled_share=round(float(d.model.mean()), 3),
                       space_r=round(float(np.corrcoef(g.observed, g.modelled)[0, 1]), 2) if g.modelled.std() > 0 else np.nan)
            for z in ("Cape Point", "communities", "road"):
                n = d[d.zone == z].outcome.value_counts()
                row[f"{z}_hit_miss_fa"] = f"{n.get('hit', 0)}/{n.get('miss', 0)}/{n.get('false_alarm', 0)}"
            meta = __import__("json").loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
            row["dune_ceiling_recorded"] = meta.get("identity", {}).get("dune_ceiling", "")
            row["rmse_interior_m_yr"] = float(meta["skill"]["rmse_interior_m_yr"])
            rows.append(row)
            print(f"  scored {start} {v}", flush=True)
    t = pd.DataFrame(rows)
    pdom = pd.concat(per_dom, ignore_index=True)
    tables = U.EXP_DIR / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    t.to_csv(tables / "scores.csv", index=False)
    pdom.to_csv(tables / "overwash_by_domain.csv", index=False)
    with pd.option_context("display.width", 320, "display.max_columns", 40):
        print(t.to_string(index=False))
    figure(pdom)


# Observed against modelled overwash frequency along the island, one panel per period
def figure(pdom):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import DOMAIN_AXIS_LABEL, INK, apply_style, figsize, save, town_bands
    apply_style()
    periods = list(dict.fromkeys(pdom.period))
    fig, axes = plt.subplots(len(periods), 1, figsize=figsize("double", height=2.6 * len(periods)), sharex=True)
    cmap = plt.get_cmap("viridis")
    vals = sorted(pdom.ceiling_m_navd88.unique())
    for ax, per in zip(np.atleast_1d(axes), periods):
        p = pdom[pdom.period == per]
        obs = p[p.ceiling_m_navd88 == vals[0]].sort_values("gis")
        ax.bar(obs.gis, obs.observed, color="0.8", width=0.9, label="observed")
        for i, v in enumerate(vals):
            m = p[p.ceiling_m_navd88 == v].sort_values("gis")
            ax.plot(m.gis, m.modelled, color=cmap(i / max(len(vals) - 1, 1)), lw=1.2 if v != ADOPTED_M else 2.2,
                    label=f"{v:.2f} m" + (" (adopted)" if v == ADOPTED_M else ""))
        ax.set_ylabel("share of images\nwith overwash")
        ax.set_ylim(0, 1.05)
        ax.set_title(per.replace("_", "–"), loc="left", color=INK)
        town_bands(ax)
    np.atleast_1d(axes)[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    np.atleast_1d(axes)[0].legend(ncol=4, fontsize=7, frameon=False, loc="upper right", title="dune ceiling, m NAVD88")
    out = U.EXP_DIR / "figures" / "overwash_by_domain_dune_ceiling_range.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    save(fig, out, close=True)
    print(f"  figure: {out}")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2])
        return
    adopt()
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        from concurrent.futures import ThreadPoolExecutor
        jobs = [(s, v) for s in U.PERIODS for v in CEILINGS if U.run_dir(s, v) is None]
        with ThreadPoolExecutor(U.WORKERS) as ex:
            list(ex.map(U.launch, jobs))
    else:
        score()


if __name__ == "__main__":
    main()
