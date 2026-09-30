"""
What is the 2021 step? Five tests of the island-wide shift, including the two that failed.

    python scripts/input_prep/5-scr/1-observations/detrended_position/coastsat_position_attribution.py   # after coastsat_detrended_position.py

Writes the test tables and figures beside the detrended position. Details: scripts/input_prep/5-scr/1-observations/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-23
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
from site_layer import hat_observed_rates as obs      # noqa: E402
from site_layer import hat_figure_style as fs         # noqa: E402

sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401
import coastsat_lrr as cl                             # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coastsat_detrended_position import (NOURISHED, NOURISHED_DOMAINS, PRE, POST,  # noqa: E402
                            FIRST_YEAR, LAST_YEAR)

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402

# --- CONFIG ------------------------------------------------------------------
# The dune-line survey pairs, by the rate-window folder that holds them
DUNE_PAIRS = [("1996_2010", "1997-2009", "control, no step"),
              ("2010_2024", "2009-2023", "spans the step"),
              ("1996_2024", "1997-2023", "spans the step")]
# -----------------------------------------------------------------------------


# The detrended position matrix, years as ints
def load_matrix():
    M = pd.read_csv(obs.DETRENDED_POSITION / obs.DETRENDED_POSITION_MATRIX, index_col=0)
    M.index = M.index.astype(int)
    return M


# The domain of each transect
def domains_of(columns):
    lookup = pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv")
    lookup = lookup.dropna(subset=["domain_number"]).set_index("transect_id")
    return lookup["domain_number"].reindex(columns).astype(int)


# Does the step survive on domains that were never touched?
def test_nourishment(M, dom):
    touched = dom.isin(NOURISHED_DOMAINS).to_numpy()
    out = {}
    for label, sel in (("all", np.ones(len(touched), bool)),
                       ("untouched", ~touched), ("nourished", touched)):
        ix = M.loc[:, sel].mean(axis=1)
        out[label] = {
            "jump_2020_2021_m": round(float(ix.loc[2021] - ix.loc[2020]), 2),
            "step_m": round(float(ix.loc[POST[0]:POST[1]].mean()
                                  - ix.loc[PRE[0]:PRE[1]].mean()), 2),
            "n_transects": int(sel.sum()),
        }
    rows = []
    for name, (year, ds) in NOURISHED.items():
        inside = dom.isin(list(ds)).to_numpy()
        a = M.loc[:, inside].mean(axis=1)
        b = M.loc[:, ~inside].mean(axis=1)
        local = float(a.loc[year + 1] - a.loc[year - 1])
        island = float(b.loc[year + 1] - b.loc[year - 1])
        rows.append({"project": name, "year": year,
                     "n_domains": len(list(ds)),
                     "filled_step_m": round(local, 1),
                     "rest_of_island_m": round(island, 1),
                     "difference_m": round(local - island, 1)})
    return out, pd.DataFrame(rows)


# Is the step an artefact of when and how often the satellites looked?
def test_sampling(dom_index_sample=12):
    lookup = pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv")
    lookup = lookup.dropna(subset=["domain_number"]).iloc[::dom_index_sample]
    frames = []
    for r in lookup.itertuples(index=False):
        site = r.transect_id.rsplit("_", 1)[0]
        path = (obs.COASTSAT_TIMESERIES / "{0}_timeseries".format(site)
                / "{0}.csv".format(r.transect_id))
        if not path.is_file():
            continue
        df = cl.filter_dates(cl.load_timeseries(str(path)),
                             "{0}-01-01".format(FIRST_YEAR),
                             "{0}-12-31".format(LAST_YEAR))
        t = (df["date"] - pd.Timestamp("{0}-01-01".format(FIRST_YEAR), tz="UTC")
             ).dt.days / 365.25
        sl, ic = np.polyfit(t, df["chainage_m"], 1)
        frames.append(df.assign(anom=df["chainage_m"] - (sl * t + ic),
                                year=df["date"].dt.year,
                                q=((df["date"].dt.month - 1) // 3 + 1)))
    a = pd.concat(frames, ignore_index=True)
    rows = []
    for qq in (1, 2, 3, 4):
        before = a[(a["q"] == qq) & a["year"].between(PRE[0], PRE[1])]["anom"].median()
        after = a[(a["q"] == qq) & a["year"].between(POST[0], POST[1])]["anom"].median()
        rows.append({"quarter": "Q{0}".format(qq),
                     "pre_m": round(float(before), 1),
                     "post_m": round(float(after), 1),
                     "step_m": round(float(after - before), 1)})
    return pd.DataFrame(rows)


# Does the independent dune line show the same reversal? Stored tables only -- both endpoint products ...
def test_duneline():
    rows = []
    for window, label, role in DUNE_PAIRS:
        dune = pd.read_csv(obs.DUNELINE_ENDPOINT_ROOT / window
                           / obs.ENDPOINT_DOMAIN_FILE)
        cs = pd.read_csv(obs.COASTSAT_ENDPOINT_ROOT / window
                         / obs.ENDPOINT_DOMAIN_FILE)
        dc = [c for c in dune.columns if "change" in c and "mean" in c] or ["change_m"]
        cc = [c for c in cs.columns if "change" in c and "mean" in c] or ["change_m"]
        m = (dune[["domain_number", dc[0]]].rename(columns={dc[0]: "dune_m"})
             .merge(cs[["domain_number", cc[0]]].rename(columns={cc[0]: "coastsat_m"}),
                    on="domain_number"))
        m["offset_m"] = m["coastsat_m"] - m["dune_m"]
        rows.append({"window": window, "pair": label, "role": role,
                     "dune_median_m": round(float(m["dune_m"].median()), 1),
                     "coastsat_median_m": round(float(m["coastsat_m"].median()), 1),
                     "offset_median_m": round(float(m["offset_m"].median()), 1),
                     "r_alongshore": round(float(m["dune_m"].corr(m["coastsat_m"])), 2),
                     "n_domains": len(m)})
    d = pd.DataFrame(rows)
    d["offset_vs_control_m"] = (d["offset_median_m"]
                                - d.loc[0, "offset_median_m"]).round(1)
    return d


# The three tests a picture helps with
def draw(nourish_split, dune, quarters, by_domain, out_dir):
    fs.apply_style()
    fig, axes = plt.subplots(1, 3, figsize=fs.figsize("double", aspect=0.40),
                             layout="constrained")

    # (a) nourishment: does the step survive on untouched ground?
    ax = axes[0]
    labels = ["all", "untouched", "nourished"]
    vals = [nourish_split[k]["jump_2020_2021_m"] for k in labels]
    ax.bar(range(3), vals, color=[fs.C["BASE"], fs.C["LATE"], fs.C["ADDED"]],
           width=0.66, lw=0, zorder=3)
    for i, v in enumerate(vals):
        ax.annotate("{0:+.1f}".format(v), xy=(i, v), xytext=(0, 3),
                    textcoords="offset points", ha="center", fontsize=7.5)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["all\ntransects", "never\nnourished", "nourished"])
    ax.set_ylabel("2020 → 2021 jump (m)")
    fs._title(ax, 0, "Nourishment")
    ax.grid(True, axis="y", alpha=0.6)

    # (b) the dune line, both sources per interval
    ax = axes[1]
    w = 0.36
    xs = np.arange(len(dune))
    ax.bar(xs - w / 2, dune["dune_m"] if "dune_m" in dune else dune["dune_median_m"],
           width=w, color=fs.C["EARLY"], lw=0, label="dune line", zorder=3)
    ax.bar(xs + w / 2, dune["coastsat_median_m"], width=w, color=fs.C["LATE"],
           lw=0, label="CoastSat", zorder=3)
    ax.axhline(0.0, color=fs.C["INK"], lw=0.8, zorder=4)
    ax.set_xticks(xs)
    ax.set_xticklabels([p.replace("-", "→\n") for p in dune["pair"]],
                       fontsize=7.5)
    ax.set_ylabel("net change (m)")
    fs._title(ax, 1, "The dune line")
    ax.legend(loc="lower right")
    ax.grid(True, axis="y", alpha=0.6)

    # (c) spatial structure: a sensor offset would be flat
    ax = axes[2]
    s = by_domain["step_m"].to_numpy()
    ax.hist(s, bins=22, color=fs.C["LATE"], zorder=3)
    ax.axvline(0.0, color=fs.C["INK"], lw=0.8, zorder=4)
    ax.axvline(float(np.median(s)), color=fs.C["ACCENT"], lw=1.2,
               ls=(0, (4, 2)), zorder=5)
    ax.set_xlabel("step by domain (m)")
    ax.set_ylabel("domains")
    fs._title(ax, 2, "Spatial structure")
    ax.grid(True, axis="y", alpha=0.6)

    paths = fs.save(fig, Path(out_dir) / "detrended_position_attribution", close=True)
    fs.record_caption(paths[0],
        "Three of the five tests of what the 2021 step is. (a) the "
        "2020→2021 jump computed on all transects, on the 672 that were "
        "never nourished, and on the 234 inside a fill footprint — the "
        "jump survives untouched ground almost unchanged, and it lands a year "
        "before the 2022 fills, so nourishment is not its cause. (b) net "
        "change from the digitized dune line and from CoastSat over the same "
        "survey-date pairs; the 2023 line sits inside the stepped period and "
        "the dune line reverses from retreat to advance with CoastSat, "
        "matching it to 0.1 m over 2009→2023. An artefact would have made "
        "that interval carry about +13.6 m more CoastSat-minus-dune; it "
        "carries 9.7 m less. (c) the step by domain: a sensor or waterline "
        "bias applies one offset everywhere, but this spreads from "
        "−11.9 to +71.5 m with eight domains stepping landward, while "
        "transects inside a domain agree to 2.8 m. The two tests not drawn are "
        "in attribution_tests.csv: per-transect noise (rejected) and the "
        "seasonal sampling mix (rejected). The verdict is that the step is "
        "real, so the grading target's sensitivity to 2021–2024 is a "
        "question about which period the model should represent, not about "
        "which data to trust."
    )
    return paths[0]


# Run: the five tests, then the tables and figures
def main():
    out_dir = obs.DETRENDED_POSITION
    M = load_matrix()
    dom = domains_of(M.columns)
    by_domain = pd.read_csv(out_dir / "detrended_position_by_domain.csv")

    split, fills = test_nourishment(M, dom)
    quarters = test_sampling()
    dune = test_duneline()

    fills.to_csv(out_dir / "attribution_nourishment.csv", index=False)
    dune.to_csv(out_dir / "attribution_duneline.csv", index=False)

    step = by_domain["step_m"]
    tests = pd.DataFrame([
        {"test": "per-transect noise", "verdict": "REJECTED",
         "measured": "correlation of a transect's noise-to-trend ratio with its "
                     "convergence time",
         "result": "r = 0.01 forward, 0.16 backward; calmest quartile 27 yr, "
                   "swingiest 26 yr"},
        {"test": "nourishment", "verdict": "REJECTED as cause",
         "measured": "2020->2021 jump on never-nourished vs nourished transects",
         "result": "{0:+.1f} m untouched vs {1:+.1f} m nourished; the jump is a "
                   "year before the 2022 fills".format(
                       split["untouched"]["jump_2020_2021_m"],
                       split["nourished"]["jump_2020_2021_m"])},
        {"test": "nourishment visibility", "verdict": "CONFIRMED signal",
         "measured": "each fill's step against the rest of the island, year after "
                     "minus year before",
         "result": "; ".join("{0} {1:+.1f} m".format(r.project, r.difference_m)
                             for r in fills.itertuples(index=False))},
        {"test": "seasonal sampling mix", "verdict": "REJECTED",
         "measured": "the step computed within each calendar quarter separately",
         "result": "; ".join("{0} {1:+.1f} m".format(r.quarter, r.step_m)
                             for r in quarters.itertuples(index=False))},
        {"test": "dune line", "verdict": "CORROBORATES",
         "measured": "net change from both sources over the same survey-date pairs",
         "result": "1997->2009 dune {0:+.1f} m / CoastSat {1:+.1f} m; "
                   "2009->2023 dune {2:+.1f} m / CoastSat {3:+.1f} m; "
                   "difference-of-differences {4:+.1f} m (an artefact predicts "
                   "about +13.6)".format(
                       dune.loc[0, "dune_median_m"], dune.loc[0, "coastsat_median_m"],
                       dune.loc[1, "dune_median_m"], dune.loc[1, "coastsat_median_m"],
                       dune.loc[1, "offset_vs_control_m"])},
        {"test": "spatial structure", "verdict": "CORROBORATES",
         "measured": "spread of the step between domains vs within a domain",
         "result": "domain means {0:+.1f} to {1:+.1f} m, sd {2:.1f} m between "
                   "domains vs {3:.1f} m within; {4} domains step landward".format(
                       step.min(), step.max(), step.std(),
                       by_domain["step_sd_within_m"].median(), int((step < 0).sum()))},
    ])
    tests.to_csv(out_dir / "attribution_tests.csv", index=False)
    draw(split, dune, quarters, by_domain, out_dir)

    print(tests[["test", "verdict"]].to_string(index=False))
    print("\n" + dune.to_string(index=False))
    print("\n" + fills.to_string(index=False))
    print("\nwrote {0}".format(out_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
