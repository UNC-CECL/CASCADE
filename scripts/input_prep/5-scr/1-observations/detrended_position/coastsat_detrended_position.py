"""
coastsat_detrended_position.py
==============================================================================
IS THERE ONE SIGNAL THE WHOLE ISLAND SHARES, ON TOP OF EACH TRANSECT'S TREND?

Detrend every CoastSat transect against its OWN 1996-2024 fit, reduce each to
one value a year (the annual median position), and average across all 906. Any
signal that survives that averaging is common to the island, because anything
local is incoherent between transects and cancels.

WHY IT WAS BUILT (2026-09-23). `3-rates/coastsat/window_convergence/` found
that no window shorter than about 25 years recovers the long-term rate, and
that the answer barely varies between transects -- a transect with a fast clean
trend needs as long as a slow noisy one. A per-transect explanation cannot
produce a per-transect-invariant answer, so the cause had to be shared. This is
the search for it.

WHAT IT FINDS. One coherent excursion: roughly flat 1996-2004, a landward sag
through 2005-2020 bottoming at -7.4 m, then +16.8 m in a single year into 2021,
held through 2024. It is only 17% of the mean transect variance, but it is the
COHERENT part, and coherent is what moves an OLS slope systematically:

    window        bias the departure alone puts on the fitted rate
    1996-2010     -0.367 m/yr      (the graded window, ~37% of a median rate)
    2010-2024     +0.875 m/yr      (the chain's second leg, opposite sign)
    1996-2024      0.000 m/yr      (zero by construction -- it is detrended
                                    against this window)

That is the answer to "why do the windows disagree": not noise, a shared
multi-decadal excursion that an OLS slope cannot separate from the trend until
the window spans the whole of it. `coastsat_position_attribution.py` beside this file
then asks what the excursion IS.

WHAT THIS IS NOT. Not a rate product and not a model input. Nothing is graded
against it; it exists to explain a property of the rate products, which is why
it lives under 1-observations rather than 3-rates.

Inputs
------
    transect_domain_lookup.csv      2-transect-frame/, via hat_observed_rates
    CoastSat time-series CSVs       1-observations/coastsat_timeseries/

Outputs  (hat_observed_rates.DETRENDED_POSITION)
--------------------------------------------
    annual_medians_detrended.csv    year x transect, m. The matrix everything
                                    else here reads, so the detrending is done
                                    once and cannot drift between scripts
    detrended_position_by_year.csv        a row per year: the index, the counts, and
                                    the nourished / untouched split
    detrended_position_by_domain.csv    a row per domain: its share of the step,
                                    and how well it tracks the index
    detrended_position.png              the index, the alongshore step, the record

Usage
-----
    python .../coastsat_detrended_position.py
==============================================================================
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

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                   # noqa: E402

# ============================================================
# CONFIG
# ============================================================

FIRST_YEAR, LAST_YEAR = 1996, 2024
MIN_YEARS = 25              # a transect needs this many years to be detrended

# The step the index turns out to contain, as two periods to difference. Named
# here rather than buried so the number quoted in the README and the figure is
# the same one.
PRE = (2016, 2020)
POST = (2021, 2024)

# NOURISHMENT, as the hindcast receives it (hatteras_site_config
# .HATTERAS_NOURISHMENT_PROJECTS, drawn in 4-mgmt-forcing/nourishment/). Held
# here so the index can be split by it; `coastsat_position_attribution.py` does the test.
NOURISHED = {"Rodanthe 2014": (2014, range(84, 90)),
             "Buxton 2022": (2022, range(6, 16)),
             "Avon 2022": (2022, range(21, 29))}
NOURISHED_DOMAINS = sorted({d for _y, ds in NOURISHED.values() for d in ds})


# ============================================================
# THE MATRIX
# ============================================================

def build_matrix():
    """Detrended annual median position, year x transect, in metres.

    Each transect is detrended against ITS OWN 1996-2024 OLS, so what is left
    is departure from that transect's long-term behaviour and nothing else. A
    transect eroding at 3 m/yr and one accreting at 1 m/yr both come out
    centred on zero, which is what lets them be averaged.
    """
    lookup = pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv")
    lookup = lookup.dropna(subset=["domain_number"]).sort_values(
        ["domain_number", "transect_id"])
    years = np.arange(FIRST_YEAR, LAST_YEAR + 1)

    cols, counts, domains = {}, {}, {}
    for r in lookup.itertuples(index=False):
        site = r.transect_id.rsplit("_", 1)[0]
        path = (obs.COASTSAT_TIMESERIES / "{0}_timeseries".format(site)
                / "{0}.csv".format(r.transect_id))
        if not path.is_file():
            continue
        df = cl.filter_dates(cl.load_timeseries(str(path)),
                             "{0}-01-01".format(FIRST_YEAR),
                             "{0}-12-31".format(LAST_YEAR))
        by_year = df.groupby(df["date"].dt.year)["chainage_m"]
        med = by_year.median().reindex(years)
        if med.notna().sum() < MIN_YEARS:
            continue
        ok = med.notna().to_numpy()
        slope, inter = np.polyfit(years[ok].astype(float), med.to_numpy()[ok], 1)
        cols[r.transect_id] = med.to_numpy() - (slope * years + inter)
        counts[r.transect_id] = by_year.size().reindex(years).to_numpy()
        domains[r.transect_id] = int(r.domain_number)

    M = pd.DataFrame(cols, index=years)
    M.index.name = "year"
    N = pd.DataFrame(counts, index=years)
    return M, N, pd.Series(domains, name="domain_number")


def index_table(M, N, domains):
    """One row per year: the island index, the split, and the sampling."""
    dom = domains.reindex(M.columns).to_numpy()
    touched = np.isin(dom, NOURISHED_DOMAINS)
    out = pd.DataFrame({
        "island_mean_m": M.mean(axis=1).round(2),
        "untouched_mean_m": M.loc[:, ~touched].mean(axis=1).round(2),
        "nourished_mean_m": M.loc[:, touched].mean(axis=1).round(2),
        "spread_m": M.std(axis=1).round(2),
        "n_transects": M.notna().sum(axis=1),
        "obs_per_transect": (N.sum(axis=1) / N.shape[1]).round(1),
    })
    out.index.name = "year"
    return out.reset_index()


def domain_table(M, domains, index):
    """One row per domain: its step, and how well it tracks the island."""
    dom = domains.reindex(M.columns)
    step = (M.loc[POST[0]:POST[1]].mean() - M.loc[PRE[0]:PRE[1]].mean())
    corr = M.apply(lambda c: c.corr(index))
    d = pd.DataFrame({"domain_number": dom, "step_m": step, "r_with_index": corr})
    g = d.groupby("domain_number")
    out = pd.DataFrame({
        "n_transects": g.size(),
        "step_m": g["step_m"].mean().round(2),
        "step_sd_within_m": g["step_m"].std().round(2),
        "r_with_index": g["r_with_index"].median().round(3),
    }).reset_index()
    out["nourished"] = out["domain_number"].isin(NOURISHED_DOMAINS)
    return out


# ============================================================
# THE FIGURE
# ============================================================

def draw(index, by_domain, out_dir):
    """Three panels: what the signal is, where it is, and how well sampled."""
    fs.apply_style()
    fig, axes = plt.subplots(3, 1, figsize=fs.figsize("double", height=7.8),
                             layout="constrained")
    yrs = index["year"].to_numpy()

    # (a) the index itself
    ax = axes[0]
    v = index["island_mean_m"].to_numpy()
    ax.bar(yrs, v, width=0.78, lw=0,
           color=[fs.C["LATE"] if x >= 0 else fs.C["EARLY"] for x in v], zorder=3)
    ax.axhline(0.0, color=fs.C["INK"], lw=0.8, zorder=4)
    # One label per YEAR, not per project: Buxton and Avon are both 2022 and
    # two rotated labels on the same bar overprint into a smear.
    by_year = {}
    for name, (fy, _ds) in NOURISHED.items():
        by_year.setdefault(fy, []).append(name.rsplit(" ", 1)[0])
    for fy, names in by_year.items():
        ax.annotate(" + ".join(names) + " {0}".format(fy),
                    xy=(fy, ax.get_ylim()[0]), xytext=(0, 3),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=6, color=fs.C["INK_MUTED"], rotation=90)
    ax.set_ylabel("detrended position (m)")
    fs._title(ax, 0, "Detrended position, island mean: every transect detrended against "
                     "its own {0}–{1} fit".format(FIRST_YEAR, LAST_YEAR))
    ax.grid(True, axis="y", alpha=0.6)

    # (b) where the step lands
    ax = axes[1]
    x = by_domain["domain_number"].to_numpy()
    s = by_domain["step_m"].to_numpy()
    ax.bar(x, s, width=0.85, lw=0,
           color=[fs.C["ADDED"] if t else fs.C["LATE"]
                  for t in by_domain["nourished"]], zorder=3)
    ax.axhline(0.0, color=fs.C["INK"], lw=0.8, zorder=4)
    ax.axhline(float(np.median(s)), color=fs.C["ACCENT"], lw=1.0,
               ls=(0, (4, 2)), zorder=5)
    ax.set_ylabel("step (m)")
    ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    fs._title(ax, 1, "The {0}–{1} step, by domain (mean {2}–{3} minus "
                     "mean {4}–{5})".format(POST[0], POST[1], POST[0],
                                                 POST[1], PRE[0], PRE[1]))
    fs.town_bands(ax, label=True)
    ax.grid(True, axis="y", alpha=0.6)
    ax.set_xlim(x.min() - 0.5, x.max() + 0.5)

    # (c) the record behind it
    ax = axes[2]
    ax.bar(yrs, index["obs_per_transect"], width=0.78, lw=0,
           color=fs.C["BASE"], zorder=3)
    ax.set_ylabel("observations per\ntransect per year")
    ax.set_xlabel("year")
    fs._title(ax, 2, "How densely the record is sampled")
    ax.grid(True, axis="y", alpha=0.6)

    for ax in (axes[0], axes[2]):
        ax.set_xlim(yrs.min() - 0.7, yrs.max() + 0.7)

    fig.legend(handles=[
        Line2D([], [], color=fs.C["LATE"], lw=6, label="seaward / not nourished"),
        Line2D([], [], color=fs.C["EARLY"], lw=6, label="landward"),
        Line2D([], [], color=fs.C["ADDED"], lw=6, label="a nourished domain"),
        Line2D([], [], color=fs.C["ACCENT"], lw=1.0, ls=(0, (4, 2)),
               label="island median step"),
    ], loc="outside lower center", ncol=4)

    paths = fs.save(fig, Path(out_dir) / "detrended_position", close=True)
    step = index.set_index("year")
    jump = step.loc[POST[0], "island_mean_m"] - step.loc[PRE[1], "island_mean_m"]
    med = float(np.median(s))
    fs.record_caption(paths[0],
        "One signal the whole island shares. (a) every CoastSat transect "
        "detrended against its own {0}–{1} fit, reduced to an annual "
        "median and averaged over {2} transects, so anything local cancels and "
        "what is left is common to the island: roughly flat to 2004, a landward "
        "sag through 2005–2020, then {3:+.1f} m in the single year to "
        "{4}, held to {1}. Nourishment years are labelled. (b) that step by "
        "domain, the mean over {5}–{6} minus the mean over {7}–{8}; "
        "amber domains were nourished, the dashed line is the island median of "
        "{9:+.1f} m. It is smooth alongshore and reach-dependent, not a "
        "constant offset. (c) the sampling density behind it — the record "
        "thickens sharply in the same year the step appears, which is why "
        "`coastsat_position_attribution.py` tests whether the step is the shoreline or "
        "the satellites. The departure carries a bias onto any fitted rate: "
        "-0.37 m/yr on 1996–2010 and +0.88 m/yr on 2010–2024, against "
        "zero on {0}–{1} by construction."
        .format(FIRST_YEAR, LAST_YEAR, int(by_domain["n_transects"].sum()),
                jump, POST[0], POST[0], POST[1], PRE[0], PRE[1], med))
    return paths[0]


README = """# 1-observations/detrended_position — the signal the whole island shares

Every CoastSat transect detrended against its **own** {first}–{last} fit,
reduced to one annual median, and averaged over all {n} of them. Anything local
is incoherent between transects and cancels; what survives is common to the
island.

## Why it exists

`3-rates/coastsat/window_convergence/` found that no window shorter than about
25 years recovers the long-term rate, and that the answer **barely varies
between transects** — a fast, clean transect needs as long as a slow, noisy
one (correlation between a transect's noise-to-trend ratio and its convergence
time: 0.01). A per-transect cause cannot produce a per-transect-invariant
answer, so the cause had to be shared. This is the search for it.

## What it found

One coherent excursion: roughly flat {first}–2004, a landward sag through
2005–2020 bottoming at −7.4 m, then **{jump:+.1f} m in the single year
to {post0}**, held through {last}. It is only 17% of the mean transect variance
— but it is the *coherent* part, and coherence is what moves an OLS slope
systematically while incoherent scatter averages out inside the window.

The bias it alone puts on a fitted rate:

| window | bias from the anomaly | against a median rate of ~1 m/yr |
|---|---|---|
| 1996–2010 | −0.367 m/yr | 37% of it, erosional |
| 2010–2024 | +0.875 m/yr | 87% of it, the other way |
| {first}–{last} | 0.000 m/yr | zero by construction |

That is the answer to *why the windows disagree*: not noise, but a shared
multi-decadal excursion an OLS slope cannot separate from the trend until the
window spans the whole of it.

```
annual_medians_detrended.csv   year x transect, m. The matrix every other
                               script here reads, so the detrending happens
                               once and cannot drift
detrended_position_by_year.csv       a row per year: the index, the nourished /
                               untouched split, the sampling density
detrended_position_by_domain.csv   a row per domain: its share of the step, and
                               how well it tracks the index
detrended_position.png             the index, the alongshore step, the record
attribution_*                  written by coastsat_position_attribution.py — what
                               the step IS
```

## What this is not

Not a rate product and not a model input. Nothing is graded against it. It
lives under `1-observations/` because its subject is the observed record
itself, not a rate fitted from it.

Producers: `scripts/input_prep/5-scr/1-observations/detrended_position/`
(`coastsat_detrended_position.py` builds it, `coastsat_position_attribution.py` tests it). Built
{today} by interview (Hannah).
"""


def main():
    out_dir = obs.DETRENDED_POSITION
    out_dir.mkdir(parents=True, exist_ok=True)

    print("building the detrended annual-median matrix ...")
    M, N, domains = build_matrix()
    print("  {0} transects x {1} years".format(M.shape[1], M.shape[0]))

    M.round(3).to_csv(out_dir / "annual_medians_detrended.csv")
    index = index_table(M, N, domains)
    index.to_csv(out_dir / "detrended_position_by_year.csv", index=False)
    by_domain = domain_table(M, domains, M.mean(axis=1))
    by_domain.to_csv(out_dir / "detrended_position_by_domain.csv", index=False)
    draw(index, by_domain, out_dir)

    # .loc throughout: a bare [2021:2024] on an integer index is POSITIONAL
    # in pandas and silently returns nothing, which printed the step as NaN.
    ix = index.set_index("year")["island_mean_m"]
    jump = ix.loc[POST[0]] - ix.loc[PRE[1]]
    (out_dir / "README.md").write_text(README.format(
        first=FIRST_YEAR, last=LAST_YEAR, n=M.shape[1], jump=jump,
        post0=POST[0], today=pd.Timestamp.today().strftime("%Y-%m-%d")),
        encoding="utf-8")

    print("\nthe index, m:")
    print(index[["year", "island_mean_m", "obs_per_transect"]].to_string(index=False))
    print("\nstep {0}-{1} minus {2}-{3}: {4:+.1f} m island-wide"
          .format(POST[0], POST[1], PRE[0], PRE[1],
                  ix.loc[POST[0]:POST[1]].mean()
                  - ix.loc[PRE[0]:PRE[1]].mean()))
    print("wrote {0}".format(out_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
