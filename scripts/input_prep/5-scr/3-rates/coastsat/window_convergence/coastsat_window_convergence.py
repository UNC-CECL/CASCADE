"""
coastsat_window_convergence.py
==============================================================================
WHICH WINDOWS RECOVER THE LONG-TERM RATE, AND WHICH ARE TOO SHORT?

The model is graded on the LRR over 1996-2010. Is that enough record for the
rate to represent the shoreline, or is it still an artefact of where the record
happens to stop? And if it is not, what window would be?

A NOTE ON COUNTING YEARS. A window here is counted in CALENDAR YEARS OF
RECORD, both ends included, because that is what the fit consumes: 1996-2010
is fifteen years of satellite positions. Rule 2 of ORGANIZATION.md counts the
same window as fourteen SIMULATED years, because the model spends 1996..2009
and the end year is a boundary. Both are right about different things; the
tables and figures here are always the fit's count.

TWO SWEEPS, ONE FROM EACH END (Hannah, 2026-09-23). Both fit the same OLS the
target uses (`coastsat_lrr.compute_lrr`) on a family of NESTED windows, and
both converge on the same reference, the 1996-2024 rate -- but they approach
it from opposite sides, and each answers a different question:

    forward_from_1996/   the START is pinned at 1996 and the END walks out:
                         1996-2000, 1996-2001 ... 1996-2024.
                         HOW MUCH RECORD DO YOU NEED from the start of the
                         chain before the rate settles?

    backward_from_2024/  the END is pinned at 2024 and the START walks back:
                         2020-2024, 2019-2024 ... 1996-2024.
                         HOW LATE CAN A WINDOW BEGIN and still recover the
                         long-term rate?

The pair brackets the answer. Forward gives a window 1996-YYYY; backward gives
a window YYYY-2024. In both the reference is the LONGEST window of that sweep,
which is 1996-2024 either way, fitted in the same loop as every other window
so it cannot drift from a stored product.

The moving year 2010 is marked in both, and in both it is a real model window:
forward that is 1996-2010, the graded window, and backward it is 2010-2024,
the second leg of the canonical chain.

TWO SCALES, in two folders under each direction:

    a-eight_sites/  one transect at the middle of each of eight evenly spaced
                    domains. The readable case: every position, every fit, one
                    panel per site. Eight transects cannot speak for an
                    island, but they show WHAT is happening.
    b-every_transect/  every CoastSat transect on the island, ~906 of them, the
                    same sweep. Answers whether the eight were representative:
                    the convergence year as an alongshore profile, and the
                    spread within each domain.

WHAT THIS CAN AND CANNOT SAY
    The windows are NESTED, so a curve converges on the reference BECAUSE the
    reference is its endpoint. "Does 1996-2010 match 1996-2024" therefore has
    no yes/no answer, and it does not need one -- what the sweep gives is the
    shape of the approach and the window at which the rate stops leaving a
    tolerance. Read it as "this window recovers the long-term rate", never as
    a match test.

    It also cannot separate a rate that was WRONG from one that was merely
    EARLY. A site whose shoreline genuinely changed behaviour in 2010 and a
    site that simply needed more observations draw the same curve. Disjoint
    windows are what separate those -- 3-rates/coastsat/5yr_bins/ is the
    product that asks when the change happened. Running BOTH directions is
    the cheapest partial answer: a real change of behaviour shows as a
    forward sweep and a backward sweep that disagree about where the good
    window is.

THREE TOLERANCES, reported side by side rather than chosen, so the sensitivity
to the choice is visible in the table:

    ci    within the 1996-2024 fit's own 95% confidence half-width. Scales
          with how well constrained each site is; nothing arbitrary to defend.
    abs   within +/-0.25 m/yr. The same band everywhere, so sites compare;
          ~7 m of shoreline over 28 yr, small against rates of 2-3 m/yr.
    rel   within +/-20% of the reference. Sensible where the rate is fast,
          unreachable where it is near zero -- which the table will show.

And two readings of each, because they differ exactly where the answer is
interesting:

    first entry    the shortest window inside the band. May be a crossing the
                   curve then leaves again.
    stable entry   the shortest window after which every LONGER window is
                   also inside. This is the convergence window; the gap
                   between the two is a site that wandered back out.

Inputs
------
    transect_domain_lookup.csv      2-transect-frame/, via hat_observed_rates
    CoastSat time-series CSVs       1-observations/coastsat_timeseries/
    coastsat_lrr.compute_lrr        5-scr/lib/, via scr_paths

Outputs  (hat_observed_rates.window_convergence_dir(direction, anchor))
-----------------------------------------------------------------------
    2-settling_window/<direction>_from_<anchor>/
        a-eight_sites/
            shoreline_position_window_fits_*.png   the record: positions, the
                                                   annual median, six fits
            window_convergence_*.png               the sweep as an ERROR
            window_convergence_transects.csv       a row per site per window
            convergence_summary.csv                a row per site
        b-every_transect/
            convergence_alongshore_*.png           the island profile
            domain_convergence_summary.csv         a row per GIS domain
            convergence_summary_all_transects.csv  a row per transect
            window_convergence_transects_all.csv   the full sweep
        c-domain_means/
            the transect sweep grouped to the 90 domains (no refit)
    README.md beside each, supporting/ for PDFs and CAPTIONS.md
    A record other than 1996-2024 files under experiments/record_cut_<end>/.

Usage
-----
    python .../coastsat_window_convergence.py                 both, both scales
    python .../coastsat_window_convergence.py --direction forward --scale sites
    (--abs-tol, --rel-tol, --domains to vary it)
==============================================================================
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Rule 5: find the root by searching upward.
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
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# ============================================================
# CONFIG
# ============================================================

# THE RECORD THE SWEEP MAY SEE, and the two pins. The longest window of either
# sweep is this one, so both directions converge on the same number, and every
# error in the product is measured against it.
#
# TRUNCATING IT IS AN EXPERIMENT, NOT A SETTING. The island-wide detrended
# position anomaly steps +17 m between 2020 and 2021 and holds there, in the
# same year the satellite record goes from 12.8 to 30.4 observations per
# transect per year. Those four years sit at the highest-leverage end of the
# 1996-2024 fit, so if the step is an artefact of the sensor mix rather than
# the shoreline, the reference itself is biased -- and the reference is the
# model's grading target. `--ref-end 2020` refits everything on the record
# before the step, and the two record spans file side by side under
# `experiments/record_cut_<end>/`, apart from the main result in
# `2-settling_window/`, so they can be differenced rather than confused.
REF_START, REF_END = 1996, 2024

# The shortest window either sweep fits. Below five years an OLS through ~90
# positions is describing a storm cycle, not a trend.
MIN_WINDOW_YEARS = 5

# EIGHT EVENLY SPACED DOMAINS (Hannah, 2026-09-23). The 90 domains run south
# (1, Cape Point) to north (90, Pea Island); these are every 12th or 13th,
# both ends included. Not chosen for behaviour -- an even spread has no
# selection argument to defend, and the alongshore gradient is covered.
SITE_DOMAINS = [1, 13, 26, 38, 51, 64, 77, 90]

# ONE LITERAL TRANSECT PER SITE, not the domain mean (Hannah, 2026-09-23):
# the noisiest case, and so an upper bound on how long convergence takes. The
# transect is the MIDDLE one of its domain by alongshore order -- a centroid
# proxy that needs no geometry, and each domain holds 9-13 transects.
MIN_OBS = 10                # a window with fewer positions gets no fit

ABS_TOL_M_YR = 0.25         # the "abs" band
REL_TOL = 0.20              # the "rel" band, as a fraction of the reference

# SEVEN TOLERANCES, AND NONE OF THEM IS THE ANSWER (Hannah, 2026-09-23, on
# being shown how strict the first three were: "would it be possible to make it
# more forgiving so it is not an exact match but is close?").
#
# The first three ask a window to land inside a band set by the REFERENCE, and
# the reference is the most precise fit in the record -- a median 95% half-width
# of 0.151 m/yr over 906 transects. A 15-year window is itself uncertain to
# about 0.44 m/yr, so those three judge an imprecise estimate as though it were
# as precise as the thing it is being compared to. That asymmetry, not the
# shoreline, is why almost nothing passed them.
#
# `overlap` is the principled loosening and the one to reach for first: a window
# passes when its OWN 95% interval reaches the reference's, i.e. when the two
# fits are not statistically distinguishable. It is a FUNNEL, not a band -- wide
# where the record is short and narrowing onto the reference -- which is exactly
# the shape the question deserves and the reason it cannot be expressed as a
# number. The rest are arbitrary widths, stated openly as such so the answer's
# sensitivity to the choice is on the page rather than in a decision no one
# wrote down.
#
# ONE THAT WAS TRIED AND REJECTED: the nested-sample standard error,
# sqrt(u_window^2 - u_ref^2), the Hausman variance of the difference between an
# efficient estimator and a less efficient one. It is correct for nested windows
# and more forgiving where the record is short -- but its band COLLAPSES TO ZERO
# as the window approaches the full record, and "stable entry" is decided at the
# long end, so it comes out STRICTER than everything here (1996-2024 in both
# directions). Kept in this comment because a plausible-sounding statistic that
# fails for a structural reason is worth one paragraph.
#
# Each entry: (tag, label for figures and prose, test(diff, unc_window,
# unc_ref, ref_lrr) -> bool array).
CRITERIA = [
    ("ci", "the reference fit's 95% band",
     lambda d, uw, ur, ref: np.abs(d) <= ur),
    ("abs", "±0.25 m/yr",
     lambda d, uw, ur, ref: np.abs(d) <= ABS_TOL_M_YR),
    ("rel", "±20%",
     lambda d, uw, ur, ref: np.abs(d) <= REL_TOL * np.abs(ref)),
    ("ci3x", "3× the reference band",
     lambda d, uw, ur, ref: np.abs(d) <= 3.0 * ur),
    ("abs50", "±0.50 m/yr",
     lambda d, uw, ur, ref: np.abs(d) <= 0.50),
    ("abs100", "±1.00 m/yr",
     lambda d, uw, ur, ref: np.abs(d) <= 1.00),
    ("overlap", "CI overlap",
     lambda d, uw, ur, ref: np.abs(d) <= uw + ur),
]
CRITERION_LABEL = dict((tag, label) for tag, label, _ in CRITERIA)

# The five Hannah asked to see side by side, in the order they loosen.
HEADLINE_TAGS = ["ci", "ci3x", "abs50", "overlap", "abs100"]

# THE ONE THE PRODUCT ANSWERS WITH (Hannah, 2026-09-23). Every other tolerance
# is still scored and tabled; this is the one the figures draw, the panel
# titles name and the READMEs lead with.
#
# CI overlap, because it is the only criterion here that is not a number
# somebody chose. It asks whether a window's rate is STATISTICALLY
# DISTINGUISHABLE from the long-term rate -- the two 95% intervals meet -- so a
# five-year fit is judged against what a five-year fit can actually resolve
# (+/-1.8 m/yr) and a twenty-five-year fit against what that can (+/-0.2). The
# strict bands ask both to land inside +/-0.15, which is the reference's own
# precision and no one else's.
#
# It is worth knowing it barely changes the answer: the island median moves
# from 27 to 26 years forward and 25 to 22 backward. That is the finding, not a
# disappointment -- the windows disagree because the shoreline changed, not
# because the fits were noisy, so no defensible loosening rescues a 15-year
# window. The strict bands stay in the tables so that claim can be checked.
HEADLINE_TAG = "overlap"
HEADLINE_LABEL = "CI overlap"

# THE MARKED YEAR, and it is the same number in both directions because in
# both it names a real model window: forward it is the end of 1996-2010, the
# window the model is graded on; backward it is the start of 2010-2024, the
# second leg of the canonical chain. See [[cascade-canonical-periods]].
MARKED_YEAR = 2010

# THE FIGURES DO NOT SINGLE OUT THE MODEL WINDOW for now (Hannah, 2026-09-29):
# she is choosing which window to use, so nothing is drawn in amber -- no
# 15-year line on the years-needed figure, no amber fit or dashed continuation
# on the eight-site figure; the marked window is drawn like any other. The
# tables and README text still score it. Set True to bring the highlight back.
HIGHLIGHT_MODEL_WINDOW = False

# SIX WINDOWS ARE DRAWN, NOT TWENTY-FIVE (Hannah, 2026-09-23: "too many lines
# to distinguish anything"). The sweep still FITS every window -- the tables
# carry all of them -- but a panel that draws them all is a smear in which the
# two that matter, the marked window and the reference, are lost among
# twenty-three that differ from their neighbour by one year of record.
DRAWN_YEARS = [2000, 2005, 2010, 2015, 2020, 2024]


def windows_for(direction):
    """The nested family, SHORTEST FIRST, so the reference is always last.

    Sorting by length rather than by year is what lets one convergence walk
    serve both directions: "stable entry" always means "the shortest window
    after which every longer one is also inside the band".

    Returns [(start, end, moving_year), ...]. The moving year is the end that
    is not pinned -- the end year going forward, the start year going back.
    """
    if direction == "forward":
        years = range(REF_START + MIN_WINDOW_YEARS - 1, REF_END + 1)
        return [(REF_START, y, y) for y in years]
    if direction == "backward":
        years = range(REF_END - MIN_WINDOW_YEARS + 1, REF_START - 1, -1)
        return [(y, REF_END, y) for y in years]
    raise ValueError("direction is 'forward' or 'backward', not {0!r}".format(direction))


def pinned_year(direction):
    return REF_START if direction == "forward" else REF_END


def window_label(direction, moving_year):
    """`1996-2010` going forward, `2010-2024` going back.

    The year is rounded because callers pass medians as well as single
    windows, and a median over 906 transects arrives as 2022.0 -- which then
    printed "1996-2022.0" into a README (found 2026-09-23).
    """
    year = int(round(float(moving_year)))
    if direction == "forward":
        return "{0}–{1}".format(REF_START, year)
    return "{0}–{1}".format(year, REF_END)


def moving_axis_label(direction):
    if direction == "forward":
        return "window END year  (every window starts {0})".format(REF_START)
    return "window START year  (every window ends {0})".format(REF_END)


# ============================================================
# THE SWEEP
# ============================================================

def load_lookup():
    """The transect-to-domain table, from the 2-transect-frame resolver."""
    return pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv")


def pick_transects(domains):
    """The middle transect of each domain, by alongshore (sorted id) order."""
    lookup = load_lookup()
    picks = []
    for d in domains:
        ids = sorted(lookup.loc[lookup["domain_number"] == float(d), "transect_id"])
        if not ids:
            raise ValueError(
                "domain {0} has no transects in the lookup; domains on file "
                "are {1:.0f}-{2:.0f}".format(d, lookup["domain_number"].min(),
                                             lookup["domain_number"].max()))
        picks.append((int(d), ids[len(ids) // 2]))
    return picks


def all_transects():
    """Every transect in the lookup that has a domain, alongshore order."""
    lookup = load_lookup().dropna(subset=["domain_number"])
    lookup = lookup.sort_values(["domain_number", "transect_id"])
    return [(int(r.domain_number), r.transect_id)
            for r in lookup.itertuples(index=False)]


def timeseries_path(transect_id):
    """`usa_NC_0034_0054` -> its CSV under coastsat_timeseries/."""
    site = transect_id.rsplit("_", 1)[0]
    path = obs.COASTSAT_TIMESERIES / "{0}_timeseries".format(site) / "{0}.csv".format(transect_id)
    if not path.is_file():
        raise FileNotFoundError(
            "no chainage series for {0}: expected {1}".format(transect_id, path))
    return path


def sweep_one(transect_id, domain, windows):
    """Every window of one family for one transect, as a list of row dicts.

    The series is loaded and sorted once and each window is a SLICE of it, so
    a 906-transect run is a couple of minutes rather than an afternoon.
    `compute_lrr` still does every fit, so the estimator is the target's.

    The reference row (the longest window) is fitted in the same loop as the
    rest, so it cannot disagree with them.
    """
    df = cl.load_timeseries(str(timeseries_path(transect_id)))
    # Dropped to naive UTC purely so searchsorted has a datetime64 array to
    # bisect: a tz-aware column comes back as objects and compares against
    # nothing. The column itself, which compute_lrr reads, stays tz-aware.
    dates = df["date"].dt.tz_localize(None).to_numpy()
    rows = []
    for start, end, moving in windows:
        lo = np.searchsorted(dates, np.datetime64("{0}-01-01T00:00:00".format(start)),
                             side="left")
        hi = np.searchsorted(dates, np.datetime64("{0}-01-01T00:00:00".format(end + 1)),
                             side="left")
        clipped = df.iloc[lo:hi]
        if len(clipped) < MIN_OBS:
            res = cl._empty_lrr(len(clipped))
            intercept = np.nan
        else:
            res = cl.compute_lrr(clipped)
            # The intercept the slope implies, so the straight line can be
            # DRAWN without a second fit disagreeing with the stored one. x is
            # years since the window's FIRST observation, as compute_lrr sets it.
            t0 = clipped["date"].min()
            x = (clipped["date"] - t0).dt.total_seconds().to_numpy() / (86400.0 * 365.25)
            intercept = float(clipped["chainage_m"].to_numpy().mean()
                              - res["lrr_m_yr"] * x.mean())
        rows.append({
            "domain_number": domain,
            # The UNIT the sweep is about. A transect names itself; a domain
            # mean names its domain. Everything downstream -- scoring, the
            # convergence walk, the figures -- keys on this, so the same code
            # serves a transect and a 10-transect mean without a branch.
            "unit_id": transect_id,
            "transect_id": transect_id,
            "start_year": start,
            "end_year": end,
            "moving_year": moving,
            "window": "{0}_{1}".format(start, end),
            "n_years": end - start + 1,
            "lrr_m_yr": res["lrr_m_yr"],
            "intercept_m": round(intercept, 3) if intercept == intercept else np.nan,
            "unc_m_yr": res["unc_m_yr"],
            "r_squared": res["r_squared"],
            "p_value": res["p_value"],
            "n_obs": res["n_obs"],
            "first_obs": res["start_date"],
            "last_obs": res["end_date"],
        })
    return rows


def score(rows):
    """Add the reference, the difference and the three in-band flags in place.

    The reference is the LAST row -- the longest window of the family, which
    is 1996-2024 in either direction.
    """
    ref = rows[-1]
    ref_lrr, ref_unc = ref["lrr_m_yr"], ref["unc_m_yr"]
    for r in rows:
        diff = r["lrr_m_yr"] - ref_lrr
        r["ref_lrr_m_yr"] = ref_lrr
        r["ref_unc_m_yr"] = ref_unc
        r["diff_m_yr"] = round(diff, 4)
        r["abs_diff_m_yr"] = round(abs(diff), 4)
        ok = not np.isnan(diff)
        for tag, _label, test in CRITERIA:
            passes = test(np.array([diff]), np.array([r["unc_m_yr"]]),
                          np.array([ref_unc]), np.array([ref_lrr]))
            r["in_" + tag] = bool(ok and bool(passes[0]))
    return rows


def convergence_entry(sub, flag):
    """(first entry, stable entry) for one transect under one flag.

    `sub` is ordered SHORTEST WINDOW FIRST. First entry is the shortest window
    inside the band. Stable entry is the shortest window after which every
    longer one is also inside -- the convergence window. Both are reported as
    the MOVING year, so forward gives an end year and backward a start year.
    The longest window is always inside (it IS the reference), so both exist.
    """
    moving = sub["moving_year"].to_numpy()
    inside = sub[flag].to_numpy(dtype=bool)
    first = int(moving[inside][0]) if inside.any() else None
    stable = int(moving[-1])
    for i in range(len(moving) - 1, -1, -1):
        if not inside[i]:
            break
        stable = int(moving[i])
    return first, stable


def summarise(sweep, direction):
    """One row per transect: reference, the marked window, the entries."""
    out = []
    for (domain, tid), sub in sweep.groupby(["domain_number", "unit_id"], sort=False):
        sub = sub.sort_values("n_years")
        ref = sub.iloc[-1]
        marked = sub[sub["moving_year"] == MARKED_YEAR]
        marked = marked.iloc[0] if len(marked) else None
        row = {
            "domain_number": domain,
            "unit_id": tid,
            "transect_id": sub["transect_id"].iloc[0],
            "direction": direction,
            "ref_window": ref["window"],
            "ref_lrr_m_yr": ref["lrr_m_yr"],
            "ref_unc_m_yr": ref["unc_m_yr"],
            "ref_n_obs": ref["n_obs"],
        }
        if marked is not None:
            row.update({
                "marked_window": marked["window"],
                "lrr_marked_m_yr": marked["lrr_m_yr"],
                "unc_marked_m_yr": marked["unc_m_yr"],
                "diff_marked_m_yr": marked["diff_m_yr"],
            })
            for tag, _label, _test in CRITERIA:
                row["marked_in_" + tag] = marked["in_" + tag]
        for tag, _label, _test in CRITERIA:
            first, stable = convergence_entry(sub, "in_" + tag)
            row["first_entry_" + tag] = first
            row["stable_entry_" + tag] = stable
            row["years_needed_" + tag] = (
                stable - REF_START + 1 if direction == "forward"
                else REF_END - stable + 1)
        out.append(row)
    return pd.DataFrame(out)


def summarise_domains(summary):
    """One row per GIS domain: the median and spread across its transects.

    The median, not the mean: a single transect that never settles would drag
    a mean to the end of the record and say the whole domain did.
    """
    g = summary.groupby("domain_number")
    out = pd.DataFrame({
        "n_transects": g.size(),
        "ref_lrr_median_m_yr": g["ref_lrr_m_yr"].median().round(4),
        "diff_marked_median_m_yr": g["diff_marked_m_yr"].median().round(4),
        "diff_marked_q25_m_yr": g["diff_marked_m_yr"].quantile(0.25).round(4),
        "diff_marked_q75_m_yr": g["diff_marked_m_yr"].quantile(0.75).round(4),
    })
    for tag, _label, _test in CRITERIA:
        out["marked_in_" + tag + "_frac"] = g["marked_in_" + tag].mean().round(4)
        out["years_needed_" + tag + "_median"] = g["years_needed_" + tag].median()
        out["years_needed_" + tag + "_q25"] = g["years_needed_" + tag].quantile(0.25)
        out["years_needed_" + tag + "_q75"] = g["years_needed_" + tag].quantile(0.75)
        out["stable_entry_" + tag + "_median"] = g["stable_entry_" + tag].median()
    return out.reset_index()


def run_sweep(picks, windows, announce=False):
    """The sweep over a list of (domain, transect_id), as one DataFrame."""
    rows = []
    for i, (domain, tid) in enumerate(picks, start=1):
        site_rows = score(sweep_one(tid, domain, windows))
        if announce:
            ref = site_rows[-1]
            print("  GIS {0:>2}  {1}  {2} = {3:+.3f} +/- {4:.3f} m/yr  ({5} obs)"
                  .format(domain, tid, ref["window"], ref["lrr_m_yr"],
                          ref["unc_m_yr"], ref["n_obs"]))
        elif i % 150 == 0 or i == len(picks):
            print("    {0}/{1} transects".format(i, len(picks)))
        rows.extend(site_rows)
    return pd.DataFrame(rows)


def domain_mean_sweep(sweep):
    """The transect sweep aggregated to the model's own unit: the domain MEAN.

    WHY A MEAN OF SLOPES, NOT A FIT THROUGH POOLED POSITIONS. This has to be
    the same construction as the grading target, and `coastsat_domain_lrr.py`
    builds that by fitting each transect and averaging the slopes. Pooling the
    positions first would be a different number, because each transect's
    chainage sits on its own arbitrary origin.

    No refitting happens here. Every window of every transect is already in
    `sweep`; this is a groupby.

    TWO UNCERTAINTIES, and the honest one is the wider. `unc_m_yr` is the MEAN
    of the transects' 95% half-widths -- the typical fit uncertainty in the
    domain. `unc_if_independent_m_yr` is what the half-width of the mean would
    be if the ~10 transects were independent samples, which they are not:
    they are 10-metre-spaced views of the same shoreline and move together.
    Propagating as if they were would divide the band by about sqrt(10) and
    manufacture a much later convergence window out of an assumption. The
    scoring uses the mean, which is conservative -- it makes convergence look
    EARLIER, not later -- and the independent figure is carried as a column so
    the size of that choice is visible rather than buried.
    """
    keys = ["domain_number", "start_year", "end_year", "moving_year",
            "n_years", "window"]
    g = sweep.groupby(keys, sort=True)
    out = g.agg(
        lrr_m_yr=("lrr_m_yr", "mean"),
        lrr_median_m_yr=("lrr_m_yr", "median"),
        lrr_std_m_yr=("lrr_m_yr", "std"),
        unc_m_yr=("unc_m_yr", "mean"),
        r_squared=("r_squared", "mean"),
        n_obs=("n_obs", "sum"),
        n_transects=("lrr_m_yr", "count"),
        first_obs=("first_obs", "min"),
        last_obs=("last_obs", "max"),
    ).reset_index()
    sq = g["unc_m_yr"].apply(lambda u: float(np.sqrt(np.nansum(u.to_numpy() ** 2))))
    out["unc_if_independent_m_yr"] = (sq.to_numpy() / out["n_transects"].to_numpy()).round(4)
    for col in ("lrr_m_yr", "lrr_median_m_yr", "lrr_std_m_yr", "unc_m_yr", "r_squared"):
        out[col] = out[col].round(4)
    out["unit_id"] = ["GIS {0:.0f} mean".format(d) for d in out["domain_number"]]
    out["transect_id"] = out["unit_id"]
    return out.sort_values(["domain_number", "n_years"]).reset_index(drop=True)


def score_units(frame):
    """Score every unit of a frame, unit by unit, shortest window first."""
    rows = []
    for _, sub in frame.groupby("unit_id", sort=False):
        rows.extend(score(sub.sort_values("n_years").to_dict("records")))
    return pd.DataFrame(rows)


# ============================================================
# THE FIGURES -- the eight sites
# ============================================================

def _annual_median(series):
    """The year-by-year median position: the trajectory under the scatter.

    A CoastSat transect gets ~18 positions a year and their spread is tens of
    metres, so the raw cloud hides the very shape the fits are arguing about.
    The median is per CALENDAR year, matching the window convention.
    """
    by_year = series.groupby(series["date"].dt.year)["chainage_m"].median()
    return by_year.index.to_numpy(), by_year.to_numpy()


def draw_fits(sweep, summary, out_dir, direction):
    """The record itself: position against time, with six windows' fits.

    The companion to the error figure, and the thing it is derived FROM. y is
    shoreline position in metres, x is time across the whole record; each
    straight line is one nested window's OLS drawn over the span it was fitted
    on, so the spread between them IS the disagreement the sweep measures.

    Going forward the marked window is also continued as a DASHED line to the
    end of the record: not a claim about the future, but the plainest way to
    see what the marked window's record would have had you believe about 2024.
    """
    fs.apply_style()
    sites = list(summary.itertuples(index=False))
    nrow, ncol = 4, 2
    fig, axes = plt.subplots(nrow, ncol, figsize=fs.figsize("double", height=8.6),
                             sharex=True, layout="constrained")
    flat = axes.ravel()
    ramp = LinearSegmentedColormap.from_list("windows", list(fs.SMOOTH_RAMP))
    moving_all = sorted(sweep["moving_year"].unique())
    # The reference is always drawn. Its moving year is REF_END going forward
    # but REF_START going back, which DRAWN_YEARS does not hold, so until
    # 2026-09-29 the backward figure listed the reference and never drew it.
    ref_moving = REF_END if direction == "forward" else REF_START
    drawn = sorted(set(y for y in DRAWN_YEARS if y in moving_all) | {ref_moving})
    x0 = pd.Timestamp("{0}-01-01".format(REF_START), tz="UTC")
    x1 = pd.Timestamp("{0}-12-31".format(REF_END), tz="UTC")

    for i, site in enumerate(sites):
        ax = flat[i]
        sub = sweep[sweep["transect_id"] == site.transect_id]
        series = cl.filter_dates(cl.load_timeseries(str(timeseries_path(site.transect_id))),
                                 "{0}-01-01".format(REF_START),
                                 "{0}-12-31".format(REF_END))

        ax.scatter(series["date"], series["chainage_m"], s=2.2,
                   color=fs.C["BASE"], alpha=0.16, lw=0, zorder=1)

        # The trajectory the fits are arguing about.
        yrs, med = _annual_median(series)
        ax.plot([pd.Timestamp("{0}-07-01".format(y), tz="UTC") for y in yrs], med,
                color=fs.C["INK_MUTED"], lw=1.0, marker="o", ms=2.4,
                mfc="white", mew=0.6, zorder=4)

        for row in sub.itertuples(index=False):
            if row.moving_year not in drawn or row.intercept_m != row.intercept_m:
                continue
            t0 = pd.Timestamp(row.first_obs, tz="UTC")
            t1 = pd.Timestamp(row.last_obs, tz="UTC")
            span = (t1 - t0).total_seconds() / (86400.0 * 365.25)
            y0, y1 = row.intercept_m, row.intercept_m + row.lrr_m_yr * span
            if row.moving_year == MARKED_YEAR and HIGHLIGHT_MODEL_WINDOW:
                colour, lw, z = fs.C["ADDED"], 2.0, 8
            elif row.n_years == REF_END - REF_START + 1:
                colour, lw, z = fs.C["ACCENT"], 2.0, 9
            else:
                frac = (row.moving_year - min(drawn)) / max(max(drawn) - min(drawn), 1)
                if direction == "backward":
                    frac = 1.0 - frac
                colour, lw, z = ramp(0.20 + 0.75 * frac), 1.5, 5
            ax.plot([t0, t1], [y0, y1], color=colour, lw=lw, zorder=z,
                    solid_capstyle="round")

            # Forward only: where the marked window's record would put 2024.
            if (row.moving_year == MARKED_YEAR and direction == "forward"
                    and HIGHLIGHT_MODEL_WINDOW):
                far = (x1 - t0).total_seconds() / (86400.0 * 365.25)
                ax.plot([t1, x1], [y1, y0 + row.lrr_m_yr * far],
                        color=fs.C["ADDED"], lw=1.2, ls=(0, (3, 2.2)), zorder=7)

            # The moving year at the line's free end, so the ramp needs no key.
            if direction == "forward":
                ax.annotate(str(row.moving_year), xy=(t1, y1), xytext=(-2, 3),
                            textcoords="offset points", ha="right", va="bottom",
                            fontsize=6.2, color=colour,
                            path_effects=fs._halo(2.0), zorder=12)
            else:
                ax.annotate(str(row.moving_year), xy=(t0, y0), xytext=(2, 3),
                            textcoords="offset points", ha="left", va="bottom",
                            fontsize=6.2, color=colour,
                            path_effects=fs._halo(2.0), zorder=12)

        fs._title(ax, i, "GIS {0} · {1}".format(
            site.domain_number, site.transect_id.replace("usa_NC_", "")))
        ax.grid(True, axis="y", alpha=0.6)
        if i % ncol == 0:
            ax.set_ylabel("shoreline position (m)")
        if i >= len(sites) - ncol:
            ax.set_xlabel("year")
        ax.set_xlim(x0, x1)

    for ax in flat[len(sites):]:
        ax.set_visible(False)

    others = ", ".join(str(y) for y in drawn
                       if y != pinned_year(direction)
                       and (y != MARKED_YEAR or not HIGHLIGHT_MODEL_WINDOW)
                       and y != (REF_END if direction == "forward" else REF_START))
    handles = [
        Line2D([], [], color="none", marker="o", ms=3.5, mfc=fs.C["BASE"],
               mec="none", label="CoastSat position"),
        Line2D([], [], color=fs.C["INK_MUTED"], lw=1.0, marker="o", ms=2.4,
               mfc="white", mew=0.6, label="annual median"),
        Line2D([], [], color=ramp(0.45), lw=1.5,
               label="windows at {0}".format(others)),
        Line2D([], [], color=fs.C["ACCENT"], lw=2.0,
               label="{0}–{1}, the reference".format(REF_START, REF_END)),
    ]
    if HIGHLIGHT_MODEL_WINDOW:
        handles.insert(3, Line2D([], [], color=fs.C["ADDED"], lw=2.0,
                                 label="{0}, the marked window".format(
                                     window_label(direction, MARKED_YEAR))))
    if direction == "forward" and HIGHLIGHT_MODEL_WINDOW:
        handles.insert(4, Line2D([], [], color=fs.C["ADDED"], lw=1.2,
                                 ls=(0, (3, 2.2)),
                                 label="that fit continued to {0}".format(REF_END)))
    fig.legend(handles=handles, loc="outside lower center", ncol=3)

    stem = "shoreline_position_window_fits_{0}_from_{1}".format(
        direction, pinned_year(direction))
    paths = fs.save(fig, Path(out_dir) / stem, close=True)
    worst = summary.loc[summary["diff_marked_m_yr"].abs().idxmax()]
    pinned = "starts" if direction == "forward" else "ends"
    extra = ("The amber dash continues the marked fit to {0} — what "
             "{1} calendar years of record would have had you believe about {0}. "
             .format(REF_END, MARKED_YEAR - REF_START + 1)
             if direction == "forward" else "")
    if not HIGHLIGHT_MODEL_WINDOW:
        fs.record_caption(paths[0],
            "The record the sweep is fitted on. Every CoastSat shoreline "
            "position at one transect in each of eight evenly spaced domains, "
            "{0} to {1}, with the annual median through them and some of the "
            "{2} nested fits drawn as straight lines, each spanning the window "
            "it was fitted on and labelled with its moving year. Every window "
            "{3} {4}. Purple is {0}–{1}, the reference every window converges "
            "on. Where the lines separate, the rate still depends on where the "
            "window is cut. Position is CoastSat chainage, seaward positive, on "
            "an origin that is arbitrary per transect: only the SLOPE compares "
            "between panels, and each panel is autoscaled."
            .format(REF_START, REF_END, len(moving_all), pinned,
                    pinned_year(direction)))
        return paths[0]
    fs.record_caption(paths[0],
        "The record the sweep is fitted on. Every CoastSat shoreline position "
        "at one transect in each of eight evenly spaced domains, {0} to {1}, "
        "with the annual median through them and six of the {2} nested fits "
        "drawn as the straight lines they are. Every window {3} {4}; each line "
        "spans the window it was fitted on and is labelled with its moving "
        "year. Amber is {5}, the marked window; {6}purple is {0}–{1}, the "
        "reference every window converges on. Where the lines separate, the "
        "rate still depends on where the window ends. The intermediate windows "
        "are fitted and tabled but not drawn. Position is CoastSat chainage, "
        "seaward positive, on an origin that is arbitrary per transect: only "
        "the SLOPE compares between panels, and each panel is autoscaled. The "
        "widest disagreement is GIS {7:.0f}, where the marked window gives "
        "{8:+.2f} m/yr against the reference's {9:+.2f}."
        .format(REF_START, REF_END, len(moving_all), pinned,
                pinned_year(direction), window_label(direction, MARKED_YEAR),
                extra, worst["domain_number"], worst["lrr_marked_m_yr"],
                worst["ref_lrr_m_yr"]))
    return paths[0]


def alongshore_x(picks):
    """Each transect's x on the domain axis: its domain, with the domain's
    transects spread evenly across [d - 0.5, d + 0.5] in alongshore order, so
    906 transects and the village bands share one axis."""
    frame = pd.DataFrame(picks, columns=["domain_number", "transect_id"])
    rank = frame.groupby("domain_number").cumcount()
    n = frame.groupby("domain_number")["transect_id"].transform("size")
    frame["x_domain"] = frame["domain_number"] - 0.5 + (rank + 0.5) / n
    return frame


# THREE PLAIN THRESHOLDS, NOT SEVEN TOLERANCES (Hannah, 2026-09-29, on the old
# figures: too many of them, an abstract "rate minus reference" axis, and
# tolerance jargon). The figure asks one thing -- how many years of record does
# each transect need before its rate stays within X m/yr of the 1996-2024 rate
# -- at three values of X, in m/yr and nothing else. All seven tolerances are
# still scored in the tables. Strictest last, so it stacks on top.
YEARS_NEEDED_THRESHOLDS = [
    ("abs100", "±1.0 m/yr", "#6baed6"),
    ("abs50", "±0.5 m/yr", "#9ecae1"),
    ("abs", "±0.25 m/yr", "#deebf7"),
]
YEARS_NEEDED_STEM = "years_needed_alongshore"


def draw_years_needed(ref_start, ref_end):
    """THE figure of the settling sweep: years needed, per transect, alongshore.

    One panel per direction, read from the stored every-transect summaries, so
    it can be redrawn without refitting. The three thresholds are NESTED -- a
    rate that stays within 0.25 m/yr also stays within 0.5 -- so they stack as
    shaded bands rather than crossing as lines. The model's 15-year window is
    the amber line: wherever a band reaches above it, 15 years was not enough
    at that threshold.
    """
    fs.apply_style()
    # x IS THE TRANSECT COUNT, 1-906 south to north (Hannah, 2026-09-29): every
    # transect one unit wide, nothing stretched to fit a domain. The domains
    # are only labels, on a second axis under panel (b).
    order = pd.DataFrame(all_transects(), columns=["domain_number", "transect_id"])
    order["x"] = np.arange(1, len(order) + 1)
    xmap = order.set_index("transect_id")["x"]
    n_tr = len(order)
    first = order.groupby("domain_number")["x"].min()
    last = order.groupby("domain_number")["x"].max()
    centre = (first + last) / 2.0
    try:
        from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS
        spans = {name: (first[lo], last[hi])
                 for name, (lo, hi) in HATTERAS_ANNOTATIONS.town_spans.items()
                 if lo in first.index and hi in last.index}
    except ImportError:
        spans = {}
    fig, axes = plt.subplots(2, 1, figsize=fs.figsize("double", height=6.4),
                             sharex=True, layout="constrained")
    marked_years = MARKED_YEAR - ref_start + 1
    out_dir = None
    notes = []
    for i, direction in enumerate(("forward", "backward")):
        pinned = ref_start if direction == "forward" else ref_end
        root = obs.window_convergence_dir(direction, pinned, ref_start, ref_end)
        out_dir = root.parent
        s = pd.read_csv(root / obs.SETTLING_SCALE_DIRS["all"]
                        / "convergence_summary_all_transects.csv")
        s["x"] = s["transect_id"].map(xmap)
        s = s.sort_values("x")
        years = (marked_years if direction == "forward"
                 else ref_end - MARKED_YEAR + 1)
        model = ("{0}–{1}".format(ref_start, MARKED_YEAR) if direction == "forward"
                 else "{0}–{1}".format(MARKED_YEAR, ref_end))
        ax = axes[i]
        below = np.zeros(len(s))
        fracs = []
        for tag, name, colour in YEARS_NEEDED_THRESHOLDS:
            y = s["years_needed_" + tag].to_numpy(float)
            ax.fill_between(s["x"], below, y, step="mid", color=colour, lw=0,
                            label="within " + name)
            fracs.append((name, 100.0 * np.mean(y <= years)))
            below = y
        box = dict(facecolor="white", edgecolor="none", pad=1.5, alpha=0.9)
        note = "{0} yr is enough at: ".format(years) + " · ".join(
            "{0:.0f}% ({1})".format(f, n) for n, f in fracs) + " of transects"
        notes.append((direction, note))
        if HIGHLIGHT_MODEL_WINDOW:
            ax.axhline(years, color=fs.C["ADDED"], lw=1.8, zorder=5,
                       path_effects=fs._halo(3.5))
            ax.text(n_tr - 5, years + 0.5,
                    "model window, {0} ({1} yr)".format(model, years),
                    ha="right", va="bottom", color=fs.C["ADDED"], fontsize=7.5,
                    fontweight="bold", zorder=6, bbox=box)
            ax.text(6, years + 0.5, note, ha="left", va="bottom", fontsize=7,
                    color=fs.INK, zorder=6, bbox=box)
        ax.set_ylim(0, ref_end - ref_start + 1.5)
        ax.set_xlim(0.5, n_tr + 0.5)
        ax.set_ylabel("years of record needed")
        ax.grid(True, axis="y", alpha=0.6)
        fs.town_bands(ax, label=(i == 0), spans=spans)
        fs._title(ax, i, "Windows starting {0} (the end moves later)".format(ref_start)
                  if direction == "forward"
                  else "Windows ending {0} (the start moves earlier)".format(ref_end))
    axes[-1].set_xlabel("transect (1–{0}, south → north)".format(n_tr))
    axes[-1].set_xticks([1] + list(range(100, n_tr + 1, 100)))
    # GIS domain labels under the transect axis: piecewise-linear between the
    # domains' centre transects, so each label sits over its own transects.
    dom = axes[-1].secondary_xaxis(
        -0.28, functions=(
            lambda x: np.interp(x, centre.to_numpy(), centre.index.to_numpy(float)),
            lambda d: np.interp(d, centre.index.to_numpy(float), centre.to_numpy())))
    dom.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    dom.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc="outside lower center", ncol=3,
               title="years needed before the rate stays")

    paths = fs.save(fig, out_dir / YEARS_NEEDED_STEM, close=True)
    fs.record_caption(paths[0],
        "How many years of CoastSat record each of the island's {n} transects "
        "needs before its shoreline change rate (OLS, the target's estimator) "
        "comes within a threshold of the {rs}–{re} rate and stays there for "
        "every longer window. (a) windows starting in {rs}, the end moving "
        "later; (b) windows ending in {re}, the start moving earlier. The "
        "bands stack from the loosest threshold (±1.0 m/yr, darkest, bottom) "
        "to the strictest (±0.25 m/yr, lightest, top): the top of each band is "
        "the years that threshold needs. {full} years means only the full "
        "record matches. {model}The x-axis "
        "counts transects south to north, each one unit wide; GIS domain "
        "numbers are given beneath it. Village spans are shaded.".format(n=len(xmap), rs=ref_start, re=ref_end,
                                   full=ref_end - ref_start + 1,
                                   model=("The amber line is the model's "
                                          "15-year window; where a band rises "
                                          "above it, 15 years was not enough at "
                                          "that threshold. The percentages are "
                                          "the share of transects where it was. "
                                          if HIGHLIGHT_MODEL_WINDOW else "")))
    return paths[0], notes


DOMAINS_README = """# {folder}/c-domain_means — the unit the model is graded on

The grading target is the domain MEAN of its transect rates
(`coastsat_domain_lrr.py` fits each transect and averages the slopes), so this
is the scale the answer actually has to be given at. `../a-eight_sites/` and
`../b-every_transect/` work on single transects. Averaging was expected to
settle sooner; it does not, so the disagreement between windows is real
shoreline behaviour rather than per-transect scatter.

No refitting happens here: every window of every transect is already in
`../b-every_transect/window_convergence_transects_all.csv`, and this is that table
grouped to the {n_domains} domains.

```
window_convergence_domains.csv       a row per domain per window
convergence_summary_domains.csv      a row per domain
```

Tables only. The figure is `../../years_needed_alongshore.png`.

## Two uncertainties, and the scoring uses the wider

`unc_m_yr` is the MEAN of the domain's transects' 95% half-widths.
`unc_if_independent_m_yr` is what the half-width of the mean would be if those
~10 transects were independent samples. They are not -- they are 10-metre-spaced
views of the same shoreline and move together -- so propagating them that way
divides the band by about √10 and manufactures a later convergence window
out of an assumption. Scoring uses the mean, which is **conservative**: it makes
convergence look earlier, not later. The independent figure is a column so the
size of the choice is visible.

## What this run found

{findings}
"""


# ============================================================
# READMEs
# ============================================================

DIRECTION_README = """# {folder} — {headline}

**When does each location settle on the long-term rate?** Fitted on the
CoastSat record **{ref_start}–{ref_end}**. {record_note}

{intro}

Both directions converge on the same reference, the {ref_start}–{ref_end}
rate, and that reference is the longest window of each sweep — fitted in
the same loop as every other window, so it cannot drift from a stored product.
The marked year {marked} is a real model window in both: forward it is
{ref_start}–{marked}, the window the model is graded on, and backward it is
{marked}–{ref_end}, the second leg of the canonical chain.

The figure for this sweep is `../years_needed_alongshore.png` (both
directions, one panel each).

```
a-eight_sites/     eight transects: shoreline position through time with the
                   fits drawn on it (the one figure here), plus tables
b-every_transect/  tables for all ~906 transects (the figure's source)
c-domain_means/    tables for the transects averaged to the 90 domains
```

The numbers below are the domain means under seven tolerances ("settles" = the
window's rate enters the tolerance and stays inside it for every longer window).
The figure uses only the three plain ones, ±0.25, ±0.5 and ±1.0 m/yr.

## The window this sweep gives

**{answer}**

{findings}

Producer:
`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_convergence.py`
`--direction {direction}`. Built {today} by interview (Hannah).
"""

SETTLING_README = """# {title}

**How many years of record does each place need before its shoreline change
rate stays within X m/yr of the {rs}–{re} rate?** One figure answers it:

![years needed](years_needed_alongshore.png)

`years_needed_alongshore.png` has one line of numbers per transect (906) along
the island, south to north. (a) is windows starting {rs}, with the end moving
later, and (b) is windows ending {re}, with the start moving earlier. The shaded
bands are three thresholds, ±1.0, ±0.5 and ±0.25 m/yr: the top of each band is
the years that threshold needs.{amber}

{notes}

"Stays within" means the window's rate is inside the threshold and every
longer window's rate is too. The first time a rate passes through the band
doesn't count, because most transects swing through it and back out.
{full} years means only the full record matches.

```
forward_from_{rs}/     windows {rs}–{rs1} … {rs}–{re}
backward_from_{re}/    windows {re4}–{re} … {rs}–{re}
    a-eight_sites/     shoreline position through time at 8 transects, with
                       the fits drawn on it (supporting figure) + tables
    b-every_transect/  per-transect tables (the figure's source)
    c-domain_means/    the same, averaged to the 90 model domains (tables)
```

Redraw without refitting:
`python scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_convergence.py --figure-only{flag}`

**Retired 2026-09-29** (Hannah: too many figures, an abstract "rate minus
reference" axis, and tolerance jargon). The figures that went are
`window_convergence_*` (eight sites and domains, the error curves),
`convergence_alongshore_*`, `tolerance_comparison_*` and
`domain_mean_vs_transect_*`, in both directions. Their tables are all kept.
PNGs are gitignored, but each one's PDF was committed in 3c6de274, under the old
layout: `window_convergence/record_<start>_<end>/<direction>/{{sites,all_transects,domain_means}}/supporting/`.
To see one: `git show 3c6de274:data/hatteras_init/5-scr/3-rates/coastsat/window_convergence/record_1996_2024/forward_from_1996/all_transects/supporting/convergence_alongshore_forward_from_1996.pdf > old.pdf`.
"""

RECORD_NOTE_FULL = (
    "The whole-profile companion question is in `../../1-rate_profiles/`.")
RECORD_NOTE_CUT = (
    "**This is an experiment**: the record is cut to {ref_start}–{ref_end}, "
    "so every window is scored against the {ref_start}–{ref_end} rate, not "
    "1996–2024. It asks whether the answer depends on the 2021 step. The main "
    "result is `../../../2-settling_window/`.")

SITES_README = """# {folder}/a-eight_sites — eight transects, in full

Eight evenly spaced domains (GIS {domain_list}), one transect each: the middle
one of its domain by alongshore order. The readable scale — every position,
every fit, one panel per site. `../b-every_transect/` runs the same sweep on all
of them and says whether these eight were representative.

```
shoreline_position_window_fits_{stem}.png
                                   the record: positions, the annual median,
                                   and six of the {n_windows} fits
window_convergence_transects.csv   a row per site per window
convergence_summary.csv            a row per site: reference, the marked
                                   window, six convergence entries
supporting/                        the PDFs and CAPTIONS.md
```

Six windows are drawn of the {n_windows} fitted. Drawing them all is a smear in
which the marked window and the reference — the two that matter — are
lost among the rest, each differing from its neighbour by one year of record
(Hannah, {today}).

## What this run found

{findings}
"""

ALL_README = """# {folder}/b-every_transect — the whole island

The same nested sweep on every CoastSat transect in the lookup
({n_transects} of them, {n_fits} fits), aggregated to the {n_domains} model
domains. It exists to answer one question about `../a-eight_sites/`: were the eight
representative, or did an even spread of eight happen to pick the unsettled
ones?

The figure drawn from these tables is `../../years_needed_alongshore.png`.

```
domain_convergence_summary.csv          a row per GIS domain: medians and
                                        quartiles across its ~10 transects
convergence_summary_all_transects.csv   a row per transect
window_convergence_transects_all.csv    the full sweep, {n_fits} rows
```

`years_needed_abs`, `_abs50` and `_abs100` are the ±0.25, ±0.5 and ±1.0 m/yr
columns the figure draws.

The domain number is the **median** across its transects, never the mean: one
transect that never settles would drag a mean to the end of the record and
report that the whole domain behaved that way.

## What this run found

{findings}
"""


def findings_text(summary, direction, scale="sites"):
    """The bullets a README carries, written from the run's own numbers."""
    lines = []
    med_ci = summary["years_needed_" + HEADLINE_TAG].median()
    lo = summary["years_needed_" + HEADLINE_TAG].min()
    hi = summary["years_needed_" + HEADLINE_TAG].max()
    n = len(summary)
    unit = {"sites": "sites", "domains": "domain means"}.get(scale, "transects")
    med_entry = summary["stable_entry_" + HEADLINE_TAG].median()
    lines.append(
        "**Headline (" + HEADLINE_LABEL + ").** The rate settles after a "
        "median of {0:.0f} years "
        "(range {1:.0f}–{2:.0f} across the {3} {4}), i.e. the window {5}."
        .format(med_ci, lo, hi, n, unit, window_label(direction, med_entry)))
    for tag, label, _test in CRITERIA:
        n_ok = int(summary["marked_in_" + tag].sum())
        lines.append(
            "Against {0}, {1} of {2} {3} ({4:.0f}%) have their {5} rate already "
            "inside the band; the median convergence window is {6}."
            .format(label, n_ok, n, unit, 100.0 * n_ok / n,
                    window_label(direction, MARKED_YEAR),
                    window_label(direction, summary["stable_entry_" + tag].median())))
    wandered = int((summary["first_entry_" + HEADLINE_TAG]
                    != summary["stable_entry_" + HEADLINE_TAG]).sum())
    lines.append(
        "{0} of {1} {2} ({3:.0f}%) enter the band and leave it again before "
        "settling, so the FIRST crossing is not the convergence window."
        .format(wandered, n, unit, 100.0 * wandered / n))
    flips = int((np.sign(summary["lrr_marked_m_yr"])
                 != np.sign(summary["ref_lrr_m_yr"])).sum())
    lines.append(
        "{0} of {1} {2} ({3:.0f}%) change SIGN between the {4} window and the "
        "reference — erosional over one and accretional over the other."
        .format(flips, n, unit, 100.0 * flips / n,
                window_label(direction, MARKED_YEAR)))
    worst = summary.loc[summary["diff_marked_m_yr"].abs().idxmax()]
    # The transect id only adds something when the unit IS a transect: at the
    # domain scale it is "GIS 36 mean" beside "GIS 36", which reads as a typo.
    named = ("GIS {0:.0f}".format(worst["domain_number"]) if scale == "domains"
             else "GIS {0:.0f} ({1})".format(worst["domain_number"],
                                             worst["transect_id"]))
    lines.append(
        "The largest {0} error is {1} at {2:+.2f} m/yr against a reference of "
        "{3:+.2f} m/yr."
        .format(window_label(direction, MARKED_YEAR), named,
                worst["diff_marked_m_yr"], worst["ref_lrr_m_yr"]))
    return "\n".join("- " + ln for ln in lines)


DIRECTION_TEXT = {
    "forward": (
        "how much record do you need from 1996?",
        "The START is pinned at {ref_start} and the END walks out, one year at a "
        "time: {ref_start}–{first} through {ref_start}–{ref_end}. It "
        "answers how much record the chain needs before the fitted rate stops "
        "depending on where it is cut off."),
    "backward": (
        "how late can a window begin?",
        "The END is pinned at {ref_end} and the START walks back, one year at a "
        "time: {last}–{ref_end} through {ref_start}–{ref_end}. It "
        "answers how recent a window can be and still recover the long-term "
        "rate — the other bracket on the same question."),
}


# ============================================================
# MAIN
# ============================================================

def one_direction(direction, domains_for_sites, scale, today):
    """Run and write one direction at whatever scales were asked for.

    The transect sweep is run ONCE and serves both the all_transects product
    and the domain means, which are a groupby of it rather than a refit.
    """
    windows = windows_for(direction)
    root = obs.window_convergence_dir(direction, pinned_year(direction),
                                      REF_START, REF_END)
    folder, stem = root.name, "{0}_from_{1}".format(direction, pinned_year(direction))
    headline = "Run a wider --scale for the island-wide answer."
    findings = ""
    summary_all = None

    if scale in ("sites", "every"):
        picks = pick_transects(domains_for_sites)
        print("{0} / sites: {1} windows x {2} transects = {3} fits"
              .format(direction, len(windows), len(picks), len(windows) * len(picks)))
        sweep = run_sweep(picks, windows, announce=True)
        summary = summarise(sweep, direction)

        out = root / obs.SETTLING_SCALE_DIRS["sites"]
        out.mkdir(parents=True, exist_ok=True)
        sweep.to_csv(out / obs.WINDOW_CONVERGENCE_SWEEP_FILE, index=False)
        summary.to_csv(out / obs.WINDOW_CONVERGENCE_SUMMARY_FILE, index=False)
        draw_fits(sweep, summary, out, direction)
        (out / "README.md").write_text(SITES_README.format(
            folder=folder, today=today, stem=stem, n_windows=len(windows),
            domain_list=", ".join(str(d) for d in domains_for_sites),
            findings=findings_text(summary, direction, "sites")), encoding="utf-8")
        print(summary[["domain_number", "ref_lrr_m_yr", "lrr_marked_m_yr",
                       "diff_marked_m_yr",
                       "stable_entry_" + HEADLINE_TAG]].to_string(index=False))
        print("wrote {0}\n".format(out))

    if scale in ("all", "domains", "every"):
        every = all_transects()
        print("{0} / every transect: {1} windows x {2} transects = {3} fits"
              .format(direction, len(windows), len(every), len(windows) * len(every)))
        sweep_all = run_sweep(every, windows)
        summary_all = summarise(sweep_all, direction)
        by_domain = summarise_domains(summary_all)

        if scale in ("all", "every"):
            out = root / obs.SETTLING_SCALE_DIRS["all"]
            out.mkdir(parents=True, exist_ok=True)
            sweep_all.to_csv(out / "window_convergence_transects_all.csv", index=False)
            summary_all.to_csv(out / "convergence_summary_all_transects.csv", index=False)
            by_domain.to_csv(out / "domain_convergence_summary.csv", index=False)
            findings = findings_text(summary_all, direction, "all")
            (out / "README.md").write_text(ALL_README.format(
                folder=folder, stem=stem, n_transects=len(summary_all),
                n_domains=len(by_domain), n_fits=len(sweep_all),
                findings=findings), encoding="utf-8")
            print("wrote {0}\n".format(out))
            headline = ("single transects settle at {0} (median over {1}, "
                        "{4}); {2:.0f}% have their {3} rate in that band"
                        .format(window_label(
                                    direction,
                                    summary_all["stable_entry_" + HEADLINE_TAG].median()),
                                len(summary_all),
                                100.0 * summary_all["marked_in_" + HEADLINE_TAG].mean(),
                                window_label(direction, MARKED_YEAR),
                                HEADLINE_LABEL))

        if scale in ("domains", "every"):
            dom_sweep = score_units(domain_mean_sweep(sweep_all))
            dom_summary = summarise(dom_sweep, direction)

            out = root / obs.SETTLING_SCALE_DIRS["domains"]
            out.mkdir(parents=True, exist_ok=True)
            dom_sweep.to_csv(out / "window_convergence_domains.csv", index=False)
            dom_summary.to_csv(out / "convergence_summary_domains.csv", index=False)
            dom_findings = findings_text(dom_summary, direction, "domains")
            (out / "README.md").write_text(DOMAINS_README.format(
                folder=folder, stem=stem, n_domains=len(dom_summary),
                findings=dom_findings), encoding="utf-8")
            print("wrote {0}\n".format(out))
            findings = dom_findings
            headline = ("DOMAIN MEANS settle at {0} ({5}, median over {1} "
                        "domains); {2:.0f}% have their {3} rate in that band. "
                        "Single transects settle at {4}"
                        .format(window_label(
                                    direction,
                                    dom_summary["stable_entry_" + HEADLINE_TAG].median()),
                                len(dom_summary),
                                100.0 * dom_summary["marked_in_" + HEADLINE_TAG].mean(),
                                window_label(direction, MARKED_YEAR),
                                window_label(
                                    direction,
                                    summary_all["stable_entry_" + HEADLINE_TAG].median()),
                                HEADLINE_LABEL))

    title, intro = DIRECTION_TEXT[direction]
    root.mkdir(parents=True, exist_ok=True)
    (root / "README.md").write_text(DIRECTION_README.format(
        folder=folder, headline=title, direction=direction, today=today,
        ref_start=REF_START, ref_end=REF_END, marked=MARKED_YEAR,
        record_note=RECORD_NOTE_FULL if (REF_START, REF_END) == (1996, 2024)
        else RECORD_NOTE_CUT.format(ref_start=REF_START, ref_end=REF_END),
        first=REF_START + MIN_WINDOW_YEARS - 1,
        last=REF_END - MIN_WINDOW_YEARS + 1,
        intro=intro.format(ref_start=REF_START, ref_end=REF_END,
                           first=REF_START + MIN_WINDOW_YEARS - 1,
                           last=REF_END - MIN_WINDOW_YEARS + 1),
        answer=headline, findings=findings), encoding="utf-8")
    return headline


def main(argv=None):
    global ABS_TOL_M_YR, REL_TOL, REF_START, REF_END

    ap = argparse.ArgumentParser(
        description="Two nested window sweeps: which windows recover the long-term rate?")
    ap.add_argument("--direction", choices=("forward", "backward", "both"),
                    default="both")
    ap.add_argument("--scale",
                    choices=("sites", "all", "domains", "every"), default="every",
                    help="sites: eight transects drawn in full. all: every "
                         "transect. domains: the domain MEAN, the unit the "
                         "model is graded on. every: all three.")
    ap.add_argument("--domains", type=int, nargs="+", default=SITE_DOMAINS,
                    help="the domains the `sites` scale draws, one transect each")
    ap.add_argument("--ref-start", type=int, default=REF_START,
                    help="first year of the record the sweep may see")
    ap.add_argument("--ref-end", type=int, default=REF_END,
                    help="last year of the record; the reference window is "
                         "<ref-start>-<ref-end> and both sweeps converge on it")
    ap.add_argument("--figure-only", action="store_true",
                    help="redraw years_needed_alongshore from the stored "
                         "tables without refitting")
    ap.add_argument("--abs-tol", type=float, default=ABS_TOL_M_YR)
    ap.add_argument("--rel-tol", type=float, default=REL_TOL)
    args = ap.parse_args(argv)

    ABS_TOL_M_YR, REL_TOL = args.abs_tol, args.rel_tol
    REF_START, REF_END = args.ref_start, args.ref_end
    if REF_END - REF_START + 1 < 2 * MIN_WINDOW_YEARS:
        raise SystemExit(
            "a record of {0} years cannot hold two sweeps of at least {1}"
            .format(REF_END - REF_START + 1, MIN_WINDOW_YEARS))
    print("record {0}-{1}; reference window {0}-{1}\n".format(REF_START, REF_END))
    today = pd.Timestamp.today().strftime("%Y-%m-%d")
    directions = (("forward", "backward") if args.direction == "both"
                  else (args.direction,))

    answers = {}
    if not args.figure_only:
        for d in directions:
            answers[d] = one_direction(d, args.domains, args.scale, today)
    # Needs both directions' every-transect tables, from this run or on disk.
    have = all((obs.window_convergence_dir(
                    d, REF_START if d == "forward" else REF_END, REF_START, REF_END)
                / obs.SETTLING_SCALE_DIRS["all"]
                / "convergence_summary_all_transects.csv").is_file()
               for d in ("forward", "backward"))
    if have:
        path, notes = draw_years_needed(REF_START, REF_END)
        print("wrote {0}".format(path))
        for d, note in notes:
            print("  {0:>8}: {1}".format(d, note))
        full = (REF_START, REF_END) == (1996, 2024)
        (path.parent / "README.md").write_text(SETTLING_README.format(
            title=("2-settling_window — how many years does each place need?"
                   if full else
                   "experiments/{0} — the settling sweep on a record cut to "
                   "{1}–{2}".format(path.parent.name, REF_START, REF_END)),
            rs=REF_START, re=REF_END, rs1=REF_START + MIN_WINDOW_YEARS - 1,
            re4=REF_END - MIN_WINDOW_YEARS + 1, full=REF_END - REF_START + 1,
            flag="" if full else " --ref-end {0}".format(REF_END),
            amber=(" The amber line is the model's 15-year window."
                   if HIGHLIGHT_MODEL_WINDOW else ""),
            notes="\n".join("- **{0}**: {1}".format(
                "(a) forward" if d == "forward" else "(b) backward", n)
                for d, n in notes)), encoding="utf-8")

    print("\n" + "=" * 70)
    for d, text in answers.items():
        print("{0:>8}: {1}".format(d, text))
    print("=" * 70)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
