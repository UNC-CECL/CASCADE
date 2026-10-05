"""
How close is each nested window's alongshore rate profile to 1996-2024? r, bias and RMSE.

    python scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_window_r_bias_rmse.py

Reads the per-transect window fits that coastsat_window_profiles.py writes (no
refitting), scores every window against 1996-2024 with a 95% domain bootstrap and
draws the r figure and the bias + RMSE figure, both directions together and one
per direction. Details: scripts/input_prep/5-scr/3-rates/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-02
"""

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

# The sibling sweep, for the reference years and the pinned year of each direction
sys.path.insert(0, str(Path(__file__).resolve().parent))
import coastsat_window_convergence as wc              # noqa: E402

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402


# --- CONFIG ------------------------------------------------------------------
REF_START, REF_END = wc.REF_START, wc.REF_END

# The r levels, the domain bootstrap; the model windows live in coastsat_window_convergence.MODEL_WINDOWS
# r levels marked on the r figures (Hannah, 2026-10-01): conventional, not data-picked
R_LEVELS = (0.5, 0.75, 0.9)
N_BOOT = 1000
BOOT_SEED = 20261001
R_YMIN = -0.3
# Bias and RMSE axes cut so the model window reads; the shortest windows run past
BIAS_Y_HALF = 2.5
RMSE_Y_MAX = 6.0
# -----------------------------------------------------------------------------


# Years in the model window for a direction, per direction since 2026-10-02 (Hannah)
def model_years(direction):
    w = wc.model_window(direction)
    return w[1] - w[0] + 1


# A window as 'start-end'
def window_label(start, end):
    return "{0}–{1}".format(start, end)


# Window against reference by window length: r (shape), bias (level), RMSE (size of the miss)

# Transect x window-length matrix of rates for a direction, from its transects CSV
def _rate_matrix(direction):
    t = pd.read_csv(obs.window_profiles_dir(direction, wc.pinned_year(direction), REF_START, REF_END)
                    / "window_profiles_transects.csv")
    m = t.pivot_table(index=["transect_id", "domain_number"], columns="n_years",
                      values="lrr_m_yr")
    return m.sort_index()


# r, mean difference and RMSE of each column of m against ref, NaN-safe
def _scores(m, ref):
    out = np.full((3, m.shape[1]), np.nan)
    for j in range(m.shape[1]):
        ok = np.isfinite(m[:, j]) & np.isfinite(ref)
        if ok.sum() > 2:
            diff = m[ok, j] - ref[ok]
            out[0, j] = np.corrcoef(m[ok, j], ref[ok])[0, 1]
            out[1, j] = diff.mean()
            out[2, j] = np.sqrt((diff ** 2).mean())
    return out


# The three scores per window length and direction, with block-bootstrap 95% intervals
def window_scores():
    fwd, bwd = _rate_matrix("forward"), _rate_matrix("backward")
    fwd, bwd = fwd.align(bwd, join="inner")
    lengths = np.array(fwd.columns, dtype=int)
    ref = fwd[REF_END - REF_START + 1].to_numpy()
    mats = {"forward": fwd.to_numpy(), "backward": bwd.to_numpy()}
    # Domains are the resampled blocks: neighbouring transects are not independent
    domains = fwd.index.get_level_values("domain_number").to_numpy()
    blocks = [np.flatnonzero(domains == d) for d in np.unique(domains)]
    rng = np.random.default_rng(BOOT_SEED)
    boot = {k: [] for k in mats}
    for _ in range(N_BOOT):
        idx = np.concatenate([blocks[i] for i in rng.integers(0, len(blocks), len(blocks))])
        for k, m in mats.items():
            boot[k].append(_scores(m[idx], ref[idx]))
    rows = []
    for k, m in mats.items():
        est = _scores(m, ref)
        lo, hi = np.nanpercentile(np.array(boot[k]), [2.5, 97.5], axis=0)
        for j, L in enumerate(lengths):
            row = {"direction": k, "n_years": int(L)}
            for i, name in enumerate(("r", "bias_m_yr", "rmse_m_yr")):
                row[name] = round(float(est[i, j]), 4)
                row[name + "_lo95"] = round(float(lo[i, j]), 4)
                row[name + "_hi95"] = round(float(hi[i, j]), 4)
            row["n_transects"] = len(fwd)
            rows.append(row)
    return pd.DataFrame(rows)


# First window length from which r stays at or above a level
def first_lasting(sub, level):
    above = (sub["r"] >= level).to_numpy()
    for i in range(len(above)):
        if above[i:].all():
            return int(sub["n_years"].iloc[i])
    return None


# Colour, marker, legend label and the window family for a direction
def _r_style(direction):
    if direction == "forward":
        return (fs.C["EARLY"], "o", "start pinned at {0}".format(REF_START),
                "{0}–{1} … {0}–{2}".format(REF_START, REF_START + 1, REF_END))
    return (fs.C["LATE"], "s", "end pinned at {0}".format(REF_END),
            "{1}–{0} … {2}–{0}".format(REF_END, REF_END - 1, REF_START))


# Where a label sits from its point, per panel and direction, clear of the other curve
LABEL_PLACE = {
    ("r", "forward"): "below-right", ("r", "backward"): "above-left",
    ("bias_m_yr", "forward"): "below-right", ("bias_m_yr", "backward"): "above-right",
    ("rmse_m_yr", "forward"): "below-left", ("rmse_m_yr", "backward"): "above-right",
}


# Model-window labels lifted clear on a leader line, where the neighbourhood is crowded
MODEL_LABEL_LIFT = {("r", "backward"): (-30, 19)}
# On 1996-2026 the forward 20-yr label lands on the 19-yr r-level label
if REF_END == 2026:
    MODEL_LABEL_LIFT[("r", "forward")] = (12, -42)


# A circled, labelled point; `lift` moves the label further out on a leader line
def _mark(ax, L, v, col, place, text, lift=None):
    ax.plot(L, v, marker="o", ms=6.5, mfc="none", mec=col, mew=1.0, zorder=4)
    above, right = place.startswith("above"), place.endswith("right")
    offset = lift or (5 if right else -6, 6 if above else -9)
    arrow = dict(arrowstyle="-", color=col, lw=0.6, shrinkA=0, shrinkB=4) if lift else None
    ax.annotate(text, (L, v), xytext=offset, textcoords="offset points",
                ha="left" if right else "right", va="bottom" if above else "top",
                fontsize=8, color=col, zorder=6, path_effects=fs._halo(2.0),
                arrowprops=arrow)


# Key points (Hannah, 2026-10-01): the model window everywhere, and on r where it stays >= each level
def _label_key_points(ax, sub, col, place, col_name, fmt, lift=None, model_L=None):
    at = lambda L: float(sub.loc[sub["n_years"] == L, col_name].iloc[0])
    _mark(ax, model_L, at(model_L), col, place,
          "{0} yr: {1}".format(model_L, fmt.format(at(model_L))), lift)
    if col_name != "r":
        return
    # Every circled point says its value (Hannah, 2026-10-01: that is why it is circled)
    for L in dict.fromkeys(first_lasting(sub, lv) for lv in R_LEVELS):
        if L is not None and L != model_L:
            _mark(ax, L, at(L), col, place, "{0} yr: {1}".format(L, fmt.format(at(L))))


# Calendar windows along the top of a one-direction figure
def _calendar_axis(ax, direction):
    ticks = [5, 10, 15, 20, 25, REF_END - REF_START + 1]
    if direction == "forward":
        labels = ["{0}–{1}".format(REF_START, REF_START + L - 1) for L in ticks]
    else:
        labels = ["{0}–{1}".format(REF_END - L + 1, REF_END) for L in ticks]
    top = ax.secondary_xaxis("top")
    top.set_xticks(ticks)
    top.set_xticklabels(labels, fontsize=7.5)
    top.set_xlabel("calendar window ({0} pinned at {1})".format(
        "start" if direction == "forward" else "end", wc.pinned_year(direction)))


# Panel specs: column, y label, title, label format, y limits
SCORE_PANELS = {
    "r": ("alongshore Pearson r", "Alongshore correlation of windowed shoreline change rates with the {0} rate", "r = {0:.2f}",
          (R_YMIN, 1.1)),
    "bias_m_yr": ("bias (m/yr)", "Mean bias relative to the {0} rate", "{0:+.2f} m/yr",
                  (-BIAS_Y_HALF, BIAS_Y_HALF)),
    "rmse_m_yr": ("RMSE (m/yr)", "Root-mean-square error relative to the {0} rate",
                  "{0:.2f} m/yr", (0, RMSE_Y_MAX)),
}


# Legend text for the dashed model-window line(s)
def model_legend(directions):
    return "model window ({0})".format(", ".join(
        "{0}–{1}".format(*wc.model_window(d)) for d in directions))


# Stacked panels on one window-length axis, one per score in `names`
def draw_scores(table, directions, names, out_path):
    fs.apply_style()
    heights = {"r": 1.3, "bias_m_yr": 1.0, "rmse_m_yr": 1.0}
    fig, axes = plt.subplots(len(names), 1, sharex=True, layout="constrained",
                             figsize=fs.figsize("double", height=1.2 + 2.4 * len(names)),
                             gridspec_kw=dict(height_ratios=[heights[n] for n in names]))
    axes = np.atleast_1d(axes)
    n_ref = REF_END - REF_START + 1
    for i, (ax, name) in enumerate(zip(axes, names)):
        ylab, title, fmt, ylim = SCORE_PANELS[name]
        for d in directions:
            col, mk, label, _ = _r_style(d)
            sub = table[table["direction"] == d]
            ax.plot(sub["n_years"], sub[name], color=col, lw=1.2, marker=mk, ms=2.8,
                    zorder=3, label=label if i == 0 else None)
            ax.fill_between(sub["n_years"], sub[name + "_lo95"], sub[name + "_hi95"],
                            color=col, alpha=0.15, lw=0, zorder=1,
                            label="95% interval (domains resampled)" if i == 0 else None)
            _label_key_points(ax, sub, col, LABEL_PLACE[(name, d)], name, fmt,
                              MODEL_LABEL_LIFT.get((name, d)), model_years(d))
        # One dashed line per model-window length, one legend entry naming the windows
        for k, L in enumerate(dict.fromkeys(model_years(d) for d in directions)):
            ax.axvline(L, color=fs.INK, lw=0.7, ls=(0, (4, 2)), zorder=2,
                       label=model_legend(directions) if i == 0 and k == 0 else None)
        ax.axhline(0.0, color=fs.C["INK_MUTED"], lw=0.5, zorder=1)
        ax.set_ylim(*ylim)
        ax.set_ylabel(ylab)
        ax.grid(True, alpha=0.5)
        text = title.format(window_label(REF_START, REF_END))
        if len(names) > 1:
            fs._title(ax, i, text)
        else:
            ax.set_title(text, loc="center")
        if name == "r":
            for k, lv in enumerate(R_LEVELS):
                ax.axhline(lv, color=fs.C["INK_MUTED"], lw=0.6, ls=(0, (1, 2)), zorder=2,
                           label="r = " + ", ".join("{0:g}".format(v) for v in R_LEVELS)
                           if k == 0 else None)
    axes[-1].set_xlim(1.5, n_ref + 0.5)
    axes[-1].set_xlabel("window length L (years of record)")
    if len(directions) == 1:
        _calendar_axis(axes[0], directions[0])
    n_entries = 2 * len(directions) + 1 + ("r" in names)
    fig.legend(loc="outside lower center", frameon=False,
               ncol=n_entries if len(directions) == 1 else 3)
    paths = fs.save(fig, out_path, close=True)
    fs.record_caption(paths[0], _scores_caption(table, directions, names))
    return paths[0]


# Caption for a score figure: what each panel is, the model-window values, the method
def _scores_caption(table, directions, names):
    at = lambda sub, L, c: float(sub.loc[sub["n_years"] == L, c].iloc[0])
    ref = window_label(REF_START, REF_END)
    letters = {n: ("({0}) ".format(chr(ord("a") + i)) if len(names) > 1 else "")
               for i, n in enumerate(names)}
    defs = {
        "r": "{0}Alongshore Pearson r: whether the hotspots are in the same places "
             "(shape only).",
        "bias_m_yr": "{0}Bias: the mean over transects of window rate minus " + ref +
                     " rate; negative means the window is more erosional than the "
                     "long-term record.",
        "rmse_m_yr": "{0}RMSE: the root-mean-square of the same difference, the typical "
                     "size of the miss at one transect, sign ignored (RMSE² = bias² + "
                     "scatter²).",
    }
    parts = []
    for d in directions:
        L = model_years(d)
        _, _, label, family = _r_style(d)
        sub = table[table["direction"] == d]
        vals = []
        if "r" in names:
            vals.append("r = {0:.2f} ({1:.2f}–{2:.2f})".format(
                at(sub, L, "r"), at(sub, L, "r_lo95"), at(sub, L, "r_hi95")))
        if "bias_m_yr" in names:
            vals.append("bias {0:+.2f} m/yr ({1:+.2f} to {2:+.2f})".format(
                at(sub, L, "bias_m_yr"), at(sub, L, "bias_m_yr_lo95"),
                at(sub, L, "bias_m_yr_hi95")))
        if "rmse_m_yr" in names:
            vals.append("RMSE {0:.2f} m/yr ({1:.2f}–{2:.2f})".format(
                at(sub, L, "rmse_m_yr"), at(sub, L, "rmse_m_yr_lo95"),
                at(sub, L, "rmse_m_yr_hi95")))
        text = "{c}: windows with the {lab} ({fam}). At {L} years ({mw}, the model window) {v}.".format(
            mw="{0}–{1}".format(*wc.model_window(d)),
            c="Red" if d == "forward" else "Blue", lab=label, fam=family, L=L,
            v=", ".join(vals))
        if "r" in names:
            text += " r stays at or above " + "; ".join(
                "{0:g} from {1} years".format(lv, first_lasting(sub, lv)) for lv in R_LEVELS) + "."
        parts.append(text)
    notes = []
    if "r" in names:
        notes.append("The dotted lines are r = {0}, conventional levels, not ones the data "
                     "pick out; circled points are the model window and the first window "
                     "length from which r stays at or above each level.".format(
                         ", ".join("{0:g}".format(v) for v in R_LEVELS)))
    else:
        notes.append("Circled points are the model window.")
    cut = [("±{0:g} m/yr".format(BIAS_Y_HALF), "bias_m_yr"), ("{0:g} m/yr".format(RMSE_Y_MAX), "rmse_m_yr")]
    cut = [c for c, n in cut if n in names]
    if cut:
        notes.append("The y-axes are cut at {0}; the shortest windows run past.".format(
            " and ".join(cut)))
    if len(directions) == 1:
        notes.append("The top axis gives the calendar window.")
    return (
        "How the CoastSat shoreline change rate profile of each window (OLS, all {n} "
        "transects, domain 1 at Cape Point to 90 at Pea Island) compares with the {ref} "
        "profile, against window length. {defs} {parts} The windows are nested in the "
        "reference, so every score reaches its perfect value at {nr} years by "
        "construction. Shading is the 95% interval from {nb} bootstrap resamples of the "
        "90 domains (domains, not transects, because neighbouring transects move "
        "together). The dashed vertical line marks the model window ({mw}). {notes}".format(
            n=int(table["n_transects"].iloc[0]), ref=ref,
            defs=" ".join(defs[n].format(letters[n]) for n in names),
            parts=" ".join(parts), nr=REF_END - REF_START + 1, nb=N_BOOT,
            mw="; ".join("{0}–{1}, {2} yr".format(*wc.model_window(d), model_years(d))
                         for d in directions), notes=" ".join(notes)))


# The score figures (Hannah, 2026-10-01: r alone, bias + RMSE together), each both
# directions together in 2-r_bias_rmse/ and one per direction in its folder
SCORE_FIGURES = {
    "window_profiles_r": ("r",),
    "window_profiles_bias_rmse": ("bias_m_yr", "rmse_m_yr"),
}


def score_figures():
    table = window_scores()
    top = obs.window_scores_dir(ref_start=REF_START, ref_end=REF_END)
    top.mkdir(parents=True, exist_ok=True)
    table.to_csv(top / "window_profiles_r_bias_rmse.csv", index=False)
    for stem, names in SCORE_FIGURES.items():
        print(draw_scores(table, ("forward", "backward"), names, top / stem))
        for d in ("forward", "backward"):
            out_dir = obs.window_scores_dir(d, wc.pinned_year(d), REF_START, REF_END)
            out_dir.mkdir(parents=True, exist_ok=True)
            print(draw_scores(table, (d,), names,
                              out_dir / "{0}_{1}_from_{2}".format(stem, d, wc.pinned_year(d))))
    for d in ("forward", "backward"):
        print(table[(table["direction"] == d)
                    & (table["n_years"] == model_years(d))].T.to_string())


if __name__ == "__main__":
    score_figures()
