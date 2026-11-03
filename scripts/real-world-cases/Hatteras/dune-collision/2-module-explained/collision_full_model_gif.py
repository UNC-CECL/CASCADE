"""
One domain's dunes through the full model, year by year and storm by storm: dune growth, overwash, collision.

    python collision_full_model_gif.py record [--gis 12]   # one hindcast run (~3 min), its dunes recorded
    python collision_full_model_gif.py gif    [--gis 12]   ->  figures/collision_6_full_model.gif (+ last-frame PNG)

record runs the hindcast as hat_run.yaml sets it (1996-2009, managed, edgeBE, groin, the collision module on at
its calibrated settings), unchanged, and only listens: inside each domain's Barrier3D year it notes the chosen
domain's dune crests (50 cells) at the end of the year before (after the management modules), after Barrier3D's
dune growth, after each storm's overwash (Barrier3D), after that storm's collision lowering (the module) and
after Barrier3D's dune diffusion. The hooks are the two DuneCollisionManager methods that stand in for
Barrier3D's DuneGrowth and DiffuseDunes (cascade/dune_collision.py), wrapped to record and then do exactly what
they did. Written to the workspace twin's data/ folder; the run goes under
output/raw_runs/dune-storm-response/2026-11-03-module-explained/. gif draws the record: one frame for each
year's growth, one for each storm that changes the dunes, one for each year's management.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-03
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import collision_explained_figures as E  # noqa: E402  (paths, settings, units, the storm series)
from site_layer.hat_figure_style import C, INK_MUTED, apply_style, caption, figsize, open_frame, save  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
GIS = 12
HINDCAST = E._HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
TOES = E.WORKSPACE / "data" / "hatteras_init" / "1-barrier3d-domains" / "1996-start" / "dune-toes" / "dune_toe_cells.csv"
DATA = E.workspace_twin(E._HERE.parent) / "data"
RUN_TAG = "2026-11-03-module-explained"
FIRST_YEAR = 1996               # Barrier3D time t is calendar year FIRST_YEAR + t
MIN_CHANGE_M = 0.005            # a storm frame is drawn when some cell moves at least this much
BIG_M = 1.0                     # a step that moves some cell this much is held longer
PAUSE_FRAMES = 3                # extra frames (merged into one longer frame) on those steps
FPS = 2.0
DPI = 100
C_GROW, C_OVER, C_COLL, C_MGMT, C_MIG = C["REF"], "0.45", E.C_MOD, C["ADDED"], C["LOCATOR"]
# -----------------------------------------------------------------------------


# The record's file for a domain
def record_path(gis):
    return DATA / f"gis{gis:02d}_dune_states.csv"


# Run the hindcast unchanged, noting the domain's crests at each step of each Barrier3D year
def record(gis):
    import cascade.dune_collision as DC

    toe = pd.read_csv(TOES)
    toe = ((toe[toe.gis == gis].sort_values("cell").toe_m_navd88.to_numpy()) - E.MHW) / 10.0  # dam MHW
    rows, state = [], {"mgr": None, "t": None}

    def note(mgr, dune_domain, t, stage, storm=None, r_high=None):
        b3d = mgr._barrier3d
        crest = dune_domain[t].max(axis=1) + b3d._BermEl
        rows.append(dict(year=FIRST_YEAR + t, t=t, stage=stage, storm=storm, r_high_m=r_high,
                         **{f"c{i:02d}": v for i, v in enumerate(E.m(crest))}))

    def mine(mgr):
        if state["mgr"] is None and np.allclose(mgr._toe_elevation, toe, atol=1e-6):
            state["mgr"] = mgr
        return mgr is state["mgr"]

    orig_growth, orig_update, orig_diffuse = DC.DuneCollisionManager.dune_growth, \
        DC.DuneCollisionManager.update, DC.DuneCollisionManager.diffuse_dunes

    def dune_growth(self, DuneDomain, t):
        if mine(self):
            note(self, DuneDomain, t - 1, "end_of_year")      # last year, after the management modules
        out = orig_growth(self, DuneDomain, t)
        if mine(self):
            note(self, DuneDomain, t, "growth")
            state["t"] = t
        return out

    def update(self, barrier3d, DuneDomain, t):
        if not mine(self):
            return orig_update(self, barrier3d, DuneDomain, t)
        storms = np.flatnonzero(barrier3d._StormSeries[:, 0] == t)
        k = self._storm
        r = float(barrier3d._StormSeries[storms[k], 1]) * 10 + E.MHW if k < len(storms) else None
        note(self, DuneDomain, t, "overwash", k + 1, r)
        out = orig_update(self, barrier3d, DuneDomain, t)
        note(self, DuneDomain, t, "collision", k + 1, r)
        return out

    def diffuse_dunes(self, DuneDomain, t):
        out = orig_diffuse(self, DuneDomain, t)
        if mine(self):
            note(self, out, t, "diffusion", self._storm, None)
        return out

    # Barrier3D's dune migration with the shoreline (after BRIE, before the management modules)
    from barrier3d import Barrier3d
    orig_migrate = Barrier3d.update_dune_domain

    def update_dune_domain(self):
        if state["mgr"] is None or self is not state["mgr"]._barrier3d:
            return orig_migrate(self)
        t = self._time_index          # the year being finished; update_dune_domain may advance it
        note(state["mgr"], self._DuneDomain, t, "before_migration")
        out = orig_migrate(self)
        note(state["mgr"], self._DuneDomain, t, "migration")
        return out

    DC.DuneCollisionManager.dune_growth = dune_growth
    DC.DuneCollisionManager.update = update
    DC.DuneCollisionManager.diffuse_dunes = diffuse_dunes
    Barrier3d.update_dune_domain = update_dune_domain

    os.environ.update({"HAT_RUN_KIND": "dune-storm-response", "HAT_RUN_TAG": f"{RUN_TAG}/gis{gis:02d}",
                       "HAT_SAVE_MODEL_STATE": "0", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                       "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    import runpy
    sys.path.insert(0, str(HINDCAST.parent))
    os.chdir(HINDCAST.parent)
    sys.argv = [str(HINDCAST)]
    try:
        runpy.run_path(str(HINDCAST), run_name="__main__")
    finally:
        if state["mgr"] is not None and state["t"] is not None:
            b3d = state["mgr"]._barrier3d
            note(state["mgr"], b3d.DuneDomain, state["t"], "end_of_year")
        if not rows:
            raise RuntimeError(f"GIS {gis}: no collision manager matched its toes; nothing recorded")
        DATA.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).round(4).to_csv(record_path(gis), index=False)
        print(f"recorded {len(rows)} states of GIS {gis} -> {record_path(gis)}")


# The record as frames: each year's growth, each storm that changes the dunes, each year's management
def frames(rec):
    cols = [c for c in rec.columns if c.startswith("c")]

    def arr(y, stage):
        s = y[y.stage == stage]
        return None if s.empty else s[cols].iloc[0].to_numpy(dtype=float)

    def moved(a, b):
        return np.nanmax(np.abs(a - b)) >= MIN_CHANGE_M

    out, prev = [], None
    for year, y in rec.groupby("year", sort=True):
        if prev is None:                                     # the start: the year before the first growth
            prev = arr(y, "end_of_year")
            if arr(y, "growth") is None:
                continue
        grown = arr(y, "growth")
        out.append(dict(year=year, kind="growth", before=prev, after=grown))
        prev = grown
        n_storms = int(y.storm.max()) if y.storm.notna().any() else 0
        for k in range(1, n_storms + 1):
            s = y[y.storm == k]
            a_ow, a_co, a_df = arr(s, "overwash"), arr(s, "collision"), arr(s, "diffusion")
            if a_ow is None or a_df is None:
                continue
            if moved(a_df, prev):
                out.append(dict(year=year, kind="storm", before=prev, overwash=a_ow, collision=a_co, after=a_df,
                                storm=k, r=s[s.stage == "overwash"].r_high_m.iloc[0], n_storms=n_storms))
            prev = a_df
        # the end of the year: Barrier3D's dune migration with the shoreline, then the management modules
        mig, end = arr(y, "migration"), arr(y, "end_of_year")
        if mig is not None:
            if moved(mig, prev):
                out.append(dict(year=year, kind="migration", before=prev, after=mig))
            prev = mig
        if end is not None:
            if moved(end, prev):
                out.append(dict(year=year, kind="mgmt", before=prev, after=end))
            prev = end
    return out


# One frame: the dunes as bars, what this step did, and the crest over the years
def draw(axes, fr, start, toe, hist, years):
    """Draw one frame; return the legend entries active in it."""
    a, b = axes
    active = {"dune", "toe", "start"}
    a.clear(); b.clear()
    x = np.arange(len(start)) * 10 + 5
    now, before = fr["after"], fr["before"]
    a.bar(x, now - E.BERM, bottom=E.BERM, width=9, color=E.C_SAND, lw=0)
    a.step(x, toe, where="mid", color="0.25", lw=1, ls=":", label="dune toe")
    a.step(x, start, where="mid", color="0.6", lw=0.9, label=f"crest at the start ({years[0]})")
    if fr["kind"] == "growth":
        a.bar(x, np.maximum(now - before, 0), bottom=before, width=9, color=C_GROW, lw=0, label="dune growth (Barrier3D)")
        head = f"{fr['year']}: dune growth, once a year before the storms (Barrier3D)"
        active.add("growth")
    elif fr["kind"] in ("mgmt", "migration"):
        col = C_MGMT if fr["kind"] == "mgmt" else C_MIG
        # a rebuild resets every cell to one design height; otherwise the road manager only adds sand
        rebuild = fr["kind"] == "mgmt" and np.nanstd(now) < 0.1
        if fr["kind"] == "migration":
            what, head = ("dune line moved with the shoreline (Barrier3D)",
                          f"End of {fr['year']}: the shoreline retreated a cell, so the dune line moved landward (Barrier3D)")
        elif rebuild:
            what, head = ("road dune rebuild",
                          f"End of {fr['year']}: a dune fell below the trigger, so the road manager rebuilds every "
                          f"dune to {np.nanmedian(now):.2f} m")
        else:
            what, head = ("overwash cleared off NC 12 onto the dune",
                          f"End of {fr['year']}: the road manager puts the overwash it clears off NC 12 on the dune")
        up, down = np.maximum(now - before, 0), np.maximum(before - now, 0)
        a.bar(x, up, bottom=before, width=9, color=col, lw=0)
        if (down > 0.005).any():
            a.bar(x, down, bottom=now, width=9, facecolor="none", edgecolor=col, hatch="////", lw=0.6)
            active.add("cut")
        active.add("migration" if fr["kind"] == "migration" else "mgmt")
    else:
        ow = np.maximum(before - fr["overwash"], 0)
        co = np.maximum(fr["overwash"] - fr["collision"], 0)
        a.bar(x, ow, bottom=before - ow, width=9, color=C_OVER, lw=0, label="lost to overwash (Barrier3D)")
        a.bar(x, co, bottom=fr["collision"], width=9, color=C_COLL, lw=0, label="lost to collision (module)")
        a.axhspan(0, fr["r"], color=E.C_WATER, alpha=0.45, lw=0, zorder=0)
        a.text(497, fr["r"] + 0.12, f"storm water level {fr['r']:.2f} m", ha="right", va="bottom", fontsize=8,
               color=C["LATE"], zorder=5, bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
        n_ow, n_co = int((ow > 0.005).sum()), int((co > 0.005).sum())
        head = (f"{fr['year']}, storm {fr['storm']} of {fr['n_storms']}: water at {fr['r']:.2f} m  ·  "
                f"{n_co} dunes lowered by collision, {n_ow} by overwash")
        active.add("water")
        if n_ow:
            active.add("overwash")
        if n_co:
            active.add("collision")
    a.set(xlim=(0, 500), ylim=(0, 10), xlabel=f"Alongshore in GIS {GIS} (m, south → north)",
          ylabel="Elevation (m NAVD88)")
    a.set_title(head, fontsize=9, loc="center")
    open_frame(a)

    # (b) the domain's median crest so far, every step
    t, med, cause = zip(*hist)
    b.plot(t, med, color="0.75", lw=1.0, zorder=1)
    b.scatter(t, med, s=9, c=[CAUSE_COLOR[c] for c in cause], lw=0, zorder=2)
    b.scatter(t[-1], med[-1], s=40, c=CAUSE_COLOR[cause[-1]], edgecolors=C["INK"], lw=0.8, zorder=3)
    b.axvline(t[-1], color="0.85", lw=0.8, zorder=0)
    b.set(xlim=(years[0], years[-1] + 1), ylim=(2.5, 8), ylabel="Median crest (m)")
    b.set_xticks(range(years[0], years[-1] + 2))
    b.set_xticklabels([str(y) if y % 2 == 0 else "" for y in range(years[0], years[-1] + 2)])
    open_frame(b)
    return active


# The legend: every entry always in its place, faded when not acting in the frame
# filled column by column: the dune, the references, Barrier3D, the module and migration, the road manager
LEGEND = (("dune", "dune (berm to crest)"), ("toe", "dune toe"),
          ("start", "crest at the start"), ("water", "storm water level"),
          ("growth", "growth (Barrier3D)"), ("overwash", "overwash (Barrier3D)"),
          ("collision", "collision (module)"), ("migration", "dune migration (Barrier3D)"),
          ("mgmt", "sand added (road management)"), ("cut", "cut down (road rebuild)"))
CAUSE_COLOR = {"growth": C_GROW, "collision": C_COLL, "overwash": C_OVER, "migration": C_MIG, "mgmt": C_MGMT,
               "none": "0.6"}


def legend_handles():
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    h = {"dune": Patch(color=E.C_SAND), "growth": Patch(color=C_GROW), "collision": Patch(color=C_COLL),
         "migration": Patch(color=C_MIG), "cut": Patch(facecolor="none", edgecolor=C_MGMT, hatch="////", lw=0.6),
         "toe": Line2D([], [], color="0.25", lw=1, ls=":"), "water": Patch(color=E.C_WATER, alpha=0.6),
         "overwash": Patch(color=C_OVER), "mgmt": Patch(color=C_MGMT), "start": Line2D([], [], color="0.6", lw=0.9)}
    return [h[k] for k, _ in LEGEND], [lab for _, lab in LEGEND]


# What moved the median in a frame, for the colour of its dot
def cause(f):
    if f["kind"] != "storm":
        return f["kind"]
    ow = np.nansum(np.maximum(f["before"] - f["overwash"], 0))
    co = np.nansum(np.maximum(f["overwash"] - f["collision"], 0))
    return "overwash" if ow > co else ("collision" if co > 0 else "none")


# Draw the record as a GIF and its last frame
def gif(gis):
    rec = pd.read_csv(record_path(gis))
    fr = frames(rec)
    cols = [c for c in rec.columns if c.startswith("c")]
    start = rec[cols].iloc[0].to_numpy(dtype=float)
    toe = pd.read_csv(TOES)
    toe = toe[toe.gis == gis].sort_values("cell").toe_m_navd88.to_numpy()
    years = sorted(int(y) for y in rec.year.unique())

    apply_style()
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", 0.58), layout="constrained",
                             gridspec_kw=dict(height_ratios=(3, 1)))
    handles, labels = legend_handles()
    leg = fig.legend(handles, labels, loc="outside upper center", ncol=5, frameon=False, fontsize=7.5,
                     columnspacing=1.2, handlelength=1.6)
    out = E.OUT / "collision_6_full_model.gif"
    hist, n = [], {"growth": 0, "storm": 0, "migration": 0, "mgmt": 0}
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(out), dpi=DPI):
        for i, f in enumerate(fr):
            # position in time: growth at the year's start, storms through it, management at its end
            tt = f["year"] + {"growth": 0.0, "migration": 0.97, "mgmt": 0.99}.get(f["kind"], 0.1 + 0.85 * (f.get("storm") or 0) / max(f.get("n_storms", 1), 1))
            hist.append((tt, float(np.nanmedian(f["after"])), cause(f)))
            n[f["kind"]] += 1
            active = draw(axes, f, start, toe, hist, years)
            for (key, _), txt, hd in zip(LEGEND, leg.get_texts(), leg.legend_handles):
                on = key in active
                txt.set_alpha(1.0 if on else 0.3)
                hd.set_alpha((0.6 if key == "water" else 1.0) if on else 0.2)
            writer.grab_frame()
            if np.nanmax(np.abs(f["after"] - f["before"])) >= BIG_M:   # hold the big moments
                for _ in range(PAUSE_FRAMES):
                    writer.grab_frame()
        for _ in range(8):                  # merged into one held frame
            writer.grab_frame()

    end = fr[-1]["after"]
    grow = sum(np.nansum(np.maximum(f["after"] - f["before"], 0)) for f in fr if f["kind"] == "growth") / len(start)
    ow = sum(np.nansum(np.maximum(f["before"] - f["overwash"], 0)) for f in fr if f["kind"] == "storm") / len(start)
    co = sum(np.nansum(np.maximum(f["overwash"] - f["collision"], 0)) for f in fr if f["kind"] == "storm") / len(start)
    mg = sum(np.nansum(f["after"] - f["before"]) for f in fr if f["kind"] == "mgmt") / len(start)
    mi = sum(np.nansum(f["after"] - f["before"]) for f in fr if f["kind"] == "migration") / len(start)
    caption(fig, f"**Tests:** nothing new; the last frame of `collision_6_full_model.gif`, one domain's dunes "
            f"through the full model. **How to read:** the hindcast as `hat_run.yaml` sets it ({years[0]}-"
            f"{years[-1]}, full management, edgeBE, the groin, the collision module at Cs {E.CS:g}, W0 "
            f"{E.W0:.0f} m), run unchanged with GIS {gis}'s 50 dune cells recorded at every step of each "
            f"Barrier3D year. (a) the dunes as bars from the berm to the crest; dotted, the toe; grey line, the "
            f"crest at the start. Frames in order: each year's dune growth (green, Barrier3D's logistic growth "
            f"before the storms); each storm that changes the dunes, with its water level (blue), what Barrier3D's "
            f"overwash took (grey) and what the module's collision took (purple), then Barrier3D's dune "
            f"diffusion; at the end of a year in which the shoreline retreated a full 10 m cell, Barrier3D's dune migration (teal: the front dune row drops off and the back row is taken from the interior, so the crest can jump up or down); and the end-of-year management (orange: the road manager puts the overwash it clears off NC 12 on the dune each year, and resets every dune cell of the domain to its design height once a front-row cell falls below its trigger). (b) the domain's median "
            f"crest after each step, each dot coloured by what moved it. The legend keeps every entry in place; "
            f"those not acting in a frame are faded. Storms that change nothing are skipped ({n['storm']} storm frames, "
            f"{n['growth']} growth, {n['migration']} migration, {n['mgmt']} management). **Shows:** over the run, per cell on average, growth "
            f"adds {grow:.2f} m, overwash takes {ow:.2f} m, collision {co:.2f} m, dune migration changes "
            f"{mi:+.2f} m and management {mg:+.2f} m; the median crest goes from {np.nanmedian(start):.2f} to {np.nanmedian(end):.2f} m.")
    save(fig, E.OUT / "collision_6_full_model", close=True)
    print(f"{out}: {len(fr)} frames ({n}); per cell growth +{grow:.2f}, overwash -{ow:.2f}, "
          f"collision -{co:.2f}, migration {mi:+.2f}, management {mg:+.2f} m; median {np.nanmedian(start):.2f} -> {np.nanmedian(end):.2f} m")


# Run: record or gif
def main() -> None:
    global GIS
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("action", choices=("record", "gif"))
    ap.add_argument("--gis", type=int, default=GIS)
    args = ap.parse_args()
    GIS = args.gis
    if args.action == "record":
        record(args.gis)
    else:
        gif(args.gis)


if __name__ == "__main__":
    main()
