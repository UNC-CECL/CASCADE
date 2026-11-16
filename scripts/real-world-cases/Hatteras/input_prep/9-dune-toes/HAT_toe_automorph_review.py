"""
Check Automorph's crest and toe by hand, in one window per domain (2026-11-08): his automated crest finder (3 m floor,
HAT_toe_method_automorph.auto_crest) suggests each crest, the adapted toe follows from it and is recomputed as you
click, and you accept, move or reject each. As his own workflow, an automated pass and then a review.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_automorph_review.py [--gis 22]
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_automorph_review.py --sample 50
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_automorph_review.py --all
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_automorph_review.py --fresh
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_automorph_review.py --report
    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_automorph_review.py --preview out.png [--gis 22]

Which transects are shown: where his crest is more than AGREE_M from your crest in crest_picks_1996.json, you rejected
that crest, or there is none. The rest have their crest accepted without a look (crest_status accepted_prior) and
their toe left unchecked. --sample N adds N of those, evenly along the island: the shown set alone is the hard cases,
so a toe score needs the sample. --all shows every transect. --fresh shows every transect from scratch, with your
earlier crests neither drawn nor used and nothing accepted without a look, saved to automorph_review_fresh_1996.json
(--report needs --fresh too). A survey with no crest picks shows every transect:
--year 2009 (the 2009 lidar on the 2009 start's dune windows, with the 2009 dune line and the 2008 aerial mosaic)
has none, and saves to 1-transects-2009/2-picks/automorph_review_2009.json.

Left click sets the crest (snaps to the highest point within 3 m, and the toe recalculates); shift + left click, or
left click after t, sets the toe. Right click rejects the transect (in toe mode, only its toe); again to undo. Middle
click resets the transect. Enter accepts what is left in the window and goes north.
Moving a crest clears a toe you had set on that transect, since the toe depends on it.

Saved at once to 2-picks/automorph_review_1996.json (D:), one entry per transect: crest_status (accepted_prior,
accepted, moved, rejected), crest_x, crest_z, toe_status (accepted, moved, rejected), toe_x, toe_z, and the automatic
auto_crest_x, auto_toe_x, auto_toe_z. Distances are m landward of the shoreline. crest_picks_1996.json is only read.
The window as last seen is saved per domain to 2-picks/automorph_review_views/. --report counts the decisions and
scores each method of 5-comparison/ against the checked toes; an accepted toe equals the adapted toe, so the moved
toes and the sample are what test it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-08
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

if "--year" in sys.argv:                                         # before HAT_toe_common reads it
    os.environ["HAT_TOE_YEAR"] = sys.argv[sys.argv.index("--year") + 1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_toe_method_automorph as A  # noqa: E402
import HAT_toe_pickers as P  # noqa: E402
from HAT_toe_common import (BERM_M, CHORD, FILL, MHW, M_COMPARE, PICKS, PROFILE, TOES, WATER, WINDOW, YEAR,  # noqa: E402
                            crest_picks)

# --- CONFIG ------------------------------------------------------------------
FRESH = "--fresh" in sys.argv   # every transect, from scratch: your earlier crests are neither shown nor used
REVIEW_FILE = PICKS / f"automorph_review{'_fresh' if FRESH else ''}_{YEAR}.json"
VIEWS = PICKS / f"automorph_review{'_fresh' if FRESH else ''}_views"
AGREE_M = 5.0              # his crest within this of your crest needs no look
VIEW_MIN_M = 220           # the profile runs at least this far landward, so a taller dune behind a crest shows
SNAP_M = 3                 # a click snaps to the highest point within this
AM, AM_DARK, YOURS, RED = "#7B3F99", "#3B1566", "#E8A200", "#B03A2E"
HELP = ("left click = crest here   shift + left click (or t, then left click) = toe here   right click = reject "
        "(in toe mode, the toe)   middle click = reset\n"
        "enter = accept the rest and go north   n / → = north   p / ← = south   t = crest / toe mode   q = save and quit")
# -----------------------------------------------------------------------------
P.CREST = AM               # the map marks the automatic crest in Automorph's colour


def arrays(r, x0):
    z = np.array([np.nan if v is None else v for v in r["z"]], float)
    return x0 + np.arange(len(z), dtype=float), z, -x0


def z_at(r, x0, x):
    i = int(round(x)) - x0
    v = r["z"][i] if 0 <= i < len(r["z"]) else None
    return None if v is None else float(v)


def toe_for(r, x0, crest_x):
    """The adapted toe for a crest at crest_x: ((x, z) or None, the reason when none, the crest as the method used it)."""
    xs, z, shore = arrays(r, x0)
    a = A.automorph(xs, z, shore, pick=dict(crest_x=crest_x), picks=True, real_peaks=True)
    zs = A.smooth(z)
    pt = lambda i: None if i is None else (float(xs[i]), float(zs[i]))
    return pt(a["toe"]), a["reason"], pt(a["crest"])


def prepare(D, profiles):
    """Each transect with its automatic crest (fc) and the crest you reviewed before (prior: x, or None), by domain.
    The profiles are the unrounded ones of run(): the viewer data's 2 decimals move some peaks and toes."""
    x0, cp, by = D["x0"], crest_picks(), {}
    for r0 in D["transects"]:
        if not r0.get("z"):
            continue
        r = dict(r0)
        s, z_full, p = profiles[r["id"]][:3]
        grid = x0 + np.arange(len(r["z"]), dtype=float)
        zz = np.interp(grid, s - s[p["shore"]], z_full, left=np.nan, right=np.nan)
        r["z"] = [None if not np.isfinite(v) else float(v) for v in zz]
        xs, z, shore = arrays(r, x0)
        c = A.auto_crest(z, shore)
        r["fc"] = None if c is None else [float(xs[c]), float(A.smooth(z)[c])]
        e = {} if FRESH else cp.get(str(r["id"]), {})
        r["prior_status"], r["prior"] = e.get("status"), e.get("crest_x")
        by.setdefault(r["gis"], []).append(r)
    return by


def why_shown(r):
    """Why a transect needs a look, '' where his crest agrees with the crest you reviewed."""
    if r["fc"] is None:
        return "no automatic crest"
    if r["prior_status"] is None:
        return "no earlier review"
    if r["prior_status"] == "rejected":
        return "you rejected the crest here"
    d = r["fc"][0] - r["prior"]
    return "" if abs(d) <= AGREE_M else f"{abs(d):.0f} m {'landward' if d > 0 else 'seaward'} of your crest"


def even_sample(rows, n):
    """n transects spread evenly along the island."""
    rows = sorted(rows, key=lambda r: (r["gis"], r["id"]))
    if n >= len(rows):
        return rows
    return [rows[i] for i in np.unique(np.linspace(0, len(rows) - 1, n).round().astype(int))]


def seed_agreed(rv, by, x0):
    """Accept, without a look, his crest where it agrees with the one you reviewed. Returns how many were new."""
    new = 0
    for rs in by.values():
        for r in rs:
            k = str(r["id"])
            if k not in rv and not why_shown(r):
                rv[k] = dict(gis=r["gis"], crest_status="accepted_prior", crest_x=r["fc"][0],
                             crest_z=z_at(r, x0, r["fc"][0]), auto_crest_x=r["fc"][0],
                             updated=datetime.now().isoformat(timespec="seconds"))
                new += 1
    return new


def _base(r, rv, x0):
    e = dict(rv.get(str(r["id"]), {}))
    auto = r["fc"][0] if r["fc"] else None
    t = toe_for(r, x0, auto)[0] if auto is not None else None
    return e | dict(gis=r["gis"], auto_crest_x=auto, auto_toe_x=t[0] if t else None, auto_toe_z=t[1] if t else None,
                    updated=datetime.now().isoformat(timespec="seconds"))


def crest_of(r, rv):
    """(x, z) of the crest in use, or None."""
    e = rv.get(str(r["id"]), {})
    if e.get("crest_status") == "rejected":
        return None
    if e.get("crest_x") is not None:
        return e["crest_x"], e.get("crest_z")
    return tuple(r["fc"]) if r["fc"] else None


def toe_of(r, rv, x0):
    """((x, z) or None, the reason when none) of the toe in use."""
    e = rv.get(str(r["id"]), {})
    if e.get("toe_status") == "rejected":
        return None, "toe rejected"
    if e.get("toe_status") == "moved":
        return (e["toe_x"], e["toe_z"]), ""
    c = crest_of(r, rv)
    if c is None:
        return None, "crest rejected" if e.get("crest_status") == "rejected" else "no crest"
    return toe_for(r, x0, c[0])[:2]


def set_crest(rv, r, x0, x):
    """A clicked crest: the highest point within SNAP_M of x. Clears a toe set for the old crest."""
    xs, z, _ = arrays(r, x0)
    k = np.flatnonzero((np.abs(xs - x) <= SNAP_M) & (xs > 0) & np.isfinite(z))
    if not len(k):
        return False
    j = k[np.argmax(z[k])]
    e = _base(r, rv, x0) | dict(crest_status="moved", crest_x=float(xs[j]), crest_z=float(z[j]))
    for f in ("toe_status", "toe_x", "toe_z"):
        e.pop(f, None)
    rv[str(r["id"])] = e
    return True


def _fix_crest(e, r, rv, x0):
    """A toe decision also confirms a crest not yet decided."""
    if e.get("crest_status") in (None, "accepted_prior") and r["fc"]:
        e |= dict(crest_status="accepted", crest_x=r["fc"][0], crest_z=z_at(r, x0, r["fc"][0]))
    return e


def set_toe(rv, r, x0, x):
    z = z_at(r, x0, x)
    if z is None:
        return False
    rv[str(r["id"])] = _fix_crest(_base(r, rv, x0), r, rv, x0) | dict(toe_status="moved", toe_x=float(round(x)), toe_z=z)
    return True


def reject(rv, r, x0, toe=False):
    """Reject the transect (toe=True: only its toe); on one already rejected, undo."""
    k, e = str(r["id"]), rv.get(str(r["id"]), {})
    if toe and e.get("toe_status") == "rejected":
        for f in ("toe_status", "toe_x", "toe_z"):
            e.pop(f, None)
    elif toe:
        rv[k] = _fix_crest(_base(r, rv, x0), r, rv, x0) | dict(toe_status="rejected", toe_x=None, toe_z=None)
    elif e.get("crest_status") == "rejected":
        rv.pop(k)
    else:
        rv[k] = _base(r, rv, x0) | dict(crest_status="rejected", crest_x=None, crest_z=None, toe_status=None, toe_x=None,
                                        toe_z=None)


def accept_rest(rv, rs, x0):
    """Enter: every undecided crest and toe in the window stands as shown."""
    for r in rs:
        e = _base(r, rv, x0)
        if e.get("crest_status") == "rejected":
            continue
        e = _fix_crest(e, r, rv, x0)
        if e.get("crest_status") is None:                        # no crest to accept: leave it undecided
            continue
        if e.get("toe_status") is None:
            t = toe_for(r, x0, e["crest_x"])[0]
            e |= dict(toe_status="accepted", toe_x=t[0], toe_z=t[1]) if t else dict(toe_status="none", toe_x=None,
                                                                                    toe_z=None)
        rv[str(r["id"])] = e


def _is_checked(e):
    return e.get("toe_status") in ("accepted", "moved", "rejected", "none") or e.get("crest_status") == "rejected"


def review(gis=None, every=False, sample=0, preview=None):
    import matplotlib
    if preview:
        matplotlib.use("Agg", force=True)
    else:
        try:
            import tkinter  # noqa: F401
            matplotlib.use("TkAgg", force=True)
        except Exception:
            print("[warn] TkAgg unavailable; the review needs an interactive backend.")
    import matplotlib.pyplot as plt
    plt.switch_backend(matplotlib.get_backend())
    vd = TOES / f"viewer_data_{YEAR}.json"
    if not vd.exists():
        P.viewer_data(*P.run())
    D = json.loads(vd.read_text(encoding="utf-8"))
    x0 = D["x0"]
    print("reading the profiles (about a minute)...", flush=True)
    by = prepare(D, P.run(write=False)[2])
    rv = json.loads(REVIEW_FILE.read_text(encoding="utf-8")) if REVIEW_FILE.exists() else {}

    def save():
        REVIEW_FILE.parent.mkdir(parents=True, exist_ok=True)
        REVIEW_FILE.write_text(json.dumps(rv, indent=1), encoding="utf-8")

    if not preview and seed_agreed(rv, by, x0):
        save()
    agreed = [r for rs in by.values() for r in rs if not why_shown(r)]
    extra = {r["id"] for r in even_sample(agreed, sample)} if sample else set()
    shown = lambda g: [r for r in by.get(g, []) if every or why_shown(r) or r["id"] in extra]
    order = [g for g in range(1, 91) if shown(g)]
    print(f"{sum(len(shown(g)) for g in order)} transects to check in {len(order)} domains "
          f"({len(agreed)} agree with your crests, {len(extra)} of those in the sample)")
    start = gis if gis in order else next((g for g in order if any(not _is_checked(rv.get(str(r["id"]), {}))
                                                                   for r in shown(g))), order[0] if order else None)
    i = order.index(start) if start in order else len(order)
    while 0 <= i < len(order):
        g = order[i]
        rs = sorted(shown(g), key=lambda r: -(r.get("o") or [0, 0])[1])          # north at the top, as on the map
        n = len(rs[0]["z"])
        xs = np.arange(x0, x0 + n)
        xend = max(VIEW_MIN_M, P.view_end(rs), *[(r["prior"] or 0) + 30 for r in rs])
        vis = xs <= xend
        zs = [z for r in rs for z, v in zip(r["z"], vis) if z is not None and v]
        ylo, yhi = np.floor(min([-0.5] + zs)), np.ceil(max([4.0] + zs) + 0.5)
        fig = plt.figure(figsize=(18, 2.3 * max(len(rs), 2) + 1.6))
        gsp = fig.add_gridspec(len(rs), 3, width_ratios=(1, 1, 2.3))
        ax_map = fig.add_subplot(gsp[:, 0])
        ax_photo = fig.add_subplot(gsp[:, 1], sharex=ax_map, sharey=ax_map)
        axes = [fig.add_subplot(gsp[k, 2]) for k in range(len(rs))]
        state = dict(go=None, mode="crest", laid_out=False)
        mp = P.domain_map(rs, x0, xend)
        sh = lambda r: mp["off"].get(r["id"], 0.0) if mp else 0.0
        xlim = (mp["c"][0], mp["c"][-1]) if mp else (xs[0], xend)

        def map_toes(r):
            return (toe_of(r, rv, x0)[0], None)

        def draw():
            flags = {str(r["id"]): dict(status="picked") for r in rs
                     if rv.get(str(r["id"]), {}).get("crest_status") == "moved"}
            kept = {ax: (ax.get_xlim(), ax.get_ylim()) for ax in [ax_map, ax_photo, *axes]} if state["laid_out"] else {}
            P.draw_map(ax_map, mp, rs, lambda r: crest_of(r, rv), map_toes, flags, xend)
            P.draw_map(ax_photo, mp, rs, lambda r: crest_of(r, rv), map_toes, flags, xend, photo=True)
            for ax, r in zip(axes, rs):
                ax.clear()
                e, o = rv.get(str(r["id"]), {}), sh(r)
                at = lambda q: (q[0] + o, q[1])
                z = np.array([np.nan if v is None else v for v in r["z"]], float)
                ax.fill_between(xs + o, z, ylo, color=FILL, lw=0)
                ax.plot(xs + o, z, color=PROFILE, lw=1.2)
                if r.get("w"):
                    ax.axvspan(r["w"][0] + o, r["w"][1] + o, color=WINDOW, alpha=0.12, lw=0)
                ax.plot(o, z[-x0], "o", ms=5, color=WATER, zorder=5)
                for xd in (mp["dlx"].get(r["id"], []) if mp else []):
                    ax.axvline(xd, color=P.DUNELINE_COLOR, lw=1.2, ls=(0, (5, 3)), zorder=4)
                if r["prior"] is not None:
                    pz = z_at(r, x0, r["prior"])
                    if pz is not None:
                        ax.plot(r["prior"] + o, pz + 0.45, "v", ms=7, mfc="none", mec=YOURS, mew=1.3, zorder=6)
                moved = e.get("crest_status") == "moved"
                if r["fc"]:
                    ax.plot(*at(r["fc"]), "^", ms=9, color=AM, alpha=0.35 if moved else 1, zorder=7)
                c, (t, why) = crest_of(r, rv), toe_of(r, rv, x0)
                if moved and c:
                    ax.plot(*at(c), "^", ms=10, color=YOURS, mec="white", mew=1, zorder=8)
                if c and c[1] is not None and t:
                    ax.plot([o, c[0] + o], [z[-x0], c[1]], color=CHORD, lw=0.7, ls=(0, (4, 2)))
                if t:
                    ax.plot(*at(t), "s", ms=8, color=YOURS if e.get("toe_status") == "moved" else AM_DARK, mec="white",
                            mew=0.7, zorder=9)
                ax.axhline(MHW, color=WATER, lw=0.7)
                ax.axhline(BERM_M, color="0.45", lw=0.6, ls=(0, (2, 2)))
                ax.set_ylim(ylo, yhi)
                ax.set_xlim(xlim[1], xlim[0])                            # the ocean on the right, as on the map
                f = lambda p: "—" if p is None or p[1] is None else f"{p[1]:.2f} m"
                ax.set_title(f"Transect {r['id']}   crest {f(c)}   toe {f(t)}" + (f" ({why})" if why and not t else ""),
                             loc="left", fontsize=9.5)
                note = "" if FRESH else why_shown(r) or ("sample" if r["id"] in extra else "")
                if note:
                    ax.text(0.99, 0.94, note, transform=ax.transAxes, ha="right", va="top", fontsize=8, color=RED)
                cs, ts = e.get("crest_status") or "unchecked", e.get("toe_status") or "unchecked"
                bad = cs == "rejected" or ts == "rejected"
                ax.set_title(f"crest {cs.upper()}  ·  toe {ts.upper()}", loc="right", fontsize=9, fontweight="bold",
                             color=RED if bad else AM_DARK if _is_checked(e) else "0.4")
                ax.set_facecolor("#FBEBEB" if bad else "white")
                ax.set_ylabel("m NAVD88", fontsize=8)
                ax.tick_params(labelsize=8)
            for ax in axes[:-1]:
                ax.tick_params(labelbottom=False)
            axes[-1].set_xlabel("m landward of the most seaward shoreline (as the maps; dot = this transect's shoreline)",
                                fontsize=9)
            done = sum(_is_checked(rv.get(str(r["id"]), {})) for r in rs)
            fig.suptitle(f"GIS {g}   CLICK MODE: {state['mode'].upper()}   ({done}/{len(rs)} checked; domain {i + 1} of "
                         f"{len(order)})\n{HELP}\npurple triangle = Automorph's crest, yellow = yours    "
                         + ("" if FRESH else "open yellow triangle above = the crest you reviewed earlier    ")
                         + "dark square = toe, yellow = yours    teal "
                         f"band = the model's dune window    pink dashed = {P.DUNELINE_YEAR} dune line",
                         fontsize=9.5, x=0.01, ha="left")
            for ax, (xl, yl) in kept.items():                            # a redraw keeps your zoom
                ax.set_xlim(xl)
                ax.set_ylim(yl)
            fig.canvas.draw_idle()

        def leave(go):
            """The window as last seen goes to VIEWS, once, on the way out."""
            VIEWS.mkdir(parents=True, exist_ok=True)
            fig.savefig(VIEWS / f"automorph_review_gis{g:02d}.png", dpi=110)
            state["go"] = go
            plt.close(fig)

        def changed():
            save()
            draw()

        def on_click(ev):
            if ev.inaxes not in axes or ev.xdata is None or (toolbar and toolbar.mode):   # zooming, not picking
                return
            r = rs[axes.index(ev.inaxes)]
            toe_mode = state["mode"] == "toe" or ev.key == "shift"
            x = ev.xdata - sh(r)
            if ev.button == 1 and (set_toe if toe_mode else set_crest)(rv, r, x0, x):
                changed()
            elif ev.button == 3:
                reject(rv, r, x0, toe=toe_mode)
                changed()
            elif ev.button == 2 and str(r["id"]) in rv:
                rv.pop(str(r["id"]))
                changed()

        def on_key(ev):
            if ev.key in ("enter", "return"):
                accept_rest(rv, rs, x0)
                save()
                leave(+1)
            elif ev.key == "t":
                state["mode"] = "toe" if state["mode"] == "crest" else "crest"
                draw()
            elif ev.key in ("n", "right"):
                leave(+1)
            elif ev.key in ("p", "left"):
                leave(-1)
            elif ev.key in ("q", "escape"):
                leave(0)

        toolbar = getattr(fig.canvas.manager, "toolbar", None)
        if not preview:                                              # full screen, before the layout is set for it
            try:
                fig.canvas.manager.window.state("zoomed")
                fig.canvas.manager.window.update()
            except Exception as e:
                print(f"[warn] could not maximize the window: {e}")
        draw()
        fig.tight_layout(rect=(0, 0, 1, 1 - 1.0 / fig.get_size_inches()[1]))
        if preview:
            fig.savefig(preview, dpi=110)
            print(f"wrote {preview} (GIS {g}, {len(rs)} transects; nothing saved to the review file)")
            plt.close(fig)
            return
        state["laid_out"] = True
        draw()
        fig.canvas.mpl_connect("button_press_event", on_click)
        fig.canvas.mpl_connect("key_press_event", on_key)
        plt.show()
        if not state["go"]:
            break
        i += state["go"]
    save()
    print(f"saved {REVIEW_FILE}")
    report(rv)


def report(rv=None):
    rv = rv if rv is not None else json.loads(REVIEW_FILE.read_text(encoding="utf-8"))
    e = pd.DataFrame(rv).T
    print("crests:", e.crest_status.value_counts().to_dict())
    print("toes:  ", e.get("toe_status", pd.Series(dtype=object)).fillna("unchecked").value_counts().to_dict())
    d = e[e.crest_status.isin(["accepted", "moved"])]
    if len(d):
        print(f"of the {len(d)} crests looked at and kept, Automorph's stood on "
              f"{100 * (d.crest_status == 'accepted').mean():.0f}%")
    if "toe_status" not in e:
        return
    k = e[e.toe_status.isin(["accepted", "moved"])].copy()
    f = M_COMPARE / f"toe_methods_transects_{YEAR}.csv"
    if not len(k) or not f.exists():
        return
    k.index = k.index.astype(int)
    c = pd.read_csv(f).set_index("transect")
    hand = pd.to_numeric(k.toe_z)
    moved = k.toe_status == "moved"
    print(f"against your {len(k)} checked toes ({int(moved.sum())} moved); accepted toes equal the adapted toe:")
    for m in c.columns:
        if m == "gis" or m.endswith(("_x", "_crest")):
            continue
        dz = (c[m].reindex(hand.index) - hand).dropna()
        if len(dz):
            print(f"   {m:20s} n {len(dz):3d}  RMSE {np.sqrt((dz ** 2).mean()):.2f} m  mean {dz.mean():+.2f} m")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gis", type=int, default=None, help="the domain to start at")
    ap.add_argument("--year", type=int, choices=(1996, 2009), default=1996, help="the survey (default 1996)")
    ap.add_argument("--all", action="store_true", help="every transect")
    ap.add_argument("--fresh", action="store_true", help="every transect from scratch, in its own file (read above)")
    ap.add_argument("--sample", type=int, default=0, help="add this many transects whose crest agreed")
    ap.add_argument("--report", action="store_true", help="count the decisions and score the methods")
    ap.add_argument("--preview", default=None, help="save the first window to this PNG and stop")
    a = ap.parse_args()
    if a.report:
        return report()
    review(a.gis, a.all or a.fresh, a.sample, a.preview)


if __name__ == "__main__":
    main()
