#!/usr/bin/env python3
r"""
HAT_imagery_review_gui.py
==============================================================================
The aerial-imagery review of the 1984 footprint as a window you work in,
instead of 62 PNGs and a spreadsheet. Same data, same overlays, same sheet as
HAT_imagery_review_1984.py - this only changes how the judgement is entered,
and lets the reviewer MEASURE the one thing the verdict rests on.

WHAT THE WINDOW SHOWS
    Left: the 1984 and 1997 photographs of one domain side by side (a third
    panel, the 1 m lidar, on request), sharing one view so pan and zoom in
    either moves both (the matplotlib toolbar below them; the mouse wheel
    zooms about the pointer). Under them the brightness strip. Every overlay
    - the two digitized dune lines, NC-12, interior row 0, the road rows,
    v3's rows, the seaward alternative, the domain box, your picks - is a
    checkbox, so the photograph can be looked at bare and the lines brought
    back.
    BLINK MODE puts both photographs in ONE panel and flips between them on
    Space (or B), or on a timer. The eye catches movement between two frames
    of the same view far better than side by side, which matters for the
    one-cell offsets most of the footprint is made of.
    Right: the domain's numbers from the footprint table, the pick buttons,
    the verdict form (the six columns of imagery_review_1984.csv, with the
    vocabulary as drop-downs), Save, Prev / Next, and Summarize.

THE PICKS (measured by the reviewer, not by code)
    Three features per year, each one click on the photograph after its button:
        toe    the seaward vegetation line, beach sand -> dune vegetation: the
               feature the digitized dune lines trace
        back   the landward edge of the dune band, hummocky dune sand -> flat
               vegetated interior
        road   the seaward edge of the pavement AS IT IS IN THAT PHOTOGRAPH
    Each click is taken on the profile nearest it and stored as a position and
    as metres landward of interior row 0 (negative = seaward). From them, per
    year, the three bands the placement question is about:
        dune_band     back - toe        the dune field
        back_to_road  road - back       the strip between the dune and the road
        toe_to_road   road - toe        the whole crest-to-road space
    and once both years are picked, their change 1997 - 1984, the shift of each
    feature (positive where the 1984 feature lay seaward, the footprint's sign),
    and - taking N as given (Hannah, 2026-09-09) - where the lost or gained
    width sat:
        lost_dune_band     = -(d_dune_band)
        lost_back_to_road  = -(d_back_to_road)
        lost_behind_road   = N x 10 m - the two above      (the remainder)
    `band_suggests` names the largest share (dune_band / back_to_road /
    behind_road / none when N is 0) and is DERIVED; the verdict is yours.
    Nothing is snapped: the number is where you clicked, on the photograph as
    georeferenced (stated accuracy 1.2 m). A domain without a model road still
    takes the road pick if a road is visible; the bands that need it stay blank
    otherwise.

WHAT IT WRITES
    On Save (Ctrl+S, or "Save & next"): the six verdict columns and the
    measured columns of the domain on screen, plus reviewed_by / reviewed_at,
    into 2-domain-reconstruction-1984/3-placement/imagery-review/imagery_review_1984.csv. Nothing else in the sheet
    is touched, and HAT_imagery_review_1984.py keeps these columns when it
    re-runs. "Summarize" runs HAT_imagery_review_summary.py on the sheet as
    it stands.

KEYS   Right / Left  next / previous domain      Ctrl+S  save
       Space or B    flip the year in blink mode   Esc     cancel a pick
       1 2 3         pick toe / back / road on the year shown (blink mode)
       (keys are ignored while the cursor is in a text box)

PERFORMANCE
    A domain takes ~5-10 s to read from the drive the first time. The arrays
    are cached under ~/.cascade/imagery_review_cache/ (outside the repo, keyed
    on domain, year, window and resolution), and the next domain in the list
    is read in the background while you look at the current one. Delete the
    cache folder to force a re-read (e.g. after changing the merge rule in
    the batch script).

USAGE
    python HAT_imagery_review_gui.py                 # every domain in the sheet
    python HAT_imagery_review_gui.py --domains 85,63 # a subset
    python HAT_imagery_review_gui.py --years 1984,1997,1996
    python HAT_imagery_review_gui.py --smoke         # open, draw one, screenshot, close
==============================================================================
"""
from __future__ import annotations

import argparse
import getpass
import sys
import threading
import tkinter as tk
import warnings
from tkinter import ttk
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

sys.path.insert(0, str(Path(__file__).resolve().parent))
import HAT_imagery_review_1984 as R  # noqa: E402   the batch script: data, overlays, sheet
off = R.off

CACHE_DIR = Path.home() / ".cascade" / "imagery_review_cache"
YEAR_COLOUR = {1984: off.LINE_STYLE[1984]["color"], 1997: off.LINE_STYLE[1997]["color"]}
PICK_YEARS = R.PICK_YEARS          # the two the sheet has columns for
PICK_KINDS = R.PICK_KINDS
PICK_LABEL = {"toe": "toe (sand -> dune)", "back": "back of dune (dune -> interior)", "road": "road edge"}
PICK_MARKER = {"toe": "o", "back": "s", "road": "D"}
BLINK_MS = 700
OVERLAYS = [                       # key, label, default on
    ("line84", "1984 dune line (digitized)", True),
    ("line97", "1997 dune line (digitized)", True),
    ("roads", "NC-12 1984 (dashed) and 2004 (solid)", True),
    ("row0", "interior row 0 (crest + 1)", True),
    ("roadrows", "road rows: measured (dotted), placed in v3 (dark)", True),
    ("v3", "v3 rows: behind / in front of NC-12", True),
    ("sea", "alternative: the same rows at the dune", True),
    ("box", "domain box", True),
    ("picks", "your picks (toe, back of dune, road)", True),
    ("stripmarks", "lines and rows on the strip", True),
]


def _num(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


# =============================================================================
# DATA
# =============================================================================

class Data:
    """Everything the window needs, loaded once; per-domain arrays cached."""

    def __init__(self, years: list[int]):
        self.years = years
        self.tab = pd.read_csv(R.FOOTPRINT_CSV).set_index("domain")
        print("loading the map layers ...")
        self.gdf = off.load_domains()
        self.lines = off.load_lines(self.gdf.crs)
        self.roads = off.m.load_roads(self.gdf.crs, clip_to=self.gdf.union_all())
        print("indexing the imagery ...")
        self.imagery = {y: R.Imagery(y, self.gdf.crs, strict=False) for y in years}   # cache carries a lost drive
        self.frames = R.ProfileFrames()
        self._cache: dict[int, dict] = {}
        self._lock = threading.Lock()
        CACHE_DIR.mkdir(parents=True, exist_ok=True)

    def sheet(self) -> pd.DataFrame:
        if R.SHEET.is_file():
            df = pd.read_csv(R.SHEET, dtype={c: str for c in R.VERDICT_COLS + R.MEASURED_COLS})
            return df.set_index("domain")
        return pd.DataFrame(columns=["role"] + R.VERDICT_COLS + R.MEASURED_COLS)

    def geometry(self, d: int) -> dict:
        t = self.tab.loc[d]
        pr = self.frames.get(d)
        pl = R.placements(t)
        bounds = self.gdf[self.gdf["domain_id"].astype(int) == d].total_bounds
        win = R.domain_window(pr, pl, bounds)
        return {"t": t, "pr": pr, "pl": pl, "win": win, "bounds": bounds}

    def load(self, d: int) -> dict:
        """Photographs (and their strips) for one domain, from memory, disk or the drive."""
        with self._lock:
            if d in self._cache:
                return self._cache[d]
        g = self.geometry(d)
        x_lo, x_hi, y_lo, y_hi = g["win"]
        x_med = float(g["pr"]["interior_x"].median())
        r_min = -int(R.SEAWARD_STRIP_M // R.CELL_M)
        r_max = int(np.ceil((x_med - x_lo) / R.CELL_M)) + 1
        photos, strips, cover = {}, {}, {}
        for y, im in self.imagery.items():
            key = CACHE_DIR / (f"GIS{d}_{y}_{R.DISPLAY_RES_M:g}m_{x_lo:.0f}_{x_hi:.0f}_{y_lo:.0f}_{y_hi:.0f}.npz")
            if key.is_file():
                img = np.load(key)["img"]
            elif im.files:
                img = im.read(x_lo, x_hi, y_lo, y_hi, R.DISPLAY_RES_M, self.gdf.crs)
                np.savez_compressed(key, img=img)
            else:                                       # drive off and nothing cached: a blank panel, not a crash
                W, H = int(round((x_hi - x_lo) / R.DISPLAY_RES_M)), int(round((y_hi - y_lo) / R.DISPLAY_RES_M))
                img = np.zeros((H, W, 3), dtype=np.uint8)
            L = R.luminance(img)
            ok = np.isfinite(L)
            cover[y] = float(ok.mean())
            if ok.any():
                lo, hi = np.nanpercentile(L, R.STRETCH_PCT)
                shown = np.clip((img.astype(np.float32) - lo) / max(hi - lo, 1.0), 0, 1)
                shown[~ok] = 0.93
                cells, S = R.strip(L, x_lo, y_hi, R.DISPLAY_RES_M, g["pr"], r_min, r_max)
                strips[y] = (cells, (S - lo) / max(hi - lo, 1.0))
            else:
                shown = None
            photos[y] = shown
        out = {**g, "photos": photos, "strips": strips, "cover": cover, "lidar": None}
        with self._lock:
            self._cache[d] = out
            if len(self._cache) > 8:                      # keep memory bounded
                for k in list(self._cache)[:-8]:
                    if k != d:
                        del self._cache[k]
        return out

    def lidar(self, d: int):
        rec = self.load(d)
        if rec["lidar"] is None:
            rec["lidar"] = off.load_1m(self.gdf, [d])
        return rec["lidar"]


# =============================================================================
# THE MEASUREMENT BEHIND A PICK
# =============================================================================

def measure_pick(rec: dict, year: int, kind: str, x: float, y: float) -> dict:
    """One click -> the columns for that feature and year, on the nearest profile."""
    pr = rec["pr"]
    i = int(np.argmin(np.abs(pr["interior_y"].to_numpy() - y)))
    p = pr.iloc[i]
    from_row0 = float(p["interior_x"]) - x                 # + = landward of row 0
    k = f"{kind}{str(year)[2:]}"
    return {f"{k}_x": round(x, 1), f"{k}_y": round(y, 1), f"{k}_profile": int(p["profile"]),
            f"{k}_from_row0_m": round(from_row0, 1)}


def derived(m: dict, n_cells: int) -> dict:
    """Bands per year, their change, and where the width sat given N."""
    out = {}
    pos = {}
    for y in PICK_YEARS:
        yy = str(y)[2:]
        pos[y] = {k: _num(m.get(f"{k}{yy}_from_row0_m")) for k in PICK_KINDS}
        t, b, r = pos[y]["toe"], pos[y]["back"], pos[y]["road"]
        out[f"dune_band{yy}_m"] = round(b - t, 1) if np.isfinite(t) and np.isfinite(b) else ""
        out[f"back_to_road{yy}_m"] = round(r - b, 1) if np.isfinite(b) and np.isfinite(r) else ""
        out[f"toe_to_road{yy}_m"] = round(r - t, 1) if np.isfinite(t) and np.isfinite(r) else ""
    y0, y1 = PICK_YEARS
    for k in PICK_KINDS:
        a, b = pos[y0][k], pos[y1][k]
        out[f"{k}_shift_m"] = round(b - a, 1) if np.isfinite(a) and np.isfinite(b) else ""
    out["toe_shift_cells"] = (int(np.trunc(out["toe_shift_m"] / R.CELL_M))
                              if out["toe_shift_m"] != "" else "")
    for band in ("dune_band", "back_to_road", "toe_to_road"):
        a, b = _num(out[f"{band}{str(y0)[2:]}_m"]), _num(out[f"{band}{str(y1)[2:]}_m"])
        out[f"d_{band}_m"] = round(b - a, 1) if np.isfinite(a) and np.isfinite(b) else ""
    # where the width sat, taking N as given: lost = 1984 - 1997 (positive where
    # the 1984 island was wider there); the remainder is behind the road
    total = n_cells * R.CELL_M
    ld = -_num(out["d_dune_band_m"])
    lb = -_num(out["d_back_to_road_m"])
    out["lost_dune_band_m"] = round(ld, 1) if np.isfinite(ld) else ""
    out["lost_back_to_road_m"] = round(lb, 1) if np.isfinite(lb) else ""
    out["lost_behind_road_m"] = round(total - ld - lb, 1) if np.isfinite(ld) and np.isfinite(lb) else ""
    if n_cells == 0:
        out["band_suggests"] = "none" if np.isfinite(ld) else ""
    elif np.isfinite(ld) and np.isfinite(lb):
        shares = {"dune_band": ld, "back_to_road": lb, "behind_road": total - ld - lb}
        sgn = 1 if n_cells > 0 else -1                 # a removal domain GAINED width: the most negative loss
        out["band_suggests"] = max(shares, key=lambda k_: sgn * shares[k_])
    else:
        out["band_suggests"] = ""
    return out


# =============================================================================
# THE WINDOW
# =============================================================================

class App:
    def __init__(self, data: Data, ids: list[int], smoke: bool = False):
        self.data = data
        self.ids = ids
        self.i = 0
        self.smoke = smoke
        self.pick_target: tuple[int, str] | None = None  # (year, kind) being picked
        self.picks: dict[tuple[int, str], dict] = {}     # (year, kind) -> the columns of that pick
        self.blink_year: int | None = None
        self.blink_images: dict[int, object] = {}
        self._blink_job = None
        self.axes: list = []
        self.all_axes: list = []
        self.axs = None

        self.root = tk.Tk()
        self.show_lidar = tk.BooleanVar(value=False)
        self.blink = tk.BooleanVar(value=False)
        self.auto_blink = tk.BooleanVar(value=False)
        self.root.title("HAT imagery review - the 1984 footprint against the photographs")
        self.root.geometry("1760x1000")
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

        left = ttk.Frame(self.root)
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        right = ttk.Frame(self.root, width=380, padding=8)
        right.pack(side=tk.RIGHT, fill=tk.Y)
        right.pack_propagate(False)

        off.apply_style()
        self.fig = Figure(figsize=(12.5, 9.5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=left)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, left, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.canvas.mpl_connect("button_press_event", self.on_click)

        self.overlay_vars: dict[str, tk.BooleanVar] = {}
        self.artists: dict[str, list] = {}
        self.build_right(right)
        self.root.bind("<Right>", lambda e: self._key(lambda: self.step(+1)))
        self.root.bind("<Left>", lambda e: self._key(lambda: self.step(-1)))
        self.root.bind("<space>", lambda e: self._key(self.flip))
        self.root.bind("<b>", lambda e: self._key(self.flip))
        self.root.bind("<Escape>", lambda e: self.cancel_pick())
        for k_, kind in enumerate(PICK_KINDS, 1):
            self.root.bind(str(k_), lambda e, kind_=kind: self._key(
                lambda: self.start_pick(self.blink_year or PICK_YEARS[0], kind_)))
        self.root.bind("<Control-s>", lambda e: self.save())
        self.status("ready")
        self.show(0)

    # ---- the right-hand panel ------------------------------------------------
    def build_right(self, f):
        nav = ttk.Frame(f)
        nav.pack(fill=tk.X)
        ttk.Button(nav, text="< Prev", command=lambda: self.step(-1)).pack(side=tk.LEFT)
        self.combo = ttk.Combobox(nav, state="readonly", width=28)
        self.combo.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
        self.combo.bind("<<ComboboxSelected>>", lambda e: self.show(self.combo.current()))
        ttk.Button(nav, text="Next >", command=lambda: self.step(+1)).pack(side=tk.LEFT)
        self.refresh_combo()

        self.info = tk.Text(f, height=11, width=46, wrap="word", relief="flat",
                            font=("Consolas", 9), background="#f4f4f4")
        self.info.pack(fill=tk.X, pady=(8, 4))
        self.info.configure(state="disabled")

        ttk.Label(f, text="View", font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(4, 0))
        ttk.Checkbutton(f, text="blink mode: one panel, Space flips the year",
                        variable=self.blink, command=lambda: self.show(self.i)).pack(anchor="w")
        ttk.Checkbutton(f, text=f"auto-blink every {BLINK_MS / 1000:.1f} s",
                        variable=self.auto_blink, command=self._auto_blink).pack(anchor="w")
        ttk.Checkbutton(f, text="lidar panel (1 m, the surface the rows are cut from)",
                        variable=self.show_lidar, command=lambda: self.show(self.i)).pack(anchor="w")

        ttk.Label(f, text="Overlays", font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(6, 0))
        for key, label, on in OVERLAYS:
            v = tk.BooleanVar(value=on)
            self.overlay_vars[key] = v
            ttk.Checkbutton(f, text=label, variable=v, command=lambda k_=key: self.toggle(k_)).pack(anchor="w")

        ttk.Label(f, text="Measure: three features per year, one click each",
                  font=("Segoe UI", 10, "bold"), wraplength=360).pack(anchor="w", pady=(8, 0))
        pk = ttk.Frame(f)
        pk.pack(fill=tk.X)
        for c_, kind in enumerate(PICK_KINDS, 1):
            ttk.Label(pk, text=PICK_LABEL[kind], font=("Segoe UI", 8), wraplength=105,
                      justify="center").grid(row=0, column=c_, padx=2)
        for r_, y in enumerate(PICK_YEARS, 1):
            ttk.Label(pk, text=str(y), foreground=YEAR_COLOUR.get(y, "0.2"),
                      font=("Segoe UI", 9, "bold")).grid(row=r_, column=0, padx=(0, 4))
            for c_, kind in enumerate(PICK_KINDS, 1):
                ttk.Button(pk, text=f"pick {kind}", width=10,
                           command=lambda y_=y, k_=kind: self.start_pick(y_, k_)).grid(row=r_, column=c_, padx=2, pady=1)
        ttk.Button(pk, text="Clear", command=self.clear_picks).grid(row=1, column=4, rowspan=2, padx=(6, 0))
        self.measure_var = tk.StringVar(value="no picks yet")
        ttk.Label(f, textvariable=self.measure_var, wraplength=360, font=("Consolas", 9)).pack(anchor="w", pady=(2, 0))

        ttk.Label(f, text="Verdict", font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(8, 0))
        self.fields: dict[str, ttk.Combobox] = {}
        for c in R.VERDICT_COLS:
            if c == "notes":
                continue
            row = ttk.Frame(f)
            row.pack(fill=tk.X, pady=1)
            ttk.Label(row, text=c, width=18).pack(side=tk.LEFT)
            cb = ttk.Combobox(row, values=[""] + [v.strip() for v in R.VERDICT_VOCAB[c].split("|")],
                              state="readonly", width=20)
            cb.pack(side=tk.LEFT, fill=tk.X, expand=True)
            self.fields[c] = cb
        ttk.Label(f, text="notes").pack(anchor="w", pady=(4, 0))
        self.notes = tk.Text(f, height=3, width=46, wrap="word", font=("Segoe UI", 9))
        self.notes.pack(fill=tk.X)

        btns = ttk.Frame(f)
        btns.pack(fill=tk.X, pady=(8, 0))
        ttk.Button(btns, text="Save  (Ctrl+S)", command=self.save).pack(side=tk.LEFT)
        ttk.Button(btns, text="Save & next", command=lambda: (self.save(), self.step(+1))).pack(side=tk.LEFT, padx=6)
        ttk.Button(btns, text="Summarize", command=self.summarize).pack(side=tk.RIGHT)
        self.status_var = tk.StringVar()
        ttk.Label(f, textvariable=self.status_var, foreground="#444", wraplength=360).pack(anchor="w", pady=(8, 0))
        ttk.Label(f, text=f"sheet: {R.SHEET.relative_to(R.SCOPE_DIR.parent)}", foreground="#888",
                  wraplength=360).pack(anchor="w", side=tk.BOTTOM)

    def refresh_combo(self):
        sh = self.data.sheet()
        labels = []
        for d in self.ids:
            t = self.data.tab.loc[d]
            n = int(t["n_cells"])
            role = sh.loc[d, "role"] if d in sh.index and "role" in sh else ("control" if n == 0 else "changed")
            done = d in sh.index and any(isinstance(sh.loc[d, c], str) and sh.loc[d, c].strip()
                                         for c in R.VERDICT_COLS if c in sh)
            picked = d in sh.index and "lost_behind_road_m" in sh \
                and isinstance(sh.loc[d, "lost_behind_road_m"], str) and sh.loc[d, "lost_behind_road_m"].strip()
            labels.append(f"GIS {d:>2}   {n:+d} rows   {role:8s} {'done' if done else '    '}{' measured' if picked else ''}")
        self.combo["values"] = labels
        if labels:
            self.combo.current(min(self.i, len(labels) - 1))

    def status(self, msg: str):
        self.status_var.set(f"{datetime.now():%H:%M:%S}  {msg}")
        self.root.update_idletasks()

    def _key(self, fn):
        """Keys act unless the cursor is in a text box or drop-down."""
        if isinstance(self.root.focus_get(), (tk.Text, ttk.Combobox)):
            return
        fn()

    # ---- navigation ---------------------------------------------------------
    def step(self, k: int):
        j = self.i + k
        if 0 <= j < len(self.ids):
            self.show(j)

    def show(self, j: int):
        self.i = j
        d = self.ids[j]
        self.combo.current(j)
        self.cancel_pick(quiet=True)
        self.status(f"loading GIS {d} ...")
        try:
            rec = self.data.load(d)
        except Exception as e:                      # a bad frame must not kill the session
            self.status(f"GIS {d}: {e}")
            return
        self.fill_form(d)                            # picks first: draw() draws them
        self.draw(d, rec)
        self.fill_info(d, rec)
        self.status(f"GIS {d} shown")
        nxt = self.ids[j + 1] if j + 1 < len(self.ids) else None
        if nxt is not None:
            threading.Thread(target=self._prefetch, args=(nxt,), daemon=True).start()
        if self.smoke:
            self.root.after(1500, self._smoke_done)

    def _prefetch(self, d: int):
        try:
            self.data.load(d)
        except Exception:
            pass

    # ---- drawing ------------------------------------------------------------
    def draw(self, d: int, rec: dict):
        self.fig.clear()
        self.artists = {k: [] for k, _, _ in OVERLAYS}
        self.blink_images = {}
        pr, pl, t = rec["pr"], rec["pl"], rec["t"]
        n = int(t["n_cells"])
        x_lo, x_hi, y_lo, y_hi = rec["win"]
        x_med = float(pr["interior_x"].median())
        years = list(self.data.years)
        blink = self.blink.get()
        n_photo = 1 if blink else len(years)
        n_pan = n_photo + (1 if self.show_lidar.get() else 0)
        gs = self.fig.add_gridspec(2, n_pan, height_ratios=[3.2, 1.0], hspace=0.18, wspace=0.04,
                                   left=0.04, right=0.99, top=0.95, bottom=0.08)
        axes = []
        for k in range(n_pan):
            ax = self.fig.add_subplot(gs[0, k], sharex=axes[0] if axes else None, sharey=axes[0] if axes else None)
            axes.append(ax)
        axs = self.fig.add_subplot(gs[1, :])
        gd = self.data.gdf[self.data.gdf["domain_id"].astype(int) == d]
        box = gd.geometry.iloc[0]
        drawn = off.clip_for_drawing(self.data.lines, box.buffer(30.0))
        A = self.artists

        def snapshot(ax):
            return list(ax.lines) + list(ax.collections) + list(ax.patches)

        def new_artists(ax, before):
            return [a for a in snapshot(ax) if a not in before]

        def overlays(ax, on_photo):
            before = snapshot(ax)
            gd.boundary.plot(ax=ax, color="0.3", linewidth=0.8, zorder=4)
            A["box"] += new_artists(ax, before)
            for yr, key in ((1997, "line97"), (1984, "line84")):     # 1984 on top, as draw_lines
                before = snapshot(ax)
                if len(drawn[yr]):
                    z = 6 + off.LINE_ORDER.index(yr) * 2
                    drawn[yr].plot(ax=ax, linestyle="-", alpha=0.9, zorder=z, **off.LINE_CASING[yr])
                    drawn[yr].plot(ax=ax, zorder=z + 1, **off.LINE_STYLE[yr])
                A[key] += new_artists(ax, before)
            before = snapshot(ax)
            off.m.draw_roads(ax, self.data.roads, scale=0.9)
            A["roads"] += new_artists(ax, before)
            A["row0"] += ax.plot(pr["interior_x"], pr["interior_y"], color="white" if on_photo else R.INK,
                                 lw=1.2, ls=(0, (4, 2)), zorder=6)
            if pl["pav"] is not None:
                A["roadrows"].append(ax.add_patch(plt.Polygon(
                    R._band_poly(pr, *pl["pav"]), closed=True, facecolor="none",
                    edgecolor="white" if on_photo else R.C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=7)))
            if pl["road"] is not None:
                for rr in pl["road"]:
                    A["roadrows"].append(ax.add_patch(plt.Polygon(
                        R._band_poly(pr, rr, rr), closed=True, facecolor=R.C_ROAD, edgecolor="none",
                        alpha=0.35, zorder=7)))
            if pl["v3"] is not None:
                A["v3"].append(ax.add_patch(plt.Polygon(
                    R._band_poly(pr, *pl["v3"]), closed=True, facecolor="none",
                    edgecolor=R.C_ADD if n > 0 else R.C_REM, lw=2.0, zorder=8,
                    hatch=None if on_photo else "////")))
            if pl["sea"] is not None:
                A["sea"].append(ax.add_patch(plt.Polygon(
                    R._band_poly(pr, *pl["sea"]), closed=True, facecolor="none", edgecolor=R.C_SEA_ALT,
                    lw=1.8, ls=(0, (5, 2)), zorder=8)))
            ax.set_xlim(x_lo, x_hi)
            ax.set_ylim(y_lo, y_hi)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])

        def photo(ax, y, visible=True):
            shown = rec["photos"].get(y)
            if shown is not None:
                return ax.imshow(shown, extent=(x_lo, x_hi, y_lo, y_hi), origin="upper", zorder=1,
                                 interpolation="bilinear", visible=visible)
            ax.set_facecolor("0.93")
            return ax.text(0.5, 0.5, f"no {y} photograph here", transform=ax.transAxes,
                           ha="center", va="center", visible=visible)

        if blink:
            ax = axes[0]
            for y in years:
                self.blink_images[y] = photo(ax, y, visible=False)
            if self.blink_year not in years:
                self.blink_year = years[0]
            overlays(ax, on_photo=True)
            off._scalebar(ax, length_m=100.0)
            off._north_arrow(ax, x=0.94, y=0.12)
        else:
            for k, y in enumerate(years):
                ax = axes[k]
                photo(ax, y)
                overlays(ax, on_photo=True)
                im = self.data.imagery[y]
                ax.set_title(f"({chr(97 + k)}) {im.date}   cover {rec['cover'].get(y, 0):.0%}", loc="left")
                if k == 0:
                    off._scalebar(ax, length_m=100.0)
                    off._north_arrow(ax, x=0.94, y=0.12)
        if self.show_lidar.get():
            ax = axes[-1]
            arr, extent = self.data.lidar(d)
            if arr is not None:
                off._hillshade(ax, arr, extent, res=1.0)
            overlays(ax, on_photo=False)
            ax.set_title(f"({chr(97 + n_photo)}) 1 m lidar (1996 foredune, 2009 backdune)", loc="left")

        # the strip
        for y in years:
            if y not in rec["strips"]:
                continue
            cells, S = rec["strips"][y]
            x = x_med - cells * R.CELL_M
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                med = np.nanmedian(S, axis=0)
                p25, p75 = np.nanpercentile(S, 25, axis=0), np.nanpercentile(S, 75, axis=0)
            col = YEAR_COLOUR.get(y, "0.3")
            axs.fill_between(x, p25, p75, color=col, alpha=0.15, lw=0)
            axs.plot(x, med, color=col, lw=1.6, label=str(self.data.imagery[y].date))
        x_l84 = float((pr["interior_x"] - pr["r_line84"] * R.CELL_M).median())
        x_l97 = float((pr["interior_x"] - pr["r_line97"] * R.CELL_M).median())
        A["stripmarks"] += [axs.axvline(x_l84, color=YEAR_COLOUR[1984], lw=1.4, zorder=3),
                            axs.axvline(x_l97, color=YEAR_COLOUR[1997], lw=1.4, zorder=3),
                            axs.axvline(x_med, color=R.INK, lw=1.0, ls=(0, (4, 2)), zorder=3)]

        def span(rr, **kw):
            if rr is not None:
                A["stripmarks"].append(axs.axvspan(x_med - rr[1] * R.CELL_M - R.CELL_M / 2,
                                                   x_med - rr[0] * R.CELL_M + R.CELL_M / 2, **kw))
        span(pl["pav"], facecolor="none", edgecolor=R.C_ROAD_OLD, lw=1.0, ls=(0, (2, 2)), zorder=2)
        if pl["road"] is not None:
            for rr in pl["road"]:
                span((rr, rr), facecolor=R.C_ROAD, alpha=0.3, lw=0, zorder=1)
        span(pl["v3"], facecolor=R.C_ADD if n > 0 else R.C_REM, alpha=0.18, lw=0, zorder=1)
        span(pl["sea"], facecolor="none", edgecolor=R.C_SEA_ALT, lw=1.6, ls=(0, (5, 2)), zorder=2)
        axs.set_xlim(x_lo, x_hi)
        axs.set_ylim(-0.05, 1.05)
        axs.set_ylabel("brightness")
        axs.set_xlabel("easting (m)")
        sec = axs.secondary_xaxis("top", functions=(lambda xx: (x_med - xx), lambda mm: x_med - mm))
        sec.set_xlabel("m landward of interior row 0", fontsize=8)
        axs.grid(True, axis="x", alpha=0.4)
        axs.legend(loc="upper left", fontsize=8)
        self.axes, self.axs = axes[:n_photo], axs
        self.all_axes = axes
        if blink:
            self._show_blink_year()
        self.draw_picks()
        for key in self.artists:
            self._set_visible(key, self.overlay_vars[key].get())
        self.canvas.draw_idle()

    def draw_picks(self):
        """Your picks, on every map panel and on the strip: colour by year, marker by feature."""
        for a in self.artists.get("picks", []):
            try:
                a.remove()
            except Exception:
                pass
        self.artists["picks"] = []
        for (y, kind), m in self.picks.items():
            k = f"{kind}{str(y)[2:]}"
            x, yv = _num(m.get(f"{k}_x")), _num(m.get(f"{k}_y"))
            if not (np.isfinite(x) and np.isfinite(yv)):
                continue
            col = YEAR_COLOUR.get(y, "0.2")
            up = y == PICK_YEARS[0]                    # 1984 labelled above, 1997 below
            dy = 22 + 13 * PICK_KINDS.index(kind)      # staggered by feature, so close picks still read
            yv_lab = yv + dy if up else yv - dy
            for ax in self.all_axes:
                self.artists["picks"] += ax.plot([x, x], [yv - 20, yv + 20], color="white", lw=4.0, zorder=11,
                                                 solid_capstyle="butt")
                self.artists["picks"] += ax.plot([x, x], [yv - 20, yv + 20], color=col, lw=2.0, zorder=12,
                                                 solid_capstyle="butt")
                self.artists["picks"] += ax.plot([x], [yv], marker=PICK_MARKER[kind], ms=8, mfc="none", mec=col,
                                                 mew=2.0, zorder=12, ls="none")
                self.artists["picks"].append(ax.annotate(f"{kind} {y}", (x, yv_lab), ha="center",
                                                         va="bottom" if up else "top",
                                                         fontsize=8, color=col, fontweight="bold", zorder=12,
                                                         bbox=dict(facecolor="white", alpha=0.8, edgecolor="none",
                                                                   boxstyle="square,pad=0.15")))
            if self.axs is not None:
                self.artists["picks"].append(self.axs.axvline(x, color=col, lw=1.6, ls=(0, (1, 1.5)), zorder=4))
        self._set_visible("picks", self.overlay_vars["picks"].get())

    def _set_visible(self, key: str, on: bool):
        for a in self.artists.get(key, []):
            try:
                a.set_visible(on)
            except Exception:
                pass

    def toggle(self, key: str):
        self._set_visible(key, self.overlay_vars[key].get())
        self.canvas.draw_idle()

    def on_scroll(self, event):
        """Wheel zoom about the pointer on the map panels (all move together)."""
        ax = event.inaxes
        if ax is None or ax not in self.all_axes or event.xdata is None:
            return
        f = 0.8 if event.button == "up" else 1.25
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_xlim(event.xdata - (event.xdata - x0) * f, event.xdata + (x1 - event.xdata) * f)
        ax.set_ylim(event.ydata - (event.ydata - y0) * f, event.ydata + (y1 - event.ydata) * f)
        self.canvas.draw_idle()

    # ---- blink ----------------------------------------------------------------
    def _show_blink_year(self):
        for y, im in self.blink_images.items():
            im.set_visible(y == self.blink_year)
        if self.blink.get() and self.axes and self.blink_year is not None:
            imy = self.data.imagery[self.blink_year]
            rec = self.data.load(self.ids[self.i])
            self.axes[0].set_title(f"{imy.date}   cover {rec['cover'].get(self.blink_year, 0):.0%}"
                                   f"    (Space flips)", loc="left", color=YEAR_COLOUR.get(self.blink_year, "0.2"))
        self.canvas.draw_idle()

    def flip(self):
        if not self.blink.get() or not self.blink_images:
            return
        years = list(self.blink_images)
        self.blink_year = years[(years.index(self.blink_year) + 1) % len(years)]
        self._show_blink_year()

    def _auto_blink(self):
        if self._blink_job is not None:
            self.root.after_cancel(self._blink_job)
            self._blink_job = None
        if self.auto_blink.get():
            if not self.blink.get():
                self.blink.set(True)
                self.show(self.i)
            self._tick()

    def _tick(self):
        if not self.auto_blink.get():
            return
        self.flip()
        self._blink_job = self.root.after(BLINK_MS, self._tick)

    # ---- picks ------------------------------------------------------------------
    def start_pick(self, year: int, kind: str):
        if self.toolbar.mode:                # leave pan/zoom so the click reaches us
            if str(self.toolbar.mode).lower().startswith("pan"):
                self.toolbar.pan()
            else:
                self.toolbar.zoom()
        self.pick_target = (year, kind)
        if self.blink.get():
            self.blink_year = year
            self._show_blink_year()
        self.canvas.get_tk_widget().configure(cursor="crosshair")
        self.status(f"click the {PICK_LABEL[kind]} on the {year} photograph (Esc cancels)")

    def cancel_pick(self, quiet: bool = False):
        self.pick_target = None
        try:
            self.canvas.get_tk_widget().configure(cursor="")
        except tk.TclError:
            pass
        if not quiet:
            self.status("pick cancelled")

    def on_click(self, event):
        if self.pick_target is None or event.button != 1 or event.inaxes is None:
            return
        if event.inaxes not in self.all_axes or event.xdata is None:
            return
        if self.toolbar.mode:                # pan/zoom drag, not a pick
            return
        rec = self.data.load(self.ids[self.i])
        year, kind = self.pick_target
        self.picks[(year, kind)] = measure_pick(rec, year, kind, float(event.xdata), float(event.ydata))
        self.cancel_pick(quiet=True)
        self.draw_picks()
        self.canvas.draw_idle()
        self.update_measure()
        prof = self.picks[(year, kind)][f"{kind}{str(year)[2:]}_profile"]
        self.status(f"{kind} {year} picked on profile {prof} - Save to keep it")

    def clear_picks(self):
        self.picks = {}
        self.draw_picks()
        self.canvas.draw_idle()
        self.update_measure()
        self.status("picks cleared (Save to write the blanks)")

    def measured(self) -> dict:
        m = {}
        for v in self.picks.values():
            m.update(v)
        m.update(derived(m, int(self.data.tab.loc[self.ids[self.i], "n_cells"])))
        return m

    def update_measure(self):
        m = self.measured()
        n = int(self.data.tab.loc[self.ids[self.i], "n_cells"])
        lines = []
        for y in PICK_YEARS:
            yy = str(y)[2:]
            have = [k for k in PICK_KINDS if f"{k}{yy}_x" in m]
            if not have:
                continue
            parts = [f"{k} {m[f'{k}{yy}_from_row0_m']:+.0f}" for k in have]
            bands = [f"{lab} {m[f'{b}{yy}_m']:.0f}" for b, lab in
                     (("dune_band", "dune"), ("back_to_road", "dune->road"), ("toe_to_road", "toe->road"))
                     if m[f"{b}{yy}_m"] != ""]
            lines.append(f"{y}: " + ", ".join(parts) + " m from row 0"
                         + ("; bands " + ", ".join(bands) + " m" if bands else ""))
        shifts = [f"{k} {m[f'{k}_shift_m']:+.0f}" for k in PICK_KINDS if m[f"{k}_shift_m"] != ""]
        if shifts:
            lines.append("shift 1984->1997 (+ = 1984 seaward): " + ", ".join(shifts) + " m"
                         + (f"; toe = {m['toe_shift_cells']:+d} cells vs N {n:+d}" if m["toe_shift_cells"] != "" else ""))
        if m["lost_behind_road_m"] != "":
            lines.append(f"1984 width, given N x 10 = {n * R.CELL_M:+.0f} m: dune band {m['lost_dune_band_m']:+.0f}, "
                         f"dune->road {m['lost_back_to_road_m']:+.0f}, behind road {m['lost_behind_road_m']:+.0f} m")
            lines.append(f"  -> bands suggest: {m['band_suggests']}   (derived; the verdict is yours)")
        self.measure_var.set("\n".join(lines) if lines else "no picks yet")

    # ---- the form -------------------------------------------------------------
    def fill_info(self, d: int, rec: dict):
        t = rec["t"]
        n = int(t["n_cells"])
        where = {"road": "behind NC-12 as placed" if n > 0 else "in front of NC-12 today",
                 "crest": "behind the crest row"}.get(str(t.get("insert_anchor", "")), "-")
        has_road = np.isfinite(_num(t.get("setback_v2_m", np.nan)))
        lines = [f"GIS {d}   {t['action']}   N = {n:+d} rows",
                 f"shift 1984->1997: {t['shift_m_median']:+.0f} m  (p10 {t['shift_m_p10']:+.0f}, p90 {t['shift_m_p90']:+.0f})",
                 f"v3 places them: {where}",
                 f"   {t.get('rows_behind_road', '') or '-'}",
                 (f"setback: {t['setback_v2_m']:.0f} m today -> {t['setback_new_m']:.0f} m for 1984"
                  if has_road else "no model road in this domain"),
                 f"flags: {t['flags'] if isinstance(t['flags'], str) else '-'}",
                 "photo cover: " + ", ".join(f"{y} {c:.0%}" for y, c in rec["cover"].items()),
                 "",
                 "+ = the 1984 line seaward of the 1997 line = rows added.",
                 "Red/blue outline: v3's rows. Purple dashed: the same rows at",
                 "the dune. Dotted white: today's pavement rows. Dark: as placed.",
                 "Picks: o toe, square back of dune, diamond road; red 1984, blue 1997."]
        self.info.configure(state="normal")
        self.info.delete("1.0", tk.END)
        self.info.insert("1.0", "\n".join(lines))
        self.info.configure(state="disabled")

    def fill_form(self, d: int):
        sh = self.data.sheet()
        for c, cb in self.fields.items():
            v = sh.loc[d, c] if d in sh.index and c in sh else ""
            cb.set(v if isinstance(v, str) else "")
        self.notes.delete("1.0", tk.END)
        v = sh.loc[d, "notes"] if d in sh.index and "notes" in sh else ""
        self.notes.insert("1.0", v if isinstance(v, str) else "")
        self.picks = {}
        for y in PICK_YEARS:
            for kind in PICK_KINDS:
                k = f"{kind}{str(y)[2:]}"
                cols = [f"{k}_{f}" for f in ("x", "y", "profile", "from_row0_m")]
                if d in sh.index and all(c in sh for c in cols) and np.isfinite(_num(sh.loc[d, f"{k}_x"])):
                    self.picks[(y, kind)] = {c: (int(_num(sh.loc[d, c])) if c.endswith("_profile")
                                                 else _num(sh.loc[d, c])) for c in cols}
        self.update_measure()

    def save(self):
        d = self.ids[self.i]
        if not R.SHEET.is_file():
            self.status("no sheet on disk - run HAT_imagery_review_1984.py first")
            return
        df = pd.read_csv(R.SHEET, dtype={c: str for c in R.VERDICT_COLS + R.MEASURED_COLS}).set_index("domain")
        if d not in df.index:
            self.status(f"GIS {d} is not in the sheet - run the batch script on it first")
            return
        for c in R.VERDICT_COLS + R.MEASURED_COLS:
            if c not in df:
                df[c] = ""
        for c, cb in self.fields.items():
            df.loc[d, c] = cb.get()
        df.loc[d, "notes"] = self.notes.get("1.0", tk.END).strip().replace("\n", " / ")
        m = self.measured()
        for c in R.MEASURED_COLS:
            if c in ("reviewed_by", "reviewed_at"):
                continue
            v = m.get(c, "")
            df.loc[d, c] = "" if (v == "" or (isinstance(v, float) and not np.isfinite(v))) else str(v)
        df.loc[d, "reviewed_by"] = getpass.getuser()
        df.loc[d, "reviewed_at"] = f"{datetime.now():%Y-%m-%d %H:%M}"
        df.to_csv(R.SHEET)
        self.refresh_combo()
        self.status(f"GIS {d} saved to {R.SHEET.name}")

    def summarize(self):
        try:
            import HAT_imagery_review_summary as S
            fig, rep = S.run()
            self.status(f"summary written: {fig.name}, {rep.name}")
        except Exception as e:
            self.status(f"summary failed: {e}")

    # ---- lifecycle ------------------------------------------------------------
    def _smoke_done(self):
        out = CACHE_DIR.parent / "imagery_review_gui_smoke.png"       # outside the repo
        self.fig.savefig(out, dpi=100)
        print(f"smoke: drew GIS {self.ids[self.i]}, saved {out}")
        self.on_close()

    def on_close(self):
        if self._blink_job is not None:
            self.root.after_cancel(self._blink_job)
        self.root.quit()
        self.root.destroy()

    def run(self):
        self.root.mainloop()


def main() -> None:
    ap = argparse.ArgumentParser(description="the imagery review as a window")
    ap.add_argument("--domains", default="", help="comma-separated GIS ids (default: the sheet's domains)")
    ap.add_argument("--years", default=",".join(str(y) for y in R.DEFAULT_YEARS))
    ap.add_argument("--smoke", action="store_true", help="open, draw the first domain, screenshot, close")
    a = ap.parse_args()
    years = [int(y) for y in a.years.split(",") if y.strip()]
    data = Data(years)
    if a.domains:
        ids = [int(x) for x in a.domains.split(",")]
    elif R.SHEET.is_file():
        ids = sorted(int(d) for d in pd.read_csv(R.SHEET)["domain"])
    else:
        ids = sorted(int(d) for d in data.tab.index if int(data.tab.loc[d, "n_cells"]) != 0)
    App(data, ids, smoke=a.smoke).run()


if __name__ == "__main__":
    main()
