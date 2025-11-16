# === cascade_npz_explorer.py (organized + registry + helpers) ===============
# Purpose:
#   - Load a CASCADE .npz run (which contains a serialized Cascade object)
#   - Introspect arrays/DataFrames within the object (and nested modules)
#   - Save:
#       * <RUN_TAG>__metadata.json
#       * <RUN_TAG>__cascade_contents.csv
#       * <RUN_TAG>__likely_outputs.csv
#       * <RUN_TAG>__explorer.log
#   - Append/maintain a cross-run registry CSV at OUT_ROOT/_run_registry.csv
#
# Usage:
#   1) Set NPZ_DIR and RUN_NAME (no extension).
#   2) Fill the RUN_TAG parts below (ISL, YEARS, SCENARIO, MGMT, NOTES, VERSION)
#   3) Run this file.
#
# Notes:
#   - RUN_TAG format: <ISL>_<YEARS>_<SCENARIO>[_<MGMT>][_NOTES][_vNN]
#   - YEARS should be "YYYY-YYYY" (Windows-safe)
#   - OUT_DIR: <OUT_ROOT>/<RUN_TAG>/
# ============================================================================

import os
import re
import json
import datetime as dt
from collections import deque

import numpy as np
import pandas as pd

# ---------------- USER INPUTS -----------------------------------------------
# Where your .npz lives (the file produced by your CASCADE run)
NPZ_DIR  = r"C:\Users\hanna\PycharmProjects\CASCADE\output\Hindcast_1978_1997_buffers"
RUN_NAME = "Hindcast_1978_1997_buffers"  # no file extension

# RUN_TAG parts (these build the folder/filenames for analysis artifacts)
ISL      = "HAT"            # e.g., HAT (Hatteras), OCR (Ocracoke)
YEARS    = "1978-1997"      # Windows-safe hyphen
SCENARIO = "hindcast"       # e.g., hindcast, adaptmgmt, bnourish, slrX
MGMT     = "nooffset"               # optional descriptor (e.g., roadsetback50m)
NOTES    = ""               # optional short token (e.g., S2, highwaves, dm_units)
VERSION  = ""               # optional version like "v01"
# Where to write analysis outputs
OUT_ROOT = r"C:\Users\hanna\PycharmProjects\CASCADE\analysis"
# ---------------------------------------------------------------------------


# === Helpers for consistent naming/paths ====================================
def _clean_token(x: str) -> str:
    if not x:
        return ""
    # Allow letters, numbers, hyphen; convert spaces to nothing; underscore later as separators
    x = x.strip().replace(" ", "")
    return re.sub(r"[^A-Za-z0-9\-]", "", x)

def build_run_tag(isl, years, scenario, mgmt="", notes="", version=""):
    parts = [_clean_token(isl), _clean_token(years), _clean_token(scenario)]
    if mgmt:   parts.append(_clean_token(mgmt))
    if notes:  parts.append(_clean_token(notes))
    if version: parts.append(_clean_token(version))
    return "_".join([p for p in parts if p])

RUN_TAG = build_run_tag(ISL, YEARS, SCENARIO, MGMT, NOTES, VERSION)

# Set up output folder
os.makedirs(OUT_ROOT, exist_ok=True)
OUT_DIR = os.path.join(OUT_ROOT, RUN_TAG)
os.makedirs(OUT_DIR, exist_ok=True)

# Optional: if you ever want timestamped subfolders
# OUT_DIR = os.path.join(OUT_DIR, dt.datetime.now().strftime("%Y%m%d-%H%M%S"))
# os.makedirs(OUT_DIR, exist_ok=True)


# === Logging (tee to console + file) ========================================
log_path = os.path.join(OUT_DIR, f"{RUN_TAG}__explorer.log")
def logprint(*args, **kwargs):
    msg = " ".join(str(a) for a in args)
    print(msg, **kwargs)
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(msg + ("\n" if not msg.endswith("\n") else ""))


# === Load the .npz and Cascade object =======================================
npz_path = os.path.join(NPZ_DIR, RUN_NAME + ".npz")
if not os.path.exists(npz_path):
    raise FileNotFoundError(f"Could not find npz at: {npz_path}")

npz = np.load(npz_path, allow_pickle=True)
if "cascade" not in npz:
    raise RuntimeError("No 'cascade' key found. Keys present: " + ", ".join(npz.keys()))
cascade = npz["cascade"][0]  # CASCADE object


# === Metadata capture ========================================================
def safe_get(obj, name, default=None):
    try:
        return getattr(obj, name)
    except Exception:
        return default

meta_keys = [
    "_nt","_ny","_number_of_communities","_sea_level_rise_rate",
    "_wave_height","_wave_period","_wave_asymmetry","_wave_angle_high_fraction",
    "_enable_shoreline_offset","_road_setback","_road_width"
]

logprint("\n=== CASCADE metadata (best-effort) ===")
metadata = {}
for k in meta_keys:
    v = safe_get(cascade, k, None)
    if isinstance(v, np.ndarray):
        metadata[k] = {"shape": tuple(v.shape), "dtype": str(v.dtype)}
        logprint(f"{k}: shape={v.shape} dtype={v.dtype}")
    else:
        metadata[k] = v
        logprint(f"{k}: {v}")

modules = [name for name in [
    "brie","barrier3d","roadways","roadway_management_module",
    "beach_nourishment_module","chom","module_lists"
] if hasattr(cascade, name)]
metadata["modules_detected"] = modules
logprint("\nModules detected:", modules)

# Save machine-readable metadata
with open(os.path.join(OUT_DIR, f"{RUN_TAG}__metadata.json"), "w", encoding="utf-8") as f:
    json.dump({
        "npz_path": npz_path,
        "run_tag": RUN_TAG,
        "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
        **metadata
    }, f, indent=2)


# === Object graph traversal to list arrays/frames ===========================
def is_simple(x):
    return isinstance(x, (int, float, str, bool, type(None), np.number, np.bool_))

def iter_public_attrs(obj):
    for name in dir(obj):
        if name.startswith("__"):
            continue
        yield name  # include private like _foo — CASCADE uses many

def summarize_ndarray(arr: np.ndarray):
    return {
        "type": "ndarray",
        "shape": tuple(arr.shape),
        "dtype": str(arr.dtype),
        "nbytes": int(arr.nbytes),
    }

def summarize_df(df: pd.DataFrame):
    return {
        "type": "DataFrame",
        "shape": (df.shape[0], df.shape[1]),
        "columns": list(map(str, df.columns[:12])),
    }

seen_ids = set()
queue = deque([((), cascade)])
summaries = []
MAX_NODES = 5000
MAX_DEPTH = 5

while queue and len(seen_ids) < MAX_NODES:
    path, obj = queue.popleft()
    oid = id(obj)
    if oid in seen_ids:
        continue
    seen_ids.add(oid)

    # Collect arrays / dataframes
    try:
        if isinstance(obj, np.ndarray):
            summaries.append({"path": ".".join(path), "summary": summarize_ndarray(obj)})
            continue
        if isinstance(obj, pd.DataFrame):
            summaries.append({"path": ".".join(path), "summary": summarize_df(obj)})
            continue
    except Exception:
        pass

    # Recurse into dicts / lists / tuples
    if isinstance(obj, dict):
        for k, v in list(obj.items())[:1000]:
            if len(path) < MAX_DEPTH:
                queue.append((path + (str(k),), v))
        continue
    if isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj[:1000]):
            if len(path) < MAX_DEPTH:
                queue.append((path + (f"[{i}]",), v))
        continue

    # Explore attributes for custom classes/modules
    if is_simple(obj):
        continue
    try:
        for name in iter_public_attrs(obj):
            if len(path) >= MAX_DEPTH:
                break
            try:
                val = getattr(obj, name)
            except Exception:
                continue
            queue.append((path + (name,), val))
    except Exception:
        continue


# === Write summaries =========================================================
def human_row(row):
    p = row["path"] or "(root)"
    s = row["summary"]
    if s["type"] == "ndarray":
        return f"{p:80s}  ndarray  shape={s['shape']}  dtype={s['dtype']}  ~{s['nbytes']:,} bytes"
    else:
        return f"{p:80s}  DataFrame  shape={s['shape']}  cols(sample)={s['columns']}"

logprint("\n=== Found arrays / dataframes in CASCADE object (truncated) ===")
for row in summaries[:2000]:
    logprint(human_row(row))

# Full contents CSV
contents_csv = os.path.join(OUT_DIR, f"{RUN_TAG}__cascade_contents.csv")
rows_out = []
for row in summaries:
    s = row["summary"]
    kind = s["type"]
    shape = s["shape"]
    extra = s.get("dtype", None) if kind == "ndarray" else ",".join(map(str, s.get("columns", [])))
    rows_out.append([row["path"], kind, shape, extra])
pd.DataFrame(rows_out, columns=["path","kind","shape","extra"]).to_csv(contents_csv, index=False)
logprint(f"\nSaved summary to {contents_csv}")

# Likely outputs CSV (filters on common names)
patterns = re.compile(r"(shoreline|dune|crest|elev|ele|height|width|overwash|storm|x_|y_|position|road|nourish)", re.I)
hits = [r for r in summaries if patterns.search(r["path"])]
hits_csv = os.path.join(OUT_DIR, f"{RUN_TAG}__likely_outputs.csv")
pd.DataFrame(
    [[r["path"], r["summary"]["type"], r["summary"].get("shape")] for r in hits],
    columns=["path","kind","shape"]
).to_csv(hits_csv, index=False)
logprint(f"Saved likely outputs list to {hits_csv}")


# === Master registry (append one row per analysis) ==========================
registry_csv = os.path.join(OUT_ROOT, "_run_registry.csv")
row = {
    "run_tag": RUN_TAG,
    "npz_path": npz_path,
    "out_dir": OUT_DIR,
    "analyzed_at": dt.datetime.now().isoformat(timespec="seconds"),
    "nt": metadata.get("_nt"),
    "ny": metadata.get("_ny"),
    "modules": ";".join(metadata.get("modules_detected", [])),
}
pd.DataFrame([row]).to_csv(
    registry_csv, mode="a", header=not os.path.exists(registry_csv), index=False
)
logprint(f"Appended registry row -> {registry_csv}")


# === Optional helpers for later plotting/inspection =========================
def get_by_path(root, path_str):
    """
    Access attribute/index path like: 'brie.shoreline_position[0]'
    Returns the referenced object (e.g., np.ndarray), or raises if not found.
    """
    tokens = re.findall(r"[A-Za-z_]\w*|\[\d+\]", path_str)
    cur = root
    for t in tokens:
        if t.startswith("["):
            cur = cur[int(t[1:-1])]
        else:
            cur = getattr(cur, t)
    return cur

def quick_plot(arr, title="", save_as=None):
    """
    Quick 1D/2D plot. If save_as is provided (filename stem), saves PNG in OUT_DIR.
    """
    import matplotlib.pyplot as plt
    if not isinstance(arr, np.ndarray):
        arr = np.asarray(arr)
    if arr.ndim == 1:
        plt.figure(figsize=(9, 4))
        plt.plot(arr)
        plt.xlabel("Index / time")
        plt.ylabel("Value")
        plt.title(title or "1D series")
        plt.tight_layout()
    elif arr.ndim == 2:
        plt.figure(figsize=(9, 4))
        plt.imshow(arr, aspect="auto", origin="lower")
        plt.colorbar(label="Value")
        plt.xlabel("Alongshore / index")
        plt.ylabel("Cross-shore / time (depends on field)")
        plt.title(title or "2D field")
        plt.tight_layout()
    else:
        logprint("Array has ndim > 2; please slice before plotting.")
        return

    if save_as:
        out_png = os.path.join(OUT_DIR, f"{RUN_TAG}__{save_as}.png")
        plt.savefig(out_png, dpi=300)
        logprint(f"Saved figure: {out_png}")
    plt.show()


logprint("\nOutput directory:", OUT_DIR)
logprint("Done.")
# ============================================================================

