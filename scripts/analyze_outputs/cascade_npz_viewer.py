# =============================================================================
# CASCADE NPZ Introspector (Explorer v1)
# - Lists all numpy arrays inside the saved 'Cascade' object
# - Saves a CSV catalog of array paths, shapes, and dtypes
# - Provides a helper to extract a time×domain matrix once you pick an attribute
# =============================================================================
import os
import csv
import numpy as np
from types import SimpleNamespace

# --- UPDATE THESE ---
NPZ_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\output\HAT_1978_1997_Natural_State\HAT_1978_1997_Natural_State.npz"
CATALOG_CSV = os.path.join(os.path.dirname(NPZ_PATH), "cascade_array_catalog.csv")

# -------------------- Load the Cascade object --------------------
data = np.load(NPZ_PATH, allow_pickle=True)
assert "cascade" in data.files, f"No 'cascade' key found. Keys: {data.files}"
cascade_obj = data["cascade"][0]
print("Loaded Cascade object:", type(cascade_obj))

# -------------------- Introspection helpers --------------------
import inspect
from collections import deque

def is_numpy_array(x):
    return isinstance(x, np.ndarray)

def safe_dir(obj):
    try:
        return dir(obj)
    except Exception:
        return []

def safe_getattr(obj, name):
    try:
        return getattr(obj, name)
    except Exception:
        return SimpleNamespace(__unreadable__=True)

def walk_numpy_arrays(root, max_depth=4, max_items_per_node=200):
    """
    Breadth-first walk that records every numpy array found inside `root`.
    Produces (path, shape, dtype) rows.
    """
    visited_ids = set()
    rows = []

    def enqueue(obj, path, depth):
        oid = id(obj)
        if oid in visited_ids:
            return
        visited_ids.add(oid)
        q.append((obj, path, depth))

    q = deque()
    enqueue(root, "cascade", 0)

    while q:
        obj, path, depth = q.popleft()
        if depth > max_depth:
            continue

        # record arrays
        if is_numpy_array(obj):
            rows.append((path, obj.shape, str(obj.dtype)))
            continue

        # containers
        if isinstance(obj, (list, tuple)):
            for idx, item in enumerate(obj[:max_items_per_node]):
                enqueue(item, f"{path}[{idx}]", depth + 1)
            continue

        if isinstance(obj, dict):
            for k, v in list(obj.items())[:max_items_per_node]:
                enqueue(v, f"{path}['{k}']", depth + 1)
            continue

        # objects: iterate attributes
        # skip builtins, callables, modules, etc.
        for name in safe_dir(obj):
            if name.startswith("__") and name.endswith("__"):
                continue
            # keep attribute list manageable
            if len(name) > 80:
                continue
            val = safe_getattr(obj, name)
            if isinstance(val, (int, float, str, bool, type(None))):
                continue
            if inspect.ismodule(val) or inspect.isfunction(val) or inspect.ismethod(val):
                continue
            enqueue(val, f"{path}.{name}", depth + 1)

    return rows

# -------------------- Run the walk & print/save catalog --------------------
rows = walk_numpy_arrays(cascade_obj, max_depth=5)
rows_sorted = sorted(rows, key=lambda r: (r[0], r[1]))

print("\nFound numpy arrays inside the Cascade object:")
for path, shape, dtype in rows_sorted[:200]:
    print(f"- {path:80s}  shape={shape}  dtype={dtype}")
if len(rows_sorted) > 200:
    print(f"... ({len(rows_sorted)-200} more not shown)")

# Save to CSV
with open(CATALOG_CSV, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["path", "shape", "dtype"])
    for path, shape, dtype in rows_sorted:
        w.writerow([path, str(shape), dtype])
print(f"\nCatalog saved: {CATALOG_CSV}")

# -------------------- Optional: extract time×domain matrix --------------------
def extract_time_by_domain(cascade_obj, per_domain_attr_path, domain_count=None):
    """
    Extract a 2D array of shape (time, domain) from a per-domain attribute.
    `per_domain_attr_path` is something like:
      "barrier3d[].Shoreline"   or   "nourishments[].beach_width"
    Use '[]' as the domain placeholder. We’ll expand it for domain indices.

    Returns (T, D) ndarray.

    Steps:
    - We evaluate each domain's attribute and expect a 1D time series (or something
      that can be squeezed to 1D time).
    - We stack them as columns to build (time, domain).
    """
    # figure out domain count
    if domain_count is None:
        try:
            domain_count = len(cascade_obj.barrier3d)
        except Exception as e:
            raise RuntimeError(f"Cannot infer domain count: {e}")

    # helper to resolve an attribute path for a specific domain index
    def resolve_for_domain(obj, path_with_brackets, d):
        # replace [] with [d]
        path = path_with_brackets.replace("[]", f"[{d}]")
        # walk the path dotted / indexed
        cur = obj
        tokens = []
        # Split on dots but keep [index] and ['key'] segments attached
        # We'll manually parse indexing after each token
        for chunk in path.split("."):
            tokens.append(chunk)

        for tok in tokens:
            # handle indexing inside token if present
            # First, attribute name before any [ ... ]
            if "[" in tok:
                attr = tok.split("[", 1)[0]
                if attr:
                    cur = getattr(cur, attr)
                rest = tok[len(attr):]
            else:
                attr = tok
                rest = ""

                if attr:
                    cur = getattr(cur, attr)

            # Now apply any chained indices like [3]['key'][0]
            while rest.startswith("["):
                # find matching ]
                r = rest.find("]")
                key = rest[1:r]
                rest = rest[r+1:]
                if key.startswith("'") and key.endswith("'"):
                    key = key[1:-1]
                    cur = cur[key]
                elif key.isdigit():
                    cur = cur[int(key)]
                else:
                    raise ValueError(f"Unsupported index in token: [{key}]")

        return cur

    series_list = []
    for d in range(domain_count):
        arr = resolve_for_domain(cascade_obj, per_domain_attr_path, d)
        arr = np.asarray(arr)
        # squeeze to 1D time if possible
        if arr.ndim > 1:
            arr = np.squeeze(arr)
            if arr.ndim != 1:
                raise ValueError(f"Domain {d}: expected 1D time series, got shape {arr.shape}")
        series_list.append(arr)

    # ensure equal length time across domains
    lengths = [len(a) for a in series_list]
    if len(set(lengths)) != 1:
        raise ValueError(f"Time lengths differ across domains: {set(lengths)}")
    T = lengths[0]
    D = len(series_list)
    M = np.column_stack(series_list)  # (T, D)
    return M

# -------------------- How to use extract_time_by_domain --------------------
# 1) Run this script once to generate the catalog CSV.
# 2) Open the CSV and find a per-domain attribute path you care about.
#    Common guesses (may vary by branch):
#      - "cascade.barrier3d[].Shoreline" (or ".shoreline")
#      - "cascade.barrier3d[].BermElTime" (if berm elevation time series exists)
#      - "cascade.nourishments[].beach_width"
#
# 3) Then, uncomment and try an extraction like:
#
# try:
#     M = extract_time_by_domain(cascade_obj, "barrier3d[].Shoreline")
#     print("Extracted matrix shape (time, domain):", M.shape)
#     # Example: save to CSV next to the npz
#     out_csv = os.path.join(os.path.dirname(NPZ_PATH), "Shoreline_time_by_domain.csv")
#     import pandas as pd
#     pd.DataFrame(M, columns=[f"domain_{i}" for i in range(M.shape[1])]).assign(timestep=np.arange(M.shape[0])).to_csv(out_csv, index=False)
#     print("Saved:", out_csv)
# except Exception as e:
#     print("Extraction failed:", e)

