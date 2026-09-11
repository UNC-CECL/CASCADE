#!/usr/bin/env python3
"""
fix_cascade_yaml.py
===================
One-shot repair for a CASCADE parameter YAML that has been "poisoned" with
NumPy scalar objects, e.g.:

    yaml.constructor.ConstructorError: could not determine a constructor for the
    tag 'tag:yaml.org,2002:python/object/apply:numpy._core.multiarray.scalar'

That happens when a numpy value (np.float64(...), etc.) got serialized into the
file instead of a plain number. CASCADE reads the file with yaml.full_load,
which refuses to reconstruct those objects, so every run crashes on load.

This script loads the file with a loader that CAN rebuild the numpy objects,
converts every numpy scalar/array back to a plain Python float/int/list, and
rewrites the file with plain numbers only. A .bak copy is made first.

USAGE
-----
Just run it. Edit PARAM_FILE below if your path differs.
    python fix_cascade_yaml.py
"""

import os
import shutil
import numpy as np
import yaml

# Path to the parameter YAML CASCADE is choking on (from your traceback).
PARAM_FILE = r"/data/hatteras_init/Hatteras-CASCADE-parameters.yaml"


def to_plain(obj):
    """Recursively convert numpy scalars/arrays to plain Python types."""
    if isinstance(obj, dict):
        return {k: to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_plain(v) for v in obj]
    if isinstance(obj, np.generic):     # np.float64, np.int64, ... -> python scalar
        return obj.item()
    if isinstance(obj, np.ndarray):
        return to_plain(obj.tolist())
    return obj


def main():
    if not os.path.isfile(PARAM_FILE):
        raise SystemExit(f"File not found: {PARAM_FILE}")

    # 1. Back up the original first (never overwrite without a copy).
    backup = PARAM_FILE + ".bak"
    shutil.copyfile(PARAM_FILE, backup)
    print(f"Backup written: {backup}")

    # 2. Load with UnsafeLoader -- unlike full_load, it can reconstruct the
    #    numpy scalar objects (numpy must be importable, which it is here).
    with open(PARAM_FILE, "r") as f:
        doc = yaml.load(f, Loader=yaml.UnsafeLoader)

    # 3. Strip every numpy type back to a plain Python number.
    doc_clean = to_plain(doc)

    # 4. Rewrite with safe_dump so only plain scalars are emitted. Key order is
    #    preserved; CASCADE reads by key, so order/comments don't matter to it.
    with open(PARAM_FILE, "w") as f:
        yaml.safe_dump(doc_clean, f, default_flow_style=False, sort_keys=False)

    # 5. Prove it: full_load (what CASCADE uses) must now succeed.
    with open(PARAM_FILE, "r") as f:
        yaml.full_load(f)

    print(f"Repaired: {PARAM_FILE}")
    print("full_load now succeeds -- CASCADE can read the file again.")
    print("If anything looks off, the original is preserved at the .bak path.")


if __name__ == "__main__":
    main()
