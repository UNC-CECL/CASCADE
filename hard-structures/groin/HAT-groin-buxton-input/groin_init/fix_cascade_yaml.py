#!/usr/bin/env python3
"""
Repair a CASCADE parameter YAML that numpy scalars were serialised into, so yaml.full_load reads it again.

    python fix_cascade_yaml.py

Rewrites PARAM_FILE in place with plain numbers, after a .bak copy beside it.
Needs numpy and PyYAML. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
from pathlib import Path

import os
import shutil
import numpy as np
import yaml

# --- CONFIG ------------------------------------------------------------------
_PATH_REPO = next(_p for _p in Path(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())

# The parameter YAML CASCADE fails to load
PARAM_FILE = str(_PATH_REPO / "data" / "hatteras_init" / "Hatteras-CASCADE-parameters.yaml")
# -----------------------------------------------------------------------------


# Numpy scalars and arrays, recursively, back to plain Python types
def to_plain(obj):
    if isinstance(obj, dict):
        return {k: to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_plain(v) for v in obj]
    if isinstance(obj, np.generic):     # np.float64, np.int64, ... -> python scalar
        return obj.item()
    if isinstance(obj, np.ndarray):
        return to_plain(obj.tolist())
    return obj


# Run: back up, load with numpy objects, strip them, rewrite, prove full_load works
def main():
    if not os.path.isfile(PARAM_FILE):
        raise SystemExit(f"File not found: {PARAM_FILE}")

    # Back up the original first
    backup = PARAM_FILE + ".bak"
    shutil.copyfile(PARAM_FILE, backup)
    print(f"Backup written: {backup}")

    # Load with UnsafeLoader, which can rebuild the numpy objects (README)
    with open(PARAM_FILE, "r") as f:
        doc = yaml.load(f, Loader=yaml.UnsafeLoader)

    # Strip every numpy type back to a plain Python number
    doc_clean = to_plain(doc)

    # Rewrite with safe_dump, plain scalars only, key order kept
    with open(PARAM_FILE, "w") as f:
        yaml.safe_dump(doc_clean, f, default_flow_style=False, sort_keys=False)

    # Prove it: full_load, which CASCADE uses, must now succeed
    with open(PARAM_FILE, "r") as f:
        yaml.full_load(f)

    print(f"Repaired: {PARAM_FILE}")
    print("full_load now succeeds -- CASCADE can read the file again.")
    print("If anything looks off, the original is preserved at the .bak path.")


if __name__ == "__main__":
    main()
