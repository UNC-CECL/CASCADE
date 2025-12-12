"""
Barrier3D Domain Inspector

Purpose:
    Inspect all attributes stored on a single Barrier3D domain object
    (cascade.barrier3d[domain_index]) in a CASCADE run.

What it does:
    - Loads a CASCADE object from a .npz (or pickle fallback).
    - Lets you pick a domain index.
    - Prints all public attributes on that Barrier3D segment with:
        * attribute name
        * type
        * if scalar (float/int): value
        * if array/list: shape, min, max, and sample values

Use this to:
    - Find dune crest elevation variables (e.g., z_dune, dune_crest_elev, etc.).
    - See what time series exist (TWL_TS, WL_TS, Hs_TS, _QowTS, etc.).
    - Generally understand what CASCADE is tracking for each domain.
"""

from pathlib import Path
import numpy as np
import pickle

# ==========================
# CONFIG: UPDATE THESE
# ==========================

# Path to your saved CASCADE object (.npz preferred)
PATH_TO_CASCADE_FILE = Path(
    r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
    r"\HAT_1978_1997_Natural_State\HAT_1978_1997_Natural_State.npz"
)

# Domain index to inspect (0-based index into cascade.barrier3d)
DOMAIN_INDEX = 40


# ==========================
# HELPERS
# ==========================

def load_cascade_from_file(cascade_path: Path):
    """Load a CASCADE object from a .npz with key 'cascade', or pickle as fallback."""
    if not cascade_path.exists():
        raise FileNotFoundError(f"CASCADE file not found at: {cascade_path}")

    suffix = cascade_path.suffix.lower()

    if suffix == ".npz":
        data = np.load(cascade_path, allow_pickle=True)
        if "cascade" in data:
            cascade = data["cascade"].item()
            print("Loaded CASCADE object from npz key 'cascade'.")
            return cascade
        else:
            keys = list(data.keys())
            raise KeyError(
                f"'cascade' key not found in npz file.\n"
                f"Available keys: {keys}"
            )

    # Fallback: try pickle
    with cascade_path.open("rb") as f:
        cascade = pickle.load(f)
    print("Loaded CASCADE object using pickle.load(...).")
    return cascade


def summarize_value(val, max_items: int = 5) -> str:
    """
    Create a short text summary of a value:
      - scalar: value
      - list/array: shape, min, max, a few sample values
      - other: type
    """
    try:
        if isinstance(val, (int, float, np.floating)):
            return f"scalar, value={float(val):.4f}"

        if isinstance(val, (list, tuple, np.ndarray)):
            arr = np.array(val, dtype=float)
            # handle empty or weird arrays safely
            if arr.size == 0:
                return "array/list, EMPTY"
            s = (
                f"array/list, shape={arr.shape}, "
                f"min={np.nanmin(arr):.4f}, max={np.nanmax(arr):.4f}"
            )
            # sample first few values flattened
            flat = arr.flatten()
            n_sample = min(max_items, flat.size)
            samples = ", ".join(f"{x:.4f}" for x in flat[:n_sample])
            s += f", samples=[{samples}{'...' if flat.size > n_sample else ''}]"
            return s

        # Fallback: just show type
        return f"type={type(val)}"

    except Exception as e:
        return f"<error summarizing value: {e}>"


def dune_domain_extra_info(arr: np.ndarray) -> str:
    """
    Extra diagnostic string for DuneDomain/_DuneDomain arrays
    (assumes shape ~ (time, cross_shore, field)).
    """
    try:
        dune = np.array(arr, dtype=float)
        if dune.ndim != 3:
            return "    [DuneDomain extra] Not 3D, cannot interpret as (time, x, field)."

        t0 = 0
        n_fields = dune.shape[2]
        lines = [f"    [DuneDomain extra] shape={dune.shape}"]

        if t0 >= dune.shape[0]:
            lines.append("    [DuneDomain extra] No time steps available.")
            return "\n".join(lines)

        for fi in range(n_fields):
            field = dune[t0, :, fi]
            lines.append(
                f"    - t=0, field {fi}: min={np.nanmin(field):.4f}, "
                f"max={np.nanmax(field):.4f}, mean={np.nanmean(field):.4f}"
            )

        return "\n".join(lines)

    except Exception as e:
        return f"    [DuneDomain extra] Error interpreting DuneDomain: {e}"


# ==========================
# MAIN
# ==========================

if __name__ == "__main__":
    # Load CASCADE
    cascade = load_cascade_from_file(PATH_TO_CASCADE_FILE)

    # Get barrier3d list
    barrier_models = cascade.barrier3d
    n_domains = len(barrier_models)
    print(f"\nCASCADE has {n_domains} alongshore domains (0–{n_domains - 1}).")

    if not (0 <= DOMAIN_INDEX < n_domains):
        raise IndexError(
            f"DOMAIN_INDEX {DOMAIN_INDEX} is out of range. "
            f"Must be between 0 and {n_domains - 1}."
        )

    seg = barrier_models[DOMAIN_INDEX]
    print(f"\nInspecting barrier3d[{DOMAIN_INDEX}]:\n")

    # Print attributes sorted alphabetically
    attrs = [a for a in dir(seg) if not a.startswith("__")]
    attrs.sort()

    for attr in attrs:
        try:
            val = getattr(seg, attr)
        except Exception as e:
            print(f"  {attr}: <error accessing attribute: {e}>")
            continue

        summary = summarize_value(val)
        print(f"  {attr}: {summary}")

        # Extra dune-specific diagnostics for DuneDomain / _DuneDomain
        if attr in ("DuneDomain", "_DuneDomain"):
            print(dune_domain_extra_info(val))
