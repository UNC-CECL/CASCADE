"""
Quick diagnostics for a CASCADE run:
- Check units (meters vs dam)
- Check whether overwash happens
- Compare model shoreline trend at a domain vs DSAS LRR

EDIT THE CONFIG SECTION ONLY, then run this file.
"""

import numpy as np
from cascade.cascade import Cascade

# =======================
# CONFIG – EDIT THESE
# =======================

# 1. Path to your CASCADE run directory
#    (where cascade_configuration.json and the outputs live)
RUN_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1997_Natural_State"

# 2. Buffer / domain setup
LEFT_BUFFER_DOMAINS = 15      # number of left buffer domains
FIRST_REAL_DOMAIN_ID = 1      # ID of your first real domain file (e.g., domain_1, domain_30, etc.)

# 3. Domain you want to compare to DSAS (e.g., "domain 80")
REAL_DOMAIN_ID_TO_CHECK = 80

# 4. DSAS long-term rate for that domain (m/yr) – fill this in from your DSAS table
DSAS_LRR_FOR_DOMAIN = -2.5    # <-- CHANGE THIS to your actual DSAS LRR value

# =======================
# HELPER FUNCTIONS
# =======================

def real_domain_id_to_index(real_domain_id: int) -> int:
    """
    Convert a real domain ID (like 80) to the index in cascade.barrier3d
    assuming barrier3d[0:LEFT_BUFFER_DOMAINS] are buffer domains.
    """
    return LEFT_BUFFER_DOMAINS + (real_domain_id - FIRST_REAL_DOMAIN_ID)


def get_first_available_attr(obj, candidates):
    """
    Try a list of possible attribute names and return the first one that exists,
    or None if none exist.
    """
    for name in candidates:
        if hasattr(obj, name):
            return name
    return None


# =======================
# LOAD CASCADE RUN
# =======================

print("Loading CASCADE run from:", RUN_DIR)
cas = Cascade(RUN_DIR)
barrier3d = cas.barrier3d
print(f"Loaded {len(barrier3d)} barrier3d domains.\n")

# =======================
# 1. UNITS CHECK (METERS VS DAM)
# =======================

b0 = barrier3d[0]  # just look at the first domain as an example

# Try to find elevation array
elev_attr = get_first_available_attr(b0, ["eta", "z", "elev", "zeta"])
if elev_attr is None:
    raise AttributeError("Could not find an elevation attribute (eta/z/elev/zeta) on barrier3d[0].")

z = getattr(b0, elev_attr)

dx = getattr(b0, "dx", None)

print("=== UNITS CHECK ===")
print(f"Cross-shore grid spacing dx: {dx}")
print(f"Elevation attribute: {elev_attr}")
print(f"Elevation min/max: {np.nanmin(z):.2f}, {np.nanmax(z):.2f} (same units as dx)")

print("\nInterpretation:")
print("- If dx is ~10 and elevations look like 0–10 or so, everything is almost certainly in METERS.")
print("- If elevations looked like 0–1 with dx ~1, that might suggest decameters, but that would be unusual here.\n")

# =======================
# 2. OVERWASH FLUX CHECK
# =======================

# Pick a central domain to inspect. We'll use the same one as the DSAS comparison.
domain_index = real_domain_id_to_index(REAL_DOMAIN_ID_TO_CHECK)
b_center = barrier3d[domain_index]

# Try to find an overwash flux attribute
ow_attr_candidates = [
    "q_ow_TS", "Q_ow_TS", "Q_overwash_TS", "q_overwash_TS",
    "overwash_flux_TS", "Qow_TS"
]
ow_attr = get_first_available_attr(b_center, ow_attr_candidates)

print("=== OVERWASH CHECK ===")
if ow_attr is None:
    print("Could not automatically find an overwash flux time series.")
    print("Attributes on the domain containing 'ow' or 'overwash':")
    ow_like = [a for a in dir(b_center) if ("ow" in a.lower()) or ("overwash" in a.lower())]
    print(ow_like)
    print("\nYou can pick the correct one from the list above and add it to ow_attr_candidates.")
else:
    ow = getattr(b_center, ow_attr)
    ow_min = np.nanmin(ow)
    ow_max = np.nanmax(ow)
    any_ow = np.any(ow > 0)

    print(f"Overwash flux attribute: {ow_attr}")
    print(f"Overwash flux min/max: {ow_min:.6f}, {ow_max:.6f}")
    print(f"Any overwash > 0? {any_ow}")

    print("\nInterpretation:")
    if not any_ow or np.isclose(ow_max, 0.0):
        print("- Overwash flux is essentially zero for this domain in this run.")
        print("- That supports the idea that storms are not overtopping the dunes.")
    else:
        print("- There IS overwash in this run (overwash flux > 0 at some times).")
        print("- We could dig deeper into when/how much if needed.\n")

# Try to inspect dune crest and water level if they exist
dune_attr = get_first_available_attr(b_center, ["z_dune_TS", "zc_TS", "z_crest_TS"])
wl_attr = get_first_available_attr(b_center, ["eta_TS", "wl_TS", "waterlevel_TS", "z_water_TS"])

print("\n=== DUNE VS WATER LEVEL (if available) ===")
if dune_attr:
    dune = getattr(b_center, dune_attr)
    print(f"Dune crest attribute: {dune_attr}")
    print(f"Dune crest min/max: {np.nanmin(dune):.2f}, {np.nanmax(dune):.2f}")
else:
    print("No dune crest time series found (z_dune_TS / zc_TS / z_crest_TS).")

if wl_attr:
    wl = getattr(b_center, wl_attr)
    print(f"Water level attribute: {wl_attr}")
    print(f"Water level min/max: {np.nanmin(wl):.2f}, {np.nanmax(wl):.2f}")
else:
    print("No water level time series found (eta_TS / wl_TS / waterlevel_TS / z_water_TS).")

print("\nInterpretation:")
print("- If dune crest min is well ABOVE water level max, storms likely don't overtop.")
print("- If water level occasionally approaches or exceeds dune crest, overwash SHOULD be happening.\n")

# =======================
# 3. SHORELINE TREND VS DSAS AT REAL_DOMAIN_ID_TO_CHECK
# =======================

print("=== SHORELINE TREND VS DSAS ===")
print(f"Real domain ID to check: {REAL_DOMAIN_ID_TO_CHECK}")
print(f"Index in barrier3d: {domain_index}")

b = b_center  # same object

# Find time and shoreline position attributes
t_attr = get_first_available_attr(b, ["t_TS", "time_TS", "t"])
xs_attr = get_first_available_attr(b, ["x_s_TS", "xs_TS", "shoreline_x_TS"])

if t_attr is None or xs_attr is None:
    print("Could not find time or shoreline position series.")
    print("Time-like attributes:", [a for a in dir(b) if "t_TS" in a or "time" in a])
    print("Shoreline-like attributes:", [a for a in dir(b) if "x_s" in a or "shoreline" in a])
else:
    t = getattr(b, t_attr)
    x = getattr(b, xs_attr)

    # Ensure arrays are 1D and same length
    t = np.asarray(t).flatten()
    x = np.asarray(x).flatten()

    # Remove NaNs (if any)
    valid = np.isfinite(t) & np.isfinite(x)
    t_valid = t[valid]
    x_valid = x[valid]

    if len(t_valid) < 2:
        print("Not enough valid points to fit a trend.")
    else:
        coeffs = np.polyfit(t_valid, x_valid, 1)
        slope_m_per_yr = coeffs[0]

        print(f"Model shoreline trend at domain {REAL_DOMAIN_ID_TO_CHECK}: {slope_m_per_yr:.3f} m/yr")
        print(f"DSAS LRR at this domain: {DSAS_LRR_FOR_DOMAIN:.3f} m/yr")
        diff = slope_m_per_yr - DSAS_LRR_FOR_DOMAIN
        print(f"Model - DSAS difference: {diff:.3f} m/yr")

        print("\nInterpretation:")
        if slope_m_per_yr < 0:
            print("- Model shoreline is on average EROSIONAL (retreating).")
        else:
            print("- Model shoreline is on average ACCRETIONAL (advancing).")

        if DSAS_LRR_FOR_DOMAIN < 0:
            print("- Observed DSAS rate is EROSIONAL.")
        else:
            print("- Observed DSAS rate is ACCRETIONAL.")

        print(f"- The model is {'more' if abs(slope_m_per_yr) > abs(DSAS_LRR_FOR_DOMAIN) else 'less'} erosional than observed (in magnitude).")
        print()
