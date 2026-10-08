"""
Why the 1996 no-groin run erases the Buxton fillet: every matrix run, then alongshore diffusion alone.

    python hard-structures/groin/2-module-tests/3-real-planform/diagnose_no_groin_relaxation.py

Run from the repo root: its paths are relative. Part 1 reads every adopted
matrix run's shoreline matrix (edge forcing x management); part 2 runs
groin_stability_option_a's emulator with no groin (no Barrier3D). Prints only;
the 2026-09-29 output is diagnose_no_groin_relaxation.txt.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
import numpy as np, glob, os

# --- CONFIG ------------------------------------------------------------------
UP, DOWN = 20, 19                          # GIS 6 / GIS 5 in the padded array
# -----------------------------------------------------------------------------


# Gap change over the run, OLS slope x run length
def tc(s):
    t=np.arange(s.size); return np.polyfit(t,s,1)[0]*(s.size-1)


# Part 1: D5-D6 change across every adopted matrix run
for per in ("1996_2010","2010_2024"):
    print(f"\n== {per}")
    for d in sorted(glob.glob(f"output/raw_runs/matrix/{per}/*/*/")):
        n=os.path.basename(d.rstrip("/\\"))
        f=os.path.join(d,n+"_shoreline_matrix.npy")
        if not os.path.exists(f): continue
        x=np.load(f); g=x[:,DOWN]-x[:,UP]
        print(f"{n[14:]:48s} gap0 {g[0]:+7.1f}  total {tc(g):+6.1f}  endpt {g[-1]-g[0]:+6.1f} | dx D5 {x[-1,DOWN]-x[0,DOWN]:+6.1f} D6 {x[-1,UP]-x[0,UP]:+6.1f} D4 {x[-1,18]-x[0,18]:+6.1f} D7 {x[-1,21]-x[0,21]:+6.1f}")

# Part 2: the same with alongshore diffusion alone (the emulator, no Barrier3D)
import sys
sys.path.insert(0, r"hard-structures/groin/2-module-tests/3-real-planform")
sys.path.insert(0, r"hard-structures/groin/2-module-tests/2-solver-audit")
import groin_stability_option_a as g
for start in (1996, 2010):
    x0 = g.load_planform(start)
    b = g.solve(x0, g.CLIMATES["optionA"], 14, lambda t: 0.0)
    gap = b.gap_m.values
    print(start, "emulator alongshore-only: gap0 %.1f  gap14 %.1f  change %.1f | D5 %+.1f D6 %+.1f" % (
        x0[19]-x0[20], gap[-1], gap[-1]-(x0[19]-x0[20]), b.x_down.iloc[-1]-x0[19], b.x_up.iloc[-1]-x0[20]))
