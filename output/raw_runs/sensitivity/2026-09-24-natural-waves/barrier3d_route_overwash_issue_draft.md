# DRAFT — not posted. For review before opening an issue on UNC-CECL/Barrier3D.

---

**Title:** `route_overwash`: row and column swapped in the subaerial check (line 1092) — reads the wrong cells, and out of bounds on narrow domains

## Summary

In `Barrier3D.route_overwash`, the test that decides whether overwash flux keeps
moving landward or is deposited into the bay indexes `Elevation` with its row and
column swapped in the second clause:

```python
# barrier3d/barrier3d.py, route_overwash (line 1092 in our copy)
if Elevation[TS, d, i] > SL or np.sum(np.greater(Elevation[TS, i, d + 1: d + 10], SL)) > 0:
```

Everywhere else in the routine `d` is the cross-shore row (axis 1) and `i` the
alongshore column (axis 2), e.g. two lines below:

```python
SedFluxIn[TS, d + 1, i] += Qs2
```

The clause appears intended to ask whether any of the next nine cells landward of
`(d, i)` is subaerial, i.e.

```python
np.sum(np.greater(Elevation[TS, d + 1: d + 10, i], SL)) > 0
```

`git blame` attributes the line to commit b11b8808 (2024-03-01).

## Effects

1. **Wrong cells, silently.** When `i` is less than the number of rows, the clause
   reads row `i`, columns `d+1 … d+9`: a strip of the domain unrelated to the flow
   path. The subaerial/subaqueous branch is then decided on the wrong cells, which
   changes whether sediment is routed onward or decayed into the bay.
2. **Out-of-bounds reads.** When `i` is not less than the number of rows (a
   domain narrower cross-shore than it is wide alongshore, e.g. an interior of 45
   rows with `BarrierLength = 50`), `Elevation[TS, i, …]` is out of bounds on
   axis 1. `route_overwash` is `numba`-jitted with bounds checking off, so this
   returns arbitrary memory and occasionally crashes the process with an access
   violation / segmentation fault and no Python traceback.

## Reproduction

With bounds checking on, the out-of-bounds case raises instead of reading garbage:

```
NUMBA_BOUNDSCHECK=1 python <any run with overwash on a domain with fewer rows than BarrierLength>
...
  File ".../barrier3d/barrier3d.py", line 1586, in update
    ) = self.route_overwash(
IndexError: index is out of bounds
```

With the JIT disabled (`NUMBA_DISABLE_JIT=1`) the message is exact:

```
  File ".../barrier3d/barrier3d.py", line 1092, in route_overwash
    if Elevation[TS, d, i] > SL or np.sum(np.greater(Elevation[TS, i, d + 1: d + 10], SL)) > 0:
IndexError: index 45 is out of bounds for axis 1 with size 45
```

We hit it in CASCADE hindcasts (BRIE-coupled, 90 domains of 50 alongshore cells,
interiors 29–189 rows) in both a default and a modified configuration; in
unchecked runs it surfaced only as occasional silent crashes, the earliest
observed years after the first out-of-bounds read.

## Suggested fix

```diff
-                    if Elevation[TS, d, i] > SL or np.sum(np.greater(Elevation[TS, i, d + 1: d + 10], SL)) > 0:
+                    if Elevation[TS, d, i] > SL or np.sum(np.greater(Elevation[TS, d + 1: d + 10, i], SL)) > 0:
```

(`d + 1: d + 10` near the bay edge is clipped by slicing, as before.)

## Effect of the fix, measured

Eight CASCADE hindcasts (14 years, 90 domains) re-run with the fix against
their unpatched twins:

- no out-of-bounds read remains (two runs completed under `NUMBA_BOUNDSCHECK=1`);
- a run that crashed unpatched in year 13 completes;
- skill scores against observed shoreline change moved in the third decimal
  (RMSE 1.141 → 1.143, 4.874 → 4.870, …);
- per-domain shoreline-change rates moved by at most 0.11 m/yr in six of the
  seven comparable pairs; in the run with the most overwash, 14 of 88 domains moved by more than
  0.1 m/yr and one by 2.2 m/yr.

So the wrong-cell reads mostly returned the same answer; the practical cost was
the crashes.
