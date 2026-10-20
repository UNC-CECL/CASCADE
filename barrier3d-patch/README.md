# Barrier3D patch for the Hatteras hindcast

The hindcast runs on published Barrier3D (UNC-CECL/Barrier3D, `master` at
`ce36866`) with one change, kept here so that it travels with this branch. It
has not yet been proposed to the Barrier3D repository.

`barrier3d-hatteras.patch` changes `barrier3d/barrier3d.py` only:

| | where | status |
|---|---|---|
| 1 | `route_overwash`, the subaerial test | **applied.** It read `Elevation[TS, i, d + 1: d + 10]`, with the axes swapped, for `[TS, d + 1: d + 10, i]`. It read the wrong cells, and past the array whenever `i` ≥ rows, which crashed runs with long storms (Isabel 2003) |
| 2 | `DuneGaps` | comment only: the last gap loses its last cell, and a lone overtopped cell gives no gap |
| 3 | gap discharge, `Discharge[:, 0, start:stop]` | comment only: `stop` is inclusive, so each gap's last cell gets no water |
| 4 | `C = 0  # Initialize` in the overwash loop | comment only: it overwrites the inundation momentum constant (added in b11b880, 2024) |

Each change is marked `HAH:` in the code. The full record, with evidence, is
`HATTERAS_FIXES.md` on the local Barrier3D branch `hatteras/adopted`; the
CASCADE side is [`HATTERAS_CASCADE_CHANGES.md`](../HATTERAS_CASCADE_CHANGES.md).

## Install

From the folder that holds `CASCADE/`:

```
git clone https://github.com/UNC-CECL/Barrier3D.git
cd Barrier3D
git checkout ce36866
git apply ../CASCADE/barrier3d-patch/barrier3d-hatteras.patch
pip install -e .
```

`git apply --check` first shows whether the patch still fits; it was checked
against `ce36866` on 2026-10-20 and reproduces `hatteras/adopted` exactly.

## Remaking the patch

After changing Barrier3D on `hatteras/adopted`:

```
git -C ../Barrier3D diff ce36866 HEAD -- barrier3d/barrier3d.py > barrier3d-patch/barrier3d-hatteras.patch
```
