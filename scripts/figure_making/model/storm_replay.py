"""
storm_replay.py
==============================================================================
Replay storms through a saved Barrier3D domain's own `update()` and read the
routing arrays out of it: water discharge, sediment flux in and out, and the
elevation of every cell for every routing step, which the model computes and
throws away. Used by model_mechanics_figures.py and overwash_routing_figures.py.

The model code is not copied or modified. A line trace on the one update()
frame snapshots its locals at the first statement after each storm's routing
loop, and stops the update once the last storm is captured.

Two facts about Barrier3D that the replay depends on, both checked against a
saved run to 0.0 m (overwash_routing_figures.storm_routing_check):
  - SeaLevel() lowers DuneDomain[t - 1] IN PLACE at the start of update t, so
    a saved object's dune slice t - 1 has already lost year t's sea-level rise
    and must have it put back. The interior is lowered into a new array, so
    DomainTS[t - 1] is the true starting grid.
  - The dune crest that decides which cells overwash (DuneDomainCrest,
    barrier3d.py:1358) is computed ONCE per year, after dune growth and
    before the first storm. Every storm that year is tested against it and
    routed over it, although each storm also lowers the dune it erodes.
"""

from __future__ import annotations

import copy
import functools
import inspect
import sys

import numpy as np
from barrier3d import Barrier3d

DAM = 10.0


class _StormsDone(Exception):
    pass


CAPTURE_TEXT = "InteriorUpdate = Elevation[-1, 1:, :]"

# THREE SUSPECTED DEFECTS in how Barrier3D starts overwash, found with this
# replay on 2026-09-27. None is fixed in barrier3d.py; a VARIANT applies a fix
# to an in-memory copy of the class only, so their effect can be measured.
#   "gaps"     DuneGaps() (2020) drops the last overtopped cell of the last
#              gap, and returns nothing at all when one cell is overtopped.
#   "slice"    update() sets gap water with Discharge[:, 0, start:stop], but
#              stop is inclusive: every gap loses its last cell of water, and
#              a one-cell gap gets none (2020).
#   "momentum" update() computes the inundation momentum constant
#              C = Cx * AvgSlope before the gap loop, then resets C = 0 inside
#              it (b11b880, 2024 Numba refactor), so inundation transport
#              Ki * (Q * (S + C))**mm always runs with C = 0.
FIXES = ("gaps", "slice", "momentum")
_SOURCE_FIXES = {
    "slice": ("Discharge[:, 0, start:stop] = Qdune", "Discharge[:, 0, start:stop + 1] = Qdune"),
    "momentum": ("C = 0  # Initialize", "pass  # C kept (replay variant)"),
}


def _dune_gaps_fixed(self, DuneDomain, Dow, bermel, Rhigh):
    """Contiguous runs of overtopped cells, every cell kept."""
    gaps = []
    if not len(Dow):
        return gaps
    runs, start = [], Dow[0]
    for a, b in zip(Dow[:-1], Dow[1:]):
        if b - a != 1:
            runs.append((start, a))
            start = b
    runs.append((start, Dow[-1]))
    for a, b in runs:
        x = DuneDomain[a:b + 1]
        hmean = sum(x) / float(len(x))
        gaps.append([a, b, Rhigh - (hmean + bermel), Rhigh / (hmean + bermel)])
    return gaps


@functools.lru_cache(maxsize=None)
def model_class(fixes=()):
    """(class, update code object, capture line). fixes=() is Barrier3d
    itself; otherwise a subclass whose update() is compiled from the model's
    own source with the named text substitutions."""
    if not fixes:
        lines, first = inspect.getsourcelines(Barrier3d.update)
        for i, text in enumerate(lines):
            if CAPTURE_TEXT in text:
                return Barrier3d, Barrier3d.update.__code__, first + i
        raise RuntimeError("Barrier3d.update has changed: capture line not found")
    import textwrap
    import barrier3d.barrier3d as b3dmod
    src = textwrap.dedent(inspect.getsource(Barrier3d.update))
    for f in fixes:
        if f in _SOURCE_FIXES:
            old, new = _SOURCE_FIXES[f]
            if src.count(old) != 1:
                raise RuntimeError(f"fix {f!r}: expected one {old!r} in update()")
            src = src.replace(old, new)
    ns = {}
    exec(compile(src, f"<Barrier3d.update {'+'.join(fixes)}>", "exec"), vars(b3dmod), ns)
    attrs = {"update": ns["update"]}
    if "gaps" in fixes:
        attrs["DuneGaps"] = _dune_gaps_fixed
    cls = type(f"Barrier3dFix_{'_'.join(fixes)}", (Barrier3d,), attrs)
    line = next(i + 1 for i, text in enumerate(src.splitlines()) if CAPTURE_TEXT in text)
    return cls, ns["update"].__code__, line


def replay(b_saved, t, storms, fixes=()):
    """Run `storms` (rows of Rhigh m MHW, Rlow m MHW, period s, duration h)
    through model year `t` of a saved Barrier3D domain, starting from the
    grid the run saved entering that year. Returns one dict per storm."""
    b = copy.deepcopy(b_saved)
    cls, code, line = model_class(tuple(sorted(fixes)))
    b.__class__ = cls
    b._time_index = t
    b._InteriorDomain = np.array(b.DomainTS[t - 1], dtype=float).copy()
    # SeaLevel() lowers DuneDomain[t - 1] IN PLACE at the start of update t
    # (barrier3d.py:19), so the saved object's slice t - 1 has already had
    # year t's sea-level rise taken off. Put it back, or the replay lowers the
    # dune twice. The interior is lowered into a new array, so DomainTS[t - 1]
    # is the true starting grid. (A 4 mm slip here moved individual cells by
    # up to 0.5 m while the storm's total overwash changed by 1%: the routing
    # is cell-scale sensitive, the volumes are not.)
    b._DuneDomain[t - 1] = b._DuneDomain[t - 1] + b._RSLR[t]
    b._StormSeries = np.array([[t, rh / DAM, rl / DAM, per, dur] for rh, rl, per, dur in storms], dtype=float)
    got = []

    def local(frame, event, arg):
        if event == "line" and frame.f_lineno == line:
            L = frame.f_locals
            got.append(dict(
                elevation=np.array(L["Elevation"]) * DAM,          # (hours*substep, rows, cols) m MHW
                discharge=np.array(L["Discharge"]) * 1000.0,       # m3/hr through each cell
                flux_in=np.array(L["SedFluxIn"]) * 1000.0,         # m3 per step
                flux_out=np.array(L["SedFluxOut"]) * 1000.0,
                inundation=int(L["inundation"]), substep=int(L["substep"]),
                gaps=list(L["gaps"]), rhigh=float(L["Rhigh"][L["n"]]) * DAM,
                rlow=float(L["Rlow"][L["n"]]) * DAM,
                crest_pre=(np.array(L["Dunes_prestorm"]) + b.BermEl) * DAM,
                owloss_cum=float(L["OWloss"]) * 1000.0 / b._BarrierLength / DAM,   # m3/m, year so far
            ))
            prev = got[-2]["owloss_cum"] if len(got) > 1 else 0.0
            got[-1]["owloss"] = got[-1]["owloss_cum"] - prev                     # this storm
            if len(got) == len(storms):
                raise _StormsDone
        return local

    def glob(frame, event, arg):
        return local if (event == "call" and frame.f_code is code) else None

    sys.settrace(glob)
    try:
        b.update()
    except _StormsDone:
        pass
    finally:
        sys.settrace(None)
    return got, b


def regime(s):
    if not s["gaps"]:
        return "collision"
    return "inundation" if s["inundation"] else "run-up"
