"""
Replay storms through a saved Barrier3D domain's own update() and read the routing arrays it discards.

    from storm_replay import replay, regime

Used by model_mechanics_figures.py and overwash_routing_figures.py. The model
code is not copied; `defects=` puts a fixed Barrier3D defect back in memory.
Details: scripts/figure_making/model/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import copy
import functools
import inspect
import sys

import numpy as np
from barrier3d import Barrier3d

DAM = 10.0


# Raised to stop update() once the last storm is captured
class _StormsDone(Exception):
    pass


CAPTURE_TEXT = "InteriorUpdate = Elevation[-1, 1:, :]"

# The three overwash defects, fixed in hatteras/adopted; `defects=` puts them back (README)
FIXES = ("gaps", "slice", "momentum")
DEFECTS = FIXES
_SOURCE_FIXES = {
    "slice": ("Discharge[:, 0, start:stop] = Qdune", "Discharge[:, 0, start:stop + 1] = Qdune"),
    "momentum": ("C = 0  # Initialize", "pass  # C kept (replay variant)"),
}
# Reverse substitutions: hatteras/adopted text back to upstream
_INUNDATION_IF = "if inundation == 1:  # Inundation regime"
_SOURCE_DEFECTS = {
    "slice": ("Discharge[:, 0, start:stop + 1] = Qdune", "Discharge[:, 0, start:stop] = Qdune"),
    "momentum": (_INUNDATION_IF, "C = 0  # upstream reset (replay variant)\n{indent}" + _INUNDATION_IF),
}


# DuneGaps with every overtopped cell kept
def _dune_gaps_fixed(self, DuneDomain, Dow, bermel, Rhigh):
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


# DuneGaps as upstream Barrier3D has it (drops the last cell, and lone cells)
def _dune_gaps_upstream(self, DuneDomain, Dow, bermel, Rhigh):
    gaps = []
    start = 0
    i = start
    while i < (len(Dow) - 1):
        adjacent = Dow[i + 1] - Dow[i]
        if adjacent == 1:
            i = i + 1
        else:
            stop = i
            x = DuneDomain[Dow[start] : (Dow[stop] + 1)]
            Hmean = sum(x) / float(len(x))
            gaps.append([Dow[start], Dow[stop], Rhigh - (Hmean + bermel), Rhigh / (Hmean + bermel)])
            start = stop + 1
            i = start
    if i > 0:
        stop = i - 1
        x = DuneDomain[Dow[start] : (Dow[stop] + 1)]
        if len(x) > 0:
            Hmean = sum(x) / float(len(x))
            gaps.append([Dow[start], Dow[stop], Rhigh - (Hmean + bermel), Rhigh / (Hmean + bermel)])
    return gaps


# (class, update code, capture line): Barrier3d, or a subclass compiled with fixes/defects
@functools.lru_cache(maxsize=None)
def model_class(fixes=(), defects=()):
    if set(fixes) & set(defects):
        raise ValueError(f"both fixed and restored: {set(fixes) & set(defects)}")
    if not fixes and not defects:
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
            if src.count(new.split("  #")[0]) == 1 or (f == "momentum" and old not in src):
                continue                    # already fixed in this Barrier3D
            if src.count(old) != 1:
                raise RuntimeError(f"fix {f!r}: expected one {old!r} in update()")
            src = src.replace(old, new)
    for f in defects:
        if f in _SOURCE_DEFECTS:
            old, new = _SOURCE_DEFECTS[f]
            if src.count(old) != 1:
                raise RuntimeError(f"defect {f!r}: expected one {old!r} in update() "
                                   "(is ../Barrier3D on hatteras/adopted?)")
            line = next(ln for ln in src.splitlines() if old in ln)
            new = new.format(indent=line[: len(line) - len(line.lstrip())])
            src = src.replace(old, new)
    tag = "+".join([f"fix-{f}" for f in fixes] + [f"upstream-{f}" for f in defects])
    ns = {}
    exec(compile(src, f"<Barrier3d.update {tag}>", "exec"), vars(b3dmod), ns)
    attrs = {"update": ns["update"]}
    if "gaps" in fixes:
        attrs["DuneGaps"] = _dune_gaps_fixed
    if "gaps" in defects:
        attrs["DuneGaps"] = _dune_gaps_upstream
    cls = type(f"Barrier3dVariant_{tag.replace('+', '_').replace('-', '_')}", (Barrier3d,), attrs)
    line = next(i + 1 for i, text in enumerate(src.splitlines()) if CAPTURE_TEXT in text)
    return cls, ns["update"].__code__, line


# Run storms (Rhigh, Rlow m MHW, period s, duration h) through model year t of a saved domain
def replay(b_saved, t, storms, fixes=(), defects=()):
    b = copy.deepcopy(b_saved)
    cls, code, line = model_class(tuple(sorted(fixes)), tuple(sorted(defects)))
    b.__class__ = cls
    b._time_index = t
    b._InteriorDomain = np.array(b.DomainTS[t - 1], dtype=float).copy()
    # Put back year t's sea-level rise, which SeaLevel() took off the saved dune in place (README)
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


# The overwash regime name of a replayed storm
def regime(s):
    if not s["gaps"]:
        return "collision"
    return "inundation" if s["inundation"] else "run-up"
