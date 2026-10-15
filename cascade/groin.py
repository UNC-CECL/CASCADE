
"""The Buxton groin for CASCADE: a share of the wave-driven drift held back, without BRIE leaking it away.

One callback, :class:`CorrectedGroinCallback`, attached to a :class:`Cascade` instance via
``cascade._groin_callback`` and invoked once per model year, immediately before BRIE's
alongshore-transport solve. It reads BRIE's shoreline and its own wave-climate tables, adds
to ``x_s_dt``, and changes nothing inside BRIE or Barrier3D.

Two numbers describe it at Buxton: a strength, the share of the arriving drift the groin
holds back (0.25), and the year it stops working (2017, the Buxton fill). Both were fitted on
one continuous 1967-2025 run against the wet/dry photo dates and annual CoastSat
(hard-structures/groin/4-hindcast/2-blocking-1996-2025/2026-10-13-pregroin-planform/).

The earlier groins this module held, a fixed -M/+M dipole and a blocking groin that cancelled
a share b of BRIE's diffusive exchange (each with a deterioration floor f), were removed on
2026-10-14. Their runs and fits are in hard-structures/groin/; git history keeps the code.

Conventions (verified against ``brie_coupler``)
-----------------------------------------------
Units
    ``x_s_dt`` is passed to BRIE in metres (``brie_coupler`` converts back to decametres).
Sign
    BRIE / Barrier3D ``x_s`` increases landward: a negative ``x_s_dt`` is a seaward advance.
Indexing
    The runner resolves 1-based GIS domain IDs to padded array indices via its own
    ``gis_to_pad`` translator and passes those indices in, so the padding width is defined in
    exactly one place (the runner), never here.

Author
------
Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

__all__ = ["CorrectedGroinCallback", "GroinOutsideTestedRange"]


class GroinOutsideTestedRange(ValueError):
    """The run has left the conditions the groin module was tested in; its results would not hold."""


def _check_pads(updrift_pad, downdrift_pad, n_domains):
    """Raise ValueError unless the two pads are adjacent and in range."""
    if abs(updrift_pad - downdrift_pad) != 1:
        raise ValueError(
            "groin domains must be adjacent (share the blocked boundary); "
            f"got updrift_pad={updrift_pad}, downdrift_pad={downdrift_pad}."
        )
    if not (0 <= updrift_pad < n_domains and 0 <= downdrift_pad < n_domains):
        raise ValueError(
            f"groin pads out of range for n_domains={n_domains}: "
            f"updrift_pad={updrift_pad}, downdrift_pad={downdrift_pad}."
        )



class CorrectedGroinCallback:
    """A groin that holds back a share of the wave-driven drift, without BRIE leaking it away.

    Two parts each year, both added to ``x_s_dt`` before BRIE's alongshore solve:

    1. **The drift groin.** BRIE's net alongshore drift for the wave climate,
       read from its own tables at the shoreline angle where the drift arrives
       at the face between the two pads, is converted to shoreline change
       (``q * dt / (dy * (h_b_crit + d_sf))``); a share ``s`` of it is moved
       from the downdrift pad to the updrift pad.
    2. **The leak correction.** BRIE's diffusive step would smooth the step the
       groin makes and spread the held sand onto the rest of the island. Within
       ``half_width`` cells of the face, the correction adds back
       ``2 (r_i - r_{i-1}) (x_i - x_{i-1})`` for this shoreline minus the same
       term for the groin-free shoreline, so it acts on the groin's own bends
       only, never on the coast's natural ones. The groin-free shoreline is this
       shoreline less the groin's tracked footprint: each year BRIE's alongshore
       step for last year (``brie.py`` lines 1293-1326) is repeated with and
       without the groin's own contribution, and the difference is the
       footprint. With ``s`` at 0 the footprint and the correction are exactly
       zero, so the run is the no-groin run.

    **The tested range.** The module is a simplified groin, tested on straight
    coasts and on Buxton (hard-structures/groin/3-module-tests/6-stress/). Each
    year it is installed it stops the run with :class:`GroinOutsideTestedRange`
    when the run leaves that range: the wave climate changes during the run;
    the groin makes or loses more than ``SAND_TOLERANCE_M`` (or
    ``SAND_TOLERANCE_FRACTION`` of its footprint) of summed shoreline, which
    happens where BRIE's diffusion at the groin is strong; or the groin face is
    steeper than ``MAX_FACE_ANGLE_DEG``, past which BRIE's diffusivity is clipped
    to zero there and the fillet grows without limit. It also refuses, when
    built, a correction zone that runs past either end of the domains.

    The strength is set per year by ``strength_history``, a list of
    ``(from_year, strength)`` pairs: the strength holds from each year until
    the next. A groin working at 0.25 until a 2017 fill is
    ``[(1969, 0.25), (2017, 0.0)]``. Fitted on Buxton over 1967-2025 in
    hard-structures/groin/4-hindcast/2-blocking-1996-2025/2026-10-13-pregroin-planform/.

    Parameters
    ----------
    updrift_pad, downdrift_pad : int
        Padded indices of the two domains sharing the face. Adjacent.
    strength_history : sequence of (int, float)
        From each year on, the share of the arriving drift held back, in [0, 1].
    start_year, install_year, n_domains
        As :class:`GroinCallback`.
    half_width : int
        Cells either side of the face the leak correction acts on (10: 5 leaks on
        a grid-aligned coast, hard-structures/groin/3-module-tests/5-long-term/).
    """
    kind: str = "corrected"
    # The tested range (hard-structures/groin/3-module-tests/6-stress/): outside it the run stops
    SAND_TOLERANCE_M: float = 10.0         # the groin may make or lose this much summed shoreline ...
    SAND_TOLERANCE_FRACTION: float = 0.05  # ... or this share of its footprint's gross, whichever is larger
    MAX_FACE_ANGLE_DEG: float = 45.0       # the groin face's steepest tested angle (Buxton's: 16-29 deg)

    def __init__(
        self,
        updrift_pad: int,
        downdrift_pad: int,
        strength_history,
        start_year: int,
        install_year: int,
        n_domains: int,
        half_width: int = 10,
    ) -> None:
        _check_pads(updrift_pad, downdrift_pad, n_domains)
        self._check_zone(min(updrift_pad, downdrift_pad), max(updrift_pad, downdrift_pad), half_width, n_domains)
        history = sorted((int(y), float(s)) for y, s in strength_history)
        if not history or any(not (0.0 <= s <= 1.0) for _, s in history):
            raise ValueError("strength_history needs (year, strength) pairs with strength in [0, 1].")
        self.updrift_pad: int = int(updrift_pad)
        self.downdrift_pad: int = int(downdrift_pad)
        self.strength_history = history
        self.start_year: int = int(start_year)
        self.install_year: int = int(install_year)
        self.n_domains: int = int(n_domains)
        self.half_width: int = int(half_width)
        self._lo: int = min(self.updrift_pad, self.downdrift_pad)
        self._hi: int = max(self.updrift_pad, self.downdrift_pad)

        self.year_TS: List[int] = []
        self.active_TS: List[bool] = []
        self.blocking_applied_TS: List[float] = []
        self.applied_dx_updrift_TS: List[float] = []
        self.applied_dx_downdrift_TS: List[float] = []
        self.trapping_rate_applied_TS: List[float] = []
        self.face_offset_m_TS: List[float] = []
        self.net_drift_m3_yr_TS: List[float] = []
        self.approach_angle_deg_TS: List[float] = []
        self.leak_correction_m_TS: List[float] = []
        self.footprint_max_m_TS: List[float] = []
        self._call_count: int = 0
        self._coast_diff0 = None
        self._q_net = None
        self._a = None
        self._prev = None

    def _effective_trapping_rate(self, year: int) -> float:
        """The strength in ``year``: 0 before install, else the history's value."""
        if year < self.install_year:
            return 0.0
        value = 0.0
        for y, s in self.strength_history:
            if year >= y:
                value = s
        return value

    @staticmethod
    def _drift_table(brie):
        """BRIE's net alongshore drift for its wave climate, per shoreline-angle bin."""
        step = 180.0 / (brie._wave_climl - 1)
        pdf = brie._angles.pdf(brie._angle_array) * np.deg2rad(step)
        conv = np.convolve(pdf, brie._coast_qs, mode="full")
        npad = len(brie._coast_qs) - 1
        first = npad - npad // 2
        return -conv[first:first + len(pdf)]

    @staticmethod
    def _diffusion_number(brie, x):
        """BRIE's diffusion number r = D dt / (2 dy^2) at each cell's forward face, clipped at zero.

        The angle index is clipped to the last entry of ``_coast_diff`` (``_wave_climl - 1``).
        BRIE's own step (brie.py:1297) clips one past it, which fails at faces steeper than ~89.5 deg.
        """
        theta = 180.0 * np.arctan2(np.roll(x, -1) - x, brie._dy) / np.pi
        k = np.clip(np.round(90 - theta).astype(int), 1, brie._wave_climl - 1)
        return np.maximum(0.0, brie._coast_diff[k] * brie._dt / 2 / brie._dy ** 2)

    @classmethod
    def _leak(cls, brie, x):
        """The leak term 2 (r_i - r_{i-1}) (x_i - x_{i-1}) for shoreline x, from BRIE's diffusivity."""
        step = np.roll(x, -1) - x
        r = cls._diffusion_number(brie, x)
        return 2 * (r - np.roll(r, 1)) * np.roll(step, 1)

    @classmethod
    def _brie_step(cls, brie, x, dx):
        """BRIE's implicit alongshore step (brie.py:1293-1326) for shoreline x and change dx."""
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import spsolve
        r = cls._diffusion_number(brie, x)
        dv = np.r_[-r[-1], -r[1:], 1 + 2 * r, -r[0:-1], -r[0]]
        rhs = x + r * (np.roll(x, -1) - 2 * x + np.roll(x, 1)) + dx
        return spsolve(csr_matrix((dv, (brie._di, brie._dj))), rhs)

    def _check_zone(self, lo, hi, half_width, n_domains):
        """Raise unless the correction zone fits inside the array: it does not wrap round the ends."""
        if lo + 1 - half_width < 0 or hi + half_width - 1 > n_domains - 1:
            raise GroinOutsideTestedRange(
                f"the correction zone ({half_width} cells either side of the face {lo}|{hi}) runs past "
                f"the end of the {n_domains} domains; it does not wrap, so the groin would leak sand there.")

    def _check_tested_range(self, brie, x_s, near):
        """Stop the run once it leaves the conditions the module was tested in (6-stress/README.md).

        Three checks, each year the groin is installed: BRIE's wave climate is the one the run
        started with (the drift table is read once); the groin's footprint (this shoreline less
        the groin-free one) sums to within SAND_TOLERANCE_M or SAND_TOLERANCE_FRACTION of its
        gross, so the groin is not making or losing sand; and the groin face is no steeper than
        MAX_FACE_ANGLE_DEG, past which BRIE's diffusion stops carrying sand back across it.
        """
        if self._coast_diff0 is None:
            self._coast_diff0 = np.array(brie._coast_diff, dtype=float, copy=True)
        elif not np.array_equal(brie._coast_diff, self._coast_diff0):
            raise GroinOutsideTestedRange(
                "the wave climate changed during the run; the groin module is tested only with one "
                "wave climate for the whole run (hard-structures/groin/3-module-tests/6-stress/).")
        net = float(np.sum(self._a))
        gross = float(np.sum(np.abs(self._a)))
        if abs(net) > max(self.SAND_TOLERANCE_M, self.SAND_TOLERANCE_FRACTION * gross):
            raise GroinOutsideTestedRange(
                f"the groin has {'made' if net < 0 else 'lost'} {abs(net):.1f} m of summed shoreline "
                f"(footprint gross {gross:.1f} m), more than the module's tolerance; BRIE's diffusion at the "
                "groin is outside the range the module was tested in (hard-structures/groin/3-module-tests/6-stress/).")
        face = abs(180.0 * np.arctan2(x_s[self._hi] - x_s[self._lo], brie._dy) / np.pi)
        if face > self.MAX_FACE_ANGLE_DEG:
            raise GroinOutsideTestedRange(
                f"the groin face is at {face:.1f} deg, steeper than {self.MAX_FACE_ANGLE_DEG:g} deg; past this "
                "BRIE's diffusion no longer carries sand back across the groin and the fillet grows without "
                "limit (hard-structures/groin/3-module-tests/6-stress/).")

    def __call__(self, cascade, x_s_dt):
        """Move a share of the arriving drift across the face, and cancel BRIE's leak of it.

        Parameters
        ----------
        cascade : Cascade
            The calling CASCADE instance; BRIE is read from ``cascade._brie_coupler._brie``.
        x_s_dt : numpy array of float
            Per-domain shoreline change (metres, landward +) before the alongshore
            solve. Modified in place and returned.
        """
        year = self.start_year + self._call_count
        self._call_count += 1
        active = year >= self.install_year
        s = self._effective_trapping_rate(year) if active else 0.0

        brie = cascade._brie_coupler._brie
        x_s = np.asarray(brie.x_s, dtype=float)
        ny, lo, hi = len(x_s), self._lo, self._hi
        if self._q_net is None:
            self._q_net = self._drift_table(brie)
            self._a = np.zeros(ny)
        if self._prev is not None:
            # Last year's step with and without the groin's own contribution; everything else cancels
            xp, dxp, gp, ap = self._prev
            self._a = self._brie_step(brie, xp, dxp) - self._brie_step(brie, xp - ap, dxp - gp)

        def angle(i):
            return 180.0 * np.arctan2(x_s[(i + 1) % ny] - x_s[i], brie._dy) / np.pi

        def q_net(theta):
            return float(self._q_net[int(np.clip(np.round(90 - theta), 1, brie._wave_climl - 1))])

        # The drift arriving at the face, from whichever side it comes
        toward_lower = q_net(angle(hi)) < 0
        approach = angle(hi) if toward_lower else angle(lo - 1)
        qv = q_net(approach)
        q = s * qv * brie._dt / (brie._dy * (brie._h_b_crit + brie._d_sf))
        dx_lo, dx_hi = (-q, q) if (active and s > 0.0) else (0.0, 0.0)
        x_s_dt[lo] += dx_lo
        x_s_dt[hi] += dx_hi

        # The leak correction on the groin's own bends: this shoreline's term less the groin-free shoreline's
        x_none = x_s - self._a
        near = np.abs(np.arange(ny) - (lo + 0.5)) <= self.half_width
        corr = np.where(near, self._leak(brie, x_s) - self._leak(brie, x_none), 0.0)
        x_s_dt += corr
        if active:
            self._check_tested_range(brie, x_s, near)
        g = corr.copy()
        g[lo] += dx_lo
        g[hi] += dx_hi
        self._prev = (x_s.copy(), np.asarray(x_s_dt, dtype=float).copy(), g, self._a.copy())

        dx_up, dx_down = (dx_lo, dx_hi) if self.updrift_pad == lo else (dx_hi, dx_lo)
        self.year_TS.append(year)
        self.active_TS.append(bool(active))
        self.blocking_applied_TS.append(s)
        self.applied_dx_updrift_TS.append(dx_up)
        self.applied_dx_downdrift_TS.append(dx_down)
        self.trapping_rate_applied_TS.append(abs(dx_up))
        self.face_offset_m_TS.append(float(x_s[hi] - x_s[lo]))
        self.net_drift_m3_yr_TS.append(qv)
        self.approach_angle_deg_TS.append(approach)
        self.leak_correction_m_TS.append(float(corr.sum()))
        self.footprint_max_m_TS.append(float(np.abs(self._a).max()))
        return x_s_dt

    @property
    def mean_trapping_rate_m_yr(self) -> float:
        """Mean shoreline change moved across the face over the active years, m/yr."""
        active = [m for m, a in zip(self.trapping_rate_applied_TS, self.active_TS) if a]
        return float(np.mean(active)) if active else float("nan")

    def summary(self) -> Dict[str, object]:
        """Return a compact run-metadata dictionary for logging."""
        return dict(
            kind=self.kind,
            updrift_pad=self.updrift_pad,
            downdrift_pad=self.downdrift_pad,
            strength_history=[list(p) for p in self.strength_history],
            install_year=self.install_year,
            start_year=self.start_year,
            half_width=self.half_width,
            years_active=int(np.sum(self.active_TS)),
            mean_trapping_rate_m_yr=self.mean_trapping_rate_m_yr,
            cumulative_updrift_m=float(np.sum(self.applied_dx_updrift_TS)),
            cumulative_downdrift_m=float(np.sum(self.applied_dx_downdrift_TS)),
        )

    def diagnostics_frame(self) -> List[Dict[str, object]]:
        """Return the per-year record: the strength, the shoreline change moved, the drift, the correction."""
        cum_up = np.cumsum(self.applied_dx_updrift_TS) if self.applied_dx_updrift_TS else []
        cum_down = np.cumsum(self.applied_dx_downdrift_TS) if self.applied_dx_downdrift_TS else []
        rows: List[Dict[str, object]] = []
        for i, year in enumerate(self.year_TS):
            rows.append(dict(
                model_year=year,
                groin_active=self.active_TS[i],
                strength_applied=self.blocking_applied_TS[i],
                trapping_rate_applied_m_yr=self.trapping_rate_applied_TS[i],
                applied_dx_updrift_m=self.applied_dx_updrift_TS[i],
                applied_dx_downdrift_m=self.applied_dx_downdrift_TS[i],
                cumulative_updrift_m=float(cum_up[i]),
                cumulative_downdrift_m=float(cum_down[i]),
                face_offset_m=self.face_offset_m_TS[i],
                net_drift_m3_yr=self.net_drift_m3_yr_TS[i],
                approach_angle_deg=self.approach_angle_deg_TS[i],
                leak_correction_m=self.leak_correction_m_TS[i],
                footprint_max_m=self.footprint_max_m_TS[i],
            ))
        return rows
