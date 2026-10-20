
"""The wave-driven groin: holds back a share of BRIE's alongshore drift, without BRIE leaking it away.

A CASCADE callback (``cascade._groin_callback``), called once per model year just before BRIE's
alongshore solve; it only adds to ``x_s_dt``. Units are metres, ``x_s`` is landward positive, and
pads are padded array indices from the runner. How it works, its tested range and its Buxton fit:
hard-structures/groin/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

__all__ = ["GroinCallback", "GroinOutsideTestedRange"]


class GroinOutsideTestedRange(ValueError):
    """The run has left the conditions the groin module was tested in; its results would not hold."""


# Shoreline angle of each cell's forward face, deg
def _face_angles(brie, x):
    return 180.0 * np.arctan2(np.roll(x, -1) - x, brie._dy) / np.pi


# Index into BRIE's per-angle tables for a shoreline angle, clipped to the table
def _table_index(brie, theta):
    return np.clip(np.round(90 - theta).astype(int), 1, brie._wave_climl - 1)


class GroinCallback:
    """A groin that holds back a share of the drift arriving at its face, and cancels BRIE's leak of it.

    Each year it adds two parts to ``x_s_dt``: the drift groin (a share ``s`` of BRIE's net drift
    moved from the downdrift to the updrift pad) and the leak correction (BRIE's diffusive term at
    the groin's own bends, cancelled within ``half_width`` cells of the face). It stops the run with
    :class:`GroinOutsideTestedRange` outside its tested range.

    Parameters
    ----------
    updrift_pad, downdrift_pad : int
        Padded indices of the two domains sharing the face. Adjacent.
    strength_history : sequence of (int, float)
        ``(from_year, strength)`` pairs; each strength, in [0, 1], holds until the next. The
        first year is the install year; strength is 0 before it.
    start_year : int
        Model year of the first call.
    n_domains : int
        Number of padded domains.
    half_width : int
        Cells either side of the face the leak correction acts on.
    """
    kind: str = "corrected"
    # The tested range: outside it the run stops (README "The module's bounds")
    SAND_TOLERANCE_M: float = 10.0         # summed shoreline the groin may make or lose, m ...
    SAND_TOLERANCE_FRACTION: float = 0.05  # ... or this share of its footprint's gross, whichever is larger
    MAX_FACE_ANGLE_DEG: float = 45.0       # steepest tested groin face, deg

    def __init__(
        self,
        updrift_pad: int,
        downdrift_pad: int,
        strength_history,
        start_year: int,
        n_domains: int,
        half_width: int = 10,
    ) -> None:
        history = sorted((int(y), float(s)) for y, s in strength_history)
        if not history or any(not (0.0 <= s <= 1.0) for _, s in history):
            raise ValueError("strength_history needs (year, strength) pairs with strength in [0, 1].")
        self.updrift_pad: int = int(updrift_pad)
        self.downdrift_pad: int = int(downdrift_pad)
        self.strength_history = history
        self.install_year: int = history[0][0]
        self.start_year: int = int(start_year)
        self.n_domains: int = int(n_domains)
        self.half_width: int = int(half_width)
        self._lo: int = min(self.updrift_pad, self.downdrift_pad)
        self._hi: int = max(self.updrift_pad, self.downdrift_pad)
        self._check_zone(self._lo, self._hi, self.half_width, self.n_domains)

        self.record: List[Dict[str, object]] = []   # one row per model year
        self._call_count: int = 0
        self._coast_diff0 = None
        self._q_net = None
        self._a = None
        self._prev = None

    def strength(self, year: int) -> float:
        """The strength in ``year``: the history's value, 0 before its first year."""
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
        """BRIE's diffusion number r = D dt / (2 dy^2) at each cell's forward face, clipped at zero."""
        k = _table_index(brie, _face_angles(brie, x))
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
        """Raise unless the pads are adjacent and the correction zone fits inside the array (it does not wrap)."""
        if hi - lo != 1:
            raise ValueError(f"groin domains must be adjacent (share the blocked boundary); got pads {lo} and {hi}.")
        if lo + 1 - half_width < 0 or hi + half_width - 1 > n_domains - 1:
            raise GroinOutsideTestedRange(
                f"the correction zone ({half_width} cells either side of the face {lo}|{hi}) runs past "
                f"the end of the {n_domains} domains; it does not wrap, so the groin would leak sand there.")

    def _check_tested_range(self, brie, x_s):
        """Raise GroinOutsideTestedRange if the wave climate changed, sand is unbalanced, or the face is too steep."""
        # Wave climate: the drift table is read once, so it must not change
        if self._coast_diff0 is None:
            self._coast_diff0 = np.array(brie._coast_diff, dtype=float, copy=True)
        elif not np.array_equal(brie._coast_diff, self._coast_diff0):
            raise GroinOutsideTestedRange(
                "the wave climate changed during the run; the groin module is tested only with one "
                "wave climate for the whole run (hard-structures/groin/3-module-tests/6-stress/).")
        # Sand balance: the footprint should sum to about zero
        net = float(np.sum(self._a))
        gross = float(np.sum(np.abs(self._a)))
        if abs(net) > max(self.SAND_TOLERANCE_M, self.SAND_TOLERANCE_FRACTION * gross):
            raise GroinOutsideTestedRange(
                f"the groin has {'made' if net < 0 else 'lost'} {abs(net):.1f} m of summed shoreline "
                f"(footprint gross {gross:.1f} m), more than the module's tolerance; BRIE's diffusion at the "
                "groin is outside the range the module was tested in (hard-structures/groin/3-module-tests/6-stress/).")
        # Face angle: past this BRIE's diffusivity is clipped to zero
        face = abs(_face_angles(brie, x_s)[self._lo])
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
        s = self.strength(year)

        brie = cascade._brie_coupler._brie
        x_s = np.asarray(brie.x_s, dtype=float)
        ny, lo, hi = len(x_s), self._lo, self._hi
        if self._q_net is None:
            self._q_net = self._drift_table(brie)
        if self._a is None:
            self._a = np.zeros(ny)
        if self._prev is not None:
            # Footprint: last year's BRIE step with and without the groin's own contribution
            xp, dxp, gp, ap = self._prev
            self._a = self._brie_step(brie, xp, dxp) - self._brie_step(brie, xp - ap, dxp - gp)

        # Drift groin: a share s of the drift arriving at the face, from whichever side, as shoreline change
        angles = _face_angles(brie, x_s)
        q_net = self._q_net[_table_index(brie, angles)]
        side = hi if q_net[hi] < 0 else lo - 1
        qv = float(q_net[side])
        q = s * qv * brie._dt / (brie._dy * (brie._h_b_crit + brie._d_sf))  # m3/yr to m of shoreline
        dx_lo, dx_hi = (-q, q) if s > 0.0 else (0.0, 0.0)
        x_s_dt[lo] += dx_lo
        x_s_dt[hi] += dx_hi

        # Leak correction near the face: this shoreline's leak term less the groin-free shoreline's
        near = np.abs(np.arange(ny) - (lo + 0.5)) <= self.half_width
        corr = np.where(near, self._leak(brie, x_s) - self._leak(brie, x_s - self._a), 0.0)
        x_s_dt += corr
        if active:
            self._check_tested_range(brie, x_s)

        # Keep this year's state and the groin's own contribution for next year's footprint
        g = corr.copy()
        g[lo] += dx_lo
        g[hi] += dx_hi
        self._prev = (x_s.copy(), np.asarray(x_s_dt, dtype=float).copy(), g, self._a)

        # Record the year
        dx_up, dx_down = (dx_lo, dx_hi) if self.updrift_pad == lo else (dx_hi, dx_lo)
        self.record.append(dict(
            model_year=year,
            groin_active=active,
            strength_applied=s,
            applied_dx_updrift_m=dx_up,
            applied_dx_downdrift_m=dx_down,
            face_offset_m=float(x_s[hi] - x_s[lo]),
            net_drift_m3_yr=qv,
            approach_angle_deg=float(angles[side]),
            leak_correction_m=float(corr.sum()),
            footprint_max_m=float(np.abs(self._a).max()),
        ))
        return x_s_dt

    @property
    def mean_trapping_rate_m_yr(self) -> float:
        """Mean shoreline change moved across the face over the active years, m/yr."""
        moved = [abs(r["applied_dx_updrift_m"]) for r in self.record if r["groin_active"]]
        return float(np.mean(moved)) if moved else float("nan")

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
            years_active=sum(r["groin_active"] for r in self.record),
            mean_trapping_rate_m_yr=self.mean_trapping_rate_m_yr,
            cumulative_updrift_m=float(np.sum([r["applied_dx_updrift_m"] for r in self.record])),
            cumulative_downdrift_m=float(np.sum([r["applied_dx_downdrift_m"] for r in self.record])),
        )

    def diagnostics_frame(self) -> List[Dict[str, object]]:
        """Return the per-year record, with the shoreline change moved so far."""
        cum_up = np.cumsum([r["applied_dx_updrift_m"] for r in self.record])
        cum_down = np.cumsum([r["applied_dx_downdrift_m"] for r in self.record])
        return [dict(model_year=r["model_year"],
                     groin_active=r["groin_active"],
                     strength_applied=r["strength_applied"],
                     trapping_rate_applied_m_yr=abs(r["applied_dx_updrift_m"]),
                     applied_dx_updrift_m=r["applied_dx_updrift_m"],
                     applied_dx_downdrift_m=r["applied_dx_downdrift_m"],
                     cumulative_updrift_m=float(cum_up[i]),
                     cumulative_downdrift_m=float(cum_down[i]),
                     face_offset_m=r["face_offset_m"],
                     net_drift_m3_yr=r["net_drift_m3_yr"],
                     approach_angle_deg=r["approach_angle_deg"],
                     leak_correction_m=r["leak_correction_m"],
                     footprint_max_m=r["footprint_max_m"])
                for i, r in enumerate(self.record)]
