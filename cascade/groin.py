
"""The wave-driven groin: holds back a share of BRIE's alongshore drift, without BRIE leaking it away.

A CASCADE callback (``cascade._groin_callback``), called once per model year just before BRIE's
alongshore solve; it only adds to ``x_s_dt``. Units are metres, ``x_s`` is landward positive, and
pads are padded array indices from the runner. How it works, its tested range and its Buxton fit:
hard-structures/groin/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""

from __future__ import annotations

from typing import Dict, List, NamedTuple

import numpy as np

__all__ = ["GroinCallback", "GroinOutsideTestedRange"]


class GroinOutsideTestedRange(ValueError):
    """The run has left the conditions the groin module was tested in; its results would not hold."""


# Last year's state, kept to work out this year's footprint (each an array over the padded domains)
class _LastYear(NamedTuple):
    x_s_m: np.ndarray          # shoreline at the start of the year
    dx_m: np.ndarray           # the year's shoreline change passed to BRIE, groin included
    groin_dx_m: np.ndarray     # the groin's own part of that change
    footprint_m: np.ndarray    # the footprint used that year


# Shoreline angle of each cell's forward face, deg
def _face_angles(brie, x_s_m):
    return 180.0 * np.arctan2(np.roll(x_s_m, -1) - x_s_m, brie._dy) / np.pi


# Index into BRIE's per-angle tables for a shoreline angle, clipped to the table
def _table_index(brie, angle_deg):
    return np.clip(np.round(90 - angle_deg).astype(int), 1, brie._wave_climl - 1)


class GroinCallback:
    """A groin that holds back a share of the drift arriving at its face, and cancels BRIE's leak of it.

    Each year it adds two parts to ``x_s_dt``: the drift groin (a share of BRIE's net drift, the
    strength, moved from the downdrift to the updrift pad) and the leak correction (BRIE's diffusive
    term at the groin's own bends, cancelled within ``half_width`` cells of the face). It stops the
    run with :class:`GroinOutsideTestedRange` outside its tested range.

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
        history = sorted((int(year), float(strength)) for year, strength in strength_history)
        if not history or any(not (0.0 <= strength <= 1.0) for _, strength in history):
            raise ValueError("strength_history needs (year, strength) pairs with strength in [0, 1].")
        self.updrift_pad: int = int(updrift_pad)
        self.downdrift_pad: int = int(downdrift_pad)
        self.strength_history = history
        self.install_year: int = history[0][0]
        self.start_year: int = int(start_year)
        self.n_domains: int = int(n_domains)
        self.half_width: int = int(half_width)
        self._lo_pad: int = min(self.updrift_pad, self.downdrift_pad)
        self._hi_pad: int = max(self.updrift_pad, self.downdrift_pad)
        self._check_zone(self._lo_pad, self._hi_pad, self.half_width, self.n_domains)

        self.record: List[Dict[str, object]] = []   # one row per model year
        self._call_count: int = 0
        self._coast_diff_at_start = None       # BRIE's diffusivity table in the first year
        self._drift_table_m3_yr = None         # BRIE's net drift per shoreline-angle bin, read once
        self._footprint_m = None               # how far the groin has moved each domain's shoreline
        self._last_year = None                 # _LastYear, for the footprint

    def strength(self, year: int) -> float:
        """The strength in ``year``: the history's value, 0 before its first year."""
        value = 0.0
        for from_year, strength in self.strength_history:
            if year >= from_year:
                value = strength
        return value

    @staticmethod
    def _drift_table(brie):
        """BRIE's net alongshore drift for its wave climate, per shoreline-angle bin."""
        bin_width_deg = 180.0 / (brie._wave_climl - 1)
        angle_pdf = brie._angles.pdf(brie._angle_array) * np.deg2rad(bin_width_deg)
        convolved = np.convolve(angle_pdf, brie._coast_qs, mode="full")
        n_pad = len(brie._coast_qs) - 1
        first_bin = n_pad - n_pad // 2
        return -convolved[first_bin:first_bin + len(angle_pdf)]

    @staticmethod
    def _diffusion_number(brie, x_s_m):
        """BRIE's diffusion number r = D dt / (2 dy^2) at each cell's forward face, clipped at zero."""
        table_rows = _table_index(brie, _face_angles(brie, x_s_m))
        return np.maximum(0.0, brie._coast_diff[table_rows] * brie._dt / 2 / brie._dy ** 2)

    @classmethod
    def _leak(cls, brie, x_s_m):
        """The leak term 2 (r_i - r_{i-1}) (x_i - x_{i-1}) for shoreline x, from BRIE's diffusivity."""
        step_to_next_m = np.roll(x_s_m, -1) - x_s_m
        diffusion_number = cls._diffusion_number(brie, x_s_m)
        return 2 * (diffusion_number - np.roll(diffusion_number, 1)) * np.roll(step_to_next_m, 1)

    @classmethod
    def _brie_step(cls, brie, x_s_m, dx_m):
        """BRIE's implicit alongshore step (brie.py:1293-1326) for shoreline x_s_m and change dx_m."""
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import spsolve
        diffusion_number = cls._diffusion_number(brie, x_s_m)
        diagonals = np.r_[-diffusion_number[-1], -diffusion_number[1:], 1 + 2 * diffusion_number,
                          -diffusion_number[0:-1], -diffusion_number[0]]
        right_hand_side = (x_s_m + diffusion_number * (np.roll(x_s_m, -1) - 2 * x_s_m + np.roll(x_s_m, 1))
                           + dx_m)
        return spsolve(csr_matrix((diagonals, (brie._di, brie._dj))), right_hand_side)

    def _check_zone(self, lo_pad, hi_pad, half_width, n_domains):
        """Raise unless the pads are adjacent and the correction zone fits inside the array (it does not wrap)."""
        if hi_pad - lo_pad != 1:
            raise ValueError(
                f"groin domains must be adjacent (share the blocked boundary); got pads {lo_pad} and {hi_pad}.")
        if lo_pad + 1 - half_width < 0 or hi_pad + half_width - 1 > n_domains - 1:
            raise GroinOutsideTestedRange(
                f"the correction zone ({half_width} cells either side of the face {lo_pad}|{hi_pad}) runs past "
                f"the end of the {n_domains} domains; it does not wrap, so the groin would leak sand there.")

    def _check_tested_range(self, brie, x_s_m):
        """Raise GroinOutsideTestedRange if the wave climate changed, sand is unbalanced, or the face is too steep."""
        # Wave climate: the drift table is read once, so it must not change
        if self._coast_diff_at_start is None:
            self._coast_diff_at_start = np.array(brie._coast_diff, dtype=float, copy=True)
        elif not np.array_equal(brie._coast_diff, self._coast_diff_at_start):
            raise GroinOutsideTestedRange(
                "the wave climate changed during the run; the groin module is tested only with one "
                "wave climate for the whole run (hard-structures/groin/3-module-tests/6-stress/).")
        # Sand balance: the footprint should sum to about zero
        net_m = float(np.sum(self._footprint_m))
        gross_m = float(np.sum(np.abs(self._footprint_m)))
        if abs(net_m) > max(self.SAND_TOLERANCE_M, self.SAND_TOLERANCE_FRACTION * gross_m):
            raise GroinOutsideTestedRange(
                f"the groin has {'made' if net_m < 0 else 'lost'} {abs(net_m):.1f} m of summed shoreline "
                f"(footprint gross {gross_m:.1f} m), more than the module's tolerance; BRIE's diffusion at the "
                "groin is outside the range the module was tested in (hard-structures/groin/3-module-tests/6-stress/).")
        # Face angle: past this BRIE's diffusivity is clipped to zero
        face_angle_deg = abs(_face_angles(brie, x_s_m)[self._lo_pad])
        if face_angle_deg > self.MAX_FACE_ANGLE_DEG:
            raise GroinOutsideTestedRange(
                f"the groin face is at {face_angle_deg:.1f} deg, steeper than {self.MAX_FACE_ANGLE_DEG:g} deg; past "
                "this BRIE's diffusion no longer carries sand back across the groin and the fillet grows without "
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
        strength = self.strength(year)

        brie = cascade._brie_coupler._brie
        x_s_m = np.asarray(brie.x_s, dtype=float)
        n_cells, lo_pad, hi_pad = len(x_s_m), self._lo_pad, self._hi_pad
        if self._drift_table_m3_yr is None:
            self._drift_table_m3_yr = self._drift_table(brie)
        if self._footprint_m is None:
            self._footprint_m = np.zeros(n_cells)
        if self._last_year is not None:
            # Footprint: last year's BRIE step with and without the groin's own contribution
            last = self._last_year
            self._footprint_m = (self._brie_step(brie, last.x_s_m, last.dx_m)
                                 - self._brie_step(brie, last.x_s_m - last.footprint_m,
                                                   last.dx_m - last.groin_dx_m))

        # Drift groin: a share of the drift arriving at the face, from whichever side, as shoreline change
        face_angles_deg = _face_angles(brie, x_s_m)
        net_drift_m3_yr = self._drift_table_m3_yr[_table_index(brie, face_angles_deg)]
        arrival_face = hi_pad if net_drift_m3_yr[hi_pad] < 0 else lo_pad - 1
        arriving_drift_m3_yr = float(net_drift_m3_yr[arrival_face])
        trapped_m = (strength * arriving_drift_m3_yr * brie._dt
                     / (brie._dy * (brie._h_b_crit + brie._d_sf)))         # m3/yr to m of shoreline
        trap_lo_m, trap_hi_m = (-trapped_m, trapped_m) if strength > 0.0 else (0.0, 0.0)
        x_s_dt[lo_pad] += trap_lo_m
        x_s_dt[hi_pad] += trap_hi_m

        # Leak correction near the face: this shoreline's leak term less the groin-free shoreline's
        in_correction_zone = np.abs(np.arange(n_cells) - (lo_pad + 0.5)) <= self.half_width
        leak_correction_m = np.where(in_correction_zone,
                                     self._leak(brie, x_s_m) - self._leak(brie, x_s_m - self._footprint_m), 0.0)
        x_s_dt += leak_correction_m
        if active:
            self._check_tested_range(brie, x_s_m)

        # Keep this year's state and the groin's own contribution for next year's footprint
        groin_dx_m = leak_correction_m.copy()
        groin_dx_m[lo_pad] += trap_lo_m
        groin_dx_m[hi_pad] += trap_hi_m
        self._last_year = _LastYear(x_s_m.copy(), np.asarray(x_s_dt, dtype=float).copy(), groin_dx_m,
                                    self._footprint_m)

        # Record the year
        trap_updrift_m, trap_downdrift_m = ((trap_lo_m, trap_hi_m) if self.updrift_pad == lo_pad
                                            else (trap_hi_m, trap_lo_m))
        self.record.append(dict(
            model_year=year,
            groin_active=active,
            strength_applied=strength,
            applied_dx_updrift_m=trap_updrift_m,
            applied_dx_downdrift_m=trap_downdrift_m,
            face_offset_m=float(x_s_m[hi_pad] - x_s_m[lo_pad]),
            net_drift_m3_yr=arriving_drift_m3_yr,
            approach_angle_deg=float(face_angles_deg[arrival_face]),
            leak_correction_m=float(leak_correction_m.sum()),
            footprint_max_m=float(np.abs(self._footprint_m).max()),
        ))
        return x_s_dt

    @property
    def mean_trapping_rate_m_yr(self) -> float:
        """Mean shoreline change moved across the face over the active years, m/yr."""
        moved_m = [abs(row["applied_dx_updrift_m"]) for row in self.record if row["groin_active"]]
        return float(np.mean(moved_m)) if moved_m else float("nan")

    def summary(self) -> Dict[str, object]:
        """Return a compact run-metadata dictionary for logging."""
        return dict(
            kind=self.kind,
            updrift_pad=self.updrift_pad,
            downdrift_pad=self.downdrift_pad,
            strength_history=[list(pair) for pair in self.strength_history],
            install_year=self.install_year,
            start_year=self.start_year,
            half_width=self.half_width,
            years_active=sum(row["groin_active"] for row in self.record),
            mean_trapping_rate_m_yr=self.mean_trapping_rate_m_yr,
            cumulative_updrift_m=float(np.sum([row["applied_dx_updrift_m"] for row in self.record])),
            cumulative_downdrift_m=float(np.sum([row["applied_dx_downdrift_m"] for row in self.record])),
        )

    def diagnostics_frame(self) -> List[Dict[str, object]]:
        """Return the per-year record, with the shoreline change moved so far."""
        cumulative_updrift_m = np.cumsum([row["applied_dx_updrift_m"] for row in self.record])
        cumulative_downdrift_m = np.cumsum([row["applied_dx_downdrift_m"] for row in self.record])
        return [dict(model_year=row["model_year"],
                     groin_active=row["groin_active"],
                     strength_applied=row["strength_applied"],
                     trapping_rate_applied_m_yr=abs(row["applied_dx_updrift_m"]),
                     applied_dx_updrift_m=row["applied_dx_updrift_m"],
                     applied_dx_downdrift_m=row["applied_dx_downdrift_m"],
                     cumulative_updrift_m=float(cumulative_updrift_m[year_index]),
                     cumulative_downdrift_m=float(cumulative_downdrift_m[year_index]),
                     face_offset_m=row["face_offset_m"],
                     net_drift_m3_yr=row["net_drift_m3_yr"],
                     approach_angle_deg=row["approach_angle_deg"],
                     leak_correction_m=row["leak_correction_m"],
                     footprint_max_m=row["footprint_max_m"])
                for year_index, row in enumerate(self.record)]
