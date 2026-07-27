"""
SteadyCache: generalized "is the stored flamelet solution still valid?" gate.

This is THE cost lever for coupling a flamelet solver into an outer flow
solver. Advancing a flamelet (a full Strang diffusion/chemistry step across
every Z-node) is the expensive operation; if the outer flow's operating
conditions (pressure, scalar dissipation, boundary stream states) have not
moved meaningfully since the last advance, the previously converged flamelet
solution is still a valid representative state and re-advancing is wasted
work that can simply be skipped.

Extracted and generalized from the rocket RIF solver's cache-gate logic
(src/hybrid_rocket/physics/chemistry/rif/manager.py, `_cache_gate_check` and
the `cache_tol_*` fields). Two gates from the source were dropped here on
purpose:

  - the "burning" precondition (T_max physical bound + is_burning +
    respark-margin check) -- that gate exists only because the rocket case
    has an ignition/extinction S-curve with a lower (quenched) branch that
    must never be cached as if it were the flame branch. A steady heat
    exchanger with no ignition/extinction has no such second branch, so
    there is nothing for this gate to protect against.

This generalized cache instead assumes: if a valid flamelet solution exists
at all (i.e. at least one advance has been recorded), and the operating
inputs have not drifted outside tolerance, it is safe to reuse.

Regime contrast (documented, not hypothetical): in the rocket's transient
combustion, chi_st jitters advance-to-advance (median ~6% relative, mean
~8%) which trips the ~5% `tol_chi` band on the majority of hot-burn
advances -- the cache hit rate collapsed to ~0.3% there. A steady or
steadily-ramping heat-exchanger flow has none of that acoustic/ignition
jitter; chi_st, p, and the boundary stream states move slowly and
monotonically (if at all) between calls, so the SAME tolerance bands that
starved the rocket cache are expected to produce a HIGH hit rate here --
this is the favorable regime for this technique.
"""
from __future__ import annotations

from typing import Optional


class SteadyCache:
    """
    Tolerance-gated reuse decision for a single flamelet's steady solution.

    Parameters
    ----------
    tol_dT : float
        Max allowed |T_max change| from the last advance (K). A large last-
        advance dT means the flamelet itself was still relaxing -- not yet
        steady -- so caching would freeze a transient, not a steady state.
    tol_p : float
        Max allowed relative pressure change since the last advance.
    tol_chi : float
        Max allowed relative chi_st change since the last advance.
    tol_T_fuel : float
        Max allowed absolute fuel (or whichever boundary you key on)
        temperature change since the last advance (K).
    p_ema_tau, chi_ema_tau : float
        Optional EMA time constants (same units as the `dt_elapsed` you pass
        to `key_p`/`key_chi`) for low-pass-filtering the cache KEYS only.
        0.0 (default) = use the instantaneous value as the key, unfiltered.
        This mirrors the source's `cache_key_filtered_p`/`cache_key_filtered_chi`
        opt-in: filtering only changes the skip/recompute DECISION, never the
        value the flamelet actually integrates with on a miss.
    """

    def __init__(
        self,
        tol_dT: float = 2.0,
        tol_p: float = 0.01,
        tol_chi: float = 0.05,
        tol_T_fuel: float = 5.0,
        p_ema_tau: float = 0.0,
        chi_ema_tau: float = 0.0,
    ):
        self.tol_dT = float(tol_dT)
        self.tol_p = float(tol_p)
        self.tol_chi = float(tol_chi)
        self.tol_T_fuel = float(tol_T_fuel)
        self.p_ema_tau = float(p_ema_tau)
        self.chi_ema_tau = float(chi_ema_tau)

        self._last: Optional[dict] = None  # {"p","chi","T_fuel","dT"}
        self._p_ema: Optional[float] = None
        self._chi_ema: Optional[float] = None

        self.n_hit = 0
        self.n_miss = 0
        self.miss_reasons: dict = {}

    # ------------------------------------------------------------------
    # Optional EMA-filtered cache keys
    # ------------------------------------------------------------------

    def key_p(self, p: float, dt_elapsed: float = 0.0) -> float:
        """Cache KEY for pressure -- instantaneous unless p_ema_tau > 0."""
        if self.p_ema_tau <= 0.0:
            return float(p)
        if self._p_ema is None or dt_elapsed <= 0.0:
            self._p_ema = float(p)
        else:
            import numpy as np
            a = 1.0 - float(np.exp(-dt_elapsed / self.p_ema_tau))
            self._p_ema = (1.0 - a) * self._p_ema + a * float(p)
        return float(self._p_ema)

    def key_chi(self, chi_st: float, dt_elapsed: float = 0.0) -> float:
        """Cache KEY for chi_st -- instantaneous unless chi_ema_tau > 0."""
        if self.chi_ema_tau <= 0.0:
            return float(chi_st)
        if self._chi_ema is None or dt_elapsed <= 0.0:
            self._chi_ema = float(chi_st)
        else:
            import numpy as np
            a = 1.0 - float(np.exp(-dt_elapsed / self.chi_ema_tau))
            self._chi_ema = (1.0 - a) * self._chi_ema + a * float(chi_st)
        return float(self._chi_ema)

    # ------------------------------------------------------------------
    # Gate evaluation
    # ------------------------------------------------------------------

    def check(self, p_key: float, chi_key: float, T_fuel: float) -> str:
        """
        Evaluate cache validity against the last recorded advance.

        Returns "" (HIT -- safe to skip re-advancing) or the name of the
        first failing gate: "no_prior_advance", "dT", "p", "chi", "T_fuel"
        (MISS -- must re-advance).
        """
        la = self._last
        if la is None:
            reason = "no_prior_advance"
        elif la["dT"] > self.tol_dT:
            reason = "dT"
        elif abs(p_key - la["p"]) > self.tol_p * max(la["p"], 1.0):
            reason = "p"
        elif abs(chi_key - la["chi"]) > self.tol_chi * max(abs(la["chi"]), 1e-12):
            reason = "chi"
        elif abs(T_fuel - la["T_fuel"]) > self.tol_T_fuel:
            reason = "T_fuel"
        else:
            reason = ""

        if reason == "":
            self.n_hit += 1
        else:
            self.n_miss += 1
            self.miss_reasons[reason] = self.miss_reasons.get(reason, 0) + 1
        return reason

    def is_hit(self, p_key: float, chi_key: float, T_fuel: float) -> bool:
        return self.check(p_key, chi_key, T_fuel) == ""

    def record_advance(self, p_key: float, chi_key: float, T_fuel: float,
                        dT: float) -> None:
        """Call after actually performing a real flamelet advance (a MISS was
        resolved by advancing), to update the state future `check()` calls
        compare against."""
        self._last = {"p": p_key, "chi": chi_key, "T_fuel": T_fuel, "dT": dT}

    @property
    def hit_rate(self) -> float:
        n = self.n_hit + self.n_miss
        return self.n_hit / n if n > 0 else 0.0
