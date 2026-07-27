"""
CoolingPFR: regime-2 companion to Flamelet for a diesel/O2 (or similar)
heat exchanger whose flow has TWO physically distinct zones:

  Regime 1 (flamelet.py, `Flamelet`):
      The mixing/reaction zone where the two feed streams (fuel, oxidizer)
      first meet and burn -- a non-premixed diffusion flame, correctly
      modeled in mixture-fraction space (this is what flamelet.py is for).

  Regime 2 (this module, `CoolingPFR`):
      Downstream of the flame, the two streams are gone -- there is now ONE
      already-mixed (mostly burnt) hot-gas stream flowing through the
      heat-exchanger channel, cooling as it gives up energy to the walls.
      There is no second stream and no meaningful mixture fraction here, so
      the flamelet-in-Z abstraction does NOT apply. The correct cheap
      finite-rate idiom for a single well-mixed stream advancing along a
      channel is a plug-flow reactor (PFR): a 0-D constant-pressure reactor
      marched along the 1-D axial coordinate, with an imposed wall heat-loss
      term.

Coupling between the two regimes: take the flamelet's burnt/product state
at (or near) its stoichiometric node -- `Flamelet.T_at_Z(Z_st)` /
`Flamelet.Y_at_Z(Z_st, ...)`, or the full channel-mean-Z composition if the
channel mixes flamelet output further -- as the PFR's inlet (T_in, Y_in).
See example_cooling.py for a worked hand-off.

Physics
-------
Per unit length, constant mass flow mdot [kg/s] through a channel of
diameter D (cross-section area A = pi D^2/4, wetted perimeter P = pi D):

    mdot * dh/dx = -q_wall(x),   q_wall(x) = h_conv * P * (T(x) - T_wall(x))

where h is the gas specific enthalpy (finite-rate chemistry changes h at
fixed p only through reaction is NOT what's assumed here -- reaction is
adiabatic at each infinitesimal step, so the ONLY thing that changes h
between the reactor's chemistry sub-step and the wall sub-step is the
imposed heat loss). This is applied via a Strang-like operator split per
axial step dx, exactly mirroring the flamelet module's split style:

    1. Chemistry sub-step: advance the constant-pressure reactor
       (Cantera IdealGasConstPressureReactor + ReactorNet, adiabatic) over
       dt = dx / u(x), u(x) = mdot / (rho(x) * A).
    2. Wall heat-loss sub-step: remove Q_seg = q_wall(x) * dx [W] from the
       gas by setting its specific enthalpy directly
       (`gas.HPY = h - Q_seg/mdot, p, Y`), which is exact (no operator-split
       error beyond the reactor's own chemistry substep) and makes the
       energy balance close to reactor tolerance:

           mdot * (h_in - h_out) == sum(Q_seg)  (== Q_wall_total)

This module depends on ONLY numpy and cantera (imports `_resolve_mechanism`
from .flamelet, an intra-package helper -- no hybrid_rocket import).
"""
from __future__ import annotations

from typing import Callable, Optional, Union

import numpy as np

try:
    from .flamelet import _resolve_mechanism
    from .steady_cache import SteadyCache
except ImportError:  # executed as a top-level script, not as part of the package
    from flamelet import _resolve_mechanism
    from steady_cache import SteadyCache


class CoolingPFR:
    """
    Plug-flow reactor with imposed wall heat loss, for the downstream
    already-mixed hot-gas leg of a two-regime (flamelet + cooling channel)
    heat exchanger.

    Parameters
    ----------
    mechanism : str
        Cantera mechanism YAML (same mechanism the upstream Flamelet used,
        so species indices/names line up for the hand-off).
    cvode_rtol, cvode_atol : float
        Cantera ReactorNet tolerances for the per-segment chemistry sub-step.
    """

    def __init__(
        self,
        mechanism: str,
        cvode_rtol: float = 1.0e-9,
        cvode_atol: float = 1.0e-15,
    ):
        try:
            import cantera as ct
        except ImportError:
            raise ImportError("Cantera is required for CoolingPFR.")
        self._ct = ct
        mech_path = _resolve_mechanism(mechanism)
        self._gas = ct.Solution(mech_path)
        self.species_names = list(self._gas.species_names)
        self._cvode_rtol = float(cvode_rtol)
        self._cvode_atol = float(cvode_atol)
        self._last_result: Optional[dict] = None

    def march(
        self,
        T_in: float,
        Y_in: np.ndarray,
        p: float,
        mdot: float,
        length: float,
        n_steps: int,
        diameter: float,
        h_conv: float,
        T_wall: Union[float, Callable[[float], float]],
    ) -> dict:
        """
        March finite-rate chemistry + wall heat loss along a channel of the
        given `length` [m], in `n_steps` equal axial segments.

        Parameters
        ----------
        T_in, Y_in, p : inlet gas state (K, mass-fraction array, Pa).
        mdot : mass flow rate [kg/s] (assumed constant along the channel).
        diameter : channel hydraulic diameter [m] (sets cross-section area
            and wetted perimeter for a circular channel; for a non-circular
            channel pass the hydraulic diameter and interpret `h_conv` as
            already referenced to that perimeter).
        h_conv : convective heat-transfer coefficient to the wall [W/m^2/K].
        T_wall : wall temperature [K], constant or callable(x) -> T_wall(x).

        Returns
        -------
        dict with keys:
            x (n_steps+1,), T (n_steps+1,), Y (n_steps+1, n_species),
            Q_wall_total_W, h_in, h_out,
            energy_balance_residual_W = mdot*(h_in-h_out) - Q_wall_total
            (should be ~0 to reactor tolerance; see tests/test_flamelet_kit.py).
        """
        ct = self._ct
        Y_in = np.asarray(Y_in, dtype=float)
        gas = self._gas
        gas.TPY = float(T_in), float(p), {
            sp: float(y) for sp, y in zip(self.species_names, Y_in) if y > 0
        }
        r = ct.IdealGasConstPressureReactor(gas, clone=True, energy="on")
        net = ct.ReactorNet([r])
        net.rtol = self._cvode_rtol
        net.atol = self._cvode_atol
        try:
            ph = r.phase
        except AttributeError:
            ph = r.thermo

        dx = float(length) / int(n_steps)
        cross_area = np.pi * diameter ** 2 / 4.0
        perimeter = np.pi * diameter

        n_sp = len(self.species_names)
        x_arr = np.zeros(n_steps + 1)
        T_arr = np.zeros(n_steps + 1)
        Y_arr = np.zeros((n_steps + 1, n_sp))
        T_arr[0] = T_in
        Y_arr[0] = Y_in

        h_in = float(ph.enthalpy_mass)
        Q_wall_total = 0.0
        x = 0.0
        for i in range(int(n_steps)):
            Tw = float(T_wall(x)) if callable(T_wall) else float(T_wall)

            u = mdot / max(float(ph.density) * cross_area, 1e-12)
            dt = dx / max(u, 1e-9)
            net.initial_time = 0.0
            net.advance(dt)

            T_chem = float(ph.T)
            q_wall = h_conv * perimeter * (T_chem - Tw)  # W/m
            Q_seg = q_wall * dx  # W

            h_before = float(ph.enthalpy_mass)
            h_after = h_before - Q_seg / max(mdot, 1e-12)
            Y_now = np.array(ph.Y, dtype=float)
            ph.HPY = h_after, float(p), Y_now
            r.syncState()

            Q_wall_total += Q_seg
            x += dx
            x_arr[i + 1] = x
            T_arr[i + 1] = float(ph.T)
            Y_arr[i + 1] = np.array(ph.Y, dtype=float)

        h_out = float(ph.enthalpy_mass)
        result = {
            "x": x_arr,
            "T": T_arr,
            "Y": Y_arr,
            "Q_wall_total_W": Q_wall_total,
            "h_in": h_in,
            "h_out": h_out,
            "mdot": float(mdot),
            "energy_balance_residual_W": mdot * (h_in - h_out) - Q_wall_total,
        }
        self._last_result = result
        return result

    def march_cached(
        self,
        cache: SteadyCache,
        T_in: float,
        Y_in: np.ndarray,
        p: float,
        mdot: float,
        **march_kwargs,
    ) -> dict:
        """
        Steady-cache-gated march: reuses `SteadyCache` exactly as flamelet.py
        does, re-purposing its three generic tolerance slots for the PFR's
        own operating inputs:

            SteadyCache "p"      slot <- channel pressure p
            SteadyCache "chi"    slot <- mass flow rate mdot
            SteadyCache "T_fuel" slot <- inlet temperature T_in

        (SteadyCache's field names are historical from the flamelet source;
        the tolerances themselves are just three independent relative/
        absolute bands, so reusing the class for a different physical
        triple is intentional -- it is the same cost lever: if the inlet
        state and flow rate have not moved outside tolerance since the last
        march, reuse the stored profile instead of re-marching.)

        Returns the (possibly cached) result dict from `march()`.
        """
        p_key = cache.key_p(p)
        chi_key = cache.key_chi(mdot)
        if cache.is_hit(p_key, chi_key, T_in) and self._last_result is not None:
            return self._last_result

        result = self.march(T_in, Y_in, p, mdot, **march_kwargs)
        dT = (0.0 if self._last_result is None
              else abs(result["T"][-1] - self._last_result["T"][-1]))
        cache.record_advance(p_key, chi_key, T_in, dT)
        return result
