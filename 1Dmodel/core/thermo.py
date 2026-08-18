"""Fluid-agnostic thermodynamic state backends — Stage A of the FV core rework.

See docs/solver_design/FV_CORE_REWORK_PLAN.md section 3.2. This module is the
canonical home for ``ThermoState`` (the former ``CoolantState`` in
``physics/liquid_flow/dispatch.py``) and the (p,h)/(T,p) state-construction
functions that back it. ``physics/liquid_flow/dispatch.py`` re-exports these
names unchanged so every existing caller (main_solve.py,
main_solve_shellntube.py, hx_adapters.py, sanity_checks.py, transient_core/*)
keeps working without modification — this is a pure relocation, not a
rewrite: the CoolProp call sequences below are byte-for-byte the calls that
used to live in dispatch.py, and are what makes the Stage A bit-identical
acceptance gate achievable.

Three backends implement the same ``ThermoBackend`` protocol:

- ``RealFluidBackend`` — (p,h)-primary real-fluid CoolProp state, covering
  single-phase, two-phase (dome), and supercritical for any CoolProp fluid.
  This is what ``coolant_model="equilibrium_liquid"`` already uses; it is now
  also the generic fluid-agnostic path for GHe/N2/water/etc.
- ``IdealGasBackend`` — the legacy fast (T,p) path used by
  ``coolant_model="single_phase_coolprop"`` today (raw ``CP.PropsSI`` calls,
  no cached ``AbstractState``). Kept only so the existing helium bit-identical
  regression can be reproduced during migration; see the "Known gaps" note in
  the design doc about deleting it once the migration is complete.
- ``ReactingGasBackend`` — thin adapter around the existing FPV/equilibrium/
  frozen manifold "gas state provider" callables already used by
  ``transient_core/adapters_shelltube.py`` and the FPV calls in
  ``main_solve.py``. Combustion-product state is keyed on
  ``(h_removed, progress_variable, p)``, not ``(p, h)`` or ``(T, p)`` — it
  intentionally does NOT implement the same ``state_ph``/``state_pT`` protocol
  as the other two backends. Forcing that shape would hide the real
  difference (no per-node Cantera calls; manifold lookup only) behind a
  misleading interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Protocol

import CoolProp.CoolProp as CP

from hps_combustor.physics.liquid_flow.coolprop_state_cache import (
    coolprop_fluid_string,
    get_cached_state,
)
from hps_combustor.physics.liquid_flow.correlations import saturation_state
from hps_combustor.physics.liquid_flow import regime as _regime


@dataclass(frozen=True)
class ThermoState:
    """Fluid-agnostic single-node thermodynamic state.

    Formerly ``CoolantState`` in ``physics/liquid_flow/dispatch.py`` — moved
    here verbatim (same fields, same semantics) as the Stage A relocation.
    ``dispatch.py`` keeps ``CoolantState`` as an alias of this class so
    existing ``isinstance``/field-access callers are unaffected.
    """

    fluid: str
    model: str
    p_Pa: float
    T_K: float
    h_J_kg: float
    rho_kg_m3: float
    mu_Pa_s: float
    k_W_m_K: float
    cp_J_kg_K: float
    Pr: float
    quality: float
    void_fraction: float
    phase: str
    # Supercritical-regime fields (defaulted so the subcritical constructors
    # above are untouched). ``is_supercritical`` gates the supercritical
    # closure branch; ``p_reduced`` = p/p_crit; ``T_pc_K`` is the
    # pseudo-critical temperature at this pressure (None subcritically).
    is_supercritical: bool = False
    p_reduced: float | None = None
    T_pc_K: float | None = None


def coolant_state_from_Tp(fluid: str, T_K: float, p_Pa: float) -> ThermoState:
    """Single-phase CoolProp state from temperature and pressure.

    Relocated verbatim from ``physics/liquid_flow/dispatch.py`` — same five
    independent ``PropsSI`` calls, same order. Each high-level ``PropsSI``
    call is a pure function of its arguments (no cross-call state), so the
    call order does not affect the numeric result; it is preserved anyway to
    keep this a mechanical move rather than a rewrite.
    """
    rho = CP.PropsSI("D", "T", T_K, "P", p_Pa, fluid)
    mu = CP.PropsSI("V", "T", T_K, "P", p_Pa, fluid)
    k = CP.PropsSI("L", "T", T_K, "P", p_Pa, fluid)
    cp = CP.PropsSI("C", "T", T_K, "P", p_Pa, fluid)
    h = CP.PropsSI("H", "T", T_K, "P", p_Pa, fluid)
    return ThermoState(
        fluid=fluid,
        model="single_phase_coolprop",
        p_Pa=float(p_Pa),
        T_K=float(T_K),
        h_J_kg=float(h),
        rho_kg_m3=float(rho),
        mu_Pa_s=float(mu),
        k_W_m_K=float(k),
        cp_J_kg_K=float(cp),
        Pr=float(cp * mu / k),
        quality=float("nan"),
        void_fraction=0.0,
        phase="single_phase",
    )


def coolant_state_from_ph(
    fluid: str, p_Pa: float, h_J_kg: float, model: str, backend: str = "HEOS"
) -> ThermoState:
    """Coolant/fluid state closure from pressure and enthalpy.

    Relocated verbatim from ``physics/liquid_flow/dispatch.py``. ``backend``
    opts into a faster, interpolated CoolProp property backend ("TTSE"/
    "BICUBIC" — see coolprop_state_cache.py) for the ``equilibrium_liquid``
    path only; the returned ``ThermoState.fluid`` is always the plain fluid
    name regardless (tagging stays an internal computation detail).
    """
    if model == "single_phase_coolprop":
        T = CP.PropsSI("T", "P", p_Pa, "H", h_J_kg, fluid)
        return coolant_state_from_Tp(fluid, T, p_Pa)
    if model != "equilibrium_liquid":
        raise ValueError(f"unknown coolant model: {model!r}")

    fluid_cp = coolprop_fluid_string(fluid, backend)
    # real_fluid_state_ph: dome-based below p_crit (bit-identical to the
    # former equilibrium_state_ph call), single-phase real-EOS above it (no
    # crash).
    eq = _regime.real_fluid_state_ph(fluid_cp, p_Pa, h_J_kg)
    supercritical = _regime.is_supercritical(fluid_cp, p_Pa)
    if eq.phase == "two_phase":
        sat = saturation_state(fluid_cp, p_Pa)
        x = min(max(eq.quality, 0.0), 1.0)
        mu = (1.0 - x) * sat.mu_l_Pa_s + x * sat.mu_v_Pa_s
        k = (1.0 - x) * sat.k_l_W_m_K + x * sat.k_v_W_m_K
        cp = (1.0 - x) * sat.cp_l_J_kg_K + x * sat.cp_v_J_kg_K
    else:
        flashed = get_cached_state(fluid_cp).flash_ph(p_Pa, h_J_kg)
        mu = flashed.viscosity()
        k = flashed.conductivity()
        cp = flashed.cpmass()
    p_crit = get_cached_state(fluid_cp).p_crit_Pa
    T_pc = _regime.pseudo_critical_temperature(fluid_cp, p_Pa) if supercritical else None
    return ThermoState(
        fluid=fluid,
        model=model,
        p_Pa=float(p_Pa),
        T_K=float(eq.T_K),
        h_J_kg=float(h_J_kg),
        rho_kg_m3=float(eq.rho_kg_m3),
        mu_Pa_s=float(mu),
        k_W_m_K=float(k),
        cp_J_kg_K=float(cp),
        Pr=float(cp * mu / k),
        quality=float(eq.quality),
        void_fraction=float(eq.void_fraction),
        phase=eq.phase,
        is_supercritical=bool(supercritical),
        p_reduced=float(p_Pa / p_crit),
        T_pc_K=None if T_pc is None else float(T_pc),
    )


def coolant_inlet_state(coolant_prop) -> ThermoState:
    """Build the inlet state from a ``coolantProp``-like dataclass instance.

    Relocated verbatim from ``physics/liquid_flow/dispatch.py``.
    """
    model = getattr(coolant_prop, "coolant_model", "single_phase_coolprop")
    fluid = getattr(coolant_prop, "coolant", "Helium")
    T_in = float(coolant_prop.T_in)
    p_in = float(coolant_prop.p_in)
    if model == "single_phase_coolprop":
        return coolant_state_from_Tp(fluid, T_in, p_in)
    backend = getattr(coolant_prop, "liquid_property_backend", "HEOS")
    h_in = CP.PropsSI("H", "T", T_in, "P", p_in, fluid)
    return coolant_state_from_ph(fluid, p_in, h_in, model, backend=backend)


class ThermoBackend(Protocol):
    """Common interface for pure-fluid property backends.

    Reacting-gas state is deliberately NOT part of this protocol — see
    ``ReactingGasBackend`` below.
    """

    def state_ph(self, fluid: str, p_Pa: float, h_J_kg: float) -> ThermoState: ...

    def state_pT(self, fluid: str, T_K: float, p_Pa: float) -> ThermoState: ...

    def p_crit(self, fluid: str) -> float: ...


class RealFluidBackend:
    """(p,h)-primary real-fluid backend — CoolProp HEOS/TTSE/BICUBIC.

    Covers single-phase, two-phase, and supercritical for any CoolProp fluid
    (GHe, N2/LN2, water/steam, O2, ...). This is the fluid-agnostic backend:
    no fluid name appears in its logic, only in the arguments it is called
    with.
    """

    def __init__(self, property_backend: str = "HEOS") -> None:
        self._property_backend = property_backend

    def state_ph(self, fluid: str, p_Pa: float, h_J_kg: float) -> ThermoState:
        return coolant_state_from_ph(
            fluid, p_Pa, h_J_kg, model="equilibrium_liquid", backend=self._property_backend
        )

    def state_pT(self, fluid: str, T_K: float, p_Pa: float) -> ThermoState:
        return coolant_state_from_Tp(fluid, T_K, p_Pa)

    def p_crit(self, fluid: str) -> float:
        fluid_cp = coolprop_fluid_string(fluid, self._property_backend)
        return get_cached_state(fluid_cp).p_crit_Pa


class IdealGasBackend:
    """Legacy (T,p) raw-``PropsSI`` fast path — ``single_phase_coolprop``.

    Exists to reproduce the existing helium bit-identical regression during
    migration (see docs/solver_design/FV_CORE_REWORK_PLAN.md Stage A/G). Also
    exposes single-property getters matching the exact per-property call
    granularity ``main_solve.py``/``main_solve_shellntube.py`` use inline
    today (e.g. only density at one point in the march, only enthalpy at
    another) — swapping a call site to one of these getters changes nothing
    about which CoolProp evaluations happen or with what arguments.
    """

    def density(self, fluid: str, T_K: float, p_Pa: float) -> float:
        return float(CP.PropsSI("D", "T", T_K, "P", p_Pa, fluid))

    def enthalpy(self, fluid: str, T_K: float, p_Pa: float) -> float:
        return float(CP.PropsSI("H", "T", T_K, "P", p_Pa, fluid))

    def viscosity(self, fluid: str, T_K: float, p_Pa: float) -> float:
        return float(CP.PropsSI("V", "T", T_K, "P", p_Pa, fluid))

    def conductivity(self, fluid: str, T_K: float, p_Pa: float) -> float:
        return float(CP.PropsSI("L", "T", T_K, "P", p_Pa, fluid))

    def cp(self, fluid: str, T_K: float, p_Pa: float) -> float:
        return float(CP.PropsSI("C", "T", T_K, "P", p_Pa, fluid))

    def cv(self, fluid: str, T_K: float, p_Pa: float) -> float:
        return float(CP.PropsSI("CVMASS", "T", T_K, "P", p_Pa, fluid))

    def compressibility(self, fluid: str, T_K: float, p_Pa: float) -> float:
        return float(CP.PropsSI("Z", "T", T_K, "P", p_Pa, fluid))

    def molar_mass(self, fluid: str) -> float:
        return float(CP.PropsSI("MOLAR_MASS", fluid))

    def state_ph(self, fluid: str, p_Pa: float, h_J_kg: float) -> ThermoState:
        return coolant_state_from_ph(fluid, p_Pa, h_J_kg, model="single_phase_coolprop")

    def state_pT(self, fluid: str, T_K: float, p_Pa: float) -> ThermoState:
        return coolant_state_from_Tp(fluid, T_K, p_Pa)

    def p_crit(self, fluid: str) -> float:
        return float(CP.PropsSI("PCRIT", fluid))


# ``ShellTubeGasState``-shaped dict/object from a
# transient_core.adapters_shelltube-style provider: callable(h_removed,
# progress_variable, node_index) -> gas state.
GasStateProvider = Callable[[float, float, int], object]


class ReactingGasBackend:
    """Thin adapter around an existing FPV/equilibrium/frozen gas-state
    provider (see ``transient_core/adapters_shelltube.py``:
    ``fpv_gas_state_provider``, ``equilibrium_gas_state_provider``,
    ``oxygen_gas_state_provider``).

    Deliberately NOT a ``ThermoBackend``: combustion-product state is keyed
    on ``(h_removed, progress_variable, node_index)`` from a tabulated
    manifold, never on ``(p, h)`` via a per-node Cantera call (forbidden on
    the march hot path — see CLAUDE.md / docs/context/TRANSIENT_STATUS.md).
    This class exists so ``core/residual.py`` can hold one object per stream
    regardless of whether that stream is a pure fluid or a reacting gas,
    without pretending the two have the same state parameterization.
    """

    def __init__(self, provider: GasStateProvider) -> None:
        self._provider = provider

    def state(self, h_removed: float, progress_variable: float, node_index: int = 0):
        return self._provider(h_removed, progress_variable, node_index)
