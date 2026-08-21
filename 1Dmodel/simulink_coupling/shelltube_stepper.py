"""Decoupled, single-communication-step wrapper around the shell-and-tube
transient core, for embedding this HX as a component in an external
multiphysics loop (Simulink, an FMU master algorithm, a MATLAB System object,
a custom co-simulation bridge, ...).

This module does not modify anything outside `simulink_coupling/`. It:

- constructs the existing, validated `shellntube_transient_solver`
  (`..main_solve_shellntube_transient`) exactly as `main_transient.py` does, so
  geometry/materials/chemistry setup is 100% reused, not re-implemented;
- re-sequences the existing per-step physics primitives from
  `..transient_core.adapters_shelltube` (the same functions
  `_run_shelltube_transient_core_mass_energy` calls internally on every loop
  iteration) into a single external call, `step(dt, boundary)`, so an external
  caller can own time-stepping instead of this code owning a fixed schedule.

See README.md in this folder for the full I/O contract, and
SIMULINK_PLUGIN_GUIDE.md for how to wire this into Simulink via the FMU
wrapper in `fmu_wrapper.py`.

Known constraint inherited from the existing solver (not introduced here): the
compressible coolant path this stepper drives is Helium-only — several of the
private helpers reused below call CoolProp with a hardcoded `"Helium"` fluid
string. This stepper does not attempt to generalize that; see README.md.

**Standalone use**: this file uses only absolute `hps_combustor.*` imports
(never `from ..x import y`) specifically so it can be copied - alone, or as
part of this whole `simulink_coupling/` folder - to a machine that has never
seen the rest of this repository. If a real `hps_combustor` install isn't on
`sys.path`, the bootstrap below falls back to the vendored copy in
`_vendor/hps_combustor/` shipped alongside this file (see
`vendor_dependencies.py` and README.md's "Standalone deployment" section).
"""

from __future__ import annotations

import sys as _sys
from dataclasses import dataclass
from pathlib import Path as _Path

try:
    import hps_combustor as _hps_combustor  # noqa: F401
except ImportError:
    _vendor_dir = _Path(__file__).resolve().parent / "_vendor"
    if str(_vendor_dir) not in _sys.path:
        _sys.path.insert(0, str(_vendor_dir))
    import hps_combustor as _hps_combustor  # noqa: F401

import numpy as np
from CoolProp.CoolProp import PropsSI

from hps_combustor.main_solve_shellntube_transient import shellntube_transient_solver
from hps_combustor.transient_core import (
    build_shelltube_core_geometry_from_solver,
    coolprop_fluid_properties,
    coolprop_state_from_mass_energy,
    enforce_density_bounds,
    enforce_internal_energy_floor,
    equilibrium_gas_state_provider,
    fpv_gas_state_provider,
    initial_mass_energy_from_TP,
    oxygen_gas_state_provider,
    semi_implicit_wall_compressible_coolant_step,
    shelltube_step_inputs,
)
from hps_combustor.transient_core.adapters_shelltube import (
    _coolant_mass_energy_from_TP_profile,
    _coolprop_fluid_properties_at_profile,
    _limit_closed_valve_outlet_discharge_to_backpressure,
    _limit_face_mdot_for_inventory,
    _shelltube_boundary_pressure_profile,
    _shelltube_face_inertance,
    _shelltube_face_resistance,
    _shelltube_initial_pressure_profile,
    _shelltube_low_mach_lumped_faces,
    _shelltube_nominal_pressure_drop,
    _shelltube_quasi_steady_faces,
)

_FLUID = "Helium"  # matches the hardcoded fluid in the reused private helpers


@dataclass(frozen=True)
class BoundaryInputs:
    """Live boundary conditions supplied by the external solver for one step.

    See README.md for the full physical meaning of each field and which ones
    are load-bearing in which `coolant_momentum_model`.
    """

    mdot_coolant: float
    p_coolant_in: float
    p_coolant_out: float
    T_coolant_in: float
    mdot_hot_total: float
    ignited: bool = True
    T_lox_in: float = 300.0


@dataclass(frozen=True)
class StepOutputs:
    """Results returned from one `ShellTubeTransientStepper.step()` call."""

    T_wall: np.ndarray
    T_coolant: np.ndarray
    T_coolant_outlet: float
    p_coolant_outlet: float
    dp_coolant_hydraulic_Pa: float
    T_gas_outlet: float
    duty_W: float
    face_mdot: np.ndarray
    energy_residual_J: float
    mass_residual_kg: float


class ShellTubeTransientStepper:
    """Single-communication-step interface to the shell-and-tube transient core.

    Construct once per run with the same config dataclasses used by
    `main_transient.py` (`1Dmodel/input_data.py`), then call `step(dt,
    boundary)` once per external communication step. Internal state (wall
    temperature, coolant mass/energy, face mass flow memory) persists across
    calls automatically; use `get_state()`/`set_state()` for external
    checkpoint/rollback (e.g. an FMI master algorithm's `set_fmu_state`).
    """

    def __init__(
        self,
        coolantProp,
        hotgasProp,
        shellTubeProp,
        numericalProp,
        system_requirements,
        transientProp,
        corrCoeffs=None,
        N_axial=80,
        flow_config="co",
        mdot_coolant_reference=None,
        p_coolant_out_initial=None,
    ):
        solver = shellntube_transient_solver(
            coolantProp,
            hotgasProp,
            shellTubeProp,
            numericalProp,
            system_requirements,
            transientProp,
            corrCoeffs=corrCoeffs,
            N_axial=N_axial,
            flow_config=flow_config,
        )
        self._solver = solver
        tp = transientProp

        self._momentum_model = (
            "low_mach"
            if getattr(tp, "coolant_momentum_model", "quasi_steady") == "low_mach"
            else "quasi_steady"
        )
        self._mdot_floor = float(getattr(tp, "transient_coolant_mdot_floor", 1e-9))

        self.geometry = build_shelltube_core_geometry_from_solver(solver)
        self.bell_geometry = solver.geom
        self.corrCoeffs = solver.corrCoeffs
        self.wall_density = float(solver.rho_t)
        self.wall_conductivity_at_T = solver.k_t
        cp_t = solver.cp_t
        # Matches main_solve_shellntube_transient.py's solve_transient_core():
        # wall_cp is called with a Kelvin array and must return an array;
        # the material lookup (`cp_t`) itself expects degC and one value at a time.
        self.wall_cp = lambda T: np.array(
            [cp_t(float(Ti) - 273.15) for Ti in np.asarray(T)]
        )
        self.inside_tube_choice = solver.stp.inside_tube_choice
        self.nusselt_selector = solver.stp.Nusselt_tube
        self.tube_roughness = solver.stp.tube_roughness
        self.corrugation_thickness = solver.stp.corrugation_thickness
        self.corrugation_pitch = solver.stp.corrugation_pitch
        self.flow_direction = self.geometry.grid.flow_direction

        if solver._chem_mode == "finite_rate":
            self._combustion_provider, self._combustion_progress0 = fpv_gas_state_provider(
                solver._fpv
            )
        else:
            self._combustion_provider, self._combustion_progress0 = (
                equilibrium_gas_state_provider(solver._eqm)
            )

        self.coolant_properties_at = coolprop_fluid_properties(coolantProp.coolant)

        n = self.geometry.grid.n_cells
        T_wall0 = np.full(n, float(tp.T_wall_initial))
        T_coolant_initial = (
            float(coolantProp.T_in)
            if getattr(tp, "T_coolant_initial", None) is None
            else float(tp.T_coolant_initial)
        )
        T_coolant0 = np.full(n, T_coolant_initial)

        p_in0 = float(coolantProp.p_in)
        mdot0 = max(float(coolantProp.mass_flow_c), 0.0)
        # Hydraulic resistance/inertance is calibrated once, here, against a
        # single reference mass flow. The existing schedule-driven solver gets
        # to look ahead at the whole run's schedule and uses its maximum; a
        # live step-by-step caller (Simulink) cannot offer that look-ahead, so
        # this must be supplied as a design-time expectation instead (default:
        # the nominal `coolantProp.mass_flow_c`). If the live boundary trajectory
        # later commands flow well above this value, resistance/inertance stay
        # calibrated at the design point rather than the true peak - see README.
        mdot_ref_design = (
            float(mdot_coolant_reference) if mdot_coolant_reference is not None else mdot0
        )
        self._mdot_reference = max(mdot_ref_design, self._mdot_floor)
        self._dp_nominal = _shelltube_nominal_pressure_drop(
            self.geometry,
            self.bell_geometry,
            mdot_shell=self._mdot_reference,
            T_coolant=T_coolant0,
            p_coolant=p_in0,
            coolant_properties_at=self.coolant_properties_at,
            corrCoeffs=self.corrCoeffs,
            mdot_floor=self._mdot_floor,
        )
        # The initial pressure profile is an assumed cold-start estimate
        # (nominal design dp) unless the caller knows the true t=0 outlet
        # pressure - e.g. co-simulation starting from a known steady
        # operating point - in which case pass it via `p_coolant_out_initial`
        # for an exact match instead of an estimate. See README.md.
        p_out0 = (
            float(p_coolant_out_initial)
            if p_coolant_out_initial is not None
            else max(p_in0 - self._dp_nominal, 1.0e3)
        )
        if self._momentum_model == "low_mach":
            p_initial = _shelltube_boundary_pressure_profile(
                self.geometry,
                inlet_pressure=p_in0,
                outlet_pressure=p_out0,
                flow_direction=self.flow_direction,
            )
        else:
            p_initial = _shelltube_initial_pressure_profile(
                self.geometry,
                inlet_pressure=p_in0,
                pressure_drop=self._dp_nominal,
                flow_direction=self.flow_direction,
            )
        mass0, U0 = initial_mass_energy_from_TP(
            T_coolant0, p_initial, self.geometry.grid.coolant_volume, _FLUID
        )
        state0 = coolprop_state_from_mass_energy(
            mass0, U0, self.geometry.grid.coolant_volume, _FLUID
        )

        self.T_wall = T_wall0
        self.coolant_mass = mass0
        self.coolant_U = U0
        self.face_mdot = np.zeros(n + 1)
        self._resistance = _shelltube_face_resistance(
            self._dp_nominal, state0.density, self._mdot_reference, n_faces=n
        )
        self._inertance = _shelltube_face_inertance(self.geometry)

    # ------------------------------------------------------------------
    def get_state(self) -> dict:
        """Return a plain-array snapshot of continuous state for checkpointing."""

        return dict(
            T_wall=self.T_wall.copy(),
            coolant_mass=self.coolant_mass.copy(),
            coolant_U=self.coolant_U.copy(),
            face_mdot=self.face_mdot.copy(),
        )

    def set_state(self, state: dict) -> None:
        """Restore a snapshot previously returned by `get_state()`."""

        self.T_wall = np.array(state["T_wall"], dtype=float)
        self.coolant_mass = np.array(state["coolant_mass"], dtype=float)
        self.coolant_U = np.array(state["coolant_U"], dtype=float)
        self.face_mdot = np.array(state["face_mdot"], dtype=float)

    # ------------------------------------------------------------------
    def step(self, dt: float, boundary: BoundaryInputs) -> StepOutputs:
        """Advance internal state by `dt` seconds under `boundary` conditions.

        Mirrors, as a single call, exactly what one iteration of
        `_run_shelltube_transient_core_mass_energy`'s internal loop does
        (`1Dmodel/transient_core/adapters_shelltube.py`), except boundary
        values come directly from the caller instead of interpolated schedules.
        """

        grid = self.geometry.grid
        flow = self.flow_direction

        state = coolprop_state_from_mass_energy(
            self.coolant_mass, self.coolant_U, grid.coolant_volume, _FLUID
        )

        mdot_cmd = max(float(boundary.mdot_coolant), 0.0)
        p_inlet = float(boundary.p_coolant_in)
        p_outlet = max(float(boundary.p_coolant_out), 1.0e3)

        if self._momentum_model == "low_mach":
            p_transport = _shelltube_boundary_pressure_profile(
                self.geometry,
                inlet_pressure=p_inlet,
                outlet_pressure=p_outlet,
                flow_direction=flow,
            )
            faces = _shelltube_low_mach_lumped_faces(
                self.geometry,
                self.bell_geometry,
                self.face_mdot,
                state.temperature,
                p_transport,
                dt=dt,
                inlet_pressure=p_inlet,
                outlet_pressure=p_outlet,
                flow_direction=flow,
                mdot_reference=self._mdot_reference,
                corrCoeffs=self.corrCoeffs,
                mdot_floor=self._mdot_floor,
            )
            mdot_cap = max(mdot_cmd, self._mdot_floor)
        else:
            faces = _shelltube_quasi_steady_faces(
                state.pressure,
                state.density,
                self._resistance,
                mdot_inlet=mdot_cmd,
                outlet_pressure=p_outlet,
                flow_direction=flow,
                mdot_floor=self._mdot_floor,
            )
            mdot_cap = max(2.0 * self._mdot_reference, self._mdot_floor)
        faces = np.clip(faces, -mdot_cap, mdot_cap)

        if self._momentum_model != "low_mach" and mdot_cmd <= self._mdot_floor:
            outlet_face = -1 if flow == 1 else 0
            outlet_cell = grid.outlet_index
            outlet_drive = float(state.pressure[outlet_cell]) - float(p_outlet)
            if flow == -1:
                outlet_drive = float(p_outlet) - float(state.pressure[outlet_cell])
            if outlet_drive > 0.0 and abs(faces[outlet_face]) <= self._mdot_floor:
                faces[outlet_face] = 0.5 * self.face_mdot[outlet_face]
            faces = _limit_closed_valve_outlet_discharge_to_backpressure(
                self.geometry,
                state,
                self.coolant_mass,
                faces,
                dt=dt,
                outlet_pressure=p_outlet,
                flow_direction=flow,
            )
            faces = _limit_face_mdot_for_inventory(
                self.coolant_mass,
                faces,
                dt,
                internal_energy_J=self.coolant_U,
                specific_enthalpy_J_kg=state.specific_enthalpy_J_kg,
            )

        p_transport = _shelltube_boundary_pressure_profile(
            self.geometry,
            inlet_pressure=p_inlet,
            outlet_pressure=p_outlet,
            flow_direction=flow,
        )
        coolant_props_for_film = _coolprop_fluid_properties_at_profile(
            state.temperature, p_transport, _FLUID
        )

        if boundary.ignited:
            provider = self._combustion_provider
            progress_initial = self._combustion_progress0
        else:
            provider, progress_initial = oxygen_gas_state_provider(
                T_inlet=float(boundary.T_lox_in),
                pressure=max(float(self._solver.hotgasProp.p0), 1.0e4),
            )
        mdot_hot_t = max(float(boundary.mdot_hot_total), 0.0)
        T_inlet = float(boundary.T_coolant_in)

        assembled = shelltube_step_inputs(
            self.geometry,
            self.bell_geometry,
            Tbar_wall=self.T_wall,
            T_coolant=state.temperature,
            mdot_coolant=max(float(np.mean(np.abs(faces))), self._mdot_floor),
            T_coolant_inlet=T_inlet,
            p_coolant=p_inlet,
            mdot_hot_total=mdot_hot_t,
            gas_state_at=provider,
            coolant_properties_at=lambda _T, _p, props=coolant_props_for_film: props,
            wall_density=self.wall_density,
            wall_cp=self.wall_cp,
            wall_conductivity_at_T=self.wall_conductivity_at_T,
            inside_tube_choice=self.inside_tube_choice,
            nusselt_selector=self.nusselt_selector,
            tube_roughness=self.tube_roughness,
            corrCoeffs=self.corrCoeffs,
            corrugation_thickness=self.corrugation_thickness,
            corrugation_pitch=self.corrugation_pitch,
            progress_initial=progress_initial,
            flow_direction=flow,
        )

        h_inlet = float(PropsSI("H", "T", T_inlet, "P", max(p_inlet, 1.0e3), _FLUID))
        thermal_step = semi_implicit_wall_compressible_coolant_step(
            self.T_wall,
            assembled.wall_heat_capacity_J_K,
            self.coolant_mass,
            self.coolant_U,
            state.temperature,
            state.specific_enthalpy_J_kg,
            faces,
            assembled.wall_coolant_inputs.hot_heat_W,
            assembled.wall_coolant_inputs.wall_to_coolant_conductance_W_per_K,
            dt,
            inlet_enthalpy_J_kg=h_inlet,
            outlet_backflow_enthalpy_J_kg=h_inlet,
            mass_floor=1.0e-12,
        )
        m_candidate, U_candidate = enforce_density_bounds(
            np.maximum(thermal_step.coolant.mass_new, 1.0e-12),
            thermal_step.coolant.internal_energy_new_J,
            grid.coolant_volume,
        )
        if self._momentum_model == "low_mach":
            provisional_U = enforce_internal_energy_floor(
                m_candidate, U_candidate, grid.coolant_volume, _FLUID, clip=False
            )
            provisional_state = coolprop_state_from_mass_energy(
                m_candidate, provisional_U, grid.coolant_volume, _FLUID
            )
            p_projected = _shelltube_boundary_pressure_profile(
                self.geometry,
                inlet_pressure=p_inlet,
                outlet_pressure=p_outlet,
                flow_direction=flow,
            )
            m_new, U_new = _coolant_mass_energy_from_TP_profile(
                provisional_state.temperature, p_projected, grid.coolant_volume, _FLUID
            )
        else:
            m_new = m_candidate
            U_new = enforce_internal_energy_floor(
                m_new, U_candidate, grid.coolant_volume, _FLUID, clip=False
            )
        new_state = coolprop_state_from_mass_energy(
            m_new, U_new, grid.coolant_volume, _FLUID
        )

        self.T_wall = thermal_step.T_wall_new
        self.coolant_mass = m_new
        self.coolant_U = U_new
        self.face_mdot = faces

        duty_W = float(thermal_step.hot_heat_added_J) / dt if dt > 0.0 else 0.0
        # Bell-Delaware whole-exchanger dp estimate from THIS step's actual mass
        # flux/properties - a genuine physics prediction, independent of
        # momentum_model. Same aggregation convention as
        # `_shelltube_nominal_pressure_drop`/the steady solver's dp_c_bar: mean
        # magnitude across cells (each cell's Bell-Delaware value is already a
        # whole-exchanger-scale estimate, not a per-cell segment to sum).
        # Contrast with `p_coolant_outlet`: that is a *state reconstruction*
        # that, in low_mach mode, is mathematically forced to equal
        # `boundary.p_coolant_out` every step (see README.md) - it is not an
        # independent prediction. This field is. For closing a loop with an
        # external system model, use this (or p_coolant_in - this), not
        # p_coolant_outlet.
        dp_coolant_hydraulic_Pa = float(
            np.nanmean(np.abs(assembled.shell_film.dp_shell_Pa))
        )
        return StepOutputs(
            T_wall=self.T_wall.copy(),
            T_coolant=new_state.temperature.copy(),
            T_coolant_outlet=float(new_state.temperature[grid.outlet_index]),
            p_coolant_outlet=float(new_state.pressure[grid.outlet_index]),
            dp_coolant_hydraulic_Pa=dp_coolant_hydraulic_Pa,
            T_gas_outlet=float(assembled.hot_gas_march.T_gas_outlet),
            duty_W=duty_W,
            face_mdot=faces.copy(),
            energy_residual_J=float(thermal_step.total_energy_residual_J),
            mass_residual_kg=float(thermal_step.coolant.mass_residual_kg),
        )
