"""FMI 2.0 Co-Simulation wrapper around `ShellTubeTransientStepper`.

Requires `pythonfmu` (`pip install pythonfmu`) at *build* time on the machine
that packages the FMU. See SIMULINK_PLUGIN_GUIDE.md for the build command and
how to import the resulting `.fmu` into Simulink.

Design notes:

- All FMI variables are plain floats/bools (FMI 2.0 has no native array type),
  so the wall temperature field is exposed as three scalars
  (`T_wall_hot_face`, `T_wall_mean`, `T_wall_cold_face` - actually here:
  max/mean/min over the axial grid, see the property definitions below) rather
  than the full per-node array. Use `get_state()`/`set_state()` on the
  underlying stepper directly (not through FMI) if the full field is needed.
- Construction-time parameters (geometry, chemistry mode, momentum model, the
  hot-gas manifold's O/F and chamber pressure) are NOT exposed as FMI
  variables here - they are fixed when the FMU is built, mirroring the
  `ShellTubeTransientStepper` constraint documented in README.md. To change
  them, edit `_build_config()` below and rebuild the FMU; this is a deliberate
  simplification, not an oversight - varying them live would require
  rebuilding the (disk-cached) combustion manifold inside `do_step`, which is
  out of scope for a first version.
"""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

import numpy as np
from pythonfmu import Boolean, Fmi2Causality, Fmi2Slave, Fmi2Variability, Real

# pythonfmu loads this file as a standalone top-level module (both when
# `build_fmu.py` introspects it and inside the packaged FMU at simulation
# time), so neither a relative import (`from .shelltube_stepper import ...`)
# nor the dotted `hps_combustor.simulink_coupling.shelltube_stepper` path
# (which requires `hps_combustor` itself to already expose this folder as a
# subpackage - not true for the vendored standalone deployment, see
# shelltube_stepper.py's own docstring) can be relied on here. Import the
# sibling module directly instead, by adding its directory to `sys.path`.
#
# Where that sibling actually lives differs by context: when this folder is
# used directly (dev repo, or a standalone copy on disk), shelltube_stepper.py
# sits right next to this file. Inside a *built* FMU, `build_fmu.py` passes
# this folder to pythonfmu as a `project_files` entry, which nests it one
# level deeper as `resources/simulink_coupling/` alongside the top-level
# `resources/fmu_wrapper.py` - so check both locations.
_here = _Path(__file__).resolve().parent
for _candidate in (_here, _here / "simulink_coupling"):
    if (_candidate / "shelltube_stepper.py").is_file() and str(_candidate) not in _sys.path:
        _sys.path.insert(0, str(_candidate))
from shelltube_stepper import BoundaryInputs, ShellTubeTransientStepper  # noqa: E402


def _build_config():
    """Return the config dataclasses used to construct the stepper.

    Edit this function to point at your own combustor/coolant/geometry
    configuration before building the FMU - see SIMULINK_PLUGIN_GUIDE.md.
    """

    from hps_combustor.input_data import (
        CorrelationCoefficients,
        coolantProp,
        hotgasProp,
        numericalProp,
        shellTubeProp,
        system_requirements,
        transientProp,
    )

    cp = coolantProp()
    hp = hotgasProp()
    stp = shellTubeProp()
    nup = numericalProp()
    sr = system_requirements()
    tp = transientProp()
    tp.chemistry_transient = "finite_rate"
    tp.coolant_momentum_model = "low_mach"
    return cp, hp, stp, nup, sr, tp, CorrelationCoefficients()


# N_axial=20 here (not the project's usual 80-node production grid) is a
# deliberately gentle demo choice: at full nominal duty, a finer grid puts
# proportionally more heat into a physically smaller first cell per unit
# time, which can outrun a coarse fixed communication step from a cold
# start (same grid/timestep stiffness interaction already documented in
# docs/context/TRANSIENT_STATUS.md for the legacy schedule-driven solver).
# Raise N_axial only alongside a correspondingly small communication step
# and/or a slower hot-gas ramp-in from Simulink - see SIMULINK_PLUGIN_GUIDE.md.


class ShellTubeTransientFmu(Fmi2Slave):
    """Fmi2Slave exposing `ShellTubeTransientStepper` as an FMU."""

    author = "hps_combustor"
    description = (
        "Shell-and-tube transient HX co-simulation slave "
        "(1Dmodel/simulink_coupling)."
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # --- FMI inputs (set by the master algorithm before each do_step) --
        self.mdot_coolant = 0.075
        self.p_coolant_in = 80.0e5
        self.p_coolant_out = 78.0e5
        self.T_coolant_in = 303.15
        self.mdot_hot_total = 0.1
        self.ignited = True
        self.T_lox_in = 300.0

        # --- FMI outputs (written by do_step) ------------------------------
        self.T_coolant_outlet = 303.15
        self.p_coolant_outlet = 78.0e5
        self.dp_coolant_hydraulic_Pa = 0.0
        self.p_coolant_outlet_predicted = 78.0e5
        self.T_gas_outlet = 300.0
        self.duty_W = 0.0
        self.T_wall_max = 293.15
        self.T_wall_mean = 293.15
        self.T_wall_min = 293.15
        self.mdot_coolant_achieved = 0.0
        self.energy_residual_J = 0.0
        self.mass_residual_kg = 0.0

        self.register_variable(Real("mdot_coolant", causality=Fmi2Causality.input))
        self.register_variable(Real("p_coolant_in", causality=Fmi2Causality.input))
        self.register_variable(Real("p_coolant_out", causality=Fmi2Causality.input))
        self.register_variable(Real("T_coolant_in", causality=Fmi2Causality.input))
        self.register_variable(Real("mdot_hot_total", causality=Fmi2Causality.input))
        self.register_variable(Boolean("ignited", causality=Fmi2Causality.input))
        self.register_variable(Real("T_lox_in", causality=Fmi2Causality.input))

        self.register_variable(Real("T_coolant_outlet", causality=Fmi2Causality.output))
        self.register_variable(Real("p_coolant_outlet", causality=Fmi2Causality.output))
        self.register_variable(Real("dp_coolant_hydraulic_Pa", causality=Fmi2Causality.output))
        self.register_variable(Real("p_coolant_outlet_predicted", causality=Fmi2Causality.output))
        self.register_variable(Real("T_gas_outlet", causality=Fmi2Causality.output))
        self.register_variable(Real("duty_W", causality=Fmi2Causality.output))
        self.register_variable(Real("T_wall_max", causality=Fmi2Causality.output))
        self.register_variable(Real("T_wall_mean", causality=Fmi2Causality.output))
        self.register_variable(Real("T_wall_min", causality=Fmi2Causality.output))
        self.register_variable(
            Real("mdot_coolant_achieved", causality=Fmi2Causality.output)
        )
        self.register_variable(
            Real("energy_residual_J", causality=Fmi2Causality.output, variability=Fmi2Variability.discrete)
        )
        self.register_variable(
            Real("mass_residual_kg", causality=Fmi2Causality.output, variability=Fmi2Variability.discrete)
        )

        self._stepper: ShellTubeTransientStepper | None = None

    def setup_experiment(self, start_time: float, stop_time=None, tolerance=None) -> None:
        cp, hp, stp, nup, sr, tp, corr = _build_config()
        self._stepper = ShellTubeTransientStepper(
            cp, hp, stp, nup, sr, tp,
            corrCoeffs=corr, N_axial=20, flow_config="co",
            p_coolant_out_initial=self.p_coolant_out,
        )

    def do_step(self, current_time: float, step_size: float) -> bool:
        if self._stepper is None:
            # Some master algorithms skip setup_experiment; build lazily.
            cp, hp, stp, nup, sr, tp, corr = _build_config()
            self._stepper = ShellTubeTransientStepper(
                cp, hp, stp, nup, sr, tp,
                corrCoeffs=corr, N_axial=20, flow_config="co",
                p_coolant_out_initial=self.p_coolant_out,
            )

        boundary = BoundaryInputs(
            mdot_coolant=float(self.mdot_coolant),
            p_coolant_in=float(self.p_coolant_in),
            p_coolant_out=float(self.p_coolant_out),
            T_coolant_in=float(self.T_coolant_in),
            mdot_hot_total=float(self.mdot_hot_total),
            ignited=bool(self.ignited),
            T_lox_in=float(self.T_lox_in),
        )
        out = self._stepper.step(float(step_size), boundary)

        self.T_coolant_outlet = out.T_coolant_outlet
        self.p_coolant_outlet = out.p_coolant_outlet
        self.dp_coolant_hydraulic_Pa = out.dp_coolant_hydraulic_Pa
        # A genuinely predicted outlet pressure (inlet minus the model's own
        # computed friction loss) - unlike p_coolant_outlet, this is safe to
        # feed back into p_coolant_out for a real closed loop, since it is not
        # forced to equal p_coolant_out by construction. See README.md.
        self.p_coolant_outlet_predicted = float(self.p_coolant_in) - out.dp_coolant_hydraulic_Pa
        self.T_gas_outlet = out.T_gas_outlet
        self.duty_W = out.duty_W
        self.T_wall_max = float(np.max(out.T_wall))
        self.T_wall_mean = float(np.mean(out.T_wall))
        self.T_wall_min = float(np.min(out.T_wall))
        self.mdot_coolant_achieved = float(np.mean(np.abs(out.face_mdot)))
        self.energy_residual_J = out.energy_residual_J
        self.mass_residual_kg = out.mass_residual_kg
        return True
