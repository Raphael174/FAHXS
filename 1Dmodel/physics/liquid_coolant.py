"""Deprecated location. Use ``hps_combustor.physics.liquid_flow.correlations``
(properties/HTC/friction) and ``hps_combustor.physics.liquid_flow.chf`` (CHF
lookup) instead. Kept as a re-export shim during the Phase 1 restructure in
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md.
"""

import warnings

from hps_combustor.physics.liquid_flow.chf import (
    GROENEVELD_2006_MASS_FLUXES,
    GROENEVELD_2006_PRESSURES_MPA,
    GROENEVELD_2006_QUALITIES,
    groeneveld_2006_chf,
    interpolate_chf_table,
    load_groeneveld_2006_lut,
    local_chf_diameter_correction,
)
from hps_combustor.physics.liquid_flow.correlations import (
    EquilibriumState,
    SaturationState,
    chisholm_two_phase_multiplier,
    darcy_friction_smooth_pipe,
    equilibrium_state_ph,
    gungor_winterton_boiling_htc,
    homogeneous_acceleration_pressure_gradient,
    homogeneous_void_fraction,
    liquid_single_phase_nusselt,
    martinelli_parameter_laminar_liquid_turbulent_vapor,
    martinelli_parameter_tt,
    muller_steinhagen_heck_friction_gradient,
    saturation_state,
    thermodynamic_quality,
    yu2002_modified_anl_boiling_htc,
    yu2002_small_channel_pressure_multiplier,
)

warnings.warn(
    "hps_combustor.physics.liquid_coolant has moved to "
    "hps_combustor.physics.liquid_flow.correlations / .chf; update imports.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "GROENEVELD_2006_MASS_FLUXES",
    "GROENEVELD_2006_PRESSURES_MPA",
    "GROENEVELD_2006_QUALITIES",
    "EquilibriumState",
    "SaturationState",
    "chisholm_two_phase_multiplier",
    "darcy_friction_smooth_pipe",
    "equilibrium_state_ph",
    "groeneveld_2006_chf",
    "gungor_winterton_boiling_htc",
    "homogeneous_acceleration_pressure_gradient",
    "homogeneous_void_fraction",
    "interpolate_chf_table",
    "liquid_single_phase_nusselt",
    "load_groeneveld_2006_lut",
    "local_chf_diameter_correction",
    "martinelli_parameter_laminar_liquid_turbulent_vapor",
    "martinelli_parameter_tt",
    "muller_steinhagen_heck_friction_gradient",
    "saturation_state",
    "thermodynamic_quality",
    "yu2002_modified_anl_boiling_htc",
    "yu2002_small_channel_pressure_multiplier",
]
