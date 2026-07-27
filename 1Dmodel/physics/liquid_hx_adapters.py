"""Deprecated location. Use ``hps_combustor.physics.liquid_flow.hx_adapters``
instead. Kept as a re-export shim during the Phase 1 restructure in
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md.
"""

import warnings

from hps_combustor.physics.liquid_flow.hx_adapters import (
    solve_helical_coil_liquid_from_data_master,
    solve_helical_coil_liquid_from_duty,
    solve_shelltube_shellside_liquid_from_duty,
    solve_shelltube_shellside_liquid_from_tube_result,
)

warnings.warn(
    "hps_combustor.physics.liquid_hx_adapters has moved to "
    "hps_combustor.physics.liquid_flow.hx_adapters; update imports.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "solve_helical_coil_liquid_from_data_master",
    "solve_helical_coil_liquid_from_duty",
    "solve_shelltube_shellside_liquid_from_duty",
    "solve_shelltube_shellside_liquid_from_tube_result",
]
