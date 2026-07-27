"""Deprecated location. Use ``hps_combustor.physics.gas_flow.governing_equations``
instead (compressible ideal-gas quasi-1D equations — NOT related to the liquid
``p,h`` equations in ``physics/liquid_flow/governing_equations.py``). Kept as a
re-export shim during the Phase 1 restructure in
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md.
"""

import warnings

from hps_combustor.physics.gas_flow.governing_equations import (
    dp__dx_IdealGas,
    dp__dx_IdealGas_logical,
    drho__dx_IdealGas,
    drho__dx_IdealGas_logical,
    dT__dx_IdealGas,
    dT_g__dx,
    dU__dx_IdealGas,
)

warnings.warn(
    "hps_combustor.physics.governing_equations has moved to "
    "hps_combustor.physics.gas_flow.governing_equations; update imports.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "dp__dx_IdealGas",
    "dp__dx_IdealGas_logical",
    "drho__dx_IdealGas",
    "drho__dx_IdealGas_logical",
    "dT__dx_IdealGas",
    "dT_g__dx",
    "dU__dx_IdealGas",
]
