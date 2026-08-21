"""Regenerate `_vendor/hps_combustor/` - the private, self-contained copy of
exactly the `hps_combustor` files `shelltube_stepper.py`/`fmu_wrapper.py`
actually import, plus the combustion mechanism files `choose_fuel()` can
reference. This is what lets someone copy *only* the `simulink_coupling`
folder - nothing else from this repo - and have it still work.

The file list below was derived empirically: constructed a
`ShellTubeTransientStepper`, ran a couple of `step()` calls, and inspected
`sys.modules` for everything under `hps_combustor.*` that actually got
loaded. Re-run this discovery (see the comment at the bottom of this file)
and regenerate whenever `shelltube_stepper.py`'s imports change, rather than
hand-editing the list - a stale vendor copy that silently diverges from what
the stepper actually needs is worse than no vendor copy at all.

Usage (from anywhere):

    python 1Dmodel/simulink_coupling/vendor_dependencies.py

Regenerates `_vendor/hps_combustor/` from scratch (removes any prior copy
first) using whatever `1Dmodel/` currently contains.
"""

from __future__ import annotations

import shutil
from pathlib import Path

# Relative to `1Dmodel/`. Every hps_combustor.* module actually imported by
# shelltube_stepper.py + fmu_wrapper.py, as traced via sys.modules - see the
# module docstring.
_PYTHON_FILES = [
    "__init__.py",
    "input_data.py",
    "main_solve_shellntube.py",
    "main_solve_shellntube_transient.py",
    "mechanical/__init__.py",
    "mechanical/geometry/__init__.py",
    "mechanical/geometry/helix_geometry.py",
    "mechanical/geometry/shelltube_geometry.py",
    "mechanical/loads.py",
    "mechanical/material_specs/__init__.py",
    "mechanical/material_specs/material_temperature_strength.py",
    "physics/__init__.py",
    "physics/bell_delaware.py",
    "physics/combustion_chemistry/__init__.py",
    "physics/combustion_chemistry/combustion_gas.py",
    "physics/combustion_chemistry/fpv_manifold.py",
    "physics/combustion_chemistry/gas_manifold.py",
    "physics/friction_correlations.py",
    "physics/heat_conduction.py",
    "physics/heat_transfer_correlations.py",
    "physics/liquid_flow/__init__.py",
    "physics/liquid_flow/chf.py",
    "physics/liquid_flow/coolprop_state_cache.py",
    "physics/liquid_flow/correlations.py",
    "physics/liquid_flow/dispatch.py",
    "physics/liquid_flow/governing_equations.py",
    "physics/liquid_flow/hx_adapters.py",
    "physics/liquid_flow/regime.py",
    "physics/liquid_flow/registry.py",
    "physics/liquid_flow/supercritical.py",
    "physics/radiation_model/__init__.py",
    "physics/radiation_model/radiation_equations.py",
    "transient_core/__init__.py",
    "transient_core/adapters_helical.py",
    "transient_core/adapters_shelltube.py",
    "transient_core/compressible_coolant.py",
    "transient_core/coolant_fv.py",
    "transient_core/diagnostics.py",
    "transient_core/grid.py",
    "transient_core/integrator.py",
    "transient_core/progress.py",
    "transient_core/schedules.py",
    "transient_core/state.py",
    "transient_core/wall_compressible_coolant.py",
    "transient_core/wall_coolant.py",
]

# Cantera mechanism files `choose_fuel()` (combustion_gas.py) can select
# among for the fuels it supports. Only ~1.3 MB total - cheap enough to
# include all of them so editing `_build_config()`'s `hotgasProp.fuel` still
# works standalone, rather than vendoring only the current default's file.
_MECHANISM_FILES = [
    "physics/combustion_chemistry/A2highT.yaml",
    "physics/combustion_chemistry/H2-O2_Burke2012.yaml",
    "physics/combustion_chemistry/llnl_gasoline_323.yaml",
    "physics/combustion_chemistry/RenKokjohn_surrogate.yaml",
]


def regenerate_vendor() -> Path:
    model_dir = Path(__file__).resolve().parents[1]  # .../1Dmodel
    vendor_root = Path(__file__).resolve().parent / "_vendor" / "hps_combustor"

    if vendor_root.exists():
        shutil.rmtree(vendor_root)

    for rel in _PYTHON_FILES + _MECHANISM_FILES:
        src = model_dir / rel
        dst = vendor_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

    return vendor_root


if __name__ == "__main__":
    path = regenerate_vendor()
    print(f"Regenerated vendor copy at: {path}")
    print(
        "To re-derive the file list itself (only needed if "
        "shelltube_stepper.py's imports change), run:\n"
        "  python -c \"\n"
        "import sys\n"
        "from hps_combustor.input_data import *\n"
        "from hps_combustor.simulink_coupling import ShellTubeTransientStepper, BoundaryInputs\n"
        "cp,hp,stp,npr,sr,tp = coolantProp(),hotgasProp(),shellTubeProp(),numericalProp(),system_requirements(),transientProp()\n"
        "tp.coolant_momentum_model='low_mach'\n"
        "s = ShellTubeTransientStepper(cp,hp,stp,npr,sr,tp,corrCoeffs=CorrelationCoefficients(),N_axial=6,flow_config='co')\n"
        "s.step(0.01, BoundaryInputs(cp.mass_flow_c,cp.p_in,cp.p_in-2e5,cp.T_in,0.02,True))\n"
        "print('\\\\n'.join(sorted(m for m in sys.modules if m.startswith('hps_combustor.'))))\n"
        "\""
    )
