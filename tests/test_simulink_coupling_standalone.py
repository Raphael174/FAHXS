"""Proves `1Dmodel/simulink_coupling/` genuinely works standalone - copied by
itself, with no access to the rest of this repo or a separately installed
`hps_combustor` package.

This is the regression test for the actual complaint that motivated
vendoring: a copy of just this folder was handed to someone without the rest
of the repo, and it failed on `import hps_combustor`. If someone edits
`shelltube_stepper.py`/`fmu_wrapper.py` to import something not covered by
`vendor_dependencies.py`'s file list, this test fails with an `ImportError`
inside the subprocess, catching the drift immediately instead of silently
shipping a folder that only works in this dev environment.

Runs the isolated check in a subprocess (not in-process sys.path/meta_path
surgery in the main test process) so it can't leak import state into the
rest of the test suite, and so it reflects what a real separate Python
process on another machine would see.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

SIMULINK_COUPLING_DIR = Path(__file__).resolve().parents[1] / "1Dmodel" / "simulink_coupling"

_ISOLATED_SCRIPT = textwrap.dedent(
    """
    import sys

    # Simulate a machine that never installed this repo's hps_combustor:
    # remove the dev repo's editable-install meta path finder so `import
    # hps_combustor` can only succeed via the vendored copy shipped inside
    # the (copied) simulink_coupling folder itself.
    sys.meta_path = [
        f for f in sys.meta_path if "editable" not in getattr(f, "__name__", "").lower()
    ]
    sys.path.insert(0, sys.argv[1])

    try:
        import hps_combustor  # noqa: F401
        print("FAIL: hps_combustor importable without the vendor fallback")
        raise SystemExit(1)
    except ImportError:
        pass

    from shelltube_stepper import ShellTubeTransientStepper, BoundaryInputs
    import hps_combustor as hc
    from hps_combustor.input_data import (
        CorrelationCoefficients, coolantProp, hotgasProp, numericalProp,
        shellTubeProp, system_requirements, transientProp,
    )

    assert "_vendor" in hc.__file__, f"hps_combustor resolved outside _vendor: {hc.__file__}"

    cp, hp = coolantProp(), hotgasProp()
    stp, npr = shellTubeProp(), numericalProp()
    sr, tp = system_requirements(), transientProp()
    tp.coolant_momentum_model = "low_mach"
    stepper = ShellTubeTransientStepper(
        cp, hp, stp, npr, sr, tp,
        corrCoeffs=CorrelationCoefficients(), N_axial=6, flow_config="co",
    )
    boundary = BoundaryInputs(
        mdot_coolant=cp.mass_flow_c, p_coolant_in=cp.p_in,
        p_coolant_out=cp.p_in - 2e5, T_coolant_in=cp.T_in,
        mdot_hot_total=0.02, ignited=True,
    )
    out = stepper.step(0.01, boundary)
    assert out.T_coolant_outlet > 0.0
    print("STANDALONE_OK")
    """
)


def test_simulink_coupling_folder_is_standalone(tmp_path):
    copied = tmp_path / "simulink_coupling"
    shutil.copytree(
        SIMULINK_COUPLING_DIR,
        copied,
        ignore=shutil.ignore_patterns("__pycache__"),
    )

    script_path = tmp_path / "isolated_check.py"
    script_path.write_text(_ISOLATED_SCRIPT, encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(script_path), str(copied)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert "STANDALONE_OK" in result.stdout, (
        f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
    )
    assert result.returncode == 0, f"stderr:\n{result.stderr}"
