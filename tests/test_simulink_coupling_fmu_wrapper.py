"""Tests for the FMI slave wrapper.

`test_fmu_slave_steps_without_error` only checks that `ShellTubeTransientFmu`
can be instantiated and stepped through the plain `pythonfmu.Fmi2Slave`
Python API in this dev environment - it catches wiring mistakes between the
FMI variable definitions and `ShellTubeTransientStepper`, but it does NOT
prove the *packaged* `.fmu` file is self-contained (pythonfmu bundles only
what `project_files` explicitly lists; running the unpacked script in this
dev repo can silently succeed via the real `hps_combustor` install even if
the actual .fmu is missing a file).

`test_packaged_fmu_runs_standalone` is the real proof for that: builds the
FMU with `build_fmu.py`, unzips it exactly as an FMI master would, and runs
`fmu_wrapper.py` from that extracted `resources/` folder in a subprocess with
this repo's editable `hps_combustor` install hidden - the same technique as
`test_simulink_coupling_standalone.py`. This is the regression test for the
bug where `build_fmu.py` bundled only `fmu_wrapper.py` itself and the FMU
silently depended on `shelltube_stepper.py`/`_vendor/` happening to sit next
to it on disk in dev, rather than actually being inside the FMU zip.
"""

import subprocess
import sys
import textwrap
import zipfile
from pathlib import Path

import pytest

pytest.importorskip("pythonfmu")

from hps_combustor.simulink_coupling.fmu_wrapper import ShellTubeTransientFmu

SIMULINK_COUPLING_DIR = Path(__file__).resolve().parents[1] / "1Dmodel" / "simulink_coupling"

_ISOLATED_RUN_SCRIPT = textwrap.dedent(
    """
    import sys

    sys.meta_path = [
        f for f in sys.meta_path if "editable" not in getattr(f, "__name__", "").lower()
    ]
    sys.path.insert(0, sys.argv[1])

    import fmu_wrapper

    slave = fmu_wrapper.ShellTubeTransientFmu(instance_name="pkgTest")
    slave.setup_experiment(0.0)
    slave.p_coolant_in = 80.0e5
    slave.p_coolant_out = 78.0e5
    slave.mdot_coolant = 0.075
    slave.T_coolant_in = 303.15
    slave.ignited = True
    t = 0.0
    for _ in range(3):
        slave.mdot_hot_total = min(0.005 + 0.095 * (t / 2.5), 0.10)
        assert slave.do_step(t, 0.01) is True
        t += 0.01

    import hps_combustor as hc
    assert "_vendor" in hc.__file__, f"resolved outside _vendor: {hc.__file__}"
    print("PACKAGED_FMU_STANDALONE_OK")
    """
)


def test_fmu_slave_steps_without_error():
    slave = ShellTubeTransientFmu(instance_name="testSlave")
    slave.setup_experiment(0.0)

    slave.p_coolant_in = 80.0e5
    slave.p_coolant_out = 78.0e5
    slave.mdot_coolant = 0.075
    slave.T_coolant_in = 303.15
    slave.ignited = True

    t = 0.0
    for _ in range(3):
        slave.mdot_hot_total = min(0.005 + 0.095 * (t / 2.5), 0.10)
        assert slave.do_step(t, 0.01) is True
        assert slave.T_coolant_outlet > 0.0
        assert slave.p_coolant_outlet > 0.0
        assert abs(slave.mass_residual_kg) < 1e-6
        t += 0.01


def test_packaged_fmu_runs_standalone(tmp_path):
    from pythonfmu import FmuBuilder

    script = SIMULINK_COUPLING_DIR / "fmu_wrapper.py"
    fmu_path = FmuBuilder.build_FMU(script, dest=tmp_path, project_files=[script.parent])

    extracted = tmp_path / "extracted"
    with zipfile.ZipFile(fmu_path) as zf:
        zf.extractall(extracted)

    script_path = tmp_path / "isolated_run.py"
    script_path.write_text(_ISOLATED_RUN_SCRIPT, encoding="utf-8")

    result = subprocess.run(
        [sys.executable, str(script_path), str(extracted / "resources")],
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert "PACKAGED_FMU_STANDALONE_OK" in result.stdout, (
        f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
    )
    assert result.returncode == 0, f"stderr:\n{result.stderr}"
