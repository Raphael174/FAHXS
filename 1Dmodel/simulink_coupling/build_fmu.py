"""Build the shell-and-tube transient HX FMU.

Usage (from the repository root, with `pythonfmu` installed in the active
Python environment):

    python -m hps_combustor.simulink_coupling.build_fmu [output_dir]

Produces `ShellTubeTransientFmu.fmu` in `output_dir` (default: current
directory). See SIMULINK_PLUGIN_GUIDE.md for how to import it into Simulink.

The FMU bundles `shelltube_stepper.py` and `_vendor/hps_combustor/` alongside
`fmu_wrapper.py` (via `project_files`) so the built `.fmu` is standalone the
same way this whole folder is: it does NOT require this repo's
`hps_combustor` to be separately installed on whatever machine runs it. It
still requires Python plus `numpy`/`scipy`/`CoolProp`/`cantera` on that
machine - those are not (and cannot practically be) bundled inside an FMU
zip. See README.md's "Standalone deployment" section.
"""

from __future__ import annotations

import sys
from pathlib import Path

from pythonfmu import FmuBuilder


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    output_dir = Path(argv[0]) if argv else Path.cwd()
    output_dir.mkdir(parents=True, exist_ok=True)

    script = Path(__file__).with_name("fmu_wrapper.py")
    fmu_path = FmuBuilder.build_FMU(
        script,
        dest=output_dir,
        # Bundles every sibling file/folder next to fmu_wrapper.py
        # (shelltube_stepper.py, _vendor/hps_combustor/, etc.) into the
        # FMU's resources/ - pythonfmu special-cases "the script's own
        # parent directory" to mean exactly that.
        project_files=[script.parent],
    )
    print(f"Built FMU: {fmu_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
