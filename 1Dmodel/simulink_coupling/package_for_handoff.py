"""Zip a self-contained handoff package for deploying this solver (including
`simulink_coupling/`) to another machine, e.g. one that runs Simulink.

This does NOT make the code dependency-free - see README.md's "What
'independent module' actually means here". It packages the real dependency
closure (the whole `1Dmodel/` package plus `pyproject.toml`) so the result is
installable with `pip install -e .` on the target machine, without needing
this repo's git history, docs, tests, or accumulated run outputs.

Usage (from the repository root):

    python 1Dmodel/simulink_coupling/package_for_handoff.py path/to/output.zip
"""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path

_EXCLUDE_DIR_NAMES = {"__pycache__", ".pytest_cache"}
_EXCLUDE_SUFFIXES = {".pyc"}


def _should_include(path: Path) -> bool:
    if any(part in _EXCLUDE_DIR_NAMES for part in path.parts):
        return False
    if path.suffix in _EXCLUDE_SUFFIXES:
        return False
    return True


def build_handoff_zip(output_path: Path) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    model_dir = repo_root / "1Dmodel"
    pyproject = repo_root / "pyproject.toml"

    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in model_dir.rglob("*"):
            if path.is_file() and _should_include(path):
                zf.write(path, path.relative_to(repo_root))
        zf.write(pyproject, pyproject.relative_to(repo_root))

    return output_path


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 1:
        print("usage: package_for_handoff.py <output.zip>", file=sys.stderr)
        return 2
    result = build_handoff_zip(Path(argv[0]))
    print(f"Wrote handoff package: {result}")
    print("On the target machine: unzip, then `pip install -e .` "
          "(plus Cantera/CoolProp) before building/running the FMU.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
