"""Build a transient HTML dashboard from saved run data.

Usage:
  hps-dashboard path/to/transient_timeseries.npz
  hps-dashboard path/to/run_folder
  hps-dashboard path/to/run_archive.zip
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .model_data_process.data_plotting_transient import TransientDashboard
from .result_package import load_transient_time_series


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(
        description="Generate a self-contained transient dashboard HTML file from saved data."
    )
    parser.add_argument(
        "input",
        help="Path to transient_timeseries.npz, a packaged run folder, or a run .zip archive.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output HTML path. Defaults beside the input data or archive.",
    )
    args = parser.parse_args(argv)

    input_path = Path(args.input)
    time_series = load_transient_time_series(input_path)
    output_path = Path(args.output) if args.output else _default_output_path(input_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    dashboard = TransientDashboard(
        time_series,
        meta={"source": str(input_path), "output": str(output_path)},
    )
    dashboard.to_html(output_path)
    print(f"Transient dashboard written: {output_path}")
    return 0


def _default_output_path(input_path: Path):
    if input_path.suffix.lower() == ".npz":
        return input_path.with_name("transient_dashboard.html")
    if input_path.suffix.lower() == ".zip":
        return input_path.with_suffix(".html")
    return input_path / "transient_dashboard.html"


if __name__ == "__main__":
    main()
