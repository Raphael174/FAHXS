"""Decoupled Simulink/multiphysics co-simulation layer for the shell-and-tube
transient HX solver.

Nothing outside this folder is modified to provide this capability — see
README.md for the I/O contract and known limitations, and
SIMULINK_PLUGIN_GUIDE.md for how to build and import the FMU into Simulink.
"""

from .shelltube_stepper import BoundaryInputs, ShellTubeTransientStepper, StepOutputs

__all__ = ["BoundaryInputs", "ShellTubeTransientStepper", "StepOutputs"]
