"""Fluid-agnostic quasi-1D FV core (new, staged migration).

See docs/solver_design/FV_CORE_REWORK_PLAN.md. This package is built up one
stage at a time alongside the existing maintained solvers; nothing here is
wired into the legacy solvers until each stage's acceptance gate passes.
"""
