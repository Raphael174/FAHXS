# Methodology

## Two-regime decomposition

A diesel/O2 (or similar two-feed-stream) heat exchanger has two physically
distinct zones, and this kit deliberately uses a different reduced model for
each:

1. **Mixing/flame zone** (`flamelet.py`): fuel and oxidizer arrive as two
   separate streams and burn in a thin non-premixed flame. A conserved
   scalar -- the mixture fraction Z -- exists here (it is materially
   conserved: it is 0 in the oxidizer stream, 1 in the fuel stream, and
   linear in element mass fractions everywhere in between). Composition and
   temperature are then well-approximated as functions of Z alone (plus one
   scalar, chi_st, controlling how fast Z mixes) -- this is the flamelet
   assumption.
2. **Cooling channel** (`cooling_pfr.py`): downstream of the flame the two
   streams no longer exist as separate entities -- there is one already-
   mixed hot-gas stream. Z is no longer a useful coordinate (it's uniform,
   or nearly so, across the channel); what matters now is how that single
   stream's composition and temperature evolve along the channel axis as it
   gives up heat to the walls. This is a plug-flow-reactor (PFR) problem,
   not a flamelet problem.

The coupling is one-directional and simple: regime 1's burnt/product state
(evaluated at or near Z_st, or wherever the channel's bulk composition sits)
becomes regime 2's inlet condition.

## Regime 1: the Peters unsteady flamelet

Governing equations in mixture-fraction space Z in [0,1], unity Lewis number,
low Mach:

    rho dY_k/dt = rho (chi/2) d^2Y_k/dZ^2  +  wdot_k
    rho cp dT/dt = rho cp (chi/2) d^2T/dZ^2
                   + rho (chi/2) [dcp/dZ + sum_k cp_k dY_k/dZ] dT/dZ
                   - sum_k h_k wdot_k  +  dp/dt

- Z=0 boundary: pure oxidizer stream (T_ox, Y_ox), Dirichlet.
- Z=1 boundary: pure fuel stream (T_fuel, Y_fuel), Dirichlet.
- `chi` = chi(Z): the scalar dissipation rate, physically the local rate at
  which Z-gradients are destroyed by molecular/turbulent mixing (units 1/s).
  It is the ONE quantity that couples this 0-D-in-physical-space equation set
  to the actual flow's mixing intensity -- everything else is intrinsic
  chemistry/thermodynamics. `chi(Z)` is parameterized by its value at the
  stoichiometric surface, `chi_st`, via the Peters counterflow form:

      chi(Z) = chi_st * exp(2*erfcinv(2*Z_st)^2 - 2*erfcinv(2*Z)^2)

  This profile's global maximum sits at Z=0.5 (the counterflow's geometric
  mixing-layer center; erfcinv(1)=0), independent of where Z_st is -- NOT at
  Z_st. chi_st is an anchor value: chi(Z_st) == chi_st by construction, not
  the profile's peak. (See `flamelet.py::_chi_profile` and the
  corresponding test in `tests/test_flamelet_kit.py`.)

  **What chi_st physically means for a heat-exchanger user**: it is the
  inverse of a local mixing time scale at the flame surface. In the source
  rocket application it came from a wall-shear-derived turbulent-strain
  closure (not portable). For a heat exchanger, estimate it instead from the
  local flow's strain field or residence time -- see `ADAPTATION_GUIDE.md`
  for concrete routes.

- Z_st (stoichiometric mixture fraction) is computed once from the two feed
  stream compositions via the Bilger (1990) coupling function using
  elemental mass fractions of C, H, O:

      beta(Y) = sum_k (Y_k / W_k) * (2*n_C(k) + 0.5*n_H(k) - n_O(k))
      Z_st = -beta(Y_ox) / (beta(Y_fuel) - beta(Y_ox))

- The Z grid is a two-sided tanh (Roberts) stretching, clustered at Z_st
  from both sides, so the (typically thin) reaction zone is well resolved
  without wasting nodes on the smooth boundary regions.

### Numerics: Strang splitting

Each `step(dt)` does:

    CN diffusion half-step (dt/2)  ->  Cantera chemistry (dt)  ->  CN diffusion half-step (dt/2)

- **Diffusion half-step**: Crank-Nicolson on the non-uniform Z grid, solved
  by the Thomas algorithm (tridiagonal, interior nodes only; boundaries are
  Dirichlet and held fixed). Species are diffused independently (with an
  optional per-species `diff_mask` for non-unity-Lewis species); temperature
  additionally carries the cp-gradient ("Pitsch-Peters") convection term,
  which is treated FULLY IMPLICIT (frozen coefficient, folded into the same
  tridiagonal solve) -- this is a validated stability fix from the source
  solver: an explicit treatment of this term overshoots on the Z_st-clustered
  grid and drives a numerical limit cycle.
- **Chemistry full-step**: each interior Z-node is an independent 0-D
  constant-pressure Cantera reactor (`IdealGasConstPressureReactor` +
  `ReactorNet`, CVODE-integrated) advanced by the full `dt`. Reactors are
  persistent across steps (created once by `init_mixing`) so state carries
  over exactly; boundary nodes are never advanced (they are Dirichlet).

## Regime 2: the cooling plug-flow reactor

A single already-mixed stream, constant p, constant mdot, marched along
length x with an imposed wall heat-loss term:

    mdot * dh/dx = -q_wall(x),   q_wall(x) = h_conv * P_wetted * (T(x) - T_wall(x))

Per axial segment dx, an analogous split is used:

1. Chemistry sub-step: advance the (adiabatic) constant-pressure reactor
   over dt = dx / u(x), u(x) = mdot / (rho(x) * A_cross).
2. Wall heat-loss sub-step: remove Q_seg = q_wall(x)*dx [W] from the gas by
   setting its specific enthalpy directly (`gas.HPY = h - Q_seg/mdot, p, Y`)
   -- exact by construction, so the energy balance closes to reactor
   tolerance: `mdot*(h_in - h_out) == Q_wall_total`.

No mixture fraction, no second stream, no diffusion operator -- this is
intentionally the simpler PFR idiom, not a flamelet.

## Steady-cache rationale

Both regimes' expensive operation is "run the finite-rate solve" (a flamelet
Strang step across all Z-nodes, or a PFR march along the channel). If the
governing inputs (pressure, chi_st or mdot, boundary/inlet temperature) have
not moved outside tolerance since the last solve, the previously computed
solution is still valid and re-solving is wasted work -- skip it.

This is exactly the reuse gate from the source rocket RIF manager
(`_cache_gate_check`), generalized and with the ignition-specific "burning"
precondition removed (see `ADAPTATION_GUIDE.md`). The regime contrast that
motivates keeping this in a STEADY kit: in the rocket's transient
combustion, chi_st jitter (median ~6% relative, mean ~8% advance-to-advance)
tripped the ~5% tolerance band on most advances, collapsing the cache hit
rate to ~0.3%. A steady or steadily-ramping heat-exchanger flow has none of
that acoustic/ignition jitter -- p, chi_st/mdot, and the boundary states move
slowly and (if at all) monotonically between calls, so the identical
tolerance bands that starved the rocket cache are expected to give a HIGH
hit rate here. This is the favorable regime for the technique.
