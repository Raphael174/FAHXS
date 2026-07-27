# Steady Liquid/Boiling Flow: Required Changes (From Scratch)

Starting point: a quasi-1D steady solver whose governing equations, state
variables, and correlations are all derived for/valid on an ideal gas.
Target: any HX design in the same architecture running a genuine liquid
coolant that boils. This is the build list, independent of any prior
implementation attempt.

## 1. State variables and EOS

Replace `(T, p)` with `(p, h)` as the coolant march state. `T` and `p` are not
independent inside the two-phase dome (`T = T_sat(p)`), so a `(T, p)` state
cannot represent a two-phase point at all, and an energy update expressed as
`cp·ΔT` has no latent-heat sink — the coolant would heat straight through the
boiling point with no phase change. Every property lookup (density, T,
quality, transport properties) must come from a real-fluid `(p, h)` flash
(e.g. CoolProp), not the ideal-gas relations currently in place. This is the
precondition everything else depends on — do it first, and verify it against
the existing gas baseline before touching anything else (non-regression is
the cheapest test available and isolates refactor bugs from new-physics
bugs).

## 2. Momentum equation

`-dp/dz = friction + acceleration`. An ideal-gas friction-only pressure drop
is insufficient once boiling starts: the accelerational term,
`G²·d(1/ρ)/dz`, arises because density collapses by orders of magnitude
through the dome, and can dominate the friction term near boiling onset. Both
terms must be carried; dropping acceleration is a real, not negligible,
simplification. Two-phase friction should use a correlation meant for it
(Müller-Steinhagen-Heck or Friedel), not the single-phase Darcy/Colebrook
form outside the dome.

## 3. Heat-transfer/friction correlation chain

A single correlation cannot span the whole state space. Required chain,
each valid only in its regime:

1. Subcooled liquid — Gnielinski/Dittus-Boelter-type single-phase.
2. Onset of nucleate boiling — a real ONB criterion (Bergles-Rohsenow or
   similar), not a bare quality threshold.
3. Saturated flow boiling — Chen or Gungor-Winterton.
4. **CHF check** — a lookup (Groeneveld 2006 LUT is the standard workhorse)
   that must actively override the boiling HTC once exceeded, distinguishing
   DNB (low quality) from dryout (high quality) — a diagnostic-only CHF
   number that never feeds back into the HTC silently under-predicts wall
   temperature past the real limit.
5. Post-CHF/post-dryout — a degraded HTC closure (a conservative
   single-phase-vapor-at-full-flux approximation is a defensible minimum; a
   real dispersed-flow/mist correlation is better).
6. Superheated vapor — single-phase, real-fluid properties.
The transition between adjacent regimes (especially 1→3) must be smoothed
over a small window, not hard-switched — a hard switch produces a spurious
step-change in HTC (and therefore wall temperature) at a single node that is
a numerical artifact, not physics.

## 4. Sound speed / compressibility

`c = √(γRT)` is an ideal-gas relation and does not apply here. Use the real
EOS sound speed outside the dome and Wood's equation (void-fraction weighted
mixture sound speed) inside it. This is not cosmetic: two-phase mixture sound
speed can collapse to a few hundred m/s or less, so choking is a real risk at
velocities that would be entirely unremarkable for a single-phase gas or
liquid — any Mach-number-based diagnostic needs this to mean anything for the
coolant side.

## 5. Numerical robustness

- **EOS validity bounds**: property back-ends have hard temperature/pressure
  limits (e.g. ~3000 K ceilings are typical). Very low coolant flow relative
  to duty can drive superheated coolant past this ceiling. The march must
  detect and fail gracefully (hold last valid state, flag the run), not
  crash — this is a real, reachable operating point, not a corner case to
  ignore.
- **Boundary conditions**: a single-phase `(T, p)` guess cannot seed a
  two-phase starting state for any configuration with a prescribed *outlet*
  (e.g. counter-flow). Use a shooting method on the physically-known inlet
  state (enthalpy, not temperature) instead of guessing the far boundary.

## 6. Validation

- Non-regression against the existing gas-coolant baseline (do this before
  and after every phase above).
- Energy balance closure (`Σ wall duty` vs. `ṁ·Δh`) across the full march —
  should agree to numerical precision; this is the strongest single defect
  detector for a state-variable refactor.
- Validation against a literature benchmark with independent wall-temperature
  and CHF-location data for the target fluid/geometry class (e.g. Bennett et
  al. 1967 for water in heated tubes) — a small-scale correlation-fit
  comparison alone is not sufficient evidence for design use.

## Sequencing note

Items 1-2 are a from-scratch refactor and carry the most schedule risk
(saturation-boundary numerics, in particular Jacobian/derivative
discontinuities if any implicit scheme is used). Item 3 depends on 1. Items 4
and 5 can be added incrementally once 1-3 are stable. Item 6 should run
continuously, not as a final gate — regime-switching logic makes it easy to
fix one regime and silently break another.
