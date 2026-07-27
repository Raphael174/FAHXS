# Steady Liquid/Boiling Flow: Required Changes (Any Design)

Updates `water_coolant_conversion_plan.md` against what is actually implemented
as of 2026-07-16. Scope: the changes needed to take *any* steady HX design in
this codebase from ideal-gas-derived quasi-1D Navier-Stokes to a valid liquid/
boiling coolant, not just this combustor's specific geometry.

## 1. State variables and EOS — **required, done**

Switch the coolant march from `(T, p)` to `(p, h)`. `T` and `p` are not
independent inside the two-phase dome, so any `(T, p)` energy update (plain
`cp·ΔT`) has no latent-heat sink — this is not a refinement, it is the
precondition for the coolant not heating a boiling fluid straight through
1000+ K. Implemented for both maintained solvers via
`physics/liquid_flow/dispatch.py::evaluate_coolant_closure` /
`correlations.py::equilibrium_state_ph` (CoolProp `(P,H)` flash).

## 2. Momentum equation — **required, done**

`-dp/dz = friction + acceleration`. Friction alone (Müller-Steinhagen-Heck
two-phase / Colebrook-type single-phase) is not sufficient: the accelerational
term, `G²·d(1/ρ)/dz`, can reach **~90% of the total gradient right at boiling
onset** (measured this session) because density collapses through the dome.
**Done** for the helical solver (one-node-lagged, matching the existing
boiling-HTC lag pattern). **Done (2026-07-16)** for the shell-and-tube
solver's shell side too: `_shell_h_at`/`_shell_side_march` now march an actual
(one-sweep-lagged) pressure profile instead of holding `p_in` constant —
Müller-Steinhagen-Heck's friction gradient (already computed for the boiling
HTC call, previously discarded) inside the two-phase dome, and the lumped
whole-bundle Bell-Delaware `dp_shell` (evaluated with local density/viscosity,
apportioned evenly per axial node) outside it. This is still an approximation
— apportioning assumes each axial slice carries an equal share of the total
drop, not a true baffle-by-baffle discretization — but it replaces a flat
assumption with a real, density-aware pressure field that feeds back into
local Tsat(p) and the CHF/HTC closures. Gas mode is unaffected (still holds
`p_in`; out of scope, no phase change to make pressure-driven Tsat shifts
matter there).

## 3. Heat-transfer/friction correlation stack — **required, done**

Subcooled liquid (Gnielinski) → **blended** over a small quality window (not a
hard switch — a hard switch produces a real, spurious ~5x HTC jump in one
node) → saturated boiling (Gungor-Winterton) → **CHF-margin gate** (Groeneveld
2006 LUT) that **must** override the boiling HTC once exceeded, not just
report a number → post-CHF closure (implemented as a conservative
100%-vapor-at-full-flux approximation — real dispersed-flow correlations
would improve accuracy, not correctness) → single-phase vapor. All four
regimes are required; skipping the CHF→HTC feedback silently under-predicts
wall temperature past dryout.

## 4. Sound speed / choking — **required, done**

Ideal-gas `c = √(γRT)` does not apply to a liquid or two-phase mixture. Real
EOS sound speed (CoolProp) outside the dome; **Wood's equation** (void-fraction
weighted, not quality-weighted) inside it — mixture sound speed can collapse
to a few hundred m/s even far from either pure-phase value. Needed to flag
two-phase choking, a real risk at modest velocities once density craters.

## 5. Numerical robustness — **required, done**

- EOS-validity ceiling: very low coolant flow relative to duty can superheat
  past CoolProp's backend temperature limit (~3000 K for water). Must fail
  gracefully (freeze last valid state), not crash.
- Bracket/shooting search for counter-flow (or any prescribed-outlet)
  configurations: a single `(T,p)` guess cannot seed a two-phase state.
  Physically-anchored shooting on hot-end enthalpy is required, not optional,
  for liquid coolant in that configuration.

## 6. Validation — **required, partially done**

- Non-regression vs. the existing gas baseline (cheapest, catches
  implementation bugs before new-physics bugs). **Done**, currently blocked by
  an unrelated pre-existing baseline drift (flagged, not fixed). Re-confirmed
  2026-07-16 after the §2/§7 additions (ONB, CHF regime, shell pressure march):
  helical and shell-and-tube, gas and liquid mode, all ran end-to-end with no
  new crashes and unchanged gas-mode results; the known shell-and-tube water
  edge-of-validity oscillation (50 g/s, ±7 K after 25 sweeps) reproduced
  identically, confirming it is unrelated to this session's changes.
- Energy balance closure across the full march. **Done** (machine-precision
  agreement demonstrated).
- **Not done**: validation against a classic literature wall-temperature/CHF
  benchmark (e.g. Bennett et al. 1967). Current validation (Yu 2002) is a
  narrower small-channel dataset — a materiality gap for absolute-accuracy
  confidence, not a checkbox. This requires digitized external reference data
  this session did not have access to; it remains a genuinely open task, not
  something to approximate from memory.

## 7. Outstanding for a fully general implementation

| Item | Why it matters | Effort | Status |
|---|---|---|---|
| TTSE/bicubic CoolProp tabulation | Raw `PropsSI` calls are ~100x slower than ideal-gas; this is the actual bottleneck behind slow shooting-method convergence | Medium — must be opt-in, validated near the dome where interpolation error is largest | **Done (2026-07-16)**. Two layers: (1) `physics/liquid_flow/coolprop_state_cache.py` reuses a persistent low-level `CP.AbstractState` per (backend, fluid) instead of re-parsing the fluid string on every `PropsSI` call — a pure performance refactor (same HEOS equation of state, confirmed bit-identical `Q_tot` on both solvers after the change), which alone cuts per-call cost roughly 10-15x. (2) `coolantProp.liquid_property_backend` ("HEOS" default / "TTSE" / "BICUBIC") opts into CoolProp's tabulated interpolation backends, measured a further ~10-12x faster per property call (see `validation/liquid_ttse_backend_validation.py` — max relative error ~2e-4 in T/rho and ~1e-9 in quality across a dense grid at the saturation boundaries, 30/80/150 bar). Note: CoolProp's high-level `PropsSI("BICUBIC::Fluid", ...)` string form does **not** work in this environment (raises "cannot be used in the high-level interface") — only the low-level `AbstractState` API supports these backends, which is why the cache module exists rather than a simple string tag. Also found and fixed a real CoolProp quirk: only the *first* `AbstractState` instance of a given tabulated backend+fluid pair can compute `surface_tension()`; every other coexisting instance of that backend raises a spurious "only defined within the two-phase region" even at an exact saturation point — routed through a dedicated HEOS probe instead (unaffected, and surface tension isn't on the hot path anyway). |
| DNB vs. dryout distinction | Currently one undifferentiated CHF check; low-quality and high-quality CHF are different physical mechanisms | Small | **Done (2026-07-16)**: `chf.chf_regime()` classifies by quality at the worst-margin node (fixed x=0.1 cutoff, not a fluid/pressure-dependent boundary — a labeling aid, not new lookup physics), surfaced in both solvers' printed diagnostics and `LiquidMarchSanityReport.chf_regime_at_min_margin` |
| Real ONB criterion (Bergles-Rohsenow) | Current blend is numerical smoothing, not physics | Small | **Done (2026-07-16)**: `correlations.bergles_rohsenow_onb_wall_superheat()` gates an active warning (both solvers) when the estimated wall superheat in the subcooled region exceeds the ONB threshold at the local heat flux/pressure — a diagnostic gate, not a restructuring of the subcooled-region HTC itself (still the blended single-phase closure) |
| Shell-side cross-flow boiling enhancement | Gungor-Winterton is a tube-flow correlation; shell-side cross-flow boiling is unmodeled and likely non-conservative | Medium-large, needs literature closure | Not done — deliberately skipped rather than fabricate an unvalidated correlation extension; still a real gap |
| Shell-side coolant pressure profile | See §2 | Small-medium | **Done (2026-07-16)** — see §2 |

**Bottom line**: items 1-5 are the non-negotiable precondition set — without
them the solver is not just less accurate, it is wrong (unphysical
temperatures, silent CHF blindness, or a crash) once the coolant boils. Item 6
is what lets you trust the numbers. Item 7 is where accuracy improves from
"qualitatively right" to "quantitatively validated for design" — as of
2026-07-16 the one remaining open row is shell-side cross-flow boiling
enhancement (a real physics gap requiring literature grounding this session
did not attempt); TTSE/BICUBIC tabulation is now done.
