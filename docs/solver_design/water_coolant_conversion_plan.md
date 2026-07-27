# Heat Exchanger Model: Helium → Water Coolant Conversion

**Context:** Quasi-1D perfect-gas Navier-Stokes solver for hot-gas/coolant heat exchanger, currently validated for helium. Objective: extend to liquid water coolant, including phase change (subcooled → nucleate boiling → saturated boiling → dryout → superheated steam).

---

## 1. What Stays the Same

The quasi-1D finite-volume conservation laws are already fluid-agnostic:

- Mass: `∂(ρA)/∂t + ∂(ρuA)/∂x = 0`
- Momentum: `∂(ρuA)/∂t + ∂((ρu² + p)A)/∂x = p·dA/dx − τ_w·P_w`
- Energy: `∂(ρEA)/∂t + ∂(ρuHA)/∂x = q″·P_h`

Nothing here assumes perfect gas. The perfect-gas assumption lives entirely in **(1) the equation-of-state closure, (2) the choice of state variables, and (3) the numerics** — that's where the work is.

---

## 2. State Variables: Switch to (p, h)

Temperature is a poor state variable for water: inside the saturation dome, T and p are not independent (`T = T_sat(p)`), so a (p, T) flash is degenerate.

**Use pressure–enthalpy (p, h) instead.** Everything is single-valued in this pair:

- Density: `ρ = ρ(p, h)`
- Temperature: `T = T(p, h)`
- Vapor quality: `x = (h − h_l(p)) / h_fg(p)`

CoolProp supports PH flash natively (`PropsSI('D','P',p,'H',h,'Water')`) across subcooled liquid, two-phase mixture, and superheated steam — and this formulation works unchanged for helium too, making the solver genuinely coolant-agnostic.

---

## 3. Physical Model Hierarchy

1. **Homogeneous Equilibrium Model (HEM)** — both phases share velocity and are at saturation; treated as a pseudo-single fluid `ρ_mix(p, h)`. The three conservation equations survive intact. **Recommended starting point.**
2. **Drift-flux / slip models** — vapor moves faster than liquid (slip ratio S ≠ 1); adds one algebraic closure (e.g., Zuber-Findlay) for better void fraction and pressure drop accuracy.
3. **Two-fluid six-equation models** (RELAP/CATHARE-style) — separate mass/momentum/energy per phase. Likely overkill for a design-level tool.

HEM's main weakness: assumes thermal equilibrium, missing real subcooled boiling (vapor at wall, subcooled core). Acceptable for sizing; wall correlations compensate.

---

## 4. Assumptions That Break

- **Compressibility doesn't disappear, it gets stranger.** Pure liquid sound speed ~1500 m/s → very stiff explicit acoustic CFL (go implicit, low-Mach preconditioning, or quasi-steady marching). In the two-phase region, mixture sound speed collapses (Wood's formula) to tens of m/s — two-phase choking becomes a real risk even at modest velocities.
- **Pressure drop budget shifts.** Density drops ~1000× on evaporation; the accelerational pressure-drop term (negligible for gas) can dominate frictional losses. Must stay in the momentum equation.
- **Wall coupling becomes nonlinear and discontinuous.** Heat transfer coefficient depends nonlinearly on wall superheat; above CHF it drops by an order of magnitude as a temperature jump, not a smooth degradation. If the hot-gas side imposes heat flux, local CHF exceedance causes a wall temperature excursion the solver must handle.

---

## 5. Correlation Stack (replaces single Dittus-Boelter/Gnielinski correlation)

| Regime | Correlation |
|---|---|
| Subcooled liquid | Gnielinski (unchanged) |
| Onset of nucleate boiling (ONB) | Bergles-Rohsenow or Davis-Anderson |
| Saturated flow boiling | Chen (classic), or Gungor-Winterton / Kandlikar / Steiner-Taborek (more accurate) |
| Critical heat flux / dryout | Groeneveld 2006 look-up table (workhorse); Katto-Ohno as alternative. Distinguish DNB (low quality, high flux) from dryout (high quality, annular film depletion) |
| Post-dryout (optional) | Groeneveld film boiling correlations |
| Two-phase friction | Friedel or Müller-Steinhagen & Heck (preferred over Lockhart-Martinelli) |

---

## 6. Numerical Pitfalls

- Property derivatives (`∂ρ/∂p`, `∂ρ/∂h`) are discontinuous at the saturation boundary — breaks Newton iteration/Jacobian convergence.
- CoolProp PH flash calls are ~100× slower than ideal-gas evaluation — enable TTSE or bicubic tabulation (`HEOS&TTSE::Water`).
- Stiffness from sound-speed contrast strongly favors implicit or semi-implicit schemes (why RELAP/CATHARE are semi-implicit) if going transient.

---

## 7. Implementation Roadmap & Effort Estimate

*Assumes side-of-desk work (~half-time), building on the existing validated helium code.*

### Phase 1 — (p, h) refactor + CoolProp closure
**Effort: 1–2 weeks**

Scope: swap state variables, replace ideal-gas closure, add TTSE tabulation.

**Validation:**
- Non-regression vs. existing helium results (exit temperature, pressure drop, wall temperature profile) — should match to <0.1%. Cheapest, most important test; isolates implementation bugs from new physics.
- Single-phase liquid water case in a straight pipe vs. Gnielinski + Colebrook hand calc — pressure drop verifiable to a few percent.
- Exit temperature vs. simple energy balance (`ṁ·cp·ΔT = Q`) — should close to <1%.

### Phase 2 — HEM two-phase properties + flow
**Effort: 2–3 weeks**

Scope: mixture density, two-phase friction multiplier, accelerational pressure drop, mixture sound speed for CFL/choking checks. Most calendar risk lives here (saturation-boundary numerics).

**Validation:**
- Adiabatic two-phase pressure drop vs. published steam-water benchmark data. Note: Friedel itself carries ~30% scatter vs. data — acceptance criterion is "within correlation's own uncertainty," not a tight percentage.
- Energy balance closure across the dome: heat subcooled → superheated, verify exit enthalpy matches integrated wall heat input exactly, and exit steam temperature is consistent with `h_exit` at `p_exit`.
- Quality profile: analytical check — for uniform heat flux, `x(z)` should be linear in the saturated region.

### Phase 3 — Boiling heat transfer chain
**Effort: 2–3 weeks**

Scope: ONB criterion, Chen or Gungor-Winterton, regime-switching logic (blend correlations over a small quality window to avoid solver oscillation at transitions).

**Validation:**
- Wall temperature profiles vs. classic uniformly-heated tube data — **Bennett et al. (1967, AERE Harwell)** steam-water dataset is the canonical reference (matches this problem shape closely). Expect agreement within ±20–30% on h → a few K to tens of K on wall superheat.
- Exit steam temperature for a superheated-exit case: dominated by energy balance + single-phase vapor correlation, should be tight (<2–3%) — good integral check.

### Phase 4 — CHF / dryout
**Effort: 1–2 weeks**

Scope: Groeneveld 2006 look-up table (interpolation + diameter/geometry correction factors), DNB-vs-dryout flagging.

**Validation:**
- CHF location vs. Bennett dataset (reports measured dryout position) — check both critical quality and axial dryout location. Groeneveld table RMS error ~7–8% on CHF against its database; target ~10% of channel length on location.
- Parametric sanity sweeps: CHF should decrease with increasing exit quality, increase with mass flux (at low quality), decrease with pressure above ~7 MPa. Cheap check that catches sign/indexing errors.

### Phase 5 (optional) — Post-dryout film boiling
**Effort: 1–2 weeks**

Only needed if the system must survive past CHF rather than flag-and-stop. For a design code, treating CHF as a hard constraint is defensible and this phase can be skipped initially.

---

## 8. Summary Timeline

| Phase | Scope | Effort |
|---|---|---|
| 1 | (p,h) refactor + CoolProp closure | 1–2 weeks |
| 2 | HEM two-phase flow | 2–3 weeks |
| 3 | Boiling heat transfer chain | 2–3 weeks |
| 4 | CHF/dryout | 1–2 weeks |
| 5 (optional) | Post-dryout film boiling | 1–2 weeks |

**Total: ~6–10 weeks part-time (~3–5 weeks full-time equivalent)** for a validated HEM code with CHF checking.

**Key risk driver:** not coding volume, but Phase 2 saturation-boundary numerics — can add 1–2 weeks if the current scheme is fully explicit and a switch to implicit/preconditioned formulation is forced (decide this **before** Phase 1, since it affects the refactor architecture).

---

## 9. Recommendation

Build the validation suite (helium regression, water single-phase, Bennett two-phase) as an **automated test suite from day one**, re-run after every phase. Regime-switching logic makes it easy to fix one regime and silently break another.

**Sequencing:** (p,h) refactor + CoolProp closure → verify exact helium reproduction → HEM two-phase properties → boiling correlation chain → CHF checking (post-processed constraint before becoming a solver mode). Each step independently testable against the previous validated baseline.
