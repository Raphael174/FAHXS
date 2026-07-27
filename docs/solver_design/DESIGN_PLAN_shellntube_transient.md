# Design Plan — Baffled Shell-and-Tube Config + Transient Solver

Current-status note: this file is a design history and architecture rationale.
The implementation has advanced beyond the original plan: shell-and-tube
transient now supports co-flow and counter-flow, finite-rate FPV is the default
chemistry path with cached manifolds, pre-ignition GOX chilldown is modeled with
CoolProp Oxygen, and production transient runs use a linearly-implicit
`fixed_step` wall integrator. For current truth, read
`docs/TECHNICAL_REFERENCE.md`, `docs/context/TRANSIENT_STATUS.md`, and
`CLAUDE.md`.

Status: **design approved for implementation** (written 2026-07-07, architecture pass).
Implementation target: work packages WP1→WP4 below, in order. Each WP is independently
testable and leaves `main_solve.py` (helical config) fully working.

---

## 0. Context

The existing solver (`main_solve.py`) is a steady quasi-1D marching solver for the
shell-and-helical-tube combustor: hot combustion gas on the shell side (Cantera),
Helium inside the coil (CoolProp), cylindrical wall conduction + WSGGM radiation,
per-node `fsolve` for wall temperatures.

Two features are added:

- **Feature A** — a TEMA-style baffled shell-and-tube HX config (single segmental
  baffles, tubes in windows — geometry mirrors the EchTherm reference case): hot
  combustion gas flows **axially inside straight tubes**, Helium flows on the
  **shell side in zig-zag crossflow around baffles**.
- **Feature B** — a transient (dynamic) solver for start-up and shut-down, reusing
  the steady per-node physics.

### Reference operating point (from EchTherm case + current input_data)

| Quantity | Value | Consequence |
|---|---|---|
| Tubes | ~235 × Ø5 mm OD, 0.75 mm wall (ID 3.5 mm), L ≈ 235 mm, Inconel 718 | L/D ≈ 67 → thermal entrance matters |
| Shell | ID 110 mm, 15 single-segmental baffles, 20 % cut, spacing ≈ 12 mm, triangular pitch p/D = 1.3, tubes in windows | Bell-Delaware applies directly |
| Tube-side gas | ~100 g/s total → 0.43 g/s/tube, T ≈ 2000→800 K, p ~ 1–5 bar | **Re ≈ 2 000–5 000 per tube: laminar-transitional!** Pr ≈ 0.7 |
| Shell-side He | 150 g/s, 90 bar, 120 K inlet | Re_Do ≈ 10⁵ (up to 10⁶), Pr ≈ 0.66, Mach ≪ 0.1 |
| Wall | Bi = h·t/k ≈ 0.05–0.2 | radially lumpable per node |
| Time scales | τ_wall = ρ c t / h ≈ 10–60 s vs gas residence ≈ 5 ms | quasi-steady fluid during transients |

The per-tube Re estimate is the single most important design fact: **the tube-side
correlation set must cover laminar + transitional + turbulent with smooth blending**,
not just a turbulent correlation.

---

## 1. Methodology assessment (answers to the open questions)

### 1.1 Is the quasi-1D compressible ("1D Navier-Stokes") formulation right for the tubes?

**Yes — reuse it for the hot gas inside the tubes.** The equations in
`physics/governing_equations.py` (dT/dx, dU/dx, dp/dx, dρ/dx) are fluid-agnostic
quasi-1D single-duct equations. Today they are applied to the He coil; in the new
config they are applied to a **single representative tube** (all 235 tubes are
identical in this 1D idealization): per-tube ṁ = ṁ_g / N_tubes, dQ per node
multiplied by N_tubes when handed to the shell side. Gas properties keep coming
from Cantera with the existing `remove_energy()` equilibrium/frozen mechanism.
Tube-side gas can reach Mach 0.1–0.3 at the hot end, so keeping the compressible
momentum/continuity terms is genuinely useful there (unlike for He).

Radiation inside a Ø3.5 mm tube is optically thin: mean beam length ≈ 0.9·D ≈ 3 mm
vs ~30+ mm in the combustor annulus → gas emissivity is negligible.
**Default `radiation_ON = False` for this config**; if enabled, use Le = 0.9·D_tube_inner.

### 1.2 What about a "1D Navier-Stokes" for the shell-side He?

**No — and this is a physics decision, not a shortcut.** The shell fluid's actual
path is a 3D zig-zag; there is no well-defined 1D streamwise coordinate along which
a momentum ODE is meaningful. The correct altitude is exactly what Bell-Delaware
was built for:

- **Energy: 1D in the axial coordinate.** The He bulk temperature T_s(x) advances
  axially node-by-node on the *same grid* as the tube-side march
  (dT_s/dx = ±dQ_total/dx / (ṁ_c · cp_c)). This makes counter/co-flow, the wall
  coupling, and later the transient solver all share one spatial grid.
- **Momentum: 0D empirical per baffle compartment.** Δp from the Bell-Delaware
  pressure-drop method (crossflow zones + window zones + end zones), distributed
  onto the axial grid pro-rata so p_s(x) is available at every node.
- **Density/velocity: algebraic.** He Mach is ≪ 0.1 everywhere; ρ = ρ(T,p) from
  CoolProp at each node, characteristic crossflow velocity from continuity through
  S_m (the Bell crossflow area). No dU/dx ODE needed.

### 1.3 Counter-flow handling must change

The current helical solver avoids the two-point boundary-value problem by starting
the He march from a *prescribed outlet state* (`T_out`, `p_out`) — fine for
calibration, but not predictive. For the new config implement a **predictive
sweep iteration** (both inlet states known, opposite ends):

1. Initialize T_s(x) as a linear profile between He inlet T and a guess outlet.
2. March the tube-side gas from its inlet, solving the wall + dQ per node against
   the current T_s(x) profile.
3. March the shell-side He from *its* inlet using the frozen dQ(x) field.
4. Under-relax (ω ≈ 0.5) the T_s(x) update; repeat 2–4 until max|ΔT_s| < 0.05 K.

This converges in ~5–15 sweeps for gas–gas HX (capacity ratio far from 1 helps).
Keep the existing prescribed-outlet mode available behind the existing
`flow_config` switch for the helical config; the new solver gets
`solve_mode = "sweep"`.

---

## 2. Feature A — correlation set

All new correlations go in `physics/heat_transfer_correlations.py` /
`physics/friction_correlations.py` with the same pattern as today: pure functions
+ a dispatcher + knobs in `CorrelationCoefficients`.

### 2.1 Tube side — combustion gas in straight circular tube (Re 500–20 000, Pr ≈ 0.7)

New dispatcher `dispatch_nu_tube_straight(selector, ...)`:

- **Laminar (Re < 2300)** — Gnielinski/VDI laminar entrance composite
  (VDI Heat Atlas, section G1), constant-T-wall form:
  `Nu = [3.66³ + 0.7³ + (1.615·(Re·Pr·D/L)^(1/3) − 0.7)³]^(1/3)`
  At L/D ≈ 67 and Re·Pr·D/L ≈ 20–25 this gives Nu ≈ 5–6, i.e. **+40–60 % over the
  fully-developed 3.66 — the entrance term is not optional here.** Use local x for
  L when marching (Leveque-type local Nu), with a floor at the fully-developed value.
- **Turbulent (Re > 4000)** — Gnielinski, already implemented
  (`compute_Nusselt_Gnielinski`), f from Colebrook (exists) or Konakov.
- **Transitional (2300 ≤ Re ≤ 4000)** — Gnielinski (2013) linear blending:
  `γ = (Re − 2300)/(4000 − 2300)`, `Nu = (1−γ)·Nu_lam(2300) + γ·Nu_turb(4000)`.
  Never leave a discontinuity — the solver marches through this band as the gas cools.
- **Variable-property correction** — for a *gas being cooled*, Kays & Crawford /
  McEligot give exponent ≈ 0: `Nu/Nu_cp = (T_w/T_b)^n`, `n_tube_gas = 0.0` default,
  exposed as a knob (the hot-end T_b/T_w ≈ 2.5 is far outside correlation databases;
  this is a prime calibration knob against test data, exactly like `Nusselt_correction`
  today).
- **Friction**: laminar 64/Re (Darcy), Colebrook for turbulent, same γ blending in
  transition. Add to `friction_correlations.py` as `dispatch_friction_tube_straight`.

### 2.2 Shell side — He in baffled crossflow (Re 10³–10⁶, Pr ≈ 0.66)

**Method: Bell-Delaware** (Shah & Sekulić, *Fundamentals of Heat Exchanger Design*,
ch. 9 has every equation and the curve-fit constants — implement from there; also
Taborek in HEDH). Structure:

`h_shell = h_ideal · Jc · Jl · Jb · Js · Jr`

- **h_ideal — ideal tube bank**: Žukauskas (1972) staggered-bank correlation,
  piecewise in Re (constants per Re decade, covers Re up to 2·10⁶),
  `Nu = C·Re^m·Pr^0.36·(Pr/Pr_w)^0.25`. Re based on tube OD and **max (gap)
  velocity** through S_m. Pr = 0.66 is marginally below the nominal 0.7 floor —
  acceptable; expose prefactor `zukauskas_C_factor` as calibration knob.
- **Jc** — baffle-window correction from fraction of tubes in crossflow F_c
  (function of baffle cut; "tubes in windows" layout → F_c from geometry).
- **Jl** — leakage (tube-to-baffle-hole area S_tb + shell-to-baffle area S_sb vs S_m).
  Clearances come straight from the EchTherm-style inputs (diametral clearances).
- **Jb** — bundle bypass (F_sbp = bypass area / S_m, sealing-strip pairs N_ss).
- **Js** — unequal inlet/outlet baffle spacing (front/rear end lengths).
- **Jr** — laminar correction, only Re < 100 → implement as no-op with a warning.
- **Δp**: Bell-Delaware analog — ideal-bank Euler number (Žukauskas or
  Gaddis–Gnielinski), zones: `Δp = Δp_cross·(N_b−1)·R_b·R_l + Δp_window·N_b·R_l
  + Δp_end·R_b·R_s`, plus optional nozzle K-factor losses.

All J/R factors are pure functions of geometry ratios + Re → put them in a new
`physics/bell_delaware.py` so they're unit-testable in isolation against the
worked example in Shah & Sekulić (there is a complete numerical example in the
book — use it as the acceptance test).

### 2.3 Geometry module

New `mechanical/geometry/shelltube_geometry.py`:

- Inputs mirror the EchTherm GEOMETRY screen (see `shellTubeProp` below).
- Outputs: S_m (crossflow area at bundle centerline), S_w (window flow area),
  S_tb, S_sb (leakage areas), F_sbp (bypass fraction), F_c (fraction of tubes in
  crossflow), N_tcc / N_tcw (tube rows crossed per crossflow / window zone),
  D_otl (outer tube limit), window hydraulic diameter, per-compartment → axial-grid
  mapping, and tube-count sanity check from shell ID + pitch + layout.

---

## 3. Feature A — solver architecture

New class `shellntube_solver` in a new file `main_solve_shellntube.py`
(do **not** graft into `main_solver`'s while-loop; share the physics modules).

Per outer sweep-iteration, per axial node i (grid = tube axial coordinate, uniform dx):

1. Tube-side (representative tube): Cantera props → Re, Pr → `dispatch_nu_tube_straight`
   → h_t; friction f_t.
2. Shell-side: CoolProp He props at (T_s(x_i), p_s(x_i)) → Bell-Delaware h_s
   (h_s is per-compartment piecewise-constant; map compartment → node range).
3. Wall: reuse `OneDimensionalSteadyConduction_ShellnHelicalTube` — it is already
   a generic cylindrical-wall gas/coolant node. Hot fluid is now *inside*, so the
   radial direction flips: hot side = inner perimeter π·D_in, cold side = outer
   π·D_out. **Add a `hot_side` = `"inner"`|`"outer"` switch to the class rather
   than duplicating it** (the log-term and area assignments swap; ~15-line change,
   keep default `"outer"` so the helical config is untouched).
4. dQ_node (per tube) → ×N_tubes → into shell-side energy balance;
   tube-side advance via existing `governing_equations` + `remove_energy()`.
5. Stress: `stress_pressure_tube` is internal-pressure hoop; here pressure is
   *external* (90 bar shell, ~2 bar tube) → tube is in **external-pressure
   collapse**, a different failure mode. Add `stress_external_pressure_tube` +
   a simple collapse-pressure check (Roark) in `mechanical/loads.py`.

Convergence + `_check_global` energy balance identical in spirit to today.

### New input dataclasses (`input_data.py`)

```python
@dataclass
class shellTubeProp:
    # tubes
    N_tubes: int = 235
    D_tube_outer: float = 5e-3
    thickness_tube_wall: float = 0.75e-3
    L_tube: float = 235e-3
    layout: str = "triangular30"     # "triangular30" | "square90" | "rotated45"
    pitch_ratio: float = 1.3
    # shell
    D_shell_inner: float = 110e-3
    # baffles (single segmental, tubes in windows)
    N_baffles: int = 15
    baffle_cut: float = 0.20          # window opening fraction of D_shell
    baffle_thickness: float = 3e-3
    L_front_end: float = 100e-3
    L_rear_end: float = 10e-3
    # clearances (diametral) + bypass
    clearance_tube_baffle: float = 1e-3
    clearance_baffle_shell: float = 1e-3
    clearance_bundle_shell: float = 0e-3
    N_sealing_strip_pairs: int = 0
    # allocation + materials + correlation selectors
    tube_side_fluid: str = "hotgas"   # "hotgas" | "coolant" (implement hotgas first)
    Nusselt_tube: str = "gnielinski_blended"
    Nusselt_shell_baffled: str = "bell_delaware"
    material_tube: str = "INCO718"
    material_shell: str = "INCO718"
```

`CorrelationCoefficients` additions: `n_tube_gas = 0.0`, `zukauskas_C_factor = 1.0`,
`bell_Jl_factor = 1.0`, `bell_Jb_factor = 1.0`, `tube_laminar_Nu_fd = 3.66`,
`Re_transition_lo = 2300`, `Re_transition_hi = 4000`.

### Validation (WP2 acceptance)

1. **Bell-Delaware unit tests** vs the Shah & Sekulić worked example (J factors,
   S_m, Δp to within round-off of the book values).
2. **ε-NTU cross-check**: constant-property run (frozen chemistry, constant cp)
   must reproduce the analytic counterflow ε-NTU curve within ~2 %.
3. **EchTherm reproduction**: run the exact screenshot case; compare Q, outlet T's,
   both Δp's. Expect ±20 % (Bell-Delaware's honest accuracy) — document the deltas.
4. Energy balance ΣdQ vs ṁΔh on both streams < 1 %.

---

## 4. Feature B — transient solver

### 4.1 Architecture decision: quasi-steady fluid, lumped wall-energy ODE

**Third and final correction (2026-07-07, numerically verified with an honest
three-way comparison — this supersedes BOTH earlier drafts in this section).**
The architecture is: **one lumped wall-energy ODE state per axial node
(`dT̄_w,i/dt`), with the two face temperatures reconstructed each step from a
quadratic quasi-static profile.** This is the *original* plan's approach — it
was wrongly "falsified" in the second draft by a strawman test, and is here
reinstated with proper validation.

**What went wrong in the second draft.** That pass compared the true wall PDE
against a **zero-thermal-mass steady resistance network** — i.e. a model where
the wall temperature *teleports* to its new steady value the instant `h_c`
changes, with no thermal inertia at all. That model is indeed off by ~440 K
during the ramp — but it was never the proposed architecture. The proposed
architecture always integrated the wall's stored energy
(`ρc·δ·dT̄/dt = q_hot − q_cold`); only the *radial profile shape* (not the
energy) was taken quasi-static. Testing inertia-less teleportation and
concluding "we must resolve the wall thickness" was a category error.

**The honest three-way test** (script: `doc/check_wall_quasi_static_validity.py`,
grid/timestep-converged; INCO718, δ=0.85 mm, `h_g`=300 W/m²K, `h_c` ramped
300→75 000 W/m²K over 2.5 s):

| Model | max face-T error vs PDE truth | captures transient wall-ΔT peak? |
|---|---|---|
| **A** — resolved radial PDE (truth) | — (reference) | yes: ΔT peaks **68.9 K at t≈0.34 s**, 5× the 13.8 K steady value |
| **B** — lumped energy ODE + quadratic profile reconstruction | **1.8 K** (flux err ~1–2 %) | **yes — 69.0 K, essentially exact** |
| **C** — zero-mass steady network (the strawman) | 438 K | no — misses it entirely |

**Conclusion: the lumped-ODE + quadratic-reconstruction model (B) tracks the
fully-resolved PDE to <2 K through the entire fast ramp, including the
transient wall-ΔT overshoot that peaks at ~69 K** (the thermal-shock signal a
transient solver exists to catch). No radial sub-grid is needed for thin
coil/tube walls. The second draft's radial-FD mandate is withdrawn.

**The wall-energy ODE per axial node:**

```
ρ_w c_p,w δ · dT̄_w,i/dt = q_hot,i − q_cold,i
    q_hot,i  = h_g,eff·(T_g,i − T_wg,i)      (radiation folded into h_g,eff as today)
    q_cold,i = h_c,i·(T_wc,i − T_c,i)
```

`T̄_w,i` (thickness-mean wall temperature) is the ONLY time-integrated state
per node. The two face temperatures needed for the fluxes and for stress
post-processing are reconstructed each step from a **quadratic quasi-static
profile** consistent with both Robin BCs and the current mean — justified
because the wall's radial *profile-shape* relaxation (τ_diff,δ = δ²/(π²α) ≈
13 ms) is fast compared with the ramp, even though its *stored energy*
(τ_wall = ρc·δ/(h_g+h_c), see table) is not.

Profile: `T(x) = T_wg − (q_hot/k)x + (q_hot−q_cold)·x²/(2kδ)`, with
`q_hot = h_g,eff(T_g−T_wg)`, `q_cold = h_c(T_wc−T_c)`. Two constraints —
(i) thickness-mean equals `T̄_w`, (ii) the cold-face flux-balance identity —
give a **2×2 linear solve for `(T_wg, T_wc)`** each step. With radiation folded
into `h_g,eff` the system is weakly nonlinear (`h_g,eff` depends on `T_wg`
through `q_rad`), so wrap the 2×2 in **2–3 fixed-point iterations**.

⚠️ **Implementation warning (learned the hard way this session):** the 2×2
face-reconstruction algebra is easy to get subtly wrong — a sign error in the
cold-face constraint produced physically plausible-looking but diverging
garbage (wall ΔT growing to >1000 K) before it was caught. **Ship a
steady-state self-consistency assertion**: at fixed `h`, the reconstructed
`(T_wg, T_wc)` must reproduce the steady 3-resistance-network faces to machine
precision. This test is already in the guardrail script and passes for the
corrected algebra.

**Fallback clause (when lumped is NOT valid).** Lumped-B is validated for
*thin* walls (coil 0.85 mm, tube 0.75 mm). It would break if a future config
has a thick thermal mass — e.g. if the 8 mm shell wall or a thick tube-sheet
is given structural thermal inertia, δ²/(π²α) ≈ 11 s and the profile-shape is
no longer quasi-static — or for sub-100 ms ramps. **Guardrail:** the audit
(§4.2) runs the check script's B-vs-A comparison for the active geometry/ramp
at t=0 and, if B's error exceeds a threshold (say 5 K or 5 %), switches *that
component's* wall to a small implicit radial sub-grid (the second-draft
approach, kept in reserve, not deleted). Thin walls never trigger it.

**Time scales, corrected** (τ_wall = ρc·δ/(h_g+h_c) is not a fixed "~10 s" —
it *sweeps three orders of magnitude during the ramp*; verified numerically):

| Scale | Value | Role |
|---|---|---|
| τ_residence, hot gas | ~1–5 ms | gas refreshes the duct — quasi-steady OK |
| τ_diff,δ (profile shape) | ~13 ms (δ²/π²α) | radial profile relaxes — justifies quasi-static *shape* |
| **τ_wall(t)** (stored energy) | **5.05 s → 0.04 s** as `h_c` goes 300→75 000 | the wall settles *during* the ramp; this is the integrated ODE |
| τ_ramp (He 0→full) | 2–3 s | dominant transient of interest |
| τ_residence, He (see §4.6) | 10.6 s → 0.11 s as ṁ goes 1 %→100 % | **exceeds τ_wall early in ramp → cold-side quasi-steady-fluid caveat, §4.6** |

- `τ_residence,gas ≪ τ_wall` throughout ⇒ hot-side quasi-steady-fluid safe.
- `τ_diff,δ ≪ τ_ramp` ⇒ quasi-static *profile shape* safe (validated: B≈A to 2 K).
- The wall's *stored energy* is NOT quasi-static — it is the integrated state,
  correctly handled by the ODE. (This is the distinction both earlier drafts
  blurred.)

**Required refactor (WP3a), the only invasive change:** factor the *flux
evaluation* out of the steady conduction node so it can be called in two modes:
the existing steady `fsolve` (heat-in = heat-out, unchanged), and a new
**`fluxes_at_Tbar(T̄_w, h_g, h_c, T_g, T_c)`** that does the 2×2 face
reconstruction and returns `(q_hot, q_cold, T_wg, T_wc)` for the transient
driver's ODE right-hand side. Same class, added method — not a separate N-node
discretization. The steady solver path is byte-for-byte unchanged.

### 4.2 Time integration

- State vector: `T̄_w,i` (one lumped wall-energy value per axial node, §4.1).
- RHS(t, `T̄_w`): one full quasi-steady fluid pass (the sweep of §1.3) using the
  reconstructed face temperatures as the wall BC → per-node `h_g,i(t)`,
  `h_c,i(t)`, `q_hot,i`, `q_cold,i` → `dT̄_w,i/dt = (q_hot,i − q_cold,i)/(ρc·δ)_i`.
  No wall `fsolve` in the transient path — the 2×2 face reconstruction replaces it.
- Integrator: the wall ODE is non-stiff once the fluid is quasi-steady, so
  **adaptive `scipy.solve_ivp` (`RK45`, or `BDF` if the fast-ramp segments prove
  stiff)** over the `T̄_w` vector, with `max_step` capped by the fastest-changing
  active schedule. Do **not** default to a fixed `dt` sized off the final-state
  `τ_wall` — τ_wall sweeps 5.05 s → 0.04 s during the ramp (§4.1 table), so a
  `dt` chosen from the settled value undersizes the early ramp by ~2 orders of
  magnitude. Adaptive control handles this automatically; a fixed-`dt` RK4
  fallback must size `dt` off `min_t τ_wall(t)`, not the endpoint.
- Print a **time-scale audit** at t=0: τ_residence (gas and He), τ_diff,δ,
  and `τ_wall(t)` **swept across the schedule's `h_c` range** (not a single
  value — the sweep is the whole point). Warn if (a) `max_step` exceeds ~1/10
  of the fastest schedule's local `|X/(dX/dt)|`, (b) `τ_residence,gas/τ_wall >
  1/100` anywhere (hot-side quasi-steady-fluid), or (c) `τ_residence,He > τ_wall`
  during any ramp segment (cold-side early-ramp caveat, §4.6 — flag He-outlet
  temperature as unreliable in that window rather than silently trusting it).

### 4.3 Boundary-condition schedules

```python
@dataclass
class transientProp:
    t_end: float = 120.0
    max_step: float = 0.05             # solve_ivp cap; audit tightens per schedule
    solver_method: str = "RK45"        # "RK45" | "BDF" (if fast segments stiff)
    # each schedule: list[(t, value)] linearly interpolated; None => hold steady value
    schedule_mass_flow_g: list = ...   # start-up ramp / shutdown cutoff
    schedule_mass_flow_c: list = ...   # e.g. [(0,1e-3),(2.5,0.15)] = 2.5 s He ramp
    schedule_OF: list = ...
    schedule_T_c_in: list = ...
    schedule_p_c_in: list = ...
    T_wall_initial: float = 293.15     # cold-start wall temperature (uniform)
    chemistry_transient: str = "frozen"  # see 4.4
    ignition_time: float = 0.0           # gas side inert (air/no flow) before this
    flag_He_outlet_when_residence_gt_tau: bool = True  # §4.6 cold-side caveat
```

The wall ODE needs `(ρ_w·c_p,w·δ)` per node → wall **density and specific heat**
must be available vs temperature (§4.7). Radial thermal mass uses the tube/coil
wall only (thin, validated lumped); shell/tube-sheet mass is out of scope for
v1 (would trip the §4.1 fallback — flag, don't silently lump).

Shut-down scenario = ramp ṁ_g → 0 at t_cut with He still flowing (the thermal-soak
/ thermal-shock case that motivates the feature). ṁ_g = 0 nodes: h_t → natural-
convection floor (small constant), radiation off.

### 4.4 Chemistry cost control

Equilibrium Cantera per node per RHS call ≈ 10⁵–10⁶ equilibrate calls per run —
unacceptable. Default for transients: **`"frozen"`** (composition fixed at the
inlet equilibrium of the *current* O/F, h(T) from fixed composition — cheap and
consistent). Optional `"requilibrate_every_N"` for slow ramps. The inlet
combustion state itself is recomputed only when the O/F or ṁ schedule value
changes beyond a tolerance.

### 4.5 Config-agnostic by construction — BOTH configs ship transient (required)

The transient driver (`main_solve_transient.py`) owns the `T̄_w(x, t)` state
vector and calls a config-provided
`fluid_pass(T̄_w_field) → (q_hot(x), q_cold(x), h_g(x), h_c(x), diagnostics)`.
**Deliver `fluid_pass` for BOTH heat-exchanger configs — this is a hard
requirement, not an optional extension:**

- **shell-and-helical-tube** (existing): `fluid_pass_helical` is a light
  extraction of the body of `main_solver.solver()`'s march loop — everything
  from the property lookups through `dispatch_nu_*` and the flux computation,
  but reading wall faces from the reconstruction instead of the steady `fsolve`.
  Start-up/shut-down of the *current* combustor then works immediately.
- **shell-and-tube** (WP2): `fluid_pass_shelltube` wraps that solver's per-node
  physics identically.

Both share the same wall-ODE driver, integrator, audit, schedules, and
dashboard — only `fluid_pass` differs. Acceptance requires a passing transient
run for *each* config (the helical one is also the cheapest regression check
that the extraction didn't perturb the steady physics).

### 4.6 He-side inventory & early-ramp caveat (quantified, must be flagged at runtime)

Two distinct He storage effects, with very different magnitudes:

1. **Pressure/charging transient** (mass to pressurize the shell/coil void):
   coil void ≈ 1.10 L, He ≈ 14.4 kg/m³ at 90 bar/300 K → inventory ≈ 16 g;
   charging at ≥ tens of g/s is sub-second → negligible vs τ_wall. Quasi-steady
   pressure is fine. (This is the only effect the earlier draft considered.)
2. **Thermal transport lag (the one that actually bites):** the He *residence
   time* = inventory/ṁ = **10.6 s at 1 % flow, 1.06 s at 10 %, 0.11 s at full
   flow**. Compare with τ_wall = 5.05 s → 0.04 s over the same ramp (§4.1): for
   roughly **the first ~10–20 % of the He ramp, residence time > τ_wall**, so
   the "fluid is instantaneously steady" assumption fails on the *cold* side in
   that window. The heat actually removed there is small (low ṁ, low h_c), so
   the wall/gas trajectory is barely affected — but the **predicted He-outlet
   temperature is not trustworthy** during that early window. The runtime audit
   (§4.2, warning c) must flag it; do not silently report He-outlet T there.

### Validation (WP3 acceptance)

1. **Consistency**: transient run with constant BCs from a perturbed initial wall
   field must converge (t → ∞) to the steady solver's solution, node-by-node < 0.5 K.
2. **Wall-model cross-check (the §4.1 guardrail script IS this test — run it,
   don't just cite it)**: the lumped-ODE + quadratic-reconstruction wall model
   (B) must match the fully-resolved radial PDE (A) to **< 5 K face temperature
   and < 5 % flux** across a ramp matching the actual `transientProp` schedule,
   AND the steady-state self-consistency assertion (reconstructed faces = steady
   network faces) must pass to machine precision. Current status: B vs A = 1.8 K
   on the design-point ramp — passes with margin. This gate also decides whether
   the §4.1 fallback (radial sub-grid) is needed for a given config.
3. **Both configs**: a completed transient run for helical AND shell-and-tube
   (§4.5), each with the energy audit below closing.
4. **Global energy conservation in time**: ∫(Q_hot − Q_cold)dt = ΔU_wall < 1 %,
   where ΔU_wall = Σ_i (ρc·δ·A)_i·ΔT̄_w,i (exact for the lumped model).
5. `dt`/`max_step` convergence study spanning the *fastest* scheduled ramp
   segment, not just the overall run length.

### 4.7 Temperature-dependent material properties (needed for the wall ODE)

The transient wall ODE integrates `(ρ_w c_p,w δ)·dT̄/dt`, so **wall density and
specific heat** now matter (they were irrelevant to the steady solver, which
only used k, CTE, E, yield). Current state in
`mechanical/material_specs/material_temperature_strength.py`
(`init_material_temperature_properties` returns `CTE, E, Yield, Lambda, density,
poisson`):

| Property | 316L | INCO718 | Status for transient |
|---|---|---|---|
| k(T) conductivity | ✅ interp table | ✅ interp table | ready — wall ODE uses it |
| CTE(T), E(T), Yield(T) | ✅ interp | ✅ interp | ready — stress post-proc |
| density | ⚠️ **constant scalar** (7.9e3 / 8.2e3) | ⚠️ constant scalar | acceptable — metals' ρ varies <2 % to 800 °C; keep constant, note it |
| **c_p(T) specific heat** | ❌ **absent** | ❌ **absent** | **MISSING — must add for the wall ODE** |

**Plan for the gap:**
- Add `c_p,wall(T)` functions to the material module following the exact
  interp1d-with-flat-extrapolation pattern already used for k/CTE/E (so they
  slot into `init_material_temperature_properties`'s return tuple as a new
  element). Use published tables: **316L** c_p ≈ 500 J/kg·K at RT rising to
  ~600 by 800 °C (ASM / AISI Designer's Handbook, same source as the existing
  316L data); **INCO718** c_p ≈ 435 J/kg·K RT rising to ~630 by 800 °C
  (Special Metals INCO718 datasheet). **Flag every prescribed/approximated
  value with a source comment**, matching the house style already in that file.
- Where a curve is genuinely unavailable, prescribe a constant at the
  operating-temperature midpoint and emit a one-line WARNING at solver init
  ("c_p,wall for <material> prescribed constant = X J/kg·K — no T-table"), so
  approximations are visible, never silent.
- Density: keep the constant scalar; add a code comment stating the <2 %
  justification so it's a documented choice, not an oversight.
- Gas-side (Cantera) and He-side (CoolProp) properties are already fully
  temperature/pressure-resolved — no gap there.

Return-signature change: `init_material_temperature_properties` gains a `Cp`
callable. Update its two call sites in `main_solve.py` (`self.func_..._HX`
unpacking) — additive, the steady path can ignore `Cp`.

---

## 4c. Result dashboard — extend to a dynamic (time-resolved) HTML dashboard

The existing `model_data_process/data_plotting.py::HXDashboard` produces five
static matplotlib figures (thermal, helium, combustion, mechanical, radiation)
for a *steady* run. Keep it for steady runs. For transient runs add a
**self-contained dynamic HTML dashboard** (`data_plotting_transient.py`,
`TransientDashboard`) — a single file the user opens in a browser, no server.

**Design (follow the `dataviz` skill before writing chart code):**
- **Data source**: the transient driver logs, per saved time step, the same
  per-node quantities the steady `make_solver_data` dict holds (T_g, T_c, T_wg,
  T_wc, dQ, Nu, Re, stresses, radiation flux, gas composition/props from
  Cantera, He props from CoolProp) plus the new `T̄_w(x,t)` field. Serialize to
  JSON embedded in the HTML (self-contained; no external fetch — matches the
  Artifact CSP constraints if ever published).
- **Views** (time as the interactive axis — a scrubber/slider + play control):
  1. **Axial profiles vs time** — T_g, T_c, T_wg, T_wc along x, animating with
     the time slider; overlay the steady solution as a reference ghost curve.
  2. **Wall temperature heatmap** — `T̄_w(x, t)` as an x–t image (the headline
     transient view; peak thermal-shock ΔT is read straight off it).
  3. **Outlet/scalar histories** — He-outlet T, Q_tot, max wall ΔT, max
     stress/yield vs t as line charts, with the §4.6 He-outlet-unreliable
     window shaded.
  4. **Combustion/thermo panel** — inlet T_g, O/F, key species mass fractions,
     cp_g, γ vs t (from the finite-rate/equilibrium chemistry track).
  5. **Material/mechanical panel** — per-time max stress vs temperature-derated
     yield, with a margin indicator; k(T), c_p(T) operating points annotated.
- **Tech**: inline vanilla JS + a single lightweight embedded plotting approach
  (hand-rolled SVG/Canvas or an inlined minimal lib) — **no CDN**, everything in
  one HTML file. Theme-aware (light/dark) per the dataviz palette guidance.
- Provide a `TransientDashboard(time_series_data).to_html(path)` entry point and
  a `.all()`-style convenience mirroring the steady `HXDashboard` API, so the
  transient `main_solve_transient.py` closes with one call, symmetric with how
  `HX_sizing_brief` closes the steady run.

The steady `HXDashboard` is untouched; the dynamic one is additive and reuses
the same `data_master` key names so a transient run is just "a stack of steady
snapshots + the wall-temperature time field."

---

## 4b. Feature C — chemistry via an FGM/FPV flamelet manifold (tabulated)

**Architecture chosen (2026-07-07, user direction): a Flamelet-Generated-Manifold /
Flamelet-Progress-Variable table, NOT per-node Cantera calls.** This supersedes
the earlier "per-node CoolingPFR" sketch (kept below as the trajectory generator
only, not the runtime path).

**Why FGM/FPV.** For this combustor the operating regime (<100 g/s diesel/O2 at
O/F≈2, ~440 kW extracted by the He → ~4.4 MJ/kg enthalpy removed) makes the choice
of chemistry model first-order: frozen badly under-predicts (misses recombination
heat release during the deep cooldown), HP-equilibrium is the minimum acceptable
default, and at low hot-gas flow even equilibrium breaks (residence ~ chemical
time) so finite-rate is required. An FGM/FPV manifold captures all three as slices
of ONE table and stays fast (runtime = interpolation, no Cantera in the march):

- Controlling variables: mixture fraction **Z**, progress variable **C**, enthalpy
  **h** (the enthalpy dimension is mandatory — this is a strongly non-adiabatic
  HX, unlike a classic adiabatic FGM).
- **Equilibrium** ≡ the C = C_eq(h) locus; **frozen** ≡ C held at inlet;
  **finite-rate** ≡ C transported along x by its source term ω_C(C,h,Z) with
  residence time dt = dx/U.
- The cooling channel sits at a single bulk mixture fraction Z̄ = 1/(1+O/F) ≈ 0.33,
  so the *channel* interpolates a 2-D manifold (C, h) at fixed Z̄. The full Z
  dimension is only needed for the mixing/flame inlet zone (regime 1).

`flamelet_kit` is the manifold generator: regime-1 `Flamelet` builds the Z-manifold
(non-premixed diesel/O2 flame), regime-2 `CoolingPFR` traces the (C, h) cooling
trajectories at Z̄. Importable directly (`from flamelet_kit import Flamelet,
CoolingPFR, SteadyCache` — the bootstrap puts the repo parent on sys.path). Read
`flamelet_kit/METHODOLOGY.md` + `ADAPTATION_GUIDE.md` first.

### C0 — equilibrium manifold (the C=C_eq(h) slice) — IMPLEMENTED

Layer 1, shipped in `main_solve_transient.py`: at fixed inlet O/F and p, the
equilibrium cooling path is 1-D in enthalpy-removed. Precompute it once (one
equilibrium sweep, few hundred `equilibrate` calls) → tabulate T, ρ, μ, k, cp,
yH2O, yCO2 vs h_removed; the march interpolates (no per-node Cantera). Radiation
becomes a 2-D table ε(T_eval, h_removed) since composition varies along the path.
Rebuild the table (SteadyCache-gated) only when inlet O/F or p move beyond
tolerance. This makes the *required* equilibrium default as fast as frozen
(~minutes/100 s) and IS the equilibrium edge of the full manifold — directly
extended by adding the C axis in C1.

### C1 — full (Z̄, C, h) FGM/FPV manifold + progress-variable transport (finite-rate)

**Physics being captured (why low flow needs this).** Da = τ_res/τ_chem. Near the
hot inlet, recombination is fast (τ_chem tiny) → local equilibrium. As the gas
cools (and at low hot-gas flow the enthalpy removed *per unit mass* is large, so
it cools hard and fast through the ~1600→1000 K window), τ_chem for
CO+½O₂→CO₂ and H+OH recombination blows up → the composition FREEZES at a
partially-recombined state. Equilibrium wrongly keeps releasing recombination
heat below the freeze point; frozen wrongly releases none from the start. The
FPV progress variable C, transported with a source ω_C that self-extinguishes as
T falls, reproduces the freeze-out and lands between the two bounds — the
physically correct answer.

**Manifold coordinates.** Fixed channel mixture fraction Z̄ = 1/(1+O/F) (single
already-mixed stream). Progress variable from a monotone recombination combo:

    C = (Yc − Yc_inlet) / (Yc_eq(h) − Yc_inlet),   Yc = Y_CO2 + Y_H2O − Y_CO

normalized so C=0 at the inlet burnt mix and C→1 at the local-h equilibrium.
Enthalpy h is the second axis (non-adiabatic).

**Table generation (offline, once per (O/F, p) — SteadyCache-gated).** Reuse the
`flamelet_kit` reactor idiom (REPRODUCE_SPEC.md Part B), NOT a per-node PFR at
runtime:
1. For a grid of enthalpy levels h_j spanning inlet→cold (same span as the C0
   equilibrium manifold): set a constant-pressure reactor to (h_j, p, Y_inlet)
   and integrate it in time, sampling (C, T, Y, all props) and the instantaneous
   dC/dt from the reactor's net production rates as it relaxes from the inlet
   composition toward equilibrium at that fixed h. This traces the C-path and its
   rate at that h_j.
2. Assemble onto a regular (C, h) grid: state(C,h) = {T, ρ, μ, k, cp, yH2O, yCO2}
   and ω_C(C,h) [1/s]. (Off-manifold (C,h) cells that no trajectory reaches are
   filled by nearest-valid extrapolation and flagged.)
3. Radiation ε(T_eval, C, h) — but composition depends on (C,h), so this is the
   C0 2-D ε(T_eval, h) generalized to 3-D; in practice ε is weak in C at fixed h,
   so tabulate ε(T_eval, h) at the *marched* C and accept the small error, or add
   the C axis if validation demands.

**Runtime (one extra marched scalar in `fluid_pass`).** Track C alongside h_removed:
    dC/dx = ω_C(C, h) / U_g ,   h advances from the wall flux as in C0.
Interpolate all gas props + ω_C from the (C,h) table each node — still zero
Cantera calls in the march, so speed stays at the C0 level (target <10 min/100 s
for kinetic too). Limits recovered exactly: ω_C≡0 ⇒ frozen; ω_C→∞ ⇒ C=C_eq(h) ⇒
the C0 equilibrium manifold. Inlet C set from the flame/mixing zone (regime-1
`Flamelet`, or C=0 for a fully-burnt injector feed).

**Validation.** (a) frozen and equilibrium limits reproduce the C0/frozen tables;
(b) a direct per-node `CoolingPFR.march` at one operating point matches the
tabulated-FPV channel outlet T and composition to a few %; (c) sweep hot-gas mass
flow down and confirm the freeze-out (outlet composition departs from equilibrium
below a Da threshold) appears as physically expected.

### C2 — parameter-range tabulation over (O/F, p)  [user request, do as a last step]

Extend the C0/C1 manifolds from a single (O/F, p) point to a grid so O/F and p
can vary during a transient without a rebuild, and to support the ignition
handoff (C3):
- Grid: **O/F ∈ [1.5, 3.5]** (e.g. 5 levels), **p ∈ [1 atm … 5 bar]** (e.g. 4
  levels). Equilibrium manifold becomes 3-D: state(O/F, p, h). Kinetic FPV
  becomes 4-D: state + ω_C over (O/F, p, C, h).
- Build cost is offline and one-time (a few hundred equilibrate/reactor calls per
  (O/F,p) node × 20 nodes) — parallelizable; SteadyCache still gates re-use.
- Runtime: the march reads O/F(t), p(t) from the BC schedules and quad-/quintic-
  linearly interpolates the manifold. Still zero Cantera calls in the march.
- Guard: warn/clamp if a scheduled O/F or p leaves the tabulated box.

### C3 — ignition = pilot diesel flame (no PLA mechanism)  [user decision, 2026-07-08]

**Resolved simply.** The real igniter is PLA/O₂, but rather than wire a PLA
mechanism the user represents the ignition flame as a **small pilot diesel/O₂
flame** — a low-propellant-mass-flow segment of the diesel schedule. Consequences:
- **One diesel/O₂ manifold covers the entire 100 s run** (ignition pilot → full
  injection). No fuel switching, no PLA surrogate, no second manifold, no blend.
- Because the manifold is **per-unit-mass (h in J/kg)**, the propellant mass-flow
  ramp needs NO rebuild — the pilot and full-flow states are the same table,
  reached at different spatial rates.
- The whole transient is driven purely by the **mass-flow schedules** already
  implemented: `schedule_mass_flow_g` (propellants: pilot → full diesel/O₂) and
  `schedule_mass_flow_c` (He coolant). Set the ignition window as the low-ṁ_g
  segment; ramp to full at diesel-injection time.
- Only a change in **O/F or pressure** across the schedule needs the C2 (O/F, p)
  grid; if O/F and p are ~constant, the single manifold already suffices.
- The wall ODE and He side are unaffected.

### C1-legacy (per-node PFR) — trajectory generator only, NOT the runtime path

### C1 — per-node finite-rate channel chemistry (satisfies the existing TODO(finite_rate))

Do **not** call `CoolingPFR.march()` as a one-shot (it wants T_wall(x) and constant
h_conv up front — fights the per-node wall fsolve). Instead implement a
`FiniteRateNode` (new `physics/combustion_chemistry/finite_rate_node.py`) that
embeds the CoolingPFR idiom per node, with a **persistent**
`IdealGasConstPressureReactor` + `ReactorNet` across the march:

1. chemistry sub-step: `net.advance(dt)`, `dt = dx / U_g` (adiabatic, const p);
2. heat-extraction sub-step: `gas.HPY = h − dQ_node/ṁ, p_updated, Y` (+ `r.syncState()`),
   exact by construction — mirrors both `remove_energy()` and `cooling_pfr.py`.

Behind `numericalProp.chemistry_model = "finite_rate"`; works in the helical config
immediately and in the shell-and-tube config unchanged. Loosen CVODE tolerances
from the kit's 1e-9/1e-15 to ~1e-6/1e-12 for the in-loop path after a
tolerance-convergence check.

**First deliverable of C1 (before any integration): the three-way study** —
resurrect the commented-out equilibrium/frozen/finite-rate comparison in
`combustion_gas.py` (`finite_rate_hx_solve` skeleton is already there) as a
standalone study script. It quantifies the payoff at design point, fixes the
Da-window bounds for C3, and measures per-call cost of the mechanism.
**Mechanism size is the main risk**: count species in RenKokjohn_surrogate.yaml
first; if per-node CVODE cost is prohibitive, source a reduced n-alkane surrogate
mechanism (see ADAPTATION_GUIDE §1) for the finite-rate path only.

### C2 — flamelet inlet for the mixing zone

New adapter `physics/combustion_chemistry/flamelet_inlet.py` wrapping
`flamelet_kit.Flamelet`, replacing the HP-equilibrium inlet of
`combustion_gas_solve.solve()` when `chemistry_model = "finite_rate"`. Valid model:
diesel/O2 injector is genuinely non-premixed two-stream (kit applicability §8 OK).
Three correctness points:

1. **Hand-off state = flamelet mixed over Z, not T_at_Z(Z_st).** Bulk O/F = 3 →
   Z̄ = 1/(1+OF) = 0.25 (near Z_st ≈ 0.23 for diesel/O2). Channel inlet (h, Y) =
   complete-mixing adiabatic average of the flamelet profile over Z (mass-weighted
   by the Z-PDF; complete-mixing/delta-PDF at Z̄ is the v1 default, β-PDF optional
   later).
2. **Enthalpy bookkeeping**: the current `dH_tot` heuristic (LOX gasification +
   fuel Hv) must move into the flamelet's Dirichlet boundary stream enthalpies
   (effective T_fuel/T_ox), or stay as a post-mix correction — never both, never
   dropped. Verify total-enthalpy consistency vs the equilibrium inlet at the
   frozen limit.
3. **χ_st is a physical calibration knob** (strain → flame T → CO/H2 slip → in-tube
   recombination heat release). v1 estimate: χ_st ~ U_inj / L_mix (mixing_length =
   50 mm); expose `chi_st_factor` in `CorrelationCoefficients`.

Steady solve: march the flamelet to steady once per run. Transients: gate re-solves
with `SteadyCache` as shipped (slow monotone ramps are its favorable regime).

### C3 — chemistry cost control (in payoff order)

1. **Composition lagging across sweep iterations**: after sweep k, if
   max|Y(x)_k − Y(x)_{k−1}| < tol, freeze the composition field and skip reactor
   advances in later sweeps (energy/wall balance still re-solved). Trivially correct
   at convergence; biggest steady-solver win.
2. **`SteadyCache` gating** of flamelet-inlet re-solves and (in transients)
   full-channel chemistry re-advances (re-advance every N steps or on gate trip).
3. **Da-window hybrid switching**: HP-equilibrium above ~1600 K, finite-rate through
   the freeze-out window, frozen below ~1000 K — thresholds from the C1 study
   diagnostics (local heat-release relaxation time vs residence time), not hardcoded.
4. **NN/ISAT surrogate: deferred.** Keep a chemistry-provider seam (the
   `FiniteRateNode` interface) so one can slot in later if profiling justifies it.

### Validation (WP-C acceptance)

1. Finite-rate limits: χ→∞ chemistry-off reproduces "frozen"; tightened-dt
   finite-rate at very hot end approaches "equilibrium" trend (C1 study plots).
2. Energy balance per node closes to reactor tolerance (the HPY split guarantees
   it — assert it).
3. Helical-config regression: `chemistry_model = "equilibrium"` results bit-identical
   after the refactor.

## 5. Work packages (implementation order)

| WP | Content | Touches | Risk |
|---|---|---|---|
| **WP1** | `bell_delaware.py`, `shelltube_geometry.py`, tube-side blended-Nu + friction dispatchers, `shellTubeProp`, `CorrelationCoefficients` knobs, unit tests vs book example | new files + additive edits | low |
| **WP2** | `main_solve_shellntube.py` steady sweep solver, `hot_side` switch in conduction class, external-pressure stress, ε-NTU + EchTherm validation, dashboard variant | conduction class (small), loads.py | medium |
| **WP3a** | `fluxes_at_Tbar()` method on conduction node (2×2 quadratic face reconstruction + steady-consistency assertion) + `fluid_pass` extraction for **both** configs (§4.5) + `c_p,wall(T)` added to material module (§4.7) | main_solve.py loop body, material_specs, heat_conduction.py | **highest** — do behind a flag, regression-check the helical baseline bit-identical before/after |
| **WP3b** | `main_solve_transient.py` (lumped wall-ODE state, `solve_ivp` adaptive driver, schedules, `transientProp`), time-scale audit (§4.2), **both configs' `fluid_pass`** | new files | medium |
| **WP3c** | `TransientDashboard` dynamic self-contained HTML (§4c) — axial-profile animation, T̄_w(x,t) heatmap, scalar histories, combustion + material panels | new file, additive | medium (follow `dataviz` skill) |
| **WP4** | Validation suite as `tests/` (Bell example, ε-NTU, wall-model guardrail script, both-config transient runs, energy audits) + doc updates (CLAUDE.md tables) | tests, docs | low |
| **WP-C1** | Three-way chemistry study + `FiniteRateNode` per-node PFR behind `chemistry_model="finite_rate"` (§4b) | combustion_chemistry, additive | medium (mechanism cost unknown — measure first) |
| **WP-C2** | Flamelet inlet adapter (mixed-over-Z hand-off, enthalpy bookkeeping, χ_st knob) | combustion_chemistry, additive | medium |
| **WP-C3** | Chemistry cost control: composition lagging, SteadyCache gating, Da-window hybrid | solver loops | low–medium |

The WP-C track is parallel to WP1/WP2 (chemistry plugs into both configs through
the same `remove_energy` seam). Suggested interleaving: WP-C1's study script early
(it de-risks everything downstream and needs no other WP), full C1 integration after
WP2, C2/C3 after WP3 so transient gating is testable.

Also note: `numericalProp.chemistry_model = "finite_rate"` is currently documented
as "not yet implemented (stub)" in CLAUDE.md — WP-C1 closes that.

Ground rules for the implementer:
- Never modify helical-config numerical behavior; every WP ends with a run of
  `python main_solve.py` reproducing today's summary numbers.
- Follow the existing house style: pure correlation functions with docstrings citing
  source + validity range + DOI, dispatchers with `error_factor`/`corrCoeffs`
  passthrough, per-node data recorded via `data_master` key matching.
- Every new correlation gets its calibration knob in `CorrelationCoefficients`
  with the literature default and an identifiability note.

## 6. Open decisions (flagged, defaults chosen)

1. **Fluid allocation**: plan assumes hot gas in tubes / He in shell (matches the
   stated intent and the EchTherm case). `tube_side_fluid` switch reserved; the
   reversed allocation is WP-later.
2. **Is the shell-tube HX fed by the combustor?** Assumed yes (inlet state =
   combustor exit via Cantera). A standalone-inlet mode (prescribed T/p/composition)
   is a trivial addition — worth doing in WP2 for EchTherm comparison runs.
3. **Baffle orientation vs gravity, condensation, fouling**: out of scope.
4. **Tube-side per-tube maldistribution** (plenum feeding 235 tubes): out of scope
   for 1D; note as model limitation in docs.
