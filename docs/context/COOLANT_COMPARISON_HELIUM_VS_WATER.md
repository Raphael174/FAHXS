# Helium vs. Water Coolant Comparison: Why Q Barely Moves (And Where That Stops Being True)

## Scope of validity — read this first

**Everything in this document applies to the helical-coil-in-shell geometry
(`combustorProp.HX_config == "shellnHelicalTube"`, `main_solve.py`), where the
hot combustion gas is on the shell side around the coil.** The conclusion —
that swapping the coolant fluid barely changes absorbed duty `Q` — is a
consequence of *this specific geometry's* resistance stack, not a general
property of heat exchangers or of water vs. helium as fluids. See "Scope: Where
this does NOT generalize" below before applying any of this reasoning to the
shell-and-tube solver or to rocket nozzle regenerative cooling.

## 1. The observation

Running Helium (`single_phase_coolprop`) and Water (`equilibrium_liquid`) at
identical inlet temperature, inlet pressure, mass flow, and HX geometry, the
absorbed duty differs by only a few percent — the user's own solver runs
showed Helium ≈ 263 kW vs. Water ≈ 260 kW (≈1.2% gap, finite-rate chemistry).
An independent controlled diagnostic run (co-flow, frozen chemistry for speed,
`mass_flow_c = 0.15 kg/s`, `T_in = 303.15 K`, `p_in = 80 bar`) reproduced the
same pattern with a larger absolute gap due to a modeling asymmetry explained
below:

| | Q (kW) | mean T_c | mean driving ΔT (T_g − T_c) | R_gas (K/W) | R_wall (K/W) | R_coolant (K/W) | R_coolant / R_total | mean cp_c (J/kg·K) |
|---|---|---|---|---|---|---|---|---|
| **Helium** | 199.3 | 454.2 K | 2072.5 K | 14.88 | 0.156 | 0.056 | 0.37% | 5190 |
| **Water** | 192.1 | 474.8 K | 2082.7 K | 15.01 | 0.150 | 0.382 | 2.46% | 4721 |

The naive intuition — "helium has 24% higher cp than water, so it should carry
more heat" — is not wrong in direction, but it identifies a mechanism that
turns out to be nearly irrelevant in this geometry. The real story is a
resistance-network argument, plus one modeling artifact worth knowing about.

## 2. Root cause: the hot-gas film dominates the series resistance stack

Heat transfer from hot gas to coolant crosses three resistances in series:

```
R_total = R_gas + R_wall + R_coolant
UA = 1 / R_total
Q ≈ UA · (T_g − T_c)
```

In this geometry `R_gas` is 97.5–99.6% of `R_total` for both fluids (see
table above). The coolant film is the *smallest* link in the chain by two
orders of magnitude. This matches the existing `Biot_c`/`Biot_g` diagnostics
already tracked in `data_master` (`Biot_g ≈ 0.01–0.015`, `Biot_c ≈ 2–3` in
prior runs) — the gas side is the bottleneck.

### The general sensitivity law

Differentiating `UA = 1/(R_gas + R_wall + R_coolant)`:

```
d(ln UA) / d(ln R_coolant) = − R_coolant / R_total = − (R_coolant's fraction of total resistance)
```

This says: **the leverage that changing the coolant-side resistance has on Q
is capped by the coolant's own share of the total resistance stack.** With
`R_coolant/R_total` under 2.5% here, even a full fluid swap (a large,
non-infinitesimal change in `R_coolant`, going from 0.056 to 0.382 K/W — a
~7× change) can only move Q by a few percent. This is why cp and thermal
conductivity — which set `h_c` and hence `R_coolant` — have so little
purchase on Q in this specific configuration.

## 3. Two mechanisms that produce the actual (small) gap

Since the coolant film is nearly irrelevant to Q, the residual few-percent gap
comes from two places, not from a direct cp/conductivity tradeoff at the film:

### 3.1 A modeling asymmetry (not real physics): missing coil enhancement for water

The Helium path dispatches through the Mori–Nakayama **helical-coil** Nusselt
correlation (`combustorProp.Nusselt_coil`, default `"mori1967"`), which
includes a Dean-number secondary-flow enhancement from coil curvature. The
Water path uses `physics/liquid_flow/correlations.py`'s geometry-light
straight-pipe correlations (Gnielinski-type single-phase liquid Nusselt,
Gungor-Winterton boiling) — this module is explicitly documented as
"intentionally geometry-light... before wiring the same ideas into HX-specific
solvers," and the coil-curvature enhancement for the liquid path is not yet
implemented (see `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`).

This means Helium gets a coil-curvature HTC bonus that Water currently does
not, purely because the liquid closure hasn't had the same geometry-specific
enhancement work done yet — this is model incompleteness, not a physical
difference between the fluids. It accounts for the majority of the
`R_coolant` gap between the two runs (0.056 vs. 0.382 K/W). **A real coiled
boiling flow would also get a secondary-flow enhancement; the current gap
partly flatters Helium.**

### 3.2 Coolant temperature profile (where the cp intuition actually lands)

Helium heats gradually across the full coil length. Water heats to
saturation and then plateaus at `T_sat(p)` for the back ~40% of the coil
(latent heat absorption at ~constant temperature — see
`VARIABLES_REFERENCE.md` and the boiling-onset discussion elsewhere in
project docs). This makes Water run hotter on average (mean `T_c`: 474.8 K
vs. 454.2 K), which slightly reduces the driving `(T_g − T_c)`.

Helium's higher cp does contribute to keeping it cooler on average — but as
a fraction of the ~2075 K gas-to-coolant ΔT, the ~20 K mean-temperature shift
is only about 1%. **The cp effect is real, but it acts through the
temperature profile, not through the coolant film — and it's a weak (~1%)
effect relative to the gas-side-dominance argument in §2.**

## 4. How to run an honest comparison with the current codebase

### Tier 1 — configuration you must equalize

- Same `T_in`, `p_in`, `mass_flow_c` for both fluids (swap only `coolant`/`coolant_model`).
- Same `numericalProp.chemistry_model` for both runs (`finite_rate` for a real
  answer, `frozen` for fast iteration — never mixed; the gas-side duty
  depends on it and the HX is gas-side-limited, so a chemistry-model mismatch
  would swamp any coolant difference).
- Use `combustorProp.flow_config = "co"` for the cleanest comparison. Co-flow
  starts both fluids at the physical `(T_in, p_in)` with a pure forward march
  — no boundary guessing. Counter-flow is consistent between fluids now
  (both route through a physical-inlet shooting reference by default in
  `main_steady.py`), but co-flow remains the strictly simpler controlled case.
- Set `coolantProp.liquid_chf_lut_path` for the water run (default as of this
  session). Without it, `chf_margin_c` is `NaN` everywhere and the water run
  silently keeps using optimistic nucleate-boiling HTC even past CHF — a
  dishonestly favorable result for water.
- Verify `n_nodes` (length of `data_master["T_g"]`, etc.) matches between the
  two runs. If one run hits a pressure/enthalpy floor and stops early, you are
  integrating `Q` over different physical lengths.

### Tier 2 — neutralize the one real asymmetry

To compare the two fluids on equal footing *in the subcooled/single-phase
regime*, run Helium with the same straight-pipe correlation family the water
path uses:

```python
combustorProp(Nusselt_coil="Gnielinski", ...)   # instead of the default "mori1967"
```

**Caveat: do not try to equalize away Gungor-Winterton once water boils.**
That's a genuinely different (and better) heat-transfer mode — real physics,
not an artifact — so the two fluids should only be forced onto identical
correlations where both are actually in the same physical regime
(single-phase liquid/gas).

### Tier 3 — compare the metrics that actually respond to fluid choice

Because the coolant film is <2.5% of total resistance, **`Q` is the least
sensitive output to coolant choice** — reporting only `Q_tot` will always
look like "a tie" here almost by construction. An honest comparison should
report the outputs where the fluids actually diverge:

- **Wall temperature** (`T_wg`, `T_wc` peaks) — does the coolant keep the tube
  safe?
- **Coolant-side pressure drop / pumping power** (`dp_c_tot`; for water, now
  split into friction vs. acceleration contributions via `dp_c__dx_accel` —
  see `VARIABLES_REFERENCE.md`).
- **Outlet state**: Helium exit temperature vs. Water exit quality,
  void fraction, and **CHF margin** (`chf_margin_c`).

These are where the coolant choice is actually consequential in this design;
`Q_tot` is not.

## 5. Scope: where this does NOT generalize

The conclusion "coolant fluid choice barely matters" is a consequence of
`R_gas ≫ R_coolant` **in this specific geometry** (hot combustion gas in
shell-side cross-flow around a helical coil, with the coil-side coolant film
already fast by comparison). This does not carry over automatically:

### 5.1 Shell-and-tube (`main_solve_shellntube.py`)

Different geometry, different correlation set: hot gas is inside the tubes
(`Nusselt_tube = "gnielinski_blended"` by default, or the corrugated/grooved
tube correlations), and the coolant is in baffled shell-side cross-flow
(`Nusselt_shell_baffled = "bell_delaware"`). Tube confinement raises gas-side
velocity and Nu; baffled cross-flow also raises the coolant-side `h`
substantially relative to a plain annulus. The resistance ratio is not
guaranteed to be as lopsided as the helical-coil case — it should be checked
per design point, not assumed.

More importantly: **liquid/boiling coolant coupling is not yet wired into the
shell-and-tube steady solver's coupled march** — it is postprocess-only (see
`docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`). A live,
apples-to-apples water-vs-helium comparison in the same rigorous
coupled-march sense used above is not currently possible for shell-and-tube;
that requires the Phase 3 shell-side liquid coupling work (see
`docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md`) to land first.

### 5.2 Rocket nozzle regenerative cooling

Regenerative-cooled nozzle throats operate in a fundamentally different
regime: gas-side heat transfer coefficients are typically far higher there
than in this combustor's shell-side coil geometry (very high velocity,
frequently near/at the throat's characteristic Mach-1 condition, thin
boundary layers). In that regime the gas-side resistance is often *not* the
dominant term — the coolant side becomes comparably or even more limiting,
especially right at the throat. This is precisely the regime where the
resistance-fraction law in §2 flips: `R_coolant/R_total` grows, and coolant
fluid properties (this is a large part of why hydrogen is favored as a
regenerative coolant in many real engines — its heat transfer capability is
directly consequential to nozzle survivability, not a rounding error on Q).
**Do not import the "coolant choice barely matters" conclusion into any
nozzle-cooling analysis without first checking where `R_coolant/R_total`
actually falls for that geometry and heat flux.**

## 6. Summary

| | This combustor's helical-coil-in-shell design | Where this would flip |
|---|---|---|
| Dominant resistance | Hot gas (shell-side, ~97.5–99.6% of R_total) | Coolant-limited when R_gas is reduced (better gas-side correlation/geometry, or much higher heat flux as in a nozzle throat) |
| Sensitivity of Q to coolant fluid | ~R_coolant/R_total (a few %) | Approaches unity as R_coolant/R_total grows |
| What actually differs between fluids here | Wall temperature, pumping power, CHF margin/quality — not Q | Q itself, once coolant-limited |
| Known modeling gap | Water/liquid path lacks the coil-curvature Nusselt enhancement Helium gets (Phase 3 not done) | N/A — applies wherever the liquid closure is used pre-Phase-3 |
