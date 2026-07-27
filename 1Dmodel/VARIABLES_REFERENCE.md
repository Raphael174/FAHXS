# Steady helical solver — per-node variable reference

Every field listed here is a key of `solver.data_master` (the dict returned by
`make_solver_data()` in `model_data_process/data_processing.py`), populated
once per axial march node by `main_solver` (`main_solve.py`). Each value is a
1D array of length `n_nodes` (the number of coolant-channel arc-length steps
actually marched before the loop stopped — see `_coolant_flow_continues()`).

Two coolant modes populate different subsets of these fields:

- **Gas mode** (`coolantProp.coolant_model = "single_phase_coolprop"`, e.g.
  Helium): treats the coolant as a real gas via CoolProp `PropsSI` calls,
  state = `(T, p)`.
- **Liquid/boiling mode** (`coolantProp.coolant_model = "equilibrium_liquid"`,
  e.g. Water): state = `(p, h)` (pressure, specific enthalpy) — never
  `(T, p)`, because `T` and `p` are not independent inside the two-phase
  dome. Fields tied to the ideal-gas closure specifically (compressibility
  factor `Z`, Dean/Helical numbers, coil-correlation Nusselt number/friction
  factor, `cv`/`gamma`) are explicitly set to `NaN` in this mode rather than
  computed from a relation that does not apply to a liquid/two-phase fluid.
  `Mach_c`/`c_c` (sound speed) ARE populated in liquid mode too, from a real
  fluid EOS closure (Wood's equation inside the two-phase dome, CoolProp's
  real-EOS `SPEED_OF_SOUND` outside it) — not the ideal-gas relation used for
  `c_g`.

`solver.HXDashboard` (`model_data_process/data_plotting.py`) plots most of
these; see its `boiling()` method for the liquid-only quantities.

## Compressibility / coolant gas-only diagnostics

| Field | Meaning | Units | Mode |
|---|---|---|---|
| `Z` | Coolant real-gas compressibility factor, `PropsSI('Z', ...)` | `-` | Gas only (`NaN` in liquid mode) |
| `Mach_c` | Coolant bulk Mach number, `U_c / c_c` | `-` | Both |
| `c_c` | Coolant speed of sound. Gas mode: ideal-gas relation from `cp/cv`. Liquid mode: Wood's equation homogeneous-mixture sound speed (`two_phase_sound_speed`) for `0 <= quality <= 1` (collapses sharply right after boiling onset — can drop to a few hundred m/s or lower, well below either pure-phase value), CoolProp's real-EOS `SPEED_OF_SOUND` outside the dome. Blended with the two-phase value over the same boiling-onset window as `h_c`/`dp_c__dx` (see `dispatch.py`). | m/s | Both |
| `gamma_c` | Coolant specific heat ratio `cp_c/cv_c` | `-` | Gas only (`NaN` in liquid mode) |
| `cv_c` | Coolant specific heat at constant volume | J/(kg·K) | Gas only (`NaN` in liquid mode) |
| `De` | Dean number (coil curvature effect on flow), `Re_c * sqrt(Dh/D_coil)` | `-` | Gas only (`NaN` in liquid mode) |
| `He` | Helical (Germano) number, `Re_c * sqrt(Dh/(2*Rc))` | `-` | Gas only (`NaN` in liquid mode) |
| `Nu_c` | Coolant-side Nusselt number (coil correlation, e.g. Mori & Nakayama) | `-` | Gas only (`NaN` in liquid mode) |
| `f_c` | Coolant-side Darcy friction factor (coil correlation, e.g. Ali et al. 2024) | `-` | Gas only (`NaN` in liquid mode) |
| `f_fd_c` | Legacy fully-developed friction factor field — **never populated in the current solver** (only set in the archived `main_solve_CHT.py`); always an empty series in `data_master` for both modes. Not a bug, just a dead field carried over in `_SOLVER_DATA_KEYS`. | `-` | Neither (dead field) |

## Coolant thermodynamic state

| Field | Meaning | Units | Mode |
|---|---|---|---|
| `T_c` | Coolant bulk temperature. Gas mode: an independent state variable. Liquid mode: **derived** from the `(p, h)` state via the equilibrium closure — not itself a march variable. | K | Both |
| `p_c` | Coolant pressure | Pa | Both |
| `rho_c` | Coolant density | kg/m³ | Both |
| `U_c` | Coolant bulk velocity, `mass_flow_c / (rho_c * A_ch * N_ch)` | m/s | Both |
| `Re_c` | Coolant Reynolds number | `-` | Both |
| `cp_c` | Coolant specific heat at constant pressure | J/(kg·K) | Both |
| `mu_c` | Coolant dynamic viscosity | Pa·s | Both |
| `k_c` | Coolant thermal conductivity | W/(m·K) | Both |
| `Pr_c` | Coolant Prandtl number, `cp_c * mu_c / k_c` | `-` | Both |
| `h_c` | Coolant-side convective heat transfer coefficient. Gas mode: `Nu_c * k_c / Dh_ch`. Liquid mode: the boiling HTC closure output (e.g. Gungor-Winterton), which already blends nucleate + convective boiling contributions. | W/(m²·K) | Both |
| `enthalpy_c` | Coolant specific enthalpy — the primary state variable in liquid mode (not tracked/meaningful in gas mode) | J/kg | Liquid only |
| `dh_c__dx` | Axial gradient of coolant specific enthalpy | J/(kg·m) | Liquid only |
| `quality_c` | Thermodynamic vapor quality of the coolant. `< 0`: subcooled liquid. `0–1`: two-phase (boiling). `>= 1`: fully vaporized (march stops — no post-dryout closure past this point, see `_coolant_flow_continues()`). | `-` | Liquid only |
| `void_c` | Vapor void fraction (volumetric vapor fraction, via the two-phase closure's void-fraction correlation) | `-` (0–1) | Liquid only |
| `chf_margin_c` | Margin to critical heat flux, `q''_CHF / q''`. `> 1`: safe. `< 1`: CHF exceeded (dryout risk). `NaN` if `coolantProp.liquid_chf_lut_path` was not supplied (no CHF closure available). | `-` | Liquid only |
| `dU_c__dx`, `dT_c__dx`, `drho_c__dx` | Axial gradients of coolant velocity/temperature/density from the ideal-gas 1D momentum-energy-continuity closure | per-m units of the quantity | Gas only (`NaN` in liquid mode — the liquid march advances `(p, h)`, not these) |
| `dp_c__dx` | Axial coolant pressure gradient. Gas mode: from the ideal-gas closure. Liquid mode: `-(friction) + dp_c__dx_accel` — friction from Müller-Steinhagen & Heck (or the single-phase Darcy correlation outside the dome), plus the accelerational (HEM) contribution. | Pa/m | Both (different closures) |
| `dp_c__dx_accel` | Accelerational (HEM) contribution to `dp_c__dx`, liquid mode only — `-G^2 * d(1/rho)/dz` evaluated in the flow direction, using a one-node-lagged density gradient (same lagged-closure pattern as the boiling HTC's heat-flux term; see `main_solve.py`'s liquid `_advance_state()` block). Diagnostic split-out: `dp_c__dx` already includes it, this field lets you compare the friction vs. acceleration shares. Can dominate `dp_c__dx` right at boiling onset (observed up to ~90% of the total in testing) as density collapses through the two-phase dome — see `docs/solver_design/water_coolant_conversion_plan.md` section 4. | Pa/m | Liquid only |

## Hot gas state

| Field | Meaning | Units |
|---|---|---|
| `T_g` | Hot gas bulk temperature | K |
| `p_g` | Hot gas pressure | Pa |
| `rho_g` | Hot gas density | kg/m³ |
| `U_g` | Hot gas bulk velocity, `mass_flow_g / (rho_g * Ap_cc)` | m/s |
| `Mach_g` | Hot gas Mach number, `U_g / c_g` | `-` |
| `c_g` | Hot gas speed of sound (ideal-gas relation) | m/s |
| `W_g` | Hot gas mean molecular weight | kg/kmol |
| `cp_g`, `cv_g` | Hot gas specific heats | J/(kg·K) |
| `gamma_g` | Hot gas specific heat ratio | `-` |
| `mu_g` | Hot gas dynamic viscosity | Pa·s |
| `k_g` | Hot gas thermal conductivity | W/(m·K) |
| `Re_g` | Hot gas Reynolds number (tube/passage-based) | `-` |
| `Re_sh` | Shell-side Reynolds number (hot gas around the coil, uses `Dh_ch + 2*wall thickness`) | `-` |
| `Pr_g` | Hot gas Prandtl number | `-` |
| `Nu_g` | Shell-side Nusselt number (e.g. Salimpour 2008 correlation) | `-` |
| `h_g_conv` | Hot gas convective heat transfer coefficient | W/(m²·K) |
| `h_g_rad` | Hot gas radiative heat transfer coefficient (linearized equivalent) | W/(m²·K) |
| `emissivity_g` | Gas emissivity (WSGGM radiation model, evaluated at `T_g`) | `-` |
| `absorptivity_g` | Gas absorptivity (WSGGM radiation model, evaluated at `T_wg`) | `-` |
| `X_CO2`, `X_H2O` | Mole fractions of CO2 / H2O in the hot gas (radiation-relevant species) | `-` |
| `dp_g__dx` | Hot gas axial pressure gradient (Darcy friction), `-f_g * rho_g * U_g^2 / (2*Dh_cc)` | Pa/m |

## Heat transfer / wall conduction (both modes)

| Field | Meaning | Units |
|---|---|---|
| `dQ` | Heat transferred through this node (hot gas → coolant), from the 1D conduction solve | W |
| `dh_g` | Hot gas specific enthalpy drop at this node, `dQ / mass_flow_g` | J/kg |
| `dq__dx` | Heat flux per unit length, `(T_g - T_c) * UP` | W/m |
| `q_w` | Wall heat flux (total, conv + rad) | W/m² |
| `q_w_rad` | Radiative component of the wall heat flux | W/m² |
| `UP` | Conductance per unit length (series resistance of gas-side film, wall conduction, coolant-side film, all per unit dx) | W/(m·K) |
| `UA` | Overall conductance for this node (`UP` evaluated at `dx=1`, i.e. the true node conductance) | W/K |
| `Res_g`, `Res_c`, `Res_w` | Individual thermal resistances: hot-gas film, coolant film, wall conduction | K/W |
| `Biot_g` | `Res_w / Res_g` — wall conduction resistance relative to gas-side film resistance | `-` |
| `Biot_c` | `Res_w / Res_c` — wall conduction resistance relative to coolant-side film resistance | `-` |
| `T_wg` | Hot-side (gas-facing) wall temperature | K |
| `T_wc` | Cold-side (coolant-facing) wall temperature | K |
| `k_w` | Wall material thermal conductivity, evaluated at the mean wall temperature | W/(m·K) |
| `T_c_check` | Coolant temperature back-computed from the converged wall energy balance (`T_wc_new - dQ/(h_c*A_cold)`) — a self-consistency check against the coolant-side `T_c`/state used as input to that node's conduction solve, not an independent state variable | K |
| `phi_multiplier` | Tight-coil Nusselt correction multiplier — only populated when `combustorProp.Nusselt_shell == "churchill_bernstein_tightcoil"`; empty/NaN otherwise | `-` |

## Geometry

| Field | Meaning | Units |
|---|---|---|
| `L_HX` | Cumulative axial heat-exchanger length (projection of the coolant channel onto the combustor axis) | m |
| `L_ch` | Cumulative coolant channel (arc) length along the coil | m |

## Mechanical (coil tube stress, both modes)

| Field | Meaning | Units |
|---|---|---|
| `CTE` | Coil material coefficient of thermal expansion, evaluated at mean wall temperature | 1/K |
| `Modulus` | Coil material elastic (Young's) modulus, evaluated at mean wall temperature | Pa |
| `Yield` | Coil material 0.2% yield strength, evaluated at `T_wg` (conservative, default) or mean wall temperature depending on `numericalProp.yield_at_hot_wall` | Pa |
| `stress_pressure` | Hoop stress from internal (coolant) pressure | Pa |
| `stress_thermal_inner`, `stress_thermal_outer` | Thermal (radial temperature gradient) stress at the tube inner/outer wall | Pa |
| `stress_inner`, `stress_outer` | Total stress at inner/outer wall, `stress_thermal + stress_pressure` | Pa |

## Notes

- All units are SI per project convention (see `AGENTS.md`/`CLAUDE.md`); temperatures are in Kelvin throughout `data_master` (converted to °C only in plotting, e.g. `HXDashboard`).
- For the liquid mode's excluded gas-only fields (`Z`, `gamma_c`, `cv_c`, `f_c`, `Nu_c`, `De`, `He`), the `NaN` is intentional — see `main_solve.py`'s liquid branch (`main_solver.__init__`/`_advance_state`) for the exact assignment. `Mach_c`/`c_c` are NOT in this list — they are populated with real values in both modes (see above). If you see a finite, non-NaN value for one of the still-excluded fields on a liquid run, that indicates `coolantProp.coolant_model` was not actually `"equilibrium_liquid"` for that run (check for input drift — see `docs/USER_GUIDE.md` "Known Gotchas").
- `LiquidMarchSanityReport` (`physics/liquid_flow/sanity_checks.py`, `check_liquid_march()`) includes a `mach_c_max`/`mach_choking_ok` gate (hard fail at `Mach_c >= 1.0`, informational warning above `0.5`) alongside the existing energy-balance/CHF/dryout gates, printed at the end of a liquid steady run (`solver._check_global()`).
- The shell-and-tube solver (`main_solve_shellntube.py`) tracks a related but distinct per-tube field set (`solver.tube` dict) and is not covered by this document; shell-side liquid coupling is not yet wired there (postprocess-only — see `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`).
