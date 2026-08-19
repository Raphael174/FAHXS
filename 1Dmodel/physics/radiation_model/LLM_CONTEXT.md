# physics/radiation_model/ LLM Context

## Scope

Gas radiation (H2O/CO2) for the hot-gas-to-wall flux, via a tabulated Weighted-Sum-of-
Gray-Gases (WSGGM) model. Tabulated only — no per-node spectral computation — consistent
with the project-wide "radiation must stay tabulated" numerical lesson
(`docs/context/TRANSIENT_STATUS.md` / CLAUDE.md). Consumed by
`physics/heat_conduction.py`'s `OneDimensionalSteadyConduction_ShellnHelicalTube` as an
optional parallel-path addition to the hot-side convective coefficient.

## Contents

| File | Role |
|---|---|
| `radiation_equations.py` | `qrad_net_mbl(Tg, Ts, eps_g_emit, eps_g_abs, eps_s)` — net gas→surface radiative flux via a compact two-gray-surface form with a mean-beam-length-based gas emissivity/absorptivity pair; `hrad_from_q(Tg, Ts, qpp)` — converts that flux to an equivalent radiative HTC (with a small-ΔT linearized fallback to avoid divide-by-zero). |
| `radiation_build.py` | `Ehlme2025Coeffs` (loads/evaluates the polynomial-in-mixture-ratio WSGGM coefficient set from JSON), `WSGGM_Ehlme` (evaluates gas emissivity/absorptivity from T, partial pressures, and mean beam length), `RadiativeBackendEhlme`/`make_ehlme_backend()` — the adapter matching the solver's `rad_backend(T_eval=..., p=..., yH2O=..., yCO2=..., Le=..., **_)` call signature consumed by `heat_conduction.py`. |
| `radiation_model_WSGGM_parameters.py` | One-time offline script ("TO BE RUN ONCE") that tabulates the published Ehlmé et al. (2025) 4-gray-gas WSGGM coefficients and writes them to the three JSON files below. Not imported at solver runtime. |
| `ehlme2025_mixture.json`, `ehlme2025_pure_H2O.json`, `ehlme2025_pure_CO2.json` | Generated coefficient data consumed by `Ehlme2025Coeffs.from_json()` / `make_ehlme_backend()`. |

## Model details worth knowing

- Ehlmé et al. (2025) WSGGM, valid T = 500-5000 K, mixture ratio MR = Y_H2O/Y_CO2 in
  [0.4, 4.0], pressure-pathlength p·L in [0.01, 60] atm·m (coefficients fit at ~1 atm;
  validity stated in atm·m terms only). `WSGGM_Ehlme.emissivity()` clips MR into that
  band before evaluating and clips the final emissivity into `eps_clip` (default
  `(1e-9, 0.999999)`) — it will not raise on out-of-range inputs, it silently clamps.
- The mixture/pure branch selection in `WSGGM_Ehlme.emissivity()` falls back to the
  pure-H2O or pure-CO2 coefficient set when the other species' partial pressure is ~0,
  rather than evaluating the mixture polynomial at a degenerate MR.
- `radiation_equations.qrad_net_mbl` averages emissivity-at-Tg and absorptivity-at-Ts
  into a single effective gas grayness (`eps_g_eff`) for the two-surface radiation
  network — this is a simplification, not a full spectral two-temperature treatment.

## How this plugs into the wall solver

`heat_conduction.py`'s `OneDimensionalSteadyConduction_ShellnHelicalTube` takes
`rad_enabled`, `rad_backend` (a callable matching `RadiativeBackendEhlme.__call__`'s
signature), and `rad_state` (a dict with `p`, `yH2O`, `yCO2`, `Le`). When enabled, it
computes gas emissivity at `T_g` (emission) and at the current `T_wg` estimate
(absorption), converts to `h_g_rad` via `hrad_from_q`, and adds it in parallel with the
convective `h_g` (`h_g_eff = h_g + h_g_rad`). The helical solver builds this backend from
`ehlme2025_mixture.json` when `numericalProp.radiation_ON` is true, using H2O/CO2 mole
fractions from the Cantera gas state and a mean beam length driven by geometry and
`CorrelationCoefficients.mbl_factor` (per `docs/context/PHYSICS_CONTEXT.md`).

## Test coverage

No dedicated test file targets this folder by name; it is exercised indirectly wherever
`numericalProp.radiation_ON` is exercised in the baseline-regression tests. No direct
unit test for `WSGGM_Ehlme.emissivity()`'s numeric output against a published reference
point was found in `tests/`.

## TODOs

None found stated as TODO/FIXME comments in this folder's files.

## Change history

`docs/solver_design/FV_CORE_REWORK_PLAN.md` lists "WSGGM tabulated radiation" under
`physics/radiation_model/` as reused-unchanged infrastructure for the in-progress FV core
rework (Section 2, "What already exists and must be reused, not rewritten"). No other
dated history evidenced for this folder.
