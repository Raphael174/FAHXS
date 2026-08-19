# core/hotgas/ LLM Context

## Scope

Hot-gas providers consumed by `core/residual.py` (Stage D, not built yet) and
`core/wall.py`. Two very different kinds of provider live under this
package's eventual design: `combustor.py` (wraps the EXISTING Cantera/FPV
chamber state — subsonic combustor, not built yet) and `nozzle_gas.py`
(genuinely NEW quasi-1D supersonic expansion — built). Do not confuse the
two: shell-and-tube/helical (Stages A-E) never need `nozzle_gas.py`; only
Stage F (rocket-nozzle regen configs) does.

## Contents

| File | Status | Role |
|---|---|---|
| `nozzle_gas.py` | Built (Stage F groundwork, 2026-07-31, standalone/early — zero dependency on `core/mesh.py`/`residual.py`, which didn't exist yet) | Constant-gamma frozen quasi-1D area-Mach expansion, adiabatic wall temperature, Bartz-Cornelisse film-property HTC. 20 unit tests in `tests/test_nozzle_hotgas.py`, all passing, including exact closed-form isentropic checks. |
| `combustor.py` | **Does not exist yet** | Per the design doc §3, meant to wrap the existing chamber Cantera/FPV provider (`physics/combustion_chemistry/`) behind the same provider interface `nozzle_gas.py` will expose, so `residual.py` doesn't need to know which config it's driving. |
| `prescribed.py` | **Does not exist yet** | Table-driven `h_g(x)`/`T_aw(x)`/`p(x)` override path (user-supplied CEA/RPA output), through the same provider interface. Cheap cross-check + unblocks users who already have external nozzle output. |
| `__init__.py` | — | Empty/minimal package marker |

## `nozzle_gas.py` — what's genuinely new here (no legacy equivalent exists
anywhere in the repo — verified by search when this was built: the only
prior "nozzle" hits were shell-and-tube inlet/outlet pipe connections and
`combustorProp.exhaust_diameter`)

- `area_mach_ratio`/`mach_from_area_ratio`: closed-form isentropic A/A* <->
  Mach, both branches (subsonic upstream of throat, supersonic downstream).
- `choked_mass_flux`/`throat_diameter_for_mass_flow`: `mdot = G* · A_t`,
  `G*` fixed by chamber `T0`/`gamma`/`p0`. **Use this to check
  throat-diameter/chamber-pressure/mass-flow consistency before trusting a
  user-given geometry** — caught a real inconsistency in the first design
  point exercised through this code (user's stated 120 mm throat / 50 bar /
  45 kg/s only chokes ~30.5 kg/s together; resolved by deriving throat
  diameter from mass flow, the normally-authoritative thrust/Isp-driven
  quantity, giving 145.8 mm instead — see `docs/solver_design/
  FV_CORE_REWORK_PLAN.md`'s 2026-07-31 note for the full story).
- `adiabatic_wall_temperature`: `T_aw`, NOT `T_0` and NOT `T_static`, is the
  driving temperature for hot-side flux — flagged in this module's own
  docstring as "the most common modelling error in regen cooling", with an
  asserting unit test.
- `bartz_cornelisse_htc`: the **Cornelisse (1979) Eq. 8.3-3** form
  specifically — local free-stream properties + film-property ratio
  `(rho_f/rho)(mu_f/mu)` — chosen by the user over the RPE
  (Sutton & Biblarz) stagnation-referenced form. **A second, materially
  different Bartz form (RPE Eq. 8-22: `(mu_am/mu_0)^0.2`, referenced to
  STAGNATION viscosity, not local) exists in the literature and is NOT
  interchangeable line-for-line** — confirmed from two independently
  rendered textbook page images (PDF text extraction in this repo has
  silently dropped minus signs three separate times, so anything
  correlation-coefficient-critical here was visually re-verified, not
  transcribed from OCR). If a second Bartz variant is ever added, register
  it as a SEPARATE closure, do not silently merge the two forms into a
  hybrid — see the design doc §5.2 for the full source-verification
  discipline this correlation went through.
- `NozzleGasStation`/`solve_frozen_expansion`: per-station M/T/p/h_g/T_aw
  along the contour, frozen chemistry only (constant chamber-composition
  gamma) — the validation default. An `equilibrium` (shifting-composition,
  tabulated vs. `p/p_c`, cached like the FPV manifolds) mode is designed for
  but not yet implemented; per-node Cantera calls are forbidden on this path
  regardless (CLAUDE.md / `docs/context/TRANSIENT_STATUS.md`).

## TODOs (from the design doc, not invented)

- `combustor.py`, `prescribed.py`: not started.
- No wall/coolant coupling exists yet for the nozzle configs — every `q_w`
  number this module can currently produce uses an ASSUMED uniform wall
  temperature, not a solved one. A real answer needs `core/wall.py`'s rib
  path wired up (already built) plus Stage D/E's coupled residual extended
  to this config (not started). Do not quote `nozzle_c2h4_o2_bartz_example.py`
  output as a converged design number.
- The "chamber-stagnation + c* + area-ratio (A_t/A)^0.9 + throat-curvature
  (D_t/R_curv)^0.1 + sigma correction" Bartz parameterization — the most
  commonly-cited textbook form — is explicitly NOT implemented here and
  needs its own independent source verification before it is (design doc
  §5.2 point 2); do not assume it from memory even though it's the most
  commonly quoted "the Bartz equation".
- Cryogenic/copper wall material data (CuCrZr/NARloy-Z, typical regen
  nozzle materials) is not in the tree at all; existing `ST316L` tables
  clamp flat below 27°C, unusable for an LN2-cooled nozzle wall as-is. Not
  this package's problem to fix, but any wall-temperature number downstream
  of `nozzle_gas.py` for a cryogenic case inherits this gap.

## Change history

Built 2026-07-31 in one session alongside the first real nozzle design point
exercise (C2H4/O2, O/F=2.3, 50 bar chamber) — see
`validation/nozzle_c2h4_o2_bartz_example.py` and the design doc's dated notes
for the numbers produced and their caveats. No changes since. No deep git
history available (single initial commit + ongoing uncommitted project work).
