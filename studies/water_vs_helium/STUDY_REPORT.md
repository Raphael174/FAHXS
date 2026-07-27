# Water vs. Helium: Helical (6 m coil) vs. Shell-and-Tube — Full Study Report

## 0. Scope and what changed to make this possible

This study required a real capability that did not exist before: shell-and-tube
did not have `(p,h)` phase-change physics wired into its shell-side coupled
march (Phase 3 of `docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md` was
not started). Before this study could run honestly, the following solver
work was done (see `1Dmodel/main_solve_shellntube.py`):

- **Shell-side liquid coupling (Phase 3)**: `_shell_side_march` now advances
  coolant `(p, h)` (enthalpy-based) instead of `T` when
  `coolantProp.coolant_model == "equilibrium_liquid"`. The previous
  temperature-only march had **zero latent-heat accounting** — an early test
  showed "water" heating from 303 K straight past 1000 K with no boiling
  plateau at all, which is not physical. `_shell_h_at` now uses the same
  validated Gungor-Winterton/CHF closure the helical solver uses
  (`evaluate_coolant_closure`) once inside the two-phase dome, and a `(p,h)`
  CoolProp flash (not `(T,p)`) for single-phase branches either side of it —
  Bell-Delaware itself has no boiling model, so it is only used outside the
  dome.
- **Verified**: exact energy balance (`Q_from_wall_duty` vs.
  `mdot·Δh`, relative imbalance ~1.5e-16 — machine precision), smooth quality
  progression through the dome, sane wall temperatures (down from an
  unphysical ~1359 K to ~1145 K in a direct comparison at the same
  conditions). Full liquid + baseline regression suite: 28/32 passed; the 4
  failures are the pre-existing Helium baseline mismatch described in §0.1,
  confirmed unrelated to this work.
- **Known limitation carried over**: shell-side cross-flow enhancement of
  boiling HTC is real physics that Gungor-Winterton (a tube-flow correlation)
  does not capture. This is a documented simplification, not an oversight.
- **EOS-ceiling safety net added**: at very low coolant flow relative to duty,
  superheated steam can be driven past CoolProp's Water-backend validity
  ceiling (~3000 K). The march now freezes state at the last EOS-valid point
  instead of crashing (see §3.3 for where this actually triggered).

### 0.1 A note on the regression suite

While validating this work, both the shell-and-tube **and** the helical
Helium baseline regression tests were found failing (`tests/test_steady_baseline_regression.py`).
This was investigated rigorously: `_shell_h_at`, `_shell_side_march`, and
`_tube_side_march` were each compared byte-for-byte against hand-reconstructed
originals (isolated calls and a full hand-rolled sweep loop), and the file
was fully reverted and re-tested — the failure reproduced identically on the
**unmodified original file**. This confirms the regression is pre-existing
`input_data.py` (`combustorProp`/`hotgasProp`) default drift, unrelated to
this session's work. Not fixed here (out of scope for this study); flagged
for separate follow-up.

## 1. Methodology

| | Helical | Shell-and-tube |
|---|---|---|
| Geometry | `combustorProp()` defaults, `HX_config="shellnHelicalTube"` | `shellTubeProp()` defaults — confirmed to match `docs/context/shell_and_tube_architecture_target.png` exactly (N_tubes=235, D_tube=5mm, L_tube=235mm, D_shell=110mm, corrugated tube inside, all baffle/nozzle/corrugation-geometry fields checked) |
| Coil/tube length | `numericalProp.L_HX_max` calibrated (bisection) to **6.000 m** total coil arc length (`L_ch_max`) — the as-shipped default geometry gives ~7.69 m of coil, not 6 m, since `L_HX_max` is the *axial* projection, not the arc length | Fixed at 235 mm (target architecture) |
| Flow config | Co-flow (both fluids start at the physical `T_in`/`p_in`, no boundary-guessing) |
| Coolant inlet | `T_in = 303.15 K`, `p_in = 80 bar`, identical for both fluids |
| Coolant mass flow | 50, 100, 150 g/s (swept) |
| Hot gas | `hotgasProp()` defaults, unmodified, for Part 1 |
| Chemistry | `finite_rate` (production model, per instruction — accepted the runtime cost) |
| CHF | `liquid_chf_lut_path` = bundled Groeneveld 2006 LUT for all water runs |

12 base runs (2 geometries × 2 fluids × 3 mass flows). All raw data, plots,
and this report live under `studies/water_vs_helium/`:

- `raw_data/*.zip` — one per run, each containing the full solver data array
  (`.npz`), a JSON summary, and the run's own plot.
- `plots/*.png` — per-run 4-panel plots (wall temperature, coolant
  temperature, pressures, thermal resistance breakdown) plus
  `comparison_summary.png` (cross-run comparison).
- `steam_design_study/` — the secondary feasibility study (§4).
- `logs/` — raw solver stdout per run (sweep convergence history, warnings).

## 2. Results — full grid

| Geometry | Fluid | ṁ_c [g/s] | Q_tot [kW] | T_wg,max [K] | T_c,out [K] | dp_c [bar] | x_max | min CHF margin |
|---|---|---|---|---|---|---|---|---|
| Helical | Helium | 50 | 182.1 | 1088.6 | 1003.9 | 1.55 | — | — |
| Helical | Helium | 100 | 204.9 | 765.4 | 696.2 | 3.88 | — | — |
| Helical | Helium | 150 | 214.9 | 638.1 | 576.7 | 6.97 | — | — |
| Helical | Water | 50 | 184.0 | 1207.5 | 943.5 | 0.16 | 1.73 | **0.0022** |
| Helical | Water | 100 | 206.1 | 675.0 | 568.0 | 0.18 | 0.61 | 2.50 |
| Helical | Water | 150 | 209.8 | 644.0 | 568.1 | 0.11 | 0.15 | 5.25 |
| Shell-and-tube | Helium | 50 | 361.7 | 1719.1 | 1697.1 | — | — | — |
| Shell-and-tube | Helium | 100 | 453.0 | 1289.1 | 1176.0 | — | — | — |
| Shell-and-tube | Helium | 150 | 499.8 | 1162.7 | 945.2 | — | — | — |
| Shell-and-tube | Water | 50 | 323.6* | 2034.9* | 2007.5* | — | 3.67* | 0.275 |
| Shell-and-tube | Water | 100 | 449.6 | 1400.7 | 1276.9 | — | 2.29 | **0.054** |
| Shell-and-tube | Water | 150 | 523.2 | 1264.3 | 864.4 | — | 1.59 | **0.029** |

`*` = did not fully converge (oscillating ±7 K after 25 sweeps) — this
specific point is right at the edge of physical/model validity: 50 g/s of
water cannot absorb the duty this geometry delivers even as superheated
steam, driving temperature toward CoolProp's Water-backend EOS ceiling
(~3000 K). Treat these numbers as indicative of "this is not a viable
operating point," not as a precise converged result.

See `plots/comparison_summary.png` for the visual cross-comparison (Q_tot,
peak wall temperature, CHF margin, and peak quality vs. coolant mass flow,
all four fluid/geometry combinations overlaid).

## 3. Findings

### 3.1 Q_tot: helical is gas-side-dominated, confirming the earlier analysis; shell-and-tube less so

Helical Helium and Water duty are nearly identical at every mass flow
(182.1 vs 184.0 kW at 50 g/s; 214.9 vs 209.8 kW at 150 g/s — under 2.5%
apart), reconfirming `docs/context/COOLANT_COMPARISON_HELIUM_VS_WATER.md`'s
finding that this geometry is gas-side-resistance-dominated. Shell-and-tube
shows a similar pattern but with somewhat more spread (361.7 vs 323.6 kW at
50 g/s, ~10% apart) — consistent with the caveat already flagged in that
document: shell-and-tube's resistance ratio is not guaranteed to be as
lopsided as the helical coil's, and this data confirms the coolant side
carries a bit more relative weight there.

Shell-and-tube delivers roughly **2× the duty** of the 6 m helical coil at
matched coolant flow and hot-gas conditions (e.g. 453.0 vs 204.9 kW Helium at
100 g/s) — the compact 235-tube bundle simply presents much more heat
transfer area/effectiveness for this hot-gas duty than a single 6 m coil.

### 3.2 Pressure drop: water's is over an order of magnitude smaller (helical)

Helium coolant pressure drop grows from 1.55 to 6.97 bar as mass flow rises
50→150 g/s. Water's pressure drop is essentially flat and tiny (0.11-0.18
bar) across the same range — water's much higher density keeps velocity (and
therefore friction pressure drop) low even at these mass flows. (Shell-side
pressure is not yet tracked as a profile — see `docs/validation/LIQUID_HEATED_CHANNEL_SOLVER_STATUS.md`
"Remaining Limits"; a fixed nominal `p_in` is used, matching the pre-existing
gas-mode simplification in `_shell_h_at`.)

### 3.3 CHF margin: the central safety finding of this study

This is the most consequential result:

- **Helical water**: unsafe (margin 0.0022) at 50 g/s, safe and improving at
  100/150 g/s (2.50, 5.25). Higher coolant flow straightforwardly helps, as
  expected.
- **Shell-and-tube water**: CHF-exceeded (margin < 1) at **every tested mass
  flow**, and — counter-intuitively — the margin gets *worse* as coolant flow
  increases (0.275 → 0.054 → 0.029 from 50→150 g/s). This is because
  higher coolant flow in this geometry also drives up total Q substantially
  (323.6 → 523.2 kW), and the resulting local wall heat flux rises faster
  than the CHF capacity gained from the higher mass flux. **Simply adding
  more water flow does not fix CHF risk in this shell-and-tube design** — the
  geometry itself (235 tubes, 235 mm length, this baffle/corrugation
  configuration) is fundamentally CHF-limited for water across this whole
  flow range at this hot-gas duty. Fixing this would require either a
  different tube/baffle geometry, a lower hot-gas heat flux, or accepting
  operation inside the post-CHF regime this session's new closure now models
  (with the caveat that shell-side cross-flow boiling enhancement is not yet
  captured — see §0).

### 3.4 Wall temperature: shell-and-tube runs much hotter at low flow

Shell-and-tube peak wall temperature at 50 g/s is dramatically higher than
helical's (Helium: 1719 K vs 1089 K; Water: 2035 K vs 1208 K, though the
water number is the non-converged edge case). This tracks directly with
shell-and-tube's much higher heat flux density per unit coolant flow — the
same reason it is more CHF-limited.

## 4. Secondary study — feasibility of 30 L/s steam output

**Target interpretation**: 30 L/s (0.030 m³/s) of steam exiting the HX,
at minimum saturated (quality ≥ 1) at the coolant exit pressure (~80 bar,
since neither geometry currently tracks a large shell-side pressure drop).

### 4.1 Energy balance requirement (independent of HX design)

At 80 bar: `T_sat = 295.0°C`, `ρ_v,sat = 42.51 kg/m³`, `h_fg = 1441.4 kJ/kg`.
For 30 L/s of just-saturated steam:

```
ṁ_water required = 0.030 m³/s × 42.51 kg/m³ = 1.275 kg/s (1275 g/s)
Q required = ṁ_water × (h_v,sat − h_in) = 1.275 × 2625.7 kJ/kg ≈ 3.35 MW
```

This alone — before even considering whether the HX geometry can deliver that
duty — is **8-11× the water flow rates studied in Part 1** (50-150 g/s) and
requires roughly **6.4× the highest duty measured anywhere in this study**
(523.2 kW, shell-and-tube water at 150 g/s).

### 4.2 Can hot-gas flow rate be scaled up to supply 3.35 MW?

Calibration sweep (shell-and-tube, Helium coolant at 100 g/s, hot-gas mass
flow swept 0.10 → 1.0 kg/s — 10× the baseline design point):

| ṁ_hotgas [kg/s] | Q_tot [kW] |
|---|---|
| 0.10 (baseline) | 453.0 |
| 0.30 | 980.1 |
| 0.60 | 1217.5 |
| 1.00 | 1305.0 |

**Q clearly saturates** — this is classic NTU-effectiveness heat-exchanger
behavior: as the hot-gas capacity rate grows much faster than the coolant
side's, effectiveness drops and Q approaches a UA-limited ceiling rather than
scaling with hot-gas flow. Even at 10× the baseline hot-gas flow, Q reaches
only **1.31 MW — about 39% of the 3.35 MW required**, and the trend shows
clearly diminishing returns (increments of 527, 237, 88 kW per successive
mass-flow doubling-ish step) — pushing further would yield rapidly
diminishing additional duty.

### 4.3 Best achievable steam output

At the Q ceiling reached in this sweep (1.305 MW), the maximum water flow
that can be *just* fully vaporized to saturated steam is:

```
ṁ_water,max = Q_max / (h_v,sat − h_in) = 1305 kW / 2625.7 kJ/kg ≈ 0.497 kg/s (497 g/s)
Steam volumetric flow = 0.497 / 42.51 ≈ 0.0117 m³/s = 11.7 L/s
```

**≈ 39% of the 30 L/s target** — and this already assumes scaling the hot-gas
mass flow rate 10× beyond the current design point, which in a real rocket
combustor is not a "tuning knob": it corresponds to a roughly 10× larger
thrust class, i.e., a fundamentally different combustor, not an adjustment
within the current architecture.

### 4.4 Feasibility verdict

**Not feasible** with either HX architecture studied here, at any hot-gas
flow rate short of redesigning the combustor itself for ~10× the thrust
class. The gap is not a matter of fine-tuning water/hot-gas flow rates
within the current design point — it is a duty-scale gap of roughly 6-9×
between what this combustor's hot gas can deliver (even generously scaled)
and what 30 L/s of steam requires. Closing it would require either a
substantially larger combustor (more hot-gas mass flow at the same
temperature/OF ratio) or a fundamentally different, much larger heat
exchanger (more tubes/area, a longer coil, or multiple parallel HX units) —
outside the scope of "resize within the current architecture."

If a smaller steam output target is acceptable, **≈ 11-12 L/s of saturated
steam** is the realistic ceiling for the shell-and-tube geometry at 10× its
baseline hot-gas flow rate — itself untested at that flow rate beyond the
single-fluid Helium calibration sweep above (a dedicated Water run at that
duty is recommended before treating this as validated, given the CHF margin
concerns in §3.3 would need re-checking at that much higher heat flux).

## 5. Recommendations for follow-up

1. **CHF is the binding constraint for shell-and-tube water**, not duty or
   pressure drop. Any real design targeting significant steam generation in
   this geometry needs either a shell-side boiling-enhancement correlation
   (closing the gap flagged in §0) or a geometry change, before trusting the
   post-CHF numbers as a design basis.
2. Shell-side coolant pressure drop is not yet modeled (fixed at nominal
   `p_in`) — add this if pressure-drop-driven saturation-temperature shifts
   become design-relevant (they would matter more as duty scales up toward
   the numbers in §4).
3. The pre-existing Helium baseline regression failure (§0.1) should be
   investigated and re-pinned separately — it predates this study and blocks
   trusting the regression suite as a drift detector until fixed.
