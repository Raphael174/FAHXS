# Model Calibration Methodology — Combustor-HX 1D Solver

> Code: [`optimization/calibrate.py`](../optimization/calibrate.py)  
> Parameters: [`1Dmodel/input_data.py → CorrelationCoefficients`](../1Dmodel/input_data.py)

---

## 1. What calibration is — and is not

The 1D solver uses empirical correlations as **parametric surrogate models**. Each correlation has a published functional form (e.g. Nu = a·Re^b·Pr^1/3·geometry^c) and coefficients fitted to the original author's data. When applied to this specific geometry — square chamber cross-section, diesel/O2 combustion gas, He at 50 bar — those coefficients are extrapolated well beyond their validation range.

Calibration adjusts the **prefactors** of the active correlations so that model outputs match bench measurements. This is not a validation of the correlation's physical form — it is a geometry- and fluid-specific recalibration. The calibrated coefficients are properties of *this rig*, not of the correlation in general. They do not generalize to other geometries, but they do not need to.

The key methodological point (Kennedy & O'Hagan 2001): even after calibration, there is irreducible **structural model error** δ(x) from using the wrong functional form for a square chamber, omitting secondary flow effects, etc. This δ is absorbed partially into the calibrated prefactor, but residuals of ±5–15% on heat flux should be expected.

---

## 2. Active parameters at this bench — reduced working set

The full `CorrelationCoefficients` has 12 fields. For this bench and operating range, most are either inactive, insensitive, or should be fixed on physical grounds.

### Fixed — do not calibrate

| Parameter | Fixed value | Reason |
|---|---|---|
| `emissivity_wall` | **0.95** | Coil surface is sooted/charred in current tests → near-black. Range: 0.90–0.97 for carbon soot and oxidised steel in combustion environments (Touloukian & DeWitt 1970). New/clean 316L after installation: ε ≈ 0.85. Update this value at the start of each test campaign depending on coil condition. |
| `mbl_factor` | **3.4** | Geometric constant for the square chamber. Sensitivity to ±15% variation is secondary to friction and Nu. Fix initially; revisit only if radiation residuals are large. |
| `mori_a_lo`, `mori_b_lo`, `mori_c_lo` | literature | Coil-side Nu accepted at ±20% — see below. |
| `salimpour_b`, `salimpour_c`, `kays_crawford_n` | literature | Functional-form parameters; calibrate only the prefactor (`salimpour_a`), not the exponents. Fitting exponents requires a Re sweep which is not available here. |
| `ali_c_lo`, `ali_I_split` | literature | High-I branch is always active at this bench (I ≈ 9–11). Low-I branch is never reached. |

### Helium-side Nu — accept uncertainty, do not calibrate

At 50 bar / 303 K with m_dot_c = 80–100 g/s and Dh = 7–10 mm, Re_He ≈ 700 000–900 000. Curvature parameter α = Dh/(2Rc) ≈ 0.11. The coil-side resistance is almost certainly **not the dominant thermal resistance** at these operating conditions — the shell-side (low-Re combustion gas) and wall conduction dominate. A ±20% error on Nu_He has minor effect on total UA.

Use one fixed correlation from the selector:

```python
combustorProp.Nusselt_coil = "Gnielinski"    # simplest, well-validated for turbulent pipe
# or
combustorProp.Nusselt_coil = "mori1967"      # adds curvature correction (recommended)
```

The difference between Gnielinski and Mori-1967 at these Re/Pr is small (< 10%). Pick one and keep it fixed for the campaign so results are comparable across runs.

### Calibration targets — working set

| Parameter | Default | Description | Identifiable from |
|---|---|---|---|
| **`ali_c_hi`** | 0.325 | High-curvature friction prefactor | `dp_c` — always |
| **`salimpour_a`** | 0.317 | Shell-side Nu prefactor | `T_g_out` — if measured |
| *(optional)* `mori_a_lo` | 26.2 | Coil-side Nu denominator | `T_g_out` — if measured, and partially correlated with salimpour_a |

**With nominal observables only (no T_g_out): calibrate `ali_c_hi` only.**

---

## 3. What your observables actually tell you

### Nominal observable set (always available)

| Measured quantity | How obtained | Informs |
|---|---|---|
| `dp_c = p_He_in − p_He_out` | Differential pressure transducer | Friction → `ali_c_hi` |
| `Q_He = m_dot_c·(h_He_out − h_He_in)` | CoolProp from T_He_in, T_He_out, p_He_in, p_He_out, m_dot_c | Energy balance check |
| `T_g_in` (inferred) | Cantera equilibrium at known OF, p_chamber, combustion efficiency | Gas energy budget |

### The Q_He identifiability problem in counter-flow

For a counter-flow HX with C_He >> C_g (which is the likely regime here):

```
C_He  = m_dot_c × cp_He ≈ 0.090 × 5193 ≈ 467 W/K
C_g   = m_dot_g × cp_g  ≈ 0.075 × 1400 ≈ 105 W/K   [diesel/O2 at ~1500K]
C_r   = C_g / C_He      ≈ 0.22
```

Counter-flow effectiveness: ε = f(NTU, C_r). With 27 turns and reasonable UA, NTU ≫ 1 and ε → 1/(1−C_r)·... → near unity. In this regime:

```
Q_He ≈ C_g × (T_g_in − T_He_in)
```

The helium outlet temperature is determined almost entirely by the gas energy content, not by Nu. **Q_He is insensitive to Nu in this operating regime.** This is a fundamental consequence of having many coil turns and C_He >> C_g — not a modeling bug.

If you want to see Nu sensitivity in Q_He, the NTU must be low enough that ε < 0.85, which would require either fewer turns or much lower UA. Neither is available here.

### What T_g_out would give you (if measured)

In the same regime where Q_He saturates, T_g_out continues to vary with Nu:

```
T_g_out = T_He_in + ε_min × (T_g_in − T_He_in)
```

where ε_min refers to the minimum-capacitance-stream (gas) outlet. This is sensitive to Nu because it tracks how thoroughly the gas is cooled. A 25% change in Nu_shell (salimpour_a) can move T_g_out by 50–100 K even when Q_He barely changes. This is why T_g_out is the right observable for Nu calibration.

**Without T_g_out, the campaign effectively calibrates one parameter: `ali_c_hi`.**  
That is still valuable — it reduces pressure-drop prediction uncertainty from ±15% to ±3–5%.

### Inferring T_g_in from operating conditions

The model's Cantera equilibrium gives T_flame(OF, p) at 100% combustion efficiency. Real combustion efficiency in a liquid-propellant burner is η_c\* ≈ 90–98% (based on c\* efficiency). The effective gas inlet temperature to the HX section is lower:

```
T_g_in_eff ≈ T_He_in + η_c × (T_flame_adiabatic − T_He_in)
```

For now, treat η_c as a known constant (e.g. 0.95) and update it if the energy balance `Q_He vs C_g·(T_g_in_eff − T_g_out_predicted)` shows a systematic bias. With only the nominal observable set, η_c and salimpour_a are co-mingled — you cannot separate them without T_g_out.

---

## 4. Methodology summary

### 4.1 Log-posterior and uncertainty budget

```
log p(θ | data) = log L(θ | data) + log p(θ)

log L = −½ Σ_i  [ (y_model,i(θ) − y_exp,i)² / σ_total,i² ]

σ_total² = σ_meas²  +  σ_model²          [Kennedy & O'Hagan 2001]
```

`σ_model` represents the structural mismatch (square chamber, Pr extrapolation). Default values in `CalibrationRecord`:

| Observable | σ_meas | σ_model | σ_total (approx) |
|---|---|---|---|
| dp_c | 3% (transducer) | 10% (friction correlation) | 10.4% |
| T_g_out | 50 K + probe loss | 8% of T_g_out | ~70–100 K |
| Q_He | 5% (flow + T sensors) | 10% (energy balance) | 11% |

### 4.2 Prior distributions

Prefactors (ali_c_hi, salimpour_a) use **log-normal** priors — strictly positive, multiplicative uncertainty. A log-normal with σ = 0.12 (CV = 12%) means a 68% prior credible interval of [0.325·e^(−0.12), 0.325·e^(0.12)] = [0.288, 0.367] for ali_c_hi. This prevents the optimizer from finding unphysical values.

The prior prevents overfitting with sparse data: if only one run is available and the observables are weakly informative, the MAP estimate stays near the literature value.

### 4.3 Three inference levels

| Step | Function | When | Output |
|---|---|---|---|
| 1 | `calibrate_ls()` | ≥ 1 run | Point estimate, ~1 min |
| 2 | `calibrate_map()` | ≥ 1 run | Prior-regularised point estimate, ~2 min |
| 3 | `calibrate_mcmc()` | ≥ 3 runs | Full posterior, ~30 min–2 h |

Always start with Step 1. If LS and MAP give the same answer, the data is strongly constraining. If they differ significantly, the prior is controlling the result — you need more data or a different observable.

---

## 5. Experimental design

### Controllable variables

| Variable | Range | Effect |
|---|---|---|
| m_dot_g | 50–100 g/s | T_g_in, Re_g → Nu_shell, radiation |
| OF | 2.0–2.5 | T_flame, CO2/H2O mole fractions → radiation |
| m_dot_c | 80–100 g/s | Re_He → dp_c strongly |
| p_He | variable | ρ_He → dp_c inversely |

### What is identifiable vs not from nominal observables

| Parameter | Lever | Signal in dp_c | Signal in Q_He | Signal in T_g_out |
|---|---|---|---|---|
| `ali_c_hi` | m_dot_c | **strong** (dp ∝ m_dot^1.8) | none | none |
| `salimpour_a` | m_dot_g | none | negligible | **strong** |
| `mori_a_lo` | m_dot_c | none | negligible | moderate |
| `emissivity_wall × mbl` | OF at fixed m_dot_g | none | negligible | moderate |

### Recommended test matrix

Minimum 4 runs to constrain `ali_c_hi` well and verify the energy balance.
Add 3 more if T_g_out becomes available.

**Core campaign — friction calibration (nominal observables only):**

| Test | m_dot_g (g/s) | OF | m_dot_c (g/s) | Purpose |
|---|---|---|---|---|
| **T01a** | 75 | 2.25 | 90 | Baseline |
| **T01b** | 75 | 2.25 | 90 | Repeat of T01 — **this gives σ_meas experimentally** |
| **T02** | 75 | 2.25 | 80 | Low dp → ali_c_hi lower bound |
| **T03** | 75 | 2.25 | 100 | High dp → ali_c_hi upper bound |
| **T04** | 50 | 2.25 | 90 | Low gas energy → check energy balance at low Q |
| **T05** | 100 | 2.25 | 90 | High gas energy → check saturation of Q_He |

> T01a and T01b are the most important runs. Repeatability gives your actual σ_meas, which calibrates the likelihood function. Without it, σ_meas is an assumption.

**Extension — Nu calibration (requires T_g_out):**

| Test | m_dot_g (g/s) | OF | m_dot_c (g/s) | Purpose |
|---|---|---|---|---|
| **T06** | 50 | 2.25 | 90 | Low Re_g → Nu_shell low → T_g_out high |
| **T07** | 100 | 2.25 | 90 | High Re_g → Nu_shell high → T_g_out low |
| **T08** | 75 | 2.0 | 90 | Richer mix, higher flame T → radiation |

### What m_dot_c variations buy you

In the high-I Ali branch (always active here):

```text
dp_c ∝ ali_c_hi × m_dot_c^1.8 × p_He^(−1)   [approximately]
```

80→100 g/s m_dot_c change: dp_c changes by ~(100/80)^1.8 ≈ **+47%**. That is 15× the transducer uncertainty. Even with only T02 and T03, `ali_c_hi` is very tightly constrained.

Holding m_dot_c fixed and varying p_He also changes dp_c — if you have pressure controllability, a 40→55 bar sweep at fixed m_dot_c would halve then restore dp_c and provide an independent check.

---

## 6. Workflow

### Before testing

```python
# 1. Update input_data.py to bench geometry
#    combustorProp.inner_diameter, Dh_coil, coil_gap, thickness_coil_wall,
#    gap_shell2coil, N_coils (pass count, not turns)
#    coolantProp.T_in=303.0, p_in=50e5, p_out=...
#    combustorProp.Nusselt_coil = "Gnielinski"  (or "mori1967", pick and fix)
#    CorrelationCoefficients.emissivity_wall = 0.95  (sooted coil)

# 2. Check sensitivity at baseline conditions
from optimization.calibrate import compute_sensitivities, identifiability_check, CalibrationRecord

rec_T01 = CalibrationRecord(
    m_dot_c=0.090, T_He_in=303.0, T_He_out=...,    # fill in nominal estimates
    p_He_in=50e5, p_He_out=..., m_dot_g=0.075, OF=2.25
)
S = compute_sensitivities(rec_T01, params=["ali_c_hi"])
print(S)   # expect S[dp, ali_c_hi] ≈ 0.6 (strong signal)
```

### After the first two runs (T01a, T01b)

```python
# Compute experimental repeatability
dp_a, dp_b = ..., ...    # measured values
sigma_dp_meas_actual = abs(dp_a - dp_b) / dp_a   # relative spread
# Update CalibrationRecord defaults with this σ before all subsequent calibrations
```

### After T01–T03 (friction calibration)

```python
from optimization.calibrate import calibrate_ls, calibrate_map

records = [rec_T01a, rec_T02, rec_T03]
ls  = calibrate_ls(records, params=["ali_c_hi"])
bmap = calibrate_map(records, params=["ali_c_hi"])

# If ls["x"] ≈ bmap["x"]: data is strongly constraining, take the LS value.
# If they differ by > 20%: prior is dominant, need more data or wider m_dot_c range.
print(f"ali_c_hi: LS = {ls['x'][0]:.4f}, MAP = {bmap['x'][0]:.4f}, literature = 0.325")
```

### After T01–T05 (energy balance check)

Run the calibrated model at T04 and T05 conditions and compare predicted Q_He to measured. If the model reproduces Q_He within σ_total (≈ 11%) without any Nu calibration, the energy budget is correctly closed and you have no strong evidence for Nu error. If there is a systematic bias (always over- or under-predicts Q_He), the combustion efficiency η_c assumed in T_g_in may be wrong — adjust it before attributing the residual to Nu.

### Full posterior (after 4+ runs, optional)

```python
from optimization.calibrate import calibrate_mcmc, plot_posteriors

mcmc = calibrate_mcmc(records, params=["ali_c_hi"],
                       nwalkers=20, nsteps=2000, burnin=400)
print(f"ali_c_hi: {mcmc['mean'][0]:.4f} ± {mcmc['std'][0]:.4f}")
print(f"Acceptance: {mcmc['acceptance']:.3f}  (want 0.20–0.50)")
print(f"Converged: {mcmc['converged']},  ESS: {mcmc['ess'][0]:.0f}")
plot_posteriors(mcmc)
```

---

## 7. Expected outcomes

With the nominal observable set and 4–6 runs varying m_dot_c:

| Quantity | Pre-calibration | Post-calibration |
|---|---|---|
| dp_c prediction | ±15% (literature ali_c_hi extrapolated) | ±3–5% |
| Q_He prediction | ±15% (energy balance dominated) | ±8–12% (residual: η_c uncertainty) |
| T_g_out (if measured later) | ±20% (salimpour_a not calibrated) | ±15% (reduced by knowing ali_c_hi) |

After friction calibration, the dominant residual on Q_He is **combustion efficiency η_c** (enters through T_g_in) and **shell-side Nu** (enters through UA). These require T_g_out to separate.

---

## 8. Instrumentation priorities

| Instrument | Required accuracy | Impact if missing |
|---|---|---|
| Differential pressure transducer (He) | ±0.5% FS | Calibration impossible — primary target |
| Pt100 RTD — T_He_in, T_He_out | ±0.5–1 K | Needed for Q_He (secondary observable) |
| Coriolis flowmeter — m_dot_c | ±0.5% | Needed for Q_He and Re_He; preferred over thermal for He |
| Shielded thermocouple — T_g_out | ±30 K (aspirated probe) | Enables Nu calibration; bare probe unusable (±80 K radiation error) |

---

## 9. References

1. Kennedy, M.C. & O'Hagan, A. (2001). Bayesian calibration of computer models. *J. R. Statist. Soc. B* 63(3), 425–464.
2. Salimpour, M.R. (2008). Heat transfer coefficients of shell-and-coiled tube heat exchangers. *Int. J. Therm. Sci.* 47, 1027–1033.
3. Ali, M. et al. (2024). Experimental study of pressure drop in helically coiled tubes. *Exp. Therm. Fluid Sci.* 154, 111126.
4. Mori, Y. & Nakayama, W. (1967). Study on forced convective heat transfer in curved pipes. *Int. J. Heat Mass Transfer* 10(5), 681–695.
5. Touloukian, Y.S. & DeWitt, D.P. (1970). *Thermophysical Properties of Matter, Vol. 7*. Plenum. *(emissivity of sooted surfaces)*
6. Hottel, H.C. & Sarofim, A.F. (1967). *Radiative Heat Transfer*. McGraw-Hill. *(mean beam length)*
7. Foreman-Mackey, D. et al. (2013). emcee: the MCMC hammer. *PASP* 125, 306–312.
