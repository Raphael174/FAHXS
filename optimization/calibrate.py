"""
Parameter calibration for the Combustor-HX_V2 1D solver.

THEORETICAL FRAMEWORK — Kennedy & O'Hagan (2001) Bayesian calibration
======================================================================
The model output y_model(θ) is connected to the true observable y_exp by:

    y_exp = y_model(θ_true) + ε_model + ε_meas

where θ is the calibration parameter vector, ε_meas is measurement noise,
and ε_model is model-form error (structural discrepancy).  The combined
uncertainty is modelled as:

    σ_total² = σ_meas² + σ_model²

with σ_meas from instrument specs and σ_model estimated as a fixed fraction
of the nominal model output (typically 5–15 % for engineering correlations).
The log-likelihood becomes:

    log L(θ) = -½ Σ_i [(y_model,i(θ) - y_exp,i)² / σ_total,i²]
                + const

OBSERVABLES AND IDENTIFIABILITY
================================
dp_c   [Pa]   helium pressure drop = p_He_in - p_He_out
T_g_out [K]   hot gas outlet temperature (thermocouple / IR — "bonus" measurement)
Q_He   [W]    total heat to helium = m_dot_c*(h_He_out - h_He_in)

Sensitivity analysis at nominal design point (counter-flow, He 120→750 K, dp=17 bar):

    dp_c      → ali_c_hi (friction):  +54 % coeff → +34 % dp.  Cleanly identifiable.
    T_g_out   → salimpour_a, mori_a_lo:  gas exit T is sensitive to Nu (how fast heat
                is extracted per turn).  Identifiable if T_g_out is measured.
    Q_He      → NOT sensitive to Nu in counter-flow at the nominal operating point.
                Reason: the HX is gas-energy-limited (T_c_final ≈ 13 K < T_He_in).
                Q_He is still useful as a consistency check.
    NOTE: in co-flow, Q_He is Nu-sensitive because the HX never saturates.

Practical identifiability (sparse data):
    With dp_c only            → ali_c_hi                           1 param
    With dp_c + T_g_out       → ali_c_hi + salimpour_a (or mori)   2 params
    With dp_c + T_g_out + Q_He (varying m_dot_c or geometry)       3+ params

Use compute_sensitivities() + identifiability_check() before attempting multi-
parameter calibration with a single experimental record.

PRIOR DISTRIBUTIONS
===================
The PRIOR_SPEC dict below encodes physics-informed priors.  Distribution choice:

  log-normal  — for strictly-positive prefactors with multiplicative uncertainty
                (a ±CV error in the physical coefficient maps to ±CV in log-space).
                Parameterised by (mu, sigma) of the underlying normal, where
                mu = ln(x_lit) and sigma = CV (e.g. 0.25 → 25 % coefficient of variation).

  normal      — for signed exponents where CV has an additive meaning.
                Parameterised by (mu, sigma) in parameter space.

  truncated   — for bounded parameters (emissivity ∈ (0, 1]).

Each entry: {dist, mu, sigma, bounds, source}

Literature CV references
------------------------
salimpour_a  : CV ≈ 25 %  — Salimpour (2008) reports ±17 % fit on held-out data,
               but his Pr range (4–15) does not cover Pr ~ 0.65 for combustion gas
               (extrapolation ×2 in Pr).  Conservative: 25 %.
               Source: Salimpour, M.R. (2008) Int. J. Therm. Sci. 47, 1027–1033.

salimpour_b  : CV ≈ 5 %   — Re exponent is well-constrained by curve-fit.
salimpour_c  : σ ≈ 0.04   — small exponent on pitch/D_o ratio; normal prior.
kays_crawford_n : σ ≈ 0.05 — Kays & Crawford (1993) T-ratio correction exponent.
               Typical published range 0.11 (Sieder-Tate gases) – 0.36 (liquids).

mori_a_lo    : CV ≈ 15 %  — Mori & Nakayama (1967) Pr-range is Pr < 1 which
               covers He (Pr ≈ 0.67).  Authors report ±10 % on original data;
               15 % accounts for geometry extrapolation.

ali_c_hi     : CV ≈ 12 %  — Ali et al. (2024) fit to their own helical tube data.
               High-I branch is active at design (I ≈ 3.3); c_hi directly controls dp.
               Source: Ali, M. et al. (2024) Exp. Therm. Fluid Sci. 154, 111126.

ali_c_lo     : CV ≈ 15 %  — low-I branch, not active at design point.
ali_I_split  : σ ≈ 0.05   — branch threshold; treated as structural parameter.

mbl_factor   : CV ≈ 15 %  — Hottel (1954) recommends 3.6 V/A for arbitrary geometry,
               3.4 with non-grey correction.  Conservative given geometric extrapolation.

emissivity_wall : σ ≈ 0.10 — Touloukian & DeWitt (1970): oxidised 316L ε = 0.70–0.88;
               burnt/carbon-coated can reach 0.92.  Truncated normal on (0, 1].

CALIBRATION WORKFLOW
====================
1. Run a forward solve with nominal parameters → identify operating regime.
2. Call compute_sensitivities() → pick the 1–2 most sensitive parameters.
3. Call identifiability_check() → verify Fisher information is well-conditioned.
4. Call calibrate_ls() for a fast point estimate.
5. Call calibrate_map() for a regularised Bayesian point estimate.
6. Call calibrate_mcmc() for full posteriors + uncertainty quantification.
7. Call posterior_predictive() to propagate parameter uncertainty to model outputs.
8. With multiple experimental records: call calibrate_sequential() for sequential
   Bayesian updating (each experiment refines the posterior).

Usage example
-------------
    from optimization.calibrate import (
        CalibrationRecord, calibrate_ls, calibrate_map, calibrate_mcmc,
        compute_sensitivities, identifiability_check, plot_posteriors,
    )
    import numpy as np

    rec = CalibrationRecord(
        m_dot_c=0.15, T_He_in=120.0, T_He_out=750.0,
        p_He_in=90e5, p_He_out=73e5,
        m_dot_g=0.070, OF=2.9,
        T_g_out=950.0,          # measured gas outlet temperature
    )

    # 1. Check sensitivity first
    S = compute_sensitivities(rec, params=["salimpour_a", "ali_c_hi"])
    print("Sensitivity matrix:\\n", S)

    # 2. Identifiability
    cond = identifiability_check(rec, params=["salimpour_a", "ali_c_hi"])
    print(f"FIM condition number: {cond:.1f}  (< 100 → well-conditioned)")

    # 3. Point estimate
    ls = calibrate_ls(rec, params=["salimpour_a", "ali_c_hi"])

    # 4. Bayesian MAP
    bmap = calibrate_map(rec, params=["salimpour_a", "ali_c_hi"])

    # 5. Full posteriors (slow — requires: pip install emcee)
    mcmc = calibrate_mcmc(rec, params=["salimpour_a", "ali_c_hi"],
                           nwalkers=24, nsteps=3000, burnin=600)
    plot_posteriors(mcmc)

References
----------
Kennedy, M.C. & O'Hagan, A. (2001) Bayesian calibration of computer models.
    J. R. Statist. Soc. B 63(3), 425–464.
Salimpour, M.R. (2008) Int. J. Therm. Sci. 47, 1027–1033.
Mori, Y. & Nakayama, W. (1967) Int. J. Heat Mass Transfer 10(5), 681–695.
Ali, M. et al. (2024) Exp. Therm. Fluid Sci. 154, 111126.
Touloukian, Y.S. & DeWitt, D.P. (1970) Thermophysical Properties of Matter, Vol. 7.
Hottel, H.C. & Sarofim, A.F. (1967) Radiative Heat Transfer. McGraw-Hill.
Kays, W.M. & Crawford, M.E. (1993) Convective Heat and Mass Transfer. McGraw-Hill.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields, replace
from typing import Sequence

import numpy as np
from scipy.optimize import least_squares, minimize
from CoolProp.CoolProp import PropsSI

from hps_combustor.input_data import (
    coolantProp, hotgasProp, combustorProp, numericalProp,
    system_requirements, CorrelationCoefficients,
)
from hps_combustor.main_solve import main_solver


# ---------------------------------------------------------------------------
# Prior specification — one entry per calibration-eligible parameter
# ---------------------------------------------------------------------------

#  Keys match CorrelationCoefficients field names.
#  dist   : "lognormal" | "normal" | "truncated_normal"
#  mu     : distribution mean (parameter space for normal; ln(x_lit) for lognormal)
#  sigma  : distribution std  (CV for lognormal, absolute σ for normal)
#  bounds : hard physical bounds (lower, upper) — None means ±∞
#  source : short citation for the prior width

_lit = {f.name: f.default for f in fields(CorrelationCoefficients)}

PRIOR_SPEC: dict[str, dict] = {
    "salimpour_a": dict(
        dist="lognormal",
        mu=math.log(_lit["salimpour_a"]),
        sigma=0.25,
        bounds=(0.05, 2.0),
        source="Salimpour (2008) ±17% in-sample; +extrapolation to Pr~0.65 → CV=25%",
    ),
    "salimpour_b": dict(
        dist="lognormal",
        mu=math.log(_lit["salimpour_b"]),
        sigma=0.05,
        bounds=(0.4, 0.9),
        source="Salimpour (2008) well-constrained Re exponent; CV=5%",
    ),
    "salimpour_c": dict(
        dist="normal",
        mu=_lit["salimpour_c"],
        sigma=0.04,
        bounds=(-0.5, 0.0),
        source="Salimpour (2008) pitch/D_o exponent; σ=0.04 (small negative exponent)",
    ),
    "kays_crawford_n": dict(
        dist="normal",
        mu=_lit["kays_crawford_n"],
        sigma=0.05,
        bounds=(0.11, 0.40),
        source="Kays & Crawford (1993) T-ratio correction 0.11–0.36; σ=0.05",
    ),
    "mori_a_lo": dict(
        dist="lognormal",
        mu=math.log(_lit["mori_a_lo"]),
        sigma=0.15,
        bounds=(10.0, 60.0),
        source="Mori & Nakayama (1967) ±10% fit + geometry extrapolation → CV=15%",
    ),
    "mori_b_lo": dict(
        dist="normal",
        mu=_lit["mori_b_lo"],
        sigma=0.015,
        bounds=(0.0, 0.20),
        source="Mori & Nakayama (1967); small correction term σ=0.015",
    ),
    "mori_c_lo": dict(
        dist="normal",
        mu=_lit["mori_c_lo"],
        sigma=0.020,
        bounds=(0.0, 0.30),
        source="Mori & Nakayama (1967); Dean-number correction σ=0.020",
    ),
    "ali_c_lo": dict(
        dist="lognormal",
        mu=math.log(_lit["ali_c_lo"]),
        sigma=0.15,
        bounds=(0.1, 1.0),
        source="Ali et al. (2024) low-I branch; CV=15% (inactive at design)",
    ),
    "ali_c_hi": dict(
        dist="lognormal",
        mu=math.log(_lit["ali_c_hi"]),
        sigma=0.12,
        bounds=(0.1, 1.0),
        source="Ali et al. (2024) high-I branch (active at design I≈3.3); CV=12%",
    ),
    "ali_I_split": dict(
        dist="normal",
        mu=_lit["ali_I_split"],
        sigma=0.05,
        bounds=(0.5, 1.5),
        source="Ali et al. (2024) branch threshold; σ=0.05 (structural param)",
    ),
    "mbl_factor": dict(
        dist="lognormal",
        mu=math.log(_lit["mbl_factor"]),
        sigma=0.15,
        bounds=(2.0, 5.0),
        source="Hottel (1954) 3.4–3.6 for non-grey correction; CV=15%",
    ),
    "emissivity_wall": dict(
        dist="truncated_normal",
        mu=_lit["emissivity_wall"],
        sigma=0.10,
        bounds=(0.3, 1.0),
        source="Touloukian & DeWitt (1970) oxidised 316L 0.70–0.88; σ=0.10",
    ),
}


def _log_prior(params: Sequence[str], x: np.ndarray) -> float:
    """Log-prior log p(θ) summed over all calibration parameters."""
    lp = 0.0
    for name, val in zip(params, x):
        spec = PRIOR_SPEC[name]
        lo, hi = spec["bounds"]
        if not (lo <= val <= hi):
            return -np.inf
        dist = spec["dist"]
        if dist == "lognormal":
            if val <= 0:
                return -np.inf
            lp += -0.5 * ((math.log(val) - spec["mu"]) / spec["sigma"]) ** 2
            lp -= math.log(val)          # Jacobian of log transform (log-normal normalisation)
        elif dist in ("normal", "truncated_normal"):
            lp += -0.5 * ((val - spec["mu"]) / spec["sigma"]) ** 2
    return lp


# ---------------------------------------------------------------------------
# Experimental record
# ---------------------------------------------------------------------------

@dataclass
class CalibrationRecord:
    """One experimental operating point.

    Minimum required observable:  dp_c (pressure drop) → calibrates ali_c_hi.
    Strongly recommended:         T_g_out → calibrates salimpour_a / Nu.
    Optional:                     Q_He_W — computed from CoolProp if None.

    Uncertainty fields (σ) represent total measurement uncertainty 1-sigma.
    sigma_model_* are structural model discrepancy fractions (ε_model / y_nominal),
    combined internally: σ_total² = σ_meas² + σ_model².
    """
    m_dot_c: float          # He mass flow [kg/s]
    T_He_in: float          # He inlet temperature [K]
    T_He_out: float         # He outlet temperature [K]
    p_He_in: float          # He inlet pressure [Pa]
    p_He_out: float         # He outlet pressure [Pa]
    m_dot_g: float          # hot-gas mass flow [kg/s]
    OF: float               # oxidiser-to-fuel mass ratio
    T_g_out: float = None   # hot gas outlet temperature [K]; None if not measured
    Q_He_W: float = None    # total heat [W]; computed from CoolProp if None

    # Measurement uncertainty (1-sigma, instrument specs)
    sigma_dp_meas_rel: float = 0.03   # relative uncertainty on dp_c (3 % — differential pressure transducer)
    sigma_Tg_meas_K: float = 50.0     # absolute uncertainty on T_g_out [K] (type-K thermocouple + probe loss)
    sigma_Q_meas_rel: float = 0.05    # relative uncertainty on Q_He (5 % — mass-flow meter + T sensors)

    # Model structural discrepancy (fraction of nominal — Kennedy-O'Hagan ε_model)
    sigma_dp_model_rel: float = 0.10  # dp correlation structural error ≈ 10 %
    sigma_Tg_model_rel: float = 0.08  # T_g_out structural error ≈ 8 %
    sigma_Q_model_rel: float = 0.10   # Q_He structural error ≈ 10 %

    def __post_init__(self):
        if self.Q_He_W is None:
            h_in  = PropsSI('H', 'T', self.T_He_in,  'P', self.p_He_in,  'Helium')
            h_out = PropsSI('H', 'T', self.T_He_out, 'P', self.p_He_out, 'Helium')
            self.Q_He_W = self.m_dot_c * abs(h_out - h_in)

    @property
    def dp_c(self) -> float:
        return self.p_He_in - self.p_He_out

    def sigma_dp(self) -> float:
        """Combined dp uncertainty (meas + model) [Pa]."""
        s_m = self.sigma_dp_meas_rel * self.dp_c
        s_d = self.sigma_dp_model_rel * self.dp_c
        return math.sqrt(s_m**2 + s_d**2)

    def sigma_Tg(self) -> float:
        """Combined T_g_out uncertainty [K]."""
        if self.T_g_out is None:
            return 1.0
        s_d = self.sigma_Tg_model_rel * self.T_g_out
        return math.sqrt(self.sigma_Tg_meas_K**2 + s_d**2)

    def sigma_Q(self) -> float:
        """Combined Q_He uncertainty [W]."""
        s_m = self.sigma_Q_meas_rel * self.Q_He_W
        s_d = self.sigma_Q_model_rel * self.Q_He_W
        return math.sqrt(s_m**2 + s_d**2)


# ---------------------------------------------------------------------------
# Solver wrapper
# ---------------------------------------------------------------------------

def _make_corr(params: Sequence[str], x: np.ndarray) -> CorrelationCoefficients:
    """Build a CorrelationCoefficients instance with the given parameter overrides."""
    overrides = dict(zip(params, x))
    return replace(CorrelationCoefficients(), **overrides)


def _run_solver(record: CalibrationRecord, corr: CorrelationCoefficients):
    """Run a full solver pass and return (Q_He_W, dp_c_Pa, T_g_out_K)."""
    cp = coolantProp(
        mass_flow_c=record.m_dot_c,
        T_in=record.T_He_in,
        T_out=record.T_He_out,
        p_in=record.p_He_in,
        p_out=record.p_He_out,
    )
    hg = hotgasProp(
        mass_flow_g=record.m_dot_g,
        mixing_ratio=record.OF,
    )
    solver = main_solver(
        coolantProp=cp,
        hotgasProp=hg,
        combustorProp=combustorProp(),
        numericalProp=numericalProp(),
        system_requirements=system_requirements(),
        corrCoeffs=corr,
    )
    solver.solver()
    solver.compute_performance()
    Q_model  = solver.Q_He * 1e3                           # kW → W
    dp_model = solver.dp_c_tot                             # Pa
    Tg_model = float(solver.data_master["T_g"][-1])        # K
    return Q_model, dp_model, Tg_model


def _log_likelihood(record: CalibrationRecord,
                    Q_m: float, dp_m: float, Tg_m: float) -> float:
    """Gaussian log-likelihood with combined measurement + model uncertainty."""
    ll = -0.5 * ((dp_m - record.dp_c) / record.sigma_dp()) ** 2
    if record.T_g_out is not None:
        ll += -0.5 * ((Tg_m - record.T_g_out) / record.sigma_Tg()) ** 2
    else:
        # Q_He as low-weight consistency check (not the calibration driver)
        ll += -0.5 * ((Q_m - record.Q_He_W) / record.sigma_Q()) ** 2
    return ll


def _log_posterior_single(params: Sequence[str], x: np.ndarray,
                           record: CalibrationRecord) -> float:
    """log p(θ | data) for one record (unnormalised)."""
    lp = _log_prior(params, x)
    if not np.isfinite(lp):
        return -np.inf
    corr = _make_corr(params, x)
    try:
        Q_m, dp_m, Tg_m = _run_solver(record, corr)
    except Exception:
        return -np.inf
    return lp + _log_likelihood(record, Q_m, dp_m, Tg_m)


def _log_posterior_multi(params: Sequence[str], x: np.ndarray,
                          records: list[CalibrationRecord]) -> float:
    """log p(θ | D) for multiple independent records."""
    lp = _log_prior(params, x)
    if not np.isfinite(lp):
        return -np.inf
    corr = _make_corr(params, x)
    ll_total = 0.0
    for rec in records:
        try:
            Q_m, dp_m, Tg_m = _run_solver(rec, corr)
        except Exception:
            return -np.inf
        ll_total += _log_likelihood(rec, Q_m, dp_m, Tg_m)
    return lp + ll_total


# ---------------------------------------------------------------------------
# Sensitivity analysis and identifiability
# ---------------------------------------------------------------------------

def compute_sensitivities(
    record: CalibrationRecord,
    params: Sequence[str] = ("salimpour_a", "ali_c_hi"),
    h_rel: float = 0.05,
) -> np.ndarray:
    """Finite-difference normalised sensitivity matrix.

    S[i, j] = (θ_j / y_i) * ∂y_i/∂θ_j  (dimensionless — logarithmic sensitivity)

    Rows: [dp_c, T_g_out (or Q_He)],  Columns: one per parameter in `params`.

    A value of 1.0 means a 1 % change in θ_j causes a 1 % change in y_i.
    """
    x0 = np.array([_lit[p] for p in params])
    corr0 = _make_corr(params, x0)
    Q0, dp0, Tg0 = _run_solver(record, corr0)
    obs0 = np.array([dp0, Tg0 if record.T_g_out is not None else Q0])

    n_obs = len(obs0)
    n_params = len(params)
    S = np.zeros((n_obs, n_params))

    for j, (name, xj) in enumerate(zip(params, x0)):
        x_plus = x0.copy()
        x_plus[j] = xj * (1 + h_rel)
        Q_p, dp_p, Tg_p = _run_solver(record, _make_corr(params, x_plus))
        obs_p = np.array([dp_p, Tg_p if record.T_g_out is not None else Q_p])

        x_minus = x0.copy()
        x_minus[j] = xj * (1 - h_rel)
        Q_m, dp_m, Tg_m = _run_solver(record, _make_corr(params, x_minus))
        obs_m = np.array([dp_m, Tg_m if record.T_g_out is not None else Q_m])

        dobs_dtheta = (obs_p - obs_m) / (2 * h_rel * xj)
        S[:, j] = (xj / obs0) * dobs_dtheta

    return S


def identifiability_check(
    record: CalibrationRecord,
    params: Sequence[str] = ("salimpour_a", "ali_c_hi"),
    h_rel: float = 0.05,
) -> dict:
    """Fisher information matrix analysis for identifiability.

    Returns a dict with:
        condition_number  — FIM condition number (< 100: well-conditioned, > 1000: poor)
        eigenvalues       — eigenvalues of FIM (small → poorly identified direction)
        sensitivity_matrix — normalised S matrix (rows=observables, cols=params)
        verdict           — human-readable string
    """
    S = compute_sensitivities(record, params, h_rel)

    # Normalise by measurement uncertainty
    x0 = np.array([_lit[p] for p in params])
    corr0 = _make_corr(params, x0)
    Q0, dp0, Tg0 = _run_solver(record, corr0)
    obs0 = np.array([dp0, Tg0 if record.T_g_out is not None else Q0])
    sigmas = np.array([
        record.sigma_dp() / dp0,
        (record.sigma_Tg() / Tg0 if record.T_g_out is not None
         else record.sigma_Q() / Q0),
    ])

    # Scaled sensitivity (inverse of σ): S_scaled_ij = S_ij / σ_i
    S_scaled = S / sigmas[:, np.newaxis]
    FIM = S_scaled.T @ S_scaled
    eigvals = np.linalg.eigvalsh(FIM)
    cond = float(eigvals.max() / (eigvals.min() + 1e-30))

    if cond < 100:
        verdict = "Well-conditioned — all parameters identifiable."
    elif cond < 1000:
        verdict = "Marginally conditioned — some parameter correlation; more data recommended."
    else:
        verdict = "Ill-conditioned — parameters nearly co-linear; reduce to fewer calibration params."

    return dict(
        condition_number=cond,
        eigenvalues=eigvals,
        sensitivity_matrix=S,
        params=list(params),
        verdict=verdict,
    )


# ---------------------------------------------------------------------------
# 1. Least-squares calibration
# ---------------------------------------------------------------------------

def calibrate_ls(
    records,
    params: Sequence[str] = ("salimpour_a", "ali_c_hi"),
    x0: np.ndarray | None = None,
) -> dict:
    """Levenberg-Marquardt weighted least-squares.

    `records` may be a single CalibrationRecord or a list of them.
    Each record contributes one residual per observable (dp_c and T_g_out or Q_He).

    Returns dict with keys: x, params, residuals, covariance, success, message.

    Covariance from J^T J (first-order approximation); std = sqrt(diag(cov)).
    """
    if isinstance(records, CalibrationRecord):
        records = [records]

    if x0 is None:
        x0 = np.array([_lit[p] for p in params])

    # Bounds from PRIOR_SPEC
    lower = np.array([PRIOR_SPEC[p]["bounds"][0] for p in params])
    upper = np.array([PRIOR_SPEC[p]["bounds"][1] for p in params])

    def residuals(x):
        corr = _make_corr(params, x)
        rs = []
        for rec in records:
            try:
                Q_m, dp_m, Tg_m = _run_solver(rec, corr)
            except Exception:
                rs.extend([1e6, 1e6])
                continue
            rs.append((dp_m - rec.dp_c) / rec.sigma_dp())
            if rec.T_g_out is not None:
                rs.append((Tg_m - rec.T_g_out) / rec.sigma_Tg())
            else:
                rs.append((Q_m - rec.Q_He_W) / rec.sigma_Q())
        return rs

    result = least_squares(residuals, x0, bounds=(lower, upper),
                           method='trf', ftol=1e-6, xtol=1e-6, gtol=1e-6, verbose=1)

    try:
        J = result.jac
        cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        cov = None

    return dict(
        x=result.x, params=list(params),
        residuals=result.fun, covariance=cov,
        success=result.success, message=result.message,
    )


# ---------------------------------------------------------------------------
# 2. Bayesian MAP calibration
# ---------------------------------------------------------------------------

def calibrate_map(
    records,
    params: Sequence[str] = ("salimpour_a", "ali_c_hi"),
    x0: np.ndarray | None = None,
    method: str = "Nelder-Mead",
) -> dict:
    """Bayesian maximum-a-posteriori estimate.

    Prior: physics-informed per-parameter distributions (PRIOR_SPEC).
    Likelihood: Gaussian with σ_total = sqrt(σ_meas² + σ_model²).

    `records` may be a single CalibrationRecord or a list.

    Returns dict with keys: x, params, neg_log_posterior, success, message.
    """
    if isinstance(records, CalibrationRecord):
        records = [records]

    if x0 is None:
        x0 = np.array([_lit[p] for p in params])

    def neg_lp(x):
        return -_log_posterior_multi(params, x, records)

    result = minimize(neg_lp, x0, method=method,
                      options=dict(xatol=1e-6, fatol=1e-6, maxiter=5000, disp=True))

    return dict(
        x=result.x, params=list(params),
        neg_log_posterior=float(result.fun),
        success=result.success, message=result.message,
    )


# ---------------------------------------------------------------------------
# 3. MCMC full posterior — requires: pip install emcee
# ---------------------------------------------------------------------------

def calibrate_mcmc(
    records,
    params: Sequence[str] = ("salimpour_a", "ali_c_hi"),
    nwalkers: int = 24,
    nsteps: int = 3000,
    burnin: int = 600,
    x0: np.ndarray | None = None,
) -> dict:
    """Full posterior sampling with emcee (ensemble MCMC).

    `records` may be a single CalibrationRecord or a list.

    Walkers are initialised in a small ball (5 % of prior σ) around the
    literature values (or x0 if provided).

    Convergence diagnostics:
      - acceptance fraction: healthy range 0.20 – 0.50
      - integrated autocorrelation time τ: want nsteps >> 50*τ
      - effective sample size ESS = (nsteps - burnin) * nwalkers / τ: want ESS > 200

    Returns dict with keys:
        flat          : (N_samples, ndim) posterior samples after burn-in
        samples       : (nsteps, nwalkers, ndim) full chain
        params        : parameter names
        mean, std     : posterior mean and std
        acceptance    : mean acceptance fraction
        tau           : integrated autocorrelation time (per parameter)
        ess           : effective sample size (per parameter)
        converged     : bool (nsteps > 50*max(tau))
    """
    try:
        import emcee
    except ImportError:
        raise ImportError("emcee not installed. Run: pip install emcee") from None

    if isinstance(records, CalibrationRecord):
        records = [records]

    ndim = len(params)
    if x0 is None:
        x0 = np.array([_lit[p] for p in params])

    # Initial ball: 5 % of prior σ around x0
    prior_scales = []
    for p in params:
        spec = PRIOR_SPEC[p]
        if spec["dist"] == "lognormal":
            prior_scales.append(x0[list(params).index(p)] * spec["sigma"] * 0.05)
        else:
            prior_scales.append(spec["sigma"] * 0.05)

    rng = np.random.default_rng(42)
    p0 = x0 + rng.normal(0, prior_scales, size=(nwalkers, ndim))

    def log_prob(x):
        return _log_posterior_multi(params, x, records)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
    sampler.run_mcmc(p0, nsteps, progress=True)

    flat = sampler.get_chain(discard=burnin, flat=True)

    try:
        tau = sampler.get_autocorr_time(discard=burnin, tol=10, quiet=True)
    except Exception:
        tau = np.full(ndim, float('nan'))

    ess = (nsteps - burnin) * nwalkers / np.where(np.isfinite(tau), tau, np.inf)
    converged = bool(np.all(np.isfinite(tau)) and (nsteps - burnin) > 50 * np.nanmax(tau))

    return dict(
        flat=flat,
        samples=sampler.get_chain(),
        params=list(params),
        mean=flat.mean(axis=0),
        std=flat.std(axis=0),
        acceptance=float(sampler.acceptance_fraction.mean()),
        tau=tau,
        ess=ess,
        converged=converged,
    )


# ---------------------------------------------------------------------------
# 4. Sequential Bayesian updating
# ---------------------------------------------------------------------------

def calibrate_sequential(
    records: list[CalibrationRecord],
    params: Sequence[str] = ("salimpour_a", "ali_c_hi"),
    nwalkers: int = 24,
    nsteps: int = 2000,
    burnin: int = 400,
) -> list[dict]:
    """Sequential Bayesian updating across multiple experimental records.

    For each experiment in `records`, a new MCMC run is performed using the
    posterior from the previous experiment as the new prior.  The prior is
    approximated as a Gaussian (mean, cov) fitted to the previous posterior
    samples — this is the Laplace approximation of the sequential update.

    This is appropriate when:
    - Records come from the same physical system (same θ_true)
    - Experiments are conducted at different operating conditions
    - Records accumulate over time (test campaign)

    Returns a list of MCMC result dicts, one per record (in order).
    Each dict has the same keys as calibrate_mcmc(), plus:
        'record_index' : index of the record used for this update.
    """
    try:
        import emcee
    except ImportError:
        raise ImportError("emcee not installed. Run: pip install emcee") from None

    ndim = len(params)
    results = []

    # Custom prior override for sequential update
    current_prior_mean = np.array([_lit[p] for p in params])
    current_prior_cov = None  # None → use PRIOR_SPEC for first record

    for idx, rec in enumerate(records):
        if current_prior_cov is None:
            def log_prob(x):
                return _log_posterior_single(params, x, rec)
        else:
            cov_inv = np.linalg.inv(current_prior_cov)

            def log_prob(x, _mu=current_prior_mean, _C_inv=cov_inv, _rec=rec):
                # Gaussian prior from previous posterior
                diff = x - _mu
                lp = -0.5 * float(diff @ _C_inv @ diff)
                lo = np.array([PRIOR_SPEC[p]["bounds"][0] for p in params])
                hi = np.array([PRIOR_SPEC[p]["bounds"][1] for p in params])
                if np.any(x < lo) or np.any(x > hi):
                    return -np.inf
                corr = _make_corr(params, x)
                try:
                    Q_m, dp_m, Tg_m = _run_solver(_rec, corr)
                except Exception:
                    return -np.inf
                return lp + _log_likelihood(_rec, Q_m, dp_m, Tg_m)

        rng = np.random.default_rng(idx + 42)
        p0 = current_prior_mean + rng.normal(
            0,
            np.sqrt(np.diag(current_prior_cov)) * 0.05 if current_prior_cov is not None
            else np.array([PRIOR_SPEC[p]["sigma"] * _lit[p] * 0.05 for p in params]),
            size=(nwalkers, ndim),
        )

        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
        sampler.run_mcmc(p0, nsteps, progress=True)
        flat = sampler.get_chain(discard=burnin, flat=True)

        # Update prior for next record: Gaussian fit to this posterior
        current_prior_mean = flat.mean(axis=0)
        current_prior_cov = np.cov(flat.T)

        res = dict(
            flat=flat,
            samples=sampler.get_chain(),
            params=list(params),
            mean=current_prior_mean.copy(),
            std=flat.std(axis=0),
            acceptance=float(sampler.acceptance_fraction.mean()),
            record_index=idx,
        )
        results.append(res)
        print(f"[calibrate_sequential] Record {idx+1}/{len(records)} done — "
              f"mean = {dict(zip(params, current_prior_mean))}")

    return results


# ---------------------------------------------------------------------------
# 5. Posterior predictive uncertainty propagation
# ---------------------------------------------------------------------------

def posterior_predictive(
    flat_samples: np.ndarray,
    records,
    params: Sequence[str],
) -> dict:
    """Propagate parameter posterior uncertainty to model output uncertainty.

    Draws N_samples parameter vectors from `flat_samples`, runs the solver
    for each, and returns the distribution of (Q_He, dp_c, T_g_out).

    `records` may be a single CalibrationRecord or a list.  If a list, the
    posterior predictive is computed for each record independently.

    Returns dict with keys per observable: mean, std, percentiles (5, 25, 50, 75, 95).
    """
    if isinstance(records, CalibrationRecord):
        records = [records]

    results_per_record = []
    for rec in records:
        Qs, dps, Tgs = [], [], []
        for x in flat_samples:
            corr = _make_corr(params, x)
            try:
                Q_m, dp_m, Tg_m = _run_solver(rec, corr)
                Qs.append(Q_m); dps.append(dp_m); Tgs.append(Tg_m)
            except Exception:
                pass
        arr = {
            "Q_He_W": np.array(Qs),
            "dp_c_Pa": np.array(dps),
            "T_g_out_K": np.array(Tgs),
        }
        summary = {}
        for key, a in arr.items():
            if len(a) == 0:
                continue
            summary[key] = dict(
                mean=float(a.mean()),
                std=float(a.std()),
                p5=float(np.percentile(a, 5)),
                p25=float(np.percentile(a, 25)),
                p50=float(np.percentile(a, 50)),
                p75=float(np.percentile(a, 75)),
                p95=float(np.percentile(a, 95)),
            )
        results_per_record.append(summary)

    return results_per_record if len(results_per_record) > 1 else results_per_record[0]


# ---------------------------------------------------------------------------
# 6. Posterior plotting (no corner dependency)
# ---------------------------------------------------------------------------

def plot_posteriors(result: dict, truths: dict | None = None):
    """Plot marginal posteriors and pairwise correlations.

    `result`  — dict returned by calibrate_mcmc() or calibrate_sequential()[-1].
    `truths`  — dict {param_name: true_value} for reference lines (optional).

    Produces two figures:
      Fig 1: marginal histograms + kernel density estimates, one panel per parameter.
      Fig 2: pairwise scatter / 2-D density, lower triangle only.
    """
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde

    flat   = result["flat"]
    params = result["params"]
    ndim   = len(params)
    lit    = {p: _lit[p] for p in params}

    # --- Figure 1: marginals ---
    ncols = min(ndim, 4)
    nrows = math.ceil(ndim / ncols)
    fig1, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows),
                              constrained_layout=True)
    axes = np.array(axes).ravel()
    for i, p in enumerate(params):
        ax = axes[i]
        samples_i = flat[:, i]
        ax.hist(samples_i, bins=50, density=True, alpha=0.6, label="posterior")
        try:
            kde = gaussian_kde(samples_i)
            xs = np.linspace(samples_i.min(), samples_i.max(), 300)
            ax.plot(xs, kde(xs), 'k-', lw=1.5)
        except Exception:
            pass
        ax.axvline(lit[p], color="C1", ls="--", lw=1.2, label="literature")
        if truths and p in truths:
            ax.axvline(truths[p], color="C2", ls=":", lw=1.5, label="true")
        spec = PRIOR_SPEC[p]
        ax.set_title(
            f"{p}\nμ={flat[:, i].mean():.4g}  σ={flat[:, i].std():.3g}",
            fontsize=9,
        )
        ax.set_xlabel(p, fontsize=8)
        ax.legend(fontsize=7)
    for ax in axes[ndim:]:
        ax.set_visible(False)
    fig1.suptitle("Posterior marginal distributions", fontsize=11)

    # --- Figure 2: pairwise correlations (lower triangle) ---
    if ndim > 1:
        fig2, axes2 = plt.subplots(ndim, ndim, figsize=(2.5 * ndim, 2.5 * ndim),
                                   constrained_layout=True)
        for i in range(ndim):
            for j in range(ndim):
                ax = axes2[i, j]
                if i == j:
                    ax.hist(flat[:, i], bins=40, density=True, alpha=0.7)
                    ax.axvline(lit[params[i]], color="C1", ls="--", lw=1)
                    ax.set_xlabel(params[i], fontsize=7)
                elif i > j:
                    try:
                        xy = np.vstack([flat[:, j], flat[:, i]])
                        kde2 = gaussian_kde(xy)
                        xg = np.linspace(flat[:, j].min(), flat[:, j].max(), 50)
                        yg = np.linspace(flat[:, i].min(), flat[:, i].max(), 50)
                        Xg, Yg = np.meshgrid(xg, yg)
                        Z = kde2(np.vstack([Xg.ravel(), Yg.ravel()])).reshape(50, 50)
                        ax.contourf(Xg, Yg, Z, levels=10, cmap="Blues")
                    except Exception:
                        ax.scatter(flat[:, j], flat[:, i], s=1, alpha=0.3)
                    ax.set_xlabel(params[j], fontsize=7)
                    ax.set_ylabel(params[i], fontsize=7)
                else:
                    ax.set_visible(False)
        fig2.suptitle("Pairwise posterior correlations", fontsize=11)

    plt.show()
    return fig1


# ---------------------------------------------------------------------------
# Quick demo (run as script)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import pprint

    # Synthetic record at default input_data conditions.
    # Replace with real test-bench measurements.
    rec = CalibrationRecord(
        m_dot_c=0.15, T_He_in=120.0, T_He_out=750.0,
        p_He_in=90e5, p_He_out=73e5,
        m_dot_g=0.070, OF=2.9,
        # T_g_out=950.0,   # uncomment when T_g_out is measured
    )

    print(f"Record: Q_He = {rec.Q_He_W/1e3:.2f} kW,  dp_c = {rec.dp_c/1e5:.2f} bar")
    print(f"σ_dp = {rec.sigma_dp()/rec.dp_c*100:.1f}%,  σ_Q = {rec.sigma_Q()/rec.Q_He_W*100:.1f}%")
    print()

    print("=== Sensitivity analysis ===")
    S = compute_sensitivities(rec, params=["salimpour_a", "ali_c_hi"])
    print(f"  S[dp, salimpour_a]  = {S[0,0]:+.3f}  (should be ~0 in counter-flow saturation regime)")
    print(f"  S[dp, ali_c_hi]     = {S[0,1]:+.3f}  (should be ~0.6)")
    print(f"  S[obs2, salimpour_a] = {S[1,0]:+.3f}")
    print(f"  S[obs2, ali_c_hi]   = {S[1,1]:+.3f}")
    print()

    print("=== Identifiability check ===")
    chk = identifiability_check(rec, params=["salimpour_a", "ali_c_hi"])
    pprint.pprint({k: v for k, v in chk.items() if k != "sensitivity_matrix"})
    print()

    print("=== Least-squares calibration ===")
    ls = calibrate_ls(rec, params=["salimpour_a", "ali_c_hi"])
    for name, val in zip(ls["params"], ls["x"]):
        print(f"  {name:20s} = {val:.4f}  (literature: {_lit[name]:.4f})")
    print(f"  Success: {ls['success']}")
    print()

    print("=== Bayesian MAP calibration ===")
    bmap = calibrate_map(rec, params=["salimpour_a", "ali_c_hi"])
    for name, val in zip(bmap["params"], bmap["x"]):
        print(f"  {name:20s} = {val:.4f}  (literature: {_lit[name]:.4f})")
    print(f"  Success: {bmap['success']}")
    print()

    # Uncomment to run MCMC (slow, requires: pip install emcee):
    # print("=== MCMC full posterior ===")
    # mcmc = calibrate_mcmc(rec, params=["salimpour_a", "ali_c_hi"],
    #                        nwalkers=24, nsteps=2000, burnin=400)
    # for name, mean, std in zip(mcmc["params"], mcmc["mean"], mcmc["std"]):
    #     print(f"  {name:20s} = {mean:.4f} ± {std:.4f}  (literature: {_lit[name]:.4f})")
    # print(f"  Acceptance: {mcmc['acceptance']:.3f},  converged: {mcmc['converged']}")
    # print(f"  ESS: {dict(zip(mcmc['params'], mcmc['ess'].astype(int)))}")
    # plot_posteriors(mcmc)
